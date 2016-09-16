import collections
import csv
import inspect
import os
import random
import re
import shutil
import subprocess
import tempfile

from decontamlib import utils
from decontamlib.fastq import FastqSplitter
from decontamlib.sam import get_mapped_reads


def FilteringTool(config):
    tool_cls = tools_available[config["method"]]
    # Proceed stepwise here to improve quality of error messages.
    tool_args = []
    for argname in tool_cls.get_argnames():
        arg = config[argname]
        tool_args.append(arg)
    return tool_cls(*tool_args)


class _FilteringTool(object):
    def __init__(self, index):
        self.index = index

    @classmethod
    def get_argnames(cls):
        return inspect.getargspec(cls.__init__)[0][1:]

    def decontaminate(self, fwd_fp, rev_fp, output_dir, organism, pct, frac):
        annotations = self.annotate(fwd_fp, rev_fp, pct, frac, output_dir)
        with FastqSplitter(fwd_fp, output_dir) as s:
            s.partition(annotations, organism)
        with FastqSplitter(rev_fp, output_dir) as s:
            s.partition(annotations, organism)
        summary_data = summarize_annotations(annotations)
        return summary_data

    def annotate(self, fwd_fp, rev_fp, pct, frac, output_dir):
        raise NotImplementedError()

    def _get_mapped_reads(self, filename, pct, frac):
        """Extracts set of qnames from SAM file."""
        mapped = set()
        for qname, is_read1, rname in get_mapped_reads(filename, pct, frac):
            if rname is not None:
                mapped.add(qname)
        return mapped

    def make_index(self):
        raise NotImplementedError()

    def index_exists(self):
        return True


class SamFile(_FilteringTool):
    def __init__(self, sam_fp):
        self.sam_fp = sam_fp

    def annotate(self, R1, R2, pct, frac, output_dir):
        mapped = self._get_mapped_reads(self.sam_fp, pct, frac)
        ids = utils.parse_read_ids(R1)
        return [(id, True if id in mapped else False) for id in ids]


class Bwa(_FilteringTool):
    def __init__(self, index, bwa_fp, num_threads, keep_sam_file):
        self.index = index
        self.bwa_fp = bwa_fp
        self.num_threads = num_threads
        self.keep_sam_file = keep_sam_file
        
    def annotate(self, R1, R2, pct, frac, output_dir):
        sam_file, stderr_file = self._run(R1, R2, output_dir)
        mapped = self._get_mapped_reads(sam_file.name, pct, frac)
        ids = utils.parse_read_ids(R1)
        return [(id, True if id in mapped else False) for id in ids]

    def _command(self, fwd_fp, rev_fp):
        return [self.bwa_fp, "mem", "-M", "-t", str(self.num_threads), self.index, fwd_fp, rev_fp]

    def _run(self, R1, R2, output_dir):
        if self.keep_sam_file:
            fwd_base_filename = os.path.splitext(os.path.basename(R1))[0]
            stdout_fp = os.path.join(output_dir, fwd_base_filename + ".sam")
            stdout_file = open(stdout_fp, "w")
        else:
            stdout_file = tempfile.NamedTemporaryFile()
        stderr_file = tempfile.NamedTemporaryFile()
        command = self._command(R1, R2)
        subprocess.check_call(command, stdout=stdout_file, stderr=stderr_file)
        return (stdout_file, stderr_file)

    def make_index(self):
        cmd = [self.bwa_fp, "index", self.index]
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT)

    def index_exists(self):
        return os.path.exists(self.index + ".amb" )


class Bowtie(Bwa):
    def __init__(self, index, bowtie2_fp, keep_sam_file):
        self.index = index
        self.bowtie2_fp = bowtie2_fp
        self.keep_sam_file = keep_sam_file

    def _command(self, fwd_fp, rev_fp):
        return [
            self.bowtie2_fp, "--local", "--very-sensitive-local",
            "-1", fwd_fp, "-2", rev_fp,
            "-x", self.index]

    def make_index(self):
        build_fp = self.bowtie2_fp + "-build"
        cmd = [build_fp, "-f", self.index, self.index]
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT)

    def index_exists(self):
        index_fp = self.index + ".1.bt2"
        return os.path.exists(index_fp)


class Random_human(_FilteringTool):
    def __init__(self, percent_human):
        assert 0.0 <= percent_human <= 100.0
        self.fraction_human = percent_human / 100.0

    def annotate(self, R1, R2, pct, frac, output_dir):
        ids = utils.parse_read_ids(R1)
        return [
            (id, True if random.random() <= self.fraction_human else False)
            for id in ids]


class None_human(_FilteringTool):
    def __init__(self):
        pass

    def annotate(self, R1, R2, pct, frac, output_dir):
        ids = utils.parse_read_ids(R1)
        return [(id, False) for id in ids]


class All_human(_FilteringTool):
    def __init__(self):
        pass

    def annotate(self, R1, R2, pct, frac, output_dir):
        ids = utils.parse_read_ids(R1)
        return [(id, True) for id in ids]


def summarize_annotations(annotations):
    return dict(collections.Counter(a for _, a in annotations))


tools_available = {
    "bwa": Bwa,
    "all_human": All_human,
    "no_human": None_human,
    "random_human": Random_human,
    "bowtie2": Bowtie,
    "samfile": SamFile,
}
