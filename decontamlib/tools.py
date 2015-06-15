import collections
import csv
import inspect
import os
import random
import re
import shutil
import subprocess
import tempfile
import utils


from decontamlib.fastq import FastqSplitter
from decontamlib.sam import get_mapped_reads


def run_command(command, error_message):
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print e.stderr


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

    def decontaminate(self, fwd_fp, rev_fp, output_dir):
        annotations = self.annotate(fwd_fp, rev_fp)
        with FastqSplitter(fwd_fp, output_dir) as s:
            s.partition(annotations)
        with FastqSplitter(rev_fp, output_dir) as s:
            s.partition(annotations)
        summary_data = summarize_annotations(annotations)
        return summary_data

    def _get_mapped_reads(self, filename):
        """Extracts set of qnames from SAM file."""
        mapped = set()
        for qname, is_read1, rname in get_mapped_reads(filename):
            if rname is not None:
                mapped.add(qname)
        return mapped


class Bwa(_FilteringTool):
    def __init__(self, index, bwa_fp):
        self.index = index
        self.bwa_fp = bwa_fp

    def annotate(self, R1, R2):
        sam_file, stderr_file = self._run(R1, R2)
        mapped = self._get_mapped_reads(sam_file.name)
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]

    def _command(self, fwd_fp, rev_fp):
        return [self.bwa_fp, "mem", "-M", self.index, fwd_fp, rev_fp]

    def _run(self, R1, R2):
        stdout_file = tempfile.NamedTemporaryFile()
        stderr_file = tempfile.NamedTemporaryFile()
        command = self._command(R1, R2)
        subprocess.check_call(command, stdout=stdout_file, stderr=stderr_file)
        return (stdout_file, stderr_file)


class Bowtie(Bwa):
    def __init__(self, index, bowtie2_fp):
        self.index = index
        self.bowtie2_fp = bowtie2_fp

    def _command(self, fwd_fp, rev_fp):
        return [
            self.bowtie2_fp, "--local", "--very-sensitive-local",
            "-1", fwd_fp, "-2", rev_fp,
            "-x", self.index]


class Random_human(_FilteringTool):
    def __init__(self, percent_human):
        assert 0.0 <= percent_human <= 100.0
        self.fraction_human = percent_human / 100.0

    def annotate(self, R1, R2):
        ids = utils.parse_read_ids(R1)
        return [
            (id, True if random.random() <= self.fraction_human else False)
            for id in ids]


class None_human(_FilteringTool):
    def __init__(self):
        pass

    def annotate(self, R1, R2):
        ids = utils.parse_read_ids(R1)
        return [(id, False) for id in ids]


class All_human(_FilteringTool):
    def __init__(self):
        pass

    def annotate(self, R1, R2):
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
}
