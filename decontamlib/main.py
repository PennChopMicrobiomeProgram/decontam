import argparse
import collections
import csv
import itertools
import json
import os.path
import subprocess
import sys

from decontamlib.version import __version__
from decontamlib.parser import (
    parse_tool_names,
    )
from decontamlib.human_filtering_tools import (
    Snap, Blat, Bwa, Bmfilter, Bmtagger,
    All_human, None_human, Random_human,
    Bowtie,
    )


tools_available = {
    "snap": Snap,
    "blat": Blat,
    "bwa": Bwa,
    "bmfilter": Bmfilter,
    "bmtagger": Bmtagger,
    "all_human": All_human,
    "no_human": None_human,
    "random_human": Random_human,
    "bowtie2": Bowtie,
}


def write_results(path, filename, results):
    """Write results to tab-separated file."""
    writer = csv.writer(open(path + filename, 'w'), delimiter="\t")
    writer.writerow(["tool", "sample", "read_id", "is_human"])
    writer.writerows(results)                  


def get_non_human_read_ids(results):
    r_id = set()
    for result in results:
        (tool_name, name_sample, read_id, is_human) = result
        if not is_human:
            r_id.add(read_id)
    return r_id


def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return itertools.izip(*args)


def parse_fastq(f):
    """ parse original fastq file and write new fastq file the filtered non-human reads.
    """
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        yield desc, seq, qual


def write_fastq(out_fastq, desc, seq, qual):
    out_fastq.write("@" + desc + "\n")
    out_fastq.write(seq + "\n")
    out_fastq.write("+\n")
    out_fastq.write(qual + "\n")
             

def filter_fastq(in_fastq, out_fastq, r_id):
    reads = parse_fastq(in_fastq)
    for desc, seq, qual in reads:
        if desc.split(" ")[0] in r_id:
            write_fastq(out_fastq, desc, seq, qual)
    in_fastq.close()
    out_fastq.close()


def filter_human_from_fastq(results, sample, path):
    """ Get non-human read ids and filter fastq file for non-human reads.
    """
    (tool_name, name_sample, read_id, is_human) = results[0]
    (sample_name, R1_fastq_file, R2_fastq_file) = sample

    #get non-human read ids.
    r_id = get_non_human_read_ids(results)
    fname_r1 = path + tool_name + "_" + sample_name + "-R1.fastq"
    fname_r2 = path + tool_name + "_" + sample_name + "-R2.fastq"
    filter_fastq(open(R1_fastq_file), open(fname_r1, "w"), r_id)
    filter_fastq(open(R2_fastq_file), open(fname_r2, "w"), r_id)


default_config = {
    "method": "bowtie2",
    "bowtie2_fp": "/home/kyle/software/bowtie2-2.2.5/bowtie2"
    }

def human_filter_main(argv=None):
    p = argparse.ArgumentParser()
    # Input
    p.add_argument(
        "--forward-reads", required=True,
        type=argparse.FileType("r"),
        help="FASTQ file of forward reads")
    p.add_argument(
        "--reverse-reads", required=True,
        type=argparse.FileType("r"),
        help="FASTQ file of reverse reads")
    p.add_argument(
        "--config-file",
        type=argparse.FileType("r"),
        help="JSON configuration file")
    # Output
    p.add_argument(
        "--summary-file", required=True,
        type=argparse.FileType("w"),
        help="long table of results.")
    p.add_argument(
        "--output-dir", required=True,
        help="Path to output directory")
    args = p.parse_args(argv)

    config = default_config.copy()
    if args.config_file:
        user_config = json.load(args.config_file)
        config.update(user_config)

    fwd_fp = args.forward_reads.name
    rev_fp = args.reverse_reads.name
    args.forward_reads.close()
    args.reverse_reads.close()

    tool_cls = tools_available[config["method"]]
    # Proceed stepwise here to improve quality of error messages.
    tool_args = []
    for argname in tool_cls.get_argnames():
        arg = config[argname]
        tool_args.append(arg)
    tool = tool_cls(*tool_args)

    if os.path.exists(args.output_dir):
        p.error("Output directory already exists")
    os.mkdir(args.output_dir)

    annotations = tool.annotate(fwd_fp, rev_fp)

    with open(fwd_fp) as f:
        partition_fastq(f, annotations, args.output_dir)
    with open(rev_fp) as f:
        partition_fastq(f, annotations, args.output_dir)
    summary_data = summarize_annotations(annotations)
    save_summary(args.summary_file, config, summary_data)


class SplitFastqWriter(object):
    suffixes = {
        True: "_human",
        False: "",
        }

    def __init__(self, input_fp, output_dir):
        self.output_dir = output_dir
        input_filename = os.path.basename(input_fp)
        self.input_root, self.input_ext = os.path.splitext(input_filename)
        self._open_files = {}

    def write(self, read, annotation):
        desc, seq, qual = read
        suffix = self.suffixes[annotation]
        output_filename = self.input_root + suffix + self.input_ext
        fp = os.path.join(self.output_dir, output_filename)
        if fp not in self._open_files:
            self._open_files[fp] = open(fp, "w")
        f = self._open_files[fp]
        write_fastq(f, desc, seq, qual)

    def close(self):
        for f in self._open_files.values():
            f.close()


def partition_fastq(f, annotations, output_dir):
    open_files = {}
    ids_to_annotation = dict(annotations)
    input_filename = os.path.basename(f.name)
    writer = SplitFastqWriter(f.name, output_dir)
    for desc, seq, qual in parse_fastq(f):
        read_id = desc.split(" ")[0]
        annotation = ids_to_annotation[read_id]
        writer.write((desc, seq, qual), annotation)
    writer.close()


def summarize_annotations(annotations):
    return dict(collections.Counter(a for _, a in annotations))


def save_summary(f, config, data):
    result = {
        "program": "decontam",
        "version": __version__,
        "config": config,
        "data": data,
        }
    json.dump(result, f)


def command_line_arguments(argv):
    parser = argparse.ArgumentParser(
        description = "Makes index files for Bowtie2, Bwa, Blat and BMTagger.")
    parser.add_argument(
        "-g", required=True, type=str,
        help="Path to the genome file")
    parser.add_argument(
        "-o", required=True, type=str,
        help="Name of the organism")
    parser.add_argument(
        "-bowtie", required=False, type=int, default=0,
        help="Create index for bowtie 1 or 0. Default:0")
    parser.add_argument(
        "-bwa", required=False, type=int, default=0,
        help="Create index for bwa 1 or 0. Default:0")
    parser.add_argument(
        "-blat", required=False, type=int, default=0,
        help="Create index for blat 1 or 0. Default:0")
    parser.add_argument(
        "-bmtagger", required=False, type=int, default=0,
        help="Create index for bmtagger 1 or 0. Default:0")
    args = parser.parse_args(argv)
    return args


def run_command(command, error_message):
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError:
        print error_message


def make_index_bowtie(genome, organism):
    command = ("bowtie2-build -f " + genome + " " + organism + ".fasta")
    run_command(command, "cannot run bowtie2. Check path to genome file.")


def make_index_bmtagger(genome, organism):
    command = ("bmtool -d " + genome + " -o " + organism + ".bitmask -A 0 -w 18")
    run_command(command, "cannot run bmtool. Check path to genome file.")


def make_index_bwa(genome, organism):
    command = ("bwa index " + genome)
    run_command(command, "cannot run bwa index. Check path to genome file.")


def make_index_blat(genome, organism):
    command = (
        "blat " + genome + " /dev/null /dev/null -tileSize=11 -makeOoc=" +
        organism + ".ooc repMatch=1024")
    run_command(command, "cannot run blat. Check path to genome file.")


def make_index_main(argv=None):
    args = command_line_arguments(argv)
    if args.bowtie:
        make_index_bowtie(args.g, args.o)
    if args.bwa:
        make_index_bwa(args.g, args.o)
    if args.bmtagger:
        make_index_bmtagger(args.g, args.o)
    if args.blat:
        make_index_blat(args.g, args.o)
