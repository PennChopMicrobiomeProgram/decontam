import argparse
import collections
import csv
import itertools
import json
import os.path
import subprocess
import sys

from decontamlib.version import __version__
from decontamlib.human_filtering_tools import FilteringTool


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
    "bowtie2_fp": "/home/kyle/software/bowtie2-2.2.5/bowtie2",
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

    tool = FilteringTool(config)

    if os.path.exists(args.output_dir):
        p.error("Output directory already exists")
    os.mkdir(args.output_dir)

    summary_data = tool.decontaminate(fwd_fp, rev_fp, args.output_dir)
    save_summary(args.summary_file, config, summary_data)


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
