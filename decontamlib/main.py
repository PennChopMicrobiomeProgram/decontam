import argparse
import collections
import csv
import itertools
import json
import os.path
import subprocess
import sys

from decontamlib.version import __version__
from decontamlib.tools import FilteringTool


default_config = {
    "method": "bowtie2",
    "bowtie2_fp": "/home/kyle/software/bowtie2-2.2.5/bowtie2",
    "bwa_fp": "/home/kyle/software/bwa-0.7.12/bwa",
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
