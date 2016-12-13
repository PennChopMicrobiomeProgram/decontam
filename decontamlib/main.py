
import argparse
import collections
import csv
import itertools
import json
import os.path
import subprocess
import sys

from decontamlib.version import __version__
from decontamlib.tools import FilteringTool, tools_available


# def get_config(user_config_file, organism):
#     config = {
#         "method": "bwa",
#         "bowtie2_fp": "bowtie2",
#         "bwa_fp":"bwa",
#         "num_threads":8
#     }

#     if user_config_file is None:
#         if organism == "human":
#             default_user_config_fp = os.path.expanduser("~/.decontam_human.json")
#         elif organism == "phix":
#             default_user_config_fp = os.path.expanduser("~/.decontam_phix.json")
#         if os.path.exists(default_user_config_fp):
#             user_config_file = open(default_user_config_fp)

#     if user_config_file is not None:
#         user_config = json.load(user_config_file)
#         config.update(user_config)
#     return config


def get_config(args, organism):
    config = {
        "method":args.method,
        "bowtie2_fp":args.bowtie2_fp,
        "bwa_fp":args.bwa_fp,
        "num_threads":args.threads,
        "index":args.index
    }
    return config

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
        "--organism", required=True,
        help="reference organism to filter from",
        choices=['human', 'phix'])
    p.add_argument(
        "--summary-file", required=True,
        type=argparse.FileType("w"),
        help="long table of results.")
    p.add_argument(
        "--output-dir", required=True,
        help="Path to output directory")
    p.add_argument(
        "--method", default="bwa",
        help="Method to use (default %(default)s)",
        choices=tools_available.keys())
    p.add_argument(
        "--index", 
        help="Path to contaminant genome indexes")
    p.add_argument(
        "--threads", default=1,
        type=int, help="Number of threads to use")                   
    p.add_argument(
        "--pct", required=False,
        type=float, default=0.5,
        help="Percent identity (default: %(default)s)")
    p.add_argument(
        "--frac", required=False,
        type=float, default=0.6,
        help="Fraction of alignment length (default: %(default)s)")
    p.add_argument(
        "--keep-sam-file", action="store_true",
        help="Write SAM file to output directory.")
    p.add_argument(
        "--sam-file", required=False,
        type=argparse.FileType("r"),
        help="File of alignments to reference database (SAM format, optional)")
    p.add_argument(
        "--bowtie2_fp", default="bowtie2",
        help="Path to bowtie2 (if using)")
    p.add_argument(
        "--bwa_fp", default="bwa",
        help="Path to bwa (if using)")

    # Output
    args = p.parse_args(argv)

    config = get_config(args, args.organism)

    fwd_fp = args.forward_reads.name
    rev_fp = args.reverse_reads.name
    args.forward_reads.close()
    args.reverse_reads.close()

    if args.sam_file is not None:
        config["method"] = "samfile"
        config["sam_fp"] = args.sam_file.name

    if args.keep_sam_file is True:
        config["keep_sam_file"] = True
    else:
        config["keep_sam_file"] = False

    tool = FilteringTool(config)
    if not tool.index_exists():
        tool.make_index()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    summary_data = tool.decontaminate(fwd_fp, rev_fp, args.output_dir,
                                      args.organism, args.pct, args.frac)
    save_summary(args.summary_file, config, summary_data)


def save_summary(f, config, data):
    result = {
        "program": "decontam",
        "version": __version__,
        "config": config,
        "data": data,
        }
    json.dump(result, f)


def make_index_main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--config-file", required=True,
        type=argparse.FileType("r"),
        help="JSON configuration file")
    p.add_argument(
        "--organism", required=True,
        help="reference organism to filter from")
    args = p.parse_args(argv)
    config = get_config(args.config_file, args.organism)
    
    tool = FilteringTool(config)

    if not tool.index_exists():
        tool.make_index()
