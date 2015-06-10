import argparse
import csv
import itertools
import os.path
import subprocess
import sys

import .parser
import .tools
import .utils


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


def human_filter_main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run human filtering tools.")
    parser.add_argument(
        "-s", "--samples", required=True,
        help="list of paired-end samples.")
    parser.add_argument(
        "-t", "--tools", required=True,
        help="human filtering tools.")
    parser.add_argument(
        "-p", "--parameters_for_tools",
        help=(
            "parameters for different tools (by default search "
            "parameters.json file in current folder.)"))
    parser.add_argument(
        "-o", "--output", default="result.dat",
        help="long table of results.")
    parser.add_argument(
        "-path", required=True,
        help="Path to output directory.")
    args = parser.parse_args(argv)

    tool_names = parser.parse_tool_names(open(args.tools)) 
    if args.parameters_for_tools:
        tool_parameters = tools.get_parameters_for_tools(args.parameters_for_tools)
    else:
        tool_parameters = tools.get_parameters_for_tools()
    tools = tools.create_tools(tool_names, tool_parameters)

    samples = parser.parse_samples(open(args.samples))

    results = []

    for sample in samples:
        (sample_name, R1_fastq_file, R2_fastq_file) = sample

        R1_read_ids = utils.parse_read_ids(R1_fastq_file)
        R2_read_ids = utils.parse_read_ids(R2_fastq_file)
        if not utils.check_all_read_ids_are_consistent(R1_read_ids, R2_read_ids):
            print "R1 and R2 file should have the same ids."
            print "R1  file  " + R1_fastq_file + " R2 file " + R2_fastq_file + " are not consistent."
            sys.exit(1)

        for tool in tools:
            human_annotation = tool.get_human_annotation(R1_fastq_file, R2_fastq_file)
            tool_name = tool.name
            results_for_tool_sample = utils.add_tool_sample(tool_name, sample_name, human_annotation)
            results += results_for_tool_sample
            filter_human_from_fastq(results_for_tool_sample, sample, args.path)
    write_results(args.path, args.output, results)


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
