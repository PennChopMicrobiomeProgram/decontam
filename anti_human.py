#!/usr/bin/python

import argparse
import os.path
import sys
import csv
import itertools
import parser
import tools
import utils

def check_file_exists_or_die(file_name):
    if not os.path.isfile(file_name): 
        print "ERROR: file " + file_name + " does not exist: check the file name and path."
        sys.exit(1)

def command_line_arguments():
    """creates command-line args.
    Returns:
        arguments as attributes
    """
    parser = argparse.ArgumentParser(description="Run human filtering tools.")
    parser.add_argument("-s", "--samples", required=True, type=str, help="list of paired-end samples.")
    parser.add_argument("-t", "--tools", required=True, type=str, help="human filtering tools.")
    parser.add_argument("-p", "--parameters_for_tools", 
        help="parameters for different tools(by default search parameters.json file in current folder.)")
    parser.add_argument("-o", "--output", type=str, default="result.dat", help="long table of results.")
    parser.add_argument("-path", required=True, type=str, help="Path to output directory.")

    args = parser.parse_args()
    check_file_exists_or_die(args.samples) 
    check_file_exists_or_die(args.tools) 
    return args

def write_results(path,filename, results):
    """ write tool, sample, read_id, is_human to tab-separated file.
    """
    writer = csv.writer(open(path + filename, 'w'), delimiter="\t")
    writer.writerow(["tool", "sample", "read_id", "is_human"])
    writer.writerows(results)                  

def get_non_human_read_ids(results):
    r_id = [] 
    for result in results:
        (tool_name, name_sample, read_id, is_human) = result
        if not is_human:
            r_id.append(read_id)

    return sorted(r_id)

'''
def write_filtered_reads_to_fastq(fastq_file, r_id, tool_name, sample_name, is_r1, path):

    filtered_reads = []
    flag = 0
    with open(fastq_file, "r") as sample:
        for line in sample:
            if line.startswith("@"):
                line = line.rstrip()
                if str(line[1:]) in r_id:
                    flag = 1
                else:
                    flag = 0

            if flag:
                filtered_reads.append(line.rstrip())

    if is_r1:
        fname = path + tool_name + "_" + sample_name + "-R1.fastq"
        print fname
    else:
        fname = path + tool_name + "_" + sample_name + "-R2.fastq"
        
    with open(fname, "w") as filter:
        for read in filtered_reads:
            filter.write(read)
            filter.write("\n")
'''

def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF
    args = [iter(iterable)] * n
    return itertools.izip(*args)


def parse_fastq(f):
    seq_by_id = {} 
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        #yield desc, seq, qual
        seq_by_id[desc] = (seq,qual)
    return seq_by_id



def filter_human_from_fastq(results, sample, path):
    (tool_name, name_sample, read_id, is_human) = results[0]
    (sample_name, R1_fastq_file, R2_fastq_file) = sample
    r_id = get_non_human_read_ids(results)
    #write_filtered_reads_to_fastq(R1_fastq_file, r_id, tool_name, sample_name, 1, path)    
    #write_filtered_reads_to_fastq(R2_fastq_file, r_id, tool_name, sample_name, 0, path)
    r1 = open(R1_fastq_file, "r")
    r2 = open(R2_fastq_file, "r")
    fname_r1 = path + tool_name + "_" + sample_name + "-R1.fastq"
    fname_r2 = path + tool_name + "_" + sample_name + "-R2.fastq"
    r1_seq = parse_fastq(r1)
    r2_seq = parse_fastq(r2)
    write_filtered_to_fastq(r1_seq, r_id, fname_r1)
    write_filtered_to_fastq(r2_seq, r_id, fname_r2)    

def write_filtered_to_fastq(r, r_id, fname):
    with open(fname, "w") as fastq:
        for key in sorted(r.keys()):
            if key in r_id:
                fastq.write("@" + key + "\n")
                fastq.write(r[key][0] + "\n")
                fastq.write("+\n")
                fastq.write(r[key][1] + "\n")
        


if __name__=="__main__":
    args = command_line_arguments()

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
