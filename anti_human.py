#!/usr/bin/python

import argparse
import os.path
import sys
import csv

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
    parser.add_argument("-o", "--output", required=True, help="long table of results.")
    args = parser.parse_args()
    check_file_exists_or_die(args.samples) 
    check_file_exists_or_die(args.tools) 
    return args

def write_results(filename, results):
    """ write tool, sample, read_id, is_human to tab-separated file.
    """
    writer = csv.writer(open(filename, 'w'), delimiter="\t")
    writer.writerow(["tool", "sample", "read_id", "is_human"])
    writer.writerows(results)
  

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
    write_results(args.output, results)
             
 
    
