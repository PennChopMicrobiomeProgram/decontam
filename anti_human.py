#!/usr/bin/python

import argparse
import os.path
import sys

def check_file_exists_or_die(file_name):
    if not os.path.isfile(file_name): 
        print "ERROR: file " + file_name + " does not exist: check the file name and path."
        sys.exit(1)

def command_line_arguments():
    """creates command-line args.
    returns arguments as attributes
    """
    parser = argparse.ArgumentParser(description="Compares performance of human filtering tools.")
    parser.add_argument("-s", "--samples", required=True, type=str, help="list of samples and corresponding annotation files")
    parser.add_argument("-t", "--tools", required=True, type=str, help="human filtering tools")
    parser.add_argument("-o", "--output", required=True, help="comparison of tools to expected results.")
    args = parser.parse_args()
    check_file_exists_or_die(args.samples) 
    check_file_exists_or_die(args.tools) 
    return args

if __name__=="__main__":
    args = command_line_arguments()

    tool_names = parse_tool_names(open(args.tools)) 
    tools = Tools.create(tool_names)

    samples = parse_samples(args.samples)

    results = []

    for sample in samples:
        sample_name = sample[0]
        fastq_file = sample[1]
        annotation_file = sample[2]

        annotation_file_handle = open(annotation_file)
        annotation = parse_annotation(annotation_file_handle)
        fastq_file_handle = open(fastq_file)
        reads = parse_reads(fastq_file_handle)
        check_all_reads_are_annotated(annotation, reads)
        for tool in tools:
            predicted_human = tool.get_human_reads(sample_name)
            tool_name = tool.get_name()
#            quality_metrics = check_tool_vs_oracle(predicted_human, annotation)
#            results.append( (sample_name, tool_name, quality_metrics) )
#
#    write_results(results)
             
 
    
