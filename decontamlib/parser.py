""" parsers for input files
"""

import re
import csv
import os

def parse_tool_names(tool_file_handle):
    """ parse file with list of tools, each tool should be on its own line:
        tool1
        tool2
        ...
        
        Returns:
             list with tool names
        Raises:
             IOError if cannot parse file
    """
    tools = []
    for line in tool_file_handle:
        tool = str.strip(line)
        tools.append(tool)
    if not tools:
        raise IOError("empty tool file")
    return tools


def check_file_exists(label, file_name):
    if not os.path.isfile(file_name): 
        raise IOError(label + "file :" + file_name + 
            "does not exist. Check file name or remove it from sample file.")


def check_number_of_rows(row, num_description_for_sample):
    if len(row) != num_description_for_sample:
        raise IOError("each sample should have " + str(num_description_for_sample) +
            " columns: sample_name R1_fastq R2_fastq") 


def parse_samples(samples_file_handle):
    """parses file with list of samples. 

    Each test case has name, R1 FASTQ and R2 FASTQ.

    Returns:
        list of tuples, each element of the list is a sample, 
        each sample is a tuple
        (1st element sample name, 2nd R1 fastq file, 3rd R2 fastq file)

    Raises: 
        IOError if cannot parse file
        IOError if cannot find fastq file
    """
    num_description_for_sample = 3 # name, forward(R1) fastq, reverse(R2) fastq
    samples = []
    reader = csv.reader(samples_file_handle, delimiter="\t")
    for row in reader:
        check_number_of_rows(row, num_description_for_sample)
        (sample_name, r1_fastq_file, r2_fastq_file) = row
        check_file_exists("FASTQ", r1_fastq_file)
        check_file_exists("FASTQ", r2_fastq_file)
        samples.append( (sample_name, r1_fastq_file, r2_fastq_file) )
    if not samples:
        raise IOError("empty sample file") 
    return samples
