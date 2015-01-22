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
        
        return list with tool names
        raises IOError if cannot parse file
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
            " columns: sample_name fastq annotation") 

def parse_samples(samples_file_handle):
    """ parses file with list of test cases. 

    Each test case has name, FASTQ and annotation(for each read from FASTQ).

    return list of tuples, each element of the list is a test case, 
    each test case is a tuple(1st element sample name, 2nd fastq file, 3rd annotation file)

    raises exception if cannot parse file
    raises exception if cannot find fastq file
    raises exception if cannot find annotation file
    """
    num_description_for_sample = 3 # name fastq annotation
    samples = []
    reader = csv.reader(samples_file_handle, delimiter="\t")
    for row in reader:
        check_number_of_rows(row, num_description_for_sample)
        sample_name = row[0]
        fastq_file = row[1]
        check_file_exists("FASTQ", fastq_file)
        annotation_file = row[2]
        check_file_exists("Annotation", fastq_file)
        samples.append( (sample_name, fastq_file, annotation_file) )
    if not samples:
        raise IOError("empty sample file") 
    return samples
