from Bio import SeqIO
from pprint import pprint


def extract_column(filehandle, column_number, skip_lines=0, sep="\t"):
    """ gets unique values for a given column (column number starting 1) """
    ids = set()
    skip = 0
    for line in filehandle:
        if skip < skip_lines:
            skip +=1
            continue       
        line = line.strip()
        columns = line.split(sep)
        ids.add(columns[column_number - 1])
    return ids

def get_column(filehandle, column_number, sep="\t"):
    id = list()
    
    for line in filehandle:
        line = line.strip()
        columns = line.split(sep)              
        id.append(columns[column_number - 1])
    filehandle.seek(0)
    return id

def get_multiple_columns(filehandle, column_numbers, sep="\t"):
    data = []
    for line in filehandle:
        line = line.strip()
        columns = line.split(sep)
        data.append([ columns[number - 1] for number in column_numbers ])

    return data

def parse_read_ids(fastq_filename):
    """extract fastq ids.
    
    Returns:
        set of ids.
    """
    records = SeqIO.index(fastq_filename, "fastq")
    ids = set()
    for record in records:
        ids.add(record)
    return ids


def check_all_read_ids_are_consistent(id1, id2):
    """ check if ids are the same.
        
    Args:
        id1, id2 sets of ids
    Returns:
        True if sets are the same 
        False if setst are diffent
    """
    return set(id1) == set(id2)


def add_tool_sample(tool, sample, human_annotation):
    """ add tool and sample column for a given human annotation.
        
    Args:
        tool string with tool name
        sample stiring with sample name
        human_annotation list of tuples: (read_id, is_human). 

    Returns:
        list of tuples: each tuple is (tool, sample, read_id, is_human) value.
    """
    number_or_rows = len(human_annotation)
    tool_rows = [ tool ] * number_or_rows
    sample_rows = [ sample ] * number_or_rows
    (read_id, is_human) = zip(*human_annotation)
    return zip(tool_rows, sample_rows, read_id, is_human)
