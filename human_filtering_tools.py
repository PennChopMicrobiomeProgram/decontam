import csv
import tempfile
import subprocess
import random
import os

import utils

def run_command(command, error_message):
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError:
        print error_message

class Bmtagger:

    name = "bmtagger"

    def __init__(self, parameters):
        if parameters is None or "bitmask" not in parameters:
            raise KeyError("parameter dictionary for Bmtagger should have key 'bitmask'.")
        self.bitmask = parameters["bitmask"] 

    def get_human_annotation(self, R1, R2):
        """ creates human read annotation by running a tool.
            Args:
                R1, R2 forward and reverese reads
            Returns:
                list of tuples, each tuple is a pair of: read_id and is_human
                is_human = 1 if read annotated as human read
                is_human = 0 if read annotated as NON human read
                (so in R sum(is_human) is number of reads annotated as human reads`)
                length of list equal to number of read in the fastq files.
        """
        output = self.run_bmtagger(R1, R2)
        read_classification = self.parse_bmtagger_output(open(output))
        return read_classification

    def run_bmtagger(self, R1, R2):
        """ run bmtagger and return filename with the output file."""
        output = tempfile.NamedTemporaryFile()
        command = ("bmfilter -1 " + R1 + " -2 " + R2 + " -q 1 -T" + 
                " -b " + self.bitmask + " -o " + output.name)
        run_command(command, "cannot run bmtagger. Check path to bitmask.")
        return output.name + ".tag"

    def encode_human(self, is_human):
        """ bmtagger encodes human read as H we want 1 and 0 otherwise."""
        return (1 if is_human == "H" else 0)
    
    def parse_bmtagger_output(self, output):
        """ annotate each read.
        Args:
            output bmtagger tag file (as filehandle)
        """
        human_annotation = []
        reader = csv.reader(output, delimiter="\t")
        reader.next() # skip header
        for row in reader:
            if len(row) != 2:
                raise IOError("cannot process bmtagger output.")
            (read_id, is_human) = row
            is_human_encoded = self.encode_human(is_human)
            human_annotation.append( (read_id, is_human_encoded) )
        return human_annotation


class Blat:

    name = "blat"
    def __init__(self, parameters):
        if parameters is None or "index" not in parameters:
            raise KeyError("parameter dictionary for blat should have key 'index'.")
        self.index = parameters["index"]

    def get_human_annotation(self, R1, R2):
        if not ".fasta" in R1:
            #Blat only accepts input in fasta format.
            R1_fasta = self.fastq_to_fasta(R1)
            R2_fasta = self.fastq_to_fasta(R2)
        
        output_R1 = self.run_blat(R1_fasta)
        output_R2 = self.run_blat(R2_fasta)
        
        mapped = self.get_mapped_reads(output_R1)
        mapped.update(self.get_mapped_reads(output_R2))
        os.remove(output_R1)
        os.remove(output_R2)
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]

    def fastq_to_fasta(self, filename):
        fasta = tempfile.NamedTemporaryFile(delete=False)
        command = ("seqtk seq -a  " + filename + " > " + fasta.name)
        run_command(command, "cannot run seqtk for blat.")
        return fasta.name

        
    def run_blat(self, R):
        output = tempfile.NamedTemporaryFile(delete=False)
        command = ("blat -minScore=50 -fastMap "+ self.index + " " + R +
                    " " +  output.name)
        run_command(command, "cannot run blat. Check path to index file.")
        return output.name

    def get_mapped_reads(self, filename):
        QNAME = tempfile.NamedTemporaryFile(delete=False)
        command = ("cut -c20-64  " + filename + " > " + QNAME.name)
        run_command(command, "cannot run cut for blat.")
        mapped = self.parse_ids(QNAME)
        os.remove(QNAME.name)
        return mapped

    def parse_ids(self, filehandle):
        ids = set()
        for line in filehandle:
            ids.add(line.strip())
        return ids
                                                                                                                                                                                        

class Bwa:
    
    name = "bwa"
    
    def __init__(self, parameters):
        if parameters is None or "index" not in parameters:
            raise KeyError("parameter dictionary for Bwa should have key 'index'.")
        self.index = parameters["index"]
        
    def get_human_annotation(self, R1, R2):
        output = self.run_bwa(R1, R2)
        mapped = self.get_mapped_reads(output)
        os.remove(output)
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]
                                      
    def run_bwa(self, R1, R2):
        output = tempfile.NamedTemporaryFile(delete=False)
        command = ("bwa mem -M " + self.index + " " + R1 + " " + R2 +
            " > " +  output.name)
        run_command(command, "cannot run bwa. Check path to index file.")
        return output.name

    def get_mapped_reads(self, filename):
        QNAME = tempfile.NamedTemporaryFile(delete=False)
        command = ("samtools view -F 4 " + filename + " | cut -f 1 > " + QNAME.name)
        run_command(command, "cannot run samtools view.")
        mapped = self.parse_ids(QNAME)
        os.remove(QNAME.name)
        return mapped

    def parse_ids(self, filehandle):
        ids = set()
        for line in filehandle:
            ids.add(line.strip())
        return ids
                                                                                    

class All_human:

    name = "all_human"

    def __init__(self, params):
        pass

    def get_human_annotation(self, R1, R2):
        ids = utils.parse_read_ids(R1)
        return [(id, 1) for id in ids]

class Bowtie:

    name = "bowtie"

    def __init__(self, parameters):
        if parameters is None or "index" not in parameters:
            raise KeyError("Bowtie requires parameter key/value 'index' : 'path to reference'.")
        self.index = parameters["index"]

    def get_human_annotation(self, R1, R2):
        output = self.run_bowtie(R1, R2)
        mapped = self.get_mapped_reads(output)
        os.remove(output)
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]

    def run_bowtie(self, R1, R2):
        output = tempfile.NamedTemporaryFile(delete=False)
        command = ("bowtie2 -1 " + R1 + " -2 " + R2 + 
                " -x " + self.index + " -S " + output.name)
        run_command(command, "cannot run bowtie2. Check path to index file.")
        return output.name

    def get_mapped_reads(self, filename):
        QNAME = tempfile.NamedTemporaryFile(delete=False)
        command = ("samtools view -F 4 " + filename + " | cut -f 1 > " + QNAME.name)
        run_command(command, "cannot run samtools view.")
        mapped = self.parse_ids(QNAME)
        os.remove(QNAME.name)
        return mapped

    def parse_ids(self, filehandle):
        ids = set()
        for line in filehandle:
            ids.add(line.strip())
        return ids

class Random_human:

    name = "random_human"

    def __init__(self, params):
        assert 0.0 <= params["percent_human"] <= 100.0
        self.fraction_human = params["percent_human"]/100.0

    def get_human_annotation(self, R1, R2):
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if random.random() <= self.fraction_human else 0) for id in ids]

class None_human:

    name = "none_human"

    def __init__(self, params):
        pass

    def get_human_annotation(self, R1, R2):
        ids = utils.parse_read_ids(R1)
        return [(id, 0) for id in ids]



