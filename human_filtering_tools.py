import csv
import tempfile
import subprocess

class Bmtagger:

    name = "bmtagger"

    def __init__(self, parameters):
        if "bitmask" not in parameters:
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
        """ run bmtagger and return filename with the output file.
        """
        output = tempfile.NamedTemporaryFile()
        command = ("bmfilter -1 " + R1 + " -2 " + R2 + " -q 1 -T" + 
                " -b " + self.bitmask + " -o " + output.name)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError:
            print "cannot run bmtagger. Check path to bitmask."
        output = output.name + ".tag"
        return output

    def encode_human(self, is_human):
        """ bmtagger encodes human read as H we want 1 and 0 otherwise.
        """
        is_human_encoded = 0
        if is_human == "H":
            is_human_encoded = 1
        return is_human_encoded
    
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


class All_human:

    name = "all_human"

    def __init__(self, params):
        pass

    def get_human_annotation(self, fastq):
        pass

class Random_human:

    name = "random_human"

    def __init__(self, params):
        pass

    def get_human_annotation(self, fastq):
        pass

class None_human:

    name = "none_human"

    def __init__(self, params):
        pass

    def get_human_annotation(self, fastq):
        pass



