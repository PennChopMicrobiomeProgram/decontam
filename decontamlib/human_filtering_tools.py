import csv
import os
import random
import re
import shutil
import subprocess
import tempfile
import utils

#compile regex to save time
hard = re.compile("\d+H")
soft = re.compile("\d+S")
all_cigar = re.compile("\D")
matches = re.compile("\d+M")


def run_command(command, error_message):
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print e.stderr


class _FilteringTool(object):
    def __init__(self, index):
        self.index = index

    def _get_mapped_reads(self, filename):
        """Extracts set of qnames from SAM file."""

        QNAME = tempfile.NamedTemporaryFile(delete=False)
        command = ("samtools view -F 4 " + filename + " > " + QNAME.name)
        run_command(command, "cannot run samtools view.")
        qname = utils.get_column(QNAME, 1)
        cigar = utils.get_column(QNAME, 6)
        mismatches = self._parse_mismatches_from_lists(
            utils.get_multiple_columns(QNAME, [12,13,14,15]))

        mapped = self._get_mapped_reads_from_cigar(qname, cigar, mismatches)
        os.remove(QNAME.name)
        return mapped

    def _get_mismatch(self, mismatch):
        return int(mismatch.split(":")[2])

    def _has_mismatch(self, x):
        return x.startswith("XM:i:")

    def _parse_mismatches_from_lists(self, cols):
        """Parse mismatch from variaous misc. column in sam file (XM:i:*)"""
        mismatches = []
        for row in cols:
            for x in row:
                if self._has_mismatch(x):
                    mismatches.append(self._get_mismatch(x))
                    break
            else:
                mismatches.append(0)
        return mismatches

    def _get_pattern_sum(self, matched_pattern):
        pattern_parsed = [s[0:len(s)-1] for s in matched_pattern]
        return sum(map(int, pattern_parsed))

    def _calculate_alignment_length(self, cigar_str):
        """Calculate alignment length from cigar string."""
        all = all_cigar.split(cigar_str)
        sum_all = sum(map(int, all[0:len(all)-1]))

        #Remove soft/hard clipped nucleotides from alignment length
        sum_soft_clip = 0
        sum_hard_clip = 0

        if  soft.match(cigar_str):
            sum_soft_clip = self._get_pattern_sum(soft.findall(cigar_str))

        elif hard.match(cigar_str):
            sum_hard_clip = self._get_pattern_sum(hard.findall(cigar_str))

        return sum_all - sum_soft_clip - sum_hard_clip

    def _calculate_identities(self, cigar_str, mismatch):
        """Calculate number of identities from cigar str.

        Equal to matches from cigar string minus mismatches (XM field in
        sam file).
        """
        sum_match = self._get_pattern_sum(matches.findall(cigar_str))
        return sum_match - mismatch

    def _get_mapped_reads_from_cigar(self, qname, cigar, mismatches):
        alignment_length = []
        identities = []

        for cigar_str, mismatch in zip(cigar, mismatches):
            alignment_length.append(
                self._calculate_alignment_length(cigar_str))
            identities.append(self._calculate_identities(cigar_str, mismatch))

        mapped = self._filter_mapped_reads(
            qname, self._calculate_pct_identity(identities, alignment_length),
            alignment_length)
        return mapped

    def _calculate_pct_identity(self, identities, alignment_length):
        return [float(iden)/alen for iden, alen in zip(identities, alignment_length)]

    def _filter_mapped_reads(
            self, qname, pct_identity, alignment_length, pct_identity_threshold=0.5,
            alignment_length_threshold=100):
        mapped = set()
        for qn, pct, alen in zip(qname, pct_identity, alignment_length):
            if pct >= pct_identity_threshold and alen >= alignment_length_threshold :
                mapped.add(qn)
        return mapped               


class Snap(_FilteringTool):
    def get_human_annotation(self, R1, R2):
        output = self._run(R1, R2)
        mapped = self._get_mapped_reads(output)
        os.remove(output)
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]

    def _run(self, R1, R2):
        output = tempfile.NamedTemporaryFile(delete=False)
        command = ("snap paired " + self.index + " " + R1 + " " + R2 +
            "  -o -sam " +  output.name)
        run_command(command, "cannot run snap. Check path to index file.")
        return output.name


class Bmfilter(_FilteringTool):
    def __init__(self, bitmask):
        self.bitmask = bitmask

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
        output = self._run(R1, R2)
        read_classifications = self._parse_bmtagger_output(open(output))
        return list(read_classifications)

    def _run(self, R1, R2):
        """ run bmtagger and return filename with the output file."""
        output = tempfile.NamedTemporaryFile()
        command = ("bmfilter -1 " + R1 + " -2 " + R2 + " -q 1 -T" + 
                " -b " + self.bitmask + " -o " + output.name)
        run_command(command, "cannot run bmtagger. Check path to bitmask.")
        return output.name + ".tag"

    @staticmethod
    def _parse_bmtagger_output(output):
        """ annotate each read.
        Args:
            output bmtagger tag file (as filehandle)
        """
        reader = csv.reader(output, delimiter="\t")
        reader.next() # skip header
        for row in reader:
            if len(row) != 2:
                raise ValueError(
                    "Expected 2 fields in BmTagger output, saw %s" % len(row))
            (read_id, annotation) = row
            is_human = (annotation == "H")
            yield (read_id, is_human)


class Bmtagger(Bmfilter):
    def __init__(self, bitmask, srprism):
        self.bitmask = bitmask
        self.srprism = srprism

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
        output = self._run_bmtagger(R1, R2)
        mapped = self._read_classification(open(output))
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]

    def _run_bmtagger(self, R1, R2):
        output = tempfile.NamedTemporaryFile(delete=False)
        tmp = tempfile.mkdtemp()
        command = ("/media/THING1/kyle/1205_PLEASE/bmtagger.sh -b " + self.bitmask + " -x " +
                   self.srprism + " -1 " + R1 + " -2 " + R2 + " -q 1 -T " + tmp + " -o " +
                   output.name)
        run_command(command, "cannot run bmtagger. Check path to bitmask.")
        shutil.rmtree(tmp)
        return output.name

    def _read_classification(self, out):
        mapped = set()
        for line in out:
            mapped.add(line.rstrip())
        return mapped
        
     
class Blat(object):
    def __init__(self, parameters):
        self.index = parameters["index"]

    def get_human_annotation(self, R1, R2):
        mapped = self._extract_blat_hits(R1)
        mapped.update(self._extract_blat_hits(R2))
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]

    def _extract_blat_hits(self, fastq_file):
        fasta = self.fastq_to_fasta(fastq_file)
        blat_psl = self.run_blat(fasta)
        blat_psl.seek(0)
        mapped = utils.extract_column(blat_psl, 10, 5)
        os.remove(blat_psl.name)
        os.remove(fasta)
        return mapped

    def _fastq_to_fasta(self, filename):
        fasta = tempfile.NamedTemporaryFile(delete=False)
        command = ("seqtk seq -a  " + filename + " > " + fasta.name)
        run_command(command, "cannot run seqtk for blat.")
        return fasta.name

    def _run_blat(self, R):
        output = tempfile.NamedTemporaryFile(delete=False)
        command = ("blat -minScore=50 -fastMap "+ self.index + " " + R +
                    " " +  output.name)
        run_command(command, "cannot run blat. Check path to index file.")
        return output


class Bwa(_FilteringTool):
    def __init__(self, parameters):
        self.index = parameters["index"]

    def get_human_annotation(self, R1, R2):
        output = self._run_bwa(R1, R2)
        mapped = self._get_mapped_reads(output)
        os.remove(output)
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]

    def _run_bwa(self, R1, R2):
        output = tempfile.NamedTemporaryFile(delete=False)
        command = ("bwa mem -M " + self.index + " " + R1 + " " + R2 +
            " > " +  output.name)
        run_command(command, "cannot run bwa. Check path to index file.")
        return output.name


class Bowtie(_FilteringTool):
    def __init__(self, index):
        self.index = index

    def get_human_annotation(self, R1, R2):
        output = self._run_bowtie(R1, R2)
        mapped = self._get_mapped_reads(output)
        os.remove(output)
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if id in mapped else 0) for id in ids]

    def _run_bowtie(self, R1, R2):
        output = tempfile.NamedTemporaryFile(delete=False)
        command = ("bowtie2 --local --very-sensitive-local -1 " + R1 + " -2 " + R2 + 
                " -x " + self.index + " -S " + output.name)
        run_command(command, "cannot run bowtie2. Check path to index file.")
        return output.name


class Random_human:
    def __init__(self, percent_human):
        assert 0.0 <= percent_human <= 100.0
        self.fraction_human = percent_human / 100.0

    def get_human_annotation(self, R1, R2):
        ids = utils.parse_read_ids(R1)
        return [(id, 1 if random.random() <= self.fraction_human else 0) for id in ids]


class None_human:
    def __init__(self):
        pass

    def get_human_annotation(self, R1, R2):
        ids = utils.parse_read_ids(R1)
        return [(id, 0) for id in ids]


class All_human:
    def __init__(self):
        pass

    def get_human_annotation(self, R1, R2):
        ids = utils.parse_read_ids(R1)
        return [(id, 1) for id in ids]
