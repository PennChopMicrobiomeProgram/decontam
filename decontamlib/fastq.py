import itertools
import os

class FastqSplitter(object):
    suffixes_human = {
        True: "_human",
        False: "",
        }
    suffixes_phix = {
        True: "_phix",
        False: "",
        }

    def __init__(self, input_fp, output_dir):
        self.output_dir = output_dir
        self.input_fp = input_fp
        input_filename = os.path.basename(input_fp)
        self.input_root, self.input_ext = os.path.splitext(input_filename)
        self._open_files = {}

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()

    def partition(self, annotations, organism):
        ids_to_annotation = dict(annotations)
        with open(self.input_fp) as f:
            for desc, seq, qual in parse_fastq(f):
                read_id = desc.split()[0]
                annotation = ids_to_annotation[read_id]
                self._write((desc, seq, qual), annotation, organism)

    def _write(self, read, annotation, organism):
        desc, seq, qual = read
        if organism == "human":
            suffix = self.suffixes_human[annotation]
        elif organism == "phix":
            suffix = self.suffixes_phix[annotation]
        output_filename = self.input_root + suffix + self.input_ext
        fp = os.path.join(self.output_dir, output_filename)
        if fp not in self._open_files:
            self._open_files[fp] = open(fp, "w")
        f = self._open_files[fp]
        write_fastq(f, desc, seq, qual)

    def close(self):
        for f in self._open_files.values():
            f.close()


def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return zip(*args)


def parse_fastq(f):
    """ parse original fastq file and write new fastq file the filtered non-host reads.
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
