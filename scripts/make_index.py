import argparse
import subprocess

def command_line_arguments():
    parser = argparse.ArgumentParser(
        description = "Makes index files for Bowtie2, Bwa, Blat and BMTagger.")
    parser.add_argument(
        "-g", required=True, type=str,
        help="Path to the genome file")
    parser.add_argument(
        "-o", required=True, type=str,
        help="Name of the organism")
    parser.add_argument(
        "-bowtie", required=False, type=int, default=0,
        help="Create index for bowtie 1 or 0. Default:0")
    parser.add_argument(
        "-bwa", required=False, type=int, default=0,
        help="Create index for bwa 1 or 0. Default:0")
    parser.add_argument(
        "-blat", required=False, type=int, default=0,
        help="Create index for blat 1 or 0. Default:0")
    parser.add_argument(
        "-bmtagger", required=False, type=int, default=0,
        help="Create index for bmtagger 1 or 0. Default:0")
    args = parser.parse_args()
    return args


def run_command(command, error_message):
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError:
        print error_message


def make_index_bowtie(genome, organism):
    command = ("bowtie2-build -f " + genome + " " + organism + ".fasta")
    run_command(command, "cannot run bowtie2. Check path to genome file.")


def make_index_bmtagger(genome, organism):
    command = ("bmtool -d " + genome + " -o " + organism + ".bitmask -A 0 -w 18")
    run_command(command, "cannot run bmtool. Check path to genome file.")


def make_index_bwa(genome, organism):
    command = ("bwa index " + genome)
    run_command(command, "cannot run bwa index. Check path to genome file.")


def make_index_blat(genome, organism):
    command = (
        "blat " + genome + " /dev/null /dev/null -tileSize=11 -makeOoc=" +
        organism + ".ooc repMatch=1024")
    run_command(command, "cannot run blat. Check path to genome file.")


if __name__=="__main__":
    args = command_line_arguments()
    if args.bowtie:
        make_index_bowtie(args.g, args.o)
    if args.bwa:
        make_index_bwa(args.g, args.o)
    if args.bmtagger:
        make_index_bmtagger(args.g, args.o)
    if args.blat:
        make_index_blat(args.g, args.o)
