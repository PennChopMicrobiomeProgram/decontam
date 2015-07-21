import json
import os
import glob
import csv
import argparse

def makeReport(illqc_summary_dir, decontam_summary_dir, output_fp):
    """
    Compiles the summary files for individual samples into one tsv file
    :param illqc_summary_dir: path to the directory where all the illqc json files are located
    :param decontam_summary_dir: path to the directory where all the decontam json files are located
    :param output_fp: filepath of the output
    """
    illqc_header = ['input', 'both kept', 'rev only', 'dropped', 'fwd only']
    decontam_header = ['true', 'false']

    with open(output_fp, 'wb') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(['Sample'] + [s.replace(" ", "_") for s in illqc_header] + ['human', 'non_human'])
        
        for file in glob.glob(os.path.join(illqc_summary_dir, "*.json")):
            file = os.path.basename(file)
            
            ill = getValues(os.path.join(illqc_summary_dir, file), illqc_header)
            de = getValues(os.path.join(decontam_summary_dir, file), decontam_header)
                
            ill.extend(de)
            ill.insert(0, file.rsplit('_', 1)[0])
            writer.writerow(ill)


def getValues(fp, headers):
    """
    Obtains the necessary values from the dictionary.
    :param fp: filepath of the summary file
    :param headers: the keys that will be extracted
    :return: list of values. None if value is not found
    """
    if os.path.isfile(fp):
        with open(fp) as f_in:
            summary = json.load(f_in)
            summary_data = summary.get('data', {})
    else:
        summary_data={}
    return [summary_data.get(header, None) for header in headers]

def main(argv=None):
    p=argparse.ArgumentParser()

    # input
    p.add_argument("--illqc-dir", required=True,
                   help="Directory for illqc summary files")
    p.add_argument("--decontam-dir", required=True,
                   help="Direcrory for decontamination summary files")

    # output
    p.add_argument("--output-fp", required=True,
                   help="Output report file")
    args=p.parse_args(argv)
    
    makeReport(args.illqc_dir, args.decontam_dir, args.output_fp)
