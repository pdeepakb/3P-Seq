#!/usr/bin/env python3

import sys
import re
import os
import subprocess
import argparse
import gzip

__author__ = 'Deepak Poduval' #updated from original 'Praveen Anand'
__date__ = '25/02/2025/' #updated from '27/10/2015'

def get_args():
    """This function parses and returns arguments passed in."""
    parser = argparse.ArgumentParser(
        description='Script is a wrapper for 3P-Seq pipeline')
    parser.add_argument(
        '-q', '--fastq', type=str, help='RAW FASTQ file of the reads', required=True)
    parser.add_argument(
        '-g', '--genome', type=str, help='Path to the complete genome file', required=False)
    parser.add_argument(
        '-c', '--config', type=str, help='Configuration file', required=False)
    parser.add_argument(
        '-o', '--output', type=str, help='Output directory to store the results', required=False)
    args = parser.parse_args()
    return args.fastq, args.genome, args.config, args.output

fastqfile, genome, config, output = get_args()

def open_fastq(filename):
    """Open a FASTQ file in text mode, using gzip if necessary."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def process_3pseq(fastqfile):
    """
    Processes the raw FASTQ file for 3P-Seq.
    Steps:
      1. Reverse complement the read sequence.
      2. Trim the terminal poly-A's to retain only 2 A's.
      3. Exclude sequences that contain 'N'.
      4. Exclude sequences shorter than 20 nucleotides.
    """
    # Determine processed filename; if gzipped, remove '.gz' appropriately.
    if fastqfile.endswith('.gz'):
        processed_filename = fastqfile.replace('.fastq.gz', '.processed.fastq')
    else:
        processed_filename = fastqfile.replace('.fastq', '.processed.fastq')

    with open(processed_filename, 'w') as processed_fastq_file, open_fastq(fastqfile) as fastqraw:
        current_name = None
        compdna = str.maketrans("ATGC", "TACG")
        trimmed_seq = ""
        second_name = ""
        for i, line in enumerate(fastqraw):
            line = line.rstrip('\n')
            if i % 4 == 0:
                current_name = line
            elif i % 4 == 1:
                sequence = line
                # Generate reverse complement and reverse the string
                revcomp = sequence.translate(compdna)[::-1]
                # Trim terminal A's to leave only two A's
                trimmed_seq = re.sub(r"A{2,}$", "AA", revcomp)
            elif i % 4 == 2:
                second_name = line
            elif i % 4 == 3:
                # Reverse the quality string so that its length matches trimmed_seq
                quality = line[::-1][:len(trimmed_seq)]
                if (trimmed_seq[-2:] == 'AA' and len(trimmed_seq) >= 20 and 'N' not in trimmed_seq):
                    processed_fastq_file.write(current_name + '\n')
                    processed_fastq_file.write(trimmed_seq + '\n')
                    processed_fastq_file.write(second_name + '\n')
                    processed_fastq_file.write(quality + '\n')

if __name__ == '__main__':
    process_3pseq(fastqfile)
