#!/usr/bin/python

import sys
import re
import os
import subprocess
import argparse

__author__ = 'Deepak Poduval' ##modified from 'Praveen Anand'
__date__ = '25/02/25' # modified from  '27/10/2015'

def get_args():
    '''This function parses and return arguments passed in'''
    # Wrapper script description
    parser = argparse.ArgumentParser(description='Script is a wrapper for 3P-Seq pipeline')
    # Add arguments
    parser.add_argument('-q', '--fastq', type=str, help='Processed FASTQ file of the reads [output from 3Pseq_iniprocess.py]', required=True)
    parser.add_argument('-g', '--genome', type=str, help='Path to the bowtie genome index file', required=True)
    parser.add_argument('-c', '--config', type=str, help='Configuration file', required=False)
    parser.add_argument('-o', '--output', type=str, help='Output file path to store the results', required=True)
    # Array for all arguments passed to script
    args = parser.parse_args()
    return args.fastq, args.genome, args.config, args.output

# Getting the arguments
fastq, genome, config, output = get_args()

# Directories and file paths
bam_file = output.replace('.bed', '.bam')
sam_file = output.replace('.bed', '.sam')
unmapped_fastq = fastq.replace('.fastq', '_unmapped.fastq')

# Fixed paths and options
BOWTIE_PATH = "/usr/local/bin/bowtie"
BOWTIE_M = 4
BOWTIE_V = 2
BOWTIE_P = 20

SAMTOOLS_PATH = "/usr/local/bin/samtools"
BEDTOOLS_PATH = "/usr/local/bin/bedtools"

# Check for tool installations
if not os.path.exists(BOWTIE_PATH):
    sys.exit("Error: Bowtie not found on your system.")
if not os.path.exists(SAMTOOLS_PATH):
    sys.exit("Error: Samtools not found on your system.")
if not os.path.exists(BEDTOOLS_PATH):
    sys.exit("Error: Bedtools not found on your system.")

def run_command(command):
    '''Run system commands'''
    try:
        print(f"Running command: {command}")
        result = subprocess.run(command, shell=True, check=True, capture_output=True)
        print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode(), file=sys.stderr)
        sys.exit(1)

def run_bowtie(fastq, genome, sam_file, unmapped_fastq):
    '''Run Bowtie'''
    cmd = (f"{BOWTIE_PATH} -q -m {BOWTIE_M} -v {BOWTIE_V} -p {BOWTIE_P} --chunkmbs 200 --un {unmapped_fastq} "
           f"{genome} {fastq} --sam {sam_file}")
    run_command(cmd)

def filter_sam(sam_file, filtered_sam_file):
    '''Filter SAM file'''
    with open(sam_file) as infile, open(filtered_sam_file, 'w') as outfile:
        for line in infile:
            if line.startswith("@") or re.search(r"MD:Z:\d+[A-Z]$", line):
                outfile.write(line)

def sam_to_bam(sam_file, bam_file):
    '''Convert SAM to BAM using Samtools'''
    cmd = f"{SAMTOOLS_PATH} view -bS {sam_file} > {bam_file}"
    run_command(cmd)

def bam_to_bed(bam_file, bed_file):
    '''Convert BAM to BED using Bedtools'''
    cmd = f"{BEDTOOLS_PATH} bamtobed -i {bam_file} > {bed_file}"
    run_command(cmd)

def count_bed_entries(bed_file, output_file):
    '''Count BED entries and write to output file'''
    range_count = {}
    with open(bed_file) as infile:
        for line in infile:
            bedcolumns = line.strip().split('\t')
            uniqrange = f"{bedcolumns[0]}_{bedcolumns[5]}_{bedcolumns[1]}_{bedcolumns[2]}"
            range_count[uniqrange] = range_count.get(uniqrange, 0) + 1

    with open(output_file, 'w') as outfile:
        for range_key, count in range_count.items():
            columns = range_key.split("_") + [str(count)]
            outfile.write("\t".join(columns) + "\n")

# Run alignment pipeline
run_bowtie(fastq, genome, sam_file, unmapped_fastq)
filtered_sam_file = sam_file.replace('.sam', '_filtered.sam')
filter_sam(sam_file, filtered_sam_file)
sam_to_bam(filtered_sam_file, bam_file)
bam_to_bed(bam_file, output)
count_bed_entries(output, output.replace('.bed', '.bedcount'))
print(f"Alignment pipeline completed. Output written to {output}")
