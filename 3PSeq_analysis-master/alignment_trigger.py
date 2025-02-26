#!/usr/bin/python

import sys
import string
import re
import os
import subprocess
import argparse

__author__ = 'Deepak Poduval' ##modified from 'Praveen Anand'
__date__ = '25/02/25' # modified from  '27/10/2015'

def get_args():
    '''This function parses and return arguments passed in'''
    # Wrapper script description
    parser = argparse.ArgumentParser(
        description='Script is a wrapper for 3P-Seq pipeline')
    # Add arguments
    parser.add_argument(
        '-q', '--fastq', type=str, help='Processed FASTQ file of the reads [output from 3Pseq_iniprocess.py]', required=False)
    parser.add_argument(
        '-g', '--genome', type=str, help='Path to the bowtie genome index file', required=True)
    parser.add_argument(
        '-o', '--output', type=str, help='Output directory to store the results', required=False)
    # Array for all arguments passed to script
    args = parser.parse_args()
    # Assign args to variables
    fastq = args.fastq
    genome = args.genome
    output = args.output
    # Return all variable values
    return fastq, genome, output

# Fixed paths and options
BOWTIE_PATH = "/usr/local/bin/bowtie"
BOWTIE_M = 4
BOWTIE_V = 2
BOWTIE_P = 20

SAMTOOLS_PATH = "/usr/local/bin/samtools"
BEDTOOLS_PATH = "/usr/local/bin/bedtools"

# Getting the arguments
processedfastqfile, genome, output = get_args()

# Check for tool installations
if os.path.exists(BOWTIE_PATH):
    print("Found bowtie installation!!")
else:
    print("Sorry bowtie not found on your system.")
    sys.exit(1)

if os.path.exists(SAMTOOLS_PATH):
    print("Found samtools installation!!")
else:
    print("Sorry samtools installation not found on your system.")
    sys.exit(1)

if os.path.exists(BEDTOOLS_PATH):
    print("Found bam2bed!!!")
else:
    print("Sorry could not find 'bamToBed' program... Please specify a proper path.")
    sys.exit(1)

def run_bowtie3P(processedfastqfile):
    bowtie_cmd = (f"{BOWTIE_PATH} -q -m {BOWTIE_M} -v {BOWTIE_V} -p {BOWTIE_P} --chunkmbs 200 --un "
                  f"{processedfastqfile.replace('.fastq', '_unmapped.fastq')} {genome} {processedfastqfile} --sam "
                  f"{processedfastqfile.replace('.fastq', '_aligned.sam')}")
    p = subprocess.Popen(bowtie_cmd, stdin=None, stdout=None, shell=True)
    p.wait()
    if p.returncode != 0:
        return "Sorry!!..Something went wrong with the bowtie run... Ensure that all variables are specified properly."

    processed_sam_file = open(processedfastqfile.replace('.fastq', '_filtered.sam'), 'w')

    with open(processedfastqfile.replace('.fastq', '_aligned.sam')) as infile:
        for line in infile:
            if line[0] == "@":
                processed_sam_file.write(line)
            if line.split('\t')[1] == "0":
                MDZtag = str(line.split('\t')[12]).rstrip('\n')
                if re.search(r"[A-Z]0$", MDZtag):
                    processed_sam_file.write(line)
            if line.split('\t')[1] == "16":
                MDZtag = str(line.split('\t')[12]).rstrip('\n')
                if re.search(r"MD:Z:0[A-Z]", MDZtag):
                    processed_sam_file.write(line)
    return "Done with bowtie alignment and processing of sam file"

def gettingbed_file(filtered_samfile):
    samtools_cmd = (f"{SAMTOOLS_PATH} view -bS {filtered_samfile} > {filtered_samfile.replace('.sam', '.bam')}")
    samtoolrun = subprocess.Popen(samtools_cmd, stdin=None, stdout=None, shell=True)
    samtoolrun.wait()

    bedtools_cmd = (f"{BEDTOOLS_PATH} bamtobed -i {filtered_samfile.replace('.sam', '.bam')} > "
                    f"{filtered_samfile.replace('.sam', '.bed')}")
    bedtools_run = subprocess.Popen(bedtools_cmd, stdin=None, stdout=None, shell=True)
    bedtools_run.wait()

    bedfile = open(filtered_samfile.replace('.sam', '.bed'))
    range_count = {}
    for line in bedfile:
        bedcolumns = line.split('\t')
        uniqrange = f"{bedcolumns[0]}_{bedcolumns[5].rstrip('\n')}_{bedcolumns[1]}_{bedcolumns[2]}"
        try:
            range_count[uniqrange] += 1
        except KeyError:
            range_count[uniqrange] = 1

    processed_bed_file = open(filtered_samfile.replace('.sam', '.bedcount'), 'w')

    for i in range_count:
        new_list = i.split("_")
        new_list.append(str(range_count[i]))
        processed_bed_file.write("\t".join(new_list) + "\n")

print(run_bowtie3P(processedfastqfile))
print(gettingbed_file(processedfastqfile.replace('.fastq', '_filtered.sam')))
