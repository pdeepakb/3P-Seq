#!/usr/bin/env python3

import sys
import gzip

def open_gff(filename):
    """
    Open a GFF file. If the file is gzipped (ends with .gz),
    it will be opened in text mode.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def process_gff(gff_file, output_handle):
    """
    Process the GFF file line by line.
    This example writes lines that are not comments to the output.
    Replace the processing logic as needed.
    """
    for line in gff_file:
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue  # skip comment lines
        output_handle.write(line + "\n")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("Usage: identify_last_exons.py <genome.gff.gz> <output.txt>")
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    with open_gff(input_file) as infile, open(output_file, 'w') as outfile:
        process_gff(infile, outfile)
