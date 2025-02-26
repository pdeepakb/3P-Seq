import pandas as pd

def combine_results(abs_file, p3seq_file, exon_file, output_file):
    """
    Combine results from ABS-Scan, 3PSeq_analysis, and last exon identification.
    """
    abs_data = pd.read_csv(abs_file, sep="\t")
    p3seq_data = pd.read_csv(p3seq_file, sep="\t")
    exon_data = pd.read_csv(exon_file, sep="\t", header=None, names=["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])
    
    # Extract transcript IDs from exon attributes
    exon_data['transcript_id'] = exon_data['attributes'].str.extract(r'transcript_id \"(.*?)\"')
    
    # Merge datasets
    combined = abs_data.merge(p3seq_data, on=["chr", "strand"], how="inner")
    combined = combined.merge(exon_data[["chr", "start", "end", "strand", "transcript_id"]], on=["chr", "strand"], how="inner")
    
    # Save combined results
    combined.to_csv(output_file, sep="\t", index=False)

# Run script
if __name__ == "__main__":
    import sys
    abs_file = sys.argv[1]
    p3seq_file = sys.argv[2]
    exon_file = sys.argv[3]
    output_file = sys.argv[4]
    combine_results(abs_file, p3seq_file, exon_file, output_file)
