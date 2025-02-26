# Snakemake pipeline for 3P-Seq Analysis

This repository contains a Snakemake pipeline for processing and analyzing sequencing data. The pipeline uses several tools, including ABS-Scan and 3PSeq Analysis, along with custom scripts for post-processing. It also takes care of sexual and asexual seperately.

## Overview
This pipeline processes 3P-Seq data to analyze polyadenylation events in *Schmidtea mediterranea*. The workflow includes:
- Quality control with `fastqc`
- Read alignment using `bowtie`
- Peak detection for 3' polyadenylation sites
- Hexamer sequence analysis and secondary PAS signal detection

## Installation and Running the Pipeline
Run the following command to install dependencies:
```bash
brew install fastqc bowtie bedtools samtools wget
pip3 install joblib scipy numpy
```
### 1. Environment Setup

Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual) installed. Then, create and activate a conda environment for running Snakemake:

```bash
# Create the environment (if not already created)
conda create -n snakemake -c bioconda -c conda-forge snakemake

# Activate the environment
conda activate snakemake
```

*(Note: If your environment is already named "snakemake", simply activate it with `conda activate snakemake`.)*

### 2. Verifying the Installation

Confirm that Snakemake is installed by checking its version:

```bash
snakemake --version
```

A version number should be printed.

### 3. Configuring the Pipeline

Review the `config.yaml` file to ensure the sample data, genome paths, and output directories match your setup. Update any paths if needed.

### 4. Testing the Snakefile

Before running the full pipeline, perform a dry-run to detect any issues:

```bash
snakemake -n
```

This command will print the planned jobs without executing them. You can also check the workflow syntax with:

```bash
snakemake --lint
```

This command will generate graphical representation of the piepline
```bash
snakemake --dag | dot -Tsvg > dag.svg
```

### 5. Running the Pipeline

After verifying that the dry-run completes without errors, start the pipeline. Specify the number of cores (adjust based on your system’s capability):

```bash
snakemake --cores 4
```

If your workflow rules use individual conda environments, run:

```bash
snakemake --use-conda --cores 4
```

### 6. Monitoring and Output

- **Logs and Output Directories:**  
  The pipeline writes output to directories specified in the `config.yaml` file (e.g., `qc/`, `data/trimmed/`, `data/aligned/`, `results/`). Check these directories to monitor the progress and results of your analysis.

- **Interrupting and Resuming:**  
  You can interrupt the pipeline with `Ctrl+C`. Resume the analysis by rerunning the same `snakemake` command; it will pick up from where it left off.

## Additional Notes

- Make sure all input files (FASTQ, genome FASTA, GFF) are in the correct locations as specified in `config.yaml`.
- Custom scripts (`identify_last_exons.py`, `combine_results.py`) and tools (`ABS-Scan-master`, `3PSeq_analysis-master`) should have the proper permissions to execute.

# making bowtie index

``bash
bowtie-build --wrapper basic-0 SmedSxl_genome_v4.0.nt.gz bowtie_index/genome
``


By following these steps, you can install, test, and run the Snakemake-based workflow successfully.
```
