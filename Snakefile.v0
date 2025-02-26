# Snakefile for the 3P-Seq analysis pipeline
# This workflow processes all samples (preproc, align) and then branches into:
# - Combined analysis (all samples)
# - Separate asexual-only analysis
# - Separate sexual-only analysis
# It further performs downstream analyses: hexamer analysis, Fischer test and 
# integrates genome annotation via last exon identification and result combination.


import os

configfile: "config.yaml"

# Sample Group Definitions
SAMPLES = list(config["samples"].keys())
ASEXUAL_SAMPLES = [s for s, d in config["samples"].items() if d["strain"] == "Asexual"]
SEXUAL_SAMPLES = [s for s, d in config["samples"].items() if d["strain"] == "Sexual"]

# Global Parameters & Script Paths
READ_CUTOFF = config.get("read_cutoff", 10)
SEED = config.get("seed", 42)
SCRIPTS = config["scripts"]
ENVS = config["envs"]
GENOME = config["genome"]
OUTPUT = config["output"]

def get_output_dir(output_file):
    return os.path.dirname(output_file)

##############################
# Final Targets (rule all)
##############################
rule all:
    input:
        # Peak files
        expand(f"{OUTPUT['peaks']}{{strain}}_peaks.fa", strain=["combined", "asexual", "sexual"]),
        expand(f"{OUTPUT['peaks']}{{strain}}_peaks.txt", strain=["combined", "asexual", "sexual"]),
        
        # Hexamer analysis files
        expand(f"{OUTPUT['hexamer']}{{strain}}_hexamer_output.txt", strain=["combined", "asexual", "sexual"]),
        
        # Fisher test files
        expand(f"{OUTPUT['fisher']}{{strain}}_fisher_output.txt", strain=["combined", "asexual", "sexual"]),
        
        # Combined results
        f"{OUTPUT['results']}combined_results.txt"

##############################
# Main Processing (All Samples)
##############################
rule preproc:
    input:
        lambda wc: config["samples"][wc.sample]["data"]
    output:
        f"{OUTPUT['preprocessed']}{{sample}}.processed.fastq"
    params:
        genome=GENOME["fasta"],
        configfile="config.yaml",
        script=SCRIPTS["preproc"],
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['preprocessed']}logs/preproc_{{sample}}.log"
    conda:
        ENVS["preproc"]
    shell:
        """
        mkdir -p {params.output_dir}
        python {params.script} -q {input} -g {params.genome} -c {params.configfile} -o {output} > {log} 2>&1
        """

rule align:
    input:
        f"{OUTPUT['preprocessed']}{{sample}}.processed.fastq"
    output:
        f"{OUTPUT['aligned']}{{sample}}.bed"
    params:
        bowtie_index = GENOME["index"],
        configfile = "config.yaml",
        script = SCRIPTS["align"],
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['aligned']}logs/align_{{sample}}.log"
    conda:
        ENVS["align"]
    shell:
        """
        mkdir -p {params.output_dir}
        python {params.script} -q {input} -g {params.bowtie_index} -c {params.configfile} -o {output} > {log} 2>&1
        """

##############################
# Merging for Peak Detection
##############################
rule merge_bed:
    input:
        lambda wildcards: expand(f"{OUTPUT['aligned']}{{sample}}.bed",
                                 sample=SAMPLES if wildcards.strain == "combined" else
                                        ASEXUAL_SAMPLES if wildcards.strain == "asexual" else SEXUAL_SAMPLES)
    output:
        f"{OUTPUT['merged']}{{strain}}.bed"
    params:
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['merged']}logs/merge_bed_{{strain}}.log"
    conda:
        ENVS["default"]
    shell:
        """
        mkdir -p {params.output_dir}
        cat {input} > {output} 2> {log}
        """

##############################
# Peak Detection
##############################
rule peaks:
    input:
        f"{OUTPUT['merged']}{{strain}}.bed"
    output:
        f"{OUTPUT['peaks']}{{strain}}_peaks.txt"
    params:
        read_cutoff = READ_CUTOFF,
        script = SCRIPTS["peaks"],
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['peaks']}logs/peaks_{{strain}}.log"
    conda:
        ENVS["peaks"]
    shell:
        """
        mkdir -p {params.output_dir}
        python {params.script} {input} {params.read_cutoff} > {output} 2> {log}
        """

##############################
# FASTA Extraction using bedtools
##############################
rule peak_fasta:
    input:
        bed=f"{OUTPUT['peaks']}{{strain}}_peaks.txt",
        genome=GENOME["fasta"]
    output:
        f"{OUTPUT['peaks']}{{strain}}_peaks.fa"
    params:
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['peaks']}logs/peak_fasta_{{strain}}.log"
    conda:
        ENVS["bedtools"]
    shell:
        """
        mkdir -p {params.output_dir}
        bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output} 2> {log}
        """

##############################
# Extended Downstream Analysis (per branch)
##############################
rule hexamer:
    input:
        f"{OUTPUT['peaks']}{{strain}}_peaks.fa"
    output:
        f"{OUTPUT['hexamer']}{{strain}}_hexamer_output.txt"
    params:
        seed = SEED,
        script = SCRIPTS["hexamer"],
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['hexamer']}logs/hexamer_{{strain}}.log"
    conda:
        ENVS["hexamer"]
    shell:
        """
        mkdir -p {params.output_dir}
        python {params.script} {input} {params.seed} > {output} 2> {log}
        """

rule fisher:
    input:
        f"{OUTPUT['peaks']}{{strain}}_peaks.fa"
    output:
        f"{OUTPUT['fisher']}{{strain}}_fisher_output.txt"
    params:
        script = SCRIPTS["fisher"],
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['fisher']}logs/fisher_{{strain}}.log"
    conda:
        ENVS["fisher"]
    shell:
        """
        mkdir -p {params.output_dir}
        python {params.script} {input} > {output} 2> {log}
        """

##############################
# Extended Annotation Analysis (combined results)
##############################
rule identify_last_exons:
    input:
        genome_gtf = GENOME["gtf"]
    output:
        f"{OUTPUT['results']}last_exons.txt"
    params:
        script = SCRIPTS["identify_last_exons"],
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['results']}logs/identify_last_exons.log"
    conda:
        ENVS["default"]
    shell:
        """
        mkdir -p {params.output_dir}
        python {params.script} {input.genome_gtf} {output} > {log} 2>&1
        """

rule combine_results:
    input:
        last_exons = f"{OUTPUT['results']}last_exons.txt",
        fischer_result = f"{OUTPUT['fisher']}combined_fisher_output.txt"
    output:
        f"{OUTPUT['results']}combined_results.txt"
    params:
        script = SCRIPTS["combine_results"],
        output_dir=lambda wildcards, output: get_output_dir(output[0])
    log:
        f"{OUTPUT['results']}logs/combine_results.log"
    conda:
        ENVS["default"]
    shell:
        """
        mkdir -p {params.output_dir}
        python {params.script} {input.last_exons} {input.fischer_result} > {output} 2> {log}
        """
