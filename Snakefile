from snakemake.utils import min_version
min_version("7")
import os
import pandas as pd

configfile: "config.yml"

samples = pd.read_table(config["metadata_file"], sep=",")
cells = pd.read_table(config["cell_dirs"], sep=",")

reference_url = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz"
    if config["species"] == "human"
    else "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz"
)

#
# Download the mouse M36 (GRCm39) full transcriptome reference genome.
#
rule download_reference_genome:
    output:
        reference_genome = "data/reference_genome/gencode.v37.transcripts.fa.gz"
    params:
        reference_genome_dir = "data/reference_genome"
        reference_url = reference_url
    shell:
        """
            mkdir -p {params.reference_genome_dir}
            wget -O {output.reference_genome} {params.reference_url}
        """

#
# Build the kallisto index for full human transcriptome using the downloaded reference genome.
#
rule build_kallisto_index:
    input:
        reference_genome = "data/reference_genome/gencode.v37.transcripts.fa.gz"
    output:
        kallisto_index = "data/reference_genome/kallisto_index.idx"
    params:
        reference_genome_dir = "data/reference_genome"
    conda: "envs/kallisto.yml"
    shell:
        """
            mkdir -p {params.reference_genome_dir}
            kallisto index -i {output.kallisto_index} {input.reference_genome}
        """

# 
# Transfer the reads one by one to minimize the load on the cluster,
# which are going to be deleted after the kallisto quant workflow.
#
# Determine what are the names of the directories inside: 
# data/reads/{sample}/rnaseq/mmu/{sample}/rawdata/
#
def scp_fastq_location(wildcards):
    return f"{config['scp_address']}:{config['scp_reads_location']}{wildcards.sample}/{wildcards.sample}_fastq.tar"

rule transfer_reads:
    output:
        reads = temp(dir("data/reads/{sample}/")),
        cell_dirs = config["cell_dirs"]
    params:
        output_dir = "data/reads",
        scp_fastq_location = lambda wildcards: scp_fastq_location(wildcards)
    shell:
        """
            mkdir -p {params.output_dir}
            scp {params.scp_fastq_location} {params.output_dir}
            tar -xvf {params.output_dir}/{wildcards.sample}_fastq.tar -C {output.reads}
            rm {params.output_dir}/{wildcards.sample}_fastq.tar
            echo "ID" > {output.cell_dirs}
            ls $(find {input.reads} -type d -name "rawdata" -print -quit) >> {output.cell_dirs}
        """

#
# Run fastq quality control.
#
rule fastqc:
    input:
        reads = dir("data/reads/{sample}")
    output:
        fastqc_report = "data/raw_fastqc/{sample}/{cell_dir}/{cell_dir}_R1_fastqc.html",
        fastqc_zip = "data/raw_fastqc/{sample}/{cell_dir}/{cell_dir}_R1_fastqc.zip"
    params:
        input_reads = "data/reads/{sample}/{cell_dir}/{cell_dir}_R1.fastq.gz",
        output_dir = "data/fastqc/{sample}/{cell_dir}"
    conda: "envs/kallisto.yml"
    threads: 5
    shell:
        """
            fastqc \
                --outdir {params.output_dir} \
                --threads {threads} \
                --quiet \
                {params.input_reads}
        """

#
# Run fastq screen to check for contamination.
#
# fastq_screen --outdir {params} \
#     --aligner {params} \
#     --conf {params} \
#     --subset {params} \
#     --force \
#     {output}

#
# Run fastp for trimming and filtering and then 
# fastqc on fastp processed reads.
#
rule fastp:
    input:
        reads = dir("data/reads/{sample}")
    output:
        fastp_report = "data/fastp_fastqc/{sample}/{cell_dir}/{cell_dir}_R1.fastq.gz",
        fastp_json = "data/fastp_fastqc/{sample}/{cell_dir}/{cell_dir}_R1.fastq.json",
        fastp_reads = "data/reads/{sample}/{cell_dir}/fastp_{cell_dir}_R1.fastq.gz",
        fastqc_zip = "data/fastp_fastqc/{sample}/{cell_dir}/{cell_dir}_R1.fastq.zip",
        fastqc_report = "data/fastp_fastqc/{sample}/{cell_dir}/{cell_dir}_R1_fastqc.html"
    params:
        input_reads = "data/reads/{sample}/{cell_dir}/{cell_dir}_R1.fastq.gz",
        output_dir = "data/fastp_fastqc/{sample}/{cell_dir}",
        fastp_parameters = config["fastp_parameters"]
    conda: "envs/kallisto.yml"
    threads: 5
    shell:
        """
            mkdir -p {params.output_dir}
            fastp \
                {params.fastp_parameters} \
                --thread {threads} \
                {params.input_reads} \
                {params.output_reads} \
                --html {output.fastp_report} \
                --json {output.fastp_json}

            fastqc \
                --outdir {params.output_dir} \
                --threads {threads} \
                --quiet \
                {input.fastp_reads}
        """

#
# Run kallisto quant workflow.
#
rule run_kallisto:
    input:
        kallisto_index = "data/reference_genome/kallisto_index.idx",
        fastp_reads = "data/reads/{sample}/{cell_dir}/fastp_{cell_dir}_R1.fastq.gz"
    output:
        abundance_h5 = "data/quant/{sample}/{cell_dir}/abundance.h5",
        abundance_tsv = "data/quant/{sample}/{cell_dir}/abundance.tsv",
        pseudoalignments = "data/quant/{sample}/{cell_dir}/pseudoalignments.bam",
        run_info = "data/quant/{sample}/{cell_dir}/run_info.json",
        quant = dir("data/quant/{sample}/{cell_dir}"),
        log = "data/logs/{sample}/{cell_dir}_kallisto.log"
    params:
        technology_type = config["technology_type"],
        bootstrap = config["bootstrap"]
    threads: 8
    conda: "envs/kallisto.yml"
    shell:
        """
            mkdir -p {output.quant}

            STATS=($(./src/read_stats.sh {input.fastp_reads}))

            kallisto quant \
                --index {input.kallisto_index} \
                --output-dir {output.quant} \
                -b {params.bootstrap} \
                --threads {threads} \
                --pseudo \
                --single \
                -l ${STATS[0]} \
                -s ${STATS[1]} \
                {input.fastp_reads} &> {output.log}
        """

#
# Run samtools to index the bam file and generate stats.
#
rule samtools:
    input:
        pseudoalignments = "data/quant/{sample}/{cell_dir}/pseudoalignments.bam"
    output:
        index = "data/quant/{sample}/{cell_dir}/pseudoalignments.bam.bai",
        stats = "data/quant/{sample}/{cell_dir}/pseudoalignments_stats.txt",
        flagstat = "data/quant/{sample}/{cell_dir}/pseudoalignments_flagstat.txt"
    conda: "envs/kallisto.yml"
    threads: 5
    shell:
        """
            samtools index {input.pseudoalignments} {output.index}
            samtools stats {input.pseudoalignments} > {output.stats}
            samtools flagstats {input.pseudoalignments} > {output.flagstat}
        """

#
# Run multiqc to generate a report from the fastqc, 
# fastp and samtools reports.
#
rule multiqc:
    input:
        expand("data/raw_fastqc/{{sample}}/{cell_dir}/{cell_dir}_R1_fastqc.zip",
               cell_dir = cells["ID"]),
        expand("data/fastp_fastqc/{{sample}}/trimmed_{cell_dir}/{cell_dir}_R1.fastq.zip",
               cell_dir = cells["ID"]),
        # expand("fastqc_screen"),
        expand("data/quant/{{sample}}/{cell_dir}/run_info.json",
               cell_dir = cells["ID"]),
        expand("data/quant/{{sample}}/{cell_dir}/pseudoalignments_stats.txt",
               cell_dir = cells["ID"]),
        expand("data/quant/{{sample}}/{cell_dir}/pseudoalignments_flagstat.txt",
               cell_dir = cells["ID"])
    output:
        "data/reports/{sample}/multiqc_report.html"
    params:
        extra = "-c ./multiqc/multiqc_config.yml"
    log: "data/logs/{sample}}/multiqc.log"
    wrapper: "v3.9.0/bio/multiqc"

#
# Generate a report from the kallisto quant workflow.
#
rule generate_report:
    input:
        quant_dirs = expand(dir("data/quant/{{sample}}/{cell_dir}"),
                            cell_dir = cells["ID"])
    output:
        report = "data/reports/{sample}.html"
    conda: "envs/kallisto.yml"
    script: "generate_report.R"

#
# Final rule to generate all of the files.
#
rule all:
    input:
        reports = expand("data/reports/{sample}.html", 
                         sample = samples["ID"])