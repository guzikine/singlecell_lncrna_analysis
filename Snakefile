from snakemake.utils import min_version
min_version("7")
import os
import pandas as pd

configfile: "config.yml"

samples = pd.read_table(config["metadata_file"], sep=",")

reference_url = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz"
    if config["species"] == "human"
    else "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz"
)

wildcard_constraints:
    sample="[^/]+"



#
# Download the mouse M36 (GRCm39) full transcriptome reference genome.
#
rule download_reference_genome:
    output:
        reference_genome = "data/reference_genome/gencode.v37.transcripts.fa.gz"
    params:
        reference_genome_dir = "data/reference_genome",
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

checkpoint transfer_reads:
    output:
        reads = directory("data/raw_reads/{sample}/")
    params:
        output_dir = "data/raw_reads",
        scp_fastq_location = lambda wildcards: scp_fastq_location(wildcards)
    shell:
        """
            mkdir -p {params.output_dir}
            cp -r /Users/karolis/Downloads/Adameyko/SS2_15_0066/ {output.reads}
        """
            # mkdir -p {params.output_dir}
            # scp {params.scp_fastq_location} {params.output_dir}
            # tar -xvf {params.output_dir}/{wildcards.sample}_fastq.tar -C {output.reads}
            # rm {params.output_dir}/{wildcards.sample}_fastq.tar



#
# Move the reads locally from the rawdata/ directory.
#
rule move_reads:
    input:
        reads = directory("data/raw_reads/{sample}/")
    output:
        cell_reads = temp("data/raw_reads/{sample}/{cell_id}_R1.fastq.gz")
    params:
        cell_reads_location = "data/raw_reads/{sample}/rnaseq/mmu/{sample}/rawdata/{cell_id}/{cell_id}_R1.fastq.gz"
    shell:
        """
            mv {params.cell_reads_location} {output.cell_reads}
        """



#
# Run fastq quality control.
#
rule fastqc:
    input:
        "data/{type}_reads/{sample}/{cell_id}_R1.fastq.gz"
    output:
        html = "data/{type}_fastqc/{sample}/{cell_id}_R1.html",
        zip = "data/{type}_fastqc/{sample}/{cell_id}_R1_fastqc.zip"
    params:
        extra = config["preprocessing"]["fastqc_parameters"]
    log:
        "logs/{type}_fastqc/{sample}/{cell_id}.log"
    threads: 5
    resources:
        mem_mb = 1024
    wrapper:
        "v6.1.0/bio/fastqc"



#
# Run fastq screen to check for contamination.
#
rule fastq_screen:
    input:
        "data/raw_reads/{sample}/{cell_id}_R1.fastq.gz"
    output:
        txt = "data/fastq_screen/{sample}/{cell_id}.fastq_screen.txt",
        png = "data/fastq_screen/{sample}/{cell_id}.fastq_screen.png"
    params:
        fastq_screen_config = "config/fastq_screen.conf",
        subset = 100000,
        aligner = 'bowtie2'
    threads: 5
    wrapper:
        "v6.1.0/bio/fastq_screen"



#
# Run fastp for trimming and filtering the reads.
#
rule fastp_se:
    input:
        sample = ["data/raw_reads/{sample}/{cell_id}_R1.fastq.gz"]
    output:
        trimmed = "data/fastp_reads/{sample}/{cell_id}_R1.fastq.gz",
        failed = "data/fastp_reads/{sample}/{cell_id}.failed.fastq.gz",
        html = "data/fastp_reads/{sample}/{cell_id}.html",
        json = "data/fastp_reads/{sample}/{cell_id}.json"
    log:
        "logs/fastp/{sample}/{cell_id}.log"
    params:
        adapters = config["preprocessing"]["fastp_adapters"],
        extra = config["preprocessing"]["fastp_parameters"]
    threads: 5
    wrapper:
        "v6.1.0/bio/fastp"



#
# Run kallisto quant workflow.
#
rule run_kallisto:
    input:
        kallisto_index = "data/reference_genome/kallisto_index.idx",
        fastp_reads = "data/fastp_reads/{sample}/{cell_id}_R1.fastq.gz"
    output:
        abundance_h5 = "data/quant/{sample}/{cell_id}/abundance.h5",
        abundance_tsv = "data/quant/{sample}/{cell_id}/abundance.tsv",
        pseudoalignments = "data/quant/{sample}/{cell_id}/pseudoalignments.bam",
        run_info = "data/quant/{sample}/{cell_id}/run_info.json",
        log = "data/quant/{sample}/{cell_id}/kallisto.log"
    params:
        technology_type = config["technology_type"],
        bootstrap = config["bootstrap"],
        quant = "data/quant/{sample}/{cell_id}"
    threads: 8
    conda: "envs/kallisto.yml"
    shell:
        """
            mkdir -p {params.quant}

            STATS=($(./src/read_stats.sh {input.fastp_reads}))

            kallisto quant \
                --index {input.kallisto_index} \
                --output-dir {params.quant} \
                -b {params.bootstrap} \
                --threads {threads} \
                --pseudo \
                --single \
                -l ${{STATS[0]}} \
                -s ${{STATS[1]}} \
                {input.fastp_reads} &> {output.log}
        """



#
# Run samtools to index the bam file and generate 
# pseudoalignment statistics.
#
rule samtools:
    input:
        pseudoalignments = "data/quant/{sample}/{cell_id}/pseudoalignments.bam"
    output:
        stats = "data/quant/{sample}/{cell_id}/pseudoalignments_stats.txt",
        flagstat = "data/quant/{sample}/{cell_id}/pseudoalignments_flagstat.txt"
    conda: "envs/kallisto.yml"
    threads: 5
    shell:
        """
            samtools stats {input.pseudoalignments} > {output.stats}
            samtools flagstats {input.pseudoalignments} > {output.flagstat}
        """



#
# Generate a multiqc report.
#
def multiqc_inputs(wildcards):
    checkpoint_output = checkpoints.transfer_reads.get(sample = wildcards.sample).output.reads
    cell_ids = os.listdir(os.path.join(checkpoint_output, f"rnaseq/mmu/{wildcards.sample}/rawdata/"))
    
    return expand([
        "data/raw_fastqc/{sample}/{cell_id}_R1_fastqc.zip",
        "data/fastp_fastqc/{sample}/{cell_id}_R1_fastqc.zip",
        #"data/fastq_screen/{sample}/{cell_id}.fastq_screen.txt",
        "data/quant/{sample}/{cell_id}/kallisto.log",
        "data/quant/{sample}/{cell_id}/pseudoalignments_stats.txt",
        "data/quant/{sample}/{cell_id}/pseudoalignments_flagstat.txt",
    ], sample = wildcards.sample, cell_id = cell_ids)

rule multiqc:
    input:
        multiqc_inputs,
        config = "config/multiqc.yml"
    output:
        "data/reports/{sample}/multiqc_report.html"
    params:
        extra = "--verbose"
    log: 
        "logs/multiqc/{sample}.log"
    wrapper: "v3.9.0/bio/multiqc"



#
# Generate a report from the kallisto quant workflow.
#
# rule generate_report:
#     input:
#         multiqc = "data/reports/{sample}/multiqc_report.html"
#     output:
#         report = "data/reports/{sample}.html"
#     conda: "envs/kallisto.yml"
#     script: "generate_report.R"

#
# Final rule to generate all of the files.
#
rule all:
    input:
        reports = expand("data/reports/{sample}/multiqc_report.html", 
                         sample = ["SS2_15_0066"]) # samples["ID"]