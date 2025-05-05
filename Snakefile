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
# which are going to be deleted after the kallisto bus workflow.
#
rule transfer_reads:
    output:
        reads = "data/reads/{sample}_R1.fastq.gz"
    params:
        output_dir = "data/reads",
        scp_address = config["scp_address"],
        scp_reads_location = config["scp_reads_location"]
    shell:
        """
            mkdir -p {params.output_dir}
            cp {input.reads} {output.reads}
        """

#
# Run kallisto bus workflow.
#
rule run_kallisto:
    input:
        kallisto_index = "data/reference_genome/kallisto_index.idx",
        reads = "data/reads/{sample}_R1.fastq.gz"
    output:
        quant = "data/quant/{sample}"
    params:
        technology_type = config["technology_type"],
        bootstrap = config["bootstrap"]
    threads: 8
    conda: "envs/kallisto.yml"
    shell:
        """
            mkdir -p {output.quant}

            STATS=($(./src/read_stats.sh {input.reads}))

            kallisto quant \
                --index {input.kallisto_index} \
                --output-dir {output.bus} \
                -b {params.bootstrap} \
                --threads {threads} \
                --pseudo \
                --single \
                -l ${STATS[0]} \
                -s ${STATS[1]} \
                {input.reads}
        """

#
# Final rule to generate all of the files.
#
rule all:
    input:
        bus = expand("data/bus/{sample}_R1.bus", 
                     sample = samples["ID"])
    output:
        tmp = "done.txt"
    shell:
        """
            touch {output.tmp}
        """