#
# Setup metadata and folder files.
# 
R --script src/metadata.R
touch folder.csv
echo "ID" >> folder.csv

#
# Get and format metadata.
#
rule get_metadata:
    output:
        metadata_file = config["metadata_file"]
    params:
        metadata_url = config["metadata_url"]
    conda: "envs/r.yml"
    script: "metadata.R"

#
# Run snakemake.
#
snakemake -c1 --threads all 