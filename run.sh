#
# Setup metadata files.
#
# Setup file path variables.
metadata_file=$(yq ".metadata_file" config.yml)

# If files already exist, skip the creation.
if [ -f "$metadata_file" ]; then
    mkdir -p $(dirname $metadata_file)
    Rscript src/metadata.R
fi

#
# Run snakemake.
#
SNAKEMAKE_JOBS=${JOBS:=1}
SNAKEMAKE_CORES=${CORES:=1}

snakemake \
    --use-conda \
    --conda-frontend mamba \
    --rerun-triggers mtime \
    --cores $SNAKEMAKE_CORES \
    --jobs $SNAKEMAKE_JOBS \
    $1
