#
# Setup metadata file.
#
metadata_file=$(yq ".metadata_file" config.yml)
if [ -f "$metadata_file" ]; then
    Rscript src/metadata.R
fi
mkdir -p $(yq ".cell_dirs" config.yml | dirname)
echo "ID" > $(yq ".cell_dirs" config.yml)

#
# Run snakemake.
#
# SNAKEMAKE_JOBS=${JOBS:=1}
# SNAKEMAKE_CORES=${CORES:=1}

# snakemake \
#     --use-conda \
#     --conda-frontend mamba \
#     --rerun-triggers mtime \
#     --cores $SNAKEMAKE_CORES \
#     --jobs $SNAKEMAKE_JOBS \
#     $1
