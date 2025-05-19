saveRDS(snakemake, ".generate_count_matrix.R.RDS")

#
# Load libraries.
#
library(data.table)
library(stringr)
library(dplyr)
library(tidyverse)
library(foreach)
library(Seurat)

#
# Load TSV files.
#
tsv_files <- snakemake@input$abundance_tsv
sample_id <- snakemake@wildcards$sample

#
# Analysis.
#
cell_ids <- sapply(tsv_files, function(x) {
  splitted <- str_split(x, "/")[[1]]
  return(splitted[length(splitted) - 1])
})

matrices <- lapply(tsv_files, fread, sep = "\t")
names(matrices) <- cell_ids

matrices <- lapply(matrices, function(x) {
  dt <- as.data.table(x) %>%
    .[, c("transcript_id", "gene_id", "ott_gene_id", "ott_transcript_id", "transcript_name", "gene_name", "transcript_length", "biotype") := 
        tstrsplit(target_id, "\\|", 
                  fixed = FALSE)]
})

common_elements <- Reduce(intersect, 
                          lapply(matrices, function(dt) dt[["transcript_name"]]))

matrices_filtered <- lapply(matrices, function(dt) {
  dt_filtered <- dt[transcript_name %in% common_elements]
  dt_filtered <- dt_filtered[match(common_elements, transcript_name)]
  return(dt_filtered[["est_counts"]])
})

count_matrix <- do.call(cbind, matrices_filtered)

rownames(count_matrix) <- common_elements
colnames(count_matrix) <- cell_ids

annotations <- matrices[[1]][ transcript_name %in% common_elements, c("transcript_id", "gene_id", "ott_gene_id", "ott_transcript_id", "transcript_name", "gene_name", "transcript_length", "biotype")]
annotations <- annotations[match(common_elements, transcript_name)]

#
# Transfer everything to Seurat.
#
seurat_obj <- CreateSeuratObject(
  counts = count_matrix,
  project = sample_id,
  min.cells = 3,
  min.features = 200
)

#
# Save data object.
#
data <- list(
  seurat_obj = seurat_obj,
  annotations = annotations
)

saveRDS(data, snakemake@output$count_matrix)