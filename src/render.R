library(rmarkdown)
library(tools)

saveRDS(snakemake, file = ".snakemake_render.R.RDS")
report_file <- snakemake@params$report

sample_id <- snakemake@wildcards$sample

report_template <- normalizePath(snakemake@params$report, mustWork = TRUE)
output_path <- paste0(getwd(), "/", snakemake@output$report)

# Print log info (optional for debugging)
cat("Rendering Rmd file:\n")
cat("  Current directory: ", getwd(), "\n")
cat("  Rmd template     : ", report_template, "\n")
cat("  Output path      : ", output_path, "\n")

output_dir <- dirname(snakemake@output$report)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

render(input = report_template,
       output_file = output_path,
       envir = new.env(parent = globalenv()))