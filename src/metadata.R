library(stringr)
library(data.table)
library(tidyverse)
library(dplyr)
library(yaml)

# Read the configuration file.
current_dir <- getwd()
config <- yaml::read_yaml(paste0(current_dir, "/config.yml"))

# Extract necessary URL parts for google sheets to access
# the SmartSeq2 tab in the spreadsheet.
url_parts <- str_split(config$metadata_url, "/")[[1]]
sheet_id <- url_parts[[6]]
gid <- str_extract(url_parts[7], "\\d+")
url_csv <- paste0("https://docs.google.com/spreadsheets/d/", sheet_id, "/export?format=csv&gid=", gid)

# Read metadata.
metadata <- read.csv(url_csv)
metadata <- metadata %>%
  as.data.table() %>%
  .[ , 1:5]
colnames(metadata) <- c("ID", "Strain", "Tracing", "Source", "Species")

# Save metadata to a CSV file.
dir.create(dirname(paste0(current_dir, "/", config$metadata_file)), 
           recursive = TRUE, 
           showWarnings = FALSE)
fwrite(metadata, paste0(current_dir, "/", config$metadata_file))

