# Reconstruct phyloseq object from BExIS CSVs
# Author: Niklas Neisse
# Date: 2025-08-13

library(tidyverse)
library(janitor)
library(phyloseq)

# Directories
dir_wd   <- file.path('C:', 'code', 'BEO_DDR')
dir_data <- file.path(dir_wd, 'data')
dir_script <- file.path(dir_wd, 'scripts')

# File paths
file_rarefied <- file.path(dir_data, 'rarefied.RDS')
file_meta     <- file.path(dir_data, 'df_bexis_meta.csv')
file_asv      <- file.path(dir_data, 'df_bexis_asv.csv')
file_distances <- file.path(dir_data, 'distances_df.csv')

if (!file.exists(file_rarefied)) {
  
  message('rarefied.RDS not found — reconstructing from CSVs...')
  
  # Load and clean
  df_meta <- read.csv(file_meta) %>% clean_names()
  df_asv  <- read.csv(file_asv) %>% clean_names()
  
  # Detect taxonomy columns automatically
  tax_cols <- df_asv %>%
    select(-id) %>%
    select(where(~ !is.numeric(.))) %>%
    names()
  
  # Taxonomy table
  df_taxa <- df_asv %>%
    select(id, all_of(tax_cols)) %>%
    column_to_rownames('id')
  
  # OTU table (samples = columns, taxa = rows)
  df_otu <- df_asv %>%
    select(-all_of(tax_cols)) %>%
    column_to_rownames('id')
  
  # Convert all sample names to uppercase
  colnames(df_otu) <- toupper(colnames(df_otu))
  if ('id' %in% names(df_meta)) {
    df_meta$id <- toupper(df_meta$id)
  }
  
  # If meta has a sample_id/id column, use it as row names
  if ('sample_id' %in% names(df_meta)) {
    rownames(df_meta) <- df_meta$sample_id
    df_meta <- df_meta %>% select(-sample_id)
  } else if ('id' %in% names(df_meta)) {
    rownames(df_meta) <- df_meta$id
    df_meta <- df_meta %>% select(-id)
  }
  
  # Ensure ordering matches OTU table
  df_meta <- df_meta[colnames(df_otu), , drop = FALSE]
  
  # Build phyloseq object
  ps <- phyloseq(
    otu_table(as.matrix(t(df_otu)), taxa_are_rows = FALSE),
    tax_table(as.matrix(df_taxa)),
    sample_data(df_meta)
  )
  
  saveRDS(ps, file_rarefied)
  message('rarefied.RDS created successfully.')
  
} else {
  message('rarefied.RDS already exists — skipping reconstruction.')
}

# --- Extra step: Run data_distance.R if distances_df.csv doesn't exist ---
if (!file.exists(file_distances)) {
  message('distances_df.csv not found — running data_distance.R...')
  source(file.path(dir_script, 'data_distance.R'))
} else {
  message('distances_df.csv already exists — skipping distance calculation.')
}
