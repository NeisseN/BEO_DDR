# Bray-Curtis distance matrix from relative abundances

# Niklas Neiße
# neisse.n@protonmail.com
# 2025, Jan. 23


# Calculating Bray-Curtis distance matrix between samples
#  from a dataframe of relative abundances (species).  

# Note: 
# 1. The samples are defined as rownames. 
# 2. Each species is one column
# 3. All entries much be numeric and not NA!
# 4. Saves the data in the dir_data as: 'prefix'_bray_df.csv


# Function ---------------------------------------------------------------------

fun_bray_curtis <- function(data, prefix, dir_data = "./") {
  
  # Step 1: Compute Bray-Curtis dissimilarity matrix
  bray_dist <- vegdist(data, method = "bray")
  
  # Step 2: Convert the distance matrix to a full matrix
  bray_matrix_full <- as.matrix(bray_dist)
  
  # Step 3: Set the upper triangle and diagonal to NA
  bray_matrix_full[upper.tri(bray_matrix_full)] <- NA 
  diag(bray_matrix_full) <- NA
  
  # Step 4: Convert the distance matrix into long format
  bray_long <- as.data.frame(as.table(bray_matrix_full))
  
  # Step 5: Rename columns for clarity
  colnames(bray_long) <- c("id", "id2", paste0(prefix, "_abund_bray"))
  
  # Step 6: Filter out rows where the distance is NA
  bray_long <- bray_long %>%
    filter(!is.na(get(paste0(prefix, "_abund_bray"))))
  
  # Optional: view first few rows
  head(bray_long)
  
  # Step 7: Save the long format data to CSV
  write.csv(bray_long,
            paste0(dir_data, prefix, "_bray_df.csv"), 
            row.names = FALSE)
  
  return(bray_long)  # Return the long format dataframe
}

