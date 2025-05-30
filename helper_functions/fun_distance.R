# Distance matrix from any set of varibales

# Niklas Neiße
# neisse.n@protonmail.com
# 2025, Jan. 23


# Calculating Distance matrix between samples
#  from a dataframe of any variables.  

# Note: 
# 1. The samples are defined as rownames. 
# 2. Each variable is one column
# 3. All entries much be numeric
#     (NA works for all `dist` fun methods)
# 4. Saves the data in the dir_data as: prefix_method_df.csv

# 5. Method: `vegdist` fun methods:
#     "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", 
#     "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", 
#     "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", 
#     "aitchison", or "robust.aitchison"
# 5.1 as well as all distm {geosphere} methods:
#     "distCosine", "distGeo", "distHaversine"

# Function ---------------------------------------------------------------------

fun_distance <- function(data, prefix, method = "manhattan", scale_data = FALSE, dir_data = "./") {
  
  # Extract the first 4 characters of the method for the filename "meth_fix"
  meth_fix <- substr(method, 1, 4)
  
  # Step 1: Optionally scale the data
  # Standardize the data if scale_data = TRUE
  if (scale_data) {
    data <- scale(data)  
  }
  
  # Step 2: Check if method starts with 'dist' (for geosphere::distm) or is 'euclidean' or 'manhattan'
  if (startsWith(method, "dist")) {
    # Extract the function name from the 'method' string (e.g., "distHaversine")
    dist_fun <- get(method, envir = asNamespace("geosphere"))
    # Assuming `data` has coordinates (e.g., latitude, longitude) in columns 1 and 2
    dist <- geosphere::distm(data, fun = dist_fun)
    rownames(dist) <- rownames(data)
    colnames(dist) <- rownames(data)
    
    # Extract the 5-9 characters of the method for the filename "meth_fix"
    meth_fix <- tolower(substr(method, 5, 9))
    
  } else if ((method == "euclidean" || method == "manhattan") && any(is.na(data))) {
    # If method is 'euclidean' or 'manhattan' and there are NAs, use dist() instead of vegdist()
    dist <- dist(data, method = method)  
  } else {
    # Otherwise, use vegdist() as standard
    dist <- vegdist(data, method = method)
  }
  
  # Step 3: Convert the distance matrix to a full matrix
  matrix_full <- as.matrix(dist)
  
  # Step 4: Set the upper triangle and diagonal to NA
  matrix_full[upper.tri(matrix_full)] <- NA 
  diag(matrix_full) <- NA
  
  # Step 5: Convert the distance matrix into long format
  matrix_long <- as.data.frame(as.table(matrix_full))
  
  # Step 6: Rename columns for clarity
  colnames(matrix_long) <- c("id", "id2", paste0(prefix, "_", meth_fix))
  
  # Step 7: Filter out rows where the distance is NA
  matrix_long <- matrix_long %>%
    filter(!is.na(get(paste0(prefix, "_", meth_fix))))
  
  # Step 8: Save the long format data to CSV
  write.csv(matrix_long,
            file.path(dir_data, paste0(prefix, "_", meth_fix, "_df.csv")), 
            row.names = FALSE)
  
  # Return the long format dataframe
  return(matrix_long)  
}