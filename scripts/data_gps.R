# GPS coordinates

# Author: Niklas Neiße*, Stephanie Jurburg
# Mail  : neisse.n@protonmail.com
# Date  : 2025.05.21


# Setup ------------------------------------------------------------------------
# # Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100)


# Packages ---------------------------------------------------------------------
# Define CRAN packages
.cran_packages <- c('tidyverse', 'geosphere', 'sp')
.inst <- .cran_packages %in% rownames(installed.packages())
if (any(!.inst)) install.packages(.cran_packages[!.inst])
sapply(.cran_packages, require, character.only = TRUE)


# Directories ------------------------------------------------------------------
dir_data <- file.path('data') # lowercase 'data'
getwd()  # just to confirm current working directory


# Data -------------------------------------------------------------------------
UTM.raw <- read_csv(file.path(dir_data, 'Potsdam23.GPS_UTM.csv'))


# Coordinate transformation ----------------------------------------------------
# Select UTM coordinates
UTM <- UTM.raw %>% select(X_coord, Y_coord)

# Coordinate transformation (UTM to WGS84 decimal degrees)
utm <- SpatialPoints(UTM, proj4string = CRS('+proj=utm +zone=32 +datum=WGS84'))
coordd <- sp::spTransform(utm, CRS('+proj=longlat +datum=WGS84'))
coordd <- as.data.frame(coordd) %>%
  rename(lon = coords.x1, lat = coords.x2) %>%
  select(lat, lon)

# Export
write.csv(coordd, row.names = FALSE, 
          file.path(dir_data, 'coordinates_decimal_degrees.csv'))


# Distance matrix --------------------------------------------------------------
# Calculate distance matrix using Haversine formula
distmat_m <- as.data.frame(geosphere::distm(coordd, fun=distHaversine)) 
# ?distm  # Display documentation for distm function

# Create metadata dataframe
meta <- data.frame(
  'QNR' = rep(c(2,5,9),18),
  'L'   = c(rep('HEG',27), rep('SEG',27)),
  'PNr'=  rep(c(
    6, 7, 8, 16, 20, 22, 25, 41, 42, 2, 3, 5, 7, 11, 14, 16, 18, 37),
    each = 3)) %>% 
  unite(L,PNr, col = 'EP_Plotid', sep = '') %>% 
  mutate(ep = EP_Plotid, qnr = QNR) %>% 
  unite(ep,qnr, col='id', sep = '_') %>%
  select(id, EP_Plotid, QNR) 

distmat_m <- distmat_m %>% 
  mutate(id = meta$id, .before = 1) # Add 'id' column to distance matrix
names(distmat_m) <- c('id', meta$id) # Rename columns of distmat_m
rownames(distmat_m) <- meta$id
write.csv(distmat_m, row.names = F,  
          file.path(dir_data, 'distance_matrix_in_meters.csv'))
distmat_m <- read.csv(file.path(dir_data, 'distance_matrix_in_meters.csv'))

distmatlong <- distmat_m %>% 
  pivot_longer(
    cols = !id, names_to = 'id2', values_to = 'dist_m') # Convert distance matrix to long format
write.csv(distmatlong, file.path(dir_data, 'distance_matrix_long_format.csv'))
