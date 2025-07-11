# Data distances

# Author: Niklas Neiße*, Stephanie Jurburg
# Mail  : neisse.n@protonmail.com
# Date  : 2025.05.21

# Distance matrix for of each sample-level variable 


# Setup ------------------------------------------------------------------------
# rm(list = ls())
set.seed(100)


# Packages ---------------------------------------------------------------------
# Define CRAN packages
.cran_packages <- c('tidyverse', 'vegan', 'corrplot', 'geosphere') 
.bioc_packages <- c('phyloseq')  

# Install missing CRAN packages
.inst <- .cran_packages %in% installed.packages() 
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst]) 
}

# Check if Bioconductor packages are installed
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  # Source Bioconductor installation script
  source('http://bioconductor.org/biocLite.R')
  # Install missing Bioconductor packages without asking for confirmation
  biocLite(.bioc_packages[!.inst], ask = F)
}

# Load required packages
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)


# Directories ------------------------------------------------------------------
getwd()
dir_wd   <- file.path('C:', 'code', 'BEO_DDR')
dir_data <- file.path(dir_wd, 'data')
dir_fun  <- file.path(dir_wd, 'helper_functions')
dir_resu <- file.path(dir_wd, 'output')


# Functions --------------------------------------------------------------------
# Distance
source(file.path(dir_fun, 'fun_distance.R'))


# Data -------------------------------------------------------------------------
# Rarefied phyloseq object
rarefied <- readRDS(file.path(dir_data, 'rarefied.RDS'))

# Meta data
rare_meta <- as.data.frame(as.matrix(rarefied@sam_data))
names(rare_meta)


# Sample coordinate data -------------------------------------------------------
# Coordinates df
coordinates_df <- rare_meta[c('lon', 'lat')] %>% 
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat))

coor_have_df <- fun_distance(
  data = coordinates_df, method = 'distHaversine',
  prefix = 'coor', dir_data = dir_data
)


# Bacterial relative abundance Bray-Curtis -------------------------------------
# Bacterial relative abundance data 
rare_abund <- as.data.frame(as.matrix(rarefied@otu_table))

# Bacterial relative abundance Bray-Curtis
bac_bray_df <- fun_distance(
  data = rare_abund, method = 'bray', 
  prefix = 'bac', dir_data = dir_data)


# Plant relative abundance Bray-Curtis -----------------------------------------
# Plant abundance df
plant_abund <- rare_meta %>% 
  dplyr::select(
    grep('^[a-zA-Z]{4}\\.[a-zA-Z]{4}$', names(rare_meta), value = TRUE)
  )

# Replace all NAs with 0, biologically: NA here means not present
plant_abund[is.na(plant_abund)] <- 0

# Convert all columns to numeric
plant_abund <- data.frame(lapply(plant_abund, as.numeric))

# Define the samples as rownames
rownames(plant_abund) <- rownames(rare_meta)

plant_bray_df <- fun_distance(
  data = plant_abund, method = 'bray', 
  prefix = 'plant', dir_data = dir_data)


# Scaling ----------------------------------------------------------------------
# A note on scaling:
# When you use one variable to explain or correlate with another, 
# the relationship (or correlation) between them does not change 
# if you scale (or standardize) the variables. 
# What changes is the interpretation of the effect size.
# Scaling is often done when comparing multiple predictors with different units, 
# as it helps normalize the effect sizes 
# and makes them more interpretative in a comparable way


# Plant-traits -----------------------------------------------------------------
# Plant-trait df
plant_traits_all_df <- rare_meta %>% 
  dplyr::select(
    c('id', 'RawBioM', 'BioMass_g', 'BioMass_g_m2', 'tdw.g',
      'h.veg_cm', 'h.reg_cm', 'LMA', 'SLA',
      'COV_Bare', 'COV_Litter', 'COV_Veg', 'COV_Moss', 
      'COV_Living', 'COV_Senesc',
      'COV_Forb', 'COV_Gras', 'COV_Leg')) %>% 
  mutate(id = as.factor(id))
plant_traits_all_df[2:ncol(plant_traits_all_df)] <- lapply(
  plant_traits_all_df[2:ncol(plant_traits_all_df)],
  as.numeric
)

corrplot(cor(plant_traits_all_df[c(2:ncol(plant_traits_all_df))]), 
         type = 'upper', method = 'square', diag = FALSE, 
         addCoef.col = 'darkgrey', number.cex = .7,
         tl.srt = 45, tl.cex = 1, tl.offset = 1)

# Here is my keep suggestion: 
# Biomass representing productivity
# 
plant_traits_df <- plant_traits_all_df %>% 
  dplyr::select(
    c('BioMass_g',
      # LMA: mass of a leaf per unit area in grams per square meter (g/m²) VS.
      # SLA: reciprocal LMA, leaf area per unit leaf mass, in square meters per gram (m²/g). 
      'SLA',
      'COV_Bare', 'COV_Litter', 'COV_Moss', 
      'COV_Senesc'))

corrplot(cor(plant_traits_df[c(2:ncol(plant_traits_df))]), 
         type = 'upper', method = 'square', diag = FALSE, 
         addCoef.col = 'darkgrey', number.cex = .7,
         tl.srt = 45, tl.cex = 1, tl.offset = 1)


# Plant trait euclidean distance
p_trait_eucl_df <- fun_distance(
  data = plant_traits_df, 
  method = 'euclidean', scale_data = T,
  prefix = 'p_trait', dir_data = dir_data
)

# Biomass manhattan
pt_biomass_manh_df <- fun_distance(
  data = plant_traits_df[c('BioMass_g')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pt_biomass', dir_data = dir_data
)

# SLA manhattan
pt_sla_manh_df <- fun_distance(
  data = plant_traits_df[c('SLA')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pt_sla', dir_data = dir_data
)

# COV_Bare manhattan
pt_cov_bare_manh_df <- fun_distance(
  data = plant_traits_df[c('COV_Bare')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pt_cov_bare', dir_data = dir_data
)

# COV_Litter manhattan
pt_cov_litter_manh_df <- fun_distance(
  data = plant_traits_df[c('COV_Litter')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pt_cov_litter', dir_data = dir_data
)

# COV_Moss manhattan
pt_cov_moss_manh_df <- fun_distance(
  data = plant_traits_df[c('COV_Moss')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pt_cov_moss', dir_data = dir_data
)

# COV_Senesc manhattan
pt_cov_senesc_manh_df <- fun_distance(
  data = plant_traits_df[c('COV_Senesc')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pt_cov_senesc', dir_data = dir_data
)


# Physio-chemical distances ----------------------------------------------------
# Physio-chemical df
physchem_df <- rare_meta %>% 
  dplyr::select(c('C_.', 'N_.', 'Soil_moisture', 'pH')) %>% 
  rename(C = 'C_.', N = 'N_.')
physchem_df[] <- lapply(physchem_df[], as.numeric)  
physchem_df <- physchem_df %>% mutate(cnr = C/N)

# Physio-chemical euclidean distance
physchem_eucl_df <- fun_distance(
  data = physchem_df, method = 'euclidean', scale_data = T,
  prefix = 'physchem', dir_data = dir_data
)

# C manhattan
pc_c_manh_df <- fun_distance(
  data = physchem_df[c('C')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pc_c', dir_data = dir_data
)

# N manhattan
pc_n_manh_df <- fun_distance(
  data = physchem_df[c('N')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pc_n', dir_data = dir_data
)

# Soil_moisture manhattan
pc_soil_moist_manh_df <- fun_distance(
  data = physchem_df[c('Soil_moisture')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pc_soil_moist', dir_data = dir_data
)

# pH manhattan
pc_ph_manh_df <- fun_distance(
  data = physchem_df[c('pH')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pc_ph', dir_data = dir_data
)

# cnr manhattan
pc_cnr_manh_df <- fun_distance(
  data = physchem_df[c('cnr')], 
  method = 'manhattan', scale_data = T,
  prefix = 'pc_cnr', dir_data = dir_data
)


# Merge distance dfs -----------------------------------------------------------
distances_df <- merge(coor_have_df, bac_bray_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, plant_bray_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, p_trait_eucl_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, physchem_eucl_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pt_biomass_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pt_cov_bare_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pt_cov_litter_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pt_cov_moss_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pt_cov_senesc_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pt_sla_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pc_c_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pc_cnr_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pc_n_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pc_ph_manh_df, by = c('id', 'id2'), all = TRUE)
distances_df <- merge(distances_df, pc_soil_moist_manh_df, by = c('id', 'id2'), all = TRUE)


# Binned distances -------------------------------------------------------------
# 3 categories of distances: within plot, within region, among region 
# Create a function to check the pattern
pattern_check <- function(id1, id2) {
  # Extract the part before the first underscore
  prefix1 <- sub('_.*', '', id1)
  prefix2 <- sub('_.*', '', id2)
  
  # Extract the first 3 letters
  first3_letters1 <- substr(id1, 1, 3)
  first3_letters2 <- substr(id2, 1, 3)
  
  # Check for the patterns and return corresponding value
  if (prefix1 == prefix2) {
    # If prefixes are the same
    return(1)  
  } else if (first3_letters1 == first3_letters2) {
    # If first 3 letters are the same
    return(2)  
  } else {
    # If both differ
    return(3)  
  }
}

distances_df <- distances_df %>%
  mutate(dist_bin = mapply(pattern_check, id, id2)) %>%
  dplyr::select(id, id2, dist_bin, everything())


# Save the data ----------------------------------------------------------------
write.csv(distances_df, file.path(dir_data, 'distances_df.csv'),
          row.names = F)

