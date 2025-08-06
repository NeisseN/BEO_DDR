# Meta data for BExIS

# Niklas Neisse
# 2025.07.04


getwd()

# Packages
library(tidyverse)
library(janitor)

# Directories
dir_wd   <- file.path('C:', 'code', 'BEO_DDR')
dir_data <- file.path(dir_wd,'data')


# Data
df_rarefied <- readRDS(file.path(dir_data, 'rarefied.RDS'))

df_meta <- df_rarefied@sam_data %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  clean_names() %>% 
  remove_rownames() %>% 
  dplyr::select(-c(
    'elevation', 'slope_d', 'slope_p', 'aspect', 'soil_texture', 'comments', 
    'op', 'l_ammonium', 'l_nitrate', 'l_phosphate', 'l_sulfate', 'l_calcium',
    'l_magnesium', 'l_potassium', 'l_phosphorus', 'l_sulfur', 'rs18', 'ts18',
    'ms18', 'rs19', 'ts19', 'ms19', 'grazing', 'mowing', 'fertilization')) %>% 
  rename(ph = p_h)
  

df_taxa <- df_rarefied@tax_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  clean_names() %>% 
  rownames_to_column(var = "id")

df_otu <- df_rarefied@otu_table %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>% 
  clean_names() %>% 
  rownames_to_column(var = "id")

df_asv <- df_taxa %>% 
  full_join(df_otu, by = 'id')


write.csv(df_meta, row.names = F,
          file.path(dir_data, 'df_bexis_meta.csv'))
write.csv(df_asv , row.names = F,
          file.path(dir_data, 'df_bexis_asv.csv'))
