# Analysis distances decay relationship linear

# Author: Niklas Neiße*, Stephanie Jurburg
# Mail  : neisse.n@protonmail.com
# Date  : 2025.05.25

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
dir_wd   <- file.path('C:', 'code', 'eddsbc')
dir_data <- file.path(dir_wd, 'data')
dir_fun  <- file.path(dir_wd, 'helper_functions')
dir_resu <- file.path(dir_wd, 'output')


# Functions --------------------------------------------------------------------
# Distance
source(file.path(dir_fun, 'fun_distance.R'))


# Data -------------------------------------------------------------------------