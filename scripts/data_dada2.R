# DADA2: Reads 2 Comunitay Analyses

# Author: Niklas Neiße*, Stephanie Jurburg
# Mail  : neisse.n@protonmail.com
# Date  : 2025.05.21

# Bioinformatic pipeline by Callahan et al. (2016)
#  handling raw microbial sequence reads 


# Setup ------------------------------------------------------------------------
# # Remove all objects in the global environment
# rm(list = ls()) 
# Set seed for reproducibility
set.seed(100) 


# Packages ---------------------------------------------------------------------
# Define CRAN packages
.cran_packages <- c('tidyverse','gridExtra','viridis','cowplot') 
# Define Bioconductor packages
.bioc_packages <- c('dada2', 'phyloseq', 'DECIPHER', 'phangorn') 

# Check if CRAN packages are installed
.inst <- .cran_packages %in% installed.packages() 
if(any(!.inst)) {
  # Install missing CRAN packages
  install.packages(.cran_packages[!.inst]) 
}

# Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install missing Bioconductor packages
.inst <- .bioc_packages %in% rownames(installed.packages())
if (any(!.inst)) {
  BiocManager::install(.bioc_packages[!.inst], ask = FALSE)
}

# Load required packages
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE) 


# Directories ------------------------------------------------------------------
getwd() 
dir_data <- file.path('data')
dir_seq  <- file.path(dir_data, 'MiSeq.16S.29082023')
# Define the path for filtered files
dir_filt <- file.path(dir_seq, 'filtered')
dir_out  <- file.path('output')


# Data -------------------------------------------------------------------------
# Get forward sequence files
fnFs <- sort(list.files(dir_seq, pattern='_L001_R1_001.fastq')) 
# Get reverse sequence files
fnRs <- sort(list.files(dir_seq, pattern='_L001_R2_001.fastq')) 

sample.names <- sapply(strsplit(fnFs,'_L001_R1_001.fastq.gz'), `[`, 1) # Extract sample names

fnFs <- file.path(dir_seq, fnFs) # Create full paths for forward sequence files
fnRs <- file.path(dir_seq, fnRs) # Create full paths for reverse sequence files

filtFs <- file.path(dir_filt, paste0(sample.names, '_1.filtered.fastq.gz')) # Define full paths for filtered forward sequences
filtRs <- file.path(dir_filt, paste0(sample.names, '_2.filtered.fastq.gz')) # Define full paths for filtered reverse sequences


# Read quality -----------------------------------------------------------------
Forwardpqp.A  <- plotQualityProfile(fnFs[1]) # Plot quality profile for the first forward sequence
Reversepqp.A  <- plotQualityProfile(fnRs[1]) # Plot quality profile for the first reverse sequence
Qualityprof.A <- cowplot::plot_grid(Forwardpqp.A, Reversepqp.A, ncol=2) # Combine quality profiles of the first sample
Forwardpqp.O  <- plotQualityProfile(fnFs[54]) # Plot quality profile for the 54th forward sequence
Reversepqp.O  <- plotQualityProfile(fnRs[54]) # Plot quality profile for the 54th reverse sequence
Qualityprof.O <- cowplot::plot_grid(Forwardpqp.O, Reversepqp.O, ncol=2) # Combine quality profiles of the 54th sample
Qualityprof   <- cowplot::plot_grid(
  Forwardpqp.A, Reversepqp.A, Forwardpqp.O, Reversepqp.O, ncol=2, nrow = 2) # Combine quality profiles of all samples
plot(Qualityprof) # Plot combined quality profiles
QualityRawPerRead <- ggplot2::last_plot() # Save the plot as 'QualityRawPerRead'
ggsave(file.path(dir_out, 'QualityRawPerRead.png'), 
       plot = QualityRawPerRead, width = 10, height = 10) # Save the plot as an image


# Filter and trim the sequences ------------------------------------------------
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs, 
  # Most important parameters!
  truncLen = c(230,220), trimLeft = 10, maxN = 0, maxEE = c(3,3), 
  truncQ = 2, rm.phix = T, compress = T, multithread = T) 

write.csv2(out, file = file.path(dir_data, 'filtered.csv'))
out <- read_delim(
  file.path(dir_data, 'filtered.csv'),
  delim = ';', escape_double = FALSE, trim_ws = TRUE)


# Error Model ------------------------------------------------------------------
#  this deviates from Callaham 2016
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
write.csv2(errF, file = file.path(dir_data, 'ErrMF'))
write.csv2(errR, file = file.path(dir_data, 'ErrMR'))

plotErrors(errF, nominalQ=TRUE) # Plot forward errors
plotErrors(errR, nominalQ=TRUE) # Plot reverse errors



# Infer ASVs -------------------------------------------------------------------
# ?dada
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)



# Merge inferred ASVs ----------------------------------------------------------
mergers <- mergePairs(
  dadaFs, filtFs, dadaRs, filtRs, verbose = T, minOverlap = 10) 

seqtab <- makeSequenceTable(mergers) # Create a sequence table
dim(seqtab) # Check dimensions of the sequence table
colnames(seqtab) # Check column names of the sequence table
row.names(seqtab) # Check row names of the sequence table



# Chimera ----------------------------------------------------------------------
# Remove chimeras from the sequence table
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method = 'consensus', multithread = T, verbose = T) 
# Identified 199 bimeras out of 12758 input sequences.

# What percentage of reads are likely chimeras?
1- (sum(seqtab.nochim)/sum(seqtab)) 
# 0.006119889
# What percentage of ASVs are likely chimeras?
1- (dim(seqtab.nochim)[2]/dim(seqtab)[2])
# 0.01559806

saveRDS(seqtab.nochim, file.path(dir_data, 'seqtab.nochim.RDS'))


# Lost reads -------------------------------------------------------------------
getN <- function(x) sum(getUniques(x)) # Define a function to count unique reads
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim)) # Create a tracking table
colnames(track) <- c('id', 'input', 'filtered',
                       'denoisedF', 'denoisedR', 'merged', 'nonchim') # Set column names
rownames(track) <- sample.names # Set row names
head(track) # Display the first few rows of the tracking table
write.table(track, 'Data/tracking.table.txt', sep='\t') # Write the tracking table to a text file
track <- as.data.frame(track) # Convert tracking table to data frame
track$input/track$nonchim # Calculate the ratio of input reads to non-chimeric reads
mean(track$input/track$nonchim) # Calculate the mean of the input-to-non-chimeric read ratio

viridis_pal()(6) # Generate a color palette
names(track) # Get column names of the tracking table
track$id <- as.factor(track$id) # Convert 'id' column to factor

# Create a plot visualizing lost reads
lost_reads <- ggplot(track) +
  geom_bar(aes(x = id, y = input, fill = 'Input'), stat = 'identity', position = 'stack', color = 'black') +
  geom_bar(aes(x = id, y = filtered, fill = 'Filtered'), stat = 'identity', position = 'stack', color = 'black') +
  geom_bar(aes(x = id, y = denoisedF, fill = 'Denoised Forward'), stat = 'identity', position = 'stack', color = 'black') +
  geom_bar(aes(x = id, y = denoisedR, fill = 'Denoised Reverse'), stat = 'identity', position = 'stack', color = 'black') +
  geom_bar(aes(x = id, y = merged, fill = 'Merge'), stat = 'identity', position = 'stack', color = 'black') +
  geom_bar(aes(x = id, y = nonchim, fill = 'Non-Chimeric'), stat = 'identity', position = 'stack', color = 'black') +
  scale_fill_manual(values = c(
    Input              = '#440154FF',
    Filtered           = '#414487FF',
    `Denoised Forward` = '#2A788EFF',
    `Denoised Reverse` = '#22A884FF',
    Merge              = '#7AD151FF',
    `Non-Chimeric`     = '#FDE725FF'), 
    breaks = c('Input', 'Filtered', 'Denoised Forward', 'Denoised Reverse', 
               'Merge', 'Non-Chimeric')) +
  labs(fill = 'Filter Process Status',
       y = 'Number of Reads',
       x = 'Sample') +
  theme_classic() +
  theme(axis.text   = element_text(size = 16),    # Increase axis text size
        axis.title  = element_text(size = 20),   # Increase axis title size
        axis.line   = element_line(linewidth = 1.1),   # Increase axis line thickness
        axis.ticks  = element_line(linewidth = 1.5),  # Increase axis tick thickness
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3, color = 'black'),
        axis.text.y = element_text(angle = 45, hjust = 1, color = 'black'),
        legend.title    = element_text(size = 18),  # Increase legend title size
        legend.text     = element_text(size = 16), # Increase legend text size
        legend.position = c(.8, .85)) +  
  scale_x_discrete(labels = 1:54) +  # Set x-axis labels from 1 to 54
  guides(fill = guide_legend(reverse = F))  # Reverse the order of items in the legend
lost_reads

ggsave(file.path(dir_out, 'LostReads.png'), 
       plot = lost_reads, width = 10, height = 9) # Save the plot as an image


# Assign taxonomy --------------------------------------------------------------
taxa <- assignTaxonomy(
  seqtab.nochim, file.path(dir_data, 'silva_nr99_v138.1_train_set.fa.gz'),
  multithread = TRUE, tryRC = TRUE) 
# Clean up this object
taxa <- as.data.frame(taxa) 
taxa[] <- lapply(taxa, as.character)

saveRDS(taxa, file.path(dir_data, 'taxtable.RDS'))

