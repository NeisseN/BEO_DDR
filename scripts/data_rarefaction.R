# Rarefaction

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
.cran_packages <- c('tidyverse','phyloseq') # 'tidyverse','vegan','ggplot2','nlme','viridis'

.inst <- .cran_packages %in% installed.packages()
# Check if CRAN packages are installed
if(any(!.inst)) {
  # Install missing CRAN packages
  install.packages(.cran_packages[!.inst]) 
}

# Load required packages
sapply(.cran_packages, require, character.only = TRUE)


# Directories ------------------------------------------------------------------
getwd()
dir_wd   <- file.path('C:', 'code', 'eddsbc')
dir_data <- file.path(dir_wd, 'data')
dir_out  <- file.path(dir_wd, 'output')


# Data -------------------------------------------------------------------------
physeq <- readRDS(file.path(dir_data, 'phyloseq.rds'))
# How many samples? & How many ASVs?
data.frame(Samples = nsamples(physeq), ASVs = ntaxa(physeq))


reads_per_sample <- data.frame(nreads = sort(sample_sums(physeq), T), 
                               samples = 1:nsamples(physeq),
                               type = "Samples")
### ReadsPerSample ###
ReadsPerSample <- ggplot() +
  geom_bar(data=reads_per_sample, 
           aes(x = samples, y = nreads),
           stat="identity", fill="#238A8DFF") +
  geom_abline(intercept = 23942, color = "#DCE318FF", linewidth = 5) +
  xlab("Ordered Samples") +
  scale_y_continuous(breaks = c(0,23945,40000,80000,120000)) +
  ylab("Number of Reads") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        axis.text.x = element_text(size = 26, color = 'black'),
        axis.text.y = element_text(size = 26, angle = 60, hjust = 1, vjust = -.3, color = 'black'),
        axis.ticks.x = element_line(size = 2, lineend = "butt"),
        axis.ticks.y = element_line(size = 2, lineend = "butt"),
        axis.ticks.length  = unit(0.4, "cm"),
        axis.line = element_line(size = 1.1))
ReadsPerSample
ggsave(file.path(dir_out, "ReadsPerSample.png"), 
       plot = ReadsPerSample, width = 10, height = 9)


# Rarefaction ------------------------------------------------------------------
# the lowest is SEG2_5 with 2220 and then it jumps to 23942
rarefied <- rarefy_even_depth(physeq, sample.size = 23942, rngseed=1) # on set.seed(1)
rarefied@sam_data$id <- as.factor(rarefied@sam_data$id)
rarefied@sam_data$loc <- as.factor(rarefied@sam_data$loc)
rarefied@sam_data$EP_Plotid <- as.factor(rarefied@sam_data$EP_Plotid)
# 213 ASVs were removed because they are no longer present in any sample after random subsampling
saveRDS(rarefied, file.path(dir_data, "rarefied.RDS"))



reads_per_ASV <- data.frame(nreads = sort(taxa_sums(physeq), T), SVs = 1:ntaxa(physeq))
reads_per_ASV.raref <- data.frame(nreads = sort(taxa_sums(rarefied), T), SVs = 1:ntaxa(rarefied))
AsvsReadsComparison <- ggplot() +
  geom_area(data = reads_per_ASV, 
            aes(x = SVs, y = nreads, fill="Original"), 
            stat="identity") + 
  geom_area(data = reads_per_ASV.raref, 
            aes(x = reads_per_ASV.raref$SVs, y = reads_per_ASV.raref$nreads, fill="Rarefied"), 
            stat="identity") +
  scale_fill_manual(values = c("Original" = "#238A8DFF", "Rarefied" = "#DCE318FF")) + 
  scale_x_log10() +
  xlab("Ordered ASVs (log)") +
  ylab("Number of Reads") +
  ylim(0, 30000) + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        axis.text.x = element_text(size = 26, color = 'black'),
        axis.text.y = element_text(size = 26, angle = 60, hjust = 1, vjust = -.3, color = 'black'),
        axis.ticks.x = element_line(size = 2, lineend = "butt"),
        axis.ticks.y = element_line(size = 2, lineend = "butt"),
        axis.ticks.length  = unit(0.4, "cm"),
        axis.line = element_line(size = 1.1),
        legend.title = element_blank(),  # Hide legend title
        legend.text = element_text(size = 26),  # Customize legend text size
        legend.key.size = unit(0.8, "cm"), # Customize legend key size)
        legend.position = c(0.8, 0.8))
AsvsReadsComparison
ggsave(file.path(dir_out, "AsvsReadsComparison.png"), 
       plot = AsvsReadsComparison, width = 10, height = 9)


rarefaction_data <- rarecurve(
  as.matrix(rarefied@otu_table@.Data), 
  step = 1000, sample.size = 23942, label = F)


