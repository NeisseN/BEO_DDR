#### Taxa Range ####
### Neisse, Niklas 
## 22052024 

#### Setting the stage ####
rm(list = ls()) # Remove all objects in the global environment

### Packages ###
.cran_packages <- c('GeoRange','tidyverse','viridis') # Define CRAN packages
# , 'cowplot' , 'BiocManager', 'iNEXT', 'gridExtra' , 'knitr','vegan', 'nlme', 'corrplot','skimr','corrplot'
.bioc_packages <- c('phyloseq') # Define Bioconductor packages
#'dada2', 'phyloseq', 'DECIPHER', 'phangorn','BiocStyle', 'vegetarian'

.inst <- .cran_packages %in% installed.packages() # Check if CRAN packages are installed
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst]) # Install missing CRAN packages
}

.inst <- .bioc_packages %in% installed.packages() # Check if Bioconductor packages are installed
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R") # Source Bioconductor installation script
  biocLite(.bioc_packages[!.inst], ask = F) # Install missing Bioconductor packages without asking for confirmation
}

sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE) # Load required packages
set.seed(100) # Set seed for reproducibility
vpal18 <- viridis_pal()(18)

#### Data ####
rarefied <- readRDS("./data/rarefied.RDS")

###  Convex Hull area calculation
## Uses the cylindrical equal area projection in order to check if the minimum convex hull wraps around the prime meridian
CHAdf <- data.frame(asv = colnames(rarefied@otu_table))
for (i in 1:ncol(rarefied@otu_table)) {
  sub.nr <- which(rarefied@otu_table[,i]>0)
  CHAdf$chull_areas[i] <- CHullArea(rarefied@sam_data$lon[sub.nr],
                                 rarefied@sam_data$lat[sub.nr])
}

### GCD 
## Calculates the maximum pairwise great circle distance from a set of decimal degree coordinates
GCDdf <- data.frame(asv = colnames(rarefied@otu_table))
for (i in 1:ncol(rarefied@otu_table)) {
  sub.nr <- which(rarefied@otu_table[,i]>0)
  GCDdf$gcd[i] <- GCD(rarefied@sam_data$lon[sub.nr],
                   rarefied@sam_data$lat[sub.nr])
}
GCDdf$gcd[is.na(GCDdf$gcd)] <- 0

### Prevalence
## Proportion of sites in which the species was recorded as present
Prevdf <- data.frame(asv = colnames(rarefied@otu_table))
for (i in 1:ncol(rarefied@otu_table)) {
  Prevdf$prev[i] <- sum(rarefied@otu_table[,i] != 0)/nrow(rarefied@otu_table[,i])
}

### Abundance
## Number of individuals devided by the rarefaction depth of about 24000 and divided by the number of samples.
Abundf <- data.frame(asv = colnames(rarefied@otu_table))
Abundf$abun <- colSums(rarefied@otu_table)/23942/nrow(rarefied@otu_table)

### Taxa
Taxadf <- data.frame(rarefied@tax_table)
Taxadf$asv <- rownames(Taxadf)

### Big Merge 
dfs <- list(Abundf, Prevdf, CHAdf, GCDdf, Taxadf) # List of data frames to merge
df <- Reduce(function(x, y) merge(x, y, by = "asv", all.x = TRUE), dfs) # Perform the merges

### Dominate Phyla
# Calculate the dominant phyla and update the df
dominant_phyla <- df %>%
  group_by(Phylum) %>%
  summarise(Sum = sum(abun)) %>%
  arrange(desc(Sum)) %>%
  slice(1:8) %>%
  pull(Phylum)

# Update df with d9phyl column and convert it to a factor with ordered levels
df <- df %>%
  mutate(d8phyl = ifelse(Phylum %in% dominant_phyla, Phylum, NA),
         d8phyl = factor(d8phyl, levels = dominant_phyla))

write.csv(df, file = 'data/TaxaRange.df')
df <- read.csv('data/TaxaRange.df')

#### Plots #### 
## Abundance by Prevalence
ggplot() +
  geom_point(data = df, aes(x = prev, y= abun)) +
  labs(x = "Rel. Abundace", y = "Prevalence") +
  theme_classic()

## Convex Hull Area by Prevalence
ggplot() +
  geom_point(data = df, aes(x = prev, y= chull_areas)) +
  labs(x = "Convex Hull Area", y = "Prevalence") +
  theme_classic()

## Convex Hull Area by Abundance 
ggplot() +
  geom_point(data = df, aes(x = abun, y= chull_areas)) +
  labs(x = "Rel. Abundace", y = "Convex Hull Area") +
  theme_classic()

ggplot() +
  geom_point(data = df, aes(x = abun, y= chull_areas, 
                            colour = factor(d8phyl, levels = dominant_phyla))) +
  scale_color_manual(values = viridis_pal()(8)) +
  labs(x = "Rel. Abundace", y = "Convex Hull Area", colour = 'Phylum') +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 10, angle = 60, vjust = .6, color = 'black'),
        axis.text.y = element_text(size = 10, angle = 60, hjust = 1, vjust = -.3, color = 'black'),
        axis.ticks.x = element_line(linewidth = 2, lineend = "butt"),
        axis.ticks.y = element_line(linewidth = 2, lineend = "butt"),
        axis.ticks.length  = unit(0.1, "cm"),
        axis.line = element_line(linewidth = 1.1),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10),  # Customize legend text size
        legend.key.size = unit(0.6, "cm")) + # Customize legend key size
  facet_wrap(~d8phyl)

## GCD by Abundance
ggplot() +
  geom_jitter(data = df, aes(x = abun, 
                             y= log(gcd + 0.001))) + # + a small value to avoid -inv
  labs(x = "Rel. Abundace", y = "GCD log transformed", colour = 'Phylum') +
  theme_classic()

ggplot() +
  geom_jitter(data = df, aes(x = abun, y= log(gcd+0.001), 
                             colour = factor(d8phyl, levels = dominant_phyla))) +
  scale_color_manual(values = viridis_pal()(8)) +
  labs(x = "Rel. Abundace", y = "GCD log transformed", colour = 'Phylum') +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 10, angle = 60, vjust = .6, color = 'black'),
        axis.text.y = element_text(size = 10, angle = 60, hjust = 1, vjust = -.3, color = 'black'),
        axis.ticks.x = element_line(linewidth = 2, lineend = "butt"),
        axis.ticks.y = element_line(linewidth = 2, lineend = "butt"),
        axis.ticks.length  = unit(0.1, "cm"),
        axis.line = element_line(linewidth = 1.1),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10),  # Customize legend text size
        legend.key.size = unit(0.6, "cm")) + # Customize legend key size
  facet_wrap(~d8phyl)


#### Anova ####

hist(df$abun)
hist(log(df$abun))
qqnorm(log(df$abun), main='Normal')
qqline(log(df$abun))

#shapiro.test(log(df$abun))
ks.test(log(df$abun), 'pnorm')

kruskal.test(df$abun ~ df$Phylum)
#pairwise.wilcox.test(df$abun, df$Phylum, p.adjust.method = "bonferroni") 
