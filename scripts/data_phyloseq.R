# Phyloseq object

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
.cran_packages <- c('tidyverse','phyloseq')

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
dir_data <- file.path('data')


# Data -------------------------------------------------------------------------
# non-chimeric sequences
seqtab.nochim <- readRDS(file.path(dir_data, 'seqtab.nochim.RDS'))
# taxonomy
taxa          <- readRDS(file.path(dir_data, 'taxtable.RDS'))

#length of the first read
nchar(colnames(seqtab.nochim)[1])
# length of all reads
hist(nchar(colnames(seqtab.nochim)))

# create a variable that contains all the sequences.
seqs <- colnames(seqtab.nochim)
# and convert the sequences into a format that phyloseq will recognize as sequences
seqs <- Biostrings::DNAStringSet(seqs)


# Data: Meta -------------------------------------------------------------------
### This is what we've got:
# 1. seqtab.nochim: 1st col  : sample ID; 
#                   12558 col: seq. reads as counts
# 2. taxa         : assigned bac's form data bank; one row one bac. 
# 3. track        : counts of reads as the clean up progresses

# 4. Meta
meta <- data.frame('QNR'= rep(c(2,5,9),18),
                   'L'  = c(rep('HEG',27),rep('SEG',27)),
                   'loc'= c(rep('HEG',27),rep('SEG',27)),
                   'PNr'= c(rep(6,3), 
                            rep(7,3), 
                            rep(8,3), 
                            rep(16,3), 
                            rep(20,3), 
                            rep(22,3), 
                            rep(25,3), 
                            rep(41,3), 
                            rep(42,3), 
                            rep(2,3), 
                            rep(3,3), 
                            rep(5,3), 
                            rep(7,3), 
                            rep(11,3), 
                            rep(14,3), 
                            rep(16,3), 
                            rep(18,3), 
                            rep(37,3))) %>% 
  unite(L, PNr, col = 'EP_Plotid', sep = '') %>% 
  mutate(ep = EP_Plotid, qnr = QNR) %>% 
  unite(ep, qnr, col = 'id', sep = '_') %>% 
  dplyr::select(id, loc, EP_Plotid, QNR)
row.names(seqtab.nochim) <- meta$id
row.names(meta) <- row.names(seqtab.nochim)


# Data: Plant ------------------------------------------------------------------
# 5.1 p.traits: 
p.traits <- read_csv(file.path(dir_data, 'Potsdam23.PlantTraits.csv')) %>%
  rename(LA = 'LA_mm²',
         SLA = 'SLA_mm²/mg') %>% 
  dplyr::select(-c(1:3,5,7,9,11,12,24)) %>% 
  group_by('EP', 'QNR') %>% 
  group_by(EP,QNR) %>% 
  summarise(tdw.g = mean(TDW_g, na.rm=T), #Total tiller dry weight
            h.veg_cm = mean(Hveg_cm, na.rm=T), # Vegetative height 
            h.reg_cm = mean(Hrep_cm, na.rm=T), # Reproductive height
            LMA = mean(`LMA_g*m-2`, na.rm=T), # Leaf mass per area g*m-2
            SLA = mean(SLA, na.rm=T))  %>% 
  unite(EP,QNR, col='id', sep = '_')

# 5.2 p.inv: 
p.inv <- read_csv(
  file.path(
    dir_data, 'SEBAS_June 2023_species inventories corrected_13.10.23.csv')) 
p.inv$QNR <- sub('^0+', '', p.inv$QNR) 
p.inv <- p.inv %>% 
  mutate(Quadrat = as.factor(Quadrat)) %>% 
  dplyr::select(c('EP','QNR','Acro','Cover')) %>% 
  spread('Acro', 'Cover', fill = NA) %>% 
  unite(EP,QNR, col='id', sep = '_')

# 5.3 p.biomass:
p.biomass <- read_csv(file.path(dir_data, 'Potsdam23.Biomass.csv')) %>% 
  dplyr::select(-c(1:4,6,8)) %>% 
  unite(EP,QNR, col='id', sep = '_') %>% 
  rename(
    RawBioM      = 'Raw_Biomass_g_0.36m²',
    BioMass_g    = 'Sophia Nicole Meyer: weight without CrispSac (1 bag 6.782g)  Biomass_g_0.36m²',
    BioMass_g_m2 = 'Sophia Nicole Meyer: column K divided by 0.36  Biomass_g_m2')

# 5.4 p.cover:
p.cover <- read_csv(file.path(dir_data, 'Potsdam23.PlotPlantCover.csv')) %>% 
  dplyr::select(-c(1:3,5,7,'Observers','Remarks')) %>% 
  unite(EP,QNR, col='id', sep = '_')
 

# Data: Physicochemical --------------------------------------------------------
# 6.1 Land use indices 2018 
#  homogenized in 20m resolution raster before processing
LUI <- read_csv(
  file.path(dir_data, 'LUI_default components set__2024-03-08.txt')) %>%
  mutate(Date = as.Date(Year, format = '%m/%d/%Y'), .before =1) %>% 
  filter(Date >= as.Date('2021-01-01')) %>% 
  dplyr::select(-c(2,3,5,6)) %>% 
  rename(EP_Plotid = `EP_PlotID`) %>% 
  group_by(EP_Plotid) %>% 
  summarise(Grazing = mean(TotalGrazing), Mowing = mean(TotalMowing), 
            Fertilization = mean(TotalFertilization))

# 6.2 Topographical Parameters; 2023
#  Elevation:	height above sea level; average over plot area
#  Slope d  :	slope; average over plot area
#  Slope p  :	slope in percent
#  Aspect   :	aspect; circular average over plot area; 360=0=North
topography <- read_csv(
  file.path(dir_data, 'EP.TopographicalParameters.csv')) %>% 
  dplyr::select(-c(1,3,4,6,8,11)) %>% 
  mutate(`EP ID` = as.factor(`EP ID`)) %>% 
  rename(EP_Plotid = `EP ID`,
         Slope_d = 'Slope d',
         Slope_p = 'Slope p')
length(unique(topography$EP_Plotid)) == length(topography$EP_Plotid)

# 6.3 Soil respiration; 2019
#  Rs: Soil respiration (g/(m^2*d))
#  Ts: Soil temperature 10 cm depth below the surface of the mineral soil	(deg C)
#  Ms: Soil volumetric water content 10 cm below the surface of the mineral soil	(%)
respiration <- read_csv(
  file.path(dir_data, 'MinSoil19.SoilRespirationFnG.csv')) %>% 
  dplyr::select(-c(2:5,7,8)) %>% 
  mutate(EP_Plotid = as.factor(EP_Plotid))
r1 <- respiration %>% 
  filter(Year == 2018)
r2 <- respiration %>% 
  filter(Year == 2019)
respira <- full_join(r1, r2, by = join_by(EP_Plotid), suffix=c('18','19')) %>% 
  dplyr::select(-c('Year18', 'Year19')) 
rm(r1,r2,respiration)
length(unique(respira$EP_Plotid)) == length(respira$EP_Plotid)

# 6.4 Soil nutrient leaching; 2018
#  Annual Leaching (kg/ha)
leach <- read_csv(
  file.path(dir_data, 'MinSoil18.SoilNutrientLeachingFnG.csv')) %>%
  dplyr::select(-c(2:8))
colnames(leach) <- paste('L', colnames(leach), sep = '.') 
leach <- leach %>% 
  rename('EP_Plotid' = 'L.EP_Plotid') %>% 
  mutate(EP_Plotid = as.factor(EP_Plotid))
length(unique(leach$EP_Plotid)) == length(leach$EP_Plotid)

# 6.5 sodium bicarbonate (NaHCO3)-extractable phosphorus (P), 
#  commonly termed Olsen-P; 2021; in mg/kg
P <- read_csv(file.path(dir_data, 'SSC21.EP0-10cm.OlsenP.csv')) %>% 
  mutate(EP_Plotid = as.factor(EP_Plotid)) %>% 
  dplyr::select(-c(2:4)) %>% 
  rename(OP = 'Olsen-P')
length(unique(P$EP_Plotid)) == length(P$EP_Plotid)

# 6.6 soil analysis
soil <- read_csv(
  file.path(dir_data, 'SeBAS_SoilAnalyses_CampaignJune2023_240403.csv'),
  col_types = cols(Time_sampled = col_character(),
                   Time_frozen  = col_character())) %>% 
  rename('EP_Plotid' = 'Plot') %>% 
  mutate(EP_Plotid = gsub('^SEG0(\\d)', 'SEG\\1', EP_Plotid),
         EP_Plotid = gsub('^HEG0(\\d)', 'HEG\\1', EP_Plotid),
         id = paste(EP_Plotid, Quadrat, sep = '_')) %>% 
  dplyr::select(-c(1:8))

# 6.7 GPS 
coordd <- read_csv(file.path(dir_data, 'coordinates_decimal_degrees.csv')) %>% 
  mutate(id = meta$id)


# Merging data -----------------------------------------------------------------
meta <- meta %>% 
  left_join(coordd, by = join_by(id)) %>% 
  left_join(topography, by = join_by(EP_Plotid)) %>% 
  left_join(soil, by = join_by(id)) %>%
  left_join(P, by = join_by(EP_Plotid)) %>%
  left_join(leach, by = join_by(EP_Plotid)) %>% 
  left_join(respira, by = join_by(EP_Plotid)) %>% 
  left_join(p.biomass, by = join_by(id)) %>% 
  left_join(p.traits, by = join_by(id)) %>% 
  left_join(p.cover, by = join_by(id)) %>% 
  left_join(p.inv, by = join_by(id)) %>% 
  left_join(LUI, by = join_by(EP_Plotid))

meta_clean <- meta %>% 
  dplyr::select(id, loc, EP_Plotid, QNR) %>% 
  left_join(coordd, by = join_by(id)) %>% 
  left_join(soil, by = join_by(id)) %>% 
  left_join(p.biomass, by = join_by(id)) %>% 
  left_join(p.traits, by = join_by(id)) %>% 
  left_join(p.cover, by = join_by(id)) %>% 
  left_join(p.inv, by = join_by(id))

rm(coordd,P,leach,respira,LUI,soil,topography, 
   p.biomass, p.traits, p.cover, p.inv)
  
write.csv(meta_clean, row.names = F,
          file.path(dir_data, 'df_meta.csv'))


# Unclassified taxa ------------------------------------------------------------
taxa$Phylum[is.na(taxa$Phylum)]<-taxa$Kingdom[is.na(taxa$Phylum)]
taxa$Class[is.na(taxa$Class)]<-taxa$Phylum[is.na(taxa$Class)]
taxa$Order[is.na(taxa$Order)]<-taxa$Class[is.na(taxa$Order)]
taxa$Family[is.na(taxa$Family)]<-taxa$Order[is.na(taxa$Family)]
taxa$Genus[is.na(taxa$Genus)]<-taxa$Family[is.na(taxa$Genus)]
taxa[] <- lapply(taxa, as.factor)
taxa=as.matrix(taxa)

# remove the sequences from the taxon names
colnames(seqtab.nochim)=paste0('ASV',1:ncol(seqtab.nochim))

# make all taxon names match
names(seqs)=colnames(seqtab.nochim)
row.names(taxa)=colnames(seqtab.nochim)
row.names(meta)=row.names(seqtab.nochim)


# Phyloseq object --------------------------------------------------------------
physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F),
                   sample_data(meta),
                   tax_table(taxa), seqs)

saveRDS(physeq,   file.path(dir_data, 'phyloseq.rds'))
physeq <- readRDS(file.path(dir_data, 'phyloseq.rds'))
