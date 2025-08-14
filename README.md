# R Code for Publication: "Still, the Environment Selects: Disentangling Spatial and Environmental Effects on Soil Bacterial Communities" – SeBAS June 2023 Grassland Plots

This repository contains the **R code** used in the publication:

**Neisse et al., 2025.** *Still, the Environment Selects: Disentangling Spatial and Environmental Effects on Soil Bacterial Communities.*

The code allows full reproduction of the analyses presented in the study. To run it, you will need datasets **32156** and **32155**. Detailed setup and usage instructions are provided in the repository’s README.

---

## Project Overview

This study investigates the drivers of **soil bacterial community assembly** in managed grassland ecosystems. The main goal is to understand how **spatial distance, environmental variation, and plant communities** influence soil microbial diversity and distribution patterns, focusing specifically on the **distance-decay relationship (DDR)** of soil bacteria at local (meter) to regional (kilometer) scales.

**Sampling details:**
- Date: June 2023  
- Locations: Two Biodiversity Exploratories regions in Germany: **Hainich-Dün (Thuringia)** and **Schorfheide-Chorin (Brandenburg)**  
- Design: 18 grassland plots, each with three 1 m² subplots along a 50-meter south–north transect

**Data processing and analysis:**
- Microbial data processed with **DADA2** in R, with taxonomic assignment via **SILVA**.  
- Statistical analyses included:
  - **GLM** with Gamma distribution  
  - **Commonality analysis** to assess geographic, abiotic, plant, and vegetation influences on DDRs  
  - **Variation partitioning** and permutation tests for factor significance  
- Spatial range of ASVs quantified using **convex hull area** based on geographic distribution

---

## Data

Processed data and metadata are available in the **Biodiversity Exploratories Information System (BExIS)**:

- **Amplicon sequence variants** (ASVs) and relative abundance:  
  Neisse, N. (2025a). *Amplicon sequence variants, their taxonomic classification, and relative abundance per sample in grassland experimental plots (SeBAS June 2023, No. 32155; Version 4)* [Dataset]. [BExIS link](https://www.bexis.uni-jena.de/ddm/data/Showdata/32155)

- **Geographic, soil, and vegetation data:**  
  Neisse, N. (2025b). *Geographic location, soil physicochemical properties, plant community composition and vegetation characteristics in grassland experimental plots (SeBAS June 2023, No. 32156; Version 6)* [Dataset]. [BExIS link](https://www.bexis.uni-jena.de/ddm/data/Showdata/32156)

**Setup script:**  
Run the following to prepare your data for analysis:  

```R
source("data_bexis_from.R")
