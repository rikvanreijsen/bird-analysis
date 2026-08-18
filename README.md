# Macroevolutionary Ratcheting and Oceanographic Boundaries Shape Avian Structural Colouration: Flight Biomechanics and Marine Optical Fouling in Non-Passerines

> **Short Title:** BIOPHYSICAL BOUNDARIES OF STRUCTURAL COLOUR  
> **Repository for:** Class-wide macroevolutionary, biomechanical, and spatial analysis of non-iridescent structural colouration across non-passerine birds ($N = 3,905$).

---

## 📌 Overview

This repository contains the dataset matrices, phylogenetic trees, environmental rasters, and reproducible **R** code scripts for analyzing the macroevolutionary drivers and physical constraints of non-iridescent structural colouration (NISC) in non-passerine birds.

### Key Conceptual Findings
* **Macroevolutionary Ratchet:** Trait acquisition ($q_{01}$) systematically outpaces secondary loss ($q_{10}$) across blue (20-fold), green (29-fold), and purple (168-fold) hues ($p < 0.001$), anchored by developmental canalization.
* **Dietary Decoupling:** Structural blue ($p = 0.34$) and green ($p = 0.44$) plumage evolved independently of carotenoid-rich foraging niches.
* **Biomechanical Trade-offs:** High Hand-Wing Index ($\beta = -0.028, p < 0.001$) and body mass ($\beta = -0.227, p < 0.001$) suppress NISC due to mechanical wear on porous medullary barbs under aerodynamic load.
* **Marine Optical Fouling:** Marine habitat represents a primary filter ($\beta = -3.07, p = 0.037$), where dynamic soaring ($\text{HWI } \beta = -3.09, p < 0.01$) and sea surface temperature ($\text{SST } \beta = -2.62, p < 0.05$) drive suppression as evaporative salt micro-crystals foul the keratin-air refractive index interface ($n \approx 1.54$ vs. $n = 1.00$).

---

## 📁 Repository Structure

.
├── data/
│   ├── processed/          # Cleaned species-level phenotype and trait matrices
│   └── spatial/            # Species distribution overlays (WorldClim v2.1 & Bio-ORACLE v2.0)
├── phylo/
│   └── trees/              # Consensus phylogenetic tree files (.tre / .nexus)
├── R/
│   ├── 01_macroevolution.R # Continuous-time Markov models (ARD vs. ER) & transition rates
│   ├── 02_diet_pcap.R      # Phylogenetic logistic regressions (phylolm) for dietary niches
│   ├── 03_biomechanics.R   # HWI, body mass, and morphological PCAs
│   └── 04_oceanography.R   # Marine sub-models, SST, and surface radiation benchmarking
└── README.md

---

## 💻 Requirements & Dependencies

To reproduce the statistical analyses, spatial overlays, and figure generations, you will need **R (v4.0+)**. The script pipeline relies on the following R package libraries:

```r
# Comparative Methods & Phylogenetics
library(ape)
library(geiger)
library(phytools)
library(TreeTools)
library(phangorn)
library(caper)
library(phylolm)

# Spatial Data, GIS & Macroecology
library(sf)
library(terra)
library(raster)
library(sp)
library(rgeos)
library(maptools)
library(maps)
library(ecospat)
library(spThin)
library(rnaturalearth)
library(rnaturalearthdata)
library(sdmpredictors)

# Data Wrangling, Plotting & Formatting
library(dplyr)
library(ggplot2)
library(plotrix)
library(knitr)
```

**Note on Spatial Dependencies:** Spatial operations primarily utilize modern sf and terra frameworks. Legacy spatial dependencies (raster, sp, rgeos, maptools) are retained for backward compatibility with specific spatial matrices or archived data layers.

---

## 📜 Citation & License
This project is licensed under the MIT License for code and CC-BY 4.0 for data matrices. A permanent DOI will be issued via Zenodo upon publication.
