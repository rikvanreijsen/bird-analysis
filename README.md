# Macroevolutionary & Computational Analysis of Avian Structural Phenotypes

An automated spatial data-augmentation, phylogenetic comparative, and biophysical analysis pipeline evaluating 
the macroevolutionary dynamics of Non-Iridescent Structural Colouration (NISC) in global non-passerine lineages ($N = 4,335$).

## 📌 Project Overview
This repository contains the full analytical codebase and computational pipelines developed for my undergraduate thesis at the 
University of British Columbia (UBC). The research isolates proximate biophysical mechanics and ultimate macroevolutionary patterns 
governing angle-independent avian structural hues (blues, purples, and structural greens).

[ Lineage & Trait Data ]        [ Evolutionary Tree ]
(DataRAnalysis.csv)       (Stage2_Hackett_MCC_no_neg.tree)
│                                │
├────────────────────────────────┤
▼                                ▼
┌───────────────────┐            ┌───────────────────┐
│   phylo_tree.R    │            │ phylo_logistic.R  │
│ (Pruning & Tree   │            │ (Phylogenetic     │
│   Preparation)    │            │     Regression)   │
└───────────────────┘            └───────────────────┘

[ Spatial/Look-up Data ]        [ Spectral Data ]
(points_missing.csv)           (reflecdata.csv)
│                                │
▼ (Fills gaps in)                ▼
┌───────────────────┐            ┌───────────────────┐
│ habitat_analysis.R│            │ feather_reflec.R  │
│ (Fills unknown    │            │ (Spectrophotometry│
│   data points)    │            │    Analysis)      │
└───────────────────┘            └───────────────────┘

---

## 🔬 Core Code Repositories & Execution Logic

### Track A: Phylogenetic Comparative Framework
* **`phylo_tree.R`**: Quality control script that loads a massive consensus Maximum Clade Credibility (MCC) tree built on the Hackett 
backbone, rectifies branch-length anomalies (negative/zero lengths), drops non-target lineages, and prunes the topology down to the 
target $4,335$ non-passerine species.
* **`phylo_logistic.R`**: Executes high-dimensional Phylogenetic Logistic Regressions (using a binary logit link framework). By 
correcting error structures using the pruned phylogenetic distance matrix, this script isolates true adaptive environmental responses 
from neutral phylogenetic niche conservatism.

### Track B: Macroecological Data Augmentation & Spatial Auditing
* **`habitat_analysis.R`**: Addresses severe geographic data-gap bottlenecks by capturing breeding centroids on coastal margins and 
remote island networks that fall outside default vector layers. 
  * Integrates a spatial patch matrix (`points_missing.csv`) using a continental boundary buffering framework (`gBuffer`).
  * Applies spatial thinning via `spThin` (10 km parameter constraint) to standardize global search effort.
  * Employs an automated taxonomic frequency filter to dynamically map point densities by family and order.

### Track C: Biophysical Spectral Metrics
* **`feather_reflec.R`**: Sweeps raw reflectance matrices generated from Scanning Electron Microscopy (SEM) and objective feather 
spectrophotometry across the full tetrachromatic avian visual spectrum ($300\text{–}700\text{ nm}$). Isolates continuous phenotypic 
properties (brightness, chroma, peak wavelengths) from background interference noise.

---

## 📊 Summary of Major Empirical Findings

* **Independent Macroevolutionary Scaling:** NISC evolved independently **27 times** across deep-time non-passerine lineages, 
manifesting as a highly clustered trait exhibiting powerful phylogenetic signals. Transition rate modelling dictates that lineages 
display a higher probability of gaining NISC than losing it.
* **The Marine Environmental Driver:** Phylogenetic comparative regressions revealed **no statistical correlation** between structural 
plumage traits and terrestrial biomes (forest, grassland, wetland) or latitudinal gradients. However, a **significant, non-coincidental 
correlation emerged with marine habitats**, identifying marine systems as unique evolutionary selective forces.
* **Geospatial Hotspots:** Utilizing spatial data-gap clearing corrected significant sample-size artifacts (which initially inflated 
localized prevalence to an unrealistic $100\%$ across isolated sites). Stabilized spatial outputs definitively position **Oceania 
($26.61\%$)** and **Asia ($23.11\%$)** as global hotspots, with **South-Eastern Asia ($29.48\%$)**, **Melanesia ($29.17\%$)**, 
and **Polynesia ($29.03\%$)** acting as evolutionary epicentres.
* **Ultraviolet Convergence & Discovery:** Spectrophotometric arrays confirmed a powerful positive correlation between colour intensity 
and ultraviolet (UV) reflectance. Furthermore, microstructural analyses led to the discovery of characteristic NISC configurations within 
the family **Ardeidae** (herons and egrets)—a structural mechanism previously documented almost exclusively within **Spheniscidae** 
(penguins).

---

## 🛠️ Data-Visualization Dashboards

### Global Taxonomic Distributions (By Order)
The script automates frequency filtering to isolate the top ten most abundant non-passerine orders, compressing rarer lineages into a 
uniform background baseline (`"Other"`) to prevent categorical inflation and preserve spatial context.

### Global Ecological Niche Distributions (By Habitat)
Geospatial projection of thinned coordinates mapped directly onto categorical environmental niches (where numbers $1\text{–}6$ represent 
distinct forest, grassland, wetland, marine, and generalist matrices). This map visually isolates the distinct clustering of NISC 
phenotypes along coastal fringes, island chains, and high-latitude Antarctic corridors.

---

## 🧰 Dependencies & Setup
To reproduce the analytical pipelines, ensure the following R environment libraries are initialized:
```R
install.packages(c("ecospat", "knitr", "maps", "maptools", "raster", "rgeos", "sp", "spThin", "viridis"))
