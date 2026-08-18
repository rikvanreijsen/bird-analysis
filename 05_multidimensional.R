##### PART I: CLASS-WIDE ANALYSES (ALL NON-PASSERINES) #####
##### SECTION 1: GEOGRAPHIC RANGE FOOTPRINT & LATITUDINAL SPAN #####
# Analysis: Phylogenetic Logistic Regression (MPLE)
# Question: Do species with small ranges or narrow spans hoard NISC?

df_s1           <- na.omit(pruned_phylolm_data[, c("species_tree", "nisc", "log_range", "lat_span", "abs_lat")])
rownames(df_s1) <- df_s1$species_tree
tree_s1         <- keep.tip(pruned_phylolm_tree, df_s1$species_tree)

df_s1$scale_log_range <- scale(df_s1$log_range)
df_s1$scale_lat_span  <- scale(df_s1$lat_span)
df_s1$scale_abs_lat   <- scale(df_s1$abs_lat)

fit_s1_scaled <- phyloglm(nisc ~ scale_log_range + scale_lat_span + scale_abs_lat, data = df_s1, phy = tree_s1,
                          method = "logistic_MPLE", btol = 30)

##### SECTION 2: AERODYNAMIC BIOMECHANICS & DISPERSAL CAPABILITY #####
# Analysis: Phylogenetic Logistic Regression (MPLE)
# Question: Does flight efficiency (HWI) or body size constrain NISC?

df_s2           <- na.omit(pruned_phylolm_data[, c("species_tree", "nisc", "hand_wing_index", "log_mass")])
rownames(df_s2) <- df_s2$species_tree
tree_s2         <- keep.tip(pruned_phylolm_tree, df_s2$species_tree)

fit_s2 <- phyloglm(nisc ~ hand_wing_index + log_mass, data = df_s2, phy = tree_s2, 
                          method = "logistic_MPLE", btol = 30)

##### SECTION 3: DIETARY DECOUPLING (BLUE VS. CAROTENOID GREEN) #####
# Analysis: Comparative Phylogenetic Logistic Regressions
# Question: Is pure blue diet-independent while green relies on diet?

# Define carotenoid-rich guilds (plant/fruit/seed/nectar based)
carotenoid_guilds <- c("Frugivore", "Herbivore", "Granivore", "Nectarivore")

# Create binary predictor: 1 = carotenoid-rich, 0 = non-carotenoid
pruned_phylolm_data$carotenoid_diet <- ifelse(pruned_phylolm_data$trophic_niche %in% carotenoid_guilds, 1, 0)

# Subset complete cases
df_s3           <- na.omit(pruned_phylolm_data[, c("species_tree", "b_nisc", "g_nisc", "carotenoid_diet")])
rownames(df_s3) <- df_s3$species_tree
tree_s3         <- keep.tip(pruned_phylolm_tree, df_s3$species_tree)

# Run models with increased btol
# Model 3A: Structural Blue vs Carotenoid Diet (Hypothesis: Beta = 0, no effect)
fit_s3_blue  <- phyloglm(b_nisc ~ carotenoid_diet, data = df_s3, phy = tree_s3, 
                         method = "logistic_MPLE", btol = 30)

# Model 3B: Structural Green vs Carotenoid Diet (Hypothesis: Beta > 0, strong positive effect)
fit_s3_green <- phyloglm(g_nisc ~ carotenoid_diet, data = df_s3, phy = tree_s3, 
                         method = "logistic_MPLE", btol = 30)

rm(carotenoid_guilds)

##### SECTION 4: MICROHABITAT DENSITY & CANOPY LIGHT STRATIFICATION #####
# Analysis: Phylogenetic Logistic Regression (MPLE)
# Question: Does forest canopy closure (habitat_density) favor NISC?

df_s4 <- na.omit(pruned_phylolm_data[, c("species_tree", "nisc", "habitat_density", "primary_lifestyle")])

rownames(df_s4) <- df_s4$species_tree
tree_s4         <- keep.tip(pruned_phylolm_tree, df_s4$species_tree)

fit_s4 <- phyloglm(nisc ~ habitat_density + primary_lifestyle, data = df_s4, phy = tree_s4, 
                   method = "logistic_MPLE", btol = 10)

##### SECTION 5: MULTIDIMENSIONAL ECOMORPHOSPACE (PHYLOGENETIC PCA) #####
# Analysis: Phylogenetic PCA (BM Correlation Matrix)
# Question: Do NISC species occupy a specialized region of morphospace?

morph_cols <- c("beak_length_culmen", "beak_width", "beak_depth", "tarsus_length", "wing_length",
                "tail_length", "mass")

df_s5           <- na.omit(pruned_phylolm_data[, c("species_tree", "nisc", morph_cols)])
rownames(df_s5) <- df_s5$species_tree
tree_s5         <- keep.tip(pruned_phylolm_tree, df_s5$species_tree)

morph_matrix <- log(df_s5[, morph_cols])
ppca_result  <- phyl.pca(tree_s5, morph_matrix, method = "BM", mode = "corr")
print(ppca_result)

# Test PC1 (Overall Size) and PC2 (Shape/Proportions) against NISC
df_s5$PC1_size  <- ppca_result$S[, 1]
df_s5$PC2_shape <- ppca_result$S[, 2]

fit_s5 <- phyloglm(nisc ~ PC1_size + PC2_shape, data = df_s5, phy = tree_s5, method = "logistic_MPLE", btol = 30)

print("=== SECTION 1: RANGE FOOTPRINT & GEOGRAPHY ==="); print(summary(fit_s1_scaled))
print("=== SECTION 2: AERODYNAMIC BIOMECHANICS ==="); print(summary(fit_s2))
print("=== SECTION 3A: DIETARY DECOUPLING (BLUE NISC) ==="); print(summary(fit_s3_blue))
print("=== SECTION 3B: DIETARY DECOUPLING (GREEN NISC) ==="); print(summary(fit_s3_green))
print("=== SECTION 4: MICROHABITAT STRATIFICATION ==="); print(summary(fit_s4))
print("=== SECTION 5: ECOMORPHOSPACE (PHYLO-PCA) ==="); print(summary(fit_s5))

##### PART II: MARINE-SPECIFIC ANALYSES (SEABIRDS / HABITAT == 5) #####
# Data Preparation
pruned_marine_data <- pruned_phylolm_data %>% filter(habitat == 5)

# Download Bio-ORACLE Rasters
options(sdmpredictors_datadir = tempdir())
layer_codes <- c("BO2_tempmean_ss", "BO2_salinitymean_ss", "BO22_parmean", "BO2_ppmean_ss", "BO2_chlomean_ss")

bio_oracle_raw <- load_layers(layer_codes, rasterstack = TRUE)
marine_layers  <- rast(bio_oracle_raw)
names(marine_layers) <- c("SST", "PSU", "PAR", "PROD", "CHLA")

# Extract environmental raster values at range centroids (long_x, lat_y)
marine_coords  <- pruned_marine_data[, c("long_x", "lat_y")]
marine_points  <- vect(as.matrix(marine_coords), crs = "+proj=longlat +datum=WGS84")
extracted_vals <- extract(marine_layers, marine_points)

pruned_marine_data$SST  <- extracted_vals$SST
pruned_marine_data$PSU  <- extracted_vals$PSU
pruned_marine_data$PAR  <- extracted_vals$PAR
pruned_marine_data$PROD <- extracted_vals$PROD
pruned_marine_data$CHLA <- extracted_vals$CHLA

# Clean dataset & prune tree for marine complete cases
marine_clean_vars <- c("species_tree", "nisc", "b_nisc", "g_nisc", "PAR", "PSU", "SST", "PROD", "CHLA",
                       "hand_wing_index", "log_mass", "log_range")

pruned_marine_data <- na.omit(pruned_marine_data[, marine_clean_vars])
rownames(pruned_marine_data) <- pruned_marine_data$species_tree
marine_tree <- keep.tip(pruned_phylolm_tree, pruned_marine_data$species_tree)

rm(layer_codes, bio_oracle_raw, marine_layers, marine_coords, marine_points, extracted_vals)

##### SECTION 6: OCEANOGRAPHIC PHYSICAL OPTICS VS. TROPHIC PATHWAYS #####
# Analysis: AIC Model Comparison of Abiotic (PAR, SST) vs Biotic (PROD, CHLA)
pruned_marine_data$scale_PAR      <- as.numeric(scale(pruned_marine_data$PAR))
pruned_marine_data$scale_PSU      <- as.numeric(scale(pruned_marine_data$PSU))
pruned_marine_data$scale_SST      <- as.numeric(scale(pruned_marine_data$SST))
pruned_marine_data$scale_log_PROD <- as.numeric(scale(log(pruned_marine_data$PROD + 0.001)))
pruned_marine_data$scale_log_CHLA <- as.numeric(scale(log(pruned_marine_data$CHLA + 0.001)))

rownames(pruned_marine_data) <- pruned_marine_data$species_tree

# Abiotic Model
fit_s6_abiotic <- phyloglm(nisc ~ scale_PAR + scale_PSU + scale_SST, data = pruned_marine_data, 
                           phy = marine_tree, method = "logistic_MPLE", btol = 30)

# Biotic Model
fit_s6_biotic  <- phyloglm(nisc ~ scale_log_PROD + scale_log_CHLA, data = pruned_marine_data, 
                           phy = marine_tree, method = "logistic_MPLE", btol = 30)

# Single Benchmarks
fit_s6_par     <- phyloglm(nisc ~ scale_PAR, data = pruned_marine_data, phy = marine_tree, 
                           method = "logistic_MPLE", btol = 30)
fit_s6_sst     <- phyloglm(nisc ~ scale_SST, data = pruned_marine_data, phy = marine_tree, 
                           method = "logistic_MPLE", btol = 30)
fit_s6_prod    <- phyloglm(nisc ~ scale_log_PROD, data = pruned_marine_data, phy = marine_tree, 
                           method = "logistic_MPLE", btol = 30)

##### SECTION 7: MARINE DISPERSAL & OCEANIC LIGHTSCAPE INTERACTION #####
# Analysis: Multi-predictor Model combining HWI and PAR
# Question: Do open-ocean seabirds combine high flight mobility with low-PAR tuning?

fit_s7 <- phyloglm(nisc ~ scale(hand_wing_index) + scale_PAR + scale(log_mass), 
                   data = pruned_marine_data, phy = marine_tree, method = "logistic_MPLE", btol = 30)

##### SECTION 8: TERRESTRIAL SOLAR SENSITIVITY CONTROL (ALL TERRESTRIAL BIRDS) #####
# Analysis: Terrestrial Solar Radiation Benchmark against Marine PAR
# Question: Does terrestrial solar radiation (SOLAR) NOT drive NISC on land, 
#           confirming PAR constraints are unique to marine biomes?

pruned_terrestrial_data <- pruned_phylolm_data %>% 
  filter(habitat != 5) %>%
  filter(!is.na(long_x) & !is.na(lat_y))

solar_monthly     <- worldclim_global(var = "srad", res = 5, path = tempdir())
global_solar_grid <- mean(solar_monthly)
names(global_solar_grid) <- "SOLAR"

terr_points <- vect(
  as.matrix(pruned_terrestrial_data[, c("long_x", "lat_y")]), 
  crs = "EPSG:4326"
)

extracted_vals <- extract(global_solar_grid, terr_points)
pruned_terrestrial_data$SOLAR <- extracted_vals$SOLAR

df_s8 <- pruned_terrestrial_data %>% 
  filter(!is.na(nisc) & !is.na(SOLAR)) %>% 
  select(species_tree, nisc, SOLAR)

df_s8$scale_SOLAR <- as.numeric(scale(df_s8$SOLAR))

rownames(df_s8) <- df_s8$species_tree
tree_s8         <- keep.tip(pruned_phylolm_tree, df_s8$species_tree)

fit_s8_scaled   <- phyloglm(nisc ~ scale_SOLAR, data = df_s8, phy = tree_s8, 
                            method = "logistic_MPLE", btol = 30)

rm(solar_monthly, global_solar_grid, terr_points, extracted_vals)

##### SECTION 9: SALINITY-TEMPERATURE GRADIENTS #####
# Analysis: Bivariate/Multi-predictor Model evaluating SST and Salinity (PSU)
# Question: Do ocean thermal and osmotic conditions structure NISC occurrence in marine lineages?

fit_s9 <- phyloglm(nisc ~ scale_SST + scale_PSU, data = pruned_marine_data, phy = marine_tree, 
                   method = "logistic_MPLE", btol = 30)

##### SECTION 10: SYNTHETIC MARINE MULTIVARIATE MODEL #####
# Analysis: Full Synthetic Marine Model (PAR + SST + HWI + Range Size + Body Mass)
# Question: What are the relative contributions of environmental, biomechanical, and macroecological factors?

# Subset complete cases across all variables used in the synthetic model
synthetic_vars <- c("nisc", "scale_PAR", "scale_SST", "hand_wing_index", "log_range", "log_mass")
complete_marine_data <- pruned_marine_data[complete.cases(pruned_marine_data[, synthetic_vars]), ]

# Ensure phylogenetic tree matches complete cases
synthetic_marine_tree <- keep.tip(marine_tree, complete_marine_data$species_tree)

fit_s10 <- phyloglm(nisc ~ scale_PAR + scale_SST + scale(hand_wing_index) + scale(log_range) + scale(log_mass), 
                    data = complete_marine_data, phy = synthetic_marine_tree, method = "logistic_MPLE", 
                    btol = 30)

rm(synthetic_vars, complete_marine_data)

print("=== SECTION 6A: MULTIVARIATE ABIOTIC (PAR + PSU + SST) ==="); print(summary(fit_s6_abiotic))
print("=== SECTION 6B: MULTIVARIATE BIOTIC (PROD + CHLA) ==="); print(summary(fit_s6_biotic))
print("=== SECTION 6C: SINGLE BENCHMARK (PAR) ==="); print(summary(fit_s6_par))
print("=== SECTION 6D: SINGLE BENCHMARK (SST) ==="); print(summary(fit_s6_sst))
print("=== SECTION 6E: SINGLE BENCHMARK (PROD) ==="); print(summary(fit_s6_prod))
print("=== SECTION 7: MARINE DISPERSAL & LIGHTSCAPE (HWI + PAR + MASS) ==="); print(summary(fit_s7))
print("=== SECTION 8: TERRESTRIAL SOLAR CONTROL (SOLAR) ==="); print(summary(fit_s8_scaled))
print("=== SECTION 9: SALINITY-TEMPERATURE GRADIENTS (SST + PSU) ==="); print(summary(fit_s9))
print("=== SECTION 10: SYNTHETIC MARINE MULTIVARIATE MODEL ==="); print(summary(fit_s10))
