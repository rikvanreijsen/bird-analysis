##### PHYLOGENETIC LOGISTIC REGRESSION ANALYSIS #####
# Programs
library(phangorn); library(plotrix)

# Data Preparation
phylo_log_data <- meta_data %>%
  select(species_tree, fam_species, order_species, sc, b_nisc, g_nisc, p_nisc, nisc, min_lat, max_lat,
         long_x, lat_y, range_size, habitat, habitat_density, migration)

# Prune Data & Tree
rownames(phylo_log_data) <- phylo_log_data[,1]
row.names(phylo_log_data)
names(phylo_log_data)

pruned_phylo_log_tree <- treedata(bird_tree, phylo_log_data, sort = T)$phy
name.check(bird_tree, phylo_log_data)

pruned_phylo_log_data <- treedata(bird_tree, phylo_log_data, sort = T)$data
pruned_phylo_log_data <- data.frame(pruned_phylo_log_data)
name.check(pruned_phylo_log_tree, pruned_phylo_log_data)

lat_y_naomit <- pruned_phylo_log_data[complete.cases(pruned_phylo_log_data[, 11]), ]
lat_y_naomit$lat_y <- as.numeric(lat_y_naomit$lat_y)
lat_y_N <- subset(lat_y_naomit, lat_y > "0")
lat_y_S <- subset(lat_y_naomit, lat_y < "0")
lat_y_S$lat_y <- as.numeric(abs(lat_y_S$lat_y))

habitat_naomit <- pruned_phylo_log_data[complete.cases(pruned_phylo_log_data[, 14]), ]
habitat_naomit$habitat <- as.numeric(habitat_naomit$habitat)
habitat <- habitat_naomit

habitat_density_naomit <- pruned_phylo_log_data[complete.cases(pruned_phylo_log_data[, 15]), ]
habitat_density_naomit$habitat <- as.numeric(habitat_density_naomit$habitat)
habitat_density <- habitat_density_naomit

range_size_naomit <- pruned_phylo_log_data[complete.cases(pruned_phylo_log_data[, 13]), ]
range_size_naomit$range_size <- as.numeric(range_size_naomit$range_size)
range_size <- range_size_naomit

migration_naomit <- pruned_phylo_log_data[complete.cases(pruned_phylo_log_data[, 16]), ]
migration_naomit$migration <- as.numeric(migration_naomit$migration)
migration <- migration_naomit

rm(lat_y_naomit, range_size_naomit, habitat_naomit, habitat_density_naomit, migration_naomit)

##### Map Traits on Tree #####
# North Latitude
lat_y_N <- data.frame(as.list(lat_y_N))
rownames(lat_y_N) <- lat_y_N[, 1] 
row.names(lat_y_N)
names(lat_y_N)

bird_tree_lat_N <- treedata(pruned_phylo_log_tree, lat_y_N, sort = T, warnings = T)$phy
name.check(bird_tree_lat_N, lat_y_N)

# South Latitude
lat_y_S <- data.frame(as.list(lat_y_S))
rownames(lat_y_S) <- lat_y_S[, 1]
row.names(lat_y_S)
names(lat_y_S)

bird_tree_lat_S <- treedata(pruned_phylo_log_tree, lat_y_S, sort = T, warnings = T)$phy
name.check(bird_tree_lat_S, lat_y_S)

# Habitat
habitat <- data.frame(as.list(habitat))
rownames(habitat) <- habitat[, 1]
row.names(habitat)
names(habitat)

bird_tree_hab <- treedata(pruned_phylo_log_tree, habitat, sort = T, warnings = T)$phy
name.check(bird_tree_hab, habitat)

# Habitat Density
habitat_density <- data.frame(as.list(habitat_density))
rownames(habitat_density) <- habitat_density[, 1]
row.names(habitat_density)
names(habitat_density)

bird_tree_hab_dens <- treedata(pruned_phylo_log_tree, habitat_density, sort = T, warnings = T)$phy
name.check(bird_tree_hab_dens, habitat_density)

# Range Size
range_size <- data.frame(as.list(range_size))
rownames(range_size) <- range_size[, 1]
row.names(range_size)
names(range_size)

bird_tree_rng_sz <- treedata(pruned_phylo_log_tree, range_size, sort = T, warnings = T)$phy
name.check(bird_tree_rng_sz, range_size)

# Migration
migration <- data.frame(as.list(migration))
rownames(migration) <- migration[, 1]
row.names(migration)
names(migration)

bird_tree_mig <- treedata(pruned_phylo_log_tree, migration, sort = T, warnings = T)$phy
name.check(bird_tree_mig, migration)

##### Map Traits on Phylogeny #####
lat_N_trait    <- setNames(lat_y_N$lat_y, lat_y_N$species_tree)
lat_S_trait    <- setNames(lat_y_S$lat_y, lat_y_S$species_tree)
hab_trait      <- as.factor(setNames(habitat$habitat, habitat$species_tree))
hab_dens_trait <- as.factor(setNames(habitat_density$habitat_density, habitat_density$species_tree))
rng_sz_trait   <- setNames(range_size$range_size, range_size$species_tree)
mig_trait      <- setNames(migration$migration, migration$species_tree)

##### Data Exploration #####
# Phylogenetic Signal for Binary Variable Purvis' D (Set Zero-Length Branches to Be 1/1000000 Total Tree Length)
require(caper); require(phylolm)

# North Latitude
no_zero_tree_lat_N <- bird_tree_lat_N
no_zero_tree_lat_N$edge.length[no_zero_tree_lat_N$edge.length == 0] <- max(nodeHeights(bird_tree_lat_N))*1e-6
no_zero_tree_lat_N$edge.length

col_lat_N.D <- phylo.d(lat_y_N, no_zero_tree_lat_N, names.col = species_tree, binvar = nisc, permut = 1000,
                       rnd.bias = NULL)
print(col_lat_N.D)

# South Latitude
no_zero_tree_lat_S <- bird_tree_lat_S
no_zero_tree_lat_S$edge.length[no_zero_tree_lat_S$edge.length == 0] <- max(nodeHeights(bird_tree_lat_S))*1e-6
no_zero_tree_lat_S$edge.length

col_lat_S.D <- phylo.d(lat_y_S, no_zero_tree_lat_S, names.col = species_tree, binvar = nisc, permut = 1000,
                       rnd.bias = NULL)
print(col_lat_S.D)

# Habitat
no_zero_tree_hab <- bird_tree_hab
no_zero_tree_hab$edge.length[no_zero_tree_hab$edge.length == 0] <- max(nodeHeights(bird_tree_hab))*1e-6
no_zero_tree_hab$edge.length

col_hab.D <- phylo.d(habitat, no_zero_tree_hab, names.col = species_tree, 
                     binvar = nisc, permut = 1000, rnd.bias = NULL)
print(col_hab.D)

# Habitat Density
no_zero_tree_hab_dens <- bird_tree_hab_dens
no_zero_tree_hab_dens$edge.length[no_zero_tree_hab_dens$edge.length == 0] <- 
  max(nodeHeights(bird_tree_hab_dens))*1e-6
no_zero_tree_hab_dens$edge.length

col_hab_dens.D <- phylo.d(habitat_density, no_zero_tree_hab_dens, names.col = species_tree, binvar = nisc, 
                          permut = 1000, rnd.bias = NULL)
print(col_hab_dens.D)

# Range Size
no_zero_tree_rng_sz <- bird_tree_rng_sz
no_zero_tree_rng_sz$edge.length[no_zero_tree_rng_sz$edge.length == 0] <- max(nodeHeights(bird_tree_rng_sz))*1e-6
no_zero_tree_rng_sz$edge.length

col_rng_sz.D <- phylo.d(range_size, no_zero_tree_rng_sz, names.col = species_tree, binvar = nisc, permut = 1000,
                        rnd.bias = NULL)
print(col_rng_sz.D)

# Migration
no_zero_tree_mig <- bird_tree_mig
no_zero_tree_mig$edge.length[no_zero_tree_mig$edge.length == 0] <- max(nodeHeights(bird_tree_mig))*1e-6
no_zero_tree_mig$edge.length

col_mig.D <- phylo.d(migration, no_zero_tree_mig, names.col = species_tree, binvar = nisc, permut = 1000,
                     rnd.bias = NULL)
print(col_mig.D)

rm(bird_tree_lat_N, bird_tree_lat_S, bird_tree_hab, bird_tree_hab_dens, bird_tree_rng_sz, bird_tree_mig)

##### Phylogenetic Logistic Regression Model #####
# North Latitude
model_lat_N = phyloglm(nisc ~ lat_N_trait, phy = no_zero_tree_lat_N, data = lat_y_N, 
                       method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_lat_N)

# South Latitude
model_lat_S = phyloglm(nisc ~ lat_S_trait, phy = no_zero_tree_lat_S, data = lat_y_S,
                       method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_lat_S)

# Habitat
model_hab = phyloglm(nisc ~ -1 + hab_trait, phy = no_zero_tree_hab, data = habitat,
                     method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_hab)

# Habitat Density
model_hab_dens = phyloglm(nisc ~ -1 + hab_dens_trait, phy = no_zero_tree_hab_dens, data = habitat_density,
                     method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_hab_dens)

# Range Size
model_rng_sz = phyloglm(nisc ~ rng_sz_trait, phy = no_zero_tree_rng_sz, data = range_size,
                     method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_rng_sz)

# Migration
model_mig = phyloglm(nisc ~ mig_trait, phy = no_zero_tree_mig, data = migration,
                        method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_mig)

rm(no_zero_tree_lat_N, no_zero_tree_lat_S, no_zero_tree_hab, no_zero_tree_hab_dens, no_zero_tree_rng_sz,
   no_zero_tree_mig, lat_N_trait, lat_S_trait, hab_trait, hab_dens_trait, rng_sz_trait, mig_trait)
