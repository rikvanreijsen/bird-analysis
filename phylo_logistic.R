# Reset Commands
closeAllConnections()
rm(list = ls())

# Programs
require(ape); require(geiger); require(ggplot2); require(phangorn); require(phytools); require(plotrix)

##### DATA PREPARATION #####
setwd("~/Desktop/R/Data")
bird_tree <- read.tree("Stage2_Hackett_MCC_no_neg.tree")
is.ultrametric(bird_tree)
bird_tree
data_overall <- read.csv2("~/Desktop/R/Data/DataRAnalysis.csv")
data_overall <- dataframe(as.list(data_overall))

# Prune Data & Tree
rownames(data_overall) <- data_overall[,1]
row.names(data_overall)
names(data_overall)

pruned_tree_overall <- treedata(bird_tree, data_overall, sort = T)$phy
name.check(bird_tree, data_overall)

pruned_data_overall <- treedata(bird_tree, data_overall, sort = T)$data
pruned_data_overall <- data.frame(pruned_data_overall)
name.check(pruned_tree_overall, pruned_data_overall)

lat_y_naomit <- pruned_data_overall[complete.cases(pruned_data_overall[,7]),]
lat_y_naomit$lat_y <- as.numeric(lat_y_naomit$lat_y)
lat_y_N <- subset(lat_y_naomit, lat_y > "0")
lat_y_S <- subset(lat_y_naomit, lat_y < "0")
lat_y_S$lat_y <- as.numeric(abs(lat_y_S$lat_y))
habitat_naomit <- pruned_data_overall[complete.cases(pruned_data_overall[,5]),]
habitat_naomit$habitat <- as.numeric(habitat_naomit$habitat)

# Map Traits on Tree
# North Latitude
lat_y_N <- data.frame(as.list(lat_y_N))
rownames(lat_y_N) <- lat_y_N[,1]
row.names(lat_y_N)
names(lat_y_N)

bird_tree_lat_N <- treedata(pruned_tree_overall, lat_y_N, sort = T, warnings = T)$phy
name.check(bird_tree_lat_N, lat_y_N)

# South Latitude
lat_y_S <- data.frame(as.list(lat_y_S))
rownames(lat_y_S) <- lat_y_S[,1]
row.names(lat_y_S)
names(lat_y_S)

bird_tree_lat_S <- treedata(pruned_tree_overall, lat_y_S, sort = T, warnings = T)$phy
name.check(bird_tree_lat_S, lat_y_S)

# Habitat
habitat_naomit <- data.frame(as.list(habitat_naomit))
rownames(habitat_naomit) <- habitat_naomit[,1]
row.names(habitat_naomit)
names(habitat_naomit)

bird_tree_hab <- treedata(pruned_-tree_overall, habitat_naomit, sort = T, warnings = T)$phy
name.check(bird_tree_hab, habitat_naomit)

# Map Traits on Phylogeny
lat_N_trait <- setNames(lat_y_N$lat_y, lat_y_N$species_tree)
lat_S_trait <- setNames(Lat_y_S$lat_y, lat_y_S$species_tree)
hab_trait <- as.factor(setNames(habitat_naomit$habitat, habitat_naomit$species_tree))

##### DATA EXPLORATION #####
# Phylogenetic Signal for Binary Variable Purvis' D
require(caper); require(phylolm)

col_lat_N.D <- phylo.d(lat_y_N, bird_tree_lat_N, names.col = species_tree, binvar = nisc, permut = 1000, rnd.bias = NULL)
col_lat_S.D <- phylo.d(lat_y_S, bird_tree_lat_S, names.col = species_tree, binvar = nisc, permut = 1000, rnd.bias = NULL)
col_hab.D <- phylo.d(habitat_naomit, bird_tree_hab, names.col = species_tree, binvar = nisc, permut = 1000, rnd.bias = NULL)

# Set Zero-Length Branches to Be 1/1000000 Total Tree Length
# North Latitude
no_zero_tree_lat_N <- bird_tree_lat_N
no_zero_tree_lat_N$edge.length[no_zero_tree_lat_N$edge.length == 0] <- max(nodeHeights(bird_tree_lat_N))*1e-6
no_zero_tree_lat_N$edge.length

col_lat_N.D <- phylo.d(lat_y_N, no_zero_tree_lat_N, names.col = species_tree, binvar = nisc, permut = 1000, rnd.bias = NULL)
print(col_lat_N.D)
plot(col_lat_N.D)

# South Latitude
no_zero_tree_lat_S <- bird_tree_lat_S
no_zero_tree_lat_S$edge.length[no_zero_tree_lat_S$edge.length == 0] <- max(nodeHeights(bird_tree_lat_S))*1e-6
no_zero_tree_lat_S$edge.length

col_lat_S.D <- phylo.d(lat_y_S, no_zero_tree_lat_S, names.col = species_tree, binvar = nisc, permut = 1000, rnd.bias = NULL)
print(col_lat_S.D)
plot(col_lat_S.D)

# Habitat
no_zero_tree_hab <- bird_tree_hab
no_zero_tree_hab$edge.length[no_zero_tree_hab$edge.length == 0] <- max(nodeHeights(bird_tree_hab))*1e-6
no_zero_length_tree_hab$edge.length

col_hab.D <- phylo.d(habitat_naomit, no_zero_tree_hab, names.col = species_tree, binvar = nisc, permut = 1000, rnd.bias = NULL)
print(col_hab.D)
plot(col_hab.D)

##### Creating the Plots #####
# Phylogenetic Logistic Regression Model North Latitude
model_lat_N = phyloglm(nisc ~ lat_N_trait, phy = no_zero_tree_lat_N, data = lat_y_N, method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_lat_N)
plot(model_lat_N)

plot(jitter(lat_y_N$lat_y), lat_y_N$nisc, pch = 1, col = 'gray80', xlab = "North Latitude", ylab = "Probability of NISC")

x <- lat_y_N$lat_y
cc1 <- coef(model_lat_N)
curve(plogis(cc1[1]+cc1[2]*x), col = 'indianred2', lwd = 3, add = T)

# Phylogenetic Logistic Regression Model South Latitude
model_lat_S = phyloglm(nisc ~ lat_S_trait, phy = no_zero_tree_lat_S, data = lat_y_S, method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_lat_S)
plot(model_lat_S)

plot(jitter(lat_y_S$lat_y), lat_y_S$nisc, pch = 1, col = 'gray80', xlab = "South Latitude", ylab = "Probability of NISC")

y = lat_y_S$lat_y
cc2 <- coef(model_lat_S)
curve(plogis(cc2[1] + cc2[2]*y), col = 'indianred2', lwd = 3, add = T)

# Phylogenetic Logistic Regression Model Habitat
model_hab = phyloglm(nisc ~ -1 + hab_trait, phy = no_zero_tree_hab, data = habitat_naomit, method = c("logistic_MPLE"), log.alpha.bound = 5, btol = 30)
summary(model_hab)
plot(model_hab)

plot(jitter(habitat_naomit$habitat), habitat_naomit$nisc, pch = 1, col = 'gray80', xlab = "Habitat", ylab = "Probability of NISC")

z = habitat_naomit$habitat
cc3 <- coef(model_hab)
curve(plogis(cc3[1]+cc3[2]+cc3[3]+cc3[4]+cc3[5]+cc3[6]*z), col = 'indianred2', lwd = 3, add = T)
