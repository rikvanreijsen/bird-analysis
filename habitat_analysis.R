# Reset Commands
closeAllConnections()
rm(list = ls())

# Programs
require(ecospat); require(knitr); require(maps); require(maptools)
require(raster); require(rgeos); require(sp); require(spThin)

# Data Preparation
setwd("~/Desktop/R/Data")
data_overall <- read.csv2("~/Desktop/R/Data/DataRAnalysis.csv")
occ_points_overall <- subset(data_overall, select = -c(fam_species, order_species, habitat))
occ_points_overall <- occ_points_overall[complete.cases(occ_points_overall[,3]),]
data(wrld_simpl)

##### DISTRIBUTION OF SC #####
# Match Distribution of Occurrence Records with Absence and Presence of SC to World Map
occ_ab_sc <- subset(occ_points_overall, sc == "0")
occ_ab_thin_sc <- thin(occ_ab_sc, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                       spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                       write.log.file = F, locs.thinned.list.return = T)

occ_pr_sc <- subset(occ_points_overall, sc == "1")
occ_pr_thin_sc <- thin(occ_pr_sc, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                       spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                       write.log.file = F, locs.thinned.list.return = T)

plot(wrld_simpl, asp = 1)
points(occ_ab_thin_sc[[1]], col = "gray80", pch = 4, cex = 0.5)
points(occ_pr_thin_sc[[1]], col = "indianred2", pch = 3, cex = 0.5)

##### DISTRIBUTION OF NISC #####
# Match Distribution of Occurrence Records with Absence and Presence of NISC to World Map
occ_ab_nisc <- subset(occ_points_overall, nisc == "0")
occ_ab_thin_nisc <- thin(occ_ab_nisc, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                         spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                         write.log.file = F, locs.thinned.list.return = T)

occ_pr_nisc <- subset(occ_points_overall, nisc == "1")
occ_pr_thin_nisc <- thin(occ_pr_nisc, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                         spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                         write.log.file = F, locs.thinned.list.return = T)

plot(wrld_simpl, asp = 1)
points(occ_ab_thin_nisc[[1]], col = "gray80", pch = 4, cex = 0.5)
points(occ_pr_thin_nisc[[1]], col = "indianred2", pch = 3, cex = 0.5)

##### DISTRIBUTION OF BLUE NISC #####
# Match Distribution of Occurrence Records with Absence and Presence of Blue NISC to World Map
occ_ab_b <- subset(occ_points_overall, b_nisc == "0")
occ_ab_thin_b <- thin(occ_ab_b, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                      spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                      write.log.file = F, locs.thinned.list.return = T)

occ_pr_b <- subset(occ_points_overall, b_nisc == "1")
occ_pr_thin_b <- thin(occ_pr_b, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                      spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                      write.log.file = F, locs.thinned.list.return = T)

plot(wrld_simpl, asp = 1)
points(occ_ab_thin_b[[1]], col = "gray80", pch = 4, cex = 0.5)
points(occ_pr_thin_b[[1]], col = "midnightblue", pch = 3, cex = 0.5)

##### DISTRIBUTION OF GREEN NISC #####
# Match Distribution of Occurrence Records with Absence and Presence of Green NISC to World Map
occ_ab_g <- subset(occ_points_overall, g_nisc == "0")
occ_ab_thin_g <- thin(occ_ab_g, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                      spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                      write.log.file = F, locs.thinned.list.return = T)

occ_pr_g <- subset(occ_points_overall, g_nisc == "1")
occ_pr_thin_g <- thin(occ_pr_g, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                      spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                      write.log.file = F, locs.thinned.list.return = T)

plot(wrld_simpl, asp = 1)
points(occ_ab_thin_g[[1]], col = "gray80", pch = 4, cex = 0.5)
points(occ_pr_thin_g[[1]], col = "seagreen4", pch = 3, cex = 0.5)

##### DISTRIBUTION OF PURPLE NISC #####
# Match Distribution of Occurrence Records with Absence and Presence of Purple NISC to World Map
occ_ab_p <- subset(occ_points_overall, p_nisc == "0")
occ_ab_thin_p <- thin(occ_ab_p, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                      spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                      write.log.file = F, locs.thinned.list.return = T)

occ_pr_p <- subset(occ_points_overall, p_nisc == "1")
occ_pr_thin_p <- thin(occ_pr_p, verbose = F, lat.col = "lat_y", long.col = "long_x", 
                      spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                      write.log.file = F, locs.thinned.list.return = T)

plot(wrld_simpl, asp = 1)
points(occ_ab_thin_p[[1]], col = "gray80", pch = 4, cex = 0.5)
points(occ_pr_thin_p[[1]], col = "purple4", pch = 3, cex = 0.5)
