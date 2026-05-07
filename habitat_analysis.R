# Reset Commands
closeAllConnections()
rm(list = ls())

# Programs
require(ecospat); require(knitr); require(maps); require(maptools);
require(raster); require(rgeos); require(sp); require(spThin)

# Data Preparation
setwd("~/Desktop/R/Data")
data_overall <- read.csv2("~/Desktop/R/Data/DataRAnalysis.csv")
occ_points_overall <- subset(data_overall, select = -c(fam_species, order_species, habitat))
occ_points_overall_naomit <- occ_points_overall[complete.cases(occ_points_overall[,3]),]
data(wrld_simpl)

##### Distribution of Structural Coloration #####
# Match Distribution of Occurrence Records with Absence and Presence of Structural Coloration to World Map
occ_points_absence_sc <- subset(occ_points_overall_naomit, sc == "0")
occ_points_absence_thin_sc <- thin(occ_points_absence_sc, verbose = F, lat.col = "lat_y", long.col = "long_x", spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, write.log.file = F, locs.thinned.list.return = T)

occ_points_presence_sc <- subset(occ_points_overall_naomit, sc == "1")
occ_points_presence_thin_sc <- thin(occ_points_presence_sc, verbose = F, lat.col = "lat_y", long.col = "long_x", spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, write.log.file = F, locs.thinned.list.return = T)

plot(wrld_simpl)
points(occ_points_absence_thin_sc[[1]], col = "gray80", pch = 16, cex = 0.5)
points(occ_points<presence_thin_sc[[1]], col = "indianred2", pch = 16, cex = 0.5)

