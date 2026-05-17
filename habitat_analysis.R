# Reset Commands
closeAllConnections()
rm(list = ls())

# Programs
library(ecospat); library(knitr); library(maps); library(maptools)
library(raster); library(rgeos); library(sp); library(spThin)

# Data Preparation
setwd("~/Desktop/R/Data")
data_overall <- read.csv2("~/Desktop/R/Data/DataRAnalysis.csv")
occ_points_overall <- subset(data_overall, select = -c(fam_species, order_species, habitat))
occ_points_overall <- occ_points_overall[complete.cases(occ_points_overall[, 3]), ]
data(wrld_simpl)

##### NISC DISTRIBUTION ANALYSIS #####

# Define Categories and Colors
nisc_cols <- c("b_nisc", "g_nisc", "p_nisc")
occ_points_overall$nisc_count <- rowSums(occ_points_overall[, nisc_cols], na.rm = T)

occ_points_overall$category <- "absence"
occ_points_overall$category[occ_points_overall$sc == 1] <- "sc"
occ_points_overall$category[occ_points_overall$b_nisc == 1] <- "blue"
occ_points_overall$category[occ_points_overall$g_nisc == 1] <- "green"
occ_points_overall$category[occ_points_overall$p_nisc == 1] <- "purple"
occ_points_overall$category[occ_points_overall$nisc_count > 1] <- "multi"

col_map <- c("absence" = "gray80", "sc" = "gold", "blue" = "royalblue", 
                "green" = "seagreen", "purple" = "slateblue2", "multi" = "sienna2")

# Global Thinning (10 km to standardize search effort)
occ_points_thin <- thin(occ_points_overall, lat.col = "lat_y", long.col = "long_x", 
                        spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                        write.log.file = F, locs.thinned.list.return = T)


occ_points_thin <- occ_points_thin[[1]]
occ_points_thin <- merge(occ_points_thin, occ_points_overall, by.x = c("Longitude", "Latitude"),
                         by.y = c("long_x", "lat_y"))

# Convert to Spatial Object
pts <- SpatialPointsDataFrame(coords = occ_points_thin[, c("Longitude", "Latitude")], 
                              data = occ_points_thin, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts <- spTransform(pts, CRS(proj4string(wrld_simpl)))

slot(pts, "proj4string") <- slot(wrld_simpl, "proj4string")

region_lookup <- c("2" = "Africa", "19" = "Americas", "150" = "Europe", 
                   "142" = "Asia", "9" = "Oceania", "0" = "Antarctica")

##### LOOP BY REGION #####
for (cont_code in unique(wrld_simpl$REGION)) {
  cont_name <- region_lookup[as.character(cont_code)]
  if(is.na(cont_name)) cont_name <- paste("Region", cont_code)
  
  cont_map <- wrld_simpl[wrld_simpl$REGION == cont_code, ]
  
  # Buffer to capture coastal breeding centroids
  cont_map     <- gBuffer(cont_map, width = 0.1, byid = F)
  check_inside <- over(pts, cont_map)
  cont_pts     <- pts[!is.na(check_inside), ]
  
  if(length(cont_pts) == 0) next
  
  ext   <- extent(cont_pts)
  x_buf <- (ext@xmax - ext@xmin) * 0.1
  y_buf <- (ext@ymax - ext@ymin) * 0.1
  
  # Plotting
  plot(cont_map, main = paste("NISC Breeding Prevalence in", cont_name), border = "gray60", 
       xlim = c(ext@xmin - x_buf, ext@xmax + x_buf), ylim = c(ext@ymin - y_buf, ext@ymax + y_buf),
       asp = 1)
  
  # Legend/visual hint: focus is on colored NISC points (blue, green, purple, multi)
  points(cont_pts, col = col_map[as.character(cont_pts$category)], 
         pch = ifelse(cont_pts$category == "absence", 4, 
                      ifelse(cont_pts$category == "sc", 8, 16)), cex = 0.7)
}

##### GEOGRAPHIC TABLES & PERCENTAGE CALCULATION #####
slot(pts, "proj4string")  <- slot(wrld_simpl, "proj4string")
geo_data                  <- over(pts, wrld_simpl)
occ_points_thin$country   <- geo_data$NAME
occ_points_thin$region    <- as.character(geo_data$REGION)
occ_points_thin$subregion <- as.character(geo_data$SUBREGION)

sub_names <- c("15" = "Northern Africa", "11" = "Western Africa", "17" = "Middle Africa", 
               "14" = "Eastern Africa", "18" = "Southern Africa", "143" = "Central Asia", 
               "145" = "Western Asia", "34" = "Southern Asia", "35" = "South-Eastern Asia", 
               "30" = "Eastern Asia", "154" = "Northern Europe", "155" = "Western Europe", 
               "151" = "Eastern Europe", "39" = "Southern Europe", "21" = "Northern America", 
               "29" = "Caribbean", "13" = "Central America", "5" = "South America", 
               "53" = "Australia and New Zealand", "54" = "Melanesia", "57" = "Micronesia", 
               "61" = "Polynesia")

occ_points_thin$region    <- region_lookup[occ_points_thin$region]
occ_points_thin$subregion <- sub_names[occ_points_thin$subregion]

# Tables (Absence, SC, blue, green, purple, multi)
# Indices in table: 1 = Absence, 2 = SC, 3-5 = NISC categories
country_occ   <- as.data.frame.matrix(table(occ_points_thin$country, occ_points_thin$category)[, c(1, 5, 2:4)])
region_occ    <- as.data.frame.matrix(table(occ_points_thin$region, occ_points_thin$category)[, c(1, 5, 2:4)])
subregion_occ <- as.data.frame.matrix(table(occ_points_thin$subregion, occ_points_thin$category)[, c(1, 5, 2:4)])

# Exclude SC from tables: 1 = Absence, 2 = Blue NISC, 3 = Green NISC, 4 = Multi
country_occ   <- subset(country_occ[, -2])
region_occ    <- subset(region_occ[, -2])
subregion_occ <- subset(subregion_occ[, -2])

# --- Calculate Proportional Breeding Prevalence ---

# 1. Total thinned sites (Demoninator)
country_occ$total_sites   <- rowSums(country_occ)
region_occ$total_sites    <- rowSums(region_occ)
subregion_occ$total_sites <- rowSums(subregion_occ)

# 2. Total NISC counts (Numerator)
# Sum of NISC color columns + multi
country_occ$nisc_total_count   <- rowSums(country_occ[, 2:4])
region_occ$nisc_total_count    <- rowSums(region_occ[, 2:4])
subregion_occ$nisc_total_count <- rowSums(subregion_occ[, 2:4])

# Calculate Percentage (%)
country_occ$perc_nisc_prevalence   <- (country_occ$nisc_total_count / country_occ$total_sites) * 100
region_occ$perc_nisc_prevalence    <- (region_occ$nisc_total_count / region_occ$total_sites) * 100
subregion_occ$perc_nisc_prevalence <- (subregion_occ$nisc_total_count / subregion_occ$total_sites) * 100

# Exclude Countries, Regions, and Subregions with No Proportional Breeding Prevalence
country_occ   <- subset(country_occ, is.finite(perc_nisc_prevalence) & perc_nisc_prevalence > 0)
region_occ    <- subset(region_occ, is.finite(perc_nisc_prevalence) & perc_nisc_prevalence > 0)
subregion_occ <- subset(subregion_occ, is.finite(perc_nisc_prevalence) & perc_nisc_prevalence > 0)

# Sort table so highest occurrences are at top
country_occ   <- country_occ[order(-country_occ$perc_nisc_prevalence), ]
region_occ    <- region_occ[order(-region_occ$perc_nisc_prevalence), ]
subregion_occ <- subregion_occ[order(-subregion_occ$perc_nisc_prevalence), ]

# View top of list
head(country_occ)
head(region_occ)
head(subregion_occ)

##### INCLUDE UNKNOWN POINTS #####
points_unknown <- occ_points_thin[is.na(occ_points_thin$country), ]
plot(wrld_simpl, border = "gray60")
points(points_unknown$Longitude, points_unknown$Latitude, col = "red", pch = 4, cex = 0.7)

# Combine thinned points and new points
points_missing <- read.csv("~/Desktop/R/Data/points_missing.csv")
occ_points_com <- rbind(occ_points_thin[!is.na(occ_points_thin$country), ], points_missing)

# Geographic tables and percentage calculation
# Tables (Absence, SC, blue, green, purple, multi)
# Indices in table: 1 = Absence, 2 = SC, 3-5 = NISC categories
country_occ   <- as.data.frame.matrix(table(occ_points_com$country, occ_points_com$category)[, c(1, 5, 2:4)])
region_occ    <- as.data.frame.matrix(table(occ_points_com$region, occ_points_com$category)[, c(1, 5, 2:4)])
subregion_occ <- as.data.frame.matrix(table(occ_points_com$subregion, occ_points_com$category)[, c(1, 5, 2:4)])

# Exclude SC from tables: 1 = Absence, 2 = Blue NISC, 3 = Green NISC, 4 = Multi
country_occ   <- subset(country_occ[, -2])
region_occ    <- subset(region_occ[, -2])
subregion_occ <- subset(subregion_occ[, -2])

# --- Calculate Proportional Breeding Prevalence ---

# 1. Total thinned sites (Demoninator)
country_occ$total_sites   <- rowSums(country_occ)
region_occ$total_sites    <- rowSums(region_occ)
subregion_occ$total_sites <- rowSums(subregion_occ)

# 2. Total NISC counts (Numerator)
# Sum of NISC color columns + multi
country_occ$nisc_total_count   <- rowSums(country_occ[, 2:4])
region_occ$nisc_total_count    <- rowSums(region_occ[, 2:4])
subregion_occ$nisc_total_count <- rowSums(subregion_occ[, 2:4])

# Calculate Percentage (%)
country_occ$perc_nisc_prevalence   <- (country_occ$nisc_total_count / country_occ$total_sites) * 100
region_occ$perc_nisc_prevalence    <- (region_occ$nisc_total_count / region_occ$total_sites) * 100
subregion_occ$perc_nisc_prevalence <- (subregion_occ$nisc_total_count / subregion_occ$total_sites) * 100

# Exclude Countries, Regions, and Subregions with No Proportional Breeding Prevalence
country_occ   <- subset(country_occ, is.finite(perc_nisc_prevalence) & perc_nisc_prevalence > 0)
region_occ    <- subset(region_occ, is.finite(perc_nisc_prevalence) & perc_nisc_prevalence > 0)
subregion_occ <- subset(subregion_occ, is.finite(perc_nisc_prevalence) & perc_nisc_prevalence > 0)

# Sort table so highest occurrences are at top
country_occ   <- country_occ[order(-country_occ$perc_nisc_prevalence), ]
region_occ    <- region_occ[order(-region_occ$perc_nisc_prevalence), ]
subregion_occ <- subregion_occ[order(-subregion_occ$perc_nisc_prevalence), ]

# View top of list
head(country_occ)
head(region_occ)
head(subregion_occ)

##### DISTRIBUTION BASED ON FAMILY, ORDER, AND HABITAT #####
# Data Preparation
occ_points_char <- subset(data_overall, select = -c(sc, b_nisc, g_nisc, p_nisc, nisc))
occ_points_char <- occ_points_char[complete.cases(occ_points_char[, 5]), ]

# Global Thinning
occ_points_char_thin <- thin(occ_points_char, lat.col = "lat_y", long.col = "long_x", 
                        spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                        write.log.file = F, locs.thinned.list.return = T)

occ_points_char_thin <- occ_points_char_thin[[1]]
occ_points_char_thin <- merge(occ_points_char_thin, occ_points_char, by.x = c("Longitude", "Latitude"),
                         by.y = c("long_x", "lat_y"))

# --- Plot by Family with Gray 'Other' ---
top_fam <- names(sort(table(occ_points_char_thin$fam_species), decreasing = T))[1:10]
occ_points_char_thin$fam_plot <- ifelse(occ_points_char_thin$fam_species %in% top_fam, 
                                   as.character(occ_points_char_thin$fam_species), "Other")
# Define colors
unique_fams <- setdiff(unique(occ_points_char_thin$fam_plot), "Other")
fam_col <- setNames(viridis(length(unique_fams)), unique_fams)
fam_col["Other"] <- "gray80"

# Plotting
plot(wrld_simpl, border = "gray60", main = "Global Bird Distribution by Family")
points(occ_points_char_thin$Longitude, occ_points_char_thin$Latitude, 
       col = fam_col[occ_points_char_thin$fam_plot], pch = 4, cex = 0.7)

# --- Plot by Order with Gray 'Other' ---
top_order <- names(sort(table(occ_points_char_thin$order_species), decreasing = T))[1:10]
occ_points_char_thin$order_plot <- ifelse(occ_points_char_thin$order_species %in% top_order,
                                          as.character(occ_points_char_thin$order_species), "Other")

# Define colors
unique_order <- setdiff(unique(occ_points_char_thin$order_plot), "Other")
order_col <- setNames(viridis(length(unique_order)), unique_order)
order_col["Other"] <- "gray80"

# Plotting
plot(wrld_simpl, border = "gray60", main = "Global Bird Distribution by Order")
points(occ_points_char_thin$Longitude, occ_points_char_thin$Latitude,
       col = order_col[occ_points_char_thin$order_plot], pch = 4, cex = 0.7)

# --- Plot by Habitat ---
top_hab <- names(table(occ_points_char_thin$habitat))

# Define colors
hab_col <- c("1" = "forestgreen", "2" = "gold2", "3" = "olivedrab3", "4" = "cyan4", 
             "5" = "royalblue3", "6" = "darkorange2")

# Plotting
plot(wrld_simpl, border = "gray60", main = "Global Bird Distribution by Habitat")
points(occ_points_char_thin$Longitude, occ_points_char_thin$Latitude,
       col = hab_col[occ_points_char_thin$habitat], pch = 4, cex = 0.7)
