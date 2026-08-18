##### GEOSPATIAL ANALYSIS #####
# Programs
library(ecospat); library(knitr); library(maps); library(maptools); library(raster); library(rgeos); library(sp);
library(spThin); library(sf); library(rnaturalearth); library(rnaturalearthdata)

# Data Preparation
geospat_data <- meta_data %>%
  select(species_tree, sc, b_nisc, g_nisc, p_nisc, nisc, long_x, lat_y)

geospat_data <- geospat_data[complete.cases(geospat_data[, 8]), ]

##### NISC Distribution #####
# Define Categories and Colors
nisc_cols <- c("b_nisc", "g_nisc", "p_nisc")
geospat_data$nisc_count <- rowSums(geospat_data[, nisc_cols], na.rm = T)

geospat_data$category <- "absence"
geospat_data$category[geospat_data$sc == 1] <- "sc"
geospat_data$category[geospat_data$b_nisc == 1] <- "blue"
geospat_data$category[geospat_data$g_nisc == 1] <- "green"
geospat_data$category[geospat_data$p_nisc == 1] <- "purple"
geospat_data$category[geospat_data$nisc_count > 1] <- "multi"

# Global Thinning (10 km to standardize search effort)
geospat_data_thin <- thin(geospat_data, lat.col = "lat_y", long.col = "long_x", 
                        spec.col = "species_tree", thin.par = 10, reps = 1, write.files = F, 
                        write.log.file = F, locs.thinned.list.return = T)

geospat_data_thin <- geospat_data_thin[[1]]
geospat_data_thin <- merge(geospat_data_thin, geospat_data, by.x = c("Latitude", "Longitude"),
                         by.y = c("lat_y", "long_x"))

print(nrow(geospat_data_thin))

##### Locate Bird Species #####
# Identify Original Column Names
original_cols <- colnames(geospat_data_thin)

# Load World Map Shapefile
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Convert ALL Thinned Data to Spatial Points (WGS84 projection)
geospat_sf <- st_as_sf(
  geospat_data_thin,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# Enable Spherical/Geodesic Geometry Engine
sf_use_s2(TRUE)

# Global Spatial Join on All Thinned Points
joined_all <- st_join(geospat_sf, world, join = st_intersects)

geospat_data_thin_all <- joined_all %>%
  st_drop_geometry() %>%
  mutate(
    country   = name,
    region    = ifelse(!is.na(region_un), as.character(region_un), as.character(continent)),
    subregion = as.character(subregion)
  )

# Resolve Missing Points or False Antarctica Matches (Latitude > -60)
points_unknown_sf <- geospat_sf %>%
  filter(
    is.na(geospat_data_thin_all$country) | 
      trimws(geospat_data_thin_all$country) == "" | 
      (geospat_data_thin_all$country == "Antarctica" & Latitude > -60)
  )

if (nrow(points_unknown_sf) > 0) {
  # Find nearest land polygon using exact spherical distances
  nearest_idx   <- st_nearest_feature(points_unknown_sf, world)
  nearest_world <- world[nearest_idx, ]
  
  points_unknown_resolved <- points_unknown_sf %>%
    st_drop_geometry() %>%
    mutate(
      country   = nearest_world$name,
      region    = ifelse(!is.na(nearest_world$region_un), 
                         as.character(nearest_world$region_un), 
                         as.character(nearest_world$continent)),
      subregion = as.character(nearest_world$subregion)
    )
  
  # Separate valid points
  geospat_data_thin_known <- geospat_data_thin_all %>%
    filter(
      !is.na(country) & 
        trimws(country) != "" & 
        !(country == "Antarctica" & Latitude > -60)
    )
  
  # Combine back together
  geospat_data_thin_updated <- bind_rows(geospat_data_thin_known, points_unknown_resolved)
  
} else {
  geospat_data_thin_updated <- geospat_data_thin_all
}

# Strict Column Selection (Keep ONLY Original Columns + Country, Region, Subregion)
target_cols <- unique(c(original_cols, "name_long", "region", "subregion"))

geospat_data_thin_updated <- geospat_data_thin_updated %>% select(all_of(target_cols))
geospat_data_thin_updated <- geospat_data_thin_updated %>% rename(country = name_long)
geospat_data_thin <- geospat_data_thin_updated

##### Data Extraction #####
# 1. Total Thinned Points
total_thinned_points <- nrow(geospat_data_thin)

# 2. Reclassified / Distance-Resolved Points (Offshore / False Antarctica)
total_reclassified_points <- nrow(points_unknown_sf)

# 3. Calculate Regional Prevalence Matrix
prevalence_matrix <- geospat_data_thin %>%
  group_by(region, subregion, country) %>%
  summarise(
    total_sites = n(),
    nisc_present = sum(nisc == 1, na.rm = TRUE),
    prevalence = (nisc_present / total_sites) * 100,
    .groups = "drop"
  ) %>%
  filter(total_sites >= 5) # Filter out regions with very low sample counts (e.g. < 5 records)

# 4. Summary Metrics Across Countries/Regions
total_polygons_analyzed <- nrow(prevalence_matrix)
min_prev  <- min(prevalence_matrix$prevalence, na.rm = TRUE)
mean_prev <- mean(prevalence_matrix$prevalence, na.rm = TRUE)
max_prev  <- max(prevalence_matrix$prevalence, na.rm = TRUE)

# 5. Top 5 Hotspots
top_hotspots <- prevalence_matrix %>%
  arrange(desc(prevalence)) %>%
  select(country, region, subregion, total_sites, nisc_present, prevalence) %>%
  head(5)

# --- PRINT FINAL RESULTS ---
cat("\n================ SECTION 3 METRICS ================\n") &
cat("1. Total Thinned Spatial Points:", total_thinned_points, "\n") &
cat("2. Reclassified/Resolved Points:", total_reclassified_points, "\n") &
cat("3. Total Geographic Units Analyzed:", total_polygons_analyzed, "\n") &
cat("4. Prevalence Range (%): Min =", round(min_prev, 2), 
    "| Mean =", round(mean_prev, 2), 
    "| Max =", round(max_prev, 2), "\n\n") &
cat("Top 5 Hotspot Regions:\n") &
print(top_hotspots) &
cat("===================================================\n")

##### Geographic Tables & Percentage Calculation #####
# 1. Calculate Regional Prevalence with Cleaned Polygons & N >= 10 Threshold
prevalence_matrix_hues <- geospat_data_thin %>%
  filter(!is.na(country)) %>% # Remove offshore/unassigned spatial points
  group_by(region, subregion, country) %>%
  summarise(
    total_sites   = n(),
    
    # Occurrence counts
    nisc_count    = sum(nisc == 1, na.rm = TRUE),
    bnisc_count   = sum(b_nisc == 1, na.rm = TRUE),
    gnisc_count   = sum(g_nisc == 1, na.rm = TRUE),
    pnisc_count   = sum(p_nisc == 1, na.rm = TRUE),
    
    # Proportional Prevalence (%)
    prev_NISC     = (nisc_count / total_sites) * 100,
    prev_BNISC    = (bnisc_count / total_sites) * 100,
    prev_GNISC    = (gnisc_count / total_sites) * 100,
    prev_PNISC    = (pnisc_count / total_sites) * 100,
    .groups       = "drop"
  ) %>%
  filter(total_sites >= 10) # Stricter sampling threshold to eliminate small-N noise

# 2. Extract Top 5 Hotspots per Category
# Overall NISC
top_nisc_hotspots <- prevalence_matrix_hues %>%
  arrange(desc(prev_NISC)) %>%
  select(country, region, total_sites, nisc_count, prev_NISC) %>%
  head(5)

# Blue NISC (BNISC)
top_blue_hotspots <- prevalence_matrix_hues %>%
  arrange(desc(prev_BNISC)) %>%
  select(country, region, total_sites, bnisc_count, prev_BNISC) %>%
  head(5)

# Green NISC (GNISC)
top_green_hotspots <- prevalence_matrix_hues %>%
  arrange(desc(prev_GNISC)) %>%
  select(country, region, total_sites, gnisc_count, prev_GNISC) %>%
  head(5)

# Purple NISC (PNISC)
top_purple_hotspots <- prevalence_matrix_hues %>%
  arrange(desc(prev_PNISC)) %>%
  select(country, region, total_sites, pnisc_count, prev_PNISC) %>%
  head(5)

# 3. Output Results Cleanly
print("=== TOP OVERALL NISC HOTSPOTS (N >= 10) ==="); print(top_nisc_hotspots)
print("=== TOP BLUE NISC HOTSPOTS (N >= 10) ==="); print(top_blue_hotspots)
print("=== TOP GREEN NISC HOTSPOTS (N >= 10) ==="); print(top_green_hotspots)
print("=== TOP PURPLE NISC HOTSPOTS (N >= 10) ==="); print(top_purple_hotspots)

rm(original_cols, world, geospat_sf, joined_all, geospat_data_thin_all, points_unknown_sf, nearest_idx,
   nearest_world, points_unknown_resolved, geospat_data_thin_known, geospat_data_thin_updated)

##### MULTIDIMENSIONAL MACROECOLOGICAL, ECOMORPHOLOGICAL & OCEANOGRAPHIC ANALYSIS OF AVIAN NISC #####
# Programs
library(sdmpredictors); library(terra)

# Data Preparation
phylolm_data <- meta_data

# Prune Phylolm Data & Tree
rownames(phylolm_data) <- phylolm_data[, 1]
row.names(phylolm_data)
names(phylolm_data)

pruned_phylolm_tree <- treedata(bird_tree, phylolm_data, sort = T)$phy
name.check(bird_tree, phylolm_data)

pruned_phylolm_data <- treedata(bird_tree, phylolm_data, sort = T)$data
pruned_phylolm_data <- data.frame(pruned_phylolm_data)
name.check(pruned_phylolm_tree, pruned_phylolm_data)

# Force Numeric Conversions on Columns
num_cols <- c("sc", "b_nisc", "g_nisc", "p_nisc", "nisc", "min_lat", "max_lat", "long_x", "lat_y", "range_size",
              "habitat", "habitat_density", "migration", "beak_length_culmen", "beak_length_nares", "beak_width",
              "beak_depth", "tarsus_length", "wing_length", "kipps_distance", "secondary", "hand_wing_index",
              "tail_length", "mass")

for (col in num_cols) {
  if (col %in% colnames(pruned_phylolm_data)) {
    pruned_phylolm_data[[col]] <- as.numeric(as.character(pruned_phylolm_data[[col]]))
  }
}

# Calculated Features (Global)
pruned_phylolm_data$lat_span  <- pruned_phylolm_data$max_lat - pruned_phylolm_data$min_lat
pruned_phylolm_data$abs_lat   <- abs(pruned_phylolm_data$lat_y)
pruned_phylolm_data$log_range <- log(pruned_phylolm_data$range_size)
pruned_phylolm_data$log_mass  <- log(pruned_phylolm_data$mass)

# Factorize categorical variables
pruned_phylolm_data$trophic_lvl       <- as.factor(pruned_phylolm_data$trophic_lvl)
pruned_phylolm_data$trophic_niche     <- as.factor(pruned_phylolm_data$trophic_niche)
pruned_phylolm_data$primary_lifestyle <- as.factor(pruned_phylolm_data$primary_lifestyle)

rm(col, num_cols)
