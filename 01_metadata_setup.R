# Reset Commands
closeAllConnections()
rm(list = ls())

##### SETTING META DATA #####
# Programs
library(dplyr)

# Data Preparation
setwd("~/Desktop/R/Data")
data_r <- read.csv("DataRAnalysis.csv", stringsAsFactors = FALSE)
avonet <- read.csv("AVONET1_BirdLife.csv", stringsAsFactors = FALSE)

# Select AVONET Variables
avonet_sub <- avonet %>% select(Species1, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
                                Tarsus.Length, Wing.Length, Kipps.Distance, Secondary1, Hand.Wing.Index,
                                Tail.Length, Mass, Habitat.Density, Migration, Trophic.Level, Trophic.Niche,
                                Primary.Lifestyle, Min.Latitude, Max.Latitude, Centroid.Latitude,
                                Centroid.Longitude, Range.Size)

# Remove Old Centroids
meta_data <- data_r %>% select(-long_x, -lat_y)

# Merge AVONET Variables Into Meta Data Set
meta_data <- meta_data %>% left_join(avonet_sub, by = c("species_tree" = "Species1")) %>%
  rename(
    min_lat             = Min.Latitude,
    max_lat             = Max.Latitude,
    long_x              = Centroid.Longitude,
    lat_y               = Centroid.Latitude,
    range_size          = Range.Size,
    habitat_density     = Habitat.Density,
    migration           = Migration,
    trophic_lvl         = Trophic.Level,
    trophic_niche       = Trophic.Niche,
    primary_lifestyle   = Primary.Lifestyle,
    beak_length_culmen  = Beak.Length_Culmen,
    beak_length_nares   = Beak.Length_Nares,
    beak_width          = Beak.Width,
    beak_depth          = Beak.Depth,
    tarsus_length       = Tarsus.Length,
    wing_length         = Wing.Length,
    kipps_distance      = Kipps.Distance,
    secondary           = Secondary1,
    hand_wing_index     = Hand.Wing.Index,
    tail_length         = Tail.Length,
    mass                = Mass
  )

# Order Columns
meta_data <- meta_data %>%
  select(species_tree, fam_species, order_species, sc, b_nisc, g_nisc, p_nisc, nisc, min_lat, max_lat,
         long_x, lat_y, range_size, habitat, habitat_density, migration, trophic_lvl, trophic_niche,
         primary_lifestyle, beak_length_culmen, beak_length_nares, beak_width, beak_depth, tarsus_length,
         wing_length, kipps_distance, secondary, hand_wing_index, tail_length, mass)

meta_data <- data.frame(as.list(meta_data))

rm(avonet, avonet_sub, data_r)
