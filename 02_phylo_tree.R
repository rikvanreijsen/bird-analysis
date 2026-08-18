##### PHYLOGENETIC TREE ANALYSIS #####
# Programs
library(ape); library(geiger); library(ggplot2); library(phytools); library(TreeTools)

# Data Preparation
phylo_tree_data <- meta_data %>%
  select(species_tree, fam_species, order_species, sc, b_nisc, g_nisc, p_nisc, nisc)

bird_tree <- read.tree("Stage2_Hackett_MCC_no_neg.tree")
is.ultrametric(bird_tree)
bird_tree

# Prune Data & Tree
rownames(phylo_tree_data) <- phylo_tree_data[,1]
names(phylo_tree_data)

pruned_phylo_tree <- treedata(bird_tree, phylo_tree_data, sort = T)$phy
name.check(bird_tree, phylo_tree_data)

pruned_phylo_tree_data <- treedata(bird_tree, phylo_tree_data, sort = T)$data
pruned_phylo_tree_data <- data.frame(pruned_phylo_tree_data)
name.check(pruned_phylo_tree, pruned_phylo_tree_data)

## Overall SC ##
pruned_phylo_tree_data$sc <- as.factor(pruned_phylo_tree_data$sc)
sc_overall <- setNames(pruned_phylo_tree_data$sc, pruned_phylo_tree_data$species_tree)

# Raw Ancestral State Reconstruction
onerate_sc <- ace(sc_overall, pruned_phylo_tree, type = "discrete", model = "ER")
tworate_sc <- ace(sc_overall, pruned_phylo_tree, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
logL_er  <- onerate_sc$loglik
logL_ard <- tworate_sc$loglik
lrt_stat <- 2 * (logL_ard - logL_er)
p_val    <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)

# Extract Transition Rates
q01_gain <- tworate_sc$rates[1]
q10_loss <- tworate_sc$rates[2]

# Format Output Table
single_trait_result <- data.frame(
  Phenotype    = "SC",
  logL_ER      = round(logL_er, 3),
  logL_ARD     = round(logL_ard, 3),
  LRT_Chi2     = round(lrt_stat, 3),
  P_Value      = format.pval(p_val, digits = 3),
  q01_Gain_ARD = round(q01_gain, 5),
  q10_Loss_ARD = round(q10_loss, 5),
  Best_Model   = ifelse(p_val < 0.05, "ARD", "ER")
)

print(single_trait_result)

# How Many Species & Families Display Trait
total_species  <- length(unique(pruned_phylo_tree_data$species_tree))
total_families <- length(unique(pruned_phylo_tree_data$fam_species))

species_sc <- aggregate(as.numeric(as.character(sc)) ~ species_tree, data = pruned_phylo_tree_data,
                        FUN = max, na.rm = TRUE)
family_sc  <- aggregate(as.numeric(as.character(sc)) ~ fam_species, data = pruned_phylo_tree_data, 
                        FUN = max, na.rm = TRUE)

colnames(species_sc)[2] <- "sc"
sc_species <- sum(species_sc$sc == 1)

colnames(family_sc)[2] <- "sc"
sc_families <- sum(family_sc$sc == 1)

cat("There are",
    sc_species, "out of", total_species, "SC species (", sc_species/total_species*100, "%), distributed among",
    sc_families, "out of", total_families, "families (", sc_families/total_families*100, "%)", "\n")

##### Overall NISC #####
pruned_phylo_tree_data$nisc <- as.factor(pruned_phylo_tree_data$nisc)
nisc_overall <- setNames(pruned_phylo_tree_data$nisc, pruned_phylo_tree_data$species_tree)

# Raw Ancestral State Reconstruction
onerate_nisc <- ace(nisc_overall, pruned_phylo_tree, type = "discrete", model = "ER")
tworate_nisc <- ace(nisc_overall, pruned_phylo_tree, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
logL_er  <- onerate_nisc$loglik
logL_ard <- tworate_nisc$loglik
lrt_stat <- 2 * (logL_ard - logL_er)
p_val    <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)

# Extract Transition Rates
q01_gain <- tworate_nisc$rates[1]
q10_loss <- tworate_nisc$rates[2]

# Format Output Table
single_trait_result <- data.frame(
  Phenotype    = "NISC",
  logL_ER      = round(logL_er, 3),
  logL_ARD     = round(logL_ard, 3),
  LRT_Chi2     = round(lrt_stat, 3),
  P_Value      = format.pval(p_val, digits = 3),
  q01_Gain_ARD = round(q01_gain, 5),
  q10_Loss_ARD = round(q10_loss, 5),
  Best_Model   = ifelse(p_val < 0.05, "ARD", "ER")
)

print(single_trait_result)

# How Many Species & Families Display Trait
species_nisc <- aggregate(as.numeric(as.character(nisc)) ~ species_tree, data = pruned_phylo_tree_data,
                          FUN = max, na.rm = TRUE)
family_nisc  <- aggregate(as.numeric(as.character(nisc)) ~ fam_species, data = pruned_phylo_tree_data,
                          FUN = max, na.rm = TRUE)

colnames(species_nisc)[2] <- "nisc"
nisc_species <- sum(species_nisc$nisc == 1)

colnames(family_nisc)[2] <- "nisc"
nisc_families <- sum(family_nisc$nisc == 1)

cat("There are",
    nisc_species, "out of", total_species, "NISC species (", nisc_species/total_species*100, "%), 
    distributed among",
    nisc_families, "out of", total_families, "families (", nisc_families/total_families*100, "%)", "\n")

family_summary <- pruned_phylo_tree_data %>%
  group_by(fam_species) %>%
  summarise(
    total_species = n(),
    nisc_species  = sum(nisc == 1, na.rm = TRUE),
    prop_nisc     = nisc_species / total_species,
    bnisc_species = sum(b_nisc == 1, na.rm = TRUE),
    gnisc_species = sum(g_nisc == 1, na.rm = TRUE),
    pnisc_species = sum(p_nisc == 1, na.rm = TRUE)
  ) %>%
  arrange(desc(nisc_species))

head(family_summary, 10)

family_summary %>%
  filter(total_species >= 10) %>%
  arrange(desc(prop_nisc)) %>%
  head(10)

##### Blue NISC #####
pruned_phylo_tree_data$b_nisc <- as.factor(pruned_phylo_tree_data$b_nisc)
b_nisc_overall <- setNames(pruned_phylo_tree_data$b_nisc, pruned_phylo_tree_data$species_tree)

# Raw Ancestral State Reconstruction
onerate_b_nisc <- ace(b_nisc_overall, pruned_phylo_tree, type = "discrete", model = "ER")
tworate_b_nisc <- ace(b_nisc_overall, pruned_phylo_tree, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
logL_er  <- onerate_b_nisc$loglik
logL_ard <- tworate_b_nisc$loglik
lrt_stat <- 2 * (logL_ard - logL_er)
p_val    <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)

# Extract Transition Rates
q01_gain <- tworate_b_nisc$rates[1]
q10_loss <- tworate_b_nisc$rates[2]

# Format Output Table
single_trait_result <- data.frame(
  Phenotype    = "BNISC",
  logL_ER      = round(logL_er, 3),
  logL_ARD     = round(logL_ard, 3),
  LRT_Chi2     = round(lrt_stat, 3),
  P_Value      = format.pval(p_val, digits = 3),
  q01_Gain_ARD = round(q01_gain, 5),
  q10_Loss_ARD = round(q10_loss, 5),
  Best_Model   = ifelse(p_val < 0.05, "ARD", "ER")
)

print(single_trait_result)

# How Many Species & Families Display Trait
species_b_nisc <- aggregate(as.numeric(as.character(b_nisc)) ~ species_tree, data = pruned_phylo_tree_data,
                            FUN = max, na.rm = TRUE)
family_b_nisc  <- aggregate(as.numeric(as.character(b_nisc)) ~ fam_species, data = pruned_phylo_tree_data,
                            FUN = max, na.rm = TRUE)

colnames(species_b_nisc)[2] <- "b_nisc"
b_nisc_species <- sum(species_b_nisc$b_nisc == 1)

colnames(family_b_nisc)[2] <- "b_nisc"
b_nisc_families <- sum(family_b_nisc$b_nisc == 1)

cat("There are",
    b_nisc_species, "out of", total_species, "BNISC species (", b_nisc_species/total_species*100, "%), 
    distributed among",
    b_nisc_families, "out of", total_families, "families (", b_nisc_families/total_families*100, "%)", "\n")

##### Green NISC #####
pruned_phylo_tree_data$g_nisc <- as.factor(pruned_phylo_tree_data$g_nisc)
g_nisc_overall <- setNames(pruned_phylo_tree_data$g_nisc, pruned_phylo_tree_data$species_tree)

# Raw Ancestral State Reconstruction
onerate_g_nisc <- ace(g_nisc_overall, pruned_phylo_tree, type = "discrete", model = "ER")
tworate_g_nisc <- ace(g_nisc_overall, pruned_phylo_tree, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
logL_er  <- onerate_g_nisc$loglik
logL_ard <- tworate_g_nisc$loglik
lrt_stat <- 2 * (logL_ard - logL_er)
p_val    <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)

# Extract Transition Rates
q01_gain <- tworate_g_nisc$rates[1]
q10_loss <- tworate_g_nisc$rates[2]

# Format Output Table
single_trait_result <- data.frame(
  Phenotype    = "GNISC",
  logL_ER      = round(logL_er, 3),
  logL_ARD     = round(logL_ard, 3),
  LRT_Chi2     = round(lrt_stat, 3),
  P_Value      = format.pval(p_val, digits = 3),
  q01_Gain_ARD = round(q01_gain, 5),
  q10_Loss_ARD = round(q10_loss, 5),
  Best_Model   = ifelse(p_val < 0.05, "ARD", "ER")
)

print(single_trait_result)

# How Many Species & Families Display Trait
species_g_nisc <- aggregate(as.numeric(as.character(g_nisc)) ~ species_tree, data = pruned_phylo_tree_data,
                            FUN = max, na.rm = TRUE)
family_g_nisc  <- aggregate(as.numeric(as.character(g_nisc)) ~ fam_species, data = pruned_phylo_tree_data,
                            FUN = max, na.rm = TRUE)

colnames(species_g_nisc)[2] <- "g_nisc"
g_nisc_species <- sum(species_g_nisc$g_nisc == 1)

colnames(family_g_nisc)[2] <- "g_nisc"
g_nisc_families <- sum(family_g_nisc$g_nisc == 1)

cat("There are",
    g_nisc_species, "out of", total_species, "GNISC species (", g_nisc_species/total_species*100, "%), 
    distributed among",
    g_nisc_families, "out of", total_families, "families (", g_nisc_families/total_families*100, "%)", "\n")

##### Purple NISC #####
pruned_phylo_tree_data$p_nisc <- as.factor(pruned_phylo_tree_data$p_nisc)
p_nisc_overall <- setNames(pruned_phylo_tree_data$p_nisc, pruned_phylo_tree_data$species_tree)

# Raw Ancestral State Reconstruction
onerate_p_nisc <- ace(p_nisc_overall, pruned_phylo_tree, type = "discrete", model = "ER")
tworate_p_nisc <- ace(p_nisc_overall, pruned_phylo_tree, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
logL_er  <- onerate_p_nisc$loglik
logL_ard <- tworate_p_nisc$loglik
lrt_stat <- 2 * (logL_ard - logL_er)
p_val    <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)

# Extract Transition Rates
q01_gain <- tworate_p_nisc$rates[1]
q10_loss <- tworate_p_nisc$rates[2]

# Format Output Table
single_trait_result <- data.frame(
  Phenotype    = "PNISC",
  logL_ER      = round(logL_er, 3),
  logL_ARD     = round(logL_ard, 3),
  LRT_Chi2     = round(lrt_stat, 3),
  P_Value      = format.pval(p_val, digits = 3),
  q01_Gain_ARD = round(q01_gain, 5),
  q10_Loss_ARD = round(q10_loss, 5),
  Best_Model   = ifelse(p_val < 0.05, "ARD", "ER")
)

print(single_trait_result)

# How Many Species & Families Display Trait
species_p_nisc <- aggregate(as.numeric(as.character(p_nisc)) ~ species_tree, data = pruned_phylo_tree_data,
                            FUN = max, na.rm = TRUE)
family_p_nisc  <- aggregate(as.numeric(as.character(p_nisc)) ~ fam_species, data = pruned_phylo_tree_data,
                            FUN = max, na.rm = TRUE)

colnames(species_p_nisc)[2] <- "p_nisc"
p_nisc_species <- sum(species_p_nisc$p_nisc == 1)

colnames(family_p_nisc)[2] <- "p_nisc"
p_nisc_families <- sum(family_p_nisc$p_nisc == 1)

cat("There are",
    p_nisc_species, "out of", total_species, "PNISC species (", p_nisc_species/total_species*100, "%), 
    distributed among",
    p_nisc_families, "out of", total_families, "families (", p_nisc_families/total_families*100, "%)", "\n")

##### Proportional Representation Rate #####
trait_cols <- c("sc", "nisc", "b_nisc", "g_nisc", "p_nisc")

table_ratios <- sapply(trait_cols, function(col) {
  mean(pruned_phylo_tree_data[[col]] %in% c(1, "1", "Present", "TRUE"), na.rm = TRUE)
})

round(table_ratios, 3)

rm(family_b_nisc, family_g_nisc, family_p_nisc, family_sc,
   onerate_sc, onerate_b_nisc, onerate_g_nisc, onerate_p_nisc, onerate_nisc,
   sc_overall, b_nisc_overall, g_nisc_overall, p_nisc_overall, nisc_overall,
   trait_cols)
