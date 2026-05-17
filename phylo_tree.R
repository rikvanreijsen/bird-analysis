# Reset Commands
closeAllConnections()
rm(list = ls())

# Programs
library(ape); library(geiger); library(ggplot2); library(phytools); library(TreeTools)

##### Overall Data #####
# Data Preparation
setwd("~/Desktop/R/Data")
bird_tree <- read.tree("Stage2_Hackett_MCC_no_neg.tree")
is.ultrametric(bird_tree)
bird_tree
data_overall <- read.csv2("~/Desktop/R/Data/DataRAnalysis.csv")
data_overall <- data.frame(as.list(data_overall))

# Prune Data & Tree
rownames(data_overall) <- data_overall[,1]
row.names(data_overall)
names(data_overall)

pruned_tree_overall <- treedata(bird_tree, data_overall, sort = T)$phy
name.check(bird_tree, data_overall)

pruned_data_overall <- treedata(bird_tree, data_overall, sort = T)$data
pruned_data_overall <- data.frame(pruned_data_overall)
name.check(pruned_tree_overall, pruned_data_overall)

##### Overall SC #####
pruned_data_overall$sc <- as.factor(pruned_data_overall$sc)
sc_overall <- setNames(pruned_data_overall$sc, pruned_data_overall$species_tree)

# Raw Ancestral State Reconstruction
onerate_sc <- ace(sc_overall, pruned_tree_overall, type = "discrete", model = "ER")
tworate_sc <- ace(sc_overall, pruned_tree_overall, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
anova(onerate_sc, tworate_sc)
onerate_sc
tworate_sc

# Plotting Tree
plot.phylo(pruned_tree_overall, type = "fan", no.margin = T, show.tip.label = F, edge.width = 0.2,
           cex = 0.5, edge.color = "gray50")
tiplabels(pch = 20, col = c("gray95", "lightcoral")[as.factor(pruned_data_overall$sc)], lwd = 5)
nodelabels(pie = tworate_sc$lik.anc, piecol = setNames(c("whitesmoke", "lightcoral"), c(0, 1)), 
           cex = 0.2)

##### Overall NISC #####
pruned_data_overall$nisc <- as.factor(pruned_data_overall$nisc)
nisc_overall <- setNames(pruned_data_overall$nisc, pruned_data_overall$species_tree)

# Raw Ancestral State Reconstruction
onerate_nisc <- ace(nisc_overall, pruned_tree_overall, type = "discrete", model = "ER")
tworate_nisc <- ace(nisc_overall, pruned_tree_overall, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
anova(onerate_nisc, tworate_nisc)
onerate_nisc
tworate_nisc

# Plotting Tree
plot.phylo(pruned_tree_overall, type = "fan", no.margin = T, show.tip.label = F, edge.width = 0.2, 
           cex = 0.5, edge.color = "gray50")
tiplabels(pch = 20, col = c("gray95", "lightcoral")[as.factor(pruned_data_overall$nisc)], lwd = 5)
nodelabels(pie = tworate_nisc$lik.anc, piecol = setNames(c("whitesmoke", "lightcoral"), c(0, 1)), 
           cex = 0.2)

##### Blue NISC #####
pruned_data_overall$b_nisc <- as.factor(pruned_data_overall$b_nisc)
b_nisc_overall <- setNames(pruned_data_overall$b_nisc, pruned_data_overall$species_tree)

# Raw Ancestral State Reconstruction
onerate_b_nisc <- ace(b_nisc_overall, pruned_tree_overall, type = "discrete", model = "ER")
tworate_b_nisc <- ace(b_nisc_overall, pruned_tree_overall, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
anova(onerate_b_nisc, tworate_b_nisc)
onerate_b_nisc
tworate_b_nisc

# Plotting Tree
plot.phylo(pruned_tree_overall, type = "fan", no.margin = T, show.tip.label = F, edge.width = 0.2, 
           cex = 0.5, edge.color = "gray50")
tiplabels(pch = 20, col = c("gray95", "royalblue")[as.factor(pruned_data_overall$b_nisc)], lwd = 5)

##### Green NISC #####
pruned_data_overall$g_nisc <- as.factor(pruned_data_overall$g_nisc)
g_nisc_overall <- setNames(pruned_data_overall$g_nisc, pruned_data_overall$species_tree)

# Raw Ancestral State Reconstruction
onerate_g_nisc <- ace(g_nisc_overall, pruned_tree_overall, type = "discrete", model = "ER")
tworate_g_nisc <- ace(g_nisc_overall, pruned_tree_overall, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
anova(onerate_g_nisc, tworate_g_nisc)
onerate_g_nisc
tworate_g_nisc

# Plotting Tree
plot.phylo(pruned_tree_overall, type = "fan", no.margin = T, show.tip.label = F, edge.width = 0.2, 
           cex = 0.5, edge.color = "gray50")
tiplabels(pch = 20, col = c("gray95", "mediumseagreen")[as.factor(pruned_data_overall$g_nisc)], lwd = 5)

##### Purple NISC #####
pruned_data_overall$p_nisc <- as.factor(pruned_data_overall$p_nisc)
p_nisc_overall <- setNames(pruned_data_overall$p_nisc, pruned_data_overall$species_tree)

# Raw Ancestral State Reconstruction
onerate_p_nisc <- ace(p_nisc_overall, pruned_tree_overall, type = "discrete", model = "ER")
tworate_p_nisc <- ace(p_nisc_overall, pruned_tree_overall, type = "discrete", model = "ARD")

# Likelihood Ratio for Ancestral Character Estimation
anova(onerate_p_nisc, tworate_p_nisc)
onerate_p_nisc
tworate_p_nisc

# Plotting Tree
plot.phylo(pruned_tree_overall, type = "fan", no.margin = T, show.tip.label = F, edge.width = 0.2, 
           cex = 0.5, edge.color = "gray50")
tiplabels(pch = 20, col = c("gray95", "mediumpurple1")[as.factor(pruned_data_overall$p_nisc)], lwd = 5)
