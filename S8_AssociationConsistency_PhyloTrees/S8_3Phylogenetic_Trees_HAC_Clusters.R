
# install.packages("metacoder")
# https://github.com/grunwaldlab/metacoder
library(metacoder)
library(dplyr)

supragingival_species_lineage_df <- data.frame(read.delim("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4supragingival_species_HAC90_lineage.txt", header=FALSE))
colnames(supragingival_species_lineage_df) <- c("tax_id","lineage")

supragingival_species_lineage_df$lineage <- trimws(supragingival_species_lineage_df$lineage)
supragingival_species_lineage_df$lineage <- gsub("\\[|\\]", "", supragingival_species_lineage_df$lineage)
supragingival_species_lineage_df$lineage <- gsub("/", "_", supragingival_species_lineage_df$lineage)
supragingival_species_lineage_df$lineage <- sub("^cellular organisms;", "", supragingival_species_lineage_df$lineage)

## Add species name from lineage 
supragingival_species_lineage_df$species <- sapply(strsplit(supragingival_species_lineage_df$lineage, ";"), tail, 1)

## Parse taxonomy
supragingival_tax_obj <- parse_tax_data(
  supragingival_species_lineage_df,
  class_cols = "lineage",
  class_sep = ";"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_3Supragingival_HAC_metacoder_All.pdf", width = 22, height = 12)

heat_tree(
  supragingival_tax_obj,
  node_label = taxon_names,
  # Nodes
  node_size = n_obs,
  node_size_range = c(0.02, 0.04),
  
  # Branches
  edge_size = n_obs,
  edge_size_range = c(0.008, 0.015),
  
  # Labels
  node_label_size = n_obs,
  #node_label_size_range = c(0.02, 0.04),
  node_color = "gray60",
  edge_color = "gray70",
  layout = "kamada-kawai",
  repel_labels = TRUE,
  repel_force = 4,
  make_node_legend = FALSE,
  make_edge_legend = FALSE
)

dev.off()

####### Now plot the same thing for cluster 2 species only.
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Supragingival_PCAdf.RData")

## first separate the cluster species lineage and then follow same thing as done earlier. 
supra_sp_filt <- rownames(df_pca_supragingival[df_pca_supragingival$cluster == 2,])
supra_sp_filt <- gsub("_"," ", supra_sp_filt)
supragingival_species_lineage_df2 <- supragingival_species_lineage_df[supragingival_species_lineage_df$species %in% supra_sp_filt,]


## Parse taxonomy
supragingival_tax_obj_cluster2 <- parse_tax_data(
  supragingival_species_lineage_df2,
  class_cols = "lineage",
  class_sep = ";"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_3Supragingival_HAC_metacoder_cluster2.pdf", width = 20, height = 10)

heat_tree(
  supragingival_tax_obj_cluster2,
  node_label = taxon_names,
  # Nodes
  node_size = n_obs,
  node_size_range = c(0.02, 0.04),
  
  # Branches
  edge_size = n_obs,
  edge_size_range = c(0.008, 0.02),
  
  # Labels
  node_label_size = n_obs,
  #node_label_size_range = c(0.02, 0.04),
  node_color = "gray60",
  edge_color = "gray70",
  layout = "davidson-harel",
  repel_labels = TRUE,
  repel_force = 8,
  make_node_legend = FALSE,
  make_edge_legend = FALSE
)

dev.off()


########################
############################ Subgingival
subgingival_species_lineage_df <- data.frame(read.delim("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4subgingival_species_HAC90_lineage.txt", header=FALSE))
colnames(subgingival_species_lineage_df) <- c("tax_id","lineage")

subgingival_species_lineage_df$lineage <- trimws(subgingival_species_lineage_df$lineage)
subgingival_species_lineage_df$lineage <- gsub("\\[|\\]", "", subgingival_species_lineage_df$lineage)
subgingival_species_lineage_df$lineage <- gsub("/", "_", subgingival_species_lineage_df$lineage)
subgingival_species_lineage_df$lineage <- sub("^cellular organisms;", "", subgingival_species_lineage_df$lineage)

## Add species name from lineage 
subgingival_species_lineage_df$species <- sapply(strsplit(subgingival_species_lineage_df$lineage, ";"), tail, 1)

## Parse taxonomy
subgingival_tax_obj <- parse_tax_data(
  subgingival_species_lineage_df,
  class_cols = "lineage",
  class_sep = ";"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_3Subgingival_HAC_metacoder_All.pdf", width = 22, height = 12)

heat_tree(
  subgingival_tax_obj,
  node_label = taxon_names,
  # Nodes
  node_size = n_obs,
  node_size_range = c(0.02, 0.04),
  
  # Branches
  edge_size = n_obs,
  edge_size_range = c(0.008, 0.015),
  
  # Labels
  node_label_size = n_obs,
  #node_label_size_range = c(0.02, 0.04),
  node_color = "gray60",
  edge_color = "gray70",
  layout = "kamada-kawai",
  repel_labels = TRUE,
  repel_force = 4,
  make_node_legend = FALSE,
  make_edge_legend = FALSE
)

dev.off()

####### Now plot the same thing for cluster 2 species only.
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_PCAdf.RData")

## first separate the cluster species lineage and then follow same thing as done earlier. 
subgingival_species_lineage_df2 <- subgingival_species_lineage_df[subgingival_species_lineage_df$species %in% gsub("_"," ",rownames(df_pca_subgingival[df_pca_subgingival$cluster == 1,])),]


## Parse taxonomy
subgingival_tax_obj_cluster2 <- parse_tax_data(
  subgingival_species_lineage_df2,
  class_cols = "lineage",
  class_sep = ";"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_3Subgingival_HAC_metacoder_cluster1.pdf", width = 20, height = 10)

heat_tree(
  subgingival_tax_obj_cluster2,
  node_label = taxon_names,
  # Nodes
  node_size = n_obs,
  node_size_range = c(0.02, 0.04),
  
  # Branches
  edge_size = n_obs,
  edge_size_range = c(0.008, 0.02),
  
  # Labels
  node_label_size = n_obs,
  #node_label_size_range = c(0.02, 0.04),
  node_color = "gray60",
  edge_color = "gray70",
  layout = "davidson-harel",
  repel_labels = TRUE,
  repel_force = 8,
  make_node_legend = FALSE,
  make_edge_legend = FALSE
)

dev.off()







########################
############################ tongue
tongue_species_lineage_df <- data.frame(read.delim("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4tongue_species_HAC90_lineage.txt", header=FALSE))
colnames(tongue_species_lineage_df) <- c("tax_id","lineage")

tongue_species_lineage_df$lineage <- trimws(tongue_species_lineage_df$lineage)
tongue_species_lineage_df$lineage <- gsub("\\[|\\]", "", tongue_species_lineage_df$lineage)
tongue_species_lineage_df$lineage <- gsub("/", "_", tongue_species_lineage_df$lineage)
tongue_species_lineage_df$lineage <- sub("^cellular organisms;", "", tongue_species_lineage_df$lineage)

## Add species name from lineage 
tongue_species_lineage_df$species <- sapply(strsplit(tongue_species_lineage_df$lineage, ";"), tail, 1)

## Parse taxonomy
tongue_tax_obj <- parse_tax_data(
  tongue_species_lineage_df,
  class_cols = "lineage",
  class_sep = ";"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_3Tongue_HAC_metacoder_All.pdf", width = 22, height = 12)

heat_tree(
  tongue_tax_obj,
  node_label = taxon_names,
  # Nodes
  node_size = n_obs,
  node_size_range = c(0.02, 0.04),
  
  # Branches
  edge_size = n_obs,
  edge_size_range = c(0.008, 0.015),
  
  # Labels
  node_label_size = n_obs,
  #node_label_size_range = c(0.02, 0.04),
  node_color = "gray60",
  edge_color = "gray70",
  layout = "kamada-kawai",
  repel_labels = TRUE,
  repel_force = 4,
  make_node_legend = FALSE,
  make_edge_legend = FALSE
)

dev.off()

####### Now plot the same thing for cluster 6 species only.
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Tongue_PCAdf.RData")

## first separate the cluster species lineage and then follow same thing as done earlier. 
tongue_species_lineage_df2 <- tongue_species_lineage_df[tongue_species_lineage_df$species %in% gsub("_"," ",rownames(df_pca_tongue[df_pca_tongue$cluster == 6,])),]


## Parse taxonomy
tongue_tax_obj_cluster2 <- parse_tax_data(
  tongue_species_lineage_df2,
  class_cols = "lineage",
  class_sep = ";"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_3Tongue_HAC_metacoder_cluster6.pdf", width = 20, height = 10)

heat_tree(
  tongue_tax_obj_cluster2,
  node_label = taxon_names,
  # Nodes
  node_size = n_obs,
  node_size_range = c(0.02, 0.04),
  
  # Branches
  edge_size = n_obs,
  edge_size_range = c(0.008, 0.02),
  
  # Labels
  node_label_size = n_obs,
  #node_label_size_range = c(0.02, 0.04),
  node_color = "gray60",
  edge_color = "gray70",
  layout = "davidson-harel",
  repel_labels = TRUE,
  repel_force = 8,
  make_node_legend = FALSE,
  make_edge_legend = FALSE
)

dev.off()






save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_3Phylogenetic_Trees_HAC_Clusters_Workspace.RData")
