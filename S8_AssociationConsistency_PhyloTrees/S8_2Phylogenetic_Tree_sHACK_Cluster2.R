
install.packages("metacoder")
# https://github.com/grunwaldlab/metacoder
library(metacoder)
library(dplyr)

species_lineage_df <- data.frame(read.delim("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_sHACK_species90_lineage.txt", header=FALSE))
colnames(species_lineage_df) <- c("tax_id","lineage")

species_lineage_df$lineage <- trimws(species_lineage_df$lineage)
species_lineage_df$lineage <- gsub("\\[|\\]", "", species_lineage_df$lineage)
species_lineage_df$lineage <- gsub("/", "_", species_lineage_df$lineage)
species_lineage_df$lineage <- sub("^cellular organisms;", "", species_lineage_df$lineage)

## Add species name from lineage 
species_lineage_df$species <- sapply(strsplit(species_lineage_df$lineage, ";"), tail, 1)
rownames(species_lineage_df) <- species_lineage_df$species
# species_lineage_df$species[species_lineage_df$species == "Hoylesella shahii"] <- "Prevotella shahii"
# species_lineage_df$species[species_lineage_df$species == "Segatella oulorum"] <- "Prevotella oulorum"
# species_lineage_df$species[species_lineage_df$species == "Schaalia odontolytica"] <- "Actinomyces odontolyticus"
# species_lineage_df$species[species_lineage_df$species == "Pseudoleptotrichia goodfellowii"] <- "Leptotrichia goodfellowii"
# species_lineage_df$species[species_lineage_df$species == "Segatella salivae"] <- "Prevotella salivae"
# species_lineage_df$species[species_lineage_df$species == "Hoylesella loescheii"] <- "Prevotella loescheii"

## Parse taxonomy
tax_obj <- parse_tax_data(
  species_lineage_df,
  class_cols = "lineage",
  class_sep = ";"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_2Saliva_sHACK_metacoder_All.pdf", width = 22, height = 12)

heat_tree(
  tax_obj,
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
  node_color = "gray70",
  edge_color = "gray60",
  layout = "kamada-kawai",
  repel_labels = TRUE,
  repel_force = 4,
  make_node_legend = FALSE,
  make_edge_legend = FALSE
)

dev.off()




###############################
#################### Now plot the same thing for cluster 2 species only.
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_PCAdf.RData")

## first separate the cluster species lineage and then follow same thing as done earlier. 
cluster2_species <- rownames(df_pca_saliva[df_pca_saliva$cluster == 2, ])
cluster2_species <- gsub("_", " ", cluster2_species)

species_lineage_df2 <- species_lineage_df[
  species_lineage_df$species %in% cluster2_species,
]

## Parse taxonomy
tax_obj_cluster2 <- parse_tax_data(
  species_lineage_df2,
  class_cols = "lineage",
  class_sep = ";"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_2Saliva_sHACK_metacoder_cluster2.pdf", width = 20, height = 10)

heat_tree(
  tax_obj_cluster2,
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
  node_color = "gray70",
  edge_color = "gray60",
  layout = "davidson-harel",
  repel_labels = TRUE,
  repel_force = 8,
  make_node_legend = FALSE,
  make_edge_legend = FALSE
)

dev.off()


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_2Phylogenetic_Trees_sHACK_Cluster2_Workspace.RData")
