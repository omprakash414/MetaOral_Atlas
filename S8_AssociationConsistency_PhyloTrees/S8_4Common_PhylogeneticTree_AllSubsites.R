

S8_3Phylogenetic_Trees_HAC_Clusters_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_3Phylogenetic_Trees_HAC_Clusters_Workspace.RData", envir = S8_3Phylogenetic_Trees_HAC_Clusters_Workspace)
attach(S8_3Phylogenetic_Trees_HAC_Clusters_Workspace)
supragingival_species_lineage_df <- supragingival_species_lineage_df
supragingival_species_lineage_df2 <- supragingival_species_lineage_df2
df_pca_supragingival <- df_pca_supragingival
subgingival_species_lineage_df <- subgingival_species_lineage_df
subgingival_species_lineage_df2 <- subgingival_species_lineage_df2
df_pca_subgingival <- df_pca_subgingival
tongue_species_lineage_df <- tongue_species_lineage_df
tongue_species_lineage_df2 <- tongue_species_lineage_df2
df_pca_tongue<- df_pca_tongue
detach(S8_3Phylogenetic_Trees_HAC_Clusters_Workspace)
rm(S8_3Phylogenetic_Trees_HAC_Clusters_Workspace)


S8_2Phylogenetic_Trees_sHACK_Cluster2_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_2Phylogenetic_Trees_sHACK_Cluster2_Workspace.RData", envir = S8_2Phylogenetic_Trees_sHACK_Cluster2_Workspace)
attach(S8_2Phylogenetic_Trees_sHACK_Cluster2_Workspace)
species_lineage_df <- species_lineage_df
species_lineage_df2 <- species_lineage_df2
df_pca_saliva <- df_pca_saliva
detach(S8_2Phylogenetic_Trees_sHACK_Cluster2_Workspace)
rm(S8_2Phylogenetic_Trees_sHACK_Cluster2_Workspace)

library(metacoder)
library(dplyr)


### Add subsite column to each lineage df
species_lineage_df$subsite <- "saliva"
supragingival_species_lineage_df$subsite <- "supragingival"
subgingival_species_lineage_df$subsite <- "subgingival"
tongue_species_lineage_df$subsite <- "tongue"

### Now combine all the lineage dfs into one and then only take unique rows
combined_lineage_df <- bind_rows(species_lineage_df,supragingival_species_lineage_df,subgingival_species_lineage_df,tongue_species_lineage_df)
# Remove duplicate rows
combined_lineage_df <- distinct(combined_lineage_df)

### Create metacoder object from combined lineage df
common_tax_obj <- parse_tax_data(
  combined_lineage_df,
  class_cols = "lineage",
  class_sep = ";"
)


### Add subsite information to the metacoder object
# Extract the taxon IDs and names from the metacoder object
taxon_ids_vec <- taxon_ids(common_tax_obj)
taxon_names_vec <- taxon_names(common_tax_obj)

# Make lineage long table: one row per taxon name per species per subsite
lineage_long_df <- combined_lineage_df %>%
  dplyr::select(species, subsite, lineage) %>%
  tidyr::separate_rows(lineage, sep = ";") %>%
  dplyr::rename(taxon_name = lineage) %>%
  dplyr::mutate(taxon_name = trimws(taxon_name))

# For each node/taxon, find subsite combination
taxon_combo <- sapply(taxon_names_vec, function(tx_name) {
  print(paste("Finding subsites for taxon:", tx_name))

  subsites_found <- lineage_long_df$subsite[lineage_long_df$taxon_name == tx_name]
  paste(sort(unique(subsites_found)), collapse = "+")
})

# Number of subsites
taxon_n_subsites <- sapply(strsplit(taxon_combo, "\\+"), length)

# If empty, set to 0
taxon_n_subsites[taxon_combo == ""] <- 0

# Create taxon-level data directly in the metacoder object
common_tax_obj$data$taxon_info <- data.frame(
  taxon_id = taxon_ids_vec,
  taxon_name = taxon_names_vec,
  subsite_combination = taxon_combo,
  n_subsites = taxon_n_subsites,
  dummy = 1,
  stringsAsFactors = FALSE
)

# Color by exact subsite combination
comb_levels <- sort(unique(taxon_combo))
comb_levels
comb_cols <- c(
  "saliva" = "#1B9E77",
  "subgingival" = "#D95F02",
  "supragingival" = "#7570B3",
  "tongue" = "#E7298A",
  "saliva+subgingival" = "#66A61E",
  "saliva+supragingival" = "#E6AB02",
  "saliva+tongue" = "#A6761D",
  "subgingival+supragingival" = "#666666",
  "subgingival+supragingival+tongue" = "#1F78B4",
  "saliva+supragingival+tongue" = "#B2DF8A",
  "saliva+subgingival+supragingival" = "#FB9A99",
  "saliva+subgingival+supragingival+tongue" = "#E31A1C"
)

# Add node color to the metacoder object based on subsite combination as above
common_tax_obj$data$taxon_info$node_col <-
  comb_cols[common_tax_obj$data$taxon_info$subsite_combination]


# check if any combination is missing color
setdiff(unique(common_tax_obj$data$taxon_info$subsite_combination),names(comb_cols))

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_4Common_SubsiteCombination_metacoder.pdf",width = 60,height = 60)
heat_tree(
  common_tax_obj,
  node_label = taxon_names,

  node_size = n_obs,
  node_size_range = c(0.015, 0.035),
  
  edge_size = n_obs,
  edge_size_range = c(0.002, 0.04),
  
  node_color = node_col,
  edge_color = node_col,
  
  node_label_size = n_obs,
  node_label_size_range = c(0.014, 0.04),
  
  layout = "kamada-kawai",
  repel_labels = TRUE,
  repel_force = 1,
  make_node_legend = F,
  make_edge_legend = F
)

dev.off()






#############################
#################### Now get the legends for the subsite combinations and plot them separately
library(ggplot2)
library(dplyr)

legend_df <- data.frame(
  combination = names(comb_cols),
  color = unname(comb_cols),
  stringsAsFactors = FALSE)

legend_df$n_subsites <- sapply(strsplit(legend_df$combination, "\\+"), length)

legend_df <- legend_df %>%
  arrange(n_subsites, combination)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_4Common_SubsiteCombination_legend.pdf",width = 10,height = 6)
ggplot(legend_df, aes(x = 1, y = reorder(combination, -n_subsites))) +
  geom_tile(aes(fill = combination), width = 0.25, height = 0.7) +
  geom_text(aes(x = 1.2, label = combination), hjust = 0, size = 6) +
  scale_fill_manual(values = comb_cols) +
  xlim(0.8, 3.8) +
  theme_void() +
  theme(legend.position = "none")

dev.off()


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_4Common_PhylogeneticTree_AllSubsites_Workspace.RData")









##################################
######################## Now here we will plot the phylogenetic tree for cluster species only.
######################## Follow same exact steps as done in above in this script. Take df2 for every subsite
##################################

### adding subsite column to each lineage df
species_lineage_df2$subsite <- "saliva"
supragingival_species_lineage_df2$subsite <- "supragingival"
subgingival_species_lineage_df2$subsite <- "subgingival"
tongue_species_lineage_df2$subsite <- "tongue"

### Now combine all the lineage dfs into one and then only take unique rows
combined_lineage_df2 <- bind_rows(species_lineage_df2,supragingival_species_lineage_df2,subgingival_species_lineage_df2,tongue_species_lineage_df2)
# Remove duplicate rows
combined_lineage_df2 <- distinct(combined_lineage_df2)

### Create metacoder object from combined lineage df2
common_tax_obj2 <- parse_tax_data(
  combined_lineage_df2,
  class_cols = "lineage",
  class_sep = ";"
)

### Add subsite information to the metacoder object
# Extract the taxon IDs and names from the metacoder object
taxon_ids_vec2 <- taxon_ids(common_tax_obj2)
taxon_names_vec2 <- taxon_names(common_tax_obj2)

# Make lineage long table: one row per taxon name per species per subsite
lineage_long_df2 <- combined_lineage_df2 %>%
  dplyr::select(species, subsite, lineage) %>%
  tidyr::separate_rows(lineage, sep = ";") %>%
  dplyr::rename(taxon_name = lineage) %>%
  dplyr::mutate(taxon_name = trimws(taxon_name))

# For each node/taxon, find subsite combination
taxon_combo2 <- sapply(taxon_names_vec2, function(tx_name) {
  print(paste("Finding subsites for taxon:", tx_name))
  subsites_found <- lineage_long_df2$subsite[lineage_long_df2$taxon_name == tx_name]
  paste(sort(unique(subsites_found)), collapse = "+")
})

# Number of subsites
taxon_n_subsites2 <- sapply(strsplit(taxon_combo2, "\\+"), length)

# If empty, set to 0
taxon_n_subsites2[taxon_combo2 == ""] <- 0

# Create taxon-level data directly in the metacoder object
common_tax_obj2$data$taxon_info <- data.frame(
  taxon_id = taxon_ids_vec2,
  taxon_name = taxon_names_vec2,
  subsite_combination = taxon_combo2,
  n_subsites = taxon_n_subsites2,
  dummy = 1,
  stringsAsFactors = FALSE
)

# Color by exact subsite combination
comb_levels2 <- sort(unique(taxon_combo2))  
comb_levels2
# Add node color to the metacoder object based on subsite combination as above
common_tax_obj2$data$taxon_info$node_col <-
  comb_cols[common_tax_obj2$data$taxon_info$subsite_combination]

# check if any combination is missing color
setdiff(unique(common_tax_obj2$data$taxon_info$subsite_combination),names(comb_cols))

# Plot the phylogenetic tree for cluster species only
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_4Common_SubsiteCombination_metacoder_ClusterSpeciesOnly.pdf",width = 60,height = 60)
heat_tree(
  common_tax_obj2,
  node_label = taxon_names,
  node_size = n_obs,
  node_size_range = c(0.015, 0.035),
  edge_size = n_obs,
  edge_size_range = c(0.002, 0.035),
  node_color = node_col,
  edge_color = node_col,
  node_label_size = n_obs,
  node_label_size_range = c(0.014, 0.04),
  layout = "kamada-kawai",
  repel_labels = TRUE,
  repel_force = 1,
  make_node_legend = F,
  make_edge_legend = F
)
dev.off()



save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_4Common_PhylogeneticTree_AllSubsites_Workspace.RData")

