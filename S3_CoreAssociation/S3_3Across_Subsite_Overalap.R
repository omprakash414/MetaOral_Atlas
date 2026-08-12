######## Take top 50 core species with high CS from each subsite and generate a common heatmap from it. 

load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_CoreAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1supragingival_CoreAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1subgingival_CoreAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1buccal_palate_CoreAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1tongue_tonsil_CoreAssociationScore.RData")

# # Separate the top 50 species with high CS having value greater than 0.7 from each subsite  
# saliva_top <- saliva_CoreAssociationScore[saliva_CoreAssociationScore$CoreAssociationScore > 0.7,]
# supragingival_top <- supragingival_CoreAssociationScore[supragingival_CoreAssociationScore$CoreAssociationScore > 0.7,]
# subgingival_top <- subgingival_CoreAssociationScore[subgingival_CoreAssociationScore$CoreAssociationScore > 0.7,]
# buccal_palate_top <- buccal_palate_CoreAssociationScore[buccal_palate_CoreAssociationScore$CoreAssociationScore > 0.7,]
# tongue_tonsil_top <- tongue_tonsil_CoreAssociationScore[tongue_tonsil_CoreAssociationScore$CoreAssociationScore > 0.7,]

# Seaprate the top 50 species with high CS from each subsite
saliva_top50 <- saliva_CoreAssociationScore[order(saliva_CoreAssociationScore$CoreAssociationScore, decreasing = TRUE),][1:50,]
supragingival_top50 <- supragingival_CoreAssociationScore[order(supragingival_CoreAssociationScore$CoreAssociationScore, decreasing = TRUE),][1:50,]
subgingival_top50 <- subgingival_CoreAssociationScore[order(subgingival_CoreAssociationScore$CoreAssociationScore, decreasing = TRUE),][1:50,]
buccal_palate_top50 <- buccal_palate_CoreAssociationScore[order(buccal_palate_CoreAssociationScore$CoreAssociationScore, decreasing = TRUE),][1:50,]
tongue_tonsil_top50 <- tongue_tonsil_CoreAssociationScore[order(tongue_tonsil_CoreAssociationScore$CoreAssociationScore, decreasing = TRUE),][1:50,]


######  combine these dfs, but each of this df have specie snames in rownames and also in the column species. Also now all species are same in all dfs. So, we need to match the speceis name and add accordingly.
# create a list of all unique species from the top 50 lists
all_species <- unique(c(rownames(saliva_top50), rownames(supragingival_top50), rownames(subgingival_top50), rownames(buccal_palate_top50), rownames(tongue_tonsil_top50)))
# create an empty dataframe to store the combined data
combined_df <- data.frame(matrix(NA, nrow = length(all_species), ncol = 5))
rownames(combined_df) <- all_species
colnames(combined_df) <- c("Saliva", "Supragingival", "Subgingival", "Buccal_Palate", "Tongue_Tonsil")

## Now fill the df, some specie snames might be same across dfs, so scores of all dfs for that species should also be added.
for (species in all_species) {
  if (species %in% rownames(saliva_top50)) {
    combined_df[species, "Saliva"] <- saliva_top50[species, "CoreAssociationScore"]
  }
  if (species %in% rownames(supragingival_top50)) {
    combined_df[species, "Supragingival"] <- supragingival_top50[species, "CoreAssociationScore"]
  }
  if (species %in% rownames(subgingival_top50)) {
    combined_df[species, "Subgingival"] <- subgingival_top50[species, "CoreAssociationScore"]
  }
  if (species %in% rownames(buccal_palate_top50)) {
    combined_df[species, "Buccal_Palate"] <- buccal_palate_top50[species, "CoreAssociationScore"]
  }
  if (species %in% rownames(tongue_tonsil_top50)) {
    combined_df[species, "Tongue_Tonsil"] <- tongue_tonsil_top50[species, "CoreAssociationScore"]
  }
}

## First replace NAs with 0 and then Now get the order of species after they are clustered in the heatmap, and then arrange the df according to that order.
combined_df[is.na(combined_df)] <- 0

library(pheatmap)
# generate the heatmap and get the order of species
heatmap_result <- pheatmap(combined_df, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE)
species_order <- heatmap_result$tree_row$order
subsite_order <- heatmap_result$tree_col$order

# arrange the combined_df according to the order of species
combined_df <- combined_df[species_order, subsite_order]

write.csv(combined_df, "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_3Across_Subsite_SpeciesOveralap.csv")

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_3Across_Subsite_SpeciesOveralap_Workspace.RData")

