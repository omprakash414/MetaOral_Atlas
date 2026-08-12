############## Here we are computing the Core–Health Concordance Score from Core Association score and Health Association score.
## No stability score is there for subsites other than saliva

library(dplyr)
########### Load all the scores computed from three different Analysis 
## Core Association Score
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1subgingival_CoreAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1supragingival_CoreAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1tongue_tonsil_CoreAssociationScore.RData")


## Health Association Score
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_HealthAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_HealthAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Tongue_tonsil_HealthAssociationScore.RData")




# ############## Buccal Palate Subsite (We don't have health score for Buccal as it has only 2 study-cohorts and 2 diseases)
# buccal_palate_DiseaseAnalysis_HealthScore <- buccal_palate_DiseaseAnalysis_HealthScore[rownames(buccal_palate_CoreAssociationScore), , drop = FALSE]
# buccal_palate_CoreAssociationScore$species <- NULL

# ### Combine  dfs
# Combined_BuccalPalate_Scores <- cbind(buccal_palate_CoreAssociationScore, buccal_palate_DiseaseAnalysis_HealthScore)
# colnames(Combined_BuccalPalate_Scores) <- c("CoreAssociationScore", "HealthAssociationScore")

# ### Caclulate Combined Score

# Combined_BuccalPalate_Scores$CHC_Score <- apply(Combined_BuccalPalate_Scores[,1:2],1,mean) * (1-apply(Combined_BuccalPalate_Scores[,1:2],1,Gini))
# Combined_BuccalPalate_Scores$CHC_Score[is.nan(Combined_BuccalPalate_Scores$CHC_Score)] <- 0

# Combined_BuccalPalate_Scores$CHC_Score_RankScaled <- rank_scale(Combined_BuccalPalate_Scores$CHC_Score)
# ### Order based on the value in decreasing order
# Combined_BuccalPalate_Scores <- Combined_BuccalPalate_Scores[order(Combined_BuccalPalate_Scores$CHC_Score_RankScaled, decreasing = TRUE),]

# save(Combined_BuccalPalate_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2BuccalPalate_CombinedScores.RData")

# ### temp for supplementary
# temp_Combined_BuccalPalate_Scores <- Combined_BuccalPalate_Scores
# temp_Combined_BuccalPalate_Scores$mean_score <- apply(temp_Combined_BuccalPalate_Scores[,1:2],1,mean)
# temp_Combined_BuccalPalate_Scores$reward_Score <- (1-apply(Combined_BuccalPalate_Scores[,1:2],1,Gini))

# write.csv(temp_Combined_BuccalPalate_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2BuccalPalate_scores_supplementary.csv")



############## Subgingival Subsite
subgingival_DiseaseAnalysis_HealthScore <- subgingival_DiseaseAnalysis_HealthScore[rownames(subgingival_CoreAssociationScore), , drop = FALSE]
subgingival_CoreAssociationScore$species <- NULL

### Combine dfs
Combined_Subgingival_Scores <- cbind(subgingival_CoreAssociationScore,subgingival_DiseaseAnalysis_HealthScore)
colnames(Combined_Subgingival_Scores) <- c("CoreAssociationScore", "HealthAssociationScore")

### Calculate Combined Score
library(LaplacesDemon)
library(DescTools)

rank_scale=function(x){
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

Combined_Subgingival_Scores$CHC_Score <- apply(Combined_Subgingival_Scores[,1:2], 1, mean) * (1 - apply(Combined_Subgingival_Scores[,1:2], 1, Gini))
Combined_Subgingival_Scores$CHC_Score[is.nan(Combined_Subgingival_Scores$CHC_Score)] <- 0

Combined_Subgingival_Scores$CHC_Score_RankScaled <- rank_scale(Combined_Subgingival_Scores$CHC_Score)

### Order
Combined_Subgingival_Scores <- Combined_Subgingival_Scores[order(Combined_Subgingival_Scores$CHC_Score_RankScaled, decreasing = TRUE), ]

save(Combined_Subgingival_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2Subgingival_CombinedScores.RData")


### temp for supplementary
temp_Combined_Subgingival_Scores <- Combined_Subgingival_Scores
temp_Combined_Subgingival_Scores$mean_score <- apply(temp_Combined_Subgingival_Scores[,1:2],1,mean)
temp_Combined_Subgingival_Scores$reward_Score <- (1-apply(temp_Combined_Subgingival_Scores[,1:2],1,Gini))

write.csv(temp_Combined_Subgingival_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2Subgingival_scores_supplementary.csv")


############## Supragingival Subsite
supragingival_DiseaseAnalysis_HealthScore <- supragingival_DiseaseAnalysis_HealthScore[rownames(supragingival_CoreAssociationScore), , drop = FALSE]
supragingival_CoreAssociationScore$species <- NULL

### Combine dfs
Combined_Supragingival_Scores <- cbind(supragingival_CoreAssociationScore,supragingival_DiseaseAnalysis_HealthScore)
colnames(Combined_Supragingival_Scores) <- c("CoreAssociationScore", "HealthAssociationScore")

### Calculate Combined Score
Combined_Supragingival_Scores$CHC_Score <- apply(Combined_Supragingival_Scores[,1:2], 1, mean) * (1 - apply(Combined_Supragingival_Scores[,1:2], 1, Gini))
Combined_Supragingival_Scores$CHC_Score[is.nan(Combined_Supragingival_Scores$CHC_Score)] <- 0

Combined_Supragingival_Scores$CHC_Score_RankScaled <- rank_scale(Combined_Supragingival_Scores$CHC_Score)

### Order
Combined_Supragingival_Scores <- Combined_Supragingival_Scores[order(Combined_Supragingival_Scores$CHC_Score_RankScaled, decreasing = TRUE), ]

save(Combined_Supragingival_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2Supragingival_CombinedScores.RData")


### temp for supplementary
temp_Combined_Supragingival_Scores <- Combined_Supragingival_Scores
temp_Combined_Supragingival_Scores$mean_score <- apply(temp_Combined_Supragingival_Scores[,1:2],1,mean)
temp_Combined_Supragingival_Scores$reward_Score <- (1-apply(temp_Combined_Supragingival_Scores[,1:2],1,Gini))

write.csv(temp_Combined_Supragingival_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2Supragingival_scores_supplementary.csv")




############## Tongue–Tonsil Subsite
tongue_tonsil_DiseaseAnalysis_HealthScore <- tongue_tonsil_DiseaseAnalysis_HealthScore[rownames(tongue_tonsil_CoreAssociationScore), , drop = FALSE]
tongue_tonsil_CoreAssociationScore$species <- NULL

### Combine dfs
Combined_TongueTonsil_Scores <- cbind(tongue_tonsil_CoreAssociationScore,tongue_tonsil_DiseaseAnalysis_HealthScore)
colnames(Combined_TongueTonsil_Scores) <- c("CoreAssociationScore", "HealthAssociationScore")

### Calculate Combined Score
Combined_TongueTonsil_Scores$CHC_Score <- apply(Combined_TongueTonsil_Scores[,1:2], 1, mean) * (1 - apply(Combined_TongueTonsil_Scores[,1:2], 1, Gini))
Combined_TongueTonsil_Scores$CHC_Score[is.nan(Combined_TongueTonsil_Scores$CHC_Score)] <- 0

Combined_TongueTonsil_Scores$CHC_Score_RankScaled <- rank_scale(Combined_TongueTonsil_Scores$CHC_Score)

### Order
Combined_TongueTonsil_Scores <- Combined_TongueTonsil_Scores[order(Combined_TongueTonsil_Scores$CHC_Score_RankScaled, decreasing = TRUE), ]

save(Combined_TongueTonsil_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2TongueTonsil_CombinedScores.RData")


### temp for supplementary
temp_Combined_TongueTonsil_Scores <- Combined_TongueTonsil_Scores
temp_Combined_TongueTonsil_Scores$mean_score <- apply(temp_Combined_TongueTonsil_Scores[,1:2],1,mean)
temp_Combined_TongueTonsil_Scores$reward_Score <- (1-apply(temp_Combined_TongueTonsil_Scores[,1:2],1,Gini))

write.csv(temp_Combined_TongueTonsil_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2TongueTonsil_scores_supplementary.csv")



################ Compare all the combined scores across the four subsites (specie sin each of the subsites are different, while some are same) and then Plot a heatmap for species showing combine score of >=0.7 in atleast 1 subsite.

# Combine all the combined scores into a single dataframe
# BP <- data.frame(
#   species = rownames(Combined_BuccalPalate_Scores),
#   BuccalPalate = Combined_BuccalPalate_Scores$CHC_Score_RankScaled,
#   row.names = NULL
# )

SG <- data.frame(
  species = rownames(Combined_Subgingival_Scores),
  Subgingival = Combined_Subgingival_Scores$CHC_Score_RankScaled,
  row.names = NULL
)

SPG <- data.frame(
  species = rownames(Combined_Supragingival_Scores),
  Supragingival = Combined_Supragingival_Scores$CHC_Score_RankScaled,
  row.names = NULL
)

TT <- data.frame(
  species = rownames(Combined_TongueTonsil_Scores),
  TongueTonsil = Combined_TongueTonsil_Scores$CHC_Score_RankScaled,
  row.names = NULL
)


Combined_Scores <- SG %>%
  full_join(SPG, by = "species") %>%
  full_join(TT,  by = "species")

rownames(Combined_Scores) <- Combined_Scores$species
Combined_Scores$species <- NULL

# Filter rows having atleast one value >=0.8
Combined_Scores_filt <- Combined_Scores[apply(Combined_Scores, 1, function(x) any(x >= 0.8, na.rm = TRUE)), ]

Combined_Scores_filt2 <- Combined_Scores_filt
# Create a heatmap of the filtered combined scores across subsites
Combined_Scores_filt2[is.na(Combined_Scores_filt2)] <- 0

library(pheatmap)
# generate the heatmap and get the order of species
heatmap_result <- pheatmap(Combined_Scores_filt2, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE)
species_order <- heatmap_result$tree_row$order
subsite_order <- heatmap_result$tree_col$order

# arrange the combined_df according to the order of species
Combined_Scores_filt <- Combined_Scores_filt[species_order, subsite_order]

write.csv(Combined_Scores_filt, "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2OtherSubsites_HAC_SpeciesOveralap.csv")




save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2OtherSubsites_CombinedAssociation_Workspace.RData")

