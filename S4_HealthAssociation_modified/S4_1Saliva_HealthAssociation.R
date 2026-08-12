######## Health association score in Saliva subsite (Code is little modified version of the code used in GutHACK paper)


############ Load the control-disease data for each subsite

#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_BuccalPalate_DisCtrlCohort.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Supragingival_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_TongueTonsil_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Subgingival_DisCtrlCohort.RData")

# load the selected species for saliva subsite
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")


################## Saliva

##### Prepare the data 
## get the study specific summary (80-20)
saliva_summary_df <- make_summary_single_subsite(MetadataDf_saliva_DisCtrl)
saliva_summary_df2 <- saliva_summary_df$summary_raw
saliva_summary_df3 <- saliva_summary_df$summary_filtered # Here LinM_2020 has 0.197 and NagataN_2022 has 0.82, this are borderlined so will not exclude from the disease analysis


## Get the species profile for saliva subsite
SpDf_saliva_DisCtrl <- SpDf_saliva_DisCtrl[,colnames(SpDf_saliva_DisCtrl)%in% saliva_AssociatedSpecies]
SpDf_saliva_DisCtrl <- SpDf_saliva_DisCtrl/rowSums(SpDf_saliva_DisCtrl)

## get species profile and metadata for only studies that are balanced with control-disease samples.
MetadataDf_saliva_DisCtrl <- MetadataDf_saliva_DisCtrl[MetadataDf_saliva_DisCtrl$study_name %in% saliva_summary_df2$study_name,]
SpDf_saliva_DisCtrl <- SpDf_saliva_DisCtrl[rownames(MetadataDf_saliva_DisCtrl),]
AllControlSamples <- rownames(MetadataDf_saliva_DisCtrl[MetadataDf_saliva_DisCtrl$study_condition == "Control",])
AllDiseaseSamples <- rownames(MetadataDf_saliva_DisCtrl[MetadataDf_saliva_DisCtrl$study_condition != "Control",])

selected_studies <- saliva_summary_df2$study_name
##### Calculate the health association score for each species in saliva subsite

Saliva_DiseaseAnalysis_iterative <- healthAssociation_iterations(10,0.65,selected_studies,AllControlSamples,AllDiseaseSamples,MetadataDf_saliva_DisCtrl,SpDf_saliva_DisCtrl,saliva_AssociatedSpecies)
Saliva_DiseaseAnalysis_HealthScore <- Saliva_DiseaseAnalysis_iterative$HealthAssociationScore

save(Saliva_DiseaseAnalysis_HealthScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociationScore.RData")


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociation_Workspace.RData")




############### Now run for all the studies. No iterations. This will be further used for Study specific association of species. i.e species to study (i.e disease specific)
Saliva_DiseaseAnalysis_single <- healthAssociation_iterations(1,1,selected_studies,AllControlSamples,AllDiseaseSamples,MetadataDf_saliva_DisCtrl,SpDf_saliva_DisCtrl,saliva_AssociatedSpecies)

Saliva_HealthScore_NonIter <- Saliva_DiseaseAnalysis_single$HealthAssociationScore


save(Saliva_HealthScore_NonIter, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociationScore_NonIter.RData")
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociation_Workspace.RData")


### save study v/s species association df 
saliva_StudyWise_SpeciesAssociation <- Saliva_DiseaseAnalysis_single$df_comparison_last[[1]]

write.csv(saliva_StudyWise_SpeciesAssociation, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_StudyWise_SpeciesAssociation.csv", row.names = T)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociation_Workspace.RData")






###########################
########################### Plot the study wise association of species for saliva subsite as a heatmap. 
#### Before that we want only species that have strong signatures (so we need to reduce the species number to only those that have strong association with health or disease.)
library(dplyr)
saliva_StudyWise_SpeciesAssociation2 <- saliva_StudyWise_SpeciesAssociation
saliva_StudyWise_SpeciesAssociation2[is.na(saliva_StudyWise_SpeciesAssociation2)] <- 0

# temp_StudyWise_SpeciesAssociation <- data.frame(Negative = rowSums(saliva_StudyWise_SpeciesAssociation2 < -1,na.rm = TRUE),
#                                  Positive = rowSums(saliva_StudyWise_SpeciesAssociation2 > 1, na.rm = TRUE))

# # health association
# df_DiseaseAssociated <- temp_StudyWise_SpeciesAssociation[(temp_StudyWise_SpeciesAssociation$Positive <= 3) & (temp_StudyWise_SpeciesAssociation$Negative >= 8),]
# df_DiseaseAssociated <- df_DiseaseAssociated[order(df_DiseaseAssociated$Positive - df_DiseaseAssociated$Negative),]

# # disease association
# df_HealthAssociated <- temp_StudyWise_SpeciesAssociation[(temp_StudyWise_SpeciesAssociation$Negative <= 3) & (temp_StudyWise_SpeciesAssociation$Positive >= 8),]
# df_HealthAssociated <- df_HealthAssociated[order(df_HealthAssociated$Positive - df_HealthAssociated$Negative),]

# # Add both healthy and disease associated markers.
# temp_StudyWise_SpeciesAssociation <- rbind(df_DiseaseAssociated,df_HealthAssociated)

# Arrange the names of the species based on the health score and then take top 50 and bottom 50 species for plotting. 
Saliva_DiseaseAnalysis_HealthScore <- Saliva_DiseaseAnalysis_HealthScore[order(Saliva_DiseaseAnalysis_HealthScore$HealthAssociationScore, decreasing = TRUE),,drop = FALSE]
saliva_StudyWise_SpeciesAssociation2 <- saliva_StudyWise_SpeciesAssociation2[rownames(Saliva_DiseaseAnalysis_HealthScore),]

# filter the species and Order: disease-associated on top, health-associated bottom
heatmap_filtered <- rbind(head(saliva_StudyWise_SpeciesAssociation2, 50), tail(saliva_StudyWise_SpeciesAssociation2, 50))
# Plot the heatmap
library(pheatmap)
color_map <- c("-3"  = "#8B0000","-2"  = "#CD5C5C","-1"  = "white","0"  = "white","1" = "white","2" = "#329732","3" = "#006400")

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_heatmap_StudyWise_Species.pdf", width = 30, height = 15)
pheatmap(t(heatmap_filtered),
         color = color_map,
         fontsize_row = 8,
         fontsize_col = 8,
         cellheight = 23,
         cellwidth = 20,
         cluster_rows = T,
         cluster_cols = F,
         border_color = "black",
         treeheight_row = 0,
         treeheight_col = 0
)
dev.off()


## save this heatmap carpet as csv after getting the ordered columns (.ie study names ordered as above we are using dendogram for columns)
ph <- pheatmap(
  t(heatmap_filtered),
  cluster_rows = T,
  cluster_cols = F,
  silent = TRUE
) 
ordered_cols <- colnames(heatmap_filtered)[ph$tree_col$order]
heatmap_filtered <- heatmap_filtered[, ordered_cols]

write.csv(heatmap_filtered, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_carpet_heatmap_StudyWise_Species.csv")



### Now plot the dot and line plot showing the number of postitive and negative associations for each of the species. 
mat <- as.matrix(heatmap_filtered)
species_count_df <- data.frame(
  species = rownames(mat),
  Negative_count = rowSums(matrix(mat %in% c(-2, -3), nrow = nrow(mat))),
  Positive_count = rowSums(matrix(mat %in% c(2, 3), nrow = nrow(mat)))
)

species_count_long <- reshape(
  species_count_df,
  varying = c("Negative_count", "Positive_count"),
  v.names = "Count",
  timevar = "AssociationType",
  times = c("Negative", "Positive"),
  direction = "long"
)

species_count_long$species <- factor(
  species_count_long$species,
  levels = rownames(heatmap_filtered)
)

library(ggplot2)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_species_positive_negative_counts_dot_line_plot.pdf", width = 30, height = 8)

ggplot(
  species_count_long,
  aes(x = species, y = Count, group = AssociationType, color = AssociationType)
) +
  geom_line(linewidth = 1.2) +
  geom_point(color = "black", size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    text = element_text(size = 14)
  ) +
  labs(
    x = "Species",
    y = "Number of strong associations",
    color = "Association type"
  )

dev.off()















### Also plot the single column of Health associated score for species that are in heatmap_filtered
heatmap_filtered_HS_df <- Saliva_DiseaseAnalysis_HealthScore[rownames(heatmap_filtered),,drop = FALSE]
write.csv(heatmap_filtered_HS_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_heatmap_StudyWise_Species_HealthScores.csv")




save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociation_Workspace.RData")





############## Get the health scores in each iteration and export as csv
Iterative_comparison_df <- Saliva_DiseaseAnalysis_iterative$df_directions_all
combined_Iterative_comparison_df <- do.call(cbind, Iterative_comparison_df)

write.csv(combined_Iterative_comparison_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociationScore_AllIterations.csv", row.names = TRUE)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociation_Workspace.RData")

############## Now get the ranked score for each iteration and export as csv
Iteration_Ranked_ScoreDf <- Saliva_DiseaseAnalysis_iterative$RankedIteration_scores

write.csv(Iteration_Ranked_ScoreDf, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociationScore_RankedIterations.csv", row.names = TRUE)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociation_Workspace.RData")

############## now get the Mean ranked score (i.e Health score)
Mean_Ranked_Score <- Saliva_DiseaseAnalysis_iterative$HealthAssociationScore

write.csv(Mean_Ranked_Score, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociationScore_MeanRankedScores.csv", row.names = TRUE)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociation_Workspace.RData")
