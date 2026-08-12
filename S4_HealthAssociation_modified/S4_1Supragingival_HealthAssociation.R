
############### Supragingival Health Association Analysis ###############
############ Load the control-disease data for each subsite

#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_BuccalPalate_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_DisCtrlCohort.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Supragingival_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_TongueTonsil_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Subgingival_DisCtrlCohort.RData")

# load the selected species for supragingival subsite
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")


################## supragingival

##### Prepare the data 
## get the study specific summary (80-20)
supragingival_summary_df <- make_summary_single_subsite(MetadataDf_supragingival_DisCtrl)
supragingival_summary_df2 <- supragingival_summary_df$summary_raw
supragingival_summary_df3 <- supragingival_summary_df$summary_filtered 

## Get the species profile for saliva subsite
SpDf_supragingival_DisCtrl <- SpDf_supragingival_DisCtrl[,colnames(SpDf_supragingival_DisCtrl)%in% supragingival_AssociatedSpecies]
SpDf_supragingival_DisCtrl <- SpDf_supragingival_DisCtrl/rowSums(SpDf_supragingival_DisCtrl)

## get species profile and metadata for only studies that are balanced with control-disease samples.
MetadataDf_supragingival_DisCtrl <- MetadataDf_supragingival_DisCtrl[MetadataDf_supragingival_DisCtrl$study_name %in% supragingival_summary_df2$study_name,]
SpDf_supragingival_DisCtrl <- SpDf_supragingival_DisCtrl[rownames(MetadataDf_supragingival_DisCtrl),]
AllControlSamples <- rownames(MetadataDf_supragingival_DisCtrl[MetadataDf_supragingival_DisCtrl$study_condition == "Control",])
AllDiseaseSamples <- rownames(MetadataDf_supragingival_DisCtrl[MetadataDf_supragingival_DisCtrl$study_condition != "Control",])

selected_studies <- supragingival_summary_df2$study_name
##### Calculate the health association score for each species in saliva subsite

supragingival_DiseaseAnalysis_iterative <- healthAssociation_iterations(10,0.65,selected_studies,AllControlSamples,AllDiseaseSamples,MetadataDf_supragingival_DisCtrl,SpDf_supragingival_DisCtrl,supragingival_AssociatedSpecies)
supragingival_DiseaseAnalysis_HealthScore <- supragingival_DiseaseAnalysis_iterative$HealthAssociationScore


save(supragingival_DiseaseAnalysis_HealthScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_HealthAssociationScore.RData")




###########
############### Compare core association score with health association score for supragingival subsite and plot it
##### Import the Core Association Score and then Compare it with Health Association Score
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1supragingival_CoreAssociationScore.RData")
supragingival_CoreAssociationScore <- supragingival_CoreAssociationScore[rownames(supragingival_DiseaseAnalysis_HealthScore),]


supragingival_CoreAS_HealthAS <- bind_cols(supragingival_DiseaseAnalysis_HealthScore, supragingival_CoreAssociationScore)

supragingival_CoreAS_HealthAS$quadrant <- with(
  supragingival_CoreAS_HealthAS,
  ifelse(HealthAssociationScore >= 0.7 & CoreAssociationScore >= 0.7, "Q1",
         ifelse(HealthAssociationScore < 0.7 & CoreAssociationScore >= 0.7, "Q2",
                ifelse(HealthAssociationScore < 0.7 & CoreAssociationScore < 0.7, "Q3", "Q4")))
)


# Plot both the scores to see the concordence
library(dplyr)
library(ggplot2)
library(ggrepel)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_CoreAssociation_HealthAssociation.pdf", width = 30, height = 18)
ggplot(supragingival_CoreAS_HealthAS, aes(x = HealthAssociationScore, y = CoreAssociationScore)) +
  geom_point(color = "black", size = 2) +  # keep dots black
  geom_text_repel(aes(label = species, color = quadrant),
                  size = 7, max.overlaps = 50) +  # color only labels
  geom_vline(xintercept = 0.7, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "solid", color = "black") +
  theme_bw() +
  xlab("Health Association") +
  ylab("Core Association") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  scale_color_manual(values = c("Q1" = "purple", "Q2" = "red", "Q3" = "brown", "Q4" = "blue"))

dev.off()


write.csv(supragingival_CoreAS_HealthAS, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_HealthAssociation_CoreAssociation_Scores.csv")
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_HealthAssociation_Workspace.RData")




############### Now run for all the studies. No iterations. This will be further used for Study specific association of species. i.e species to study (i.e disease specific)
supragingival_DiseaseAnalysis_single <- healthAssociation_iterations(1,1,selected_studies,AllControlSamples,AllDiseaseSamples,MetadataDf_supragingival_DisCtrl,SpDf_supragingival_DisCtrl,supragingival_AssociatedSpecies)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_HealthAssociation_Workspace.RData")


### save study v/s species association df 
supragingival_StudyWise_SpeciesAssociation <- supragingival_DiseaseAnalysis_single$df_comparison_last[[1]]

write.csv(supragingival_StudyWise_SpeciesAssociation, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_StudyWise_SpeciesAssociation.csv", row.names = T)
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_HealthAssociation_Workspace.RData")



############## Scatter plot of HS and CS (Only show score >=0.60) 
supragingival_CoreAS_HealthAS2 <- supragingival_CoreAS_HealthAS[supragingival_CoreAS_HealthAS$HealthAssociationScore >= 0.60 & supragingival_CoreAS_HealthAS$CoreAssociationScore >= 0.60, ]

library(dplyr)
library(ggplot2)
library(ggrepel)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_CoreAssociation_HealthAssociation_Atleast_60Score.pdf", width = 20, height = 20)
ggplot(supragingival_CoreAS_HealthAS2, aes(x = HealthAssociationScore, y = CoreAssociationScore)) +
  geom_point(color = "black", size = 5) +  # keep dots black
  geom_text_repel(aes(label = rownames(supragingival_CoreAS_HealthAS2), color = quadrant),
                  size = 11, max.overlaps = 60) +  # color only labels
  geom_vline(xintercept = 0.7, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "solid", color = "black") +
  theme_bw() +
  xlab("Health Association") +
  ylab("Core Association") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  scale_color_manual(values = c("Q1" = "purple", "Q2" = "red", "Q3" = "brown", "Q4" = "blue"))

dev.off()


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_HealthAssociation_Workspace.RData")
