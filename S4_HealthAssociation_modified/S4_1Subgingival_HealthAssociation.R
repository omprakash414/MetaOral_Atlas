

############### Subgingival Health Association Analysis ###############
############ Load the control-disease data for each subsite

#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_BuccalPalate_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Supragingival_DisCtrlCohort.RData")
#load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_TongueTonsil_DisCtrlCohort.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Subgingival_DisCtrlCohort.RData")

# load the selected species for subgingival subsite
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")


################## subgingival

##### Prepare the data 
## get the study specific summary (80-20)
subgingival_summary_df <- make_summary_single_subsite(MetadataDf_subgingival_DisCtrl)
subgingival_summary_df2 <- subgingival_summary_df$summary_raw
subgingival_summary_df3 <- subgingival_summary_df$summary_filtered 


## Get the species profile for Subgingival subsite
SpDf_subgingival_DisCtrl <- SpDf_subgingival_DisCtrl[,colnames(SpDf_subgingival_DisCtrl)%in% subgingival_AssociatedSpecies]
SpDf_subgingival_DisCtrl <- SpDf_subgingival_DisCtrl/rowSums(SpDf_subgingival_DisCtrl)

## get species profile and metadata for only studies that are balanced with control-disease samples.
MetadataDf_subgingival_DisCtrl <- MetadataDf_subgingival_DisCtrl[MetadataDf_subgingival_DisCtrl$study_name %in% subgingival_summary_df2$study_name,]
SpDf_subgingival_DisCtrl <- SpDf_subgingival_DisCtrl[rownames(MetadataDf_subgingival_DisCtrl),]
AllControlSamples <- rownames(MetadataDf_subgingival_DisCtrl[MetadataDf_subgingival_DisCtrl$study_condition == "Control",])
AllDiseaseSamples <- rownames(MetadataDf_subgingival_DisCtrl[MetadataDf_subgingival_DisCtrl$study_condition != "Control",])

selected_studies <- subgingival_summary_df2$study_name
##### Calculate the health association score for each species in Subgingival subsite

subgingival_DiseaseAnalysis_iterative <- healthAssociation_iterations(10,0.65,selected_studies,AllControlSamples,AllDiseaseSamples,MetadataDf_subgingival_DisCtrl,SpDf_subgingival_DisCtrl,subgingival_AssociatedSpecies)
subgingival_DiseaseAnalysis_HealthScore <- subgingival_DiseaseAnalysis_iterative$HealthAssociationScore



save(subgingival_DiseaseAnalysis_HealthScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_HealthAssociationScore.RData")



############
################ Now compare the Core Association Score and Health Association Score for subgingival subsite and then plot it.
##### Import the Core Association Score and then Compare it with Health Association Score
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1subgingival_CoreAssociationScore.RData")
subgingival_CoreAssociationScore <- subgingival_CoreAssociationScore[rownames(subgingival_DiseaseAnalysis_HealthScore),]


subgingival_CoreAS_HealthAS <- bind_cols(subgingival_DiseaseAnalysis_HealthScore, subgingival_CoreAssociationScore)

subgingival_CoreAS_HealthAS$quadrant <- with(
  subgingival_CoreAS_HealthAS,
  ifelse(HealthAssociationScore >= 0.7 & CoreAssociationScore >= 0.7, "Q1",
         ifelse(HealthAssociationScore < 0.7 & CoreAssociationScore >= 0.7, "Q2",
                ifelse(HealthAssociationScore < 0.7 & CoreAssociationScore < 0.7, "Q3", "Q4")))
)


library(dplyr)
library(ggplot2)
library(ggrepel)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_CoreAssociation_HealthAssociation.pdf", width = 30, height = 18)
ggplot(subgingival_CoreAS_HealthAS, aes(x = HealthAssociationScore, y = CoreAssociationScore)) +
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

write.csv(subgingival_CoreAS_HealthAS, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_HealthAssociation_CoreAssociation_Scores.csv")


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_HealthAssociation_Workspace.RData")


############### Now run for all the studies. No iterations. This will be further used for Study specific association of species. i.e species to study (i.e disease specific)
subgingival_DiseaseAnalysis_single <- healthAssociation_iterations(1,1,selected_studies,AllControlSamples,AllDiseaseSamples,MetadataDf_subgingival_DisCtrl,SpDf_subgingival_DisCtrl,subgingival_AssociatedSpecies)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_HealthAssociation_Workspace.RData")



### save study v/s species association df 
subgingival_StudyWise_SpeciesAssociation <- subgingival_DiseaseAnalysis_single$df_comparison_last[[1]]

write.csv(subgingival_StudyWise_SpeciesAssociation, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_StudyWise_SpeciesAssociation.csv", row.names = T)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_HealthAssociation_Workspace.RData")



############# Now plot the scatterplot of Health vs Core score (Score >=0.60 in both)
subgingival_CoreAS_HealthAS2 <- subgingival_CoreAS_HealthAS[subgingival_CoreAS_HealthAS$HealthAssociationScore >= 0.6 & subgingival_CoreAS_HealthAS$CoreAssociationScore >= 0.6,]
library(dplyr)
library(ggplot2)
library(ggrepel)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_CoreAssociation_HealthAssociation_Atleast_60Score.pdf", width = 20, height = 20)
ggplot(subgingival_CoreAS_HealthAS2, aes(x = HealthAssociationScore, y = CoreAssociationScore)) +
  geom_point(color = "black", size = 5) +  # keep dots black
  geom_text_repel(aes(label = rownames(subgingival_CoreAS_HealthAS2), color = quadrant),
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

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_HealthAssociation_Workspace.RData")
