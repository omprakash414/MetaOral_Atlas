
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_DisCtrlCohort.RData")


load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")


source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")


#################################### Do the same things after dividing the data into three cohorts - WGS, 16S, and 16S - Exposure
##### divide the data into three:
## Exposure Studies
positive_exposure <- c("alcoholic", "smoker", "tobacco", "e-cigar")

keep_row <- sapply(MetadataDf_saliva_DisCtrl$exposure, function(x) {
  # split by comma
  tokens <- trimws(unlist(strsplit(tolower(x), "[, ;]")))
  
  # check if ANY positive exposure is present
  any(tokens %in% positive_exposure)
})

MetadataDf_saliva_DisCtrl_16s_exposure <- MetadataDf_saliva_DisCtrl[keep_row, ]

saliva_summary_df <- make_summary_single_subsite(MetadataDf_saliva_DisCtrl_16s_exposure,9,10)
saliva_summary_df2 <- saliva_summary_df$summary_raw
saliva_summary_df3 <- saliva_summary_df$summary_filtered

## FIlter the metadata and species profile wrt balanced studies.
MetadataDf_saliva_DisCtrl_16s_exposure2 <- MetadataDf_saliva_DisCtrl_16s_exposure[MetadataDf_saliva_DisCtrl_16s_exposure$study_name %in% saliva_summary_df3$study_name,]

SpDf_saliva_DisCtrl_16s_exposure<- SpDf_saliva_DisCtrl[rownames(MetadataDf_saliva_DisCtrl_16s_exposure2),colnames(SpDf_saliva_DisCtrl)%in% saliva_AssociatedSpecies]
SpDf_saliva_DisCtrl_16s_exposure <- SpDf_saliva_DisCtrl_16s_exposure/rowSums(SpDf_saliva_DisCtrl_16s_exposure)
SpDf_saliva_DisCtrl_16s_exposure <- SpDf_saliva_DisCtrl_16s_exposure[,colSums(SpDf_saliva_DisCtrl_16s_exposure)>0]
saliva_species_16s_exposure <- colnames(SpDf_saliva_DisCtrl_16s_exposure)

# separate the control and disease samples for 16s-Exposure samples, and also the studies.
AllControlSamples_16sExposure <- rownames(MetadataDf_saliva_DisCtrl_16s_exposure2[MetadataDf_saliva_DisCtrl_16s_exposure2$study_condition == "Control",])
AllDiseaseSamples_16sExposure <- rownames(MetadataDf_saliva_DisCtrl_16s_exposure2[MetadataDf_saliva_DisCtrl_16s_exposure2$study_condition != "Control",])
selected_studies_16sExposure <- saliva_summary_df3$study_name

# Run the health association function for 16s-Exposure samples.
saliva_DiseaseAnalysis_iterative_16_exposure <- healthAssociation_iterations(10,0.65,selected_studies_16sExposure,AllControlSamples_16sExposure,AllDiseaseSamples_16sExposure,MetadataDf_saliva_DisCtrl_16s_exposure2,SpDf_saliva_DisCtrl_16s_exposure,saliva_species_16s_exposure)
saliva_HealthScore_16s_exposure <- saliva_DiseaseAnalysis_iterative_16_exposure$HealthAssociation




## 16s Studies
# All the studies used in exposure if extracted from any study then its remianing samples will not be used in 16s as their count will be less than 10 (confirmed deom S0)
MetadataDf_saliva_DisCtrl_16s_only <- MetadataDf_saliva_DisCtrl[MetadataDf_saliva_DisCtrl$seq_type == "16s", ]
MetadataDf_saliva_DisCtrl_16s_only <- MetadataDf_saliva_DisCtrl_16s_only[! rownames(MetadataDf_saliva_DisCtrl_16s_only) %in% rownames(MetadataDf_saliva_DisCtrl_16s_exposure), ]
# [1] 3952   30

saliva_summary_df4 <- make_summary_single_subsite(MetadataDf_saliva_DisCtrl_16s_only)
saliva_summary_df5 <- saliva_summary_df4$summary_raw
saliva_summary_df6 <- saliva_summary_df4$summary_filtered

MetadataDf_saliva_DisCtrl_16s_only <- MetadataDf_saliva_DisCtrl_16s_only[MetadataDf_saliva_DisCtrl_16s_only$study_name %in% saliva_summary_df5$study_name,]
# [1] 3839   30

## Get the species profile for supragingival subsite
SpDf_saliva_DisCtrl_16s_only<- SpDf_saliva_DisCtrl[rownames(MetadataDf_saliva_DisCtrl_16s_only),colnames(SpDf_saliva_DisCtrl)%in% saliva_AssociatedSpecies]
SpDf_saliva_DisCtrl_16s_only <- SpDf_saliva_DisCtrl_16s_only[,colSums(SpDf_saliva_DisCtrl_16s_only)>0]
saliva_species_16s_only <- colnames(SpDf_saliva_DisCtrl_16s_only)

SpDf_saliva_DisCtrl_16s_only <- SpDf_saliva_DisCtrl_16s_only/rowSums(SpDf_saliva_DisCtrl_16s_only)

# separate the control and disease samples for 16s-Exposure samples, and also the studies.
AllControlSamples_16s_only <- rownames(MetadataDf_saliva_DisCtrl_16s_only[MetadataDf_saliva_DisCtrl_16s_only$study_condition == "Control",])
AllDiseaseSamples_16s_only <- rownames(MetadataDf_saliva_DisCtrl_16s_only[MetadataDf_saliva_DisCtrl_16s_only$study_condition != "Control",])
selected_studies_16s_only <- saliva_summary_df5$study_name

# Run the health association function for 16s-Exposure samples.
saliva_DiseaseAnalysis_iterative_16_only <- healthAssociation_iterations(10,0.65,selected_studies_16s_only,AllControlSamples_16s_only,AllDiseaseSamples_16s_only,MetadataDf_saliva_DisCtrl_16s_only,SpDf_saliva_DisCtrl_16s_only,saliva_species_16s_only)
saliva_HealthScore_16s_only <- saliva_DiseaseAnalysis_iterative_16_only$HealthAssociation





## WGS Studies
MetadataDf_saliva_DisCtrl_WGS <- MetadataDf_saliva_DisCtrl[MetadataDf_saliva_DisCtrl$seq_type == "WGS", ]
  
saliva_summary_df7 <- make_summary_single_subsite(MetadataDf_saliva_DisCtrl_WGS)
saliva_summary_df8 <- saliva_summary_df7$summary_raw
saliva_summary_df9 <- saliva_summary_df7$summary_filtered

unique(MetadataDf_saliva_DisCtrl_WGS$study_name)
unique(saliva_summary_df8$study_name)

MetadataDf_saliva_DisCtrl_WGS <- MetadataDf_saliva_DisCtrl_WGS[MetadataDf_saliva_DisCtrl_WGS$study_name %in% saliva_summary_df8$study_name,]

## Get the species profile for supragingival subsite
SpDf_saliva_DisCtrl_WGS<- SpDf_saliva_DisCtrl[rownames(MetadataDf_saliva_DisCtrl_WGS),colnames(SpDf_saliva_DisCtrl)%in% saliva_AssociatedSpecies]
SpDf_saliva_DisCtrl_WGS <- SpDf_saliva_DisCtrl_WGS[,colSums(SpDf_saliva_DisCtrl_WGS)>0]
saliva_species_WGS <- colnames(SpDf_saliva_DisCtrl_WGS) 

SpDf_saliva_DisCtrl_WGS <- SpDf_saliva_DisCtrl_WGS/rowSums(SpDf_saliva_DisCtrl_WGS)

# separate the control and disease samples for 16s-Exposure samples, and also the studies.
AllControlSamples_WGS <- rownames(MetadataDf_saliva_DisCtrl_WGS[MetadataDf_saliva_DisCtrl_WGS$study_condition == "Control",])
AllDiseaseSamples_WGS <- rownames(MetadataDf_saliva_DisCtrl_WGS[MetadataDf_saliva_DisCtrl_WGS$study_condition != "Control",])
selected_studies_WGS <- saliva_summary_df8$study_name

# Run the health association function for 16s-Exposure samples.
saliva_DiseaseAnalysis_iterative_WGS <- healthAssociation_iterations(10,0.65,selected_studies_WGS,AllControlSamples_WGS,AllDiseaseSamples_WGS,MetadataDf_saliva_DisCtrl_WGS,SpDf_saliva_DisCtrl_WGS,saliva_species_WGS)
saliva_HealthScore_WGS <- saliva_DiseaseAnalysis_iterative_WGS$HealthAssociation    



# ## Non-Exposure Studies
# MetadataDf_saliva_DisCtrl_nonExposure <- MetadataDf_saliva_DisCtrl[rownames(MetadataDf_saliva_DisCtrl) %in% c(rownames(MetadataDf_saliva_DisCtrl_WGS),rownames(MetadataDf_saliva_DisCtrl_16s_only)), ]

# saliva_summary_df10 <- make_summary_single_subsite(MetadataDf_saliva_DisCtrl_nonExposure)
# saliva_summary_df11 <- saliva_summary_df10$summary_raw
# saliva_summary_df12 <- saliva_summary_df10$summary_filtered

# unique(MetadataDf_saliva_DisCtrl_nonExposure$study_name)
# unique(saliva_summary_df11$study_name)

# MetadataDf_saliva_DisCtrl_nonExposure <- MetadataDf_saliva_DisCtrl_nonExposure[MetadataDf_saliva_DisCtrl_nonExposure$study_name %in% saliva_summary_df11$study_name,]

# ## Get the species profile for supragingival subsite
# SpDf_saliva_DisCtrl_nonExposure<- SpDf_saliva_DisCtrl[rownames(MetadataDf_saliva_DisCtrl_nonExposure),colnames(SpDf_saliva_DisCtrl)%in% saliva_AssociatedSpecies]
# SpDf_saliva_DisCtrl_nonExposure <- SpDf_saliva_DisCtrl_nonExposure[,colSums(SpDf_saliva_DisCtrl_nonExposure)>0]
# saliva_species_nonExposure <- colnames(SpDf_saliva_DisCtrl_nonExposure) 

# SpDf_saliva_DisCtrl_nonExposure <- SpDf_saliva_DisCtrl_nonExposure/rowSums(SpDf_saliva_DisCtrl_nonExposure)
# SpDf_saliva_DisCtrl_nonExposure[is.na(SpDf_saliva_DisCtrl_nonExposure)] <- 0

# # separate the control and disease samples for 16s-Exposure samples, and also the studies.
# AllControlSamples_nonExposure <- rownames(MetadataDf_saliva_DisCtrl_nonExposure[MetadataDf_saliva_DisCtrl_nonExposure$study_condition == "Control",])
# AllDiseaseSamples_nonExposure <- rownames(MetadataDf_saliva_DisCtrl_nonExposure[MetadataDf_saliva_DisCtrl_nonExposure$study_condition != "Control",])
# selected_studies_nonExposure <- saliva_summary_df11$study_name

# # Run the health association function for 16s-Exposure samples.
# saliva_DiseaseAnalysis_iterative_nonExposure <- healthAssociation_iterations(10,0.65,selected_studies_nonExposure,AllControlSamples_nonExposure,AllDiseaseSamples_nonExposure,MetadataDf_saliva_DisCtrl_nonExposure,SpDf_saliva_DisCtrl_nonExposure,saliva_species_nonExposure)
# saliva_HealthScore_nonExposure <- saliva_DiseaseAnalysis_iterative_nonExposure$HealthAssociation    







#############
################# Now combine with overall health score computed for saliva subsite and then compare it with cohort specific health scores (16s/wgs/16sExposure)
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociationScore.RData")

Saliva_DiseaseAnalysis_HealthScore <- Saliva_DiseaseAnalysis_HealthScore[order(Saliva_DiseaseAnalysis_HealthScore$HealthAssociationScore, decreasing = TRUE), , drop=FALSE]

# Add the species as a column from rownames
Saliva_DiseaseAnalysis_HealthScore$species <- rownames(Saliva_DiseaseAnalysis_HealthScore)
saliva_HealthScore_16s_exposure$species    <- rownames(saliva_HealthScore_16s_exposure)
saliva_HealthScore_16s_only$species        <- rownames(saliva_HealthScore_16s_only)
saliva_HealthScore_WGS$species             <- rownames(saliva_HealthScore_WGS)

colnames(Saliva_DiseaseAnalysis_HealthScore)[1] <- "HealthAssociationScore"

colnames(saliva_HealthScore_16s_exposure)[1] <- "HealthAssociationScore_16s_exposure"
colnames(saliva_HealthScore_16s_only)[1]     <- "HealthAssociationScore_16s"
colnames(saliva_HealthScore_WGS)[1]          <- "HealthAssociationScore_WGS"

saliva_CohortWise_HealthAssociationScore <-
  Reduce(
    function(x, y) merge(x, y, by = "species", all.x = TRUE),
    list(
      Saliva_DiseaseAnalysis_HealthScore,
      saliva_HealthScore_16s_exposure,
      saliva_HealthScore_16s_only,
      saliva_HealthScore_WGS
    )
  )

rownames(saliva_CohortWise_HealthAssociationScore) <- saliva_CohortWise_HealthAssociationScore$species
saliva_CohortWise_HealthAssociationScore$species <- NULL




#############
############### save

save(saliva_CohortWise_HealthAssociationScore,file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_CohortWise_HealthAssociationScore.RData")
write.csv(saliva_CohortWise_HealthAssociationScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_CohortWise_HealthAssociationScore.csv", row.names = TRUE, quote = FALSE)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_CohortWise_HealthAssociation_Workspace.RData")







#############
################# Count the samples used each of the cohorts for disease association analysis
## Exposure:
unique(saliva_summary_df3$study_name)
sum(saliva_summary_df3$control_count)
sum(saliva_summary_df3$disease_count)
unique(saliva_summary_df3$diseases)

# studies : 5
# total samples: 657
# control count: 287
# disease count: 370
# diseases: 
#   [1] "Pancreatic_Cancer"    "HIV"                 
#   [3] "Head_And_Neck_Cancer" "Oral_SCC_Variants"

## 16s only:
unique(saliva_summary_df5$study_name)
sum(saliva_summary_df5$control_count)
sum(saliva_summary_df5$disease_count)
unique(MetadataDf_saliva_DisCtrl_16s_only$disease)

# studies : 28
# total samples: 3839
# control count: 1780
# disease count: 2059
# diseases:
#  [1] "Control"                     "Nasopharyngeal_Carcinoma"   
#  [3] "Pancreatic_Cancer"           "Myasthenia_Gravis"          
#  [5] "Dental_Caries"               "Oral_SCC_Variants"          
#  [7] "Periodontitis"               "Polycystic_Ovary_Syndrome"  
#  [9] "Cheilitis_Granulomatosa"     "T1D"                        
# [11] "Celiac_Disease"              "Bone_Associated_Disorders"  
# [13] "Oral_Lichen_Planus"          "Chronic_Fatigue_Syndrome"   
# [15] "Other_SCC_Variants"          "Apthous_Infl_Variants"      
# [17] "GERD"                        "COVID"                      
# [19] "Hypothyroidism"              "ECC"                        
# [21] "IBD_GutInflammation"         "Gingivitis"                 
# [23] "Neurodegenerative_Disorders" "CAD_heart_related_syndromes"

## WGS only:
unique(saliva_summary_df8$study_name)
sum(saliva_summary_df8$control_count)
sum(saliva_summary_df8$disease_count)
unique(MetadataDf_saliva_DisCtrl_WGS$disease)

# studies : 5
# total samples: 845
# control count: 500
# disease count: 345
# diseases:
# [1] "Control"           "Pancreatic_Cancer" "COVID"             "Periodontitis"    
# [5] "Dental_Caries"     "Colorectal_Cancer"





#####################################
####################################
######################################

#### Compare the scores across cohorts
### corr between WGS and 16s
library(ggplot2)
library(ggrepel)
library(dplyr)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_Correlation_16s_WGS_AllSpecies.pdf",width = 5, height = 6)
ggplot(
  saliva_CohortWise_HealthAssociationScore,
  aes(x = HealthAssociationScore_WGS,y = HealthAssociationScore_16s, color = HealthAssociationScore_WGS >= 0.80 & HealthAssociationScore_16s >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "HA score WGS",
    y = "HA score 16S"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

dev.off()

cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_WGS,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s, use = "complete.obs")
#cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_WGS,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s, use = "complete.obs", method = "spearman")
#         Pearson's product-moment correlation

# data:  saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_WGS and saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s
# t = 0.55247, df = 224, p-value = 0.5812
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.09406511  0.16658659
# sample estimates:
#        cor 
# 0.03688812



#### 16s vs Overall
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_Correlation_16s_Overall_AllSpecies.pdf",
    width = 5, height = 6)
ggplot(
  saliva_CohortWise_HealthAssociationScore,
  aes(x = HealthAssociationScore_16s,y = HealthAssociationScore, color = HealthAssociationScore_16s >= 0.80 & HealthAssociationScore >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "HA score 16s",
    y = "HA score Overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )
dev.off()

cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s, use = "complete.obs", method = "pearson")
#         Pearson's product-moment correlation

# data:  saliva_CohortWise_HealthAssociationScore$HealthAssociationScore and saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s
# t = 39.743, df = 483, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8524966 0.8944535
# sample estimates:
#       cor 
# 0.8751096

format((cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s, use = "complete.obs", method = "pearson"))$p.value, scientific = TRUE, digits = 3)
# [1] "2.32e-154"


#### WGS vs Overall
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_Correlation_WGS_Overall_AllSpecies.pdf",
    width = 5, height = 6)
ggplot(
  saliva_CohortWise_HealthAssociationScore,
  aes(x = HealthAssociationScore_WGS,y = HealthAssociationScore,color = HealthAssociationScore_WGS >= 0.80 & HealthAssociationScore >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none"
  ) +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "HA score WGS",
    y = "HA score Overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )
dev.off()

cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_WGS, use = "complete.obs", method = "pearson")
#         Pearson's product-moment correlation

# data:  saliva_CohortWise_HealthAssociationScore$HealthAssociationScore and saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_WGS
# t = 6.435, df = 238, p-value = 6.728e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2715814 0.4878211
# sample estimates:
#      cor 
# 0.384972 

format((cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_WGS, use = "complete.obs", method = "pearson"))$p.value, scientific = TRUE, digits = 3)
# "6.73e-10"


#### Overall vs Exposure
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_Correlation_Exposure_Overall_AllSpecies.pdf",
    width = 5, height = 6)
ggplot(
  saliva_CohortWise_HealthAssociationScore,
  aes(x = HealthAssociationScore_16s_exposure,y = HealthAssociationScore, color = HealthAssociationScore_16s_exposure >= 0.80 & HealthAssociationScore >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none"
  ) +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "HA score Exposure",
    y = "HA score Overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

dev.off()

cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s_exposure, use = "complete.obs", method = "pearson")
#         Pearson's product-moment correlation

# data:  saliva_CohortWise_HealthAssociationScore$HealthAssociationScore and saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s_exposure
# t = 7.708, df = 467, p-value = 7.752e-14
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2531077 0.4139077
# sample estimates:
#       cor 
# 0.3359534


format((cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s_exposure, use = "complete.obs", method = "pearson"))$p.value, scientific = TRUE, digits = 3)
# [1] "7.75e-14"




#### Exposure vs Non-Exposure
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_Correlation_Exposure16s_NonExposure16s_AllSpecies.pdf",
    width = 5, height = 6)
ggplot(
  saliva_CohortWise_HealthAssociationScore,
  aes(x = HealthAssociationScore_16s_exposure,y = HealthAssociationScore_16s, color = HealthAssociationScore_16s_exposure >= 0.80 & HealthAssociationScore_16s >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none"
  ) +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "HA score Exposure",
    y = "HA score Non-Exposure"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

dev.off()

cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s_exposure, use = "complete.obs", method = "pearson")
#         Pearson's product-moment correlation

# data:  saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s and saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s_exposure
# t = 3.0122, df = 467, p-value = 0.002734
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.04811141 0.22577757
# sample estimates:
#       cor 
# 0.1380549


format((cor.test(saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s,saliva_CohortWise_HealthAssociationScore$HealthAssociationScore_16s_exposure, use = "complete.obs", method = "pearson"))$p.value, scientific = TRUE, digits = 3)
# [1] "2.73e-03"


