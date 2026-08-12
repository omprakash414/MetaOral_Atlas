############# CoreAssociation for validation cohorts.


library(dplyr)
######## Import the data (2025 and collected in 2026)
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_ValidationCohort.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_ValidationCohort2026.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/JiY_2020_SpDf_Metadata.RData")


## Combine the two validtaion species profiles and metadata and also JiY_2020_SpDf (as it has saliva samples which are not included in the previous two datasets or in discovery cohort but is present in Main Oral Repository)
SpDf_saliva_Validation <- bind_rows(SpDf_saliva_Validation, combined_SpDf, JiY_2020_SpDf)
SpDf_saliva_Validation[is.na(SpDf_saliva_Validation)] <- 0

SpDf_saliva_Validation <- SpDf_saliva_Validation[,colSums(SpDf_saliva_Validation)>0]


## Make some changes in the columns of combined_metadata and then add with previous metadata (2025)
combined_metadata$body_site_category <- "saliva_sputum_oral_wash"
combined_metadata$age <- as.character(combined_metadata$age)
combined_metadata$timepoint <- as.character(combined_metadata$timepoint)
combined_metadata$BMI <- as.character(combined_metadata$BMI)

colnames(combined_metadata)[colnames(combined_metadata) == "severity"] <- "disease_severity"

MetadataDf_saliva_validation <- bind_rows(MetadataDf_saliva_validation, combined_metadata, JiY_2020_metadata)
MetadataDf_saliva_validation <- MetadataDf_saliva_validation[!MetadataDf_saliva_validation$study_name %in% c("CabralD_2017", "PRJNA1013236_Korea","PRJNA1228548_CRC_China2", "PRJNA1293413_Periodontitis_China1"),]

###### modify ChenC_2018_metadata as it has pre-post disease samples along with controls 
library(readxl)
ChenC_2018_metadata <- data.frame(read_excel("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/ChenC_2018_metadata.xlsx"))
rownames(ChenC_2018_metadata) <- ChenC_2018_metadata$sample_id

MetadataDf_saliva_validation[rownames(ChenC_2018_metadata), "subject_id"] <- ChenC_2018_metadata$subject_id
MetadataDf_saliva_validation[rownames(ChenC_2018_metadata), "original_sample_id"] <- ChenC_2018_metadata$original_sample_id
MetadataDf_saliva_validation[rownames(ChenC_2018_metadata), "timepoint"] <- ChenC_2018_metadata$timepoint
MetadataDf_saliva_validation[rownames(ChenC_2018_metadata), "study_condition"] <- ChenC_2018_metadata$study_condition
MetadataDf_saliva_validation$disease <- ifelse(MetadataDf_saliva_validation$study_name == "ChenC_2018" & MetadataDf_saliva_validation$study_condition == "Post-Disease", "Diseased", MetadataDf_saliva_validation$disease)

## Remove T2D samples from SchmidtT_2019B as their count is 5
MetadataDf_saliva_validation <- MetadataDf_saliva_validation[!(MetadataDf_saliva_validation$study_name == "SchmidtT_2019B" & MetadataDf_saliva_validation$disease == "T2D"),]




SpDf_saliva_Validation <- SpDf_saliva_Validation[rownames(MetadataDf_saliva_validation),]
SpDf_saliva_Validation <- SpDf_saliva_Validation[,colSums(SpDf_saliva_Validation)>0]

save(SpDf_saliva_Validation,MetadataDf_saliva_validation, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1SpeciesProfile_Metadata_ValidationCohort.RData")
######## Load the functions
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")

######## Import the combined score (Core/Health/Stability) of Discovery cohort:
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedScores.RData")

selected_saliva_species <- rownames(Combined_Saliva_Scores)


######## Keep only selected species in specie sprofile and then normalize it.
SpDf_saliva_Validation <- SpDf_saliva_Validation[,colnames(SpDf_saliva_Validation) %in% selected_saliva_species]
SpDf_saliva_Validation <- SpDf_saliva_Validation[rowSums(SpDf_saliva_Validation)>0,] # No any row got removed after running this command
SpDf_saliva_Validation <- SpDf_saliva_Validation/rowSums(SpDf_saliva_Validation)

MetadataDf_saliva_validation <- MetadataDf_saliva_validation[rownames(SpDf_saliva_Validation),]
dim(MetadataDf_saliva_validation)
# 4149
######## separate the control samples metadata and species profile
control_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_condition == "Control",]
control_SpDf <- SpDf_saliva_Validation[rownames(control_metadata),]

## One study have multiple disease and they have their own control samoples
control_metadata$study_name <- ifelse(
  !is.na(control_metadata$original_disease) & control_metadata$original_disease == "Control_BC",
  paste0(control_metadata$study_name, "_BC"),
  control_metadata$study_name
)

control_metadata$study_name <- ifelse(
  !is.na(control_metadata$original_disease) & control_metadata$original_disease == "Control_CC",
  paste0(control_metadata$study_name, "_CC"),
  control_metadata$study_name
)

control_metadata$study_name <- ifelse(
  !is.na(control_metadata$original_disease) & control_metadata$original_disease == "Control_PC",
  paste0(control_metadata$study_name, "_PC"),
  control_metadata$study_name
)

selected_species <- colnames(control_SpDf)
control_SpDf$study_name <- control_metadata$study_name
study_list <- unique(control_SpDf$study_name)

######## Run the core influencer function
validationCore_Influencers <- keystoneInfluence(control_SpDf, selected_species, 0.85)

## Summarising r2 df and pr df (pval df)
summary_r2PrDf <- Summarise_r2pval_CoreKeystone(validationCore_Influencers$ENV_fit_summary, rank_scale)

r2_rankedDf <- summary_r2PrDf$r2df_ranked
r2_unrankedDf <- summary_r2PrDf$r2df_unranked
pval_Df <- summary_r2PrDf$pvalue_df


########### Calculate the prevalence of the species in each of the study (same as done in S0 but this time selected species)
prevalDf <- compute_detection(control_SpDf,selected_species,"study_name",study_list)


########### Now calculate the core influencers or core associated species in each of the study cohort (prevalance threshold as 0.85 and r2 threshold as 0.75)
ValidationCoreKeyStoneDf <- data.frame(apply(r2_rankedDf,2,function(x)(ifelse(x>=0.75,1,0))) * apply(prevalDf,2,function(x)(ifelse(x>=0.85,1,0))))



Combined_Saliva_Scores2 <- Combined_Saliva_Scores[rownames(ValidationCoreKeyStoneDf),]

ValidationCoreKeyStoneDf <- cbind(ValidationCoreKeyStoneDf,Combined_Saliva_Scores2)

########### See how many species are being overlapped in each of the study outoff 499. And then plot the plot for those species only.

detected_species_list <- list()

for(i in study_list){
    temp_spDf <- control_SpDf[control_SpDf$study_name == i,]
    temp_spDf$study_name <- NULL
    temp_spDf <- temp_spDf[,colSums(temp_spDf)>0]
    temp_species <- colnames(temp_spDf)
    detected_species_list[[i]] <- temp_species
}



########## Prepare the data for bean plot
bean_data  <- list()
pval_list  <- list()

for (st in study_list) {

  core_vals <- ValidationCoreKeyStoneDf$CoreAssociationScore[
    rownames(ValidationCoreKeyStoneDf) %in% detected_species_list[[st]] &
    ValidationCoreKeyStoneDf[[st]] == 1 &
    ValidationCoreKeyStoneDf$CoreAssociationScore != 0
  ]

  noncore_vals <- ValidationCoreKeyStoneDf$CoreAssociationScore[
    rownames(ValidationCoreKeyStoneDf) %in% detected_species_list[[st]] &
    ValidationCoreKeyStoneDf[[st]] != 1 & 
    ValidationCoreKeyStoneDf$CoreAssociationScore != 0
  ]

  ## ---- store beanplot data ----
  bean_data[[length(bean_data) + 1]] <- core_vals
  bean_data[[length(bean_data) + 1]] <- noncore_vals

  ## ---- Wilcoxon test (core vs non-core) ----
  if (length(core_vals) > 1 && length(noncore_vals) > 1) {
    wt <- wilcox.test(core_vals, noncore_vals)
    pval_list[[st]] <- wt$p.value
  } else {
    pval_list[[st]] <- NA
  }
}

######### Now plot the bean plot
library(beanplot)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_AllStudies.pdf",
    width = 15, height = 7)

par(mar = c(7, 3, 4, 1))

beanplot(
  bean_data,
  names = study_list,
  side = "both",
  what = c(1, 1, 1, 0),
  overallline = "median",
  col = list("#6cd7ec", "#f5f5f5"),
  las = 2,
  cex.axis = 0.9
)

dev.off()



####### Also plot the correlation plot where, detection percent of every species in number of studies is showna nd its Core-Associations and HACK Score is correlated. 
ValidationCoreKeyStoneDf$Detection <- rowSums(ValidationCoreKeyStoneDf[,1:21])/21


## correlation plot of Core-Associations Score and detection percentage.
library(ggplot2)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreScore_vs_detection.pdf", height = 6, width = 5)
ggplot(ValidationCoreKeyStoneDf,aes(x=Detection,y=CoreAssociationScore))+geom_point(size=2)+geom_smooth(method='lm')+ylim(0,1)+xlim(0,1)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+xlab("")+ylab("") 
dev.off()

Corr_CoreAsso <- cor.test(ValidationCoreKeyStoneDf$Detection, ValidationCoreKeyStoneDf$CoreAssociationScore, method = 'spearman')
#         Spearman's rank correlation rho

# data:  ValidationCoreKeyStoneDf$Detection and ValidationCoreKeyStoneDf$CoreAssociationScore
# S = 3681076, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8222432
# > Corr_CoreAsso$p.value
# [1]  9.562937e-124

## correlation plot of HACK Score and detection percentage.
library(ggplot2)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1HACKScore_vs_detection.pdf", height = 6, width = 5)
ggplot(ValidationCoreKeyStoneDf,aes(x=Detection,y=HAC_Score_RankScaled))+geom_point(size=2)+geom_smooth(method='lm')+ylim(0,1)+xlim(0,1)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+xlab("")+ylab("") 
dev.off()

Corr_HACK<- cor.test(ValidationCoreKeyStoneDf$Detection, ValidationCoreKeyStoneDf$HAC_Score_RankScaled, method = 'spearman')
#         Spearman's rank correlation rho

# data:  ValidationCoreKeyStoneDf$Detection and ValidationCoreKeyStoneDf$HAC_Score_RankScaled
# S = 6951313, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6643256
# > Corr_HACK$p.value
# [1] 7.915137e-65

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_Validation_Workspace.RData")



######################### Calculate the healthassociation difference using core/non-core species, similary also do for stability and overall health association
###### Health Association
bean_data_HA  <- list()
pval_list_HA <- list()

for (st in study_list) {

  core_vals <- ValidationCoreKeyStoneDf$HealthAssociationScore[
    rownames(ValidationCoreKeyStoneDf) %in% detected_species_list[[st]] &
    ValidationCoreKeyStoneDf[[st]] == 1 &
    ValidationCoreKeyStoneDf$HealthAssociationScore != 0
  ]

  noncore_vals <- ValidationCoreKeyStoneDf$HealthAssociationScore[
    rownames(ValidationCoreKeyStoneDf) %in% detected_species_list[[st]] &
    ValidationCoreKeyStoneDf[[st]] != 1 &
    ValidationCoreKeyStoneDf$HealthAssociationScore != 0
  ]

  ## ---- store beanplot data ----
  bean_data_HA[[length(bean_data_HA) + 1]] <- core_vals
  bean_data_HA[[length(bean_data_HA) + 1]] <- noncore_vals

  ## ---- Wilcoxon test (core vs non-core) ----
  if (length(core_vals) > 1 && length(noncore_vals) > 1) {
    wt <- wilcox.test(core_vals, noncore_vals)
    pval_list_HA[[st]] <- wt$p.value
  } else {
    pval_list_HA[[st]] <- NA
  }
}


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1HealthAssociation_AllStudies.pdf",
    width = 15, height = 7)

par(mar = c(7, 3, 4, 1))

beanplot(
  bean_data_HA,
  names = study_list,
  side = "both",
  what = c(1, 1, 1, 0),
  overallline = "median",
  col = list("#6cd7ec", "#f5f5f5"),
  las = 2,
  cex.axis = 0.9
)

dev.off()




###### Stability Association
bean_data_SA  <- list()
pval_list_SA <- list()

for (st in study_list) {

  core_vals <- ValidationCoreKeyStoneDf$StabilityAssociationScore[
    rownames(ValidationCoreKeyStoneDf) %in% detected_species_list[[st]] &
    ValidationCoreKeyStoneDf[[st]] == 1 &
    ValidationCoreKeyStoneDf$StabilityAssociationScore != 0
  ]

  noncore_vals <- ValidationCoreKeyStoneDf$StabilityAssociationScore[
    rownames(ValidationCoreKeyStoneDf) %in% detected_species_list[[st]] &
    ValidationCoreKeyStoneDf[[st]] != 1 &
    ValidationCoreKeyStoneDf$StabilityAssociationScore != 0
  ]

  ## ---- store beanplot data ----
  bean_data_SA[[length(bean_data_SA) + 1]] <- core_vals
  bean_data_SA[[length(bean_data_SA) + 1]] <- noncore_vals

  ## ---- Wilcoxon test (core vs non-core) ----
  if (length(core_vals) > 1 && length(noncore_vals) > 1) {
    wt <- wilcox.test(core_vals, noncore_vals)
    pval_list_SA[[st]] <- wt$p.value
  } else {
    pval_list_SA[[st]] <- NA
  }
}


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1StabilityAssociation_AllStudies.pdf",
    width = 15, height = 7)

par(mar = c(7, 3, 4, 1))

beanplot(
  bean_data_SA,
  names = study_list,
  side = "both",
  what = c(1, 1, 1, 0),
  overallline = "median",
  col = list("#6cd7ec", "#f5f5f5"),
  las = 2,
  cex.axis = 0.9
)

dev.off()





###### Overall Association (Ranked HAC Score)
bean_data_HAC  <- list()
pval_list_HAC <- list()

for (st in study_list) {

  core_vals <- ValidationCoreKeyStoneDf$HAC_Score_RankScaled[
    rownames(ValidationCoreKeyStoneDf) %in% detected_species_list[[st]] &
    ValidationCoreKeyStoneDf[[st]] == 1 &
    ValidationCoreKeyStoneDf$HAC_Score_RankScaled != 0
  ]

  noncore_vals <- ValidationCoreKeyStoneDf$HAC_Score_RankScaled[
    rownames(ValidationCoreKeyStoneDf) %in% detected_species_list[[st]] &
    ValidationCoreKeyStoneDf[[st]] != 1 &
    ValidationCoreKeyStoneDf$HAC_Score_RankScaled != 0
  ]

  ## ---- store beanplot data ----
  bean_data_HAC[[length(bean_data_HAC) + 1]] <- core_vals
  bean_data_HAC[[length(bean_data_HAC) + 1]] <- noncore_vals

  ## ---- Wilcoxon test (core vs non-core) ----
  if (length(core_vals) > 1 && length(noncore_vals) > 1) {
    wt <- wilcox.test(core_vals, noncore_vals)
    pval_list_HAC[[st]] <- wt$p.value
  } else {
    pval_list_HAC[[st]] <- NA
  }
}


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1HAC_Score_AllStudies.pdf",
    width = 15, height = 7)

par(mar = c(7, 3, 4, 1))

beanplot(
  bean_data_HAC,
  names = study_list,
  side = "both",
  what = c(1, 1, 1, 0),
  overallline = "median",
  col = list("#6cd7ec", "#f5f5f5"),
  las = 2,
  cex.axis = 0.9
)

dev.off()




ValidationCoreKeyStoneDf$count <- rowSums(ValidationCoreKeyStoneDf[,1:21])

### Plot how many studies in which the species is coming as core associated.
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_HACKScore_BoxPlot.pdf", width = 8, height = 6)
boxplot(ValidationCoreKeyStoneDf[,26] ~ cut(ValidationCoreKeyStoneDf$count,
                                  breaks = c(0, 1, 16, 12)), col = c("#e09549", "#74abe6", "#93e688"), outline = TRUE, ylab = "HACK Score (Discovery Cohort)", xlab = "No. of studies with core association taxa")
dev.off()


groups <- cut(ValidationCoreKeyStoneDf$count,
              breaks = c(0, 1, 16, 12))


library(dunn.test)
HACK_dunntest <- dunn.test(x = ValidationCoreKeyStoneDf[,26],g = groups,method = "bh")
HACK_dunntest
# $chi2
# [1] 24.53942

# $Z
# [1] -4.7723710 -3.5470444 -0.7181528

# $P
# [1] 9.103485e-07 1.947894e-04 2.363315e-01

# $P.adjusted
# [1] 2.731046e-06 2.921841e-04 2.363315e-01

# $comparisons
# [1] "(0,1] - (1,12]"   "(0,1] - (12,16]"  "(1,12] - (12,16]"


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_CoreScore_BoxPlot.pdf", width = 8, height = 6)
boxplot(ValidationCoreKeyStoneDf[,22] ~ cut(ValidationCoreKeyStoneDf$count,
                                  breaks = c(0, 1, 16, 12)), col = c("#e09549", "#74abe6", "#93e688"), outline = TRUE, ylab = "Core Score (Discovery Cohort)", xlab = "No. of studies with core association taxa")
dev.off()

Core_dunntest <- dunn.test(x = ValidationCoreKeyStoneDf[,22],g = groups,method = "bh")
Core_dunntest
# $chi2
# [1] 74.57256

# $Z
# [1] -6.490349 -8.186831 -4.835960

# $P
# [1] 4.281900e-11 1.340968e-16 6.625223e-07

# $P.adjusted
# [1] 6.422850e-11 4.022905e-16 6.625223e-07

# $comparisons
# [1] "(0,1] - (1,12]"   "(0,1] - (12,16]"  "(1,12] - (12,16]"






save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_Validation_Workspace.RData")




###### extract the 1/0 matrix of validation and sent it to sir.
ValidationCoreSp_StudyWise <- ValidationCoreKeyStoneDf[,1:21]
Validation_CoreSp_StudyWise <- data.frame(t(ValidationCoreSp_StudyWise))

save(Validation_CoreSp_StudyWise, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_Validation_BinaryMatrix.RData")
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_Validation_Workspace.RData")
