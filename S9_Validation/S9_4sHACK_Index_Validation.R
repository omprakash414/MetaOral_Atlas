
library(dplyr)
library(vegan) # for distance calculation
library(compositions) # for CLR transformation
library(ggplot2) # for plotting


S9_2HealthAssociation_Validation_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_Validation_Workspace.RData", envir = S9_2HealthAssociation_Validation_Workspace)
attach(S9_2HealthAssociation_Validation_Workspace)
ChenC_2018_SpDf <- ChenC_2018_SpDf
ChenC_2018_ctrl <- ChenC_2018_ctrl
ChenC_2018_dis <- ChenC_2018_dis
ChenC_2018_postDis <- ChenC_2018_postDis
FanX_2018a_SpDf <- FanX_2018a_SpDf
FanX_2018a_ctrl_Al <- FanX_2018a_ctrl_Al
FanX_2018a_dis_Al <- FanX_2018a_dis_Al
FanX_2018a_ctrl_NonAl <- FanX_2018a_ctrl_NonAl
FanX_2018a_dis_NonAl <- FanX_2018a_dis_NonAl
ZhangT_2020_SpDf <- ZhangT_2020_SpDf
ZhangT_2020_ctrl <- ZhangT_2020_ctrl
ZhangT_2020_dis <- ZhangT_2020_dis
ZhangT_2020_dis_rem <- ZhangT_2020_dis_rem
SchmidtT_2019B_SpDf <- SchmidtT_2019B_SpDf
SchmidtT_2019B_ctrl <- SchmidtT_2019B_ctrl
SchmidtT_2019B_dis <- SchmidtT_2019B_dis
ChenX_2025_SpDf <- ChenX_2025_SpDf
ChenX_2025_ctrl <- ChenX_2025_ctrl
ChenX_2025_dis <- ChenX_2025_dis
TakayanagiK_2023_SpDf <- TakayanagiK_2023_SpDf
TakayanagiK_2023_ctrl <- TakayanagiK_2023_ctrl
TakayanagiK_2023_dis <- TakayanagiK_2023_dis
WuZ_2025_SpDf <- WuZ_2025_SpDf
WuZ_2025_ctrl <- WuZ_2025_ctrl
WuZ_2025_dis <- WuZ_2025_dis
WuZ_2025_dis_nonmet <- WuZ_2025_dis_nonmet
ZhangL_2023_SpDf <- ZhangL_2023_SpDf
ZhangL_2023_ctrl <- ZhangL_2023_ctrl
ZhangL_2023_dis <- ZhangL_2023_dis
ChenJ_2021_SpDf <- ChenJ_2021_SpDf
ChenJ_2021_ctrl <- ChenJ_2021_ctrl
ChenJ_2021_dis_OSCC <- ChenJ_2021_dis_OSCC
ChenJ_2021_dis_OVH <- ChenJ_2021_dis_OVH

CirsteaM_2022_SpDf <- CirsteaM_2022_SpDf
CirsteaM_2022_ctrl <- CirsteaM_2022_ctrl
CirsteaM_2022_dis <- CirsteaM_2022_dis

FinkelsteinS_2025_SpDf <- FinkelsteinS_2025_SpDf
FinkelsteinS_2025_ctrl <- FinkelsteinS_2025_ctrl
FinkelsteinS_2025_dis <- FinkelsteinS_2025_dis

IglesiasA_2024_SpDf <- IglesiasA_2024_SpDf
IglesiasA_2024_ctrl <- IglesiasA_2024_ctrl
IglesiasA_2024_dis <- IglesiasA_2024_dis

LiuY_2021_SpDf <- LiuY_2021_SpDf
LiuY_2021_ctrl <- LiuY_2021_ctrl
LiuY_2021_dis <- LiuY_2021_dis

NearingJ_2023_SpDf <- NearingJ_2023_SpDf
NearingJ_2023_ctrl_BC <- NearingJ_2023_ctrl_BC
NearingJ_2023_ctrl_CC <- NearingJ_2023_ctrl_CC
NearingJ_2023_ctrl_PC <- NearingJ_2023_ctrl_PC
NearingJ_2023_dis_BC <- NearingJ_2023_dis_BC
NearingJ_2023_dis_CC <- NearingJ_2023_dis_CC
NearingJ_2023_dis_PC <- NearingJ_2023_dis_PC

RelvasM_2021_SpDf <- RelvasM_2021_SpDf
RelvasM_2021_ctrl <- RelvasM_2021_ctrl
RelvasM_2021_dis_Dg1 <- RelvasM_2021_dis_Dg1
RelvasM_2021_dis_Dg2 <- RelvasM_2021_dis_Dg2
RelvasM_2021_dis_Dg3 <- RelvasM_2021_dis_Dg3
RelvasM_2021_dis_Pg1 <- RelvasM_2021_dis_Pg1
RelvasM_2021_dis_Pg2 <- RelvasM_2021_dis_Pg2
RelvasM_2021_dis_Pg3 <- RelvasM_2021_dis_Pg3

JiY_2020_SpDf <- JiY_2020_SpDf
JiY_2020_ctrl <- JiY_2020_ctrl
JiY_2020_dis <- JiY_2020_dis
JiY_2020_post_dis <- JiY_2020_post_dis

detected_species_list <- detected_species_list
detach(S9_2HealthAssociation_Validation_Workspace)
rm(S9_2HealthAssociation_Validation_Workspace)



S9_1CoreAssociation_Validation_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_Validation_Workspace.RData", envir = S9_1CoreAssociation_Validation_Workspace)
attach(S9_1CoreAssociation_Validation_Workspace)
MetadataDf_saliva_validation <- MetadataDf_saliva_validation
SpDf_saliva_Validation <- SpDf_saliva_Validation
study_list <- study_list
Combined_Saliva_Scores <- Combined_Saliva_Scores
detach(S9_1CoreAssociation_Validation_Workspace)
rm(S9_1CoreAssociation_Validation_Workspace)

colnames(Combined_Saliva_Scores) <- c("CoreAssociationScore","HealthAssociationScore","StabilityAssociationScore","sHACK_UnRanked","sHACK_Ranked")


#sHACK_top_species <- c("Fusobacterium_periodonticum","Solobacterium_moorei","Eubacterium_sulci","Actinomyces_graevenitzii","Gemella_sanguinis","Lachnoanaerobaculum_umeaense","Campylobacter_concisus","Prevotella_oulorum","Veillonella_parvula","Veillonella_rogosae","Cardiobacterium_hominis","Capnocytophaga_gingivalis","Neisseria_elongata","Prevotella_shahii","Haemophilus_parainfluenzae","Prevotella_melaninogenica","Leptotrichia_goodfellowii","Campylobacter_showae","Capnocytophaga_leadbetteri","Prevotella_salivae","Streptococcus_infantis","Prevotella_maculosa","Leptotrichia_hongkongensis","Leptotrichia_hofstadii","Haemophilus_pittmaniae","Actinomyces_odontolyticus","Propionibacterium_propionicum","Cardiobacterium_valvarum")
Top_28_species <- c("Fusobacterium_periodonticum","Eubacterium_sulci","Lachnoanaerobaculum_umeaense","Cardiobacterium_hominis","Neisseria_elongata","Campylobacter_showae","Leptotrichia_goodfellowii","Catonella_morbi","Haemophilus_parainfluenzae","Capnocytophaga_leadbetteri","Capnocytophaga_gingivalis","Gemella_sanguinis","Haemophilus_pittmaniae","Cardiobacterium_valvarum","Propionibacterium_propionicum","Prevotella_aurantiaca","Neisseria_meningitidis","Actinomyces_massiliensis","Veillonella_rogosae","Neisseria_flavescens","Corynebacterium_durum","Streptobacillus_moniliformis","Haemophilus_sputorum","Eubacterium_yurii","Porphyromonas_catoniae","Prevotella_shahii","Neisseria_oralis","Alloprevotella_rava")
sHACK_top_species <- rownames(Combined_Saliva_Scores[rowSums(Combined_Saliva_Scores[,1:3] >= 0.80) == 3, ])

common_species <- intersect(Top_28_species, sHACK_top_species)

###### Write a function to calculate the sHACK23 Score and to calculate the dysbiosis score for microbiome samples.

rank_scale=function(x)
{
  # x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

## sHack top species score calculation
sHACK_Score <- function(data)
{
  score <- rowSums(apply(data[,intersect(colnames(data),sHACK_top_species)],2,rank_scale))
  cat("sHACK score computed\n")
  return(score)
}

# Now write a function to calculate the dysbiosis score using three different distance matrices for microbiome samples. 
## Function for distance calculation
compute_distance_matrix <- function(data, method) {
  if (method == "bray") {
    return(as.matrix(vegdist(data, method = "bray")))
    
  } else if (method == "canberra") {
    return(as.matrix(vegdist(data, method = "canberra")))
    
  } else if (method == "aitchison") {
    # CLR transform requires all-positive values
    data[data == 0] <- 1e-6
    
    clr_mat <- clr(data)
    return(as.matrix(dist(clr_mat, method = "euclidean")))
    
  } else {
    stop("Unknown distance method. Use: bray, canberra, aitchison")
  }
}


# Dysbiosis score calculation
dysbiosis_score <- function(data, control_group, method = "bray") {
  
  all_samples <- rownames(data)
  disease_samples <- setdiff(all_samples, control_group)
  
  cat("Distance:", method, "\n")
  cat("Total samples:", length(all_samples), "\n")
  cat("Control samples:", length(control_group), "\n")
  cat("Disease samples:", length(disease_samples), "\n")
  
  ## ---- Disease samples ----
  cat("Processing disease samples...\n")
  
  sub_data <- data[c(disease_samples, control_group), , drop = FALSE]
  
  dmat <- compute_distance_matrix(sub_data, method)
  diag(dmat) <- NA
  
  dmat_subset <- dmat[disease_samples, control_group, drop = FALSE]
  disease_scores <- apply(dmat_subset, 1, function(x) median(x, na.rm = TRUE))
  
  ## ---- Control samples ----
  cat("Processing control samples...\n")
  
  control_data <- data[control_group, , drop = FALSE]
  control_dmat <- compute_distance_matrix(control_data, method)
  diag(control_dmat) <- NA
  
  control_scores <- apply(control_dmat, 1, function(x) median(x, na.rm = TRUE))
  
  ## ---- Combine and reorder ----
  full_scores <- c(control_scores, disease_scores)
  full_scores <- full_scores[rownames(data)]
  
  cat("Dysbiosis score completed:", method, "\n")
  return(full_scores)
}

dysbiosis_score_controls <- function(data, method = "bray") {
  
  control_samples <- rownames(data)
  
  cat("Distance:", method, "\n")
  cat("Total control samples:", length(control_samples), "\n")
  cat("Processing control samples only...\n")
  
  ## Compute distance matrix
  dmat <- compute_distance_matrix(data, method)
  diag(dmat) <- NA
  
  ## Dysbiosis = median distance to all other controls
  control_scores <- apply(dmat, 1, function(x) median(x, na.rm = TRUE))
  
  ## Keep original order
  control_scores <- control_scores[rownames(data)]
  
  cat("Control dysbiosis score completed:", method, "\n")
  return(control_scores)
}

############## Run each study's data using the above functions.
combined_Sp_Df_validation <- SpDf_saliva_Validation

##### ChenC_2018
ChenC_2018_SpDf$study_name <- NULL
ChenC_2018_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(ChenC_2018_SpDf))
ChenC_2018_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score(data = ChenC_2018_SpDf,control_group = ChenC_2018_ctrl,method = "bray"))
ChenC_2018_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score(data = ChenC_2018_SpDf,control_group = ChenC_2018_ctrl,method = "aitchison"))
ChenC_2018_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score(data = ChenC_2018_SpDf,control_group = ChenC_2018_ctrl,method = "canberra"))


##### FanX_2018a_ctrl_Al
FanX_2018a_SpDf$study_name <- NULL
FanX_2018a_Al_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(FanX_2018a_SpDf)) 
FanX_2018a_Al_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score(data = FanX_2018a_SpDf,control_group = FanX_2018a_ctrl_Al,method = "bray"))
FanX_2018a_Al_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score(data = FanX_2018a_SpDf,control_group = FanX_2018a_ctrl_Al,method = "aitchison"))
FanX_2018a_Al_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score(data = FanX_2018a_SpDf,control_group = FanX_2018a_ctrl_Al,method = "canberra"))  

##### ZhangT_2020
ZhangT_2020_SpDf$study_name <- NULL
ZhangT_2020_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(ZhangT_2020_SpDf)) 
ZhangT_2020_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score(data = ZhangT_2020_SpDf,control_group = ZhangT_2020_ctrl,method = "bray"))
ZhangT_2020_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score(data = ZhangT_2020_SpDf,control_group = ZhangT_2020_ctrl,method = "aitchison"))  
ZhangT_2020_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score(data = ZhangT_2020_SpDf,control_group = ZhangT_2020_ctrl,method = "canberra")) 

##### SchmidtT_2019B
SchmidtT_2019B_SpDf$study_name <- NULL
SchmidtT_2019B_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(SchmidtT_2019B_SpDf)) 
SchmidtT_2019B_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score(data = SchmidtT_2019B_SpDf,control_group = SchmidtT_2019B_ctrl,method = "bray"))
SchmidtT_2019B_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score(data = SchmidtT_2019B_SpDf,control_group = SchmidtT_2019B_ctrl,method = "aitchison"))  
SchmidtT_2019B_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score(data = SchmidtT_2019B_SpDf,control_group = SchmidtT_2019B_ctrl,method = "canberra"))  

###### ChenX_2025 
ChenX_2025_SpDf$study_name <- NULL
ChenX_2025_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(ChenX_2025_SpDf))
ChenX_2025_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score(data = ChenX_2025_SpDf,control_group = ChenX_2025_ctrl,method = "bray"))
ChenX_2025_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score(data = ChenX_2025_SpDf,control_group = ChenX_2025_ctrl,method = "aitchison"))
ChenX_2025_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score(data = ChenX_2025_SpDf,control_group = ChenX_2025_ctrl,method = "canberra"))

###### TakayanagiK_2023
TakayanagiK_2023_SpDf$study_name <- NULL
TakayanagiK_2023_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(TakayanagiK_2023_SpDf))
TakayanagiK_2023_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score(data = TakayanagiK_2023_SpDf,control_group = TakayanagiK_2023_ctrl,method = "bray"))
TakayanagiK_2023_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score(data = TakayanagiK_2023_SpDf,control_group = TakayanagiK_2023_ctrl,method = "aitchison"))
TakayanagiK_2023_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score(data = TakayanagiK_2023_SpDf,control_group = TakayanagiK_2023_ctrl,method = "canberra"))  

###### WuZ_2025
WuZ_2025_SpDf$study_name <- NULL
WuZ_2025_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(WuZ_2025_SpDf))
WuZ_2025_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score(data = WuZ_2025_SpDf,control_group = WuZ_2025_ctrl,method = "bray"))
WuZ_2025_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score(data = WuZ_2025_SpDf,control_group = WuZ_2025_ctrl,method = "aitchison"))
WuZ_2025_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score(data = WuZ_2025_SpDf,control_group = WuZ_2025_ctrl,method = "canberra"))


####### ZhangL_2023
ZhangL_2023_SpDf$study_name <- NULL
ZhangL_2023_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(ZhangL_2023_SpDf))
ZhangL_2023_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score(data = ZhangL_2023_SpDf,control_group = ZhangL_2023_ctrl,method = "bray"))
ZhangL_2023_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score(data = ZhangL_2023_SpDf,control_group = ZhangL_2023_ctrl,method = "aitchison"))
ZhangL_2023_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score(data = ZhangL_2023_SpDf,control_group = ZhangL_2023_ctrl,method = "canberra"))


####### LicandroH_2023
LicandroH_2023_SpDf <- combined_Sp_Df_validation[MetadataDf_saliva_validation$study_name == "LicandroH_2023",]
LicandroH_2023_SpDf$study_name <- NULL
LicandroH_2023_sHACK_Score <- data.frame("sHACK_Score" = sHACK_Score(LicandroH_2023_SpDf))
LicandroH_2023_Dysbiosis_Bray <- data.frame("Dysbiosis_score_bray" = dysbiosis_score_controls(data = LicandroH_2023_SpDf,method = "bray"))
LicandroH_2023_Dysbio_Aitchison <- data.frame("Dysbiosis_score_aitchison" = dysbiosis_score_controls(data = LicandroH_2023_SpDf,method = "aitchison"))
LicandroH_2023_Dysbio_Canberra <- data.frame("Dysbiosis_score_canberra" = dysbiosis_score_controls(data = LicandroH_2023_SpDf,method = "canberra"))



############ Now combine the score per study and then add its metadata as control and disease. 
ChenC_2018_Scores <- cbind(ChenC_2018_sHACK_Score, ChenC_2018_Dysbiosis_Bray, ChenC_2018_Dysbio_Aitchison, ChenC_2018_Dysbio_Canberra)
ChenC_2018_Scores$study_name <- "ChenC_2018"

FanX_2018a_Scores <- cbind(FanX_2018a_sHACK_Score, FanX_2018a_Dysbiosis_Bray, FanX_2018a_Dysbio_Aitchison, FanX_2018a_Dysbio_Canberra)
FanX_2018a_Scores$study_name <- "FanX_2018a"

# ZhangT has control/disease/remission so add three of them in disease status
ZhangT_2020_Scores <- cbind(ZhangT_2020_sHACK_Score, ZhangT_2020_Dysbiosis_Bray, ZhangT_2020_Dysbio_Aitchison, ZhangT_2020_Dysbio_Canberra)
ZhangT_2020_Scores$study_name <- "ZhangT_2020"

SchmidtT_2019B_Scores <- cbind(SchmidtT_2019B_sHACK_Score, SchmidtT_2019B_Dysbiosis_Bray, SchmidtT_2019B_Dysbio_Aitchison, SchmidtT_2019B_Dysbio_Canberra)
SchmidtT_2019B_Scores$study_name <- "SchmidtT_2019B"

ChenX_2025_Scores <- cbind(ChenX_2025_sHACK_Score, ChenX_2025_Dysbiosis_Bray, ChenX_2025_Dysbio_Aitchison, ChenX_2025_Dysbio_Canberra)
ChenX_2025_Scores$study_name <- "ChenX_2025"  

TakayanagiK_2023_Scores <- cbind(TakayanagiK_2023_sHACK_Score, TakayanagiK_2023_Dysbiosis_Bray, TakayanagiK_2023_Dysbio_Aitchison, TakayanagiK_2023_Dysbio_Canberra)
TakayanagiK_2023_Scores$study_name <- "TakayanagiK_2023"  

WuZ_2025_Scores <- cbind(WuZ_2025_sHACK_Score, WuZ_2025_Dysbiosis_Bray, WuZ_2025_Dysbio_Aitchison, WuZ_2025_Dysbio_Canberra)
WuZ_2025_Scores$study_name <- "WuZ_2025"

ZhangL_2023_Scores <- cbind(ZhangL_2023_sHACK_Score, ZhangL_2023_Dysbiosis_Bray, ZhangL_2023_Dysbio_Aitchison, ZhangL_2023_Dysbio_Canberra)
ZhangL_2023_Scores$study_name <- "ZhangL_2023"

LicandroH_2023_Scores <- cbind(LicandroH_2023_sHACK_Score, LicandroH_2023_Dysbiosis_Bray, LicandroH_2023_Dysbio_Aitchison, LicandroH_2023_Dysbio_Canberra)
LicandroH_2023_Scores$study_name <- "LicandroH_2023"

## combine all the scores and also have 
combined_Scores_Validation <- bind_rows(
  ChenC_2018_Scores,
  #FanX_2018a_Scores,
  ZhangT_2020_Scores,
  SchmidtT_2019B_Scores,
  ChenX_2025_Scores,
  TakayanagiK_2023_Scores,
  WuZ_2025_Scores,
  ZhangL_2023_Scores,
  LicandroH_2023_Scores
)

# Add study_condition, disease, original_disease, timepoint, subject_id to the combined_Scores_Validation
combined_Scores_Validation$original_sample_id <- MetadataDf_saliva_validation$original_sample_id[match(rownames(combined_Scores_Validation), rownames(MetadataDf_saliva_validation))]
combined_Scores_Validation$subject_id <- MetadataDf_saliva_validation$subject_id[match(rownames(combined_Scores_Validation), rownames(MetadataDf_saliva_validation))]
combined_Scores_Validation$timepoint <- MetadataDf_saliva_validation$timepoint[match(rownames(combined_Scores_Validation), rownames(MetadataDf_saliva_validation))]
combined_Scores_Validation$study_condition <- MetadataDf_saliva_validation$study_condition[match(rownames(combined_Scores_Validation), rownames(MetadataDf_saliva_validation))]
combined_Scores_Validation$disease <- MetadataDf_saliva_validation$disease[match(rownames(combined_Scores_Validation), rownames(MetadataDf_saliva_validation))]


############ save the workspace.
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_4sHACK_Index_Validation_Workspace.RData")

########### save the combined scores, combined SpDf and combined saliva score in one RData file
save(combined_Scores_Validation, combined_Sp_Df_validation, Combined_Saliva_Scores, MetadataDf_saliva_validation,file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_4sHACK_Index_Validation_DataDf.RData")





###################### Compare the sHACK Score and Dysbiosis Scores based on the disease_condition
study_list <- unique(combined_Scores_Validation$study_name)
study1 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[1],])
study2 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[2],])
study3 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[3],])
study4 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[4],])
study5 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[5],])
study6 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[6],])
study7 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[7],])
study8 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[8],])
study9 <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_name == study_list[9],])

combined_Scores_Validation$study_condition <- ifelse(combined_Scores_Validation$study_condition == "Remission IBD_GutInflammation", "Remission",combined_Scores_Validation$study_condition)

Control <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_condition=="Control",])
Disease <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_condition=="Diseased",])
Remission <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_condition=="Remission",])
PostDisease <- rownames(combined_Scores_Validation[combined_Scores_Validation$study_condition=="Post-Disease",])


combined_Scores_Validation$SelectTaxaRankScale <- rowMeans(apply(combined_Sp_Df_validation[rownames(combined_Scores_Validation),sHACK_top_species],2,rank_scale))
combined_Scores_Validation$SelectTaxaDiversity <- diversity(apply(combined_Sp_Df_validation[rownames(combined_Scores_Validation),sHACK_top_species],2,rank_scale))
combined_Scores_Validation$NonSelectTaxaRankScale <- rowMeans(apply(combined_Sp_Df_validation[rownames(combined_Scores_Validation),setdiff(colnames(combined_Sp_Df_validation),sHACK_top_species)],2,rank_scale))
combined_Scores_Validation$NonSelectTaxaDiversity <- diversity(apply(combined_Sp_Df_validation[rownames(combined_Scores_Validation),setdiff(colnames(combined_Sp_Df_validation),sHACK_top_species)],2,rank_scale))

combined_Scores_Validation$temp_OMWI <- (log((combined_Scores_Validation$SelectTaxaRankScale+0.000001)/(combined_Scores_Validation$NonSelectTaxaRankScale+0.000001),10))

refined_Scores_Validation <- combined_Scores_Validation[!combined_Scores_Validation$study_name %in% "LicandroH_2023",]
refined_Scores_Validation$study_condition <- ifelse(refined_Scores_Validation$study_condition == "Control",1,0)

library(metafor)
compute_meta_effsize <- function(data,var1,var2,grouping_variable,grouping_list)
{
  temp_meta <- data.frame(matrix(0,length(grouping_list),7))
  colnames(temp_meta) <- c("dataset","m1i","m2i","sd1i","sd2i","n1i","n2i")
  for(i in 1:length(grouping_list))
  {
    group <- grouping_list[i]
    temp_meta[i,1] <- group
    print(group)
    data_group_Case <- data[(data[,var2]==1)&(data[,grouping_variable]==group),var1]
    data_group_Control <- data[(data[,var2]==0)&(data[,grouping_variable]==group),var1]
    temp_meta[i,2] <- mean(data_group_Case)
    temp_meta[i,3] <- mean(data_group_Control)
    temp_meta[i,4] <- sd(data_group_Case)
    temp_meta[i,5] <- sd(data_group_Control)
    temp_meta[i,6] <- length(data_group_Case)
    temp_meta[i,7] <- length(data_group_Control)
    print(temp_meta)
  }
  temp_meta <- mutate(temp_meta,study_id=grouping_list)
  rownames(temp_meta) <- grouping_list
  #temp_meta <- temp_meta %>% select(study_id, ri:ni)
  temp_meta <- escalc(measure="SMD",m1i=m1i,m2i=m2i,sd1i=sd1i,sd2i=sd2i,n1i=n1i,n2i=n2i,data=temp_meta)
  res <- rma(yi, vi, data=temp_meta)
  res$ids <- rownames(temp_meta)
  #res$slabs <- rownames(temp_meta)
  res_list <- list("df_studies"=temp_meta,"model"=res)
  return(res_list)
}

res_OMWI <- compute_meta_effsize(refined_Scores_Validation,"temp_OMWI","study_condition","study_name",unique(refined_Scores_Validation$study_name))

res_DysBray <- compute_meta_effsize(refined_Scores_Validation,"Dysbiosis_score_bray","study_condition","study_name",unique(refined_Scores_Validation$study_name))

res_DysAit <- compute_meta_effsize(refined_Scores_Validation,"Dysbiosis_score_aitchison","study_condition","study_name",unique(refined_Scores_Validation$study_name))

res_DysCan <- compute_meta_effsize(refined_Scores_Validation,"Dysbiosis_score_canberra","study_condition","study_name",unique(refined_Scores_Validation$study_name))


res_All <- data.frame(
  Score = c("WISH_Index","DysbiosisBray","DysbiosisAitchison","DysbiosisCanberra"),
  Effect = c(as.numeric(res_OMWI$model$b),as.numeric(res_DysBray$model$b),as.numeric(res_DysAit$model$b),as.numeric(res_DysCan$model$b)),
  Pval = c(res_OMWI$model$pval,res_DysBray$model$pval,res_DysAit$model$pval,res_DysCan$model$pval))

res_All$logP <- -log10(res_All$Pval)


res_All$LabelColor <- ifelse(res_All$Pval < 0.05 & res_All$Effect > 0, "Positive",
  ifelse(res_All$Pval < 0.05 & res_All$Effect < 0, "Negative", "NonSig"))

library(ggplot2)
library(ggrepel)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_4Validation_MetaAnalysis_sHACKScore_DysbiosisScore.pdf", width = 8, height = 6)
ggplot(res_All,
       aes(x = Effect, y = logP, color = LabelColor)) +
  
  geom_point(size = 3) +
  
  geom_text_repel(
    aes(label = Score),
    size = 5,
    max.overlaps = 45
  ) +
  
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), 
             color = "red", 
             linetype = "dashed", 
             linewidth = 0.6) +

  scale_color_manual(values = c(
    "Positive" = "#4575b4",
    "Negative" = "#d73027",
    "NonSig"   = "grey60"
  )) +
  
  theme_minimal() +
  labs(title = "Comparison of Microbiome Indices",
       x = "Effect Size",
       y = "-log10(P-value)") +
  
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.position = "none"
  )
dev.off()



############# Now do the meta-analysis of followup distances wrt different indices (sAHCK and All Dysbiosis Indices)
S9_3StabilityAssociation_Validation_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3StabilityAssociation_Validation_Workspace.RData", envir = S9_3StabilityAssociation_Validation_Workspace)
attach(S9_3StabilityAssociation_Validation_Workspace)
ZhangT_2020_metadata_followup <- ZhangT_2020_metadata_followup
SchmidtT_2019B_metadata_followup <- SchmidtT_2019B_metadata_followup
ChenC_2018_metadata_followup <- ChenC_2018_metadata_followup
LincardonT_2023_metadata_followup <- LincardonT_2023_metadata_followup
detach(S9_3StabilityAssociation_Validation_Workspace)
rm(S9_3StabilityAssociation_Validation_Workspace)


library(vegan)
library(compositions)

bray_followup <- function(data,follow_up)
{
  follow_up_dist <- as.data.frame(matrix(NA,nrow(data),1))
  rownames(follow_up_dist) <- rownames(follow_up)
  for(i in 1:nrow(follow_up))
  {
    T0 <- follow_up[i,1]
    T1 <- follow_up[i,2]
    if(!is.na(T1))
    {
      follow_up_dist[i,1] <- as.numeric(vegdist(data[c(T1,T0),],method="bray",na.rm=TRUE))
    }
  }
  return(follow_up_dist)
}

aitchison_followup <- function(data,follow_up)
{
  data_clr <- as.data.frame(as.matrix(clr(data+0.00001)))
  data_clr <- as.data.frame(t(apply(data_clr,1,function(x)(x-min(x)))))
  follow_up_dist <- as.data.frame(matrix(NA,nrow(data),1))
  rownames(follow_up_dist) <- rownames(follow_up)
  for(i in 1:nrow(follow_up))
  {
    T0 <- follow_up[i,1]
    T1 <- follow_up[i,2]
    if(!is.na(T1))
    {
      follow_up_dist[i,1] <- as.numeric(vegdist(data_clr[c(T1,T0),],method="euclidean",na.rm=TRUE))
    }
  }
  return(follow_up_dist)
}



ZhangT_2020_bray_follow_up <- bray_followup(ZhangT_2020_SpDf[rownames(ZhangT_2020_metadata_followup),sHACK_top_species],ZhangT_2020_metadata_followup)
ZhangT_2020_aitchison_follow_up <- aitchison_followup(ZhangT_2020_SpDf[rownames(ZhangT_2020_metadata_followup),sHACK_top_species],ZhangT_2020_metadata_followup)

SchmidtT_2019B_bray_follow_up <- bray_followup(SchmidtT_2019B_SpDf[rownames(SchmidtT_2019B_metadata_followup),sHACK_top_species],SchmidtT_2019B_metadata_followup)
SchmidtT_2019B_aitchison_follow_up <- aitchison_followup(SchmidtT_2019B_SpDf[rownames(SchmidtT_2019B_metadata_followup),sHACK_top_species],SchmidtT_2019B_metadata_followup)

ChenC_2018_bray_follow_up <- bray_followup(ChenC_2018_SpDf[rownames(ChenC_2018_metadata_followup),sHACK_top_species],ChenC_2018_metadata_followup)
ChenC_2018_aitchison_follow_up <- aitchison_followup(ChenC_2018_SpDf[rownames(ChenC_2018_metadata_followup),sHACK_top_species],ChenC_2018_metadata_followup)

LincardonT_2023_bray_follow_up <- bray_followup(combined_Sp_Df_validation[rownames(LincardonT_2023_metadata_followup),sHACK_top_species],LincardonT_2023_metadata_followup)
LincardonT_2023_aitchison_follow_up <- aitchison_followup(combined_Sp_Df_validation[rownames(LincardonT_2023_metadata_followup),sHACK_top_species],LincardonT_2023_metadata_followup)

Followup_BrayDist <- rbind(
  ZhangT_2020_bray_follow_up,
  SchmidtT_2019B_bray_follow_up,
  ChenC_2018_bray_follow_up,
  LincardonT_2023_bray_follow_up)

Followup_AitchisonDist <- rbind(
  ZhangT_2020_aitchison_follow_up,
  SchmidtT_2019B_aitchison_follow_up,
  ChenC_2018_aitchison_follow_up,
  LincardonT_2023_aitchison_follow_up)

  
combined_Scores_Validation$FollowUp_Bray <- Followup_BrayDist$V1[match(rownames(combined_Scores_Validation), rownames(Followup_BrayDist))]
combined_Scores_Validation$FollowUp_Aitchison <- Followup_AitchisonDist$V1[match(rownames(combined_Scores_Validation), rownames(Followup_AitchisonDist))]

refined_stability_Scores_Validation <- combined_Scores_Validation[!is.na(combined_Scores_Validation$FollowUp_Bray) &  !is.na(combined_Scores_Validation$FollowUp_Aitchison), ]

## Now do meta-analysis
library(metafor)
library(pcaPP)
compute_meta_corr <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(NA,length(grouping_list),3))
	colnames(temp_meta) <- c("dataset","ri","ni")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		dat1 <- data[data[,grouping_variable]==group,var1]
		dat2 <- data[data[,grouping_variable]==group,var2]
		temp_meta[i,2] <- cor.fk(dat1,dat2)
		temp_meta[i,3] <- length(dat1)
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	temp_meta <- escalc(measure="ZCOR",ri=ri,ni=ni,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	return(res)
}

StaBray_OMWI <- compute_meta_corr(refined_stability_Scores_Validation[!(is.na(refined_stability_Scores_Validation$FollowUp_Bray)),],"temp_OMWI","FollowUp_Bray","study_name",unique(refined_stability_Scores_Validation$study_name))

## LincardonT_2023 dont have dysbiosis scores so remove them and tehn do others.

# refined_stability_Scores_Validation2 <- refined_stability_Scores_Validation[refined_stability_Scores_Validation$study_name != "LicandroH_2023", ]

StaBray_OMWI <- compute_meta_corr(refined_stability_Scores_Validation[!(is.na(refined_stability_Scores_Validation$FollowUp_Bray)),],"temp_OMWI","FollowUp_Bray","study_name",unique(refined_stability_Scores_Validation$study_name))
StaBray_DysBray <- compute_meta_corr(refined_stability_Scores_Validation[!(is.na(refined_stability_Scores_Validation$FollowUp_Bray)),],"Dysbiosis_score_bray","FollowUp_Bray","study_name",unique(refined_stability_Scores_Validation$study_name))
StaBray_DysAit <- compute_meta_corr(refined_stability_Scores_Validation[!(is.na(refined_stability_Scores_Validation$FollowUp_Bray)),],"Dysbiosis_score_aitchison","FollowUp_Bray","study_name",unique(refined_stability_Scores_Validation$study_name))
StaBray_DysCan <- compute_meta_corr(refined_stability_Scores_Validation[!(is.na(refined_stability_Scores_Validation$FollowUp_Bray)),],"Dysbiosis_score_canberra","FollowUp_Bray","study_name",unique(refined_stability_Scores_Validation$study_name))


## Now create a df of this combined comparisons
Sta_All <- data.frame(
  Score = c("WISH_Index","DysbiosisBray","DysbiosisAitchison","DysbiosisCanberra"),
  Effect = c(as.numeric(StaBray_OMWI$b),as.numeric(StaBray_DysBray$b),as.numeric(StaBray_DysAit$b),as.numeric(StaBray_DysCan$b)),
  Pval = c(StaBray_OMWI$pval,StaBray_DysBray$pval,StaBray_DysAit$pval,StaBray_DysCan$pval))

Sta_All$logP <- -log10(Sta_All$Pval)


Sta_All$LabelColor <- ifelse(Sta_All$Pval < 0.05 & Sta_All$Effect > 0, "Positive",
  ifelse(Sta_All$Pval < 0.05 & Sta_All$Effect < 0, "Negative", "NonSig"))

library(ggplot2)
library(ggrepel)



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_4Validation_MetaAnalysis_Stability_sHACKScore_DysbiosisScore.pdf", width = 8, height = 6)
ggplot(Sta_All,
       aes(x = Effect, y = logP, color = LabelColor)) +
  
  geom_point(size = 3) +
  
  geom_text_repel(
    aes(label = Score),
    size = 5,
    max.overlaps = 45
  ) +
  
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), 
             color = "red", 
             linetype = "dashed", 
             linewidth = 0.6) +

  scale_color_manual(values = c(
    "Positive" = "#4575b4",
    "Negative" = "#d73027",
    "NonSig"   = "grey60"
  )) +
  
  theme_minimal() +
  labs(title = "Comparison of Microbiome Indices",
       x = "Effect Size",
       y = "-log10(P-value)") +
  
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.position = "none"
  )
dev.off()




save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_4sHACK_Index_Validation_Workspace.RData")

########### save the combined scores, combined SpDf and combined saliva score in one RData file
save(combined_Scores_Validation, combined_Sp_Df_validation, Combined_Saliva_Scores, MetadataDf_saliva_validation,file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_4sHACK_Index_Validation_DataDf.RData")




################### See if this 8 species can give the difference between disease adn control samples in each of the study.

ChenC_2018_SpDf