

########### Health Association for Validation Cohort
S9_1CoreAssociation_Validation_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_Validation_Workspace.RData", envir = S9_1CoreAssociation_Validation_Workspace)
attach(S9_1CoreAssociation_Validation_Workspace)
MetadataDf_saliva_validation <- MetadataDf_saliva_validation
SpDf_saliva_Validation <- SpDf_saliva_Validation
Combined_Saliva_Scores <- Combined_Saliva_Scores
detach(S9_1CoreAssociation_Validation_Workspace)
rm(S9_1CoreAssociation_Validation_Workspace)

all(rownames(MetadataDf_saliva_validation)==rownames(SpDf_saliva_Validation))
SpDf_saliva_Validation$study_name <- MetadataDf_saliva_validation$study_name

study_list <- unique(SpDf_saliva_Validation$study_name)

# FanX_2018a (HNC - Control - USA - 16s)
# ZhangT_2020 (IBD/remission - Control - 16s - CHN - Long)
# ChenC_2018 (Periodontitis (Pre-Post) - Control - USA - 16s)
# SchmidtT_2019B (T1D - Control - LUX - WGS - Longitudinal)
# ChenX_2025 (Hypertension - Control - 16s)
# TakayanagiK_2023 (Moyamoya Disease - Control - JPN - 16s)
# WuZ_2025  (Metastatic CRC, Non-metastatic CRC - Control - CHN - 16s)
# ZhangL_2023 (CRC Polyps - Control - CHN - 16s)
# LicandroH_2023 (All Control - USA - 16s - Longitudinal)
# "ChenJ_2021"  (OSCC/Oral_Verrucous_Hyperplasia - Control - 16s - CHN)     
# "CirsteaM_2022"  (AD - Control - 16s - CAN - seniors)
# "FinkelsteinS_2025" (GDM - Control - 16s - Israel)
# "IglesiasA_2024" (Periodontitis - Control - 16s - ESP/PRT)
# "LiuY_2021"  (T2D - Control - 16s - CHN)       
# "NearingJ_2023"  (Breast Cancer/Colon_Cancer/Prostate Cancer - Control - 16s - CAN)
# "NiuC_2020"   (All Control - CHN - 16s - Longitudinal)      
# "RelvasM_2021"  (Dental/Periodontal - Control - 16s - PRT - Severity)    
# "StahringerS_2012" (All Control - USA -16s - Longitudinal - 2012 Very old data)
# "JiY_2020" (helicobacter pylori uninfected/infected/eradicated - Control - 16s - CHN)


############### Divide the data into individual studies.


ChenC_2018_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "ChenC_2018",]
ChenC_2018_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "ChenC_2018",]
ChenC_2018_ctrl <- ChenC_2018_metadata[ChenC_2018_metadata$study_condition == "Control","sample_id"]
ChenC_2018_dis <- ChenC_2018_metadata[ChenC_2018_metadata$study_condition == "Diseased","sample_id"]
ChenC_2018_postDis <- ChenC_2018_metadata[ChenC_2018_metadata$study_condition == "Post-Disease","sample_id"]


FanX_2018a_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "FanX_2018a",]
FanX_2018a_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "FanX_2018a",]
# FanX_2018a_ctrl <- FanX_2018a_metadata[FanX_2018a_metadata$study_condition == "Control","sample_id"]
# FanX_2018a_dis <- FanX_2018a_metadata[FanX_2018a_metadata$study_condition == "Diseased","sample_id"]

FanX_2018a_ctrl_Al <- FanX_2018a_metadata[FanX_2018a_metadata$study_condition == "Control" & FanX_2018a_metadata$exposure == "hnc-control,alcoholic","sample_id"]
FanX_2018a_dis_Al <- FanX_2018a_metadata[FanX_2018a_metadata$study_condition == "Diseased" & FanX_2018a_metadata$exposure == "alcoholic","sample_id"]
# 
FanX_2018a_ctrl_NonAl <- FanX_2018a_metadata[FanX_2018a_metadata$study_condition == "Control" & FanX_2018a_metadata$exposure == "hnc-control,non-alcoholic" ,"sample_id"]
FanX_2018a_dis_NonAl <- FanX_2018a_metadata[FanX_2018a_metadata$study_condition == "Diseased" & FanX_2018a_metadata$exposure == "non-alcoholic","sample_id"]

ZhangT_2020_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "ZhangT_2020",]
ZhangT_2020_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "ZhangT_2020",]
ZhangT_2020_ctrl <- ZhangT_2020_metadata[ZhangT_2020_metadata$study_condition == "Control","sample_id"]
ZhangT_2020_dis <- ZhangT_2020_metadata[ZhangT_2020_metadata$study_condition == "Diseased","sample_id"]
ZhangT_2020_dis_rem <- ZhangT_2020_metadata[ZhangT_2020_metadata$study_condition == "Remission IBD_GutInflammation","sample_id"]


SchmidtT_2019B_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "SchmidtT_2019B",]
SchmidtT_2019B_metadata <- SchmidtT_2019B_metadata[!SchmidtT_2019B_metadata$disease %in% c("T2D"),]
SchmidtT_2019B_metadata$study_condition[SchmidtT_2019B_metadata$study_condition == "Diseased"] <- "Diseased"

SchmidtT_2019B_SpDf <- SpDf_saliva_Validation[rownames(SpDf_saliva_Validation) %in% rownames(SchmidtT_2019B_metadata),]
SchmidtT_2019B_ctrl <- SchmidtT_2019B_metadata[SchmidtT_2019B_metadata$study_condition == "Control","sample_id"]
SchmidtT_2019B_dis <- SchmidtT_2019B_metadata[SchmidtT_2019B_metadata$study_condition == "Diseased" & SchmidtT_2019B_metadata$disease == "T1D","sample_id"]


ChenX_2025_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "ChenX_2025",]
ChenX_2025_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "ChenX_2025",]
ChenX_2025_ctrl <- ChenX_2025_metadata[ChenX_2025_metadata$study_condition == "Control","sample_id"]
ChenX_2025_dis <- ChenX_2025_metadata[ChenX_2025_metadata$study_condition == "Diseased","sample_id"]

TakayanagiK_2023_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "TakayanagiK_2023",]
TakayanagiK_2023_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "TakayanagiK_2023",]
TakayanagiK_2023_ctrl <- TakayanagiK_2023_metadata[TakayanagiK_2023_metadata$study_condition == "Control","sample_id"]
TakayanagiK_2023_dis <- TakayanagiK_2023_metadata[TakayanagiK_2023_metadata$study_condition == "Diseased","sample_id"]

WuZ_2025_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "WuZ_2025",]
WuZ_2025_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "WuZ_2025",]
WuZ_2025_ctrl <- WuZ_2025_metadata[WuZ_2025_metadata$study_condition == "Control","sample_id"]
WuZ_2025_dis <- WuZ_2025_metadata[WuZ_2025_metadata$study_condition == "Diseased" & WuZ_2025_metadata$disease == "Metastatic_CRC","sample_id"]
WuZ_2025_dis_nonmet <- WuZ_2025_metadata[WuZ_2025_metadata$study_condition == "Diseased" & WuZ_2025_metadata$disease == "NonMetastatic_CRC","sample_id"]

ZhangL_2023_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "ZhangL_2023",]
ZhangL_2023_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "ZhangL_2023",]
ZhangL_2023_ctrl <- ZhangL_2023_metadata[ZhangL_2023_metadata$study_condition == "Control","sample_id"]
ZhangL_2023_dis <- ZhangL_2023_metadata[ZhangL_2023_metadata$study_condition == "Diseased","sample_id"]


ChenJ_2021_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "ChenJ_2021",]
ChenJ_2021_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "ChenJ_2021",]
ChenJ_2021_ctrl <- ChenJ_2021_metadata[ChenJ_2021_metadata$study_condition == "Control","sample_id"]
ChenJ_2021_dis_OSCC <- ChenJ_2021_metadata[ChenJ_2021_metadata$study_condition == "Diseased" & ChenJ_2021_metadata$disease == "OSCC","sample_id"]
ChenJ_2021_dis_OVH <- ChenJ_2021_metadata[ChenJ_2021_metadata$study_condition == "Diseased" & ChenJ_2021_metadata$disease == "Oral_Verrucous_Hyperplasia","sample_id"]

CirsteaM_2022_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "CirsteaM_2022",]
CirsteaM_2022_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "CirsteaM_2022",]
CirsteaM_2022_ctrl <- CirsteaM_2022_metadata[CirsteaM_2022_metadata$study_condition == "Control" & CirsteaM_2022_metadata$age_category == "senior","sample_id"]
CirsteaM_2022_dis <- CirsteaM_2022_metadata[CirsteaM_2022_metadata$study_condition == "Diseased" & CirsteaM_2022_metadata$age_category == "senior","sample_id"]

FinkelsteinS_2025_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "FinkelsteinS_2025",] 
FinkelsteinS_2025_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "FinkelsteinS_2025",] 
FinkelsteinS_2025_ctrl <- FinkelsteinS_2025_metadata[FinkelsteinS_2025_metadata$study_condition == "Control","sample_id"]
FinkelsteinS_2025_dis <- FinkelsteinS_2025_metadata[FinkelsteinS_2025_metadata$study_condition == "Diseased","sample_id"]

IglesiasA_2024_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "IglesiasA_2024",] 
IglesiasA_2024_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "IglesiasA_2024",]
IglesiasA_2024_ctrl <- IglesiasA_2024_metadata[IglesiasA_2024_metadata$study_condition == "Control","sample_id"]
IglesiasA_2024_dis <- IglesiasA_2024_metadata[IglesiasA_2024_metadata$study_condition == "Diseased","sample_id"]

LiuY_2021_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "LiuY_2021",]
LiuY_2021_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "LiuY_2021",]
LiuY_2021_ctrl <- LiuY_2021_metadata[LiuY_2021_metadata$study_condition == "Control","sample_id"]
LiuY_2021_dis <- LiuY_2021_metadata[LiuY_2021_metadata$study_condition == "Diseased","sample_id"]

NearingJ_2023_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "NearingJ_2023",] 
NearingJ_2023_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "NearingJ_2023",]
NearingJ_2023_ctrl_BC <- NearingJ_2023_metadata[NearingJ_2023_metadata$study_condition == "Control" & NearingJ_2023_metadata$original_disease == "Control_BC","sample_id"]
NearingJ_2023_ctrl_CC <- NearingJ_2023_metadata[NearingJ_2023_metadata$study_condition == "Control" & NearingJ_2023_metadata$original_disease == "Control_CC","sample_id"]
NearingJ_2023_ctrl_PC <- NearingJ_2023_metadata[NearingJ_2023_metadata$study_condition == "Control" & NearingJ_2023_metadata$original_disease == "Control_PC","sample_id"]
NearingJ_2023_dis_BC <- NearingJ_2023_metadata[NearingJ_2023_metadata$study_condition == "Diseased" & NearingJ_2023_metadata$disease == "Breast_Cancer","sample_id"]
NearingJ_2023_dis_CC <- NearingJ_2023_metadata[NearingJ_2023_metadata$study_condition == "Diseased" & NearingJ_2023_metadata$disease == "CRC","sample_id"]
NearingJ_2023_dis_PC <- NearingJ_2023_metadata[NearingJ_2023_metadata$study_condition == "Diseased" & NearingJ_2023_metadata$disease == "Prostate_Cancer","sample_id"]


RelvasM_2021_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "RelvasM_2021",] 
RelvasM_2021_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "RelvasM_2021",]
RelvasM_2021_ctrl <- RelvasM_2021_metadata[RelvasM_2021_metadata$study_condition == "Control","sample_id"]
RelvasM_2021_dis_Dg1 <- RelvasM_2021_metadata[RelvasM_2021_metadata$study_condition == "Diseased" & RelvasM_2021_metadata$Grades_Dental_Health == "Dental_Grade_1","sample_id"]
RelvasM_2021_dis_Dg2 <- RelvasM_2021_metadata[RelvasM_2021_metadata$study_condition == "Diseased" & RelvasM_2021_metadata$Grades_Dental_Health == "Dental_Grade_2","sample_id"]
RelvasM_2021_dis_Dg3 <- RelvasM_2021_metadata[RelvasM_2021_metadata$study_condition == "Diseased" & RelvasM_2021_metadata$Grades_Dental_Health == "Dental_Grade_3","sample_id"]
RelvasM_2021_dis_Pg1 <- RelvasM_2021_metadata[RelvasM_2021_metadata$study_condition == "Diseased" & RelvasM_2021_metadata$Grades_Periodontal_Health == "Perio_Grade_1","sample_id"]
RelvasM_2021_dis_Pg2 <- RelvasM_2021_metadata[RelvasM_2021_metadata$study_condition == "Diseased" & RelvasM_2021_metadata$Grades_Periodontal_Health == "Perio_Grade_2","sample_id"]
RelvasM_2021_dis_Pg3 <- RelvasM_2021_metadata[RelvasM_2021_metadata$study_condition == "Diseased" & RelvasM_2021_metadata$Grades_Periodontal_Health == "Perio_Grade_3","sample_id"]


JiY_2020_metadata <- MetadataDf_saliva_validation[MetadataDf_saliva_validation$study_name == "JiY_2020",]
JiY_2020_SpDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == "JiY_2020",]
JiY_2020_ctrl <- JiY_2020_metadata[JiY_2020_metadata$exposure == "helicobacter pylori uninfected","sample_id"]
JiY_2020_dis <- JiY_2020_metadata[JiY_2020_metadata$exposure == "helicobacter pylori infected","sample_id"]
JiY_2020_post_dis <- JiY_2020_metadata[JiY_2020_metadata$exposure == "helicobacter pylori eradicated","sample_id"]







# LicandroH_2023, NiuC_2020, StahringerS_2012 has all controls, so we will not be able to do a case-control comparison for this studies. We will skip this study for the health association analysis.

############### Calculate the direction for each of the species in each of the study.

wilcox_batch = function(x,y)
{
  p_array <- NULL;
  type_array <- NULL;
  mean1_array <- NULL;
  mean2_array <- NULL;
  x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
  y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
  z <- intersect(rownames(x),rownames(y));
  for(i in 1:length(z))
  {
    xv <- as.numeric(x[z[i], ])
    yv <- as.numeric(y[z[i], ])
    
    xv <- xv[!is.na(xv)]
    yv <- yv[!is.na(yv)]
    
    if (length(xv) >= 2 && length(yv) >= 2) {
      wt <- wilcox.test(xv, yv)
      p_array[i] <- wt$p.value
    } else {
      p_array[i] <- 1
    }
    
    type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
    mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
    i <- i + 1;
  }
  out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array,method="fdr"),mean1_array,mean2_array));
  rownames(out) <- z;
  out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
  return(t(out));
}


detected_species_list <- list()

for(i in study_list){
  temp_spDf <- SpDf_saliva_Validation[SpDf_saliva_Validation$study_name == i,]
  temp_spDf$study_name <- NULL
  temp_spDf <- temp_spDf[,colSums(temp_spDf)>0]
  temp_species <- colnames(temp_spDf)
  detected_species_list[[i]] <- temp_species
}


### ChenC_2018
wilcox_ChenC_2018 <- wilcox_batch(
  t(ChenC_2018_SpDf[ChenC_2018_dis, detected_species_list[["ChenC_2018"]]]),
  t(ChenC_2018_SpDf[ChenC_2018_ctrl, detected_species_list[["ChenC_2018"]]]))

df_Compare_ChenC_2018 <- data.frame(
  Direction = ifelse(wilcox_ChenC_2018[,1] <= 0.05, sign(wilcox_ChenC_2018[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_ChenC_2018), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_ChenC_2018), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_ChenC_2018), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_ChenC_2018), 5])

### ChenC_2018_postDis
wilcox_ChenC_2018_postDis <- wilcox_batch(
  t(ChenC_2018_SpDf[ChenC_2018_dis, detected_species_list[["ChenC_2018"]]]),
  t(ChenC_2018_SpDf[ChenC_2018_postDis, detected_species_list[["ChenC_2018"]]]))

df_Compare_ChenC_2018_postDis <- data.frame(
  Direction = ifelse(wilcox_ChenC_2018_postDis[,1] <= 0.05, sign(wilcox_ChenC_2018_postDis[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_ChenC_2018_postDis), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_ChenC_2018_postDis), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_ChenC_2018_postDis), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_ChenC_2018_postDis), 5])



### FanX_2018a
wilcox_FanX_2018a <- wilcox_batch(
  t(FanX_2018a_SpDf[FanX_2018a_dis_NonAl, detected_species_list[["FanX_2018a"]]]),
  t(FanX_2018a_SpDf[FanX_2018a_ctrl_NonAl, detected_species_list[["FanX_2018a"]]]))

df_Compare_FanX_2018a <- data.frame(
  Direction = ifelse(wilcox_FanX_2018a[,1] <= 0.05, sign(wilcox_FanX_2018a[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_FanX_2018a), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_FanX_2018a), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_FanX_2018a), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_FanX_2018a), 5])

### FanX_2018a_Al
wilcox_FanX_2018a_Al <- wilcox_batch(
  t(FanX_2018a_SpDf[FanX_2018a_dis_Al, detected_species_list[["FanX_2018a"]]]),
  t(FanX_2018a_SpDf[FanX_2018a_ctrl_Al, detected_species_list[["FanX_2018a"]]]))

df_Compare_FanX_2018a_Al <- data.frame(
  Direction = ifelse(wilcox_FanX_2018a_Al[,1] <= 0.05, sign(wilcox_FanX_2018a_Al[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_FanX_2018a_Al), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_FanX_2018a_Al), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_FanX_2018a_Al), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_FanX_2018a_Al), 5])


### ZhangT_2020
wilcox_ZhangT_2020 <- wilcox_batch(
  t(ZhangT_2020_SpDf[ZhangT_2020_dis, detected_species_list[["ZhangT_2020"]]]),
  t(ZhangT_2020_SpDf[ZhangT_2020_ctrl, detected_species_list[["ZhangT_2020"]]]))

df_Compare_ZhangT_2020 <- data.frame(
  Direction = ifelse(wilcox_ZhangT_2020[,1] <= 0.05, sign(wilcox_ZhangT_2020[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_ZhangT_2020), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_ZhangT_2020), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_ZhangT_2020), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_ZhangT_2020), 5])


### ZhangT_2020Remission
wilcox_ZhangT_2020R <- wilcox_batch(
  t(ZhangT_2020_SpDf[ZhangT_2020_dis, detected_species_list[["ZhangT_2020"]]]),
  t(ZhangT_2020_SpDf[ZhangT_2020_dis_rem, detected_species_list[["ZhangT_2020"]]]))

df_Compare_ZhangT_2020R <- data.frame(
  Direction = ifelse(wilcox_ZhangT_2020R[,1] <= 0.05, sign(wilcox_ZhangT_2020R[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_ZhangT_2020R), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_ZhangT_2020R), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_ZhangT_2020R), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_ZhangT_2020R), 5])



### SchmidtT_2019B
wilcox_SchmidtT_2019B <- wilcox_batch(
  t(SchmidtT_2019B_SpDf[SchmidtT_2019B_dis, detected_species_list[["SchmidtT_2019B"]]]),
  t(SchmidtT_2019B_SpDf[SchmidtT_2019B_ctrl, detected_species_list[["SchmidtT_2019B"]]]))

df_Compare_SchmidtT_2019B <- data.frame(
  Direction = ifelse(wilcox_SchmidtT_2019B[,1] <= 0.05, sign(wilcox_SchmidtT_2019B[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_SchmidtT_2019B), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_SchmidtT_2019B), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_SchmidtT_2019B), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_SchmidtT_2019B), 5])


## ChenX_2025
wilcox_ChenX_2025 <- wilcox_batch(
  t(ChenX_2025_SpDf[ChenX_2025_dis, detected_species_list[["ChenX_2025"]]]),
  t(ChenX_2025_SpDf[ChenX_2025_ctrl, detected_species_list[["ChenX_2025"]]]))

df_Compare_ChenX_2025 <- data.frame(
  Direction = ifelse(wilcox_ChenX_2025[,1] <= 0.05, sign(wilcox_ChenX_2025[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_ChenX_2025), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_ChenX_2025), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_ChenX_2025), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_ChenX_2025), 5])


## TakayanagiK_2023
wilcox_TakayanagiK_2023 <- wilcox_batch(
  t(TakayanagiK_2023_SpDf[TakayanagiK_2023_dis, detected_species_list[["TakayanagiK_2023"]]]),
  t(TakayanagiK_2023_SpDf[TakayanagiK_2023_ctrl, detected_species_list[["TakayanagiK_2023"]]])) 

df_Compare_TakayanagiK_2023 <- data.frame(
  Direction = ifelse(wilcox_TakayanagiK_2023[,1] <= 0.05, sign(wilcox_TakayanagiK_2023[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_TakayanagiK_2023), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_TakayanagiK_2023), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_TakayanagiK_2023), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_TakayanagiK_2023), 5])

## WuZ_2025
wilcox_WuZ_2025 <- wilcox_batch(
  t(WuZ_2025_SpDf[WuZ_2025_dis, detected_species_list[["WuZ_2025"]]]),
  t(WuZ_2025_SpDf[WuZ_2025_ctrl, detected_species_list[["WuZ_2025"]]]))

df_Compare_WuZ_2025 <- data.frame(
  Direction = ifelse(wilcox_WuZ_2025[,1] <= 0.05, sign(wilcox_WuZ_2025[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025), 5])

## WuZ_2025_nonmet
wilcox_WuZ_2025_nonmet <- wilcox_batch(
  t(WuZ_2025_SpDf[WuZ_2025_dis_nonmet, detected_species_list[["WuZ_2025"]]]),
  t(WuZ_2025_SpDf[WuZ_2025_ctrl, detected_species_list[["WuZ_2025"]]]))

df_Compare_WuZ_2025_nonmet <- data.frame(
  Direction = ifelse(wilcox_WuZ_2025_nonmet[,1] <= 0.05, sign(wilcox_WuZ_2025_nonmet[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025_nonmet), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025_nonmet), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025_nonmet), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025_nonmet), 5])  

## WuZ_2025_met_vs_nonmet
wilcox_WuZ_2025_met_nonmet <- wilcox_batch(
  t(WuZ_2025_SpDf[WuZ_2025_dis, detected_species_list[["WuZ_2025"]]]),
  t(WuZ_2025_SpDf[WuZ_2025_dis_nonmet, detected_species_list[["WuZ_2025"]]]))

df_Compare_WuZ_2025_met_nonmet <- data.frame(
  Direction = ifelse(wilcox_WuZ_2025_met_nonmet[,1] <= 0.05, sign(wilcox_WuZ_2025_met_nonmet[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025_met_nonmet), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025_met_nonmet), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025_met_nonmet), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_WuZ_2025_met_nonmet), 5])


## ZhangL_2023
wilcox_ZhangL_2023 <- wilcox_batch(
  t(ZhangL_2023_SpDf[ZhangL_2023_dis, detected_species_list[["ZhangL_2023"]]]),
  t(ZhangL_2023_SpDf[ZhangL_2023_ctrl, detected_species_list[["ZhangL_2023"]]]))

df_Compare_ZhangL_2023 <- data.frame(
  Direction = ifelse(wilcox_ZhangL_2023[,1] <= 0.05, sign(wilcox_ZhangL_2023[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_ZhangL_2023), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_ZhangL_2023), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_ZhangL_2023), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_ZhangL_2023), 5])


## ChenJ_2021 OSCC
wilcox_ChenJ_2021_OSCC <- wilcox_batch(
  t(ChenJ_2021_SpDf[ChenJ_2021_dis_OSCC, detected_species_list[["ChenJ_2021"]]]),
  t(ChenJ_2021_SpDf[ChenJ_2021_ctrl, detected_species_list[["ChenJ_2021"]]]))

df_Compare_ChenJ_2021_OSCC <- data.frame(
  Direction = ifelse(wilcox_ChenJ_2021_OSCC[,1] <= 0.05, sign(wilcox_ChenJ_2021_OSCC[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_ChenJ_2021_OSCC), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_ChenJ_2021_OSCC), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_ChenJ_2021_OSCC), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_ChenJ_2021_OSCC), 5])

## ChenJ_2021 OVH
wilcox_ChenJ_2021_OVH <- wilcox_batch(
  t(ChenJ_2021_SpDf[ChenJ_2021_dis_OVH, detected_species_list[["ChenJ_2021"]]]),
  t(ChenJ_2021_SpDf[ChenJ_2021_ctrl, detected_species_list[["ChenJ_2021"]]]))

df_Compare_ChenJ_2021_OVH <- data.frame(
  Direction = ifelse(wilcox_ChenJ_2021_OVH[,1] <= 0.05, sign(wilcox_ChenJ_2021_OVH[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_ChenJ_2021_OVH), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_ChenJ_2021_OVH), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_ChenJ_2021_OVH), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_ChenJ_2021_OVH), 5])

## CirsteaM_2022
wilcox_CirsteaM_2022 <- wilcox_batch(
  t(CirsteaM_2022_SpDf[CirsteaM_2022_dis, detected_species_list[["CirsteaM_2022"]]]),
  t(CirsteaM_2022_SpDf[CirsteaM_2022_ctrl, detected_species_list[["CirsteaM_2022"]]]))

df_Compare_CirsteaM_2022 <- data.frame(
  Direction = ifelse(wilcox_CirsteaM_2022[,1] <= 0.05, sign(wilcox_CirsteaM_2022[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_CirsteaM_2022), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_CirsteaM_2022), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_CirsteaM_2022), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_CirsteaM_2022), 5]) 

## FinkelsteinS_2025
wilcox_FinkelsteinS_2025 <- wilcox_batch(
  t(FinkelsteinS_2025_SpDf[FinkelsteinS_2025_dis, detected_species_list[["FinkelsteinS_2025"]]]),
  t(FinkelsteinS_2025_SpDf[FinkelsteinS_2025_ctrl, detected_species_list[["FinkelsteinS_2025"]]]))

df_Compare_FinkelsteinS_2025 <- data.frame(
  Direction = ifelse(wilcox_FinkelsteinS_2025[,1] <= 0.05, sign(wilcox_FinkelsteinS_2025[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_FinkelsteinS_2025), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_FinkelsteinS_2025), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_FinkelsteinS_2025), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_FinkelsteinS_2025), 5])


## IglesiasA_2024
wilcox_IglesiasA_2024 <- wilcox_batch(
  t(IglesiasA_2024_SpDf[IglesiasA_2024_dis, detected_species_list[["IglesiasA_2024"]]]),
  t(IglesiasA_2024_SpDf[IglesiasA_2024_ctrl, detected_species_list[["IglesiasA_2024"]]]))

df_Compare_IglesiasA_2024 <- data.frame(
  Direction = ifelse(wilcox_IglesiasA_2024[,1] <= 0.05, sign(wilcox_IglesiasA_2024[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_IglesiasA_2024), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_IglesiasA_2024), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_IglesiasA_2024), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_IglesiasA_2024), 5])

## NearingJ_2023_BC
wilcox_NearingJ_2023_BC <- wilcox_batch(
  t(NearingJ_2023_SpDf[NearingJ_2023_dis_BC, detected_species_list[["NearingJ_2023"]]]),
  t(NearingJ_2023_SpDf[NearingJ_2023_ctrl_BC, detected_species_list[["NearingJ_2023"]]]))

df_Compare_NearingJ_2023_BC <- data.frame(
  Direction = ifelse(wilcox_NearingJ_2023_BC[,1] <= 0.05, sign(wilcox_NearingJ_2023_BC[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_BC), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_BC), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_BC), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_BC), 5])

## NearingJ_2023_CC
wilcox_NearingJ_2023_CC <- wilcox_batch(
  t(NearingJ_2023_SpDf[NearingJ_2023_dis_CC, detected_species_list[["NearingJ_2023"]]]),
  t(NearingJ_2023_SpDf[NearingJ_2023_ctrl_CC, detected_species_list[["NearingJ_2023"]]]))

df_Compare_NearingJ_2023_CC <- data.frame(
  Direction = ifelse(wilcox_NearingJ_2023_CC[,1] <= 0.05, sign(wilcox_NearingJ_2023_CC[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_CC), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_CC), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_CC), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_CC), 5])

## NearingJ_2023_PC
wilcox_NearingJ_2023_PC <- wilcox_batch(
  t(NearingJ_2023_SpDf[NearingJ_2023_dis_PC, detected_species_list[["NearingJ_2023"]]]),
  t(NearingJ_2023_SpDf[NearingJ_2023_ctrl_PC, detected_species_list[["NearingJ_2023"]]]))  

df_Compare_NearingJ_2023_PC <- data.frame(
  Direction = ifelse(wilcox_NearingJ_2023_PC[,1] <= 0.05, sign(wilcox_NearingJ_2023_PC[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_PC), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_PC), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_PC), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_NearingJ_2023_PC), 5]) 

## RelvasM_2021_Dental_Grade_1
wilcox_RelvasM_2021_Dg1 <- wilcox_batch(
  t(RelvasM_2021_SpDf[RelvasM_2021_dis_Dg1, detected_species_list[["RelvasM_2021"]]]),
  t(RelvasM_2021_SpDf[RelvasM_2021_ctrl, detected_species_list[["RelvasM_2021"]]])) 

df_Compare_RelvasM_2021_Dg1 <- data.frame(
  Direction = ifelse(wilcox_RelvasM_2021_Dg1[,1] <= 0.05, sign(wilcox_RelvasM_2021_Dg1[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg1), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg1), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg1), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg1), 5])

## RelvasM_2021_Dental_Grade_2
wilcox_RelvasM_2021_Dg2 <- wilcox_batch(
  t(RelvasM_2021_SpDf[RelvasM_2021_dis_Dg2, detected_species_list[["RelvasM_2021"]]]),
  t(RelvasM_2021_SpDf[RelvasM_2021_ctrl, detected_species_list[["RelvasM_2021"]]]))

df_Compare_RelvasM_2021_Dg2 <- data.frame(
  Direction = ifelse(wilcox_RelvasM_2021_Dg2[,1] <= 0.05, sign(wilcox_RelvasM_2021_Dg2[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg2), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg2), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg2), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg2), 5])

## RelvasM_2021_Dental_Grade_3
wilcox_RelvasM_2021_Dg3 <- wilcox_batch(
  t(RelvasM_2021_SpDf[RelvasM_2021_dis_Dg3, detected_species_list[["RelvasM_2021"]]]),
  t(RelvasM_2021_SpDf[RelvasM_2021_ctrl, detected_species_list[["RelvasM_2021"]]])) 

df_Compare_RelvasM_2021_Dg3 <- data.frame(
  Direction = ifelse(wilcox_RelvasM_2021_Dg3[,1] <= 0.05, sign(wilcox_RelvasM_2021_Dg3[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg3), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg3), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg3), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Dg3), 5])

## RelvasM_2021_Periodontal_Grade_1
wilcox_RelvasM_2021_Pg1 <- wilcox_batch(
  t(RelvasM_2021_SpDf[RelvasM_2021_dis_Pg1, detected_species_list[["RelvasM_2021"]]]),
  t(RelvasM_2021_SpDf[RelvasM_2021_ctrl, detected_species_list[["RelvasM_2021"]]]))

df_Compare_RelvasM_2021_Pg1 <- data.frame(
  Direction = ifelse(wilcox_RelvasM_2021_Pg1[,1] <= 0.05, sign(wilcox_RelvasM_2021_Pg1[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg1), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg1), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg1), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg1), 5])

## RelvasM_2021_Periodontal_Grade_2
wilcox_RelvasM_2021_Pg2 <- wilcox_batch(
  t(RelvasM_2021_SpDf[RelvasM_2021_dis_Pg2, detected_species_list[["RelvasM_2021"]]]),
  t(RelvasM_2021_SpDf[RelvasM_2021_ctrl, detected_species_list[["RelvasM_2021"]]]))

df_Compare_RelvasM_2021_Pg2 <- data.frame(
  Direction = ifelse(wilcox_RelvasM_2021_Pg2[,1] <= 0.05, sign(wilcox_RelvasM_2021_Pg2[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg2), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg2), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg2), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg2), 5]) 

# ## RelvasM_2021_Periodontal_Grade_3
# wilcox_RelvasM_2021_Pg3 <- wilcox_batch(
#   t(RelvasM_2021_SpDf[RelvasM_2021_dis_Pg3, detected_species_list[["RelvasM_2021"]]]),
#   t(RelvasM_2021_SpDf[RelvasM_2021_ctrl, detected_species_list[["RelvasM_2021"]]]))

# df_Compare_RelvasM_2021_Pg3 <- data.frame(
#   Direction = ifelse(wilcox_RelvasM_2021_Pg3[,1] <= 0.05, sign(wilcox_RelvasM_2021_Pg3[,2]), 0),
#   Influence = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg3), 1],
#   Stability = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg3), 3],
#   Health    = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg3), 2],
#   HACKScore = Combined_Saliva_Scores[rownames(wilcox_RelvasM_2021_Pg3), 5])

## Only one samples is in RelvasM_2021_dis_Pg3 so we are skipping df_Compare_RelvasM_2021_Pg3/wilcox_RelvasM_2021_Pg3


## JiY_2020
wilcox_JiY_2020 <- wilcox_batch(
  t(JiY_2020_SpDf[JiY_2020_dis, detected_species_list[["JiY_2020"]]]),
  t(JiY_2020_SpDf[JiY_2020_ctrl, detected_species_list[["JiY_2020"]]]))

df_Compare_JiY_2020 <- data.frame(
  Direction = ifelse(wilcox_JiY_2020[,1] <= 0.05, sign(wilcox_JiY_2020[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_JiY_2020), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_JiY_2020), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_JiY_2020), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_JiY_2020), 5])


## JiY_2020 eradicated
wilcox_JiY_2020_Er <- wilcox_batch(
  t(JiY_2020_SpDf[JiY_2020_post_dis, detected_species_list[["JiY_2020"]]]),
  t(JiY_2020_SpDf[JiY_2020_ctrl, detected_species_list[["JiY_2020"]]]))

df_Compare_JiY_2020_Er <- data.frame(
  Direction = ifelse(wilcox_JiY_2020_Er[,1] <= 0.05, sign(wilcox_JiY_2020_Er[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_JiY_2020_Er), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_JiY_2020_Er), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_JiY_2020_Er), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_JiY_2020_Er), 5])


## LiuY_2021
wilcox_LiuY_2021 <- wilcox_batch(
  t(LiuY_2021_SpDf[LiuY_2021_dis, detected_species_list[["LiuY_2021"]]]),
  t(LiuY_2021_SpDf[LiuY_2021_ctrl, detected_species_list[["LiuY_2021"]]]))

df_Compare_LiuY_2021 <- data.frame(
  Direction = ifelse(wilcox_LiuY_2021[,1] <= 0.05, sign(wilcox_LiuY_2021[,2]), 0),
  Influence = Combined_Saliva_Scores[rownames(wilcox_LiuY_2021), 1],
  Stability = Combined_Saliva_Scores[rownames(wilcox_LiuY_2021), 3],
  Health    = Combined_Saliva_Scores[rownames(wilcox_LiuY_2021), 2],
  HACKScore = Combined_Saliva_Scores[rownames(wilcox_LiuY_2021), 5])






## Now plot above data for all the studies together and also add the names for each bean.
library(beanplot)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_AllStudies.pdf",width = 20, height = 7)
beanplot(
  df_Compare_ChenC_2018[df_Compare_ChenC_2018[,1]==-1,4],df_Compare_ChenC_2018[df_Compare_ChenC_2018[,1]!=-1 | df_Compare_ChenC_2018[,1]==0,4],
  df_Compare_ChenC_2018_postDis[df_Compare_ChenC_2018_postDis[,1]==-1,4],df_Compare_ChenC_2018_postDis[df_Compare_ChenC_2018_postDis[,1]!=-1 | df_Compare_ChenC_2018_postDis[,1]==0,4],
  df_Compare_FanX_2018a[df_Compare_FanX_2018a[,1]==-1,4],df_Compare_FanX_2018a[df_Compare_FanX_2018a[,1]!=-1 | df_Compare_FanX_2018a[,1]==0,4],
  df_Compare_FanX_2018a_Al[df_Compare_FanX_2018a_Al[,1]==-1,4],df_Compare_FanX_2018a_Al[df_Compare_FanX_2018a_Al[,1]!=-1 | df_Compare_FanX_2018a_Al[,1]==0,4],
  df_Compare_ZhangT_2020[df_Compare_ZhangT_2020[,1]==-1,4],df_Compare_ZhangT_2020[df_Compare_ZhangT_2020[,1]!=-1 | df_Compare_ZhangT_2020[,1]==0,4],
  df_Compare_ZhangT_2020R[df_Compare_ZhangT_2020R[,1]==-1,4],df_Compare_ZhangT_2020R[df_Compare_ZhangT_2020R[,1]!=-1 | df_Compare_ZhangT_2020R[,1]==0,4],
  df_Compare_SchmidtT_2019B[df_Compare_SchmidtT_2019B[,1]==-1,4],df_Compare_SchmidtT_2019B[df_Compare_SchmidtT_2019B[,1]!=-1 | df_Compare_SchmidtT_2019B[,1]==0,4],
  df_Compare_ChenX_2025[df_Compare_ChenX_2025[,1]==-1,4],df_Compare_ChenX_2025[df_Compare_ChenX_2025[,1]!=-1 | df_Compare_ChenX_2025[,1]==0,4],
  df_Compare_TakayanagiK_2023[df_Compare_TakayanagiK_2023[,1]==-1,4],df_Compare_TakayanagiK_2023[df_Compare_TakayanagiK_2023[,1]!=-1 | df_Compare_TakayanagiK_2023[,1]==0,4],
  df_Compare_WuZ_2025[df_Compare_WuZ_2025[,1]==-1,4],df_Compare_WuZ_2025[df_Compare_WuZ_2025[,1]!=-1 | df_Compare_WuZ_2025[,1]==0,4],
  df_Compare_WuZ_2025_nonmet[df_Compare_WuZ_2025_nonmet[,1]==-1,4],df_Compare_WuZ_2025_nonmet[df_Compare_WuZ_2025_nonmet[,1]!=-1 | df_Compare_WuZ_2025_nonmet[,1]==0,4],
  df_Compare_WuZ_2025_met_nonmet[df_Compare_WuZ_2025_met_nonmet[,1]==-1,4],df_Compare_WuZ_2025_met_nonmet[df_Compare_WuZ_2025_met_nonmet[,1]!=-1 | df_Compare_WuZ_2025_met_nonmet[,1]==0,4],
  df_Compare_ZhangL_2023[df_Compare_ZhangL_2023[,1]==-1,4],df_Compare_ZhangL_2023[df_Compare_ZhangL_2023[,1]!=-1 | df_Compare_ZhangL_2023[,1]==0,4],
  df_Compare_ChenJ_2021_OSCC[df_Compare_ChenJ_2021_OSCC[,1]==-1,4],df_Compare_ChenJ_2021_OSCC[df_Compare_ChenJ_2021_OSCC[,1]!=-1 | df_Compare_ChenJ_2021_OSCC[,1]==0,4],
  df_Compare_ChenJ_2021_OVH[df_Compare_ChenJ_2021_OVH[,1]==-1,4],df_Compare_ChenJ_2021_OVH[df_Compare_ChenJ_2021_OVH[,1]!=-1 | df_Compare_ChenJ_2021_OVH[,1]==0,4],
  df_Compare_CirsteaM_2022[df_Compare_CirsteaM_2022[,1]==-1,4],df_Compare_CirsteaM_2022[df_Compare_CirsteaM_2022[,1]!=-1 | df_Compare_CirsteaM_2022[,1]==0,4],
  df_Compare_FinkelsteinS_2025[df_Compare_FinkelsteinS_2025[,1]==-1,4],df_Compare_FinkelsteinS_2025[df_Compare_FinkelsteinS_2025[,1]!=-1 | df_Compare_FinkelsteinS_2025[,1]==0,4],
  df_Compare_IglesiasA_2024[df_Compare_IglesiasA_2024[,1]==-1,4],df_Compare_IglesiasA_2024[df_Compare_IglesiasA_2024[,1]!=-1 | df_Compare_IglesiasA_2024[,1]==0,4],
  df_Compare_NearingJ_2023_BC[df_Compare_NearingJ_2023_BC[,1]==-1,4],df_Compare_NearingJ_2023_BC[df_Compare_NearingJ_2023_BC[,1]!=-1 | df_Compare_NearingJ_2023_BC[,1]==0,4],
  df_Compare_NearingJ_2023_CC[df_Compare_NearingJ_2023_CC[,1]==-1,4],df_Compare_NearingJ_2023_CC[df_Compare_NearingJ_2023_CC[,1]!=-1 | df_Compare_NearingJ_2023_CC[,1]==0,4],
  df_Compare_NearingJ_2023_PC[df_Compare_NearingJ_2023_PC[,1]==-1,4],df_Compare_NearingJ_2023_PC[df_Compare_NearingJ_2023_PC[,1]!=-1 | df_Compare_NearingJ_2023_PC[,1]==0,4],
  df_Compare_RelvasM_2021_Dg1[df_Compare_RelvasM_2021_Dg1[,1]==-1,4],df_Compare_RelvasM_2021_Dg1[df_Compare_RelvasM_2021_Dg1[,1]!=-1 | df_Compare_RelvasM_2021_Dg1[,1]==0,4],
  df_Compare_RelvasM_2021_Dg2[df_Compare_RelvasM_2021_Dg2[,1]==-1,4],df_Compare_RelvasM_2021_Dg2[df_Compare_RelvasM_2021_Dg2[,1]!=-1 | df_Compare_RelvasM_2021_Dg2[,1]==0,4],
  df_Compare_RelvasM_2021_Dg3[df_Compare_RelvasM_2021_Dg3[,1]==-1,4],df_Compare_RelvasM_2021_Dg3[df_Compare_RelvasM_2021_Dg3[,1]!=-1 | df_Compare_RelvasM_2021_Dg3[,1]==0,4],
  df_Compare_RelvasM_2021_Pg1[df_Compare_RelvasM_2021_Pg1[,1]==-1,4],df_Compare_RelvasM_2021_Pg1[df_Compare_RelvasM_2021_Pg1[,1]!=-1 | df_Compare_RelvasM_2021_Pg1[,1]==0,4],
  df_Compare_RelvasM_2021_Pg2[df_Compare_RelvasM_2021_Pg2[,1]==-1,4],df_Compare_RelvasM_2021_Pg2[df_Compare_RelvasM_2021_Pg2[,1]!=-1 | df_Compare_RelvasM_2021_Pg2[,1]==0,4],
  df_Compare_JiY_2020[df_Compare_JiY_2020[,1]==-1,4],df_Compare_JiY_2020[df_Compare_JiY_2020[,1]!=-1 | df_Compare_JiY_2020[,1]==0,4],
  df_Compare_JiY_2020_Er[df_Compare_JiY_2020_Er[,1]==-1,4],df_Compare_JiY_2020_Er[df_Compare_JiY_2020_Er[,1]!=-1 | df_Compare_JiY_2020_Er[,1]==0,4],
  df_Compare_LiuY_2021[df_Compare_LiuY_2021[,1]==-1,4],df_Compare_LiuY_2021[df_Compare_LiuY_2021[,1]!=-1 | df_Compare_LiuY_2021[,1]==0,4],
  names = c(
    "ChenC_2018","ChenC_2018_postDis","FanX_2018a","FanX_2018a_Al",
    "ZhangT_2020","ZhangT_2020R","SchmidtT_2019B","ChenX_2025",
    "TakayanagiK_2023","WuZ_2025","WuZ_2025_nonmet","WuZ_2025_met_nonmet",
    "ZhangL_2023","ChenJ_2021_OSCC","ChenJ_2021_OVH","CirsteaM_2022",
    "FinkelsteinS_2025","IglesiasA_2024","NearingJ_2023_BC","NearingJ_2023_CC",
    "NearingJ_2023_PC","RelvasM_2021_Dg1","RelvasM_2021_Dg2","RelvasM_2021_Dg3",
    "RelvasM_2021_Pg1","RelvasM_2021_Pg2","JiY_2020","JiY_2020_Er","LiuY_2021"),
  side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"), las = 2)
dev.off()


## Now calculate the pvalue for each of the combination above and store in a df.
get_pvalue <- function(df) {
  g1 <- df[df[,1] == -1, 4]
  g2 <- df[df[,1] != -1 & df[,1] == 0, 4]
  
  # remove NA
  g1 <- g1[!is.na(g1)]
  g2 <- g2[!is.na(g2)]
  
  if(length(g1) >= 1 & length(g2) >= 1) {
    p <- wilcox.test(g1, g2)$p.value
  } else {
    p <- NA
  }
  
  return(p)
}

df_list <- list(
  ChenC_2018 = df_Compare_ChenC_2018,
  ChenC_2018_postDis = df_Compare_ChenC_2018_postDis,
  FanX_2018a = df_Compare_FanX_2018a,
  FanX_2018a_Al = df_Compare_FanX_2018a_Al,
  ZhangT_2020 = df_Compare_ZhangT_2020,
  ZhangT_2020R = df_Compare_ZhangT_2020R,
  SchmidtT_2019B = df_Compare_SchmidtT_2019B,
  ChenX_2025 = df_Compare_ChenX_2025,
  TakayanagiK_2023 = df_Compare_TakayanagiK_2023,
  WuZ_2025 = df_Compare_WuZ_2025,
  WuZ_2025_nonmet = df_Compare_WuZ_2025_nonmet,
  WuZ_2025_met_nonmet = df_Compare_WuZ_2025_met_nonmet,
  ZhangL_2023 = df_Compare_ZhangL_2023,
  ChenJ_2021_OSCC = df_Compare_ChenJ_2021_OSCC,
  ChenJ_2021_OVH = df_Compare_ChenJ_2021_OVH,
  CirsteaM_2022 = df_Compare_CirsteaM_2022,
  FinkelsteinS_2025 = df_Compare_FinkelsteinS_2025,
  IglesiasA_2024 = df_Compare_IglesiasA_2024,
  NearingJ_2023_BC = df_Compare_NearingJ_2023_BC,
  NearingJ_2023_CC = df_Compare_NearingJ_2023_CC,
  NearingJ_2023_PC = df_Compare_NearingJ_2023_PC,
  RelvasM_2021_Dg1 = df_Compare_RelvasM_2021_Dg1,
  RelvasM_2021_Dg2 = df_Compare_RelvasM_2021_Dg2,
  RelvasM_2021_Dg3 = df_Compare_RelvasM_2021_Dg3,
  RelvasM_2021_Pg1 = df_Compare_RelvasM_2021_Pg1,
  RelvasM_2021_Pg2 = df_Compare_RelvasM_2021_Pg2,
  JiY_2020 = df_Compare_JiY_2020,
  JiY_2020_Er = df_Compare_JiY_2020_Er,
  LiuY_2021 = df_Compare_LiuY_2021
)

pval_table <- data.frame(
  Study = names(df_list),
  P_value = sapply(df_list, get_pvalue)
)

pval_table$stars <- ifelse(
  is.na(pval_table$P_value),
  NA,
  as.character(cut(pval_table$P_value,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns")))
)




pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HACK_AllStudies.pdf",width = 20, height = 7)
beanplot(
  df_Compare_ChenC_2018[df_Compare_ChenC_2018[,1]==-1,5],df_Compare_ChenC_2018[df_Compare_ChenC_2018[,1]!=-1 | df_Compare_ChenC_2018[,1]==0,5],
  df_Compare_ChenC_2018_postDis[df_Compare_ChenC_2018_postDis[,1]==-1,5],df_Compare_ChenC_2018_postDis[df_Compare_ChenC_2018_postDis[,1]!=-1 | df_Compare_ChenC_2018_postDis[,1]==0,5],
  df_Compare_FanX_2018a[df_Compare_FanX_2018a[,1]==-1,5],df_Compare_FanX_2018a[df_Compare_FanX_2018a[,1]!=-1 | df_Compare_FanX_2018a[,1]==0,5],
  df_Compare_FanX_2018a_Al[df_Compare_FanX_2018a_Al[,1]==-1,5],df_Compare_FanX_2018a_Al[df_Compare_FanX_2018a_Al[,1]!=-1 | df_Compare_FanX_2018a_Al[,1]==0,5],
  df_Compare_ZhangT_2020[df_Compare_ZhangT_2020[,1]==-1,5],df_Compare_ZhangT_2020[df_Compare_ZhangT_2020[,1]!=-1 | df_Compare_ZhangT_2020[,1]==0,5],
  df_Compare_ZhangT_2020R[df_Compare_ZhangT_2020R[,1]==-1,5],df_Compare_ZhangT_2020R[df_Compare_ZhangT_2020R[,1]!=-1 | df_Compare_ZhangT_2020R[,1]==0,5],
  df_Compare_SchmidtT_2019B[df_Compare_SchmidtT_2019B[,1]==-1,5],df_Compare_SchmidtT_2019B[df_Compare_SchmidtT_2019B[,1]!=-1 | df_Compare_SchmidtT_2019B[,1]==0,5],
  df_Compare_ChenX_2025[df_Compare_ChenX_2025[,1]==-1,5],df_Compare_ChenX_2025[df_Compare_ChenX_2025[,1]!=-1 | df_Compare_ChenX_2025[,1]==0,5],
  df_Compare_TakayanagiK_2023[df_Compare_TakayanagiK_2023[,1]==-1,5],df_Compare_TakayanagiK_2023[df_Compare_TakayanagiK_2023[,1]!=-1 | df_Compare_TakayanagiK_2023[,1]==0,5],
  df_Compare_WuZ_2025[df_Compare_WuZ_2025[,1]==-1,5],df_Compare_WuZ_2025[df_Compare_WuZ_2025[,1]!=-1 | df_Compare_WuZ_2025[,1]==0,5],
  df_Compare_WuZ_2025_nonmet[df_Compare_WuZ_2025_nonmet[,1]==-1,5],df_Compare_WuZ_2025_nonmet[df_Compare_WuZ_2025_nonmet[,1]!=-1 | df_Compare_WuZ_2025_nonmet[,1]==0,5],
  df_Compare_WuZ_2025_met_nonmet[df_Compare_WuZ_2025_met_nonmet[,1]==-1,5],df_Compare_WuZ_2025_met_nonmet[df_Compare_WuZ_2025_met_nonmet[,1]!=-1 | df_Compare_WuZ_2025_met_nonmet[,1]==0,5],
  df_Compare_ZhangL_2023[df_Compare_ZhangL_2023[,1]==-1,5],df_Compare_ZhangL_2023[df_Compare_ZhangL_2023[,1]!=-1 | df_Compare_ZhangL_2023[,1]==0,5],
  df_Compare_ChenJ_2021_OSCC[df_Compare_ChenJ_2021_OSCC[,1]==-1,5],df_Compare_ChenJ_2021_OSCC[df_Compare_ChenJ_2021_OSCC[,1]!=-1 | df_Compare_ChenJ_2021_OSCC[,1]==0,5],
  df_Compare_ChenJ_2021_OVH[df_Compare_ChenJ_2021_OVH[,1]==-1,5],df_Compare_ChenJ_2021_OVH[df_Compare_ChenJ_2021_OVH[,1]!=-1 | df_Compare_ChenJ_2021_OVH[,1]==0,5],
  df_Compare_CirsteaM_2022[df_Compare_CirsteaM_2022[,1]==-1,5],df_Compare_CirsteaM_2022[df_Compare_CirsteaM_2022[,1]!=-1 | df_Compare_CirsteaM_2022[,1]==0,5],
  df_Compare_FinkelsteinS_2025[df_Compare_FinkelsteinS_2025[,1]==-1,5],df_Compare_FinkelsteinS_2025[df_Compare_FinkelsteinS_2025[,1]!=-1 | df_Compare_FinkelsteinS_2025[,1]==0,5],
  df_Compare_IglesiasA_2024[df_Compare_IglesiasA_2024[,1]==-1,5],df_Compare_IglesiasA_2024[df_Compare_IglesiasA_2024[,1]!=-1 | df_Compare_IglesiasA_2024[,1]==0,5],
  df_Compare_NearingJ_2023_BC[df_Compare_NearingJ_2023_BC[,1]==-1,5],df_Compare_NearingJ_2023_BC[df_Compare_NearingJ_2023_BC[,1]!=-1 | df_Compare_NearingJ_2023_BC[,1]==0,5],
  df_Compare_NearingJ_2023_CC[df_Compare_NearingJ_2023_CC[,1]==-1,5],df_Compare_NearingJ_2023_CC[df_Compare_NearingJ_2023_CC[,1]!=-1 | df_Compare_NearingJ_2023_CC[,1]==0,5],
  df_Compare_NearingJ_2023_PC[df_Compare_NearingJ_2023_PC[,1]==-1,5],df_Compare_NearingJ_2023_PC[df_Compare_NearingJ_2023_PC[,1]!=-1 | df_Compare_NearingJ_2023_PC[,1]==0,5],
  df_Compare_RelvasM_2021_Dg1[df_Compare_RelvasM_2021_Dg1[,1]==-1,5],df_Compare_RelvasM_2021_Dg1[df_Compare_RelvasM_2021_Dg1[,1]!=-1 | df_Compare_RelvasM_2021_Dg1[,1]==0,5],
  df_Compare_RelvasM_2021_Dg2[df_Compare_RelvasM_2021_Dg2[,1]==-1,5],df_Compare_RelvasM_2021_Dg2[df_Compare_RelvasM_2021_Dg2[,1]!=-1 | df_Compare_RelvasM_2021_Dg2[,1]==0,5],
  df_Compare_RelvasM_2021_Dg3[df_Compare_RelvasM_2021_Dg3[,1]==-1,5],df_Compare_RelvasM_2021_Dg3[df_Compare_RelvasM_2021_Dg3[,1]!=-1 | df_Compare_RelvasM_2021_Dg3[,1]==0,5],
  df_Compare_RelvasM_2021_Pg1[df_Compare_RelvasM_2021_Pg1[,1]==-1,5],df_Compare_RelvasM_2021_Pg1[df_Compare_RelvasM_2021_Pg1[,1]!=-1 | df_Compare_RelvasM_2021_Pg1[,1]==0,5],
  df_Compare_RelvasM_2021_Pg2[df_Compare_RelvasM_2021_Pg2[,1]==-1,5],df_Compare_RelvasM_2021_Pg2[df_Compare_RelvasM_2021_Pg2[,1]!=-1 | df_Compare_RelvasM_2021_Pg2[,1]==0,5],
  df_Compare_JiY_2020[df_Compare_JiY_2020[,1]==-1,5],df_Compare_JiY_2020[df_Compare_JiY_2020[,1]!=-1 | df_Compare_JiY_2020[,1]==0,5],
  df_Compare_JiY_2020_Er[df_Compare_JiY_2020_Er[,1]==-1,5],df_Compare_JiY_2020_Er[df_Compare_JiY_2020_Er[,1]!=-1 | df_Compare_JiY_2020_Er[,1]==0,5],
  df_Compare_LiuY_2021[df_Compare_LiuY_2021[,1]==-1,5],df_Compare_LiuY_2021[df_Compare_LiuY_2021[,1]!=-1 | df_Compare_LiuY_2021[,1]==0,5],
  names = c(
    "ChenC_2018","ChenC_2018_postDis","FanX_2018a","FanX_2018a_Al",
    "ZhangT_2020","ZhangT_2020R","SchmidtT_2019B","ChenX_2025",
    "TakayanagiK_2023","WuZ_2025","WuZ_2025_nonmet","WuZ_2025_met_nonmet",
    "ZhangL_2023","ChenJ_2021_OSCC","ChenJ_2021_OVH","CirsteaM_2022",
    "FinkelsteinS_2025","IglesiasA_2024","NearingJ_2023_BC","NearingJ_2023_CC",
    "NearingJ_2023_PC","RelvasM_2021_Dg1","RelvasM_2021_Dg2","RelvasM_2021_Dg3",
    "RelvasM_2021_Pg1","RelvasM_2021_Pg2","JiY_2020","JiY_2020_Er","LiuY_2021"),
  side="both",what=c(1,1,1,0),overallline="median",col=list("aquamarine","antiquewhite"), las = 2)
dev.off()

## 13, 16 both of this are showing opposite trend. 


get_pvalue2 <- function(df) {
  g1 <- df[df[,1] == -1, 5]
  g2 <- df[df[,1] != -1 & df[,1] == 0, 5]
  
  # remove NA
  g1 <- g1[!is.na(g1)]
  g2 <- g2[!is.na(g2)]
  
  if(length(g1) >= 1 & length(g2) >= 1) {
    p <- wilcox.test(g1, g2)$p.value
  } else {
    p <- NA
  }
  
  return(p)
}

pval_table2 <- data.frame(
  Study = names(df_list),
  P_value = sapply(df_list, get_pvalue2)
)

pval_table2$stars <- ifelse(
  is.na(pval_table2$P_value),
  NA,
  as.character(cut(pval_table2$P_value,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns")))
)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_Validation_Workspace.RData")




################## Number of species that are health associated and coming in how many studies. 
# create a list of all the dfs starting with df_Compare_
study_dfs <- c("df_Compare_ChenC_2018", "df_Compare_ChenC_2018_postDis", "df_Compare_FanX_2018a", "df_Compare_ZhangT_2020", "df_Compare_ZhangT_2020R", "df_Compare_SchmidtT_2019B", "df_Compare_ChenX_2025", "df_Compare_TakayanagiK_2023", "df_Compare_WuZ_2025", "df_Compare_WuZ_2025_nonmet", "df_Compare_WuZ_2025_met_nonmet", "df_Compare_ZhangL_2023",
               "df_Compare_ChenJ_2021_OSCC", "df_Compare_ChenJ_2021_OVH", "df_Compare_CirsteaM_2022", "df_Compare_FinkelsteinS_2025", "df_Compare_IglesiasA_2024", "df_Compare_NearingJ_2023_BC", "df_Compare_NearingJ_2023_CC", "df_Compare_NearingJ_2023_PC", "df_Compare_JiY_2020", "df_Compare_LiuY_2021", 
               "df_Compare_RelvasM_2021_Dg1", "df_Compare_RelvasM_2021_Dg2", "df_Compare_RelvasM_2021_Dg3", "df_Compare_RelvasM_2021_Pg1", "df_Compare_RelvasM_2021_Pg2" )


dir_mat <- data.frame(matrix(0,nrow = nrow(Combined_Saliva_Scores), ncol = length(study_dfs)))
rownames(dir_mat) <- rownames(Combined_Saliva_Scores)
colnames(dir_mat) <- study_dfs

## Now loop over each of the df from study_dfs and fill the values in dir_mat
for (s in study_dfs) {
  
  temp_df <- get(s)              # get df_Compare_*
  sp <- rownames(temp_df)        # species present in that study
  
  dir_mat[sp, s] <- temp_df[sp, "Direction"]
}


# add health score
dir_mat$HealthScore <- Combined_Saliva_Scores[rownames(dir_mat), 2]
dir_mat$HACK_Score <- Combined_Saliva_Scores[rownames(dir_mat), 5]

# count health associations
dir_mat$health_n <- apply(dir_mat[, study_dfs], 1, function(x) sum(x == -1))

dir_mat$range <- ifelse(dir_mat$health_n ==0, "<1",
                        ifelse(dir_mat$health_n >= 1 & dir_mat$health_n <=3, "[1-3]",
                        ifelse(dir_mat$health_n > 3 & dir_mat$health_n <= 6, "[4-6]", 
                        ifelse(dir_mat$health_n >= 7 & dir_mat$health_n <= 9, "[7-9]",
                        ifelse(dir_mat$health_n > 9, "[10+]", NA)))))

dir_mat$range <- factor(dir_mat$range, levels = c("<1", "[1-3]", "[4-6]", "[7-9]", "[10+]"))


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_Studies_Vs_HealthScore.pdf", height =5, width = 6)
boxplot(
  HealthScore ~ range,
  data = dir_mat,
  col = c("#f28e8e", "#ffd92f", "#a6d854","#6ca7e2","#59dd6b"),   # two different colors
  border = "black",
  ylab = "Health association score",
  xlab = "Number of studies with health association",
  main = "Health association vs reproducibility across studies",
  outline = FALSE)
dev.off()

library(dunn.test)
Dunntest_healthscore <- dunn.test(dir_mat$HealthScore, dir_mat$range, method = "bh")
Dunntest_healthscore
# $chi2
# [1] 74.30601

# $Z
#  [1] -2.6947452 -5.2673553  0.5401068 -4.4015094 -0.1631000 -1.0314341
#  [7]  2.5847535  3.2395207  6.6555816  5.2403636

# $P
#  [1] 3.522123e-03 6.920158e-08 2.945617e-01 5.375020e-06 4.352198e-01
#  [6] 1.511687e-01 4.872432e-03 5.986539e-04 1.410910e-11 8.013025e-08

# $P.adjusted
#  [1] 5.870205e-03 3.460079e-07 3.272908e-01 1.343755e-05 4.352198e-01
#  [6] 1.889608e-01 6.960617e-03 1.197308e-03 1.410910e-10 2.671008e-07

# $comparisons
#  [1] "[1-3] - [10+]" "[1-3] - [4-6]" "[10+] - [4-6]" "[1-3] - [7-9]"
#  [5] "[10+] - [7-9]" "[4-6] - [7-9]" "[1-3] - <1"    "[10+] - <1"



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_Studies_Vs_HACK_Score.pdf", height =5, width = 6)
boxplot(
  HACK_Score ~ range,
  data = dir_mat,
  col = c("#f28e8e", "#ffd92f", "#a6d854","#6ca7e2","#59dd6b"),   # two different colors
  border = "black",
  ylab = "HACK score",
  xlab = "Number of studies with health association",
  main = "Health association vs reproducibility across studies",
  outline = FALSE)
dev.off()

Dunntest_HACKscore <- dunn.test(dir_mat$HACK_Score, dir_mat$range, method = "bh")
Dunntest_HACKscore
# $chi2
# [1] 127.6555

# $Z
#  [1] -3.0846321 -5.6664258  0.7567935 -3.9567063  0.4162724 -0.4192211
#  [7]  5.9852892  4.3425575  8.8352582  5.8827566

# $P
#  [1] 1.019020e-03 7.290348e-09 2.245868e-01 3.799513e-05 3.386053e-01
#  [6] 3.375273e-01 1.080029e-09 7.041682e-06 4.993366e-19 2.017445e-09

# $P.adjusted
#  [1] 1.455743e-03 1.822587e-08 2.807335e-01 6.332522e-05 3.386053e-01
#  [6] 3.750303e-01 5.400144e-09 1.408336e-05 4.993366e-18 6.724817e-09

# $comparisons
#  [1] "[1-3] - [10+]" "[1-3] - [4-6]" "[10+] - [4-6]" "[1-3] - [7-9]"
#  [5] "[10+] - [7-9]" "[4-6] - [7-9]" "[1-3] - <1"    "[10+] - <1"


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_Validation_Workspace.RData")













######### See whether HACK score for health associated species is more or Health Score for health associated specie sis more (in comparison with disease associated species)

validation_study_list <- sub("df_Compare_","",ls(pattern="df_Compare"))

comparison_profile <- matrix(0,length(validation_study_list),2)
rownames(comparison_profile) <- validation_study_list
colnames(comparison_profile) <- c("Health","HACK")

for(i in 1:length(validation_study_list))
{
  study_name <- validation_study_list[i]
  temp_compare <- get(paste0("df_Compare_",study_name))
  dir <- sign(mean(temp_compare[temp_compare$Direction == -1,"Health"]) - mean(temp_compare[temp_compare$Direction != -1,"Health"])) 
  p_val <- wilcox.test(temp_compare[temp_compare$Direction == -1, "Health"],temp_compare[temp_compare$Direction != -1, "Health"])$p.value
  comparison_profile[study_name,"Health"] <- dir * ifelse(p_val <= 0.05,3,ifelse(p_val <= 0.1,2,1))
  
  dir <- sign(mean(temp_compare[temp_compare$Direction == -1,"HACKScore"]) - mean(temp_compare[temp_compare$Direction != -1,"HACKScore"])) 
  p_val <- wilcox.test(temp_compare[temp_compare$Direction == -1,"HACKScore"],temp_compare[temp_compare$Direction != -1, "HACKScore"])$p.value
  comparison_profile[study_name,"HACK"] <- dir * ifelse(p_val <= 0.05,3,ifelse(p_val <= 0.1,2,1))
  
}




association_consistency_validation <- as.data.frame(matrix(0,nrow(Combined_Saliva_Scores),length(validation_study_list)))
rownames(association_consistency_validation) <- rownames(Combined_Saliva_Scores)
colnames(association_consistency_validation) <- validation_study_list

for(i in 1:length(validation_study_list))
{
  study_name <- validation_study_list[i]
  temp_compare <- get(paste0("df_Compare_",study_name))
  association_consistency_validation[rownames(temp_compare),study_name] <- temp_compare$Direction
}

association_consistency_validation <- as.data.frame(apply(association_consistency_validation,2,function(x)(ifelse(is.na(x),0,x))))

association_consistency_validation$count <- apply(association_consistency_validation,1,function(x)(length(x[x==-1])))


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_HACKScore_BoxPlot.pdf", width = 8, height = 6)
boxplot(Combined_Saliva_Scores[,5] ~ cut(association_consistency_validation$count,
                                  breaks = c(0, 1, 11, 8)), col = c("#e09549", "#74abe6", "#93e688"), outline = TRUE, ylab = "HACK Score (Discovery Cohort)", xlab = "No. of studies with health association taxa")
dev.off()


groups <- cut(association_consistency_validation$count,
              breaks = c(0, 1, 11, 8))


library(dunn.test)
HACK_dunntest <- dunn.test(x = Combined_Saliva_Scores[,5],g = groups,method = "bh")
HACK_dunntest
# $chi2
# [1] 87.44089

# $Z
# [1] -8.329354 -5.882653 -2.996669

# $P
# [1] 4.064221e-17 2.018712e-09 1.364734e-03

# $P.adjusted
# [1] 1.219266e-16 3.028068e-09 1.364734e-03

# $comparisons
# [1] "(0,1] - (1,8]"  "(0,1] - (8,11]" "(1,8] - (8,11]"





pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_HealthScore_BoxPlot.pdf", width = 8, height = 6)
boxplot(Combined_Saliva_Scores[,2] ~ cut(association_consistency_validation$count,
                                  breaks = c(0, 1, 11, 8)), col = c("#e09549", "#74abe6", "#93e688"), outline = TRUE, ylab = "Health Score (Discovery Cohort)", xlab = "No. of studies with health association taxa")
dev.off()

Health_dunntest <- dunn.test(x = Combined_Saliva_Scores[,2],g = groups,method = "bh")
Health_dunntest
# $chi2
# [1] 43.12464

# $Z
# [1] -5.313912 -4.875079 -3.046620

# $P
# [1] 5.364824e-08 5.438254e-07 1.157151e-03

# $P.adjusted
# [1] 1.609447e-07 8.157381e-07 1.157151e-03

# $comparisons
# [1] "(0,1] - (1,8]"  "(0,1] - (8,11]" "(1,8] - (8,11]"




gc()
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_2HealthAssociation_Validation_Workspace.RData")


