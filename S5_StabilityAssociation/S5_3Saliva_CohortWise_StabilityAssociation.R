############ This script calculate the Longitudinal distance of a sample 



########### Divide the data into 
# positive_exposure <- c("alcoholic", "smoker", "tobacco", "e-cigar")

# keep_row <- sapply(Longitudinal_Dist_Metadata$exposure, function(x) {
#   # split by comma
#   tokens <- trimws(unlist(strsplit(tolower(x), "[, ;]")))
  
#   # check if ANY positive exposure is present
#   any(tokens %in% positive_exposure)
# })

# Longitudinal_Dist_Metadata_16s_exposure <- Longitudinal_Dist_Metadata[keep_row, ]

####### only two studies - "MakinenA_2023" (16s)  "DanckertN_2024" (WGS) are of exposure, so, we haven't did stability analysis for this separately.
####### Here MakinenA_2023" and  "DanckertN_2024"  were considered in 16s and WGS cohort respectively and not in exposure



######## Load the distance matrix calculated earlier and other needed data to calculate the stability association score
S5_1Followup_Distance_calculation_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_1Followup_Distance_calculation_Workspace.RData", envir = S5_1Followup_Distance_calculation_Workspace)
attach(S5_1Followup_Distance_calculation_Workspace)
SpDf_saliva_Long <- SpDf_saliva_Long
saliva_AssociatedSpecies <- saliva_AssociatedSpecies
std_meta_cols <- std_meta_cols
MetadataDf_saliva_Long <- MetadataDf_saliva_Long
detach(S5_1Followup_Distance_calculation_Workspace)


source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")


#############
#############
############# Calculate the distance for 16S data first 
# Separate the 16s metadata and species profile
MetadataDf_saliva_Long_16s <- MetadataDf_saliva_Long[MetadataDf_saliva_Long$seq_type == "16s",]
SpDf_saliva_Long_16s <- SpDf_saliva_Long[rownames(MetadataDf_saliva_Long_16s),]
SpDf_saliva_Long_16s <- SpDf_saliva_Long_16s[,colSums(SpDf_saliva_Long_16s)>0]
species_16s <- colnames(SpDf_saliva_Long_16s)

#### 
Longitudinal_Dist_Metadata_16s <- Longitudinal_Microbiome_Distance(MetadataDf_saliva_Long_16s, SpDf_saliva_Long_16s)

#### add the distance calculated alomg with standard metadata columns to species profile
# confirm the rownames in Longitudinal_Dist_Metadata_16s and species profile are alighed in same order
all(rownames(Longitudinal_Dist_Metadata_16s) == rownames(SpDf_saliva_Long_16s))

# Add the standard columns to species profile
SpDf_saliva_Long_withDist_16s <- bind_cols(SpDf_saliva_Long_16s, Longitudinal_Dist_Metadata_16s[,std_meta_cols])


############# Calculate the Stability score for 16s cohort
#### Calculate the stability scores using both distance metrics
## Aitchison
aitch_iters_16s <- Meta_lm_Iterrative(data = SpDf_saliva_Long_withDist_16s,species = species_16s,dist_var = "Aitchison_dist",study_var = "Study_Name", n_studies = 4)
aitch_res_16s <- Compute_Meta_Stability(iter_list = aitch_iters_16s,species = species_16s)

## Bray–Curtis
bray_iters_16s <- Meta_lm_Iterrative(data = SpDf_saliva_Long_withDist_16s,species = species_16s,dist_var = "BrayCurtis_dist",study_var = "Study_Name", n_studies = 4)
bray_res_16s <- Compute_Meta_Stability(iter_list = bray_iters_16s,species = species_16s)

## Combine the stability scores from both distance metrics
stabilityRank_16s <- cbind(aitch_res_16s$stability,bray_res_16s$stability)
colnames(stabilityRank_16s) <- c("mean_Aitchison","IQR_Aitchison","mean_BrayCurtis","IQR_BrayCurtis")


stabilityRank_16s$MeanStabilityScore <- rowMeans(stabilityRank_16s[, c(1,3)])
stabilityRank_16s <- stabilityRank_16s[order(-stabilityRank_16s$MeanStabilityScore), ]

Saliva_StabilityScore_16s <- data.frame(species = rownames(stabilityRank_16s),StabilityScore_16s = stabilityRank_16s$MeanStabilityScore)
rownames(Saliva_StabilityScore_16s) <- Saliva_StabilityScore_16s$species







#############
#############
############# Calculate the distance for WGS data first 
# Separate the 16s metadata and species profile
MetadataDf_saliva_Long_WGS <- MetadataDf_saliva_Long[MetadataDf_saliva_Long$seq_type == "WGS",]
SpDf_saliva_Long_WGS <- SpDf_saliva_Long[rownames(MetadataDf_saliva_Long_WGS),]
SpDf_saliva_Long_WGS <- SpDf_saliva_Long_WGS[,colSums(SpDf_saliva_Long_WGS)>0]
species_WGS <- colnames(SpDf_saliva_Long_WGS)


#### 
Longitudinal_Dist_Metadata_WGS <- Longitudinal_Microbiome_Distance(MetadataDf_saliva_Long_WGS, SpDf_saliva_Long_WGS)

#### add the distance calculated alomg with standard metadata columns to species profile
# confirm the rownames in Longitudinal_Dist_Metadata_WGS and species profile are alighed in same order
all(rownames(Longitudinal_Dist_Metadata_WGS) == rownames(SpDf_saliva_Long_WGS))

# Add the standard columns to species profile
SpDf_saliva_Long_withDist_WGS <- bind_cols(SpDf_saliva_Long_WGS, Longitudinal_Dist_Metadata_WGS[,std_meta_cols])



########## Calculate the Stability score for WGS cohort
######## Calculate the stability scores using both distance metrics
## Aitchison
aitch_iters_WGS <- Meta_lm_Iterrative(data = SpDf_saliva_Long_withDist_WGS,species = species_WGS,dist_var = "Aitchison_dist",study_var = "Study_Name", n_studies = 4)
aitch_res_WGS <- Compute_Meta_Stability(iter_list = aitch_iters_WGS,species = species_WGS)

## Bray–Curtis
bray_iters_WGS <- Meta_lm_Iterrative(data = SpDf_saliva_Long_withDist_WGS,species = species_WGS,dist_var = "BrayCurtis_dist",study_var = "Study_Name", n_studies = 4)
bray_res_WGS <- Compute_Meta_Stability(iter_list = bray_iters_WGS,species = species_WGS)

## Combine the stability scores from both distance metrics
stabilityRank_WGS <- cbind(aitch_res_WGS$stability,bray_res_WGS$stability)
colnames(stabilityRank_WGS) <- c("mean_Aitchison","IQR_Aitchison","mean_BrayCurtis","IQR_BrayCurtis")

stabilityRank_WGS$MeanStabilityScore <- rowMeans(stabilityRank_WGS[, c(1,3)])
stabilityRank_WGS <- stabilityRank_WGS[order(-stabilityRank_WGS$MeanStabilityScore), ]

Saliva_StabilityScore_WGS <- data.frame(species = rownames(stabilityRank_WGS),StabilityScore_WGS = stabilityRank_WGS$MeanStabilityScore)
rownames(Saliva_StabilityScore_WGS) <- Saliva_StabilityScore_WGS$species












########### Combine with overall stability score and then compare
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2Saliva_StabilityAssociationScore.RData")

Saliva_CohortWise_StabilityScore <-
  Saliva_StabilityScore %>%
  left_join(Saliva_StabilityScore_16s, by = "species") %>%
  left_join(Saliva_StabilityScore_WGS, by = "species")

rownames(Saliva_CohortWise_StabilityScore) <- Saliva_CohortWise_StabilityScore$species


save(Saliva_CohortWise_StabilityScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_3Saliva_CohortWise_StabilityAssociationScore.RData")
write.csv(Saliva_CohortWise_StabilityScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_3Saliva_CohortWise_StabilityAssociationScore.csv", row.names = T, quote = F)
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_3Saliva_CohortWise_StabilityAssociation_Workspace.RData")




#################
##################
#################

# Compare all the scores
####### 16s vs WGS
library(ggplot2)
library(ggrepel)
library(dplyr)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_3Saliva_Correlation_16s_WGS_AllSpecies.pdf",width = 5, height = 6)
ggplot(
  Saliva_CohortWise_StabilityScore %>%
    dplyr::filter(
      !is.na(StabilityScore_WGS),
      !is.na(StabilityScore_16s)
    ),
  aes(x = StabilityScore_WGS, y = StabilityScore_16s, color = StabilityScore_WGS >= 0.80 & StabilityScore_16s >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "SA score WGS",
    y = "SA score 16S"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

dev.off()

cor.test(Saliva_CohortWise_StabilityScore$StabilityScore_WGS,Saliva_CohortWise_StabilityScore$StabilityScore_16s,use = "complete.obs",method = "pearson")
#         Pearson's product-moment correlation

# data:  Saliva_CohortWise_StabilityScore$StabilityScore_WGS and Saliva_CohortWise_StabilityScore$StabilityScore_16s
# t = 4.5216, df = 220, p-value = 1.003e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1663123 0.4076218
# sample estimates:
#       cor 
# 0.2915997


#### Overall and 16s

library(ggplot2)
library(dplyr)

### Correlation between Stability score (16S vs Overall)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_3Saliva_Correlation_16s_Overall_AllSpecies.pdf",width = 5, height = 6)
ggplot(
  Saliva_CohortWise_StabilityScore %>%
    dplyr::filter(
      !is.na(StabilityScore_16s),
      !is.na(StabilityScore)
    ),
  aes(x = StabilityScore_16s, y = StabilityScore, color = StabilityScore_16s >= 0.80 & StabilityScore >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "SA score 16S",
    y = "SA score Overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

dev.off()

cor.test(Saliva_CohortWise_StabilityScore$StabilityScore_16s,Saliva_CohortWise_StabilityScore$StabilityScore,use = "complete.obs",method = "pearson")
#         Pearson's product-moment correlation

# data:  Saliva_CohortWise_StabilityScore$StabilityScore_16s and Saliva_CohortWise_StabilityScore$StabilityScore
# t = 54.746, df = 480, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9149702 0.9398000
# sample estimates:
#       cor 
# 0.9284146

format((cor.test(Saliva_CohortWise_StabilityScore$StabilityScore_16s,Saliva_CohortWise_StabilityScore$StabilityScore,use = "complete.obs",method = "pearson"))$p.value, scientific = TRUE, digits = 3)
# [1] "1.58e-208"

#### Overall vs WGS
### Correlation between Stability score (WGS vs Overall)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_3Saliva_Correlation_WGS_Overall_AllSpecies.pdf",width = 5, height = 6)
ggplot(
  Saliva_CohortWise_StabilityScore %>%
    dplyr::filter(
      !is.na(StabilityScore_WGS),
      !is.na(StabilityScore)
    ),
  aes(x = StabilityScore_WGS, y = StabilityScore, color = StabilityScore_WGS >= 0.80 & StabilityScore >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "SA score WGS",
    y = "SA score Overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

dev.off()


cor.test(Saliva_CohortWise_StabilityScore$StabilityScore_WGS,Saliva_CohortWise_StabilityScore$StabilityScore,use = "complete.obs",method = "pearson")
#         Pearson's product-moment correlation

# data:  Saliva_CohortWise_StabilityScore$StabilityScore_WGS and Saliva_CohortWise_StabilityScore$StabilityScore
# t = 13.72, df = 237, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5880825 0.7305451
# sample estimates:
#       cor 
# 0.6653276 

format((cor.test(Saliva_CohortWise_StabilityScore$StabilityScore_WGS,Saliva_CohortWise_StabilityScore$StabilityScore,use = "complete.obs",method = "pearson"))$p.value, scientific = TRUE, digits = 3)
# [1] "6.37e-32"

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_3Saliva_CohortWise_StabilityAssociation_Workspace.RData")
