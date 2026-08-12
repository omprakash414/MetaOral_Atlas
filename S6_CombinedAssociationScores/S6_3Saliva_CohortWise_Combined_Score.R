







load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_CohortWise_CoreAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_2Saliva_CohortWise_HealthAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_3Saliva_CohortWise_StabilityAssociationScore.RData")

######## Combine all three dfs into one df.
saliva_CohortWise_HealthAssociationScore$species <- rownames(saliva_CohortWise_HealthAssociationScore)
saliva_CohortWise_HealthAssociationScore <- saliva_CohortWise_HealthAssociationScore[rownames(saliva_CohortWise_CoreAssociationScore),]

Saliva_CohortWise_StabilityScore$species <- rownames(Saliva_CohortWise_StabilityScore)
Saliva_CohortWise_StabilityScore <- Saliva_CohortWise_StabilityScore[rownames(saliva_CohortWise_CoreAssociationScore),]

Saliva_combined_scores_cohortwise <- cbind(saliva_CohortWise_HealthAssociationScore, Saliva_CohortWise_StabilityScore, saliva_CohortWise_CoreAssociationScore)
Saliva_combined_scores_cohortwise <- Saliva_combined_scores_cohortwise[,c(1,2,3,4,5,7,8,9,11,12,13,14)]

######## Now count the combined HACK score for all the cohorts (16s/WGS. 16s_exposure has only health and core, as in stability we don't have any study with exposure)
# Caclulate Combined Score
library(LaplacesDemon)
library(DescTools)

rank_scale=function(x){
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}


Combined_saliva_16s_HACK_score <- Saliva_combined_scores_cohortwise[,c("HealthAssociationScore_16s","StabilityScore_16s","CoreAssociationScore_16s")]
Combined_saliva_WGS_HACK_score <- Saliva_combined_scores_cohortwise[,c("HealthAssociationScore_WGS","StabilityScore_WGS","CoreAssociationScore_WGS")]

# Take only rows that have all non-NA values first
Combined_saliva_16s_HACK_score <- Combined_saliva_16s_HACK_score[complete.cases(Combined_saliva_16s_HACK_score),]
Combined_saliva_WGS_HACK_score <- Combined_saliva_WGS_HACK_score[complete.cases(Combined_saliva_WGS_HACK_score),]

Combined_saliva_16s_HACK_score$HACKScoreUnscaled <- apply(Combined_saliva_16s_HACK_score[,1:3],1,mean) * (1-apply(Combined_saliva_16s_HACK_score[,1:3],1,Gini))
Combined_saliva_16s_HACK_score$HACScore <- rank_scale(Combined_saliva_16s_HACK_score$HACKScoreUnscaled)

Combined_saliva_WGS_HACK_score$HACKScoreUnscaled <- apply(Combined_saliva_WGS_HACK_score[,1:3],1,mean) * (1-apply(Combined_saliva_WGS_HACK_score[,1:3],1,Gini))
Combined_saliva_WGS_HACK_score$HACScore <- rank_scale(Combined_saliva_WGS_HACK_score$HACKScoreUnscaled)

# Order based on the value in decreasing order
Combined_saliva_16s_HACK_score <- Combined_saliva_16s_HACK_score[order(Combined_saliva_16s_HACK_score$HACScore, decreasing = TRUE),]
Combined_saliva_WGS_HACK_score <- Combined_saliva_WGS_HACK_score[order(Combined_saliva_WGS_HACK_score$HACScore, decreasing = TRUE),]



####### Now load the Combined HACK score counted earlier (Mixed of 16s and WGS) and then add the HACKScore column to above two dfs
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedScores.RData")

Combined_saliva_16s_HACK_score$combined_HACKScore <- Combined_Saliva_Scores$HAC_Score_RankScaled[match(rownames(Combined_saliva_16s_HACK_score), rownames(Combined_Saliva_Scores))]
Combined_saliva_WGS_HACK_score$combined_HACKScore <- Combined_Saliva_Scores$HAC_Score_RankScaled[match(rownames(Combined_saliva_WGS_HACK_score), rownames(Combined_Saliva_Scores))]


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_CohortWise_Combined_Score_Workspace.RData")



######## Now plot the correlation plots of 16s/WGS HACK score with combined HACK score and see the correlation value. 

library(ggplot2)
library(ggrepel)
library(dplyr)

### 16s HACK score vs combined HACK score
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_Correlation_16sHACK_HACKScore.pdf",width = 5,height = 6)
ggplot(
  Combined_saliva_16s_HACK_score %>%
    dplyr::filter(
      !is.na(combined_HACKScore),
      !is.na(HACScore)
    ),
  aes(
    x = HACScore,
    y = combined_HACKScore, color = HACScore >= 0.80 & combined_HACKScore >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "HACK score 16s",
    y = "HACK score Overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )
dev.off()

cor.test(
  Combined_saliva_16s_HACK_score$HACScore,
  Combined_saliva_16s_HACK_score$combined_HACKScore,
  use = "complete.obs",
  method = "pearson"
)

#         Pearson's product-moment correlation

# data:  Combined_saliva_16s_HACK_score$HACScore and Combined_saliva_16s_HACK_score$combined_HACKScore
# t = 71.857, df = 480, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9482216 0.9635265
# sample estimates:
#       cor 
# 0.9565278

format((cor.test(
  Combined_saliva_16s_HACK_score$HACScore,
  Combined_saliva_16s_HACK_score$combined_HACKScore,
  use = "complete.obs",
  method = "pearson"
))$p.value, scientific = TRUE, digits = 3)
#  "5.1e-259"




### WGS HACK score vs combined HACK score
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_Correlation_WGSHACK_HACKScore.pdf",width = 5,height = 6)
ggplot(
  Combined_saliva_WGS_HACK_score %>%
    dplyr::filter(
      !is.na(combined_HACKScore),
      !is.na(HACScore)
    ),
  aes(
    x = HACScore,
    y = combined_HACKScore, color = HACScore >= 0.80 & combined_HACKScore >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "HACK score WGS",
    y = "HACK score overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )
dev.off()

cor.test(
  Combined_saliva_WGS_HACK_score$HACScore,
  Combined_saliva_WGS_HACK_score$combined_HACKScore,
  use = "complete.obs",
  method = "pearson"
)

#         Pearson's product-moment correlation

# data:  Combined_saliva_WGS_HACK_score$HACScore and Combined_saliva_WGS_HACK_score$combined_HACKScore
# t = 9.7172, df = 227, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4435269 0.6275491
# sample estimates:
#       cor 
# 0.5420032

 format((cor.test(
  Combined_saliva_WGS_HACK_score$HACScore,
  Combined_saliva_WGS_HACK_score$combined_HACKScore,
  use = "complete.obs",
  method = "pearson"
))$p.value, scientific = TRUE, digits = 3)
# [1] "6.93e-19"






########### Now do the correlation og Core association score and Health Association Score for Exposure 16s cohort. 

Combined_saliva_Exposure_HAC_score <- Saliva_combined_scores_cohortwise[,c("HealthAssociationScore_16s_exposure","CoreAssociationScore_16s_exposure")]
Combined_saliva_Exposure_HAC_score <- Combined_saliva_Exposure_HAC_score[complete.cases(Combined_saliva_Exposure_HAC_score),]

Combined_saliva_Exposure_HAC_score$HACScoreUnscaled <- apply(Combined_saliva_Exposure_HAC_score[,1:2],1,mean) * (1-apply(Combined_saliva_Exposure_HAC_score[,1:2],1,Gini))
Combined_saliva_Exposure_HAC_score$HACScore <- rank_scale(Combined_saliva_Exposure_HAC_score$HACScoreUnscaled)

Combined_saliva_Exposure_HAC_score$combined_HACKScore <- Combined_Saliva_Scores$HAC_Score_RankScaled[match(rownames(Combined_saliva_Exposure_HAC_score), rownames(Combined_Saliva_Scores))]



########## do correlation of HACScore exposure and combined_HACKScore

### HAC score exposure (Core and Health only) vs combined HACK score
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_Correlation_ExposureHAC_HACKScore.pdf",width = 5,height = 6)
ggplot(
  Combined_saliva_Exposure_HAC_score %>%
    dplyr::filter(
      !is.na(combined_HACKScore),
      !is.na(HACScore),
      combined_HACKScore != 0,
      HACScore != 0
    ),
  aes(
    x = HACScore,
    y = combined_HACKScore, color = HACScore >= 0.80 & combined_HACKScore >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black", na.rm = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "HAC score Exposure",
    y = "HACK score overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )
dev.off()

temp_df <- Combined_saliva_Exposure_HAC_score[
  !is.na(Combined_saliva_Exposure_HAC_score$HACScore) &
  !is.na(Combined_saliva_Exposure_HAC_score$combined_HACKScore) &
  Combined_saliva_Exposure_HAC_score$HACScore != 0 &
  Combined_saliva_Exposure_HAC_score$combined_HACKScore != 0,
]

cor.test(
  temp_df$HACScore,
  temp_df$combined_HACKScore,
  method = "pearson"
)

#         Pearson's product-moment correlation

# data:  temp_df$HACScore and temp_df$combined_HACKScore
# t = 20.584, df = 381, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6744629 0.7698558
# sample estimates:
#       cor 
# 0.7256275

format((cor.test(
  Combined_saliva_Exposure_HAC_score$HACScore,
  Combined_saliva_Exposure_HAC_score$combined_HACKScore,
  use = "complete.obs",
  method = "pearson"
))$p.value, scientific = TRUE, digits = 3)
# [1] "1.61e-34"



################## Now we can also correlate the combined HACK Score with Core and health score with two lines in one plot itself
library(dplyr)
library(tidyr)

df_long <- Combined_saliva_Exposure_HAC_score %>%
  dplyr::filter(
    !is.na(combined_HACKScore),
    !is.na(HealthAssociationScore_16s_exposure),
    !is.na(CoreAssociationScore_16s_exposure)
  ) %>%
  dplyr::select(combined_HACKScore, HealthAssociationScore_16s_exposure, CoreAssociationScore_16s_exposure) %>%
  pivot_longer(
    cols = c(HealthAssociationScore_16s_exposure, CoreAssociationScore_16s_exposure),
    names_to = "ScoreType",
    values_to = "ScoreValue"
  )



library(ggplot2)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_Correlation_Health_Core_vs_combinedHACKScore.pdf",
    width = 8, height = 6)

ggplot(df_long,
       aes(x = ScoreValue,
           y = combined_HACKScore,
           color = ScoreType)) +
  geom_point(size = 3.7, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2) +
  theme_bw(base_size = 14) +
  labs(
    x = "Score value",
    y = "combined HACK score overall",
    color = "Score type"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )

dev.off()





save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_CohortWise_Combined_Score_Workspace.RData")


############ Create a single df where, a combined HACK score alogn with 16s HACK score, WGS HACK score, HAC score of exposure. Also add a column where species for which of them are non-zero
library(dplyr)
temp_all_HACKscores <- Combined_Saliva_Scores[, "HAC_Score_RankScaled", drop = FALSE]
temp_all_HACKscores$HACK_16s <- Combined_saliva_16s_HACK_score$HACScore[match(rownames(temp_all_HACKscores), rownames(Combined_saliva_16s_HACK_score))]
temp_all_HACKscores$HACK_WGS <- Combined_saliva_WGS_HACK_score$HACScore[match(rownames(temp_all_HACKscores), rownames(Combined_saliva_WGS_HACK_score))]

temp_all_HACKscores$HAC_Exposure <- Combined_saliva_Exposure_HAC_score$HACScore[match(rownames(temp_all_HACKscores), rownames(Combined_saliva_Exposure_HAC_score))]


write.csv(temp_all_HACKscores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_All_Scores_Supplementary.csv")


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_CohortWise_Combined_Score_Workspace.RData")

save(Combined_saliva_16s_HACK_score,Combined_saliva_Exposure_HAC_score,Combined_saliva_WGS_HACK_score, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_CohortWise_Combined_Score_Dfs.RData")




############## We have't done the correlation of 16s vs WGS and also plot it. 

df_16s_wgs <- temp_all_HACKscores %>%
  dplyr::filter(
    !is.na(HACK_16s),
    !is.na(HACK_WGS),
    HACK_16s != 0,
    HACK_WGS != 0
  )

cor.test(df_16s_wgs$HACK_16s, df_16s_wgs$HACK_WGS)
#         Pearson's product-moment correlation

# data:  df_16s_wgs$HACK_16s and df_16s_wgs$HACK_WGS
# t = 7.2258, df = 209, p-value = 9.175e-12
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3320672 0.5489994
# sample estimates:
#       cor 
# 0.4470832 


pdf(
  "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_Correlation_16sHACK_WGSHACK.pdf",
  width = 5,
  height = 6
)

ggplot(
  df_16s_wgs,
  aes(
    x = HACK_16s,
    y = HACK_WGS,
    color = HACK_16s >= 0.80 & HACK_WGS >= 0.80
  )
) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE" = "red"),
    guide = "none"
  ) +
  geom_point(size = 3.7, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black") +
  theme_bw(base_size = 14) +
  labs(
    x = "HACK score 16S",
    y = "HACK score WGS"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

dev.off()


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_CohortWise_Combined_Score_Workspace.RData")
