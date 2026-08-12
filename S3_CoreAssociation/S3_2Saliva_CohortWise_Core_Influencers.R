### See the Saliva Core Influencers in different cohorts (16s, WGS, exposure Cohort)


######### Import the Core Influencers results for Saliva
S3_1Saliva_Core_Detection_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_Core_Detection_Workspace.RData", envir = S3_1Saliva_Core_Detection_Workspace)
attach(S3_1Saliva_Core_Detection_Workspace)
saliva_Core_Influencers <- saliva_Core_Influencers
saliva_prevalDf <- saliva_prevalDf
saliva_r2_rankedDf <- saliva_r2_rankedDf
saliva_pval_Df <- saliva_pval_Df
saliva_CoreAssociationScore <- saliva_CoreAssociationScore
MetadataDf_saliva_control <- MetadataDf_saliva_control
detach(S3_1Saliva_Core_Detection_Workspace)
rm(S3_1Saliva_Core_Detection_Workspace)

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")
######### Extract the r2 df and Prevalance df for 16s studies only.
study_16s_list <- unique(MetadataDf_saliva_control[MetadataDf_saliva_control$seq_type == "16s",]$study_name)
study_16s_exposure_list <- study_16s_list[grepl("exposure",study_16s_list)]  ##
study_16s_only_list <- setdiff(study_16s_list,study_16s_exposure_list)


saliva_r2_rankedDf_16s <- saliva_r2_rankedDf[,study_16s_only_list]
saliva_prevalDf_16s <- saliva_prevalDf[,study_16s_only_list] 
saliva_prevalDf_16s <- saliva_prevalDf_16s[rowSums(saliva_prevalDf_16s)>0,]
saliva_r2_rankedDf_16s <- saliva_r2_rankedDf_16s[rownames(saliva_prevalDf_16s),]


saliva_r2_rankedDf_16s_exposure <- saliva_r2_rankedDf[,study_16s_exposure_list]
saliva_prevalDf_16s_exposure <- saliva_prevalDf[,study_16s_exposure_list]
saliva_prevalDf_16s_exposure <- saliva_prevalDf_16s_exposure[rowSums(saliva_prevalDf_16s_exposure)>0,]
saliva_r2_rankedDf_16s_exposure <- saliva_r2_rankedDf_16s_exposure[rownames(saliva_prevalDf_16s_exposure),]


study_WGS_list <- unique(MetadataDf_saliva_control[MetadataDf_saliva_control$seq_type == "WGS",]$study_name)
saliva_r2_rankedDf_WGS <- saliva_r2_rankedDf[,study_WGS_list]
saliva_prevalDf_WGS <- saliva_prevalDf[,study_WGS_list]
saliva_prevalDf_WGS <- saliva_prevalDf_WGS[rowSums(saliva_prevalDf_WGS)>0,]
saliva_r2_rankedDf_WGS <- saliva_r2_rankedDf_WGS[rownames(saliva_prevalDf_WGS),]


######### Now calculate the Core Score for each of the cohort
saliva_CoreKeyStoneDf_16s <- data.frame(apply(saliva_r2_rankedDf_16s,2,function(x)(ifelse(x>=0.75,1,0))) * apply(saliva_prevalDf_16s,2,function(x)(ifelse(x>=0.85,1,0))))
saliva_CoreKeyStoneDf_16s_exposure <- data.frame(apply(saliva_r2_rankedDf_16s_exposure,2,function(x)(ifelse(x>=0.75,1,0))) * apply(saliva_prevalDf_16s_exposure,2,function(x)(ifelse(x>=0.85,1,0))))
saliva_CoreKeyStoneDf_WGS <- data.frame(apply(saliva_r2_rankedDf_WGS,2,function(x)(ifelse(x>=0.75,1,0))) * apply(saliva_prevalDf_WGS,2,function(x)(ifelse(x>=0.85,1,0))))


############ Get the saliva CoreKeyStones and their study wise distribution
saliva_Core_Representation_16s <- summarize_core_keystone_detection(saliva_CoreKeyStoneDf_16s)
saliva_16s_CoreKeyStoneDf2 <- saliva_Core_Representation_16s$core_keystone_by_species

saliva_Core_Representation_16s_exposure <- summarize_core_keystone_detection(saliva_CoreKeyStoneDf_16s_exposure)
saliva_16s_exposure_CoreKeyStoneDf2 <- saliva_Core_Representation_16s_exposure$core_keystone_by_species

saliva_Core_Representation_WGS <- summarize_core_keystone_detection(saliva_CoreKeyStoneDf_WGS)
saliva_WGS_CoreKeyStoneDf2 <- saliva_Core_Representation_WGS$core_keystone_by_species


############ Add the CoreAssociation Score for each of the above dfs
saliva_16s_CoreKeyStoneDf2$CoreAssociationScore <- rank_scale((saliva_16s_CoreKeyStoneDf2$studies_detected)/length(study_16s_only_list))
saliva_16s_exposure_CoreKeyStoneDf2$CoreAssociationScore <- rank_scale((saliva_16s_exposure_CoreKeyStoneDf2$studies_detected)/length(study_16s_exposure_list))
saliva_WGS_CoreKeyStoneDf2$CoreAssociationScore <- rank_scale((saliva_WGS_CoreKeyStoneDf2$studies_detected)/length(study_WGS_list))

########### Extract the Core Association Scores
saliva_CoreAssociationScore_16s <- saliva_16s_CoreKeyStoneDf2[,tail(colnames(saliva_16s_CoreKeyStoneDf2), 2)]
saliva_CoreAssociationScore_16s <- saliva_CoreAssociationScore_16s[order(saliva_CoreAssociationScore_16s[, 2], decreasing = TRUE), ]
colnames(saliva_CoreAssociationScore_16s) <- c("species","CoreAssociationScore_16s")

saliva_CoreAssociationScore_16s_exposure <- saliva_16s_exposure_CoreKeyStoneDf2[,tail(colnames(saliva_16s_exposure_CoreKeyStoneDf2), 2)]
saliva_CoreAssociationScore_16s_exposure <- saliva_CoreAssociationScore_16s_exposure[order(saliva_CoreAssociationScore_16s_exposure[, 2], decreasing = TRUE), ]
colnames(saliva_CoreAssociationScore_16s_exposure) <- c("species","CoreAssociationScore_16s_exposure")

saliva_CoreAssociationScore_WGS <- saliva_WGS_CoreKeyStoneDf2[,tail(colnames(saliva_WGS_CoreKeyStoneDf2), 2)]
saliva_CoreAssociationScore_WGS <- saliva_CoreAssociationScore_WGS[order(saliva_CoreAssociationScore_WGS[, 2], decreasing = TRUE), ]
colnames(saliva_CoreAssociationScore_WGS) <- c("species","CoreAssociationScore_WGS")

############ Combine all the three cohort's association score along with association score of previously combined cohorts - saliva_CoreAssociationScore   

library(dplyr)

saliva_CohortWise_CoreAssociationScore <-
  saliva_CoreAssociationScore %>%
  left_join(saliva_CoreAssociationScore_16s, by = "species") %>%
  left_join(saliva_CoreAssociationScore_16s_exposure, by = "species") %>%
  left_join(saliva_CoreAssociationScore_WGS, by = "species")

rownames(saliva_CohortWise_CoreAssociationScore) <- saliva_CohortWise_CoreAssociationScore$species

save(saliva_CohortWise_CoreAssociationScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_CohortWise_CoreAssociationScore.RData")
write.csv(saliva_CohortWise_CoreAssociationScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_CohortWise_CoreAssociationScore.csv", row.names = T, quote = F)


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_CohortWise_Core_Influencers_Workspace.RData")





####### see the correlation of these CoreAssociationScores
# top17_species <- rownames(saliva_CohortWise_CoreAssociationScore[order(saliva_CohortWise_CoreAssociationScore$CoreAssociationScore,decreasing = TRUE),])[1:17]
# top67_species <- rownames(saliva_CohortWise_CoreAssociationScore[order(saliva_CohortWise_CoreAssociationScore$CoreAssociationScore,decreasing = TRUE),])[1:67]

# saliva_CohortWise_CoreAssociationScore$HAC_Top17 <- ifelse(saliva_CohortWise_CoreAssociationScore$species %in% top17_species,"1","0")

### Correlation 16s v/s WGS
library(ggplot2)
library(dplyr)
library(ggrepel)

### corr between  wgs and 16s
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_Correlation_16s_WGS.pdf",width = 5,height = 6)
ggplot(
  saliva_CohortWise_CoreAssociationScore %>%
    dplyr::filter(
      !is.na(CoreAssociationScore_WGS),
      !is.na(CoreAssociationScore_16s)),
  aes(x = CoreAssociationScore_WGS,y = CoreAssociationScore_16s, color = CoreAssociationScore_WGS >= 0.80 & CoreAssociationScore_16s >= 0.80)) +
  geom_point(size = 3.7, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black") +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  theme_bw(base_size = 14) +
  labs(x = "Core score WGS",y = "Core score 16S") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"))
dev.off()

cor.test(saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_WGS,saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s,use = "complete.obs",method = "pearson")
#         Pearson's product-moment correlation

# data:  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_WGS and saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s
# t = 12.091, df = 264, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5136004 0.6692155
# sample estimates:
#       cor 
# 0.5969947
format((cor.test(saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_WGS,saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s, use = "complete.obs"))$p.value, scientific = TRUE, digits = 20)
# 4.4559402671489488841e-27


### corr batween Overall Core-Association Score and 16s Core-Association Score
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_Correlation_16s_Overall.pdf",width = 5,height = 6)

ggplot(
  saliva_CohortWise_CoreAssociationScore %>%
    dplyr::filter(
      !is.na(CoreAssociationScore_16s),
      !is.na(CoreAssociationScore)
    ),
  aes(x = CoreAssociationScore_16s,y = CoreAssociationScore, color = CoreAssociationScore >= 0.80 & CoreAssociationScore_16s >= 0.80)) +
  geom_point(size = 3.7, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black") +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  theme_bw(base_size = 14) +
  labs(
    x = "CS  16S",
    y = "CS overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )
dev.off()
cor.test(
  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s,
  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore,
  use = "complete.obs",
  method = "pearson"
)
#         Pearson's product-moment correlation

# data:  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s and saliva_CohortWise_CoreAssociationScore$CoreAssociationScore
# t = 122.57, df = 483, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9812636 0.9868530
# sample estimates:
#       cor 
# 0.9843032

format((cor.test(saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s,saliva_CohortWise_CoreAssociationScore$CoreAssociationScore, use = "complete.obs"))$p.value, scientific = TRUE, digits = 20)
# 0e+00


### corr batween Overall Core-Association Score and WGS Core-Association Score
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_Correlation_WGS_Overall.pdf",width = 5,height = 6)

ggplot(
  saliva_CohortWise_CoreAssociationScore %>%
    dplyr::filter(
      !is.na(CoreAssociationScore_WGS),
      !is.na(CoreAssociationScore)
    ),
  aes(x = CoreAssociationScore_WGS,y = CoreAssociationScore, color = CoreAssociationScore_WGS >= 0.80 & CoreAssociationScore >= 0.80)) +
  geom_point(size = 3.7, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black") +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  theme_bw(base_size = 14) +
  labs(
    x = "CS WGS",
    y = "CS overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )
dev.off()
cor.test(
  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_WGS,
  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore,
  use = "complete.obs",
  method = "pearson"
)

#         Pearson's product-moment correlation

# data:  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_WGS and saliva_CohortWise_CoreAssociationScore$CoreAssociationScore
# t = 13.814, df = 278, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5628621 0.7026615
# sample estimates:
#       cor 
# 0.6379893

format((cor.test(saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_WGS,saliva_CohortWise_CoreAssociationScore$CoreAssociationScore, use = "complete.obs"))$p.value, scientific = TRUE, digits = 3)
# [1] "2.11e-33"


### corr between Overall Core-Association Score and Exposure
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_Correlation_Exposure_Overall.pdf",width = 5,height = 6)
ggplot(
  saliva_CohortWise_CoreAssociationScore %>%
    dplyr::filter(
      !is.na(CoreAssociationScore_16s_exposure),
      !is.na(CoreAssociationScore)
    ),
  aes(x = CoreAssociationScore_16s_exposure,y = CoreAssociationScore, color = CoreAssociationScore_16s_exposure >= 0.80 & CoreAssociationScore >= 0.80)) +
  geom_point(size = 3.7, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black") +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  theme_bw(base_size = 14) +
  labs(
    x = "CS Exposure",
    y = "CS overall"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )
dev.off()

cor.test(
  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s_exposure,
  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore,
  use = "complete.obs",
  method = "pearson"
)

#         Pearson's product-moment correlation

# data:  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s_exposure and saliva_CohortWise_CoreAssociationScore$CoreAssociationScore
# t = 34.023, df = 483, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8116621 0.8643839
# sample estimates:
#       cor 
# 0.8399945 

format((cor.test(saliva_CohortWise_CoreAssociationScore$CoreAssociationScore,saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s_exposure, use = "complete.obs"))$p.value, scientific = TRUE, digits = 3)
# [1] "2.44e-130"





### corr between 16s Core-Association Score (i.e Non-Exposure) and Exposure
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_Correlation_Exposure16s_Nonexposure16s.pdf",width = 5,height = 6)
ggplot(
  saliva_CohortWise_CoreAssociationScore %>%
    dplyr::filter(
      !is.na(CoreAssociationScore_16s_exposure),
      !is.na(CoreAssociationScore_16s)
    ),
  aes(x = CoreAssociationScore_16s_exposure,y = CoreAssociationScore_16s, color = CoreAssociationScore_16s_exposure >= 0.80 & CoreAssociationScore_16s >= 0.80)) +
  scale_color_manual(
    values = c("FALSE" = "gray40", "TRUE"  = "red"), guide = "none") +
  geom_point(size = 3.7, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 2, color = "black") +
  theme_bw(base_size = 14) +
  labs(
    x = "CS 16s Exposure",
    y = "CS 16s"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )
dev.off()

cor.test(
  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s_exposure,
  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s,
  use = "complete.obs",
  method = "pearson"
)

        Pearson's product-moment correlation

data:  saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s_exposure and saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s
t = 31.517, df = 483, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.7888367 0.8474099
sample estimates:
      cor 
0.8202622 

format((cor.test(saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s,saliva_CohortWise_CoreAssociationScore$CoreAssociationScore_16s_exposure, use = "complete.obs"))$p.value, scientific = TRUE, digits = 3)
# [1] "2.9e-119"













save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_2Saliva_CohortWise_Core_Influencers_Workspace.RData")
