
library(dplyr)
library(vegan)
library(compositions)
library(psych)
library(ggplot2)
library(dunn.test)

## Functions to use
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




## Import the validation data and then separate the data for stability analsyis and that too study wise
S9_1CoreAssociation_Validation_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_Validation_Workspace.RData", envir = S9_1CoreAssociation_Validation_Workspace)
attach(S9_1CoreAssociation_Validation_Workspace)
MetadataDf_saliva_validation <- MetadataDf_saliva_validation
SpDf_saliva_Validation <- SpDf_saliva_Validation
study_list <- study_list
Combined_Saliva_Scores <- Combined_Saliva_Scores
detach(S9_1CoreAssociation_Validation_Workspace)
rm(S9_1CoreAssociation_Validation_Workspace)

SpDf_saliva_Validation$study_name <- NULL



#################################
################################# ZhangT_2020
# Filter ZhangT_2020
ZhangT_2020_metadata <- MetadataDf_saliva_validation %>% filter(study_name == "ZhangT_2020")

# Extract T1 diseased samples
ZhangT_2020_t1 <- ZhangT_2020_metadata %>% filter(timepoint == "t1", study_condition == "Diseased") %>% select(sample_id, subject_id)

# Extract T2 samples
ZhangT_2020_t2 <- ZhangT_2020_metadata %>% filter(timepoint == "t2") %>% select(sample_id, subject_id) %>% rename(followup_sample = sample_id)

# Join T1 with T2 based on subject_id
ZhangT_2020_t1_followup <- ZhangT_2020_t1 %>% left_join(ZhangT_2020_t2, by = "subject_id")

# Now bind T2 rows separately and assign NA
ZhangT_2020_t2_withNA <- ZhangT_2020_t2 %>% mutate(sample_id = followup_sample,followup_sample = NA) %>% select(sample_id, subject_id)

# Final dataframe
ZhangT_2020_metadata_followup <- bind_rows(ZhangT_2020_t1_followup,ZhangT_2020_t2_withNA)
rownames(ZhangT_2020_metadata_followup) <- ZhangT_2020_metadata_followup$sample_id
ZhangT_2020_metadata_followup$subject_id <- NULL
colnames(ZhangT_2020_metadata_followup) <- c("T0","T1")

ZhangT_2020_SpDf <- SpDf_saliva_Validation[rownames(ZhangT_2020_metadata_followup),]
ZhangT_2020_SpDf <- ZhangT_2020_SpDf[colSums(ZhangT_2020_SpDf)>0]
ZhangT_2020_SpDf <- ZhangT_2020_SpDf[rowSums(ZhangT_2020_SpDf)>0,]

### Calculate the Bray and Aitchison distance
ZhangT_2020_bray_follow_up <- bray_followup(ZhangT_2020_SpDf,ZhangT_2020_metadata_followup)
ZhangT_2020_aitchison_follow_up <- aitchison_followup(ZhangT_2020_SpDf,ZhangT_2020_metadata_followup)

## add the followup sample ids from T1 column of ZhangT_2020_metadata_followup to new column in ZhangT_2020_SpDf using maching rownames of ZhangT_2020_SpDf and TO column in ZhangT_2020_metadata_followup
ZhangT_2020_SpDf_modified <- ZhangT_2020_SpDf
ZhangT_2020_SpDf_modified$followup_sample_id <- ZhangT_2020_metadata_followup$T1[match(rownames(ZhangT_2020_SpDf_modified),ZhangT_2020_metadata_followup$T0)]

ZhangT_2020_SpDf_modified$study_name <- "ZhangT_2020"
ZhangT_2020_SpDf_modified$timepoint <- ZhangT_2020_metadata$timepoint[match(rownames(ZhangT_2020_SpDf_modified),ZhangT_2020_metadata$sample_id)]
ZhangT_2020_SpDf_modified$study_condition <- ZhangT_2020_metadata$study_condition[match(rownames(ZhangT_2020_SpDf_modified),ZhangT_2020_metadata$sample_id)]
ZhangT_2020_SpDf_modified$bray_dist <- ZhangT_2020_bray_follow_up[match(rownames(ZhangT_2020_SpDf_modified),rownames(ZhangT_2020_bray_follow_up)),1]
ZhangT_2020_SpDf_modified$aitchison_dist <- ZhangT_2020_aitchison_follow_up[match(rownames(ZhangT_2020_SpDf_modified),rownames(ZhangT_2020_aitchison_follow_up)),1]


### Calculate the correlation
corr_ZhangT_2020_bray <- corr.test(ZhangT_2020_SpDf,ZhangT_2020_bray_follow_up[rownames(ZhangT_2020_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")
corr_ZhangT_2020_aitchison <- corr.test(ZhangT_2020_SpDf,ZhangT_2020_aitchison_follow_up[rownames(ZhangT_2020_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")

df_corr_ZhangT_2020_bray <- data.frame("R"=corr_ZhangT_2020_bray$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_bray$r)),],"P"=corr_ZhangT_2020_bray$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_bray$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_bray$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_bray$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_bray$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_bray$r)),5])
df_corr_ZhangT_2020_aitchison <- data.frame("R"=corr_ZhangT_2020_aitchison$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_aitchison$r)),],"P"=corr_ZhangT_2020_aitchison$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_aitchison$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_aitchison$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_aitchison$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_aitchison$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ZhangT_2020_aitchison$r)),5])

### add the directionnality and then select the species with negative and positive R
df_corr_ZhangT_2020_bray$dir <- ifelse(df_corr_ZhangT_2020_bray[,2]<= 0.05,2*sign(df_corr_ZhangT_2020_bray[,1]),sign(df_corr_ZhangT_2020_bray[,1]))
df_corr_ZhangT_2020_aitchison$dir <- ifelse(df_corr_ZhangT_2020_aitchison[,2]<= 0.05,2*sign(df_corr_ZhangT_2020_aitchison[,1]),sign(df_corr_ZhangT_2020_aitchison[,1]))

sig_negative_ZhangT_2020 <- union(rownames(df_corr_ZhangT_2020_aitchison[!is.na(df_corr_ZhangT_2020_aitchison[,7])&(df_corr_ZhangT_2020_aitchison[,7]==-2),]),rownames(df_corr_ZhangT_2020_bray[!is.na(df_corr_ZhangT_2020_bray[,7])&(df_corr_ZhangT_2020_bray[,7]==-2),]))
not_negative_ZhangT_2020 <- setdiff(union(rownames(df_corr_ZhangT_2020_aitchison),rownames(df_corr_ZhangT_2020_bray)),sig_negative_ZhangT_2020)


library(ggplot2)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ZhangT_2020_bray_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ZhangT_2020_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_ZhangT_2020_bray$R, df_corr_ZhangT_2020_bray$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_ZhangT_2020_bray$R and df_corr_ZhangT_2020_bray$Stability
# t = -8.7279, df = 387, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.4854002 -0.3189753
# sample estimates:
#        cor 
# -0.4055431


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ZhangT_2020_bray_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ZhangT_2020_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ZhangT_2020_bray_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ZhangT_2020_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ZhangT_2020_bray_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ZhangT_2020_bray,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ZhangT_2020_aitchison_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ZhangT_2020_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_ZhangT_2020_aitchison$R, df_corr_ZhangT_2020_aitchison$Stability)
#         Pearson's product-moment correlation

# data:  df_corr_ZhangT_2020_aitchison$R and df_corr_ZhangT_2020_aitchison$Stability
# t = -10.085, df = 387, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5315223 -0.3737246
# sample estimates:
#        cor 
# -0.4562023

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ZhangT_2020_aitchison_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ZhangT_2020_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ZhangT_2020_aitchison_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ZhangT_2020_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ZhangT_2020_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ZhangT_2020_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()





#################################
################################# SchmidtT_2019B
# Filter SchmidtT_2019B
SchmidtT_2019B_metadata <- MetadataDf_saliva_validation %>% filter(study_name == "SchmidtT_2019B")
SchmidtT_2019B_metadata <- SchmidtT_2019B_metadata[SchmidtT_2019B_metadata$disease != "T2D",]

library(dplyr)

SchmidtT_2019B_metadata_followup <-
  SchmidtT_2019B_metadata %>%
  filter(timepoint %in% c("t1", "t2", "t3")) %>%
  select(sample_id, subject_id, timepoint) %>%
  
  # Convert timepoint to ordered factor
  mutate(timepoint = factor(timepoint, levels = c("t1","t2","t3"))) %>%
  
  arrange(subject_id, timepoint) %>%
  
  group_by(subject_id) %>%
  
  # Assign next sample within subject
  mutate(followup_sample = lead(sample_id)) %>%
  
  ungroup()

SchmidtT_2019B_metadata_followup <- data.frame(SchmidtT_2019B_metadata_followup)
rownames(SchmidtT_2019B_metadata_followup) <-SchmidtT_2019B_metadata_followup$sample_id
SchmidtT_2019B_metadata_followup$subject_id <- NULL
SchmidtT_2019B_metadata_followup$timepoint <- NULL
colnames(SchmidtT_2019B_metadata_followup) <- c("T0","T1")

SchmidtT_2019B_SpDf <- SpDf_saliva_Validation[rownames(SchmidtT_2019B_metadata),]
SchmidtT_2019B_SpDf <- SchmidtT_2019B_SpDf[colSums(SchmidtT_2019B_SpDf) >0]
SchmidtT_2019B_SpDf <- SchmidtT_2019B_SpDf[rowSums(SchmidtT_2019B_SpDf)>0,]


### Calculate the Bray and Aitchison distance
SchmidtT_2019B_bray_follow_up <- bray_followup(SchmidtT_2019B_SpDf,SchmidtT_2019B_metadata_followup)
SchmidtT_2019B_aitchison_follow_up <- aitchison_followup(SchmidtT_2019B_SpDf,SchmidtT_2019B_metadata_followup)


## add the followup sample ids from T1 column of SchmidtT_2019B_metadata_followup to new column in SchmidtT_2019B_SpDf using maching rownames of SchmidtT_2019B_SpDf and TO column in SchmidtT_2019B_metadata_followup
SchmidtT_2019B_SpDf_modified <- SchmidtT_2019B_SpDf
SchmidtT_2019B_SpDf_modified$followup_sample_id <- SchmidtT_2019B_metadata_followup$T1[match(rownames(SchmidtT_2019B_SpDf_modified),SchmidtT_2019B_metadata_followup$T0)]

SchmidtT_2019B_SpDf_modified$study_name <- "SchmidtT_2019B"
SchmidtT_2019B_SpDf_modified$timepoint <- SchmidtT_2019B_metadata$timepoint[match(rownames(SchmidtT_2019B_SpDf_modified),SchmidtT_2019B_metadata$sample_id)]
SchmidtT_2019B_SpDf_modified$study_condition <- SchmidtT_2019B_metadata$study_condition[match(rownames(SchmidtT_2019B_SpDf_modified),SchmidtT_2019B_metadata$sample_id)]
SchmidtT_2019B_SpDf_modified$bray_dist <- SchmidtT_2019B_bray_follow_up[match(rownames(SchmidtT_2019B_SpDf_modified),rownames(SchmidtT_2019B_bray_follow_up)),1]
SchmidtT_2019B_SpDf_modified$aitchison_dist <- SchmidtT_2019B_aitchison_follow_up[match(rownames(SchmidtT_2019B_SpDf_modified),rownames(SchmidtT_2019B_aitchison_follow_up)),1]



### Calculate the correlation
corr_SchmidtT_2019B_bray <- corr.test(SchmidtT_2019B_SpDf,SchmidtT_2019B_bray_follow_up[rownames(SchmidtT_2019B_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")
corr_SchmidtT_2019B_aitchison <- corr.test(SchmidtT_2019B_SpDf,SchmidtT_2019B_aitchison_follow_up[rownames(SchmidtT_2019B_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")

df_corr_SchmidtT_2019B_bray <- data.frame("R"=corr_SchmidtT_2019B_bray$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_bray$r)),],"P"=corr_SchmidtT_2019B_bray$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_bray$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_bray$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_bray$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_bray$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_bray$r)),5])
df_corr_SchmidtT_2019B_aitchison <- data.frame("R"=corr_SchmidtT_2019B_aitchison$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_aitchison$r)),],"P"=corr_SchmidtT_2019B_aitchison$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_aitchison$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_aitchison$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_aitchison$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_aitchison$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_SchmidtT_2019B_aitchison$r)),5])

### add the directionnality and then select the species with negative and positive R
df_corr_SchmidtT_2019B_bray$dir <- ifelse(df_corr_SchmidtT_2019B_bray[,2]<= 0.05,2*sign(df_corr_SchmidtT_2019B_bray[,1]),sign(df_corr_SchmidtT_2019B_bray[,1]))
df_corr_SchmidtT_2019B_aitchison$dir <- ifelse(df_corr_SchmidtT_2019B_aitchison[,2]<= 0.05,2*sign(df_corr_SchmidtT_2019B_aitchison[,1]),sign(df_corr_SchmidtT_2019B_aitchison[,1]))

sig_negative_SchmidtT_2019B <- union(rownames(df_corr_SchmidtT_2019B_aitchison[!is.na(df_corr_SchmidtT_2019B_aitchison[,7])&(df_corr_SchmidtT_2019B_aitchison[,7]==-2),]),rownames(df_corr_SchmidtT_2019B_bray[!is.na(df_corr_SchmidtT_2019B_bray[,7])&(df_corr_SchmidtT_2019B_bray[,7]==-2),]))
not_negative_SchmidtT_2019B <- setdiff(union(rownames(df_corr_SchmidtT_2019B_aitchison),rownames(df_corr_SchmidtT_2019B_bray)),sig_negative_SchmidtT_2019B)


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_SchmidtT_2019B_bray_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_SchmidtT_2019B_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_SchmidtT_2019B_bray$R, df_corr_SchmidtT_2019B_bray$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_SchmidtT_2019B_bray$R and df_corr_SchmidtT_2019B_bray$Stability
# t = -7.4998, df = 224, p-value = 1.482e-12
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5465474 -0.3372151
# sample estimates:
#        cor 
# -0.4480005

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_SchmidtT_2019B_bray_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_SchmidtT_2019B_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_SchmidtT_2019B_bray_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_SchmidtT_2019B_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_SchmidtT_2019B_bray_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_SchmidtT_2019B_bray,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_SchmidtT_2019B_aitchison_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_SchmidtT_2019B_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_SchmidtT_2019B_aitchison$R, df_corr_SchmidtT_2019B_aitchison$Stability, method = "pearson") 
#         Pearson's product-moment correlation

# data:  df_corr_SchmidtT_2019B_aitchison$R and df_corr_SchmidtT_2019B_aitchison$Stability
# t = -7.1356, df = 224, p-value = 1.323e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5310364 -0.3177026
# sample estimates:
#        cor 
# -0.4303601

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_SchmidtT_2019B_aitchison_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_SchmidtT_2019B_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_SchmidtT_2019B_aitchison_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_SchmidtT_2019B_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_SchmidtT_2019B_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_SchmidtT_2019B_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()


#################################
################################# ChenC_2018
# Filter ChenC_2018
ChenC_2018_metadata <- MetadataDf_saliva_validation %>% filter(study_name == "ChenC_2018")

# Extract T1 diseased samples
ChenC_2018_t1 <- ChenC_2018_metadata %>% filter(timepoint == "Pre", study_condition == "Diseased") %>% select(sample_id, subject_id)

# Extract T2 samples
ChenC_2018_t2 <- ChenC_2018_metadata %>% filter(timepoint == "Post") %>% select(sample_id, subject_id) %>% rename(followup_sample = sample_id)

# Join T1 with T2 based on subject_id
ChenC_2018_t1_followup <- ChenC_2018_t1 %>% left_join(ChenC_2018_t2, by = "subject_id")

# Now bind T2 rows separately and assign NA
ChenC_2018_t2_withNA <- ChenC_2018_t2 %>% mutate(sample_id = followup_sample,followup_sample = NA) %>% select(sample_id, subject_id)

# Final dataframe
ChenC_2018_metadata_followup <- bind_rows(ChenC_2018_t1_followup,ChenC_2018_t2_withNA)
rownames(ChenC_2018_metadata_followup) <- ChenC_2018_metadata_followup$sample_id
ChenC_2018_metadata_followup$subject_id <- NULL
colnames(ChenC_2018_metadata_followup) <- c("T0","T1")


ChenC_2018_SpDf <- SpDf_saliva_Validation[rownames(ChenC_2018_metadata_followup),]
ChenC_2018_SpDf <- ChenC_2018_SpDf[colSums(ChenC_2018_SpDf)>0]
ChenC_2018_SpDf <- ChenC_2018_SpDf[rowSums(ChenC_2018_SpDf)>0,]

### Calculate the Bray and Aitchison distance
ChenC_2018_bray_follow_up <- bray_followup(ChenC_2018_SpDf,ChenC_2018_metadata_followup)
ChenC_2018_aitchison_follow_up <- aitchison_followup(ChenC_2018_SpDf,ChenC_2018_metadata_followup)



## add the followup sample ids from T1 column of ChenC_2018_metadata_followup to new column in ChenC_2018_SpDf using maching rownames of ChenC_2018_SpDf and TO column in ChenC_2018_metadata_followup
ChenC_2018_SpDf_modified <- ChenC_2018_SpDf
ChenC_2018_SpDf_modified$followup_sample_id <- ChenC_2018_metadata_followup$T1[match(rownames(ChenC_2018_SpDf_modified),ChenC_2018_metadata_followup$T0)]

ChenC_2018_SpDf_modified$study_name <- "ChenC_2018"
ChenC_2018_SpDf_modified$timepoint <- ChenC_2018_metadata$timepoint[match(rownames(ChenC_2018_SpDf_modified),ChenC_2018_metadata$sample_id)]
ChenC_2018_SpDf_modified$study_condition <- ChenC_2018_metadata$study_condition[match(rownames(ChenC_2018_SpDf_modified),ChenC_2018_metadata$sample_id)]
ChenC_2018_SpDf_modified$bray_dist <- ChenC_2018_bray_follow_up[match(rownames(ChenC_2018_SpDf_modified),rownames(ChenC_2018_bray_follow_up)),1]
ChenC_2018_SpDf_modified$aitchison_dist <- ChenC_2018_aitchison_follow_up[match(rownames(ChenC_2018_SpDf_modified),rownames(ChenC_2018_aitchison_follow_up)),1]




### Calculate the correlation
corr_ChenC_2018_bray <- corr.test(ChenC_2018_SpDf,ChenC_2018_bray_follow_up[rownames(ChenC_2018_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")
corr_ChenC_2018_aitchison <- corr.test(ChenC_2018_SpDf,ChenC_2018_aitchison_follow_up[rownames(ChenC_2018_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")

df_corr_ChenC_2018_bray <- data.frame("R"=corr_ChenC_2018_bray$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_bray$r)),],"P"=corr_ChenC_2018_bray$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_bray$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_bray$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_bray$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_bray$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_bray$r)),5])
df_corr_ChenC_2018_aitchison <- data.frame("R"=corr_ChenC_2018_aitchison$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_aitchison$r)),],"P"=corr_ChenC_2018_aitchison$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_aitchison$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_aitchison$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_aitchison$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_aitchison$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_ChenC_2018_aitchison$r)),5])

### add the directionnality and then select the species with negative and positive R
df_corr_ChenC_2018_bray$dir <- ifelse(df_corr_ChenC_2018_bray[,2]<= 0.05,2*sign(df_corr_ChenC_2018_bray[,1]),sign(df_corr_ChenC_2018_bray[,1]))
df_corr_ChenC_2018_aitchison$dir <- ifelse(df_corr_ChenC_2018_aitchison[,2]<= 0.05,2*sign(df_corr_ChenC_2018_aitchison[,1]),sign(df_corr_ChenC_2018_aitchison[,1]))

sig_negative_ChenC_2018 <- union(rownames(df_corr_ChenC_2018_aitchison[!is.na(df_corr_ChenC_2018_aitchison[,7])&(df_corr_ChenC_2018_aitchison[,7]==-2),]),rownames(df_corr_ChenC_2018_bray[!is.na(df_corr_ChenC_2018_bray[,7])&(df_corr_ChenC_2018_bray[,7]==-2),]))
not_negative_ChenC_2018 <- setdiff(union(rownames(df_corr_ChenC_2018_aitchison),rownames(df_corr_ChenC_2018_bray)),sig_negative_ChenC_2018)


library(ggplot2)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ChenC_2018_bray_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ChenC_2018_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_ChenC_2018_bray$R, df_corr_ChenC_2018_bray$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_ChenC_2018_bray$R and df_corr_ChenC_2018_bray$Stability
# t = -3.5139, df = 242, p-value = 0.0005269
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.33660094 -0.09743927
# sample estimates:
#        cor 
# -0.2203289 

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ChenC_2018_bray_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ChenC_2018_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ChenC_2018_bray_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ChenC_2018_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ChenC_2018_bray_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ChenC_2018_bray,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ChenC_2018_aitchison_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ChenC_2018_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()
cor.test(df_corr_ChenC_2018_aitchison$R, df_corr_ChenC_2018_aitchison$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_ChenC_2018_aitchison$R and df_corr_ChenC_2018_aitchison$Stability
# t = -1.3073, df = 242, p-value = 0.1924
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.20714785  0.04229005
# sample estimates:
#         cor 
# -0.08374052

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ChenC_2018_aitchison_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ChenC_2018_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ChenC_2018_aitchison_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ChenC_2018_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_ChenC_2018_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_ChenC_2018_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()



#################################
################################# LincardonT_2023
# Filter LincardonT_2023
LincardonT_2023_metadata <- MetadataDf_saliva_validation %>% filter(study_name == "LicandroH_2023")

# unique(LincardonT_2023_metadata$timepoint)
#  [1] "8"  "7"  "6"  "2"  "4"  "10" "5"  "9"  "3"  "1"  "11" "12"

# Join T1, T2, T2 .. T12 based on subject_id
LincardonT_2023_metadata_followup <-
  LincardonT_2023_metadata %>%
  filter(timepoint %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")) %>%
  select(sample_id, subject_id, timepoint) %>%
  
  # Convert timepoint to ordered factor
  mutate(timepoint = factor(timepoint, levels = c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
  
  arrange(subject_id, timepoint) %>%
  
  group_by(subject_id) %>%
  
  # Assign next sample within subject
  mutate(followup_sample = lead(sample_id)) %>%
  
  ungroup()


LincardonT_2023_metadata_followup <- data.frame(LincardonT_2023_metadata_followup)
rownames(LincardonT_2023_metadata_followup) <-LincardonT_2023_metadata_followup$sample_id
LincardonT_2023_metadata_followup$subject_id <- NULL
LincardonT_2023_metadata_followup$timepoint <- NULL
colnames(LincardonT_2023_metadata_followup) <- c("T0","T1")

LincardonT_2023_SpDf <- SpDf_saliva_Validation[rownames(LincardonT_2023_metadata),]
LincardonT_2023_SpDf <- LincardonT_2023_SpDf[colSums(LincardonT_2023_SpDf) >0]
LincardonT_2023_SpDf <- LincardonT_2023_SpDf[rowSums(LincardonT_2023_SpDf)>0,]


### Calculate the Bray and Aitchison distance
LincardonT_2023_bray_follow_up <- bray_followup(LincardonT_2023_SpDf,LincardonT_2023_metadata_followup)
LincardonT_2023_aitchison_follow_up <- aitchison_followup(LincardonT_2023_SpDf,LincardonT_2023_metadata_followup)


## add the followup sample ids from T1 column of LincardonT_2023_metadata_followup to new column in LincardonT_2023_SpDf using maching rownames of LincardonT_2023_SpDf and TO column in LincardonT_2023_metadata_followup
LincardonT_2023_SpDf_modified <- LincardonT_2023_SpDf
LincardonT_2023_SpDf_modified$followup_sample_id <- LincardonT_2023_metadata_followup$T1[match(rownames(LincardonT_2023_SpDf_modified),LincardonT_2023_metadata_followup$T0)]

LincardonT_2023_SpDf_modified$study_name <- "LincardonT_2023"
LincardonT_2023_SpDf_modified$timepoint <- LincardonT_2023_metadata$timepoint[match(rownames(LincardonT_2023_SpDf_modified),LincardonT_2023_metadata$sample_id)]
LincardonT_2023_SpDf_modified$study_condition <- LincardonT_2023_metadata$study_condition[match(rownames(LincardonT_2023_SpDf_modified),LincardonT_2023_metadata$sample_id)]
LincardonT_2023_SpDf_modified$bray_dist <- LincardonT_2023_bray_follow_up[match(rownames(LincardonT_2023_SpDf_modified),rownames(LincardonT_2023_bray_follow_up)),1]
LincardonT_2023_SpDf_modified$aitchison_dist <- LincardonT_2023_aitchison_follow_up[match(rownames(LincardonT_2023_SpDf_modified),rownames(LincardonT_2023_aitchison_follow_up)),1]




### Calculate the correlation
corr_LincardonT_2023_bray <- corr.test(LincardonT_2023_SpDf,LincardonT_2023_bray_follow_up[rownames(LincardonT_2023_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")
corr_LincardonT_2023_aitchison <- corr.test(LincardonT_2023_SpDf,LincardonT_2023_aitchison_follow_up[rownames(LincardonT_2023_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")

df_corr_LincardonT_2023_bray <- data.frame("R"=corr_LincardonT_2023_bray$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_bray$r)),],"P"=corr_LincardonT_2023_bray$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_bray$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_bray$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_bray$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_bray$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_bray$r)),5])
df_corr_LincardonT_2023_aitchison <- data.frame("R"=corr_LincardonT_2023_aitchison$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_aitchison$r)),],"P"=corr_LincardonT_2023_aitchison$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_aitchison$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_aitchison$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_aitchison$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_aitchison$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_LincardonT_2023_aitchison$r)),5])

### add the directionnality and then select the species with negative and positive R
df_corr_LincardonT_2023_bray$dir <- ifelse(df_corr_LincardonT_2023_bray[,2]<= 0.05,2*sign(df_corr_LincardonT_2023_bray[,1]),sign(df_corr_LincardonT_2023_bray[,1]))
df_corr_LincardonT_2023_aitchison$dir <- ifelse(df_corr_LincardonT_2023_aitchison[,2]<= 0.05,2*sign(df_corr_LincardonT_2023_aitchison[,1]),sign(df_corr_LincardonT_2023_aitchison[,1]))

sig_negative_LincardonT_2023 <- union(rownames(df_corr_LincardonT_2023_aitchison[!is.na(df_corr_LincardonT_2023_aitchison[,7])&(df_corr_LincardonT_2023_aitchison[,7]==-2),]),rownames(df_corr_LincardonT_2023_bray[!is.na(df_corr_LincardonT_2023_bray[,7])&(df_corr_LincardonT_2023_bray[,7]==-2),]))
not_negative_LincardonT_2023 <- setdiff(union(rownames(df_corr_LincardonT_2023_aitchison),rownames(df_corr_LincardonT_2023_bray)),sig_negative_LincardonT_2023)



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_LincardonT_2023_bray_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_LincardonT_2023_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_LincardonT_2023_bray$R, df_corr_LincardonT_2023_bray$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_LincardonT_2023_bray$R and df_corr_LincardonT_2023_bray$Stability
# t = -10.254, df = 162, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7121556 -0.5245602
# sample estimates:
#        cor 
# -0.6273761

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_LincardonT_2023_bray_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_LincardonT_2023_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_LincardonT_2023_bray_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_LincardonT_2023_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_LincardonT_2023_bray_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_LincardonT_2023_bray,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_LincardonT_2023_aitchison_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_LincardonT_2023_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_LincardonT_2023_aitchison$R, df_corr_LincardonT_2023_aitchison$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_LincardonT_2023_aitchison$R and df_corr_LincardonT_2023_aitchison$Stability
# t = -6.7139, df = 162, p-value = 3.038e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5784534 -0.3374408
# sample estimates:
#        cor 
# -0.4665634

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_LincardonT_2023_aitchison_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_LincardonT_2023_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_LincardonT_2023_aitchison_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_LincardonT_2023_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_LincardonT_2023_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_LincardonT_2023_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()






#################################
#################################  NiuC_2020
NiuC_2020_metadata <- MetadataDf_saliva_validation %>% filter(study_name == "NiuC_2020")

NiuC_2020_metadata_followup <-
  NiuC_2020_metadata %>%
  filter(timepoint %in% c("1", "2", "3", "4")) %>%
  select(sample_id, subject_id, timepoint) %>%
  
  # Convert timepoint to ordered factor
  mutate(timepoint = factor(timepoint, levels = c("1","2","3","4"))) %>%
  
  arrange(subject_id, timepoint) %>%
  
  group_by(subject_id) %>%
  
  # Assign next sample within subject
  mutate(followup_sample = lead(sample_id)) %>%
  
  ungroup()


NiuC_2020_metadata_followup <- data.frame(NiuC_2020_metadata_followup)
rownames(NiuC_2020_metadata_followup) <-NiuC_2020_metadata_followup$sample_id
NiuC_2020_metadata_followup$subject_id <- NULL
NiuC_2020_metadata_followup$timepoint <- NULL
colnames(NiuC_2020_metadata_followup) <- c("T0","T1")

NiuC_2020_SpDf <- SpDf_saliva_Validation[rownames(NiuC_2020_metadata),]
NiuC_2020_SpDf <- NiuC_2020_SpDf[colSums(NiuC_2020_SpDf) >0]
NiuC_2020_SpDf <- NiuC_2020_SpDf[rowSums(NiuC_2020_SpDf)>0,]


### Calculate the Bray and Aitchison distance
NiuC_2020_bray_follow_up <- bray_followup(NiuC_2020_SpDf,NiuC_2020_metadata_followup)
NiuC_2020_aitchison_follow_up <- aitchison_followup(NiuC_2020_SpDf,NiuC_2020_metadata_followup)


## add the followup sample ids from T1 column of NiuC_2020_metadata_followup to new column in NiuC_2020_SpDf using maching rownames of NiuC_2020_SpDf and TO column in NiuC_2020_metadata_followup
NiuC_2020_SpDf_modified <- NiuC_2020_SpDf
NiuC_2020_SpDf_modified$followup_sample_id <- NiuC_2020_metadata_followup$T1[match(rownames(NiuC_2020_SpDf_modified),NiuC_2020_metadata_followup$T0)]

NiuC_2020_SpDf_modified$study_name <- "NiuC_2020"
NiuC_2020_SpDf_modified$timepoint <- NiuC_2020_metadata$timepoint[match(rownames(NiuC_2020_SpDf_modified),NiuC_2020_metadata$sample_id)]
NiuC_2020_SpDf_modified$study_condition <- NiuC_2020_metadata$study_condition[match(rownames(NiuC_2020_SpDf_modified),NiuC_2020_metadata$sample_id)]
NiuC_2020_SpDf_modified$bray_dist <- NiuC_2020_bray_follow_up[match(rownames(NiuC_2020_SpDf_modified),rownames(NiuC_2020_bray_follow_up)),1]
NiuC_2020_SpDf_modified$aitchison_dist <- NiuC_2020_aitchison_follow_up[match(rownames(NiuC_2020_SpDf_modified),rownames(NiuC_2020_aitchison_follow_up)),1]





### Calculate the correlation
corr_NiuC_2020_bray <- corr.test(NiuC_2020_SpDf,NiuC_2020_bray_follow_up[rownames(NiuC_2020_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")
corr_NiuC_2020_aitchison <- corr.test(NiuC_2020_SpDf,NiuC_2020_aitchison_follow_up[rownames(NiuC_2020_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")

df_corr_NiuC_2020_bray <- data.frame("R"=corr_NiuC_2020_bray$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_bray$r)),],"P"=corr_NiuC_2020_bray$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_bray$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_bray$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_bray$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_bray$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_bray$r)),5])
df_corr_NiuC_2020_aitchison <- data.frame("R"=corr_NiuC_2020_aitchison$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_aitchison$r)),],"P"=corr_NiuC_2020_aitchison$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_aitchison$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_aitchison$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_aitchison$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_aitchison$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_NiuC_2020_aitchison$r)),5])

### add the directionnality and then select the species with negative and positive R
df_corr_NiuC_2020_bray$dir <- ifelse(df_corr_NiuC_2020_bray[,2]<= 0.05,2*sign(df_corr_NiuC_2020_bray[,1]),sign(df_corr_NiuC_2020_bray[,1]))
df_corr_NiuC_2020_aitchison$dir <- ifelse(df_corr_NiuC_2020_aitchison[,2]<= 0.05, 2*sign(df_corr_NiuC_2020_aitchison[,1]),sign(df_corr_NiuC_2020_aitchison[,1]))


sig_negative_NiuC_2020 <- union(rownames(df_corr_NiuC_2020_aitchison[!is.na(df_corr_NiuC_2020_aitchison[,7])&(df_corr_NiuC_2020_aitchison[,7]==-2),]),rownames(df_corr_NiuC_2020_bray[!is.na(df_corr_NiuC_2020_bray[,7])&(df_corr_NiuC_2020_bray[,7]==-2),]))
not_negative_NiuC_2020 <- setdiff(union(rownames(df_corr_NiuC_2020_aitchison),rownames(df_corr_NiuC_2020_bray)),sig_negative_NiuC_2020)



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_NiuC_2020_bray_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_NiuC_2020_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_NiuC_2020_bray$R, df_corr_NiuC_2020_bray$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_NiuC_2020_bray$R and df_corr_NiuC_2020_bray$Stability
# t = -1.7665, df = 394, p-value = 0.07809
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.185568511  0.009990185
# sample estimates:
#         cor 
# -0.08864331


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_NiuC_2020_bray_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_NiuC_2020_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_NiuC_2020_bray_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_NiuC_2020_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_NiuC_2020_bray_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_NiuC_2020_bray,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()






pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_NiuC_2020_aitchison_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_NiuC_2020_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_NiuC_2020_aitchison$R, df_corr_NiuC_2020_aitchison$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_NiuC_2020_aitchison$R and df_corr_NiuC_2020_aitchison$Stability
# t = -3.9876, df = 394, p-value = 7.957e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2898761 -0.1003578
# sample estimates:
#        cor 
# -0.1969562


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_NiuC_2020_aitchison_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_NiuC_2020_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_NiuC_2020_aitchison_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_NiuC_2020_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_NiuC_2020_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_NiuC_2020_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()










#################################
#################################  StahringerS_2012
StahringerS_2012_metadata <- MetadataDf_saliva_validation %>% filter(study_name == "StahringerS_2012")

StahringerS_2012_metadata_followup <-
  StahringerS_2012_metadata %>%
  filter(timepoint %in% c("1", "2", "3", "4")) %>%
  select(sample_id, subject_id, timepoint) %>%
  
  # Convert timepoint to ordered factor
  mutate(timepoint = factor(timepoint, levels = c("1","2","3","4"))) %>%
  
  arrange(subject_id, timepoint) %>%
  
  group_by(subject_id) %>%
  
  # Assign next sample within subject
  mutate(followup_sample = lead(sample_id)) %>%
  
  ungroup()



StahringerS_2012_metadata_followup <- data.frame(StahringerS_2012_metadata_followup)
rownames(StahringerS_2012_metadata_followup) <-StahringerS_2012_metadata_followup$sample_id
StahringerS_2012_metadata_followup$subject_id <- NULL
StahringerS_2012_metadata_followup$timepoint <- NULL
colnames(StahringerS_2012_metadata_followup) <- c("T0","T1")

StahringerS_2012_SpDf <- SpDf_saliva_Validation[rownames(StahringerS_2012_metadata),]
StahringerS_2012_SpDf <- StahringerS_2012_SpDf[colSums(StahringerS_2012_SpDf) >0]


### Calculate the Bray and Aitchison distance
StahringerS_2012_bray_follow_up <- bray_followup(StahringerS_2012_SpDf,StahringerS_2012_metadata_followup)
StahringerS_2012_aitchison_follow_up <- aitchison_followup(StahringerS_2012_SpDf,StahringerS_2012_metadata_followup)


## add the followup sample ids from T1 column of StahringerS_2012_metadata_followup to new column in StahringerS_2012_SpDf using maching rownames of StahringerS_2012_SpDf and TO column in StahringerS_2012_metadata_followup
StahringerS_2012_SpDf_modified <- StahringerS_2012_SpDf
StahringerS_2012_SpDf_modified$followup_sample_id <- StahringerS_2012_metadata_followup$T1[match(rownames(StahringerS_2012_SpDf_modified),StahringerS_2012_metadata_followup$T0)]

StahringerS_2012_SpDf_modified$study_name <- "StahringerS_2012"
StahringerS_2012_SpDf_modified$timepoint <- StahringerS_2012_metadata$timepoint[match(rownames(StahringerS_2012_SpDf_modified),StahringerS_2012_metadata$sample_id)]
StahringerS_2012_SpDf_modified$study_condition <- StahringerS_2012_metadata$study_condition[match(rownames(StahringerS_2012_SpDf_modified),StahringerS_2012_metadata$sample_id)]
StahringerS_2012_SpDf_modified$bray_dist <- StahringerS_2012_bray_follow_up[match(rownames(StahringerS_2012_SpDf_modified),rownames(StahringerS_2012_bray_follow_up)),1]
StahringerS_2012_SpDf_modified$aitchison_dist <- StahringerS_2012_aitchison_follow_up[match(rownames(StahringerS_2012_SpDf_modified),rownames(StahringerS_2012_aitchison_follow_up)),1]




### Calculate the correlation
corr_StahringerS_2012_bray <- corr.test(StahringerS_2012_SpDf,StahringerS_2012_bray_follow_up[rownames(StahringerS_2012_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")
corr_StahringerS_2012_aitchison <- corr.test(StahringerS_2012_SpDf,StahringerS_2012_aitchison_follow_up[rownames(StahringerS_2012_SpDf),],method="spearman",use="pairwise.complete",adjust="fdr")

df_corr_StahringerS_2012_bray <- data.frame("R"=corr_StahringerS_2012_bray$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_bray$r)),],"P"=corr_StahringerS_2012_bray$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_bray$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_bray$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_bray$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_bray$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_bray$r)),5])
df_corr_StahringerS_2012_aitchison <- data.frame("R"=corr_StahringerS_2012_aitchison$r[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_aitchison$r)),],"P"=corr_StahringerS_2012_aitchison$p[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_aitchison$r)),],"Influence"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_aitchison$r)),1],"Stability"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_aitchison$r)),3],"Health"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_aitchison$r)),2],"HACKScore"=Combined_Saliva_Scores[intersect(rownames(Combined_Saliva_Scores),rownames(corr_StahringerS_2012_aitchison$r)),5])

### add the directionnality and then select the species with negative and positive R
df_corr_StahringerS_2012_bray$dir <- ifelse(df_corr_StahringerS_2012_bray[,2]<= 0.05,2*sign(df_corr_StahringerS_2012_bray[,1]),sign(df_corr_StahringerS_2012_bray[,1]))
df_corr_StahringerS_2012_aitchison$dir <- ifelse(df_corr_StahringerS_2012_aitchison[,2]<= 0.05,2*sign(df_corr_StahringerS_2012_aitchison[,1]),sign(df_corr_StahringerS_2012_aitchison[,1]))


sig_negative_StahringerS_2012 <- union(rownames(df_corr_StahringerS_2012_aitchison[!is.na(df_corr_StahringerS_2012_aitchison[,7])&(df_corr_StahringerS_2012_aitchison[,7]==-2),]),rownames(df_corr_StahringerS_2012_bray[!is.na(df_corr_StahringerS_2012_bray[,7])&(df_corr_StahringerS_2012_bray[,7]==-2),]))
not_negative_StahringerS_2012 <- setdiff(union(rownames(df_corr_StahringerS_2012_aitchison),rownames(df_corr_StahringerS_2012_bray)),sig_negative_StahringerS_2012)



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_StahringerS_2012_bray_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_StahringerS_2012_bray,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_StahringerS_2012_bray$R, df_corr_StahringerS_2012_bray$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_StahringerS_2012_bray$R and df_corr_StahringerS_2012_bray$Stability
# t = 1.8534, df = 198, p-value = 0.06531
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.008301559  0.264537937
# sample estimates:
#       cor 
# 0.1305899


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_StahringerS_2012_bray_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_StahringerS_2012_bray,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_StahringerS_2012_bray_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_StahringerS_2012_bray,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_StahringerS_2012_bray_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_StahringerS_2012_bray,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_StahringerS_2012_aitchison_Stability_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_StahringerS_2012_aitchison,aes(x=R,y=Stability))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

cor.test(df_corr_StahringerS_2012_aitchison$R, df_corr_StahringerS_2012_aitchison$Stability, method = "pearson")
#         Pearson's product-moment correlation

# data:  df_corr_StahringerS_2012_aitchison$R and df_corr_StahringerS_2012_aitchison$Stability
# t = 1.8316, df = 198, p-value = 0.06851
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.009838633  0.263107730
# sample estimates:
#       cor 
# 0.1290786


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_StahringerS_2012_aitchison_Health_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_StahringerS_2012_aitchison,aes(x=R,y=Health))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_StahringerS_2012_aitchison_CoreInfluence_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_StahringerS_2012_aitchison,aes(x=R,y=Influence))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3corr_StahringerS_2012_aitchison_HACKScore_R.pdf", height = 4.5, width = 5)
ggplot(df_corr_StahringerS_2012_aitchison,aes(x=R,y=HACKScore))+geom_point(size = 2.5)+geom_smooth(method="lm")+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()







#### combine all the modified species profiles (added metadata columns and followup distances)
library(dplyr)
combined_Saliva_Longitudinal_SpDf_modified <- bind_rows(ZhangT_2020_SpDf_modified,SchmidtT_2019B_SpDf_modified,ChenC_2018_SpDf_modified,LincardonT_2023_SpDf_modified,StahringerS_2012_SpDf_modified,NiuC_2020_SpDf_modified)

# replace NAs with 0 in all columns except for followup_sample_id, study_name, timepoint, study_condition, bray_dist, and aitchison_dist
combined_Saliva_Longitudinal_SpDf_modified[, !names(combined_Saliva_Longitudinal_SpDf_modified) %in% c("followup_sample_id","study_name","timepoint","study_condition","bray_dist","aitchison_dist")] <-
  lapply(
    combined_Saliva_Longitudinal_SpDf_modified[, !names(combined_Saliva_Longitudinal_SpDf_modified) %in% c("followup_sample_id","study_name","timepoint","study_condition","bray_dist","aitchison_dist")],
    function(x) replace(x, is.na(x), 0))


# add new column distance_MeanRanked which is the mean of two columns of the distance after taking its rankscaled values. This has to be done study-wise
rank_scale=function(x){
  # x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

combined_Saliva_Longitudinal_SpDf_modified <- combined_Saliva_Longitudinal_SpDf_modified %>%
  group_by(study_name) %>%
  mutate(distance_MeanRanked = rowMeans(cbind(rank_scale(bray_dist), rank_scale(aitchison_dist)), na.rm = TRUE)) %>%
  ungroup()

save(combined_Saliva_Longitudinal_SpDf_modified, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3combined_Saliva_Longitudinal_SpDf_modified.RData")











##############
############## Bean plots all studies
library(beanplot)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3beanplot_HACKScore_ZhangT2020_SchmidtT2019B_ChenC_2018.pdf", height = 6, width = 5)
beanplot(Combined_Saliva_Scores[sig_negative_ZhangT_2020,5],Combined_Saliva_Scores[not_negative_ZhangT_2020,5],Combined_Saliva_Scores[sig_negative_SchmidtT_2019B,5],
          Combined_Saliva_Scores[not_negative_SchmidtT_2019B,5],Combined_Saliva_Scores[sig_negative_ChenC_2018,5],Combined_Saliva_Scores[not_negative_ChenC_2018,5], 
          Combined_Saliva_Scores[sig_negative_LincardonT_2023,5],Combined_Saliva_Scores[not_negative_LincardonT_2023,5], Combined_Saliva_Scores[sig_negative_StahringerS_2012,5],Combined_Saliva_Scores[not_negative_StahringerS_2012,5],
          Combined_Saliva_Scores[sig_negative_NiuC_2020,5],Combined_Saliva_Scores[not_negative_NiuC_2020,5],
          names = c("ZhangT_2020","SchmidtT_2019B","ChenC_2018", "LincardonT_2023","StahringerS_2012","NiuC_2020"),
          side="both",what=c(0,1,1,0),overallline="median",col=list("deepskyblue","gold"), las = 2)
dev.off()

pval_bean_hack1 <- wilcox.test(Combined_Saliva_Scores[sig_negative_ZhangT_2020,5],Combined_Saliva_Scores[not_negative_ZhangT_2020,5])
pval_bean_hack1$p.value
# [1] 4.46314e-21

pval_bean_hack2 <- wilcox.test(Combined_Saliva_Scores[sig_negative_SchmidtT_2019B,5],Combined_Saliva_Scores[not_negative_SchmidtT_2019B,5])
pval_bean_hack2$p.value
# [1] 0.0001036615

pval_bean_hack3 <- wilcox.test(Combined_Saliva_Scores[sig_negative_ChenC_2018,5],Combined_Saliva_Scores[not_negative_ChenC_2018,5])
pval_bean_hack3$p.value
# 0.008033633

pval_bean_hack4 <- wilcox.test(Combined_Saliva_Scores[sig_negative_LincardonT_2023,5],Combined_Saliva_Scores[not_negative_LincardonT_2023,5])
pval_bean_hack4$p.value
# 5.117804e-11

pval_bean_hack5 <- wilcox.test(Combined_Saliva_Scores[sig_negative_StahringerS_2012,5],Combined_Saliva_Scores[not_negative_StahringerS_2012,5])
pval_bean_hack5$p.value
# 0.8171662

pval_bean_hack6 <- wilcox.test(Combined_Saliva_Scores[sig_negative_NiuC_2020,5],Combined_Saliva_Scores[not_negative_NiuC_2020,5])
pval_bean_hack6$p.value
# 0.0001022582


library(beanplot)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3beanplot_StabilityScore_ZhangT2020_SchmidtT2019B_ChenC_2018.pdf", height = 6, width = 5)
beanplot(Combined_Saliva_Scores[sig_negative_ZhangT_2020,3],Combined_Saliva_Scores[not_negative_ZhangT_2020,3],Combined_Saliva_Scores[sig_negative_SchmidtT_2019B,3],Combined_Saliva_Scores[not_negative_SchmidtT_2019B,3],
          Combined_Saliva_Scores[sig_negative_ChenC_2018,3],Combined_Saliva_Scores[not_negative_ChenC_2018,3], Combined_Saliva_Scores[sig_negative_LincardonT_2023,3],Combined_Saliva_Scores[not_negative_LincardonT_2023,3],
          Combined_Saliva_Scores[sig_negative_StahringerS_2012,3],Combined_Saliva_Scores[not_negative_StahringerS_2012,3],Combined_Saliva_Scores[sig_negative_NiuC_2020,3],Combined_Saliva_Scores[not_negative_NiuC_2020,3],
          names = c("ZhangT_2020","SchmidtT_2019B","ChenC_2018", "LincardonT_2023","StahringerS_2012","NiuC_2020"),
          side="both",what=c(0,1,1,0),overallline="median",col=list("deepskyblue","gold"), las = 2)
dev.off()

pval_bean_stability1 <- wilcox.test(Combined_Saliva_Scores[sig_negative_ZhangT_2020,3],Combined_Saliva_Scores[not_negative_ZhangT_2020,3])
pval_bean_stability1$p.value
# [1] 5.536749e-11

pval_bean_stability2 <- wilcox.test(Combined_Saliva_Scores[sig_negative_SchmidtT_2019B,3],Combined_Saliva_Scores[not_negative_SchmidtT_2019B,3])
pval_bean_stability2$p.value
# [1] 9.079828e-05

pval_bean_stability3 <- wilcox.test(Combined_Saliva_Scores[sig_negative_ChenC_2018,3],Combined_Saliva_Scores[not_negative_ChenC_2018,3])
pval_bean_stability3$p.value
# [1] 0.2012936

pval_bean_stability4 <- wilcox.test(Combined_Saliva_Scores[sig_negative_LincardonT_2023,3],Combined_Saliva_Scores[not_negative_LincardonT_2023,3])
pval_bean_stability4$p.value
#  5.311385e-12

pval_bean_stability5 <- wilcox.test(Combined_Saliva_Scores[sig_negative_StahringerS_2012,3],Combined_Saliva_Scores[not_negative_StahringerS_2012,3])
pval_bean_stability5$p.value
# 0.008053771

pval_bean_stability6 <- wilcox.test(Combined_Saliva_Scores[sig_negative_NiuC_2020,3],Combined_Saliva_Scores[not_negative_NiuC_2020,3])
pval_bean_stability6$p.value
# 0.2157235

#############
############# Now plot a box plot of stability association of speceis coming as stability associated in number of datasets
# All species
all_species <- rownames(Combined_Saliva_Scores)

# Logical vectors
neg_Zhang   <- all_species %in% sig_negative_ZhangT_2020
neg_Schmidt <- all_species %in% sig_negative_SchmidtT_2019B
neg_ChenC   <- all_species %in% sig_negative_ChenC_2018
neg_Lin <- all_species %in% sig_negative_LincardonT_2023
neg_Stahringer <- all_species %in% sig_negative_StahringerS_2012
neg_NiuC <- all_species %in% sig_negative_NiuC_2020

# Count number of datasets where negatively associated
neg_count <- neg_Zhang + neg_Schmidt  + neg_ChenC + neg_Lin + neg_Stahringer + neg_NiuC # TRUE=1, FALSE=0
neg_count_cat <- ifelse(neg_count == 0, "0",
                  ifelse(neg_count == 1 | neg_count == 2, "1-2",
                  ifelse(neg_count == 3 | neg_count == 4, "3-4", NA)))

# Build plotting dataframe
df_plot <- data.frame(
  Species = all_species,
  Stability = Combined_Saliva_Scores[,3],
  HACKScore = Combined_Saliva_Scores[,5],
  NegCount = neg_count_cat
)

df_plot$NegCount <- factor(df_plot$NegCount, levels = c("0","1-2","3-4"))

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3StabilityAssociation_Studies_Vs_HACK_Score.pdf",
    height = 6, width = 8)

boxplot(
  HACKScore ~ NegCount,
  data = df_plot,
  col = c("#f28e8e", "#ffd92f", "#72bf9c", "#a6d854", "#e78ac3"),   # two different colors
  border = "black",
  ylab = "HACK score",
  xlab = "Number of studies with Stability association",
  outline = FALSE
)

dev.off()

library(dunn.test)
Dunntest_HACKscore <- dunn.test(df_plot$HACKScore, df_plot$NegCount, method = "bh")
Dunntest_HACKscore
# $chi2
# [1] 182.8024

# $Z
# [1] -12.735242  -6.097259  -1.260241

# $P
# [1] 1.883378e-37 5.395142e-10 1.037912e-01

# $P.adjusted
# [1] 5.650133e-37 8.092713e-10 1.037912e-01

# $comparisons
# [1] "0 - 1-2"   "0 - 3-4"   "1-2 - 3-4"


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3StabilityAssociation_Studies_Vs_Stability_Score.pdf", height =6, width = 8)
boxplot(
  Stability ~ NegCount,
  data = df_plot,
  col = c("#f28e8e", "#ffd92f", "#72bf9c", "#a6d854", "#e78ac3"),   # two different colors
  border = "black",
  ylab = "Stability score",
  xlab = "Number of studies with Stability association",
  outline = FALSE)
dev.off()

Dunntest_Stabilityscore <- dunn.test(df_plot$Stability, df_plot$NegCount, method = "bh")
Dunntest_Stabilityscore
# $chi2
# [1] 84.72855

# $Z
# [1] -8.6944987 -4.0857249 -0.7852223

# $P
# [1] 1.741822e-18 2.196971e-05 2.161616e-01

# $P.adjusted
# [1] 5.225467e-18 3.295456e-05 2.161616e-01

# $comparisons
# [1] "0 - 1-2"   "0 - 3-4"   "1-2 - 3-4"


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_3StabilityAssociation_Validation_Workspace.RData")


