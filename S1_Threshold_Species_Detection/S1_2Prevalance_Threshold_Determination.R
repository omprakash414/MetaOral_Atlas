#### Threshold that can be used for prevalance and envfit analysis to define the core association
###


library(dplyr)
library(vegan)
library(reshape2)
library(pheatmap)
library(ggplot2)

############ Load the data
S1_Threshold_Species_Detection_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Threshold_Species_Detection_Workspace.RData",envir = S1_Threshold_Species_Detection_Workspace)
attach(S1_Threshold_Species_Detection_Workspace)
SpDf_saliva_Control<- SpDf_saliva_Control
list_sample_groups_saliva <- list_sample_groups_saliva
AllSpeciesDetectionPattern_saliva <- AllSpeciesDetectionPattern_saliva

SpDf_supragingival_Control <- SpDf_supragingival_Control
list_sample_groups_supragingival<- list_sample_groups_supragingival
AllSpeciesDetectionPattern_supragingival <- AllSpeciesDetectionPattern_supragingival

SpDf_tongue_tonsil_Control <- SpDf_tongue_tonsil_Control
list_sample_groups_tonguetonsil <- list_sample_groups_tonguetonsil
AllSpeciesDetectionPattern_tonguetonsil <- AllSpeciesDetectionPattern_tonguetonsil

SpDf_subgingival_Control <- SpDf_subgingival_Control
list_sample_groups_subgingival <- list_sample_groups_subgingival
AllSpeciesDetectionPattern_subgingival <- AllSpeciesDetectionPattern_subgingival

SpDf_buccal_palate_Control <- SpDf_buccal_palate_Control
list_sample_groups_buccal <- list_sample_groups_buccal
AllSpeciesDetectionPattern_buccal <- AllSpeciesDetectionPattern_buccal

detach(S1_Threshold_Species_Detection_Workspace)
rm(S1_Threshold_Species_Detection_Workspace)
gc()

############ Load the selected species for each subsite
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")

############ Load the code library functions
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")





############ Run for Saliva
saliva_jacc_representation <- compute_detection_jaccard_representation(AllSpeciesDetectionPattern_saliva,saliva_AssociatedSpecies,SpDf_saliva_Control,list_sample_groups_saliva)

df_jaccard_saliva       <- saliva_jacc_representation$df_jaccard
df_representation_saliva <- saliva_jacc_representation$df_representation
df_patterns_saliva       <- saliva_jacc_representation$df_patterns

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_2Jaccard_Representation_Saliva_Prevalance_Threshold.pdf",width=10,height=6)
ggplot(df_patterns_saliva)+geom_boxplot(aes(y=Jaccard_Similarity,x=Threshold,group=Threshold),fill="#1B9E77",alpha=0.5)+geom_boxplot(aes(y=Representation,x=Threshold,group=Threshold),fill="#D95F02",alpha=0.5)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),panel.grid.major=element_line(linewidth=0.5,color="grey70"))
dev.off()





############ Run for Supragingival
supragingival_jacc_representation <- compute_detection_jaccard_representation(AllSpeciesDetectionPattern_supragingival,supragingival_AssociatedSpecies,SpDf_supragingival_Control,list_sample_groups_supragingival)

df_jaccard_supragingival       <- supragingival_jacc_representation$df_jaccard
df_representation_supragingival <- supragingival_jacc_representation$df_representation
df_patterns_supragingival  <- supragingival_jacc_representation$df_patterns

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_2Jaccard_Representation_Supragingival_Prevalance_Threshold.pdf",width=10,height=6)
ggplot(df_patterns_supragingival)+geom_boxplot(aes(y=Jaccard_Similarity,x=Threshold,group=Threshold),fill="#1B9E77",alpha=0.5)+geom_boxplot(aes(y=Representation,x=Threshold,group=Threshold),fill="#D95F02",alpha=0.5)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),panel.grid.major=element_line(linewidth=0.5,color="grey70"))
dev.off()





############ Run for tongue_tonsil
tongue_tonsil_jacc_representation <- compute_detection_jaccard_representation(AllSpeciesDetectionPattern_tonguetonsil,tonguetonsil_AssociatedSpecies,SpDf_tongue_tonsil_Control,list_sample_groups_tonguetonsil)

df_jaccard_tongue_tonsil       <- tongue_tonsil_jacc_representation$df_jaccard
df_representation_tongue_tonsil <- tongue_tonsil_jacc_representation$df_representation
df_patterns_tongue_tonsil  <- tongue_tonsil_jacc_representation$df_patterns

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_2Jaccard_Representation_Tongue_Tonsil_Prevalance_Threshold.pdf",width=10,height=6)
ggplot(df_patterns_tongue_tonsil)+geom_boxplot(aes(y=Jaccard_Similarity,x=Threshold,group=Threshold),fill="#1B9E77",alpha=0.5)+geom_boxplot(aes(y=Representation,x=Threshold,group=Threshold),fill="#D95F02",alpha=0.5)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),panel.grid.major=element_line(linewidth=0.5,color="grey70"))
dev.off()






############ Run for subgingival
subgingival_jacc_representation <- compute_detection_jaccard_representation(AllSpeciesDetectionPattern_subgingival,subgingival_AssociatedSpecies,SpDf_subgingival_Control,list_sample_groups_subgingival)

df_jaccard_subgingival       <- subgingival_jacc_representation$df_jaccard
df_representation_subgingival <- subgingival_jacc_representation$df_representation
df_patterns_subgingival  <- subgingival_jacc_representation$df_patterns

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_2Jaccard_Representation_Subgingival_Prevalance_Threshold.pdf",width=10,height=6)
ggplot(df_patterns_subgingival)+geom_boxplot(aes(y=Jaccard_Similarity,x=Threshold,group=Threshold),fill="#1B9E77",alpha=0.5)+geom_boxplot(aes(y=Representation,x=Threshold,group=Threshold),fill="#D95F02",alpha=0.5)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),panel.grid.major=element_line(linewidth=0.5,color="grey70"))
dev.off()






############ Run for buccal palate
buccal_palate_jacc_representation <- compute_detection_jaccard_representation(AllSpeciesDetectionPattern_buccal,buccal_AssociatedSpecies,SpDf_buccal_palate_Control,list_sample_groups_buccal)

df_jaccard_buccal_palate       <- buccal_palate_jacc_representation$df_jaccard
df_representation_buccal_palate <- buccal_palate_jacc_representation$df_representation
df_patterns_buccal_palate  <- buccal_palate_jacc_representation$df_patterns

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_2Jaccard_Representation_Buccal_Palate_Prevalance_Threshold.pdf",width=10,height=6)
ggplot(df_patterns_buccal_palate)+geom_boxplot(aes(y=Jaccard_Similarity,x=Threshold,group=Threshold),fill="#1B9E77",alpha=0.5)+geom_boxplot(aes(y=Representation,x=Threshold,group=Threshold),fill="#D95F02",alpha=0.5)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),panel.grid.major=element_line(linewidth=0.5,color="grey70"))
dev.off()




save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_2Prevalence_Threshold_Determination_Workspace.RData")










#####################################################################################
#######################################################################################


