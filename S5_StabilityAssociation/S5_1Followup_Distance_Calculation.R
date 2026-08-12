
library(dplyr)
########### Load the data to calculate the longitudinal distances (Before that we will clean the data)
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_LongCohort.RData")

## Species for Saliva
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")
rm(list = ls(pattern = "subgingival|supragingival|tongue|buccal"))

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")

########### See all the columns are okay and filter the speceis for saliva (selected in S1 - Threshold speceis)

all(nrow(MetadataDf_saliva_Long) == nrow(SpDf_saliva_Long))
all(rownames(MetadataDf_saliva_Long) == rownames(SpDf_saliva_Long))
SpDf_saliva_Long <- SpDf_saliva_Long[,intersect(colnames(SpDf_saliva_Long),saliva_AssociatedSpecies)]


########## 
Longitudinal_Dist_Metadata <- Longitudinal_Microbiome_Distance(MetadataDf_saliva_Long, SpDf_saliva_Long)

######### add the distance calculated alomg with standard metadata columns to species profile
std_meta_cols <- c("Subject_ID","Sample_ID","Time_Point","Follow_up","Time_point_type",
                   "Treatment_Condition","Timepoint_difference","Study_Name",
                   "Jaccard_dist","Aitchison_dist","Kendall_dist","BrayCurtis_dist")

# confirm the rownames in Longitudinal_Dist_Metadata and species profile are alighed in same order
all(rownames(Longitudinal_Dist_Metadata) == rownames(SpDf_saliva_Long))

# Add the standard columns to species profile
SpDf_saliva_Long_withDist <- bind_cols(SpDf_saliva_Long, Longitudinal_Dist_Metadata[,std_meta_cols])


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_1Followup_Distance_calculation_Workspace.RData")



