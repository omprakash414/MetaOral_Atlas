######## This function calculates the embeddings/associations within the species profiles


######### Load the data
S3_1Saliva_Core_Detection_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_Core_Detection_Workspace.RData", envir = S3_1Saliva_Core_Detection_Workspace)
# load("S3_1Saliva_Core_Detection_Workspace.RData", envir = S3_1Saliva_Core_Detection_Workspace)
attach(S3_1Saliva_Core_Detection_Workspace)
SpDf_saliva_Control_associated <- SpDf_saliva_Control_associated
MetadataDf_saliva_control <- MetadataDf_saliva_control
detach(S3_1Saliva_Core_Detection_Workspace)
rm(S3_1Saliva_Core_Detection_Workspace)

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")
# source("code_library_Metaoral.R")

### See if the specie sprofile is normalised
row_sums <- rowSums(
  SpDf_saliva_Control_associated[
    , colnames(SpDf_saliva_Control_associated) != "study_name"
  ],
  na.rm = TRUE
)
summary(row_sums)
# species profile is already normalised and study_name column is already in it. So, we can directly use it to calculate the associations.


######### Run the Function to get embeddings or associations
study_list <- unique(SpDf_saliva_Control_associated$study_name)
studies <- length(unique(SpDf_saliva_Control_associated$study_name))
species_list <- setdiff(colnames(SpDf_saliva_Control_associated),"study_name")

saliva_association_matrix <- compute_association_consistency_parallel(SpDf_saliva_Control_associated,species_list,study_list,studies)

######### Get the PCoA coordinates based on correlation as well as based on the distance.
pca_saliva_association <- dudi.pca(saliva_association_matrix$final_matrix,scannf=FALSE,nf=3)


save(pca_saliva_association, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1PCA_saliva_associations.RData")
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_Association_Workspace.RData")
