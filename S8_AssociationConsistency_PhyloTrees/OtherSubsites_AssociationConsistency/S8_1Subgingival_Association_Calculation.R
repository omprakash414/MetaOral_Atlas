


S3_1Subgingival_Core_Detection_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Subgingival_Core_Detection_Workspace.RData", envir = S3_1Subgingival_Core_Detection_Workspace)
# load("S3_1Subgingival_Core_Detection_Workspace.RData", envir = S3_1Subgingival_Core_Detection_Workspace)
attach(S3_1Subgingival_Core_Detection_Workspace)
SpDf_subgingival_Control_associated <- SpDf_subgingival_Control_associated
MetadataDf_subgingival_control <- MetadataDf_subgingival_control
detach(S3_1Subgingival_Core_Detection_Workspace)
rm(S3_1Subgingival_Core_Detection_Workspace)

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")
# source("code_library_Metaoral.R")

### Remove the study_name column from species profile and then normalise them and then add the study_name column
SpDf_subgingival_Control_associated$study_name <- NULL
SpDf_subgingival_Control_associated <- SpDf_subgingival_Control_associated[rowSums(SpDf_subgingival_Control_associated) > 0,]
SpDf_subgingival_Control_associated <- SpDf_subgingival_Control_associated/rowSums(SpDf_subgingival_Control_associated)
SpDf_subgingival_Control_associated[is.na(SpDf_subgingival_Control_associated)] <- 0

MetadataDf_subgingival_control <- MetadataDf_subgingival_control[rownames(SpDf_subgingival_Control_associated),]
SpDf_subgingival_Control_associated$study_name <- MetadataDf_subgingival_control$study_name


######### Run the Function to get embeddings or associations
study_list <- unique(SpDf_subgingival_Control_associated$study_name)
studies <- length(unique(SpDf_subgingival_Control_associated$study_name))
species_list <- setdiff(colnames(SpDf_subgingival_Control_associated),"study_name")

subgingival_association_matrix <- compute_association_consistency_parallel(SpDf_subgingival_Control_associated,species_list,study_list,studies)

######### Get the PCoA coordinates based on correlation as well as based on the distance.

pca_subgingival_association <- dudi.pca(subgingival_association_matrix$final_matrix,scannf=FALSE,nf=3)

save(pca_subgingival_association, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1PCA_subgingival_associations.RData")
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_Association_Workspace.RData")
