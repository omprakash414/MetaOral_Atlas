

S3_1Supragingival_Core_Detection_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Supragingival_Core_Detection_Workspace.RData", envir = S3_1Supragingival_Core_Detection_Workspace)
# load("S3_1Supragingival_Core_Detection_Workspace.RData", envir = S3_1Supragingival_Core_Detection_Workspace)
attach(S3_1Supragingival_Core_Detection_Workspace)
MetadataDf_supragingival_control <- MetadataDf_supragingival_control
SpDf_supragingival_Control_associated <- SpDf_supragingival_Control_associated
detach(S3_1Supragingival_Core_Detection_Workspace)
rm(S3_1Supragingival_Core_Detection_Workspace)


########### Prepare the data for network construction
# remove the rows that have zero rowsums in species profile
SpDf_supragingival_Control_associated <- SpDf_supragingival_Control_associated[rowSums(SpDf_supragingival_Control_associated[,colnames(SpDf_supragingival_Control_associated) != "study_name"],na.rm = TRUE) > 0,]
# normalize the data except study_name column
SpDf_supragingival_Control_associated[,colnames(SpDf_supragingival_Control_associated) != "study_name"] <- SpDf_supragingival_Control_associated[,colnames(SpDf_supragingival_Control_associated) != "study_name"] / rowSums(SpDf_supragingival_Control_associated[,colnames(SpDf_supragingival_Control_associated) != "study_name"],na.rm = TRUE)
SpDf_supragingival_Control_associated[is.na(SpDf_supragingival_Control_associated)] <- 0

species_to_check <- colnames(SpDf_supragingival_Control_associated)[colnames(SpDf_supragingival_Control_associated) != "study_name"]

# Now filter the metadata for the same 
MetadataDf_supragingival_control <- MetadataDf_supragingival_control[rownames(SpDf_supragingival_Control_associated),]


########### Start Network code on control only cohort
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")

#source("code_library_Metaoral.R")
RemNetwork_supragingival_Control <- Rem_Network2(SpDf_supragingival_Control_associated,species_to_check,"study_name",unique(SpDf_supragingival_Control_associated$study_name),species_to_check)


save(RemNetwork_supragingival_Control, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_4supragingival_ControlCohort_RemNetwork.RData")
#save(RemNetwork_supragingival_Control, file = "S7_4supragingival_ControlCohort_RemNetwork.RData")



####### Extract the dir table and then melt it to get relationship of species pairs which will be further used for plotting the network
supragingival_ControlCohort_RemNetwork_melt <- Melt_Adjacency_Matrix(RemNetwork_supragingival_Control$dir)
supragingival_ControlCohort_RemNetwork_melt_filt <- subset(supragingival_ControlCohort_RemNetwork_melt, Weight == 1)

# 12,403 edges
write.csv(supragingival_ControlCohort_RemNetwork_melt_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_4supragingival_ControlCohort_RemNetwork_melt_filt.csv")
#write.csv(supragingival_ControlCohort_RemNetwork_melt_filt, file = "S7_4supragingival_ControlCohort_RemNetwork_melt_filt.csv")

####### Detect the clusters in the network based on edge betweenness using walktrap algorithm
library(igraph)
supragingival_ControlCohort_graph <- graph_from_edgelist(as.matrix(supragingival_ControlCohort_RemNetwork_melt_filt[,1:2]))
supragingival_ControlCohort_walktrap_clusters <- cluster_walktrap(supragingival_ControlCohort_graph)
supragingival_ControlCohort_walktrap_clusters_df <- data.frame(species = names(membership(supragingival_ControlCohort_walktrap_clusters)),
                                                               ClusterID = as.numeric(membership(supragingival_ControlCohort_walktrap_clusters)),
                                                               row.names = NULL)

write.csv(supragingival_ControlCohort_walktrap_clusters_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_4supragingival_ControlCohort_walktrap_clusters_df.csv")
# write.csv(supragingival_ControlCohort_walktrap_clusters_df, file = "S7_4supragingival_ControlCohort_walktrap_clusters_df.csv")




save.image("S7_4supragingival_ControlCohort_Workspace.RData")













# ######## There are many edges in the network if we take qvalue as 0.05. so reduce it to 0.001 and rebuild the network
# qval_matrix <- data.frame(ifelse(RemNetwork_supragingival_Control$qval <= 0.0001, 1, 0))
# consistency_matrix <- data.frame(ifelse(RemNetwork_supragingival_Control$consistency >= 0.70, 1, 0))
# est_matrix <- data.frame(ifelse(RemNetwork_supragingival_Control$est > 0, 1, ifelse(RemNetwork_supragingival_Control$est < 0, -1, 0)))

# # Create dir_matrix based on conditions
# dir_matrix  <- (qval_matrix * consistency_matrix * est_matrix)

# # now melt it
# supragingival_ControlCohort_dir_matrix <- Melt_Adjacency_Matrix(dir_matrix)
# supragingival_ControlCohort_dir_matrix_filt <- subset(supragingival_ControlCohort_dir_matrix, Weight == 1 | Weight == -1)
# # edges are not 2504
# write.csv(supragingival_ControlCohort_dir_matrix_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_4supragingival_ControlCohort_RemNetwork_melt_filt_0.0001.csv")
# write.csv(supragingival_ControlCohort_dir_matrix_filt, file = "S7_4supragingival_ControlCohort_RemNetwork_melt_filt_0.0001.csv")


# # Now see the cluster numbers using walktrap algorithm
# supragingival_ControlCohort_dir_graph <- graph_from_edgelist(as.matrix(supragingival_ControlCohort_dir_matrix_filt[,1:2]))
# supragingival_ControlCohort_dir_walktrap_clusters <- cluster_walktrap(supragingival_ControlCohort_dir_graph)
# supragingival_ControlCohort_dir_walktrap_clusters_df <- data.frame(species = names(membership(supragingival_ControlCohort_dir_walktrap_clusters)),
#                                                                         ClusterID = as.numeric(membership(supragingival_ControlCohort_dir_walktrap_clusters)),
#                                                                         row.names = NULL)
# write.csv(supragingival_ControlCohort_dir_walktrap_clusters_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_4supragingival_ControlCohort_walktrap_clusters_df_0.0001.csv")



save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_4supragingival_ControlCohort_Workspace.RData")
# save.image("S7_4supragingival_ControlCohort_Workspace.RData")
