
#### Import the specie sprofile:

load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Species_Profiles_Associated_Species_SubsiteWise.RData")

rm(buccal_AssociatedSpecies,SpDf_buccal_palate_Control_associated,SpDf_subgingival_Control_associated,SpDf_supragingival_Control_associated,SpDf_tongue_tonsil_Control_associated,subgingival_AssociatedSpecies,supragingival_AssociatedSpecies,tonguetonsil_AssociatedSpecies)

#### Import the metadata
S1_Threshold_Species_Detection_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Threshold_Species_Detection_Workspace.RData", envir = S1_Threshold_Species_Detection_Workspace)
attach(S1_Threshold_Species_Detection_Workspace)

metadata_objs <- ls(envir = S1_Threshold_Species_Detection_Workspace, pattern = "^MetadataDf_saliva")
list2env(mget(metadata_objs, envir = S1_Threshold_Species_Detection_Workspace),envir = .GlobalEnv)

detach(S1_Threshold_Species_Detection_Workspace)
rm(S1_Threshold_Species_Detection_Workspace)
rm(metadata_objs)



#### Check the species profile and metadata are in same order
all(rownames(SpDf_saliva_Control_associated) == rownames(MetadataDf_saliva_control))
all(nrow(SpDf_saliva_Control_associated) == nrow(MetadataDf_saliva_control))


species_to_ckeck <- colnames(SpDf_saliva_Control_associated)


###### First normalize the profile and then add the study_name column
SpDf_saliva_Control_associated <- SpDf_saliva_Control_associated[rowSums(SpDf_saliva_Control_associated) > 0,]
SpDf_saliva_Control_associated <- SpDf_saliva_Control_associated/rowSums(SpDf_saliva_Control_associated)

MetadataDf_saliva_control <- MetadataDf_saliva_control[rownames(SpDf_saliva_Control_associated),]

SpDf_saliva_Control_associated$study_name <- MetadataDf_saliva_control$study_name

save(SpDf_saliva_Control_associated, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/Control_Associated_Data_260113.RData")
save(MetadataDf_saliva_control, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/Control_Associated_Metadata_260113.RData")

######### Start Network code on control only cohort
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")
RemNetwork_saliva_Control <- Rem_Network2(SpDf_saliva_Control_associated,species_to_ckeck,"study_name",unique(SpDf_saliva_Control_associated$study_name),species_to_ckeck)



save(RemNetwork_saliva_Control, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1saliva_ControlCohort_RemNetwork.RData")





####### Extract the dir table and then melt it to get relationship of species pairs which will be further used for plotting the network

saliva_ControlCohort_RemNetwork_melt <- Melt_Adjacency_Matrix(RemNetwork_saliva_Control$dir)
saliva_ControlCohort_RemNetwork_melt_filt <- subset(saliva_ControlCohort_RemNetwork_melt, Weight == 1 | Weight == -1)
# 12,403 edges
write.csv(saliva_ControlCohort_RemNetwork_melt_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1saliva_ControlCohort_RemNetwork_melt_filt.csv")

####### Detect the clusters in the network based on edge betweenness using walktrap algorithm
library(igraph)
saliva_ControlCohort_graph <- graph_from_edgelist(as.matrix(saliva_ControlCohort_RemNetwork_melt_filt[,1:2]))
saliva_ControlCohort_walktrap_clusters <- cluster_walktrap(saliva_ControlCohort_graph)
saliva_ControlCohort_walktrap_clusters_df <- data.frame(species = names(membership(saliva_ControlCohort_walktrap_clusters)),
                                                        ClusterID = as.numeric(membership(saliva_ControlCohort_walktrap_clusters)),
                                                        row.names = NULL)

# write.csv(saliva_ControlCohort_walktrap_clusters_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1saliva_ControlCohort_walktrap_clusters_df.csv")







######## There are many edges in the network if we take qvalue as 0.05. so reduce it to 0.001 and rebuild the network
qval_matrix <- data.frame(ifelse(RemNetwork_saliva_Control$qval <= 0.0001, 1, 0))
consistency_matrix <- data.frame(ifelse(RemNetwork_saliva_Control$consistency >= 0.70, 1, 0))
est_matrix <- data.frame(ifelse(RemNetwork_saliva_Control$est > 0, 1, ifelse(RemNetwork_saliva_Control$est < 0, -1, 0)))

# Create dir_matrix based on conditions
dir_matrix  <- (qval_matrix * consistency_matrix * est_matrix)

# now melt it
saliva_dir_matrix <- Melt_Adjacency_Matrix(dir_matrix)
saliva_dir_matrix_filt <- subset(saliva_dir_matrix, Weight == 1)
# edges are not 2420
write.csv(saliva_dir_matrix_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1saliva_ControlCohort_RemNetwork_melt_filt_0.0001.csv")


# Now see the cluster numbers using walktrap algorithm
saliva_dir_graph <- graph_from_edgelist(as.matrix(saliva_dir_matrix_filt[,1:2]))
saliva_dir_walktrap_clusters <- cluster_walktrap(saliva_dir_graph)
saliva_dir_walktrap_clusters_df <- data.frame(species = names(membership(saliva_dir_walktrap_clusters)),
                                                        ClusterID = as.numeric(membership(saliva_dir_walktrap_clusters)),
                                                        row.names = NULL)
# write.csv(saliva_dir_walktrap_clusters_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1saliva_ControlCohort_walktrap_clusters_df_0.0001.csv")








######## compare it with HAC species score and also save the HAC score of the species represented in the network using 0.0001 Qvalue threshold

# # Extract the HAC score for all the species from this network and add them to cytoscape for color gradient
# load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedScores.RData")

# HAC_score_filt <- Combined_Saliva_Scores[rownames(Combined_Saliva_Scores) %in% unique(saliva_dir_walktrap_clusters_df$species), c("HAC_Score","HAC_Score_RankScaled")]
# HAC_score_filt$HAC_group <- cut(HAC_score_filt$HAC_Score_RankScaled,breaks = seq(0, 1, 0.1),labels = paste0("group", 1:10),include.lowest = TRUE)


# write.csv(HAC_score_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1HAC_species_overlapped_ControlCohortNetwork_0.0001.csv")




save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1Network_ControlCohort_Workspace.RData")
