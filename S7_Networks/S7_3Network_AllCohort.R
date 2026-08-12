
####### Build a network from all the samples. Control as well as diseased.

####### Import the data
# Disease data
S7_2Network_DiseaseCohort_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_2Network_DiseaseCohort_Workspace.RData", envir = S7_2Network_DiseaseCohort_Workspace)
attach(S7_2Network_DiseaseCohort_Workspace)
SpDf_saliva_Dis <- SpDf_saliva_Dis
MetadataDf_saliva_Dis <- MetadataDf_saliva_Dis
detach(S7_2Network_DiseaseCohort_Workspace)
rm(S7_2Network_DiseaseCohort_Workspace)

# Control data
S7_1Network_ControlCohort_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1Network_ControlCohort_Workspace.RData", envir = S7_1Network_ControlCohort_Workspace)
attach(S7_1Network_ControlCohort_Workspace)
SpDf_saliva_Control_associated <- SpDf_saliva_Control_associated
MetadataDf_saliva_control <- MetadataDf_saliva_control
detach(S7_1Network_ControlCohort_Workspace)
rm(S7_1Network_ControlCohort_Workspace)

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")

####### Combine all those two species dfs and metadata dfs
SpDf_saliva_All <- bind_rows(SpDf_saliva_Control_associated,SpDf_saliva_Dis)
SpDf_saliva_All[is.na(SpDf_saliva_All)] <- 0
MetadataDf_saliva_All <- bind_rows(MetadataDf_saliva_control,MetadataDf_saliva_Dis)


####### Normalize the species profile


#### Check the species profile and metadata are in same order
all(rownames(SpDf_saliva_All) == rownames(MetadataDf_saliva_All))
all(nrow(SpDf_saliva_All) == nrow(MetadataDf_saliva_All))

species_to_check <- setdiff(colnames(SpDf_saliva_All),"study_name")

####### Start Network code on control only cohort
RemNetwork_saliva_All <- Rem_Network2(SpDf_saliva_All,species_to_check,"study_name",unique(SpDf_saliva_All$study_name),species_to_check)

save(RemNetwork_saliva_All, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_3saliva_AllCohort_RemNetwork.RData")






####### Extract the dir table and then melt it to get relationship of species pairs which will be further used for plotting the network

saliva_AllCohort_RemNetwork_melt <- Melt_Adjacency_Matrix(RemNetwork_saliva_All$dir)
saliva_AllCohort_RemNetwork_melt_filt <- subset(saliva_AllCohort_RemNetwork_melt, Weight == 1 | Weight == -1)

write.csv(saliva_AllCohort_RemNetwork_melt_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_3saliva_AllCohort_RemNetwork_melt_filt.csv")

####### Detect the clusters in the network based on edge betweenness using walktrap algorithm
library(igraph)
saliva_AllCohort_graph <- graph_from_edgelist(as.matrix(saliva_AllCohort_RemNetwork_melt_filt[,1:2]))
saliva_AllCohort_walktrap_clusters <- cluster_walktrap(saliva_AllCohort_graph)
saliva_AllCohort_walktrap_clusters_df <- data.frame(species = names(membership(saliva_AllCohort_walktrap_clusters)),
                                                        ClusterID = as.numeric(membership(saliva_AllCohort_walktrap_clusters)),
                                                        row.names = NULL)

write.csv(saliva_AllCohort_walktrap_clusters_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_3saliva_AllCohort_walktrap_clusters_df.csv")





######## There are many edges in the network if we take qvalue as 0.05. so reduce it to 0.001 and rebuild the network
qval_matrix <- data.frame(ifelse(RemNetwork_saliva_All$qval <= 0.0001, 1, 0))
consistency_matrix <- data.frame(ifelse(RemNetwork_saliva_All$consistency >= 0.70, 1, 0))
est_matrix <- data.frame(ifelse(RemNetwork_saliva_All$est > 0, 1, ifelse(RemNetwork_saliva_All$est < 0, -1, 0)))

# Create dir_matrix based on conditions
dir_matrix  <- (qval_matrix * consistency_matrix * est_matrix)

# now melt it
saliva_dir_matrix <- Melt_Adjacency_Matrix(dir_matrix)
saliva_dir_matrix_filt <- subset(saliva_dir_matrix, Weight == 1 | Weight == -1)
# 2919 edges
write.csv(saliva_dir_matrix_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_3saliva_AllCohort_RemNetwork_melt_filt_0.0001.csv")


# Now see the cluster numbers using walktrap algorithm
saliva_dir_graph <- graph_from_edgelist(as.matrix(saliva_dir_matrix_filt[,1:2]))
saliva_dir_walktrap_clusters <- cluster_walktrap(saliva_dir_graph)
saliva_dir_walktrap_clusters_df <- data.frame(species = names(membership(saliva_dir_walktrap_clusters)),
                                                        ClusterID = as.numeric(membership(saliva_dir_walktrap_clusters)),
                                                        row.names = NULL)
write.csv(saliva_dir_walktrap_clusters_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_3saliva_AllCohort_walktrap_clusters_df_0.0001.csv")




######## compare it with HAC species score and also save the HAC score of the species represented in the network using 0.0001 Qvalue threshold
CoreSpecies20 <- c("Fusobacterium_periodonticum","Solobacterium_moorei","Eubacterium_sulci","Lachnoanaerobaculum_umeaense","Prevotella_oulorum","Campylobacter_concisus","Actinomyces_graevenitzii","Gemella_sanguinis","Cardiobacterium_hominis","Neisseria_elongata","Veillonella_parvula","Capnocytophaga_gingivalis","Campylobacter_showae","Prevotella_maculosa","Prevotella_melaninogenica","Leptotrichia_goodfellowii","Leptotrichia_hongkongensis","Leptotrichia_hofstadii","Stomatobaculum_longum","Prevotella_salivae")
setdiff(CoreSpecies20,union(saliva_dir_matrix_filt$Node1,saliva_dir_matrix_filt$Node2))


# Extract the HAC score for all the species from this network and add them to cytoscape for color gradient
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedScores.RData")

HAC_score_filt <- Combined_Saliva_Scores[rownames(Combined_Saliva_Scores) %in% unique(saliva_dir_walktrap_clusters_df$species), c("HAC_Score","HAC_Score_RankScaled")]
HAC_score_filt$HAC_group <- cut(HAC_score_filt$HAC_Score_RankScaled,breaks = seq(0, 1, 0.1),labels = paste0("group", 1:10),include.lowest = TRUE)

write.csv(HAC_score_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_3HAC_species_overlapped_AllCohortNetwork_0.0001.csv")



save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_3Network_AllCohort_Workspace.RData")
