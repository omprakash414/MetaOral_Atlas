#### Import the disease-control samples and then filter the disease samples only.

load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_DisCtrlCohort.RData")

MetadataDf_saliva_Dis <- MetadataDf_saliva_DisCtrl[MetadataDf_saliva_DisCtrl$study_condition == "Diseased",]
SpDf_saliva_Dis <- SpDf_saliva_DisCtrl[rownames(SpDf_saliva_DisCtrl) %in% rownames(MetadataDf_saliva_Dis),]

#### Import the species to check 
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")

rm(MetadataDf_saliva_DisCtrl,SpDf_saliva_DisCtrl,SpDf_nonorm_saliva_DisCtrl,buccal_AssociatedSpecies,subgingival_AssociatedSpecies,supragingival_AssociatedSpecies,tonguetonsil_AssociatedSpecies)

SpDf_saliva_Dis <- SpDf_saliva_Dis[,colnames(SpDf_saliva_Dis)%in% saliva_AssociatedSpecies]

####### Normalize the species profile and then add study_name column to SpDf_saliva_Dis
SpDf_saliva_Dis <- SpDf_saliva_Dis[rowSums(SpDf_saliva_Dis)>0,]
SpDf_saliva_Dis <- SpDf_saliva_Dis/rowSums(SpDf_saliva_Dis)

MetadataDf_saliva_Dis <- MetadataDf_saliva_Dis[rownames(SpDf_saliva_Dis),]
SpDf_saliva_Dis$study_name <- MetadataDf_saliva_Dis$study_name


######### Start Network code on diseased only cohort
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")
RemNetwork_saliva_Disease <- Rem_Network2(SpDf_saliva_Dis,saliva_AssociatedSpecies,"study_name",unique(SpDf_saliva_Dis$study_name),saliva_AssociatedSpecies)


save(RemNetwork_saliva_Disease, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_2saliva_DisCohort_RemNetwork.RData")





####### Extract the dir table and then melt it to get relationship of species pairs which will be further used for plotting the network

saliva_DisCohort_RemNetwork_melt <- Melt_Adjacency_Matrix(RemNetwork_saliva_Disease$dir)
saliva_DisCohort_RemNetwork_melt_filt <- subset(saliva_DisCohort_RemNetwork_melt, Weight == 1 | Weight == -1)

write.csv(saliva_DisCohort_RemNetwork_melt_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_2saliva_DisCohort_RemNetwork_melt_filt.csv")

####### Detect the clusters in the network based on edge betweenness using walktrap algorithm
library(igraph)
saliva_DisCohort_graph <- graph_from_edgelist(as.matrix(saliva_DisCohort_RemNetwork_melt_filt[,1:2]))
saliva_DisCohort_walktrap_clusters <- cluster_walktrap(saliva_DisCohort_graph)
saliva_DisCohort_walktrap_clusters_df <- data.frame(species = names(membership(saliva_DisCohort_walktrap_clusters)),
                                                        ClusterID = as.numeric(membership(saliva_DisCohort_walktrap_clusters)),
                                                        row.names = NULL)

write.csv(saliva_DisCohort_walktrap_clusters_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_2saliva_DisCohort_walktrap_clusters_df.csv")


######## compare it with HAC species score and also save the HAC score of the species represented in the network using 0.0001 Qvalue threshold
CoreSpecies20 <- c("Fusobacterium_periodonticum","Solobacterium_moorei","Eubacterium_sulci","Lachnoanaerobaculum_umeaense","Prevotella_oulorum","Campylobacter_concisus","Actinomyces_graevenitzii","Gemella_sanguinis","Cardiobacterium_hominis","Neisseria_elongata","Veillonella_parvula","Capnocytophaga_gingivalis","Campylobacter_showae","Prevotella_maculosa","Prevotella_melaninogenica","Leptotrichia_goodfellowii","Leptotrichia_hongkongensis","Leptotrichia_hofstadii","Stomatobaculum_longum","Prevotella_salivae")
setdiff(CoreSpecies20,union(saliva_DisCohort_RemNetwork_melt_filt$Node1,saliva_DisCohort_RemNetwork_melt_filt$Node2))



# Extract the HAC score for all the species from this network and add them to cytoscape for color gradient
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedScores.RData")

HAC_score_filt <- Combined_Saliva_Scores[rownames(Combined_Saliva_Scores) %in% unique(saliva_DisCohort_walktrap_clusters_df$species), c("HAC_Score","HAC_Score_RankScaled")]
HAC_score_filt$HAC_group <- cut(HAC_score_filt$HAC_Score_RankScaled,breaks = seq(0, 1, 0.1),labels = paste0("group", 1:10),include.lowest = TRUE)


write.csv(HAC_score_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_2HAC_species_overlapped_DisCohortNetwork.csv")

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_2Network_DiseaseCohort_Workspace.RData")

