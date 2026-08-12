
S9_1CoreAssociation_Validation_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_1CoreAssociation_Validation_Workspace.RData")
attach(S9_1CoreAssociation_Validation_Workspace)
control_SpDf <- control_SpDf
control_metadata <- control_metadata
study_list <- study_list
detach(S9_1CoreAssociation_Validation_Workspace)
rm(S9_1CoreAssociation_Validation_Workspace)

source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")
RemNetwork_ValSaliva_Control <- Rem_Network2(control_SpDf,colnames(control_SpDf),"study_name",study_list,colnames(control_SpDf))


####### Extract the dir table and then melt it to get relationship of species pairs which will be further used for plotting the network (only co-occuring species pairs i.e 1; not -1)

ValSaliva_Network_Melt <- Melt_Adjacency_Matrix(RemNetwork_ValSaliva_Control$dir)
ValSaliva_Network_Melt_filt <- subset(ValSaliva_Network_Melt, Weight == 1)
# 1193 edges
write.csv(ValSaliva_Network_Melt_filt, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_5ValSaliva_Network_Melt_filt.csv")

####### Detect the clusters in the network based on edge betweenness using walktrap algorithm
library(igraph)
ValSaliva_Network_graph <- graph_from_edgelist(as.matrix(ValSaliva_Network_Melt_filt[,1:2]))
# ValSaliva_NetworkWalktrap_clusters <- cluster_walktrap(ValSaliva_Network_graph)
# ValSaliva_NetworkWalktrap_clusters_df <- data.frame(species = names(membership(ValSaliva_NetworkWalktrap_clusters)),
#                                                         ClusterID = as.numeric(membership(ValSaliva_NetworkWalktrap_clusters)),
#                                                         row.names = NULL)

ValSaliva_Network_graph_nodelist <- as.vector(V(ValSaliva_Network_graph)$name)

########### Now take out the properties from the network and then add the Discovery cohorts sHACK, CS, HS, SS from S8_1Saliva_NetworkNodeProperties.txt of discivery cohort
saliva_node_properties_discovery <- read.delim("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_NetworkNodeProperties.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)

# first calculate the properties of validation network. 
saliva_node_properties_Val <- data.frame(row.names=ValSaliva_Network_graph_nodelist,"degree"=degree(ValSaliva_Network_graph)[ValSaliva_Network_graph_nodelist],"betweenness"=betweenness(ValSaliva_Network_graph)[ValSaliva_Network_graph_nodelist])
saliva_node_properties_Val$species <- rownames(saliva_node_properties_Val)

# Now add the discovery cohort properties to the validation cohort properties
saliva_node_properties_discovery <- saliva_node_properties_discovery[saliva_node_properties_discovery$species %in% saliva_node_properties_Val$species,]
saliva_node_properties_Val2 <- saliva_node_properties_Val[saliva_node_properties_Val$species %in% saliva_node_properties_discovery$species,]
saliva_node_properties_Val2 <- saliva_node_properties_Val2[rownames(saliva_node_properties_discovery),]

saliva_node_properties_Val2 <- saliva_node_properties_Val2 %>% left_join(saliva_node_properties_discovery %>% dplyr::select(species, cluster, sHACK, CS, HS, SS),by = "species")
rownames(saliva_node_properties_Val2) <- saliva_node_properties_Val2$species

# now add the withgroup and acrossgroup propensity
saliva_graph_adj_matrix <- as.matrix(as_adjacency_matrix(ValSaliva_Network_graph))
saliva_graph_adj_matrix <- apply(saliva_graph_adj_matrix,2,function(x)(ifelse(x!=1,0,1)))

saliva_graph_adj_matrix <- saliva_graph_adj_matrix[rownames(saliva_node_properties_Val2),rownames(saliva_node_properties_Val2)]

saliva_node_properties_Val2$within_group_propensity <- NA
saliva_node_properties_Val2$across_group_propensity <- NA

for(i in 1:nrow(saliva_node_properties_Val2))
{
  saliva_node_properties_Val2[i,"within_group_propensity"] <- sum(saliva_graph_adj_matrix[rownames(saliva_node_properties_Val2)[i],saliva_node_properties_Val2$cluster==saliva_node_properties_Val2[rownames(saliva_node_properties_Val2)[i],"cluster"]])/(nrow(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==saliva_node_properties_Val2[rownames(saliva_node_properties_Val2)[i],"cluster"],])-1)
  
  saliva_node_properties_Val2[i,"across_group_propensity"] <- sum(saliva_graph_adj_matrix[rownames(saliva_node_properties_Val2)[i],saliva_node_properties_Val2$cluster!=saliva_node_properties_Val2[rownames(saliva_node_properties_Val2)[i],"cluster"]])/nrow(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster!=saliva_node_properties_Val2[rownames(saliva_node_properties_Val2)[i],"cluster"],])
}


write.table(saliva_node_properties_Val2,file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_5Saliva_NetworkNodeProperties_validation.txt",sep="\t", row.names = FALSE, col.names = TRUE)


ValSaliva_Network_Melt_filt2 <- ValSaliva_Network_Melt_filt[ValSaliva_Network_Melt_filt$Node1 %in% saliva_node_properties_Val2$species & ValSaliva_Network_Melt_filt$Node2 %in% saliva_node_properties_Val2$species,]
write.csv(ValSaliva_Network_Melt_filt2, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_5ValSaliva_Network_Melt_filt.csv")



save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_5Validation_NetworkControls.RData")

########## 

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_5ValidationSaliva_Network_WithinAcrossGroupPropensity_Cluster1.pdf",width=2,height=6 )
boxplot(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==1,"within_group_propensity"],saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==1,"across_group_propensity"],col=c("deepskyblue2","darkmagenta"),outline=FALSE,cex.axis=1.5)
dev.off()
wilcox.test(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==1,"within_group_propensity"],saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==1,"across_group_propensity"])
# p-value = 0.2901

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_5ValidationSaliva_Network_WithinAcrossGroupPropensity_Cluster2.pdf",width=2,height=6 )
boxplot(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==2,"within_group_propensity"],saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==2,"across_group_propensity"],col=c("deepskyblue2","darkmagenta"),outline=FALSE,cex.axis=1.5)
dev.off()
wilcox.test(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==2,"within_group_propensity"],saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==2,"across_group_propensity"])
# p-value = 2.134e-08

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_5ValidationSaliva_Network_WithinAcrossGroupPropensity_Cluster3.pdf",width=2,height=6 )
boxplot(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==3,"within_group_propensity"],saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==3,"across_group_propensity"],col=c("deepskyblue2","darkmagenta"),outline=FALSE,cex.axis=1.5)
dev.off()
wilcox.test(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==3,"within_group_propensity"],saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==3,"across_group_propensity"])$p.value
# 7.057759e-22

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S9_Validation/S9_5ValidationSaliva_Network_WithinAcrossGroupPropensity_Cluster4.pdf",width=2,height=6 )
boxplot(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==4,"within_group_propensity"],saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==4,"across_group_propensity"],col=c("deepskyblue2","darkmagenta"),outline=FALSE,cex.axis=1.5)
dev.off() 
wilcox.test(saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==4,"within_group_propensity"],saliva_node_properties_Val2[saliva_node_properties_Val2$cluster==4,"across_group_propensity"])
# p-value = 0.0205
