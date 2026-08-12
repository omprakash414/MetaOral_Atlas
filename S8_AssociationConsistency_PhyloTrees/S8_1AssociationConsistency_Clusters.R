
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/SS_DataInfo/SS2_AllMetaOral_Data_Scores.RData")


library(psych)
library(cluster)
library(clusterSim)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ade4)
library(vegan)
library(igraph)
library(dunn.test)
library(beanplot)
library(adegraphics)


compute_association_consistency <- function(data,species_list,study_list,total_studies)
{
  sp_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
  rownames(sp_matrix) <- species_list
  colnames(sp_matrix) <- species_list
  
  sn_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
  rownames(sn_matrix) <- species_list
  colnames(sn_matrix) <- species_list
  
  final_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
  rownames(final_matrix) <- species_list
  colnames(final_matrix) <- species_list
  
  for(i in 1:total_studies)
  {
    study <- study_list[i]
    print(study)
    temp_corr <- corr.test(data[data$study_name == study,species_list],method="spearman",use="pairwise.complete",adjust="fdr")
    r_matrix <- apply(temp_corr$r,2,function(x)(ifelse(is.na(x),0,x)))
    p_matrix <- apply(temp_corr$p,2,function(x)(ifelse(is.na(x),1,x)))
    dir_matrix <- apply(r_matrix,2,sign) * apply(p_matrix,2,function(x)(ifelse(x<=0.1,1,0)))
    sp_matrix <- sp_matrix + apply(dir_matrix,2,function(x)(ifelse(x==1,1,0)))
    sn_matrix <- sn_matrix + apply(dir_matrix,2,function(x)(ifelse(x== -1,1,0)))
  }
  
  final_matrix <- ((sp_matrix - sn_matrix)/total_studies)*(1-(pmin(sp_matrix,sn_matrix)+1)/(pmax(sp_matrix,sn_matrix)+1))
  
  return_list <- list("final_matrix"=final_matrix,"sp_matrix"=sp_matrix,"sn_matrix"=sn_matrix)
}



enterotyping_clustering = function(data,cluster_size)
{
  data.dist = vegdist(data,method="euclidean")
  #data.dist = as.dist(1-cor(t(data),method="kendall")/2)
  nclusters=NULL
  mean.obs.silhouette=NULL
  
  for (k in 1:cluster_size) {
    
    if (k==1) {
      nclusters[k]=NA 
      mean.obs.silhouette=NA
    } else {
      data.cluster_temp=kmeans(data.dist,k)
      nclusters[k]=index.G1(data,data.cluster_temp$cluster,  d = data.dist, centrotypes = "medoids")
      mean.obs.silhouette[k] <- mean(silhouette(data.cluster_temp$cluster, data.dist)[,3])
    }
  }
  
  returnList = list("nclusters"=nclusters,"silhouette"=mean.obs.silhouette)
  return(returnList);
}

iterative_enterotyping <- function(data,ncluster,iter)
{
  mat_ch_index <- as.data.frame(matrix(NA,iter,ncluster))
  mat_silhouette <- as.data.frame(matrix(NA,iter,ncluster))
  for(i in 1:iter)
  {
    t <- enterotyping_clustering(data,ncluster)
    print(dim(mat_ch_index))
    print(length(t$nclusters))
    mat_ch_index[i,] <- t$nclusters[1:ncluster]
    mat_silhouette[i,] <- t$silhouette[1:ncluster]
  }
  returnList <- list("ch_index"=mat_ch_index,"silhouette"=mat_silhouette)
  return(returnList)
}


print("Identifying ecologically similar taxa modules in the salivary microbiomes")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######         Saliva subsite association consistency and clustering
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

saliva_controls <- rownames(MetadataDf_All_Discovery_filt)[(MetadataDf_All_Discovery_filt$body_site_category == "saliva_sputum_oral_wash")&(MetadataDf_All_Discovery_filt$study_condition == "Control")]
SpDf_saliva_Control_associated <- SpDf_All_Discovery_filt[saliva_controls,rownames(Combined_Saliva_Scores)]
study_list <- unique(MetadataDf_All_Discovery_filt[saliva_controls,"study_name"])
total_studies <- length(study_list)
species_list <- rownames(Combined_Saliva_Scores)

data <- SpDf_saliva_Control_associated
data$study_name <- MetadataDf_All_Discovery_filt[rownames(data),"study_name"]

saliva_association_matrix <- compute_association_consistency(data,species_list,study_list,total_studies)
pca_saliva_association <- dudi.pca(saliva_association_matrix$final_matrix,scannf=FALSE,nf=3)
df_pca_saliva <- data.frame("Axis1"=pca_saliva_association$li[,1],"Axis2"=pca_saliva_association$li[,2],"Axis3"=pca_saliva_association$li[,3],"CS"=Combined_Saliva_Scores[rownames(pca_saliva_association$li),1],"HS"=Combined_Saliva_Scores[rownames(pca_saliva_association$li),2],"SS"=Combined_Saliva_Scores[rownames(pca_saliva_association$li),3],"sHACK"=Combined_Saliva_Scores[rownames(pca_saliva_association$li),5],row.names=rownames(pca_saliva_association$li))

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_PCA_sHACKClustering.pdf",width=6,height=6)
ggplot(df_pca_saliva,aes(x=Axis1,y=Axis2,colour=as.character(ifelse(sHACK>=0.90,3,ifelse(sHACK>=0.50,2,1)))))+geom_point(size=3)+scale_color_manual(values=c("3"="cornflowerblue","2"="lightblue","1"="red"))+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),legend.position="none")
dev.off()


## Find the optimal number of clusters using the CH index and silhouette index
set.seed(10);optimal_k_pca_saliva <- iterative_enterotyping(df_pca_saliva[,c(1,2)],20,10)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_CHIndex.pdf",width=6,height=5)
boxplot(optimal_k_pca_saliva$ch_index[,2],optimal_k_pca_saliva$ch_index[,3],optimal_k_pca_saliva$ch_index[,4],optimal_k_pca_saliva$ch_index[,5],optimal_k_pca_saliva$ch_index[,6],optimal_k_pca_saliva$ch_index[,7],optimal_k_pca_saliva$ch_index[,8],optimal_k_pca_saliva$ch_index[,9],optimal_k_pca_saliva$ch_index[,10],names=c("2","3","4","5","6","7","8","9","10"),cex.axis=2,outline=FALSE,ylab="CH-Index",cex.title=2)
dev.off()

print("Five optimal clusters identified.")
print("Lets have a look at these clusters.")

set.seed(10); saliva_kmeans <- kmeans(df_pca_saliva[,c(1,2)],5); df_pca_saliva$cluster <- saliva_kmeans$cluster; s.class(df_pca_saliva[,c(1,2)],as.factor(df_pca_saliva$cluster),col="red")
dev.off()

all(rownames(df_pca_saliva) == rownames(Combined_Saliva_Scores))
df_pca_saliva$sHACK <- Combined_Saliva_Scores[rownames(df_pca_saliva),5]

df_pca_saliva$cluster <- saliva_kmeans$cluster
cx <- saliva_kmeans$centers[,1]
cy <- saliva_kmeans$centers[,2]

df_pca_saliva$cx <- cx[df_pca_saliva$cluster]
df_pca_saliva$cy <- cy[df_pca_saliva$cluster]

save(df_pca_saliva, file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_PCAdf.RData")
write.csv(df_pca_saliva, file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_PCAdf.csv",row.names=TRUE,quote=FALSE)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_FreePlot.pdf",width=6,height=6)
ggplot(df_pca_saliva,aes(x=Axis1,y=Axis2))+geom_point(color="darkolivegreen",size=3)+theme_classic() + theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),legend.position="none")
dev.off()
# ggplot(df_pca_saliva,aes(x=Axis1,y=Axis2,color=as.character(ifelse(sHACK>=0.90,3,ifelse(sHACK>=0.75,2,1))))) + geom_point(size=3) + geom_segment(aes(x = cx, y = cy, xend = Axis1, yend = Axis2), color = "grey30", linewidth = 0.3, alpha = 0.8) + stat_ellipse(aes(group = cluster),color = "grey40", linewidth = 0.5, level = 0.95) + scale_color_manual(values=c("3"="cornflowerblue","2"="cyan","1"="red")) + theme_classic() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),legend.position="none") + geom_text_repel(label=ifelse(df_pca_saliva$sHACK>=0.90,rownames(df_pca_saliva),""),box.padding=0.000001,size=3.5,max.overlaps=40,col="black")

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_Clusters_sHACK_distribution.pdf",width=8,height=5)
ggplot(df_pca_saliva,aes(x=Axis1,y=Axis2,color=as.character(ifelse(sHACK>=0.90,3,ifelse(sHACK>=0.75,2,1))))) + geom_point(size=3) + geom_segment(aes(x = cx, y = cy, xend = Axis1, yend = Axis2), color = "grey30", linewidth = 0.3, alpha = 0.8) + stat_ellipse(aes(group = cluster),color = "grey40", linewidth = 0.5, level = 0.95) + scale_color_manual(values=c("3"="cornflowerblue","2"="cyan","1"="red")) + theme_classic() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),legend.position="none") + geom_text_repel(label=ifelse(df_pca_saliva$sHACK>=0.90,rownames(df_pca_saliva),""),box.padding=0.000001,size=3.5,max.overlaps=40,col="black")
dev.off()

print("Each cluster variably associates with sHACK, CS, HS, Scores")

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_Cluster_CS_Scores.pdf",width=6,height=6)
boxplot(df_pca_saliva$CS~df_pca_saliva$cluster,col=c("brown","cornflowerblue","darkolivegreen","darkgoldenrod","firebrick4"),outline=FALSE,cex.axis=2,xlab="",ylab="")
dev.off()
dunn_cluster_CS_saliva <- dunn.test(df_pca_saliva$CS,df_pca_saliva$cluster)
#   Kruskal-Wallis rank sum test

# data: x and group
# Kruskal-Wallis chi-squared = 159.246, df = 4, p-value = 0


#                     Dunn's Pairwise Comparison of x by group                    
#                                  (No adjustment)                                

# Col Mean-│
# Row Mean │          1          2          3          4
# ─────────┼────────────────────────────────────────────
#        2 │  -4.974809
#          │     0.0000*
#          │
#        3 │  -4.680428   0.469065
#          │     0.0000*    0.3195 
#          │
#        4 │   1.595996   9.054131   8.832575
#          │     0.0552     0.0000*    0.0000*
#          │
#        5 │   2.314031   8.860598   8.626602   1.215311
#          │     0.0103*    0.0000*    0.0000*    0.1121 

# α = 0.05
# Reject Ho if p ≤ α/2, where p = Pr(Z ≥ |z|)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_Cluster_HS_Scores.pdf",width=6,height=6)
boxplot(df_pca_saliva$HS~df_pca_saliva$cluster,col=c("brown","cornflowerblue","darkolivegreen","darkgoldenrod","firebrick4"),outline=FALSE,cex.axis=2,xlab="",ylab="")
dev.off()
dunn_cluster_HS_saliva <- dunn.test(df_pca_saliva$HS,df_pca_saliva$cluster)
#   Kruskal-Wallis rank sum test

# data: x and group
# Kruskal-Wallis chi-squared = 97.0871, df = 4, p-value = 0


#                     Dunn's Pairwise Comparison of x by group                    
#                                  (No adjustment)                                

# Col Mean-│
# Row Mean │          1          2          3          4
# ─────────┼────────────────────────────────────────────
#        2 │  -7.731420
#          │     0.0000*
#          │
#        3 │  -3.290105   5.613766
#          │     0.0005*    0.0000*
#          │
#        4 │  -4.324647   5.351730  -1.022127
#          │     0.0000*    0.0000*    0.1534 
#          │
#        5 │  -0.632668   8.629518   3.244691   4.649891
#          │     0.2635     0.0000*    0.0006*    0.0000*

# α = 0.05
# Reject Ho if p ≤ α/2, where p = Pr(Z ≥ |z|)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_Cluster_SS_Scores.pdf",width=6,height=6)
boxplot(df_pca_saliva$SS~df_pca_saliva$cluster,col=c("brown","cornflowerblue","darkolivegreen","darkgoldenrod","firebrick4"),outline=FALSE,cex.axis=2,xlab="",ylab="")
dev.off()
dunn_cluster_SS_saliva <- dunn.test(df_pca_saliva$SS,df_pca_saliva$cluster)
#   Kruskal-Wallis rank sum test

# data: x and group
# Kruskal-Wallis chi-squared = 85.2794, df = 4, p-value = 0


#                     Dunn's Pairwise Comparison of x by group                    
#                                  (No adjustment)                                

# Col Mean-│
# Row Mean │          1          2          3          4
# ─────────┼────────────────────────────────────────────
#        2 │  -4.631960
#          │     0.0000*
#          │
#        3 │  -2.579414   2.617402
#          │     0.0049*    0.0044*
#          │
#        4 │  -0.921218   5.345713   2.476925
#          │     0.1785     0.0000*    0.0066*
#          │
#        5 │   2.554199   8.735775   6.350130   4.779059
#          │     0.0053*    0.0000*    0.0000*    0.0000*

# α = 0.05
# Reject Ho if p ≤ α/2, where p = Pr(Z ≥ |z|)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_Cluster_sHACK_Scores.pdf",width=6,height=6)
boxplot(df_pca_saliva$sHACK~df_pca_saliva$cluster,col=c("brown","cornflowerblue","darkolivegreen","darkgoldenrod","firebrick4"),outline=FALSE,cex.axis=2,xlab="",ylab="")
dev.off()
dunn_cluster_sHACK_saliva <- dunn.test(df_pca_saliva$sHACK,df_pca_saliva$cluster)
#   Kruskal-Wallis rank sum test

# data: x and group
# Kruskal-Wallis chi-squared = 155.9252, df = 4, p-value = 0


#                     Dunn's Pairwise Comparison of x by group                    
#                                  (No adjustment)                                

# Col Mean-│
# Row Mean │          1          2          3          4
# ─────────┼────────────────────────────────────────────
#        2 │  -7.399533
#          │     0.0000*
#          │
#        3 │  -5.544407   2.435477
#          │     0.0000*    0.0074*
#          │
#        4 │  -1.594852   8.381878   5.835629
#          │     0.0554     0.0000*    0.0000*
#          │
#        5 │   1.216897   10.47446   8.316876   3.758041
#          │     0.1118     0.0000*    0.0000*    0.0001*

# α = 0.05
# Reject Ho if p ≤ α/2, where p = Pr(Z ≥ |z|)


tapply(ifelse(df_pca_saliva$sHACK>=0.90,3,ifelse(df_pca_saliva$sHACK>=0.75,2,1)),list(df_pca_saliva$cluster),table)


cluster_counts <- table(factor(df_pca_saliva[df_pca_saliva$sHACK >= 0.90, "cluster"],levels = 1:5))
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_sHACK90_Clusters_Barplot.pdf",width=4,height=6)
barplot(cluster_counts,col = "#156082",border = "black",xlab = "Cluster",ylab = "Number of species",ylim = c(0, max(cluster_counts) + 5),cex.names = 1.2,cex.axis = 1.1,cex.lab = 1.3,cex.main = 1.2)
dev.off()
## Fishers' exact test.
#Make cluster a factor so cluster 5 is included even if zero
df_pca_saliva$cluster_factor <- factor(df_pca_saliva$cluster, levels = 1:5)
## Define high sHACK
df_pca_saliva$High_sHACK <- df_pca_saliva$sHACK >= 0.90
## Define Module 2 vs other modules
df_pca_saliva$Module2 <- df_pca_saliva$cluster == 2
## 2 x 2 contingency table
saliva_fisher_table <- table(High_sHACK = df_pca_saliva$High_sHACK,Module2 = df_pca_saliva$Module2)
saliva_fisher <- fisher.test(saliva_fisher_table)
saliva_fisher$p.value # [1] 6.364959e-11
saliva_fisher$estimate # odds ratio: 8.196911





print("Lets re-examine the community using Meta-Networks")
saliva_network <- read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S7_Networks/S7_1saliva_ControlCohort_RemNetwork_melt_filt_0.0001.csv")
saliva_network <- saliva_network[saliva_network[,3] == 1,]
saliva_graph <- graph_from_edgelist(as.matrix(saliva_network[,c(1,2)]))



table(data.frame("cluster_node1"=df_pca_saliva[saliva_network[,1],"cluster"],"cluster_node2"=df_pca_saliva[saliva_network[,2],"cluster"]))

saliva_graph_nodelist <- as.vector(V(saliva_graph)$name)
saliva_node_properties <- data.frame(row.names=saliva_graph_nodelist,"cluster"=df_pca_saliva[saliva_graph_nodelist,"cluster"],"degree"=degree(saliva_graph)[saliva_graph_nodelist],"betweenness"=betweenness(saliva_graph)[saliva_graph_nodelist],"sHACK"=df_pca_saliva[saliva_graph_nodelist,"sHACK"],"CS"=Combined_Saliva_Scores[saliva_graph_nodelist,1],"HS"=Combined_Saliva_Scores[saliva_graph_nodelist,2],"SS"=Combined_Saliva_Scores[saliva_graph_nodelist,3])
saliva_node_properties2 <- saliva_node_properties
saliva_node_properties2$species <- rownames(saliva_node_properties2)
write.table(saliva_node_properties2,file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_NetworkNodeProperties.txt",sep="\t")

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_Network_DegreeDistribution_Clusters.pdf",width=6,height=6)
boxplot(saliva_node_properties$degree~saliva_node_properties$cluster)
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_Network_BetweennessDistribution_Clusters.pdf",width=6,height=6)
boxplot(saliva_node_properties$betweenness~saliva_node_properties$cluster)
dev.off()

saliva_graph_adj_matrix <- as.matrix(as_adjacency_matrix(saliva_graph))
saliva_graph_adj_matrix <- apply(saliva_graph_adj_matrix,2,function(x)(ifelse(x!=1,0,1)))

saliva_node_properties$within_group_propensity <- NA
saliva_node_properties$across_group_propensity <- NA

for(i in 1:nrow(saliva_node_properties))
{
  saliva_node_properties[i,"within_group_propensity"] <- sum(saliva_graph_adj_matrix[rownames(saliva_node_properties)[i],saliva_node_properties$cluster==saliva_node_properties[rownames(saliva_node_properties)[i],"cluster"]])/(nrow(saliva_node_properties[saliva_node_properties$cluster==saliva_node_properties[rownames(saliva_node_properties)[i],"cluster"],])-1)
  
  saliva_node_properties[i,"across_group_propensity"] <- sum(saliva_graph_adj_matrix[rownames(saliva_node_properties)[i],saliva_node_properties$cluster!=saliva_node_properties[rownames(saliva_node_properties)[i],"cluster"]])/nrow(saliva_node_properties[saliva_node_properties$cluster!=saliva_node_properties[rownames(saliva_node_properties)[i],"cluster"],])
}

print("Taxa within modules have signfiicantly more intra-module connections")
library(beanplot)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_Network_WithinAcrossGroupPropensity.pdf",width=6,height=6)
beanplot(saliva_node_properties[saliva_node_properties$cluster==1,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==1,"across_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==2,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==2,"across_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==3,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==3,"across_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==4,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==4,"across_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==5,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==5,"across_group_propensity"],saliva_node_properties[,"within_group_propensity"],saliva_node_properties[,"across_group_propensity"],side="both",what=c(1,1,1,0),col=list("deepskyblue2","darkmagenta"),outline=FALSE)
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_Network_WithinAcrossGroupPropensity_Cluster1.pdf",width=2,height=6 )
boxplot(saliva_node_properties[saliva_node_properties$cluster==1,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==1,"across_group_propensity"],col=c("deepskyblue2","darkmagenta"),outline=FALSE,cex.axis=1.5)
dev.off()
wilcox.test(saliva_node_properties[saliva_node_properties$cluster==1,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==1,"across_group_propensity"])
#         Wilcoxon rank sum test with continuity correction

# data:  saliva_node_properties[saliva_node_properties$cluster == 1, "within_group_propensity"] and saliva_node_properties[saliva_node_properties$cluster == 1, "across_group_propensity"]
# W = 1126, p-value = 0.007723
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(saliva_node_properties[saliva_node_properties$cluster ==  :
#   cannot compute exact p-value with ties

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_Network_WithinAcrossGroupPropensity_Cluster2.pdf",width=2,height=6 )
boxplot(saliva_node_properties[saliva_node_properties$cluster==2,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==2,"across_group_propensity"],col=c("deepskyblue2","darkmagenta"),outline=FALSE,cex.axis=1.5)
dev.off()
wilcox.test(saliva_node_properties[saliva_node_properties$cluster==2,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==2,"across_group_propensity"])
#         Wilcoxon rank sum test with continuity correction

# data:  saliva_node_properties[saliva_node_properties$cluster == 2, "within_group_propensity"] and saliva_node_properties[saliva_node_properties$cluster == 2, "across_group_propensity"]
# W = 4844, p-value = 9.576e-08
# alternative hypothesis: true location shift is not equal to 0

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_Network_WithinAcrossGroupPropensity_Cluster3.pdf",width=2,height=6 )
boxplot(saliva_node_properties[saliva_node_properties$cluster==3,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==3,"across_group_propensity"],col=c("deepskyblue2","darkmagenta"),outline=FALSE,cex.axis=1.5)
dev.off()
wilcox.test(saliva_node_properties[saliva_node_properties$cluster==3,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==3,"across_group_propensity"])$p.value
# [1] 1.229132e-20

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_Network_WithinAcrossGroupPropensity_Cluster4.pdf",width=2,height=6 )
boxplot(saliva_node_properties[saliva_node_properties$cluster==4,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==4,"across_group_propensity"],col=c("deepskyblue2","darkmagenta"),outline=FALSE,cex.axis=1.5)
dev.off() 
wilcox.test(saliva_node_properties[saliva_node_properties$cluster==4,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==4,"across_group_propensity"])
#         Wilcoxon rank sum test with continuity correction

# data:  saliva_node_properties[saliva_node_properties$cluster == 4, "within_group_propensity"] and saliva_node_properties[saliva_node_properties$cluster == 4, "across_group_propensity"]
# W = 1492, p-value = 0.008576
# alternative hypothesis: true location shift is not equal to 0


# wilcox.test(saliva_node_properties[saliva_node_properties$cluster==5,"within_group_propensity"],saliva_node_properties[saliva_node_properties$cluster==5,"across_group_propensity"])

# wilcox.test(saliva_node_properties[,"within_group_propensity"],saliva_node_properties[,"across_group_propensity"])

print("High sHACK microbes especially enriched in Cluster 2")
cluster1 <- rownames(df_pca_saliva[df_pca_saliva$cluster == 1,])
cluster2 <- rownames(df_pca_saliva[df_pca_saliva$cluster == 2,])
cluster3 <- rownames(df_pca_saliva[df_pca_saliva$cluster == 3,])
cluster4 <- rownames(df_pca_saliva[df_pca_saliva$cluster == 4,])
cluster5 <- rownames(df_pca_saliva[df_pca_saliva$cluster == 5,])

high_sHACK_cluster2 <- rownames(df_pca_saliva[(df_pca_saliva$cluster==2)&(df_pca_saliva$sHACK >= 0.90),])

high_sHACK_cluster3 <- rownames(df_pca_saliva[(df_pca_saliva$cluster==3)&(df_pca_saliva$sHACK >= 0.90),])






save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/260625_MetaOral_TSG_OS.RData")






# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######         Validation of the clusters in the independent validation datasets for Saliva
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Similar arrangement of taxa in the Validation datasets")
validation_controls <- rownames(validation_salivary_microbiome_scores)[validation_salivary_microbiome_scores$Study_Condition == "Control"]

data_validation <- SpDf_saliva_Validation[validation_controls,]
data_validation$study_name <- MetadataDf_saliva_validation[rownames(data_validation),"study_name"]

study_list_validation <- unique(data_validation$study_name)

total_studies_validation <- length(study_list_validation)

saliva_association_matrix_validation <- compute_association_consistency(data_validation,species_list,study_list_validation,total_studies_validation)

pca_saliva_association_validation <- dudi.pca(saliva_association_matrix_validation,scannf=FALSE,nf=3)
df_pca_saliva_validation <- data.frame("Axis1"=pca_saliva_association_validation$li[,1],"Axis2"=pca_saliva_association_validation$li[,2],"Axis3"=pca_saliva_association_validation$li[,3],row.names=rownames(pca_saliva_association_validation$li))
df_pca_saliva_validation$CS <- Combined_Saliva_Scores[rownames(df_pca_saliva_validation),"CS"]
df_pca_saliva_validation$HS <- Combined_Saliva_Scores[rownames(df_pca_saliva_validation),"HS"]
df_pca_saliva_validation$SS <- Combined_Saliva_Scores[rownames(df_pca_saliva_validation),"SS"]
df_pca_saliva_validation$sHACK <- Combined_Saliva_Scores[rownames(df_pca_saliva_validation),5]

df_pca_saliva_validation[,2] <- (-1) * df_pca_saliva_validation[,2]

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_PCA_Validation.pdf",width=6,height=6)
ggplot(df_pca_saliva_validation,aes(x=Axis1,y=Axis2,colour=as.character(ifelse(sHACK>=0.90,3,ifelse(sHACK>=0.50,2,1)))))+geom_point(size=3)+scale_color_manual(values=c("3"="cornflowerblue","2"="lightblue","1"="red"))+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.position="none")
dev.off()
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_AssociationConsistency_PCA_Validation_Clusters.pdf",width=6,height=6)
s.class(df_pca_saliva_validation[,c(1,2)],as.factor(df_pca_saliva$cluster),col=c("brown","cornflowerblue","darkolivegreen","darkgoldenrod","firebrick4"))
dev.off()

s.class(df_pca_saliva[,c(1,2)],as.factor(df_pca_saliva$cluster),col=c("brown","cornflowerblue","darkolivegreen","darkgoldenrod","firebrick4"))


print("Similar structure in validation as well")
procuste.randtest(df_pca_saliva[,c(1,2)],df_pca_saliva_validation[,c(1,2)])

print("High Association of cluster2 with health across diseases in validation cohorts is s reproduced")

#Remove Rows
ChenC_2018_remove_rows <- rownames(MetadataDf_saliva_validation[(MetadataDf_saliva_validation$study_name == "ChenC_2018")&(MetadataDf_saliva_validation$study_condition == "Post-Disease"),])

ChenJ_2021_remove_rows <- rownames(MetadataDf_saliva_validation[MetadataDf_saliva_validation$original_sample_id %in% grep("^NRVH",MetadataDf_saliva_validation[(MetadataDf_saliva_validation$study_name == "ChenJ_2021")&(MetadataDf_saliva_validation$study_condition!="Control"),"original_sample_id"],value=TRUE),])

IglesiasA_2024_remove_rows <- rownames(MetadataDf_saliva_validation[(MetadataDf_saliva_validation$study_name == "IglesiasA_2024")&(MetadataDf_saliva_validation$study_condition == "Control")&(MetadataDf_saliva_validation$exposure == "smoker"),])

SpDf_saliva_Validation <- SpDf_saliva_Validation[setdiff(rownames(SpDf_saliva_Validation),c(ChenC_2018_remove_rows,ChenJ_2021_remove_rows,IglesiasA_2024_remove_rows)),]

validation_salivary_microbiome_scores <- data.frame("sHACK"=rowMeans(apply(SpDf_saliva_Validation[,rownames(df_validation_directions[rownames(df_discovery_directions[1:20,]),])],2,rank_scale)))

validation_salivary_microbiome_scores$Study_Condition <- ifelse(MetadataDf_saliva_validation[rownames(validation_salivary_microbiome_scores),"study_condition"] == "Control","Control","Diseased")

validation_salivary_microbiome_scores$Study_Name <- MetadataDf_saliva_validation[rownames(validation_salivary_microbiome_scores),"study_name"]

validation_controls <- rownames(validation_salivary_microbiome_scores)[validation_salivary_microbiome_scores$Study_Condition == "Control"]

validation_diseased <- rownames(validation_salivary_microbiome_scores)[validation_salivary_microbiome_scores$Study_Condition != "Control"]

control_studies <- unique(validation_salivary_microbiome_scores[validation_controls,"Study_Name"])

diseased_studies <- unique(validation_salivary_microbiome_scores[validation_diseased,"Study_Name"])

matched_studies <- intersect(control_studies,diseased_studies)

df_comparison_all_taxa <- as.data.frame(matrix(NA,nrow(df_pca_saliva),length(matched_studies)))
rownames(df_comparison_all_taxa) <- rownames(df_pca_saliva)
colnames(df_comparison_all_taxa) <- matched_studies

for(i in 1:length(matched_studies))
{
  study <- matched_studies[i]
  all_study_samples <- rownames(validation_salivary_microbiome_scores)[validation_salivary_microbiome_scores$Study_Name == study]
  study_controls <- intersect(validation_controls,all_study_samples)
  study_diseased <- intersect(validation_diseased,all_study_samples)
  
  for(j in 1:nrow(df_pca_saliva))
  {
    species_name <- rownames(df_pca_saliva)[j]
    temp_cohen_d <- cohen.d(as.numeric(SpDf_saliva_Validation[study_controls,species_name]),as.numeric(SpDf_saliva_Validation[study_diseased,species_name]))$estimate
    
    temp_wilcox <- wilcox.test(as.numeric(SpDf_saliva_Validation[study_controls,species_name]),as.numeric(SpDf_saliva_Validation[study_diseased,species_name]))$p.value
    
    df_comparison_all_taxa[j,i] <- sign(temp_cohen_d) * ifelse(temp_wilcox<=0.05,3,ifelse(temp_wilcox <= 0.1,2,1))
  }
}

df_validation_directions <- data.frame("control_increased"=apply(df_comparison_all_taxa,1,function(x)(length(x[!is.na(x)&(x>2)]))),"control_depleted"=apply(df_comparison_all_taxa,1,function(x)(length(x[!is.na(x)&(x<(-2))]))))

df_validation_directions$cluster <- df_pca_saliva[rownames(df_validation_directions),"cluster"]

df_validation_directions$score <- apply(df_validation_directions,1,function(x)(((x[1]-x[2])/(15))*(1-(min(x[1],x[2])+0.00001)/(max(x[1],x[2])+0.00001))))

df_validation_directions$HS <- Combined_Saliva_Scores[rownames(df_validation_directions),2]

df_validation_directions$CS <- Combined_Saliva_Scores[rownames(df_validation_directions),1]

df_validation_directions$SS <- Combined_Saliva_Scores[rownames(df_validation_directions),3]

df_validation_directions$sHACK <- Combined_Saliva_Scores[rownames(df_validation_directions),4]


boxplot(df_validation_directions$control_increased~df_validation_directions$cluster,col=c("brown","cornflowerblue","darkolivegreen","darkgoldenrod","firebrick4"),outline=FALSE,cex.axis=2,ylab=c("","","","","",""))

dunns_validation_dataset_cluster_links <- dunn.test(df_validation_directions$control_increased,df_validation_directions$cluster)

procuste.randtest(df_pca_saliva[,c(1,2)],df_pca_saliva_validation[,c(1,2)])

validation_salivary_microbiome_scores <- data.frame("sHACK"=rowMeans(apply(SpDf_saliva_Validation[,rownames(df_validation_directions[rownames(df_discovery_directions[1:20,]),])],2,rank_scale)))

validation_salivary_microbiome_scores$Study_Condition <- ifelse(MetadataDf_saliva_validation[rownames(validation_salivary_microbiome_scores),"study_condition"] == "Control","Control","Diseased")

validation_salivary_microbiome_scores$Study_Name <- MetadataDf_saliva_validation[rownames(validation_salivary_microbiome_scores),"study_name"]

validation_salivary_microbiome_scores$Study_Condition_Binary <- ifelse(validation_salivary_microbiome_scores$Study_Condition == "Control",0,1)

sHACK <- rownames(Combined_Saliva_Scores[Combined_Saliva_Scores$HACK_Score >= 0.90,])

validation_salivary_microbiome_scores$sHACK <- rowMeans(apply(SpDf_saliva_Validation[rownames(validation_salivary_microbiome_scores),sHACK],2,rank_scale))

high_sHACK_cluster2 <- rownames(df_pca_saliva[(df_pca_saliva$sHACK>=0.90)&(df_pca_saliva$cluster == 2),])
high_sHACK_cluster3 <- rownames(df_pca_saliva[(df_pca_saliva$sHACK>=0.90)&(df_pca_saliva$cluster == 3),])

validation_salivary_microbiome_scores$high_sHACK_cluster2 <- rowMeans(apply(SpDf_saliva_Validation[rownames(validation_salivary_microbiome_scores),high_sHACK_cluster2],2,rank_scale))

validation_salivary_microbiome_scores$high_sHACK_cluster3 <- rowMeans(apply(SpDf_saliva_Validation[rownames(validation_salivary_microbiome_scores),high_sHACK_cluster3],2,rank_scale))

validation_salivary_microbiome_scores$dysbiosis_bray <- NA

validation_salivary_microbiome_scores$dysbiosis_canberra <- NA

df_comparison_sHACK_score <- as.data.frame(matrix(NA,length(matched_studies),10))
rownames(df_comparison_sHACK_score) <- matched_studies
colnames(df_comparison_sHACK_score) <- c("cohen_d_sHACK","p_val_sHACK","cohen_d_ds_bray","p_val_ds_bray","cohen_d_ds_canberra","p_val_ds_canberra","cohen_d_cluster2","p_val_cluster2","cohen_d_cluster3","p_val_cluster3")

for(i in 1:length(matched_studies))
{
  study <- matched_studies[i]
  all_study_samples <- rownames(validation_salivary_microbiome_scores)[validation_salivary_microbiome_scores$Study_Name == study]
  study_controls <- intersect(validation_controls,all_study_samples)
  study_diseased <- intersect(validation_diseased,all_study_samples)
  
  df_comparison_sHACK_score[study,"cohen_d_sHACK"] <- cohen.d(as.numeric(validation_salivary_microbiome_scores[study_controls,1]),as.numeric(validation_salivary_microbiome_scores[study_diseased,1]))$estimate
  
  df_comparison_sHACK_score[study,"p_val_sHACK"] <- wilcox.test(as.numeric(validation_salivary_microbiome_scores[study_controls,1]),as.numeric(validation_salivary_microbiome_scores[study_diseased,1]))$p.value
  
  dist_mat <- as.matrix(vegdist(SpDf_saliva_Validation[all_study_samples,],method="bray"))
  diag(dist_mat) <- NA
  df_comparison_sHACK_score[study,"cohen_d_ds_bray"] <- cohen.d(apply(dist_mat[study_controls,study_controls],1,function(x)(mean(x[!is.na(x)]))),apply(dist_mat[study_controls,study_diseased],1,function(x)(mean(x[!is.na(x)]))))$estimate
  df_comparison_sHACK_score[study,"p_val_ds_bray"] <- wilcox.test(apply(dist_mat[study_controls,study_controls],1,function(x)(mean(x[!is.na(x)]))),apply(dist_mat[study_controls,study_diseased],1,function(x)(mean(x[!is.na(x)]))))$p.value
  
  validation_salivary_microbiome_scores[all_study_samples,"dysbiosis_bray"] <- apply(dist_mat[all_study_samples,study_controls],1,function(x)(mean(x[!is.na(x)])))
  
  dist_mat <- as.matrix(vegdist(SpDf_saliva_Validation[all_study_samples,],method="canberra"))
  diag(dist_mat) <- NA
  df_comparison_sHACK_score[study,"cohen_d_ds_canberra"] <- cohen.d(apply(dist_mat[study_controls,study_controls],1,function(x)(mean(x[!is.na(x)]))),apply(dist_mat[study_controls,study_diseased],1,function(x)(mean(x[!is.na(x)]))))$estimate
  df_comparison_sHACK_score[study,"p_val_ds_canberra"] <- wilcox.test(apply(dist_mat[study_controls,study_controls],1,function(x)(mean(x[!is.na(x)]))),apply(dist_mat[study_controls,study_diseased],1,function(x)(mean(x[!is.na(x)]))))$p.value
  
  validation_salivary_microbiome_scores[all_study_samples,"dysbiosis_canberra"] <- apply(dist_mat[all_study_samples,study_controls],1,function(x)(mean(x[!is.na(x)])))
  
  df_comparison_sHACK_score[study,"cohen_d_cluster2"] <- cohen.d(as.numeric(validation_salivary_microbiome_scores[study_controls,"high_sHACK_cluster2"]),as.numeric(validation_salivary_microbiome_scores[study_diseased,"high_sHACK_cluster2"]))$estimate
  
  df_comparison_sHACK_score[study,"p_val_cluster2"] <- wilcox.test(as.numeric(validation_salivary_microbiome_scores[study_controls,"high_sHACK_cluster2"]),as.numeric(validation_salivary_microbiome_scores[study_diseased,"high_sHACK_cluster2"]))$p.value
  
  df_comparison_sHACK_score[study,"cohen_d_cluster3"] <- cohen.d(as.numeric(validation_salivary_microbiome_scores[study_controls,"high_sHACK_cluster3"]),as.numeric(validation_salivary_microbiome_scores[study_diseased,"high_sHACK_cluster3"]))$estimate
  
  df_comparison_sHACK_score[study,"p_val_cluster3"] <- wilcox.test(as.numeric(validation_salivary_microbiome_scores[study_controls,"high_sHACK_cluster3"]),as.numeric(validation_salivary_microbiome_scores[study_diseased,"high_sHACK_cluster3"]))$p.value
  
  
}

df_score_directions_validation <- as.data.frame(matrix(NA,nrow(df_comparison_sHACK_score),5))
rownames(df_score_directions_validation) <- rownames(df_comparison_sHACK_score)
colnames(df_score_directions_validation) <- c("sHACK","dysbiosis_bray","dysbiosis_canberra","high_sHACK_cluster2","high_sHACK_cluster3")

df_score_directions_validation[,"sHACK"] <- sign(df_comparison_sHACK_score[,"cohen_d_ds_bray"]) * ifelse(df_comparison_sHACK_score[,"p_val_sHACK"] <= 0.05,3,ifelse(df_comparison_sHACK_score[,"p_val_sHACK"] <= 0.1,2,1))

df_score_directions_validation[,"dysbiosis_bray"] <- sign(df_comparison_sHACK_score[,"cohen_d_ds_bray"]) * ifelse(df_comparison_sHACK_score[,"p_val_ds_bray"] <= 0.05,3,ifelse(df_comparison_sHACK_score[,"p_val_ds_bray"] <= 0.1,2,1))

df_score_directions_validation[,"dysbiosis_canberra"] <- sign(df_comparison_sHACK_score[,"cohen_d_ds_canberra"]) * ifelse(df_comparison_sHACK_score[,"p_val_ds_canberra"] <= 0.05,3,ifelse(df_comparison_sHACK_score[,"p_val_ds_canberra"] <= 0.1,2,0))

df_score_directions_validation[,"high_sHACK_cluster2"] <- sign(df_comparison_sHACK_score[,"cohen_d_cluster2"]) * ifelse(df_comparison_sHACK_score[,"p_val_cluster2"] <= 0.05,3,ifelse(df_comparison_sHACK_score[,"p_val_cluster2"] <= 0.1,2,1))

df_score_directions_validation[,"high_sHACK_cluster3"] <- sign(df_comparison_sHACK_score[,"cohen_d_cluster3"]) * ifelse(df_comparison_sHACK_score[,"p_val_cluster3"] <= 0.05,3,ifelse(df_comparison_sHACK_score[,"p_val_cluster3"] <= 0.1,2,1))

REM_saliva_scores_validation <- as.data.frame(matrix(NA,5,3))
rownames(REM_saliva_scores_validation) <- c("sHACK","dysbiosis_bray","dysbiosis_canberra","high_sHACK_cluster2","high_sHACK_cluster3")
colnames(REM_saliva_scores_validation) <- c("estimate","p_value","consistency")

temp_rem <- compute_meta_effsize(validation_salivary_microbiome_scores,"sHACK","Study_Condition_Binary","Study_Name",matched_studies)
REM_saliva_scores_validation["sHACK","estimate"] <- as.numeric(temp_rem$model$beta)
REM_saliva_scores_validation["sHACK","p_value"] <- as.numeric(temp_rem$model$pval)
REM_saliva_scores_validation["sHACK","consistency"] <- length(which(sign(as.vector(temp_rem$df_studies$yi)) == sign(as.numeric(temp_rem$model$b))))/nrow(temp_rem$df_studies)

temp_rem <- compute_meta_effsize(validation_salivary_microbiome_scores,"dysbiosis_bray","Study_Condition_Binary","Study_Name",matched_studies)
REM_saliva_scores_validation["dysbiosis_bray","estimate"] <- as.numeric(temp_rem$model$beta)
REM_saliva_scores_validation["dysbiosis_bray","p_value"] <- as.numeric(temp_rem$model$pval)
REM_saliva_scores_validation["dysbiosis_bray","consistency"] <- length(which(sign(as.vector(temp_rem$df_studies$yi)) == sign(as.numeric(temp_rem$model$b))))/nrow(temp_rem$df_studies)

temp_rem <- compute_meta_effsize(validation_salivary_microbiome_scores,"dysbiosis_canberra","Study_Condition_Binary","Study_Name",matched_studies)
REM_saliva_scores_validation["dysbiosis_canberra","estimate"] <- as.numeric(temp_rem$model$beta)
REM_saliva_scores_validation["dysbiosis_canberra","p_value"] <- as.numeric(temp_rem$model$pval)
REM_saliva_scores_validation["dysbiosis_canberra","consistency"] <- length(which(sign(as.vector(temp_rem$df_studies$yi)) == sign(as.numeric(temp_rem$model$b))))/nrow(temp_rem$df_studies)

temp_rem <- compute_meta_effsize(validation_salivary_microbiome_scores,"high_sHACK_cluster2","Study_Condition_Binary","Study_Name",matched_studies)
REM_saliva_scores_validation["high_sHACK_cluster2","estimate"] <- as.numeric(temp_rem$model$beta)
REM_saliva_scores_validation["high_sHACK_cluster2","p_value"] <- as.numeric(temp_rem$model$pval)
REM_saliva_scores_validation["high_sHACK_cluster2","consistency"] <- length(which(sign(as.vector(temp_rem$df_studies$yi)) == sign(as.numeric(temp_rem$model$b))))/nrow(temp_rem$df_studies)

temp_rem <- compute_meta_effsize(validation_salivary_microbiome_scores,"high_sHACK_cluster3","Study_Condition_Binary","Study_Name",matched_studies)
REM_saliva_scores_validation["high_sHACK_cluster3","estimate"] <- as.numeric(temp_rem$model$beta)
REM_saliva_scores_validation["high_sHACK_cluster3","p_value"] <- as.numeric(temp_rem$model$pval)
REM_saliva_scores_validation["high_sHACK_cluster3","consistency"] <- length(which(sign(as.vector(temp_rem$df_studies$yi)) == sign(as.numeric(temp_rem$model$b))))/nrow(temp_rem$df_studies)

mat <- df_score_directions_validation

heatmap.2(t(mat),density="none",trace="none",col=c("blue3","cornflowerblue","lightskyblue","orange","orange3","orange4"),colsep = 0:(nrow(mat)+1),rowsep = 0:(ncol(mat)+1),sepcolor = "black",sepwidth = c(0.05, 0.05),margins=c(15,10))

REM_saliva_scores_validation$names <- NA
REM_saliva_scores_validation["sHACK","names"] <- "high sHACK (All)"
REM_saliva_scores_validation["dysbiosis_bray","names"] <- "Dysbiosis Score (bray)"
REM_saliva_scores_validation["dysbiosis_canberra","names"] <- "Dysbiosis Score (canberra)"
REM_saliva_scores_validation["high_sHACK_cluster2","names"] <- "high sHACK (Module2)"
REM_saliva_scores_validation["high_sHACK_cluster3","names"] <- "high sHACK (Module3)"

ggplot(REM_saliva_scores_validation,aes(x=estimate,y=-log(p_value,10)))+geom_point()+theme_bw()+geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_text_repel(label=REM_saliva_scores_validation$names,size=6,col=ifelse(REM_saliva_scores_validation$estimate>0,"orange","blue"))+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+geom_hline(yintercept=-log(0.05,10),color="red")





load("~/Desktop/Manuscripts/MetaOral/S9_3combined_Saliva_Longitudinal_SpDf_modified.RData")

df_longitudinal_validation <- data.frame("sHACK"=rowMeans(apply(combined_Saliva_Longitudinal_SpDf_modified[rownames(combined_Saliva_Longitudinal_SpDf_modified),sHACK],2,rank_scale)),"high_sHACK_cluster2"=rowMeans(apply(combined_Saliva_Longitudinal_SpDf_modified[rownames(combined_Saliva_Longitudinal_SpDf_modified),high_sHACK_cluster2],2,rank_scale)),"high_sHACK_cluster3"=rowMeans(apply(combined_Saliva_Longitudinal_SpDf_modified[rownames(combined_Saliva_Longitudinal_SpDf_modified),high_sHACK_cluster3],2,rank_scale)),"bray_dist"=combined_Saliva_Longitudinal_SpDf_modified$bray_dist,"aitchison_dist"=combined_Saliva_Longitudinal_SpDf_modified$aitchison_dist,"study_name"=combined_Saliva_Longitudinal_SpDf_modified$study_name)

df_longitudinal_validation$dysbiosis_bray <- NA
df_longitudinal_validation$dysbiosis_canberra <- NA

longitudinal_study_validation <- unique(combined_Saliva_Longitudinal_SpDf_modified$study_name)

df_comparison_longitudinal <- as.data.frame(matrix(NA,length(longitudinal_study_validation),5))
rownames(df_comparison_longitudinal) <- longitudinal_study_validation
colnames(df_comparison_longitudinal) <- c("sHACK","high_sHACK_cluster2","high_sHACK_cluster3","dysbiosis_bray","dysbiosis_canberra")

for(i in 1:length(longitudinal_study_validation))
{
  study <- longitudinal_study_validation[i]
  all_study_samples <- rownames(df_longitudinal_validation)[df_longitudinal_validation$study_name == study]
  df_temp <- df_longitudinal_validation[all_study_samples,]
  temp_cor <- cor.test(df_temp$sHACK,df_temp$bray_dist,use="pairwise.complete",method="spearman")
  df_comparison_longitudinal[study,"sHACK"] <- sign(as.numeric(temp_cor$estimate)) * ifelse(temp_cor$p.value <= 0.05,3,ifelse(temp_cor$p.value <= 0.1,2,1))
  
  temp_cor <- cor.test(df_temp$high_sHACK_cluster2,df_temp$bray_dist,use="pairwise.complete",method="spearman")
  df_comparison_longitudinal[study,"high_sHACK_cluster2"] <- sign(as.numeric(temp_cor$estimate)) * ifelse(temp_cor$p.value <= 0.05,3,ifelse(temp_cor$p.value <= 0.1,2,1))
  
  temp_cor <- cor.test(df_temp$high_sHACK_cluster3,df_temp$bray_dist,use="pairwise.complete",method="spearman")
  df_comparison_longitudinal[study,"high_sHACK_cluster3"] <- sign(as.numeric(temp_cor$estimate)) * ifelse(temp_cor$p.value <= 0.05,3,ifelse(temp_cor$p.value <= 0.1,2,1))
  
  dist_mat <- as.matrix(vegdist(combined_Saliva_Longitudinal_SpDf_modified[all_study_samples,intersect(colnames(combined_Saliva_Longitudinal_SpDf_modified),rownames(Combined_Saliva_Scores))],method="canberra"))
  diag(dist_mat) <- NA
  df_longitudinal_validation[all_study_samples,"dysbiosis_canberra"] <- apply(dist_mat[all_study_samples,all_study_samples],1,function(x)(mean(x[!is.na(x)])))
  temp_cor <- cor.test(df_longitudinal_validation[all_study_samples,"dysbiosis_canberra"],df_longitudinal_validation[all_study_samples,"bray_dist"],use="pairwise.complete",method="spearman")
  df_comparison_longitudinal[study,"dysbiosis_canberra"] <- sign(as.numeric(temp_cor$estimate)) * ifelse(temp_cor$p.value <= 0.05,3,ifelse(temp_cor$p.value <= 0.1,2,1))
  
  dist_mat <- as.matrix(vegdist(combined_Saliva_Longitudinal_SpDf_modified[all_study_samples,intersect(colnames(combined_Saliva_Longitudinal_SpDf_modified),rownames(Combined_Saliva_Scores))],method="bray"))
  diag(dist_mat) <- NA
  df_longitudinal_validation[all_study_samples,"dysbiosis_bray"] <- apply(dist_mat[all_study_samples,all_study_samples],1,function(x)(mean(x[!is.na(x)])))
  temp_cor <- cor.test(df_longitudinal_validation[all_study_samples,"dysbiosis_bray"],df_longitudinal_validation[all_study_samples,"bray_dist"],use="pairwise.complete",method="spearman")
  df_comparison_longitudinal[study,"dysbiosis_bray"] <- sign(as.numeric(temp_cor$estimate)) * ifelse(temp_cor$p.value <= 0.05,3,ifelse(temp_cor$p.value <= 0.1,2,1))
}


REM_saliva_scores_longitudinal_validation <- as.data.frame(matrix(NA,5,3))
rownames(REM_saliva_scores_longitudinal_validation) <- c("sHACK","dysbiosis_bray","dysbiosis_canberra","high_sHACK_cluster2","high_sHACK_cluster3")
colnames(REM_saliva_scores_longitudinal_validation) <- c("estimate","p_value","consistency")

temp_df_longitudinal_rem <- df_longitudinal_validation[!is.na(df_longitudinal_validation$bray_dist),]

temp_rem <- compute_meta_corr(temp_df_longitudinal_rem,"sHACK","bray_dist","study_name",unique(temp_df_longitudinal_rem$study_name))
REM_saliva_scores_longitudinal_validation["sHACK","estimate"] <- as.numeric(temp_rem$beta)
REM_saliva_scores_longitudinal_validation["sHACK","p_value"] <- as.numeric(temp_rem$pval)
REM_saliva_scores_longitudinal_validation["sHACK","consistency"] <- length(which(sign(df_comparison_longitudinal[,"sHACK"])==sign(as.numeric(temp_rem$beta))))/nrow(df_comparison_longitudinal)

temp_rem <- compute_meta_corr(temp_df_longitudinal_rem,"high_sHACK_cluster2","bray_dist","study_name",unique(temp_df_longitudinal_rem$study_name))
REM_saliva_scores_longitudinal_validation["high_sHACK_cluster2","estimate"] <- as.numeric(temp_rem$beta)
REM_saliva_scores_longitudinal_validation["high_sHACK_cluster2","p_value"] <- as.numeric(temp_rem$pval)
REM_saliva_scores_longitudinal_validation["high_sHACK_cluster2","consistency"] <- length(which(sign(df_comparison_longitudinal[,"high_sHACK_cluster2"])==sign(as.numeric(temp_rem$beta))))/nrow(df_comparison_longitudinal)

temp_rem <- compute_meta_corr(temp_df_longitudinal_rem,"high_sHACK_cluster3","bray_dist","study_name",unique(temp_df_longitudinal_rem$study_name))
REM_saliva_scores_longitudinal_validation["high_sHACK_cluster3","estimate"] <- as.numeric(temp_rem$beta)
REM_saliva_scores_longitudinal_validation["high_sHACK_cluster3","p_value"] <- as.numeric(temp_rem$pval)
REM_saliva_scores_longitudinal_validation["high_sHACK_cluster3","consistency"] <- length(which(sign(df_comparison_longitudinal[,"high_sHACK_cluster3"])==sign(as.numeric(temp_rem$beta))))/nrow(df_comparison_longitudinal)

temp_rem <- compute_meta_corr(temp_df_longitudinal_rem,"dysbiosis_canberra","bray_dist","study_name",unique(temp_df_longitudinal_rem$study_name))
REM_saliva_scores_longitudinal_validation["dysbiosis_canberra","estimate"] <- as.numeric(temp_rem$beta)
REM_saliva_scores_longitudinal_validation["dysbiosis_canberra","p_value"] <- as.numeric(temp_rem$pval)
REM_saliva_scores_longitudinal_validation["dysbiosis_canberra","consistency"] <- length(which(sign(df_comparison_longitudinal[,"dysbiosis_canberra"])==sign(as.numeric(temp_rem$beta))))/nrow(df_comparison_longitudinal)

temp_rem <- compute_meta_corr(temp_df_longitudinal_rem,"dysbiosis_bray","bray_dist","study_name",unique(temp_df_longitudinal_rem$study_name))
REM_saliva_scores_longitudinal_validation["dysbiosis_bray","estimate"] <- as.numeric(temp_rem$beta)
REM_saliva_scores_longitudinal_validation["dysbiosis_bray","p_value"] <- as.numeric(temp_rem$pval)
REM_saliva_scores_longitudinal_validation["dysbiosis_bray","consistency"] <- length(which(sign(df_comparison_longitudinal[,"dysbiosis_bray"])==sign(as.numeric(temp_rem$beta))))/nrow(df_comparison_longitudinal)

longitudinal_taxa <- intersect(rownames(df_pca_saliva),colnames(combined_Saliva_Longitudinal_SpDf_modified))

df_longitudinal_comparison_all_taxa <- as.data.frame(matrix(NA,length(longitudinal_taxa),length(longitudinal_study_validation)))
rownames(df_longitudinal_comparison_all_taxa) <- longitudinal_taxa
colnames(df_longitudinal_comparison_all_taxa) <- longitudinal_study_validation

for(i in 1:length(longitudinal_study_validation))
{
  study <- longitudinal_study_validation[i]
  all_study_samples <- rownames(combined_Saliva_Longitudinal_SpDf_modified[combined_Saliva_Longitudinal_SpDf_modified$study_name == study,])
  
  for(j in 1:length(longitudinal_taxa))
  {
    species_name <- longitudinal_taxa[j]
    temp_corr <- cor.test(combined_Saliva_Longitudinal_SpDf_modified[all_study_samples,"bray_dist"],combined_Saliva_Longitudinal_SpDf_modified[all_study_samples,species_name],method="spearman",use="pairwise.complete")
    
    df_longitudinal_comparison_all_taxa[j,i] <- sign(as.numeric(temp_corr$estimate)) * ifelse(temp_corr$p.value<=0.05,3,ifelse(temp_corr$p.value <= 0.1,2,1))
  }
}

df_longitudinal_association_direction <- data.frame("stability_associated"=apply(df_longitudinal_comparison_all_taxa,1,function(x)(length(x[!is.na(x)&(x<= -2)]))),"cluster"=df_pca_saliva[longitudinal_taxa,"cluster"])

boxplot(df_longitudinal_association_direction$stability_associated~df_longitudinal_association_direction$cluster,col=c("brown","cornflowerblue","darkolivegreen","darkgoldenrod","firebrick4"),outline=FALSE,cex.axis=2,ylab=c("","","","","",""))

dunn_clusters_longitudinal <- dunn.test(df_longitudinal_association_direction$stability_associated,df_longitudinal_association_direction$cluster)







save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/260625_MetaOral_TSG_OS.RData")







# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######         Supragingival subsite association consistency and clustering
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("We extend this to other sub-sites")
print("Supragingival")

supragingival_controls <- rownames(MetadataDf_All_Discovery_filt[(MetadataDf_All_Discovery_filt$body_site_category == "supragingival_subsite")&(MetadataDf_All_Discovery_filt$study_condition == "Control"),])

SpDf_supragingival_Control_associated <- SpDf_All_Discovery_filt[supragingival_controls,rownames(Combined_Supragingival_Scores)]

study_list_supragingival <- unique(MetadataDf_All_Discovery_filt[supragingival_controls,"study_name"])
total_studies_suprgingival <- length(study_list_supragingival)
species_list_supragingival <- rownames(Combined_Supragingival_Scores)

data <- SpDf_supragingival_Control_associated
data$study_name <- MetadataDf_All_Discovery_filt[rownames(data),"study_name"]

supragingival_association_matrix <- compute_association_consistency(data,species_list_supragingival,study_list_supragingival,total_studies_suprgingival)

pca_supragingival_association <- dudi.pca(supragingival_association_matrix$final_matrix,scannf=FALSE,nf=3)
df_pca_supragingival <- data.frame("Axis1"=pca_supragingival_association$li[,1],"Axis2"=pca_supragingival_association$li[,2],"Axis3"=pca_supragingival_association$li[,3],row.names=rownames(pca_supragingival_association$li))

df_pca_supragingival <- df_pca_supragingival[rownames(Combined_Supragingival_Scores),]
all(rownames(df_pca_supragingival) == rownames(Combined_Supragingival_Scores))

df_pca_supragingival$HAC <- Combined_Supragingival_Scores$HAC_Score

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Supragingival_AssociationConsistency_PCA_HACClustering.pdf",width=6,height=6)
ggplot(df_pca_supragingival,aes(x=Axis1,y=Axis2,colour=as.character(ifelse(HAC>=0.90,3,ifelse(HAC>=0.75,2,1)))))+geom_point(size=3)+scale_color_manual(values=c("3"="cornflowerblue","2"="lightblue","1"="red"))+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.position="none")
dev.off()

library(clusterSim)
library(vegan)
optimal_k_pca_supragingival <- iterative_enterotyping(df_pca_supragingival[,c(1,2)],20,20)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Supragingival_AssociationConsistency_CHIndex.pdf",width=6,height=6)
boxplot(optimal_k_pca_supragingival$ch_index[,2],optimal_k_pca_supragingival$ch_index[,3],optimal_k_pca_supragingival$ch_index[,4],optimal_k_pca_supragingival$ch_index[,5],optimal_k_pca_supragingival$ch_index[,6],optimal_k_pca_supragingival$ch_index[,7],optimal_k_pca_supragingival$ch_index[,8],optimal_k_pca_supragingival$ch_index[,9],optimal_k_pca_supragingival$ch_index[,10],names=c("2","3","4","5","6","7","8","9","10"),cex.axis=2,outline=FALSE,ylab="CH-Index",cex.title=2)
dev.off()

set.seed(100); supragingival_kmeans <- kmeans(df_pca_supragingival[,c(1,2)],3); df_pca_supragingival$cluster <- supragingival_kmeans$cluster; s.class(df_pca_supragingival[,c(1,2)],as.factor(df_pca_supragingival$cluster),col="red")

df_pca_supragingival$cluster <- supragingival_kmeans$cluster
cx <- supragingival_kmeans$centers[,1]
cy <- supragingival_kmeans$centers[,2]

df_pca_supragingival$cx <- cx[df_pca_supragingival$cluster]
df_pca_supragingival$cy <- cy[df_pca_supragingival$cluster]

save(df_pca_supragingival, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Supragingival_PCAdf.RData")

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Supragingival_AssociationConsistency_Clusters_HAC_distribution.pdf",width=6,height=6)
ggplot(df_pca_supragingival,aes(x=Axis1,y=Axis2,color=as.character(ifelse(HAC>=0.90,3,ifelse(HAC>=0.75,2,1))))) + geom_point(size=2) + geom_segment(aes(x = cx, y = cy, xend = Axis1, yend = Axis2), color = "grey30", linewidth = 0.3, alpha = 0.8) + stat_ellipse(aes(group = cluster),color = "grey40", linewidth = 0.5, level = 0.95) + scale_color_manual(values=c("3"="cornflowerblue","2"="cyan","1"="red")) + theme_classic() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),legend.position="none") + geom_text_repel(label=ifelse(df_pca_supragingival$HAC>=0.90,rownames(df_pca_supragingival),""),box.padding=0.000001,size=3.5,max.overlaps=40,col="black")
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Supragingival_AssociationConsistency_Cluster_HAC_Scores.pdf",width=6,height=6)
boxplot(df_pca_supragingival$HAC~df_pca_supragingival$cluster,col=c("brown","cornflowerblue","darkolivegreen"),outline=FALSE,cex.axis=2,xlab="",ylab="")
dev.off()

supragingival_cluster_counts <- table(factor(df_pca_supragingival[df_pca_supragingival$HAC >= 0.90, "cluster"],levels = 1:3))
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Supragingival_HAC90_Clusters_Barplot.pdf",width=4,height=6)
barplot(supragingival_cluster_counts,col = "#156082",border = "black",xlab = "Cluster",ylab = "Number of species",ylim = c(0, max(supragingival_cluster_counts) + 5),cex.names = 1.2,cex.axis = 1.1,cex.lab = 1.3,cex.main = 1.2)
dev.off()
## Fishers' exact test.
#Make cluster a factor so cluster 2 is included even if zero
df_pca_supragingival$cluster_factor <- factor(df_pca_supragingival$cluster, levels = 1:3)
## Define high HAC as HAC >= 0.90
df_pca_supragingival$High_HAC <- df_pca_supragingival$HAC >= 0.90
## Define Module 2 vs other modules
df_pca_supragingival$Module2 <- df_pca_supragingival$cluster == 2
## 2 x 2 contingency table
supragingival_fisher_table <- table(High_HAC = df_pca_supragingival$High_HAC,Module2 = df_pca_supragingival$Module2)
supragingival_fisher <- fisher.test(supragingival_fisher_table)
supragingival_fisher$p.value # [1] 2.356991e-11
supragingival_fisher$estimate # odds ratio: 21.15184


write.csv(df_pca_supragingival,file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Supragingival_PCAdf.csv",row.names=TRUE)
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/260625_MetaOral_TSG_OS.RData")




# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######         Subgingival subsite association consistency and clustering
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Subgingival")
subgingival_controls <- rownames(MetadataDf_All_Discovery_filt[(MetadataDf_All_Discovery_filt$body_site_category == "subgingival_subsite")&(MetadataDf_All_Discovery_filt$study_condition == "Control"),])

SpDf_subgingival_Control_associated <- SpDf_All_Discovery_filt[subgingival_controls,rownames(Combined_Subgingival_Scores)]

study_list_subgingival <- unique(MetadataDf_All_Discovery_filt[subgingival_controls,"study_name"])
total_studies_subgingival <- length(study_list_subgingival)
species_list_subgingival <- rownames(Combined_Subgingival_Scores)

data <- SpDf_subgingival_Control_associated
data$study_name <- MetadataDf_All_Discovery_filt[rownames(data),"study_name"]

subgingival_association_matrix <- compute_association_consistency(data,species_list_subgingival,study_list_subgingival,total_studies_subgingival)

pca_subgingival_association <- dudi.pca(subgingival_association_matrix$final_matrix,scannf=FALSE,nf=3)
df_pca_subgingival <- data.frame("Axis1"=pca_subgingival_association$li[,1],"Axis2"=pca_subgingival_association$li[,2],"Axis3"=pca_subgingival_association$li[,3],row.names=rownames(pca_subgingival_association$li))

df_pca_subgingival <- df_pca_subgingival[rownames(Combined_Subgingival_Scores),]
all(rownames(df_pca_subgingival) == rownames(Combined_Subgingival_Scores))

df_pca_subgingival$HAC <- Combined_Subgingival_Scores$HAC_Score


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_AssociationConsistency_PCA_HACClustering.pdf",width=6,height=6)
ggplot(df_pca_subgingival,aes(x=Axis1,y=Axis2,colour=as.character(ifelse(HAC>=0.90,3,ifelse(HAC>=0.75,2,1)))))+geom_point(size=3)+scale_color_manual(values=c("3"="cornflowerblue","2"="lightblue","1"="red"))+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.position="none")
dev.off()

optimal_k_pca_subgingival <- iterative_enterotyping(df_pca_subgingival[,c(1,2)],20,20)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_AssociationConsistency_CHIndex.pdf",width=6,height=6)
boxplot(optimal_k_pca_subgingival$ch_index[,2],optimal_k_pca_subgingival$ch_index[,3],optimal_k_pca_subgingival$ch_index[,4],optimal_k_pca_subgingival$ch_index[,5],optimal_k_pca_subgingival$ch_index[,6],optimal_k_pca_subgingival$ch_index[,7],optimal_k_pca_subgingival$ch_index[,8],optimal_k_pca_subgingival$ch_index[,9],optimal_k_pca_subgingival$ch_index[,10],names=c("2","3","4","5","6","7","8","9","10"),cex.axis=2,outline=FALSE,ylab="CH-Index",cex.title=2)
dev.off()

set.seed(100); subgingival_kmeans <- kmeans(df_pca_subgingival[,c(1,2)],2); df_pca_subgingival$cluster <- subgingival_kmeans$cluster; s.class(df_pca_subgingival[,c(1,2)],as.factor(df_pca_subgingival$cluster),col="red")

df_pca_subgingival$cluster <- subgingival_kmeans$cluster
cx <- subgingival_kmeans$centers[,1]
cy <- subgingival_kmeans$centers[,2]

df_pca_subgingival$cx <- cx[df_pca_subgingival$cluster]
df_pca_subgingival$cy <- cy[df_pca_subgingival$cluster]

save(df_pca_subgingival, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_PCAdf.RData")

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_AssociationConsistency_Clusters_HAC_distribution.pdf",width=6,height=6)
ggplot(df_pca_subgingival,aes(x=Axis1,y=Axis2,color=as.character(ifelse(HAC>=0.90,3,ifelse(HAC>=0.75,2,1))))) + geom_point(size=2) + geom_segment(aes(x = cx, y = cy, xend = Axis1, yend = Axis2), color = "grey30", linewidth = 0.3, alpha = 0.8) + stat_ellipse(aes(group = cluster),color = "grey40", linewidth = 0.5, level = 0.95) + scale_color_manual(values=c("3"="cornflowerblue","2"="cyan","1"="red")) + theme_classic() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),legend.position="none")  + geom_text_repel(label=ifelse(df_pca_subgingival$HAC>=0.90,rownames(df_pca_subgingival),""),box.padding=0.000001,size=3.5,max.overlaps=40,col="black")
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_AssociationConsistency_Cluster_HAC_Scores.pdf",width=6,height=6)
boxplot(df_pca_subgingival$HAC~df_pca_subgingival$cluster,col=c("cornflowerblue","brown"),outline=FALSE,cex.axis=2,xlab="",ylab="")
dev.off()



subgingival_cluster_counts <- table(factor(df_pca_subgingival[df_pca_subgingival$HAC >= 0.90, "cluster"],levels = 1:2))
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_HAC90_Clusters_Barplot.pdf",width=4,height=6)
barplot(subgingival_cluster_counts,col = "#156082",border = "black",xlab = "Cluster",ylab = "Number of species",ylim = c(0, max(subgingival_cluster_counts) + 5),cex.names = 1.2,cex.axis = 1.1,cex.lab = 1.3,cex.main = 1.2)
dev.off()
## Fishers' exact test.
#Make cluster a factor so cluster 2 is included even if zero
df_pca_subgingival$cluster_factor <- factor(df_pca_subgingival$cluster, levels = 1:2)
## Define high HAC as HAC >= 0.90
df_pca_subgingival$High_HAC <- df_pca_subgingival$HAC >= 0.90
## Define Module 1 vs other modules
df_pca_subgingival$Module1 <- df_pca_subgingival$cluster == 1
## 2 x 2 contingency table
subgingival_fisher_table <- table(High_HAC = df_pca_subgingival$High_HAC,Module1 = df_pca_subgingival$Module1)
subgingival_fisher <- fisher.test(subgingival_fisher_table)
subgingival_fisher$p.value # [1] 0.0005041591 or 5.0e-4
subgingival_fisher$estimate # odds ratio: 15.33236


write.csv(df_pca_subgingival,file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Subgingival_PCAdf.csv",row.names=TRUE)





save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/260625_MetaOral_TSG_OS.RData")






# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######         Tongue subsite association consistency and clustering
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Tongue")
tongue_controls <- rownames(MetadataDf_All_Discovery_filt[(MetadataDf_All_Discovery_filt$body_site_category == "tongue_tonsil_subsite")&(MetadataDf_All_Discovery_filt$study_condition == "Control"),])

SpDf_tongue_Control_associated <- SpDf_All_Discovery_filt[tongue_controls,rownames(Combined_TongueTonsil_Scores)]

study_list_tongue <- unique(MetadataDf_All_Discovery_filt[tongue_controls,"study_name"])
total_studies_tongue <- length(study_list_tongue)
species_list_tongue <- rownames(Combined_TongueTonsil_Scores)

data <- SpDf_tongue_Control_associated
data$study_name <- MetadataDf_All_Discovery_filt[rownames(data),"study_name"]

tongue_association_matrix <- compute_association_consistency(data,species_list_tongue,study_list_tongue,total_studies_tongue)

 pca_tongue_association <- dudi.pca(tongue_association_matrix$final_matrix,scannf=FALSE,nf=3)
df_pca_tongue <- data.frame("Axis1"=pca_tongue_association$li[,1],"Axis2"=pca_tongue_association$li[,2],"Axis3"=pca_tongue_association$li[,3],row.names=rownames(pca_tongue_association$li))

df_pca_tongue <- df_pca_tongue[rownames(Combined_TongueTonsil_Scores),]
all(rownames(df_pca_tongue) == rownames(Combined_TongueTonsil_Scores))
df_pca_tongue$HAC <- Combined_TongueTonsil_Scores$HAC_Score

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Tongue_AssociationConsistency_PCA_HACClustering.pdf",width=6,height=6)
ggplot(df_pca_tongue,aes(x=Axis1,y=Axis2,colour=as.character(ifelse(HAC>=0.90,3,ifelse(HAC>=0.75,2,1)))))+geom_point(size=3)+scale_color_manual(values=c("3"="cornflowerblue","2"="lightblue","1"="red"))+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.position="none")
dev.off()
optimal_k_pca_tongue <- iterative_enterotyping(df_pca_tongue[,c(1,2)],20,20)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Tongue_AssociationConsistency_CHIndex.pdf",width=6,height=6)
boxplot(optimal_k_pca_tongue$ch_index[,2],optimal_k_pca_tongue$ch_index[,3],optimal_k_pca_tongue$ch_index[,4],optimal_k_pca_tongue$ch_index[,5],optimal_k_pca_tongue$ch_index[,6],optimal_k_pca_tongue$ch_index[,7],optimal_k_pca_tongue$ch_index[,8],optimal_k_pca_tongue$ch_index[,9],optimal_k_pca_tongue$ch_index[,10],names=c("2","3","4","5","6","7","8","9","10"),cex.axis=2,outline=FALSE,ylab="CH-Index",cex.title=2)
dev.off()

set.seed(100); tongue_kmeans <- kmeans(df_pca_tongue[,c(1,2)],6); df_pca_tongue$cluster <- tongue_kmeans$cluster; s.class(df_pca_tongue[,c(1,2)],as.factor(df_pca_tongue$cluster),col="red")

df_pca_tongue$cluster <- tongue_kmeans$cluster
cx <- tongue_kmeans$centers[,1]
cy <- tongue_kmeans$centers[,2]

df_pca_tongue$cx <- cx[df_pca_tongue$cluster]
df_pca_tongue$cy <- cy[df_pca_tongue$cluster]

save(df_pca_tongue, file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Tongue_PCAdf.RData")

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Tongue_AssociationConsistency_Clusters_HAC_distribution.pdf",width=6,height=6)
ggplot(df_pca_tongue,aes(x=Axis1,y=Axis2,color=as.character(ifelse(HAC>=0.90,3,ifelse(HAC>=0.75,2,1))))) + geom_point(size=2) + geom_segment(aes(x = cx, y = cy, xend = Axis1, yend = Axis2), color = "grey30", linewidth = 0.3, alpha = 0.8) + stat_ellipse(aes(group = cluster),color = "grey40", linewidth = 0.5, level = 0.95) + scale_color_manual(values=c("3"="cornflowerblue","2"="cyan","1"="red")) + theme_classic() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),legend.position="none") + geom_text_repel(label=ifelse(df_pca_tongue$HAC>=0.90,rownames(df_pca_tongue),""),box.padding=0.000001,size=3.5,max.overlaps=40,col="black")
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Tongue_AssociationConsistency_Cluster_HAC_Scores.pdf",width=6,height=6)
boxplot(df_pca_tongue$HAC~df_pca_tongue$cluster,col=c("brown","darkolivegreen","darkgoldenrod","firebrick4","indianred","cornflowerblue"),outline=FALSE,cex.axis=2,xlab="",ylab="")
dev.off()
dunn.test(df_pca_tongue$HAC,df_pca_tongue$cluster)$P
#   Kruskal-Wallis rank sum test

# data: x and group
# Kruskal-Wallis chi-squared = 34.9407, df = 5, p-value = 0


#                     Dunn's Pairwise Comparison of x by group                    
#                                  (No adjustment)                                

# Col Mean-│
# Row Mean │          1          2          3          4          5
# ─────────┼───────────────────────────────────────────────────────
#        2 │  -2.843860
#          │     0.0022*
#          │
#        3 │  -0.504352   2.451255
#          │     0.3070     0.0071*
#          │
#        4 │  -2.572284  -0.109850  -2.202458
#          │     0.0051*    0.4563     0.0138*
#          │
#        5 │  -2.977953  -0.289772  -2.603975  -0.140785
#          │     0.0015*    0.3860     0.0046*    0.4440 
#          │
#        6 │  -5.108471  -3.260425  -4.867389  -2.897488  -2.971338
#          │     0.0000*    0.0006*    0.0000*    0.0019*    0.0015*

# α = 0.05
# Reject Ho if p ≤ α/2, where p = Pr(Z ≥ |z|)
#  [1] 2.228529e-03 3.070069e-01 7.117940e-03 5.051496e-03 4.562641e-01
#  [6] 1.381648e-02 1.450898e-03 3.859955e-01 4.607469e-03 4.440197e-01
# [11] 1.623877e-07 5.562264e-04 5.654113e-07 1.880816e-03 1.482524e-03


tongue_cluster_counts <- table(factor(df_pca_tongue[df_pca_tongue$HAC >= 0.90, "cluster"],levels = 1:6))
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Tongue_HAC90_Clusters_Barplot.pdf",width=4,height=6)
barplot(tongue_cluster_counts,col = "#156082",border = "black",xlab = "Cluster",ylab = "Number of species",ylim = c(0, max(tongue_cluster_counts) + 5),cex.names = 1.2,cex.axis = 1.1,cex.lab = 1.3,cex.main = 1.2)
dev.off()
## Fishers' exact test.
#Make cluster a factor so cluster 2 is included even if zero
df_pca_tongue$cluster_factor <- factor(df_pca_tongue$cluster, levels = 1:6)
## Define high HAC as HAC >= 0.90
df_pca_tongue$High_HAC <- df_pca_tongue$HAC >= 0.90
## Define Module 6 vs other modules
df_pca_tongue$Module6 <- df_pca_tongue$cluster == 6
## 2 x 2 contingency table
tongue_fisher_table <- table(High_HAC = df_pca_tongue$High_HAC,Module6 = df_pca_tongue$Module6)
tongue_fisher <- fisher.test(tongue_fisher_table)
tongue_fisher$p.value # [1] 0.181005
tongue_fisher$estimate # odds ratio: 2.366708

write.csv(df_pca_tongue,file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Tongue_PCAdf.csv",row.names=TRUE)


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/260625_MetaOral_TSG_OS.RData")


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######         Buccal subsite association consistency and clustering
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Buccal")
buccal_controls <- rownames(MetadataDf_All_Discovery_filt[(MetadataDf_All_Discovery_filt$body_site_category == "buccal_palate_other_surface_subsite")&(MetadataDf_All_Discovery_filt$study_condition == "Control"),])

SpDf_buccal_Control_associated <- SpDf_All_Discovery_filt[buccal_controls,rownames(Combined_BuccalPalate_Scores)]

study_list_buccal <- unique(MetadataDf_All_Discovery_filt[buccal_controls,"study_name"])
total_studies_buccal <- length(study_list_buccal)
species_list_buccal <- rownames(Combined_BuccalPalate_Scores)

data <- SpDf_buccal_Control_associated
data$study_name <- MetadataDf_All_Discovery_filt[rownames(data),"study_name"]

buccal_association_matrix <- compute_association_consistency(data,species_list_buccal,study_list_buccal,total_studies_buccal)

pca_buccal_association <- dudi.pca(buccal_association_matrix$final_matrix,scannf=FALSE,nf=3)
df_pca_buccal <- data.frame("Axis1"=pca_buccal_association$li[,1],"Axis2"=pca_buccal_association$li[,2],"Axis3"=pca_buccal_association$li[,3],row.names=rownames(pca_buccal_association$li))

df_pca_buccal <- df_pca_buccal[rownames(Combined_BuccalPalate_Scores),]
all(rownames(df_pca_buccal) == rownames(Combined_BuccalPalate_Scores))
df_pca_buccal$HAC <- Combined_BuccalPalate_Scores$CoreAssociationScore

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Buccal_AssociationConsistency_PCA_HACClustering.pdf",width=6,height=6)
ggplot(df_pca_buccal,aes(x=Axis1,y=Axis2,colour=as.character(ifelse(HAC>=0.90,3,ifelse(HAC>=0.75,2,1)))))+geom_point(size=3)+scale_color_manual(values=c("3"="cornflowerblue","2"="lightblue","1"="red"))+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),legend.position="none")
dev.off()
optimal_k_pca_buccal <- iterative_enterotyping(df_pca_buccal[,c(1,2)],20,50)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Buccal_AssociationConsistency_CHIndex.pdf",width=6,height=6)
boxplot(optimal_k_pca_buccal$ch_index[,2],optimal_k_pca_buccal$ch_index[,3],optimal_k_pca_buccal$ch_index[,4],optimal_k_pca_buccal$ch_index[,5],optimal_k_pca_buccal$ch_index[,6],optimal_k_pca_buccal$ch_index[,7],optimal_k_pca_buccal$ch_index[,8],optimal_k_pca_buccal$ch_index[,9],optimal_k_pca_buccal$ch_index[,10],names=c("2","3","4","5","6","7","8","9","10"),cex.axis=2,outline=FALSE,ylab="CH-Index",cex.title=2)
dev.off()

set.seed(100); buccal_kmeans <- kmeans(df_pca_buccal[,c(1,2)],3); df_pca_buccal$cluster <- buccal_kmeans$cluster; s.class(df_pca_buccal[,c(1,2)],as.factor(df_pca_buccal$cluster),col="red")

df_pca_buccal$cluster <- buccal_kmeans$cluster
cx <- buccal_kmeans$centers[,1]
cy <- buccal_kmeans$centers[,2]

df_pca_buccal$cx <- cx[df_pca_buccal$cluster]
df_pca_buccal$cy <- cy[df_pca_buccal$cluster]

save(df_pca_buccal, file="/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Buccal_PCAdf.RData")


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Buccal_AssociationConsistency_Clusters_HAC_distribution.pdf",width=6,height=6)
ggplot(df_pca_buccal,aes(x=Axis1,y=Axis2,color=as.character(ifelse(HAC>=0.90,3,ifelse(HAC>=0.75,2,1))))) + geom_point(size=2) + geom_segment(aes(x = cx, y = cy, xend = Axis1, yend = Axis2), color = "grey30", linewidth = 0.3, alpha = 0.8) + stat_ellipse(aes(group = cluster),color = "grey40", linewidth = 0.5, level = 0.95) + scale_color_manual(values=c("3"="cornflowerblue","2"="cyan","1"="red")) + theme_classic() + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),legend.position="none")
dev.off()

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Buccal_AssociationConsistency_Cluster_HAC_Scores.pdf",width=6,height=6)
boxplot(df_pca_buccal$HAC~df_pca_buccal$cluster,col=c("brown","cornflowerblue","darkolivegreen"),outline=FALSE,cex.axis=2,xlab="",ylab="")
dev.off()


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/260625_MetaOral_TSG_OS.RData")
