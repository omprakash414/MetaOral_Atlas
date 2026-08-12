
## This script will make plots of core-species, health species, HACK/HAC species in atleast one subsite. 

###############
#################### Now plot the PCoA by taking subsite-specific study v/s species core matrix in one matrix (Also plot the heatmap).

library(dplyr)
## Import the core matrix for each subsite and combine them into one matrix
saliva_core_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_StudySpecific_CoreSpecies.csv", row.names = 1)))
supragingival_core_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1supragingival_StudySpecific_CoreSpecies.csv", row.names = 1)))
subgingival_core_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1subgingival_StudySpecific_CoreSpecies.csv", row.names = 1)))
buccal_palate_core_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1buccalpalate_StudySpecific_CoreSpecies.csv", row.names = 1)))
tongue_tonsil_core_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1tonguetonsil_StudySpecific_CoreSpecies.csv", row.names = 1)))

## Now add the subsite name to the study names in each matrix so that we can identify them later
rownames(saliva_core_matrix) <- paste("SalivaSputum", rownames(saliva_core_matrix), sep = "_")
rownames(supragingival_core_matrix) <- paste("Supragingival", rownames(supragingival_core_matrix), sep = "_")
rownames(subgingival_core_matrix) <- paste("Subgingival", rownames(subgingival_core_matrix), sep = "_")
rownames(buccal_palate_core_matrix) <- paste("BuccalPalate", rownames(buccal_palate_core_matrix), sep = "_")
rownames(tongue_tonsil_core_matrix) <- paste("TongueTonsil", rownames(tongue_tonsil_core_matrix), sep = "_")

## Now merge them into one matrix
combined_core_matrix <- bind_rows(saliva_core_matrix, supragingival_core_matrix, subgingival_core_matrix, tongue_tonsil_core_matrix)
combined_core_matrix[is.na(combined_core_matrix)] <- 0

## Now filter only species that have HACK/HAC score >= 0.80 (Already done in S6_4 script. Just import it here)
common_species <- read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4common_species_HAC_90score_AllSubsites.csv")
common_species <- common_species$x

combined_core_matrix_filtered <- combined_core_matrix[, colnames(combined_core_matrix) %in% common_species] 


## heatmap
library(pheatmap)
my_palette <- colorRampPalette(c("white","darkturquoise"))(20)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_5heatmap_CoreKeyStones_StudyWise_AllSubsites.pdf", width = 35, height = 28)
pheatmap(t(combined_core_matrix_filtered),
         color = my_palette,
         fontsize_row = 10,
         fontsize_col = 10,
         cellheight = 12,
         cellwidth = 12,
         cluster_rows = T,
         treeheight_row = 20,
         cluster_cols = T,
         treeheight_col = 20,
         border_color = "black"
)

dev.off()


## Now do PCoA
library(vegan)
library(ade4)
## distance between subsites (rows)
Euc_dist_mat <- vegdist(combined_core_matrix_filtered, method = "euclidean")

## run PCoA
pcoa_Euc <- dudi.pco(as.dist(Euc_dist_mat), scannf = FALSE, nf = 20)

## extract coordinates
pcoa_df <- as.data.frame(pcoa_Euc$li)
pcoa_df$subsite <- sub("_.*", "", rownames(pcoa_df))

pcoa_df$subsite <- as.factor(pcoa_df$subsite)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_5PCoA_CoreKeyStones_StudyWise_AllSubsites.pdf",width = 7.5,height = 5.5)
par(lwd = 2)   # increase line width
s.class(
  pcoa_df[, 1:2],
  fac = pcoa_df$subsite,
  col = c(
    "SalivaSputum" = "seagreen",
    "Supragingival" = "orangered",
    "TongueTonsil" = "orange",
    "Subgingival" = "deeppink",
    "BuccalPalate" = "blue"
  ),
  cell = 1,
  cpoint = 2.5,   # increase dot size
  xax = 1,
  yax = 2,
  clabel = 1.08
)

dev.off()
adonis2(pcoa_df[,1:2] ~ subsite, data = pcoa_df, permutations = 999, method = "euclidean")
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999

# adonis2(formula = pcoa_df[, 1:2] ~ subsite, data = pcoa_df, permutations = 999, method = "euclidean")
#           Df SumOfSqs      R2    F Pr(>F)    
# Model      3    83.43 0.16496 9.68  0.001 ***
# Residual 147   422.34 0.83504                
# Total    150   505.77 1.00000                
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_5Heatmap_PCoA_AllSubsites.RData")




###############
#################### Now plot the similar heatmap and PCoA for health species (HACK/HAC score >= 0.80) in atleast one subsite.
## load the study v/s species association df for each subsite.
Saliva_health_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_StudyWise_SpeciesAssociation.csv", row.names = 1)))
Supragingival_health_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Supragingival_StudyWise_SpeciesAssociation.csv", row.names = 1)))
Subgingival_health_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Subgingival_StudyWise_SpeciesAssociation.csv", row.names = 1)))
BuccalPalate_health_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1BuccalPalate_StudyWise_SpeciesAssociation.csv", row.names = 1)))
TongueTonsil_health_matrix <- data.frame(t(read.csv("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1TongueTonsil_StudyWise_SpeciesAssociation.csv", row.names = 1)))


## Now add the subsite name to the study names in each matrix so that we can identify them later
rownames(Saliva_health_matrix) <- paste("SalivaSputum", rownames(Saliva_health_matrix), sep = "_")
rownames(Supragingival_health_matrix) <- paste("Supragingival", rownames(Supragingival_health_matrix), sep = "_")
rownames(Subgingival_health_matrix) <- paste("Subgingival", rownames(Subgingival_health_matrix), sep = "_")
rownames(BuccalPalate_health_matrix) <- paste("BuccalPalate", rownames(BuccalPalate_health_matrix), sep = "_")
rownames(TongueTonsil_health_matrix) <- paste("TongueTonsil", rownames(TongueTonsil_health_matrix), sep = "_")

## Now merge then into one matrix
combined_health_matrix <- bind_rows(Saliva_health_matrix, Supragingival_health_matrix,Subgingival_health_matrix, TongueTonsil_health_matrix)
combined_health_matrix[is.na(combined_health_matrix)] <- 0
combined_health_matrix[combined_health_matrix == -1 | combined_health_matrix == 1] <- 0

# check  if there is no any -1 and 1 in the df (It should return FALSE)
any(combined_health_matrix == -1, na.rm = TRUE)
any(combined_health_matrix == 1, na.rm = TRUE)

# now convert -3 and -2 into -1 and 3 and 2 into 1 (to make it binary) (As during health association scores we are considering only number of +ves and number of -ves and not the actual 2 and 3)
combined_health_matrix[combined_health_matrix == -2 | combined_health_matrix == -3] <- -1
combined_health_matrix[combined_health_matrix == 2 | combined_health_matrix == 3] <- 1

## Now filter the species that have HACK/HAC score >= 0.80 (Already done in S6_4 script)
combined_health_matrix_filtered <- combined_health_matrix[,colnames(combined_health_matrix) %in% common_species]


## heatmap
library(pheatmap)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_5heatmap_HealthAssociations_StudyWise_AllSubsites.pdf", width = 35, height = 28)
pheatmap(t(combined_health_matrix_filtered),
         color = c("red", "white", "skyblue"),
         fontsize_row = 10,
         fontsize_col = 10,
         cellheight = 12,
         cellwidth = 12,
         cluster_rows = T,
         treeheight_row = 20,
         cluster_cols = T,
         treeheight_col = 20,
         border_color = "black")
dev.off()


## Now do PCoA
library(vegan)
library(ade4)
## distance between subsites (rows)
Euc_dist_mat2 <- vegdist(combined_health_matrix_filtered, method = "euclidean")

## run PCoA
pcoa_Euc2 <- dudi.pco(as.dist(Euc_dist_mat2), scannf = FALSE, nf = 20)

## extract coordinates
pcoa_df2 <- as.data.frame(pcoa_Euc2$li)
pcoa_df2$subsite <- sub("_.*", "", rownames(pcoa_df2))

pcoa_df2$subsite <- as.factor(pcoa_df2$subsite)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_5PCoA_HealthAssociations_StudyWise_AllSubsites.pdf",width = 7.5,height = 5.5)
par(lwd = 2)   # increases line width
s.class(
  pcoa_df2[, 1:2],
  fac = pcoa_df2$subsite,
  col = c(
    "SalivaSputum" = "seagreen",
    "Supragingival" = "orangered",
    "TongueTonsil" = "orange",
    "Subgingival" = "deeppink",
    "BuccalPalate" = "blue"
  ),
  cell = 1,
  cpoint = 2.5,   # increase dot size
  xax = 1,
  yax = 2,
  clabel = 1.08
)

dev.off()


adonis2(pcoa_df2[,1:2] ~ subsite, data = pcoa_df2, permutations = 999, method = "euclidean")
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 999

# adonis2(formula = pcoa_df2[, 1:2] ~ subsite, data = pcoa_df2, permutations = 999, method = "euclidean")
#          Df SumOfSqs      R2      F Pr(>F)  
# Model     3   28.946 0.12087 2.6124  0.021 *
# Residual 57  210.527 0.87913                
# Total    60  239.473 1.00000                
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_5Heatmap_PCoA_AllSubsites.RData")


###############
####################
