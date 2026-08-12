
### This script imports all the different scores for all the subsites. and also for saliva subsite there are cohort wise scores as well. Are are modified and saved in a single workspace.


load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedScores.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2Subgingival_CombinedScores.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2Supragingival_CombinedScores.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_2TongueTonsil_CombinedScores.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1buccal_palate_CoreAssociationScore.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_3Saliva_CohortWise_Combined_Score_Dfs.RData")

names(Combined_Saliva_Scores) <- c("CoreAssociationScore", "HealthAssociationScore", "StabilityAssociationScore", "HACK_Score_UnRanked", "HACK_Score")
names(Combined_Subgingival_Scores) <- c("CoreAssociationScore", "HealthAssociationScore", "HAC_Score_UnRanked", "HAC_Score")
names(Combined_Supragingival_Scores) <- c("CoreAssociationScore", "HealthAssociationScore", "HAC_Score_UnRanked", "HAC_Score")
names(Combined_TongueTonsil_Scores) <- c("CoreAssociationScore", "HealthAssociationScore", "HAC_Score_UnRanked", "HAC_Score")

Combined_BuccalPalate_Scores <- buccal_palate_CoreAssociationScore[,2,drop=F]
rm(buccal_palate_CoreAssociationScore)

names(Combined_saliva_Exposure_HAC_score) <- c("HealthAssociationScore_16s_exposure", "CoreAssociationScore_16s_exposure", "HAC_Score_UnRanked", "HAC_Score", "All_HACK_Score")
names(Combined_saliva_16s_HACK_score) <- c("HealthAssociationScore_16s", "StabilityAssociationScore_16s", "CoreAssociationScore_16s", "HACK_Score_UnRanked", "HACK_Score", "All_HACK_Score")
names(Combined_saliva_WGS_HACK_score) <- c("HealthAssociationScore_WGS", "StabilityAssociationScore_WGS", "CoreAssociationScore_WGS", "HACK_Score_UnRanked", "HACK_Score", "All_HACK_Score")

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4_AllScores_AllSubsites.RData")


Combined_Saliva_Scores
#####################
####################### Plot the heatmap for species that have HAC score >= 0.90. (For saliva subsite -  import the HAC score. Plot sHACK score as well but do not select species based on the sHACK score. For other subsites along with saliva - import the HAC score and select species based on HAC score >= 0.90 )
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_HAC_Score_noHACK.RData")

saliva_SpSelect <- saliva_HAC_Score[saliva_HAC_Score$HAC_Score_RankScaled >= 0.90,"HAC_Score_RankScaled",drop = F]
supragingival_SpSelect <- Combined_Supragingival_Scores[Combined_Supragingival_Scores$HAC_Score >= 0.90,"HAC_Score",drop = F]
subgingival_SpSelect <- Combined_Subgingival_Scores[Combined_Subgingival_Scores$HAC_Score >= 0.90, "HAC_Score", drop = F]
tongue_SpSelect <- Combined_TongueTonsil_Scores[Combined_TongueTonsil_Scores$HAC_Score >= 0.90, "HAC_Score", drop = F]

common_species <- unique(c(rownames(saliva_SpSelect), rownames(supragingival_SpSelect), rownames(subgingival_SpSelect),rownames(tongue_SpSelect)))
write.csv(common_species, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4common_species_HAC_90score_AllSubsites.csv", row.names = FALSE)


combined_HACK_HAC <- data.frame(matrix(NA,length(common_species),4))
rownames(combined_HACK_HAC) <- common_species
colnames(combined_HACK_HAC) <- c("saliva-sputum","supragingival","subgingival","tongue")

combined_HACK_HAC[rownames(saliva_SpSelect),"saliva"] <- saliva_SpSelect$HAC_Score_RankScaled
combined_HACK_HAC[rownames(supragingival_SpSelect),"supragingival"] <- supragingival_SpSelect$HAC_Score
combined_HACK_HAC[rownames(subgingival_SpSelect),"subgingival"] <- subgingival_SpSelect$HAC_Score
combined_HACK_HAC[rownames(tongue_SpSelect),"tongue"] <- tongue_SpSelect$HAC_Score

combined_HACK_HAC[is.na(combined_HACK_HAC)] <- 0

#####  Plot the heatmap
library(pheatmap)

# create a matrix of the same dimensions as combined_HACK_HAC with formatted values for display (i.e only non-zero values to be displayed.)
numbers_mat <- matrix(
  sprintf("%.2f", as.matrix(combined_HACK_HAC)),
  nrow = nrow(combined_HACK_HAC),
  ncol = ncol(combined_HACK_HAC),
  dimnames = dimnames(combined_HACK_HAC))

numbers_mat[is.na(combined_HACK_HAC)] <- ""
numbers_mat[combined_HACK_HAC == 0] <- ""

mat <- as.matrix(combined_HACK_HAC)

# Calculate the minimum and maximum values for the color scale, excluding NA and zero values
vals <- mat[!is.na(mat) & mat != 0]
min_val <- min(vals)
max_val <- max(vals)


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4combined_HAC_90score_AllSubsites.pdf", width = 10, height = 22)

pheatmap(
  mat,
  color = colorRampPalette(c("white", "indianred"))(100),
  breaks = seq(min_val, max_val, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  na_col = "grey90",
  border_color = "black",
  display_numbers = numbers_mat,
  fontsize_number = 12,
  number_color = "black",
  cellwidth = 38,
  treeheight_row = 0,
  treeheight_col = 0)

dev.off()

# exteact the carper
ph <- pheatmap(mat,cluster_rows = TRUE,cluster_cols = TRUE,silent = TRUE)
# extract row order
ph_row_order <- ph$tree_row$order
# extract column order
ph_col_order <- ph$tree_col$order
# Now you can use ph_row_order and ph_col_order to access the order of rows and columns in the heatmap. and create a new df
mat_carpet <- mat[ph_row_order, ph_col_order]

write.csv(mat_carpet, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4combined_HAC_90score_AllSubsites_Carpet.csv", row.names = TRUE)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4_AllScores_AllSubsites_workspace.RData")


###### Smilarly get the carpet for saliva subsite sHACK score only. i.e single column heatmap for saliva subsite sHACK score. (Do not select species based on sHACK score. Select species based on HAC score >= 0.90 but select the common species across all the subsites.))
# First take all the species present in the above heatmap
saliva_sHACK_SpSelect <- Combined_Saliva_Scores[rownames(mat_carpet), "HACK_Score", drop = F]
# keep all species thata re common across al subsites.

saliva_sHACK_SpSelect <- as.matrix(saliva_sHACK_SpSelect)
## Give column name
colnames(saliva_sHACK_SpSelect) <- "sHACK"
saliva_sHACK_SpSelect <- round(saliva_sHACK_SpSelect, 2)

saliva_sHACK_SpSelect2 <- saliva_sHACK_SpSelect

min_val2 <- min(saliva_sHACK_SpSelect2, na.rm = TRUE)
max_val2 <- max(saliva_sHACK_SpSelect2, na.rm = TRUE)


# now plot the heatmap for saliva subsite sHACK score only. (Do not select species based on sHACK score. Select species based on HAC score >= 0.90)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4saliva_sHACKscore_saliva_HAC90.pdf",width = 5,height = 22)

pheatmap(
  saliva_sHACK_SpSelect2,
  color = colorRampPalette(c("white", "#c7e9c0", "#238b45"))(100),
  breaks = seq(min_val2, max_val2, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  na_col = "white",
  border_color = "black",
  display_numbers = saliva_sHACK_SpSelect,
  number_format = "%.2f",
  fontsize_number = 12,
  number_color = "black",
  cellwidth = 40,
  treeheight_row = 0,
  treeheight_col = 0
)

dev.off()

################################
################### Export the species names as a txt file for which HAC>=0.90. Do this separately for each of the subsite.
supragingival_species90 <- gsub("_"," ",rownames(combined_HACK_HAC[combined_HACK_HAC$supragingival >= 0.90, ]))
subgingival_species90 <- gsub("_"," ",rownames(combined_HACK_HAC[combined_HACK_HAC$subgingival >= 0.90, ]))
tongue_species90 <- gsub("_"," ",rownames(combined_HACK_HAC[combined_HACK_HAC$tongue >= 0.90, ]))

write.table(supragingival_species90, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4supragingival_species_HAC90.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(subgingival_species90, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4subgingival_species_HAC90.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(tongue_species90, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4tongue_species_HAC90.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4_AllScores_AllSubsites_workspace.RData")

### do the correlation of HAC and sHACK for saliva subsite.
cor.test(saliva_HAC_Score[common_species,"HAC_Score_RankScaled"], Combined_Saliva_Scores[common_species,"HACK_Score"], method = "pearson")
#         Pearson's product-moment correlation

# data:  saliva_HAC_Score[common_species, "HAC_Score_RankScaled"] and Combined_Saliva_Scores[common_species, "HACK_Score"]
# t = 10.043, df = 65, p-value = 7.511e-15
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6639786 0.8590969
# sample estimates:
#       cor 
# 0.7798123



### On local laptop == Subsite-specific and subsite shared taxa  here  S6_4combined_HAC_90score_AllSubsites_Carpet.csv  is same as mat_carpet
library(eulerr)
library(polylabelr)
S6_4combined_HAC_90score_AllSubsites_Carpet <- read.csv("S6_4combined_HAC_90score_AllSubsites_Carpet.csv", row.names=1)

colSums(S6_4combined_HAC_90score_AllSubsites_Carpet!= 0, na.rm = TRUE)

saliva_sp <- rownames(S6_4combined_HAC_90score_AllSubsites_Carpet[S6_4combined_HAC_90score_AllSubsites_Carpet$SalivaSputum !=0,])
supragingival_sp <- rownames(S6_4combined_HAC_90score_AllSubsites_Carpet[S6_4combined_HAC_90score_AllSubsites_Carpet$Supragingival !=0,])
subgingival_sp <- rownames(S6_4combined_HAC_90score_AllSubsites_Carpet[S6_4combined_HAC_90score_AllSubsites_Carpet$Subgingival !=0,])
tonguetonsil_sp <- rownames(S6_4combined_HAC_90score_AllSubsites_Carpet[S6_4combined_HAC_90score_AllSubsites_Carpet$TongueTonsil !=0,])

vienn_info <- euler(list("SalivaSputum"=saliva_sp, "Supragingival"=supragingival_sp, "Subgingival"=subgingival_sp, "TongueTonsil"=tonguetonsil_sp), shape = "ellipse")

pdf("S6_4VennPlot_HAC_AllSubsites.pdf", width = 5,height = 5)
plot(vienn_info, labels = list(font = 1.2, cex = 1.4), quantities = list(cex = 1.3,font = 2))
dev.off()

print(vienn_info)














