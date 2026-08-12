

# Load the species profile and metadata for buccal_palate subsite
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_BuccalPalate_DisCtrlCohort.RData")

# Load the species list for buccal_palate subsite
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")

# Load all the functions
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")



##### Prepare the data 
## get the study specific summary (80-20)
buccal_palate_summary_df <- make_summary_single_subsite(MetadataDf_buccal_palate_DisCtrl)
buccal_palate_summary_df2 <- buccal_palate_summary_df$summary_raw
buccal_palate_summary_df3 <- buccal_palate_summary_df$summary_filtered 


## Get the species profile for buccal_palate subsite
SpDf_buccal_palate_DisCtrl <- SpDf_buccal_palate_DisCtrl[,colnames(SpDf_buccal_palate_DisCtrl)%in% buccal_AssociatedSpecies]
SpDf_buccal_palate_DisCtrl <- SpDf_buccal_palate_DisCtrl/rowSums(SpDf_buccal_palate_DisCtrl)

## get species profile and metadata for only studies that are balanced with control-disease samples.
MetadataDf_buccal_palate_DisCtrl <- MetadataDf_buccal_palate_DisCtrl[MetadataDf_buccal_palate_DisCtrl$study_name %in% buccal_palate_summary_df2$study_name,]
SpDf_buccal_palate_DisCtrl <- SpDf_buccal_palate_DisCtrl[rownames(MetadataDf_buccal_palate_DisCtrl),]
AllControlSamples <- rownames(MetadataDf_buccal_palate_DisCtrl[MetadataDf_buccal_palate_DisCtrl$study_condition == "Control",])
AllDiseaseSamples <- rownames(MetadataDf_buccal_palate_DisCtrl[MetadataDf_buccal_palate_DisCtrl$study_condition != "Control",])

selected_studies <- buccal_palate_summary_df2$study_name
##### Calculate the health association score for each species in buccal_palate subsite

Buccal_DiseaseAnalysis_single <- healthAssociation_iterations(1,1,selected_studies,AllControlSamples,AllDiseaseSamples,MetadataDf_buccal_palate_DisCtrl,SpDf_buccal_palate_DisCtrl,buccal_AssociatedSpecies)

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Buccal_palate_HealthAssociation_Workspace.RData")


### save study v/s species association df 
buccal_palate_StudyWise_SpeciesAssociation <- Buccal_DiseaseAnalysis_single$df_comparison_last[[1]]

write.csv(buccal_palate_StudyWise_SpeciesAssociation, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1BuccalPalate_StudyWise_SpeciesAssociation.csv", row.names = T)


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Buccal_palate_HealthAssociation_Workspace.RData")





SpDf_buccal_palate_DisCtrl[is.na(SpDf_buccal_palate_DisCtrl)] <- 0

########### get the normal wilcox and cohen's d vlaues for each species in each of the study.
JulianM_metadata <-MetadataDf_buccal_palate_DisCtrl[MetadataDf_buccal_palate_DisCtrl$study_name == "JulianM_2017",]
JulianM_SpDf <- SpDf_buccal_palate_DisCtrl[rownames(JulianM_metadata),]
JulianM_control_samples <- rownames(JulianM_metadata[JulianM_metadata$study_condition == "Control",])
JulianM_disease_samples <- rownames(JulianM_metadata[JulianM_metadata$study_condition != "Control",])

df_comparison <- data.frame(matrix(NA, nrow = length(buccal_AssociatedSpecies), ncol = 2))
rownames(df_comparison) <- buccal_AssociatedSpecies
colnames(df_comparison) <- c("JulianM_2017", "StehlikovaZ_2019")

for(species_name in buccal_AssociatedSpecies)
      {
        control_values <- as.numeric(JulianM_SpDf[JulianM_control_samples,species_name])
        disease_values <-as.numeric(JulianM_SpDf[JulianM_disease_samples,species_name])

        temp_cohen_d <- cohen.d(control_values,disease_values)$estimate
        temp_wilcox <- wilcox.test(control_values,disease_values)$p.value
        df_comparison[species_name,1] <- sign(temp_cohen_d) * ifelse(temp_wilcox <= 0.05,3,ifelse(temp_wilcox <= 0.10,2,1))
      }


StehlikovaZ_metadata <-MetadataDf_buccal_palate_DisCtrl[MetadataDf_buccal_palate_DisCtrl$study_name == "StehlikovaZ_2019",]
StehlikovaZ_SpDf <- SpDf_buccal_palate_DisCtrl[rownames(StehlikovaZ_metadata),]
StehlikovaZ_control_samples <- rownames(StehlikovaZ_metadata[StehlikovaZ_metadata$study_condition == "Control",])
StehlikovaZ_disease_samples <- rownames(StehlikovaZ_metadata[StehlikovaZ_metadata$study_condition != "Control",])

for(species_name in buccal_AssociatedSpecies)
      {
        control_values <- as.numeric(StehlikovaZ_SpDf[StehlikovaZ_control_samples,species_name])
        disease_values <-as.numeric(StehlikovaZ_SpDf[StehlikovaZ_disease_samples,species_name])

        temp_cohen_d <- cohen.d(control_values,disease_values)$estimate
        temp_wilcox <- wilcox.test(control_values,disease_values)$p.value
        df_comparison[species_name,2] <- sign(temp_cohen_d) * ifelse(temp_wilcox <= 0.05,3,ifelse(temp_wilcox <= 0.10,2,1))
      }


df_comparison[is.na(df_comparison)] <- 0
df_comparison_matrix <- as.matrix(df_comparison)

df_comparison_matrix[df_comparison_matrix %in% c(-1, 0, 1)] <- 0

df_comparison <- as.data.frame(df_comparison_matrix)

df_comparison <- df_comparison[!apply(df_comparison, 1, function(x) all(x == 0, na.rm = TRUE)),,drop = FALSE]


color_map <- c("-3"  = "#8B0000","-2"  = "#CD5C5C","-1"  = "white","0"  = "white","1" = "white","2" = "#329732","3" = "#006400")

library(pheatmap)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1BuccalPalate_heatmap_StudyWise_HealthSpecies.pdf", width = 40, height = 20)
pheatmap(t(df_comparison),
         color = color_map,
         fontsize_row = 8,
         fontsize_col = 8,
         cellheight = 24,
         cellwidth = 18,
         cluster_rows = T,
         cluster_cols = T,
         border_color = "black",
         treeheight_row = 0,
         treeheight_col = 0
)
dev.off()


########## Get the carpet
ph <- pheatmap(
  t(df_comparison),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  silent = TRUE
)

## Rows of t(df_comparison) = columns of df_comparison
ordered_cols <- rownames(t(df_comparison))[ph$tree_row$order]

## Columns of t(df_comparison) = rows of df_comparison
ordered_rows <- colnames(t(df_comparison))[ph$tree_col$order]

## Reorder the original data frame
df_comparison_sorted <- df_comparison[ordered_rows,ordered_cols,drop = FALSE]

write.csv(df_comparison_sorted, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1BuccalPalate_heatmap_carpet_StudyWise_HealthSpecies.csv")


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Buccal_palate_HealthAssociation_Workspace.RData")
