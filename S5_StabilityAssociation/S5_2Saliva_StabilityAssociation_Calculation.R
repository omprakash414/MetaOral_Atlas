

######## Load the distance matrix calculated earlier and other needed data to calculate the stability association score
S5_1Followup_Distance_calculation_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_1Followup_Distance_calculation_Workspace.RData", envir = S5_1Followup_Distance_calculation_Workspace)
attach(S5_1Followup_Distance_calculation_Workspace)
SpDf_saliva_Long_withDist <- SpDf_saliva_Long_withDist
saliva_AssociatedSpecies <- saliva_AssociatedSpecies
std_meta_cols <- std_meta_cols
detach(S5_1Followup_Distance_calculation_Workspace)

## Load the functions:
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")

######## Prepare the data
# It is already done in the previous script S5_1Followup_Distance_Calculation.R



######## Calculate the stability scores using both distance metrics
## Aitchison
aitch_iters <- Meta_lm_Iterrative(data = SpDf_saliva_Long_withDist,species = saliva_AssociatedSpecies,dist_var = "Aitchison_dist",study_var = "Study_Name", n_studies = 8)
aitch_res <- Compute_Meta_Stability(iter_list = aitch_iters,species = saliva_AssociatedSpecies)

## Bray–Curtis
bray_iters <- Meta_lm_Iterrative(data = SpDf_saliva_Long_withDist,species = saliva_AssociatedSpecies,dist_var = "BrayCurtis_dist",study_var = "Study_Name", n_studies = 8)
bray_res <- Compute_Meta_Stability(iter_list = bray_iters,species = saliva_AssociatedSpecies)


## Combine the stability scores from both distance metrics
stabilityRank <- cbind(aitch_res$stability,bray_res$stability)
colnames(stabilityRank) <- c("mean_Aitchison","IQR_Aitchison","mean_BrayCurtis","IQR_BrayCurtis")


stabilityRank$MeanStabilityScore <- rowMeans(stabilityRank[, c(1,3)])
stabilityRank <- stabilityRank[order(-stabilityRank$MeanStabilityScore), ]

write.csv(stabilityRank, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2Summarised_StabilityScores.csv")
#write.csv(stabilityRank, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/temp_S5_2Summarised_StabilityScores.csv")

Saliva_StabilityScore <- data.frame(species = rownames(stabilityRank),StabilityScore = stabilityRank$MeanStabilityScore)
rownames(Saliva_StabilityScore) <- Saliva_StabilityScore$species

# save
save(Saliva_StabilityScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2Saliva_StabilityAssociationScore.RData")

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2Saliva_StabilityAssociation_Calculation_Workspace.RData")



###### Export some of the files. 
bray_summarised_EffectSize <- bray_res$beta_summary
bray_summarised_Qvalue <- bray_res$qval_summary

aitch_summarised_EffectSize <- aitch_res$beta_summary
aitch_summarised_Qvalue <- aitch_res$qval_summary















##################
########################### To get the study specific stability association we need to run the same pipeline but for each study separately.
## Function to calculate the study specific stability association for each species
compute_study_species_stability <- function(data, species_list, dist_var, study_var, p_cutoff = 0.05) {
  
  studies <- unique(data[[study_var]])
  
  stability_mat <- matrix(
    0,
    nrow = length(species_list),
    ncol = length(studies)
  )
  
  rownames(stability_mat) <- species_list
  colnames(stability_mat) <- studies
  
  pval_mat <- stability_mat
  beta_mat <- stability_mat
  
  for (sp in species_list) {
    print(sp)
    
    for (study in studies) {
      
      temp_data <- data[data[[study_var]] == study, ]
      
      vec_sp <- temp_data[[sp]]
      vec_dist <- temp_data[[dist_var]]
      
      if (sum(abs(vec_sp) > 0, na.rm = TRUE) > 0 &&
          sum(abs(vec_dist) > 0, na.rm = TRUE) > 0) {
        
        tryCatch({
          
          f <- as.formula(paste0("`", sp, "` ~ ", dist_var))
          temp_rlm <- MASS::rlm(f, data = temp_data)
          
          beta <- coef(temp_rlm)[2]
          pval <- sfsmisc::f.robftest(temp_rlm, var = dist_var)$p.value
          
          beta_mat[sp, study] <- beta
          pval_mat[sp, study] <- pval
          
          ## Stable = negative association with distance and significant
          if (beta < 0 && pval <= p_cutoff) {
            stability_mat[sp, study] <- 1
          }
          
        }, error = function(e) {
          beta_mat[sp, study] <- NA
          pval_mat[sp, study] <- NA
          stability_mat[sp, study] <- 0
        })
      }
    }
  }
  
  return(list(
    stability_mat = stability_mat,
    beta_mat = beta_mat,
    pval_mat = pval_mat
  ))
}




aitch_study_stability <- compute_study_species_stability(
  data = SpDf_saliva_Long_withDist,
  species_list = saliva_AssociatedSpecies,
  dist_var = "Aitchison_dist",
  study_var = "Study_Name",
  p_cutoff = 0.10
)

bray_study_stability <- compute_study_species_stability(
  data = SpDf_saliva_Long_withDist,
  species_list = saliva_AssociatedSpecies,
  dist_var = "BrayCurtis_dist",
  study_var = "Study_Name",
  p_cutoff = 0.10
)

combined_stability_mat <- data.frame(ifelse(aitch_study_stability$stability_mat == 1 | bray_study_stability$stability_mat == 1, 1, 0))

combined_stability_mat$study_count <- rowSums(combined_stability_mat)

combined_stability_mat <- combined_stability_mat[order(-combined_stability_mat$study_count), ]

## filter only those species that are present in at least 2 studies
combined_stability_mat_filtered <- combined_stability_mat[combined_stability_mat$study_count >= 2, ]
combined_stability_mat_filtered <- combined_stability_mat_filtered[,colSums(combined_stability_mat_filtered) > 0]

## save as csv
write.csv(combined_stability_mat_filtered, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2StudyWise_StabilityAssociations.csv")


## plot the heatmap of the stability association across studies
library(pheatmap)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2Saliva_heatmap_StudyWise_StabilityAssociations.pdf", width = 30, height = 15)
pheatmap(t(combined_stability_mat_filtered[,1:12]),
         color = c("white", "#2C7FB8"),
         fontsize_row = 12,
         fontsize_col = 12,
         cellheight = 16,
         cellwidth = 16,
         cluster_rows = T,
         cluster_cols = F,
         border_color = "black",
         treeheight_row = 0,
         treeheight_col = 0
)
dev.off()


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2Saliva_StabilityAssociation_Calculation_Workspace.RData")


## now plot the line plot for the total studies in which each species is stable
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2Saliva_lineplot_StudyWise_StabilityAssociations.pdf", width = 30, height = 5)
library(ggplot2)
ggplot(
  combined_stability_mat_filtered[, c("study_count"), drop = F],
  aes(
    x = factor(rownames(combined_stability_mat_filtered), levels = rownames(combined_stability_mat_filtered)),
    y = study_count,
    group = 1
  )
) +
  geom_line() +
  geom_point(size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    x = "Species",
    y = "Number of studies"
  )
dev.off()

