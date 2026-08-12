

########### Load the required speies profile and selected species for each subsite for core detection
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Species_Profiles_Associated_Species_SubsiteWise.RData")


########### Load the metadata for all subsites
S1_Threshold_Species_Detection_Workspace <- new.env()
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Threshold_Species_Detection_Workspace.RData", envir = S1_Threshold_Species_Detection_Workspace)
attach(S1_Threshold_Species_Detection_Workspace)

metadata_objs <- ls(envir = S1_Threshold_Species_Detection_Workspace, pattern = "^MetadataDf")
list2env(mget(metadata_objs, envir = S1_Threshold_Species_Detection_Workspace),envir = .GlobalEnv)

detach(S1_Threshold_Species_Detection_Workspace)
rm(S1_Threshold_Species_Detection_Workspace)
rm(metadata_objs)


rm(list = ls(pattern = "supragingival"))
rm(list = ls(pattern = "subgingival"))
rm(list = ls(pattern = "tongue"))
rm(list = ls(pattern = "buccal"))
########### Load all the fucntions required for the core detection.
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")



########### Saliva Core (0.85)
# Check all row are in same order
all(rownames(MetadataDf_saliva_control) == rownames(SpDf_saliva_Control_associated))
SpDf_saliva_Control_associated <- SpDf_saliva_Control_associated/rowSums(SpDf_saliva_Control_associated)
SpDf_saliva_Control_associated$study_name <- MetadataDf_saliva_control$study_name

study_list <- unique(SpDf_saliva_Control_associated$study_name)



saliva_Core_Influencers <- keystoneInfluence(SpDf_saliva_Control_associated, saliva_AssociatedSpecies, 0.85)

## Summarising r2 df and pr df (pval df)
Saliva_summary_r2PrDf <- Summarise_r2pval_CoreKeystone(saliva_Core_Influencers$ENV_fit_summary, rank_scale)
saliva_r2_rankedDf <- Saliva_summary_r2PrDf$r2df_ranked
saliva_r2_unrankedDf <- Saliva_summary_r2PrDf$r2df_unranked
saliva_pval_Df <- Saliva_summary_r2PrDf$pvalue_df

write.csv(saliva_r2_rankedDf, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_R2Df_coreAssociation.csv")
write.csv(saliva_pval_Df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_PvalDf_coreAssociation.csv")

########### Calculate the prevalence of the species in each of the study (same as done in S0 but this time selected species)
SpDf_saliva_Control_associated$study_name <- NULL

## Remove the zero rowsum rows and then normalize the data.
SpDf_saliva_Control_associated <- SpDf_saliva_Control_associated[rowSums(SpDf_saliva_Control_associated)>0,]
SpDf_saliva_Control_associated <- SpDf_saliva_Control_associated/rowSums(SpDf_saliva_Control_associated)

MetadataDf_saliva_control <- MetadataDf_saliva_control[rownames(SpDf_saliva_Control_associated),]

SpDf_saliva_Control_associated$study_name <- MetadataDf_saliva_control$study_name

saliva_prevalDf <- compute_detection(SpDf_saliva_Control_associated,saliva_AssociatedSpecies,"study_name",study_list)
write.csv(saliva_prevalDf, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_PrevalenceDf_coreAssociation.csv")

########### Determining the threshold of r2 value to further define the core microbes and plot it.

r2_threshold_df <- compute_EnvFit_r2_threshold(saliva_prevalDf,saliva_r2_rankedDf,saliva_pval_Df,study_list,0.85)

## Take the mean of accuracies for each threshold across all the studies. and then plot it to see which threshold has max mean which is max accuracy which is r2 threshold to use
r2_threshold_represented  <- data.frame(Threshold = as.numeric(rownames(r2_threshold_df)), Accuracy  = rowMeans(r2_threshold_df, na.rm = TRUE))

r2_threshold_represented[which.max(r2_threshold_represented[, 2]), 1]
# its 0.75

## Plot it
library(ggplot2)
r2_threshold_represented$Threshold_f <- factor(r2_threshold_represented$Threshold)
pdf(
  "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_r2_threshold_detection.pdf",
  width = 6,
  height = 4
)

ggplot(
  r2_threshold_represented,
  aes(x = factor(Threshold), y = Accuracy)
) +
  geom_bar(
    stat = "identity",
    fill = "#4C94B0",
    width = 0.8
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = expression("Envfit " * R^2 * " rank threshold"),
    y = "Mean accuracy across studies"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 13)
  )

dev.off()


########### Now calculate the Core Association Score using prevalance threshold as 0.85 and r2 threshold as 0.75
saliva_CoreKeyStoneDf <- data.frame(apply(saliva_r2_rankedDf,2,function(x)(ifelse(x>=0.75,1,0))) * apply(saliva_prevalDf,2,function(x)(ifelse(x>=0.85,1,0))))

write.csv(saliva_CoreKeyStoneDf, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_StudySpecific_CoreSpecies.csv")

############ Get the saliva CoreKeyStones and their study wise distribution

saliva_Core_Representation <- summarize_core_keystone_detection(saliva_CoreKeyStoneDf)

saliva_CoreKeyStoneDf2 <- saliva_Core_Representation$core_keystone_by_species
saliva_CoreKeyStoneDf_cumulative <- saliva_Core_Representation$detection_summary


library(ggplot2)
## Plot these cumulative pattern
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_cumulative_CoreKeyStone_detection_pattern.pdf",width = 15, height = 8)

# Create the dot plot
ggplot(saliva_CoreKeyStoneDf_cumulative, aes(x = as.factor(studies_detected), y = cumulative_n_species)) +
  geom_bar(stat = "identity", fill = "mediumpurple2") +
  labs(x = "N Studies",
       y = "Cumulative CoreKeyStone (species)") +
  scale_y_continuous(breaks = seq(0, max(saliva_CoreKeyStoneDf_cumulative$cumulative_n_species), by = 10)) +
  geom_text(aes(y = cumulative_n_species + 10, label = cumulative_n_species), vjust = 0.3, angle = 90) +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

dev.off()


############ Calcualte the Core Association Score
saliva_CoreKeyStoneDf2$CoreAssociationScore <- rank_scale((saliva_CoreKeyStoneDf2$studies_detected)/length(study_list))

saliva_CoreAssociationScore <- saliva_CoreKeyStoneDf2[,tail(colnames(saliva_CoreKeyStoneDf2), 2)]
saliva_CoreAssociationScore <- saliva_CoreAssociationScore[order(saliva_CoreAssociationScore[, 2], decreasing = TRUE), ]

save(saliva_CoreAssociationScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_CoreAssociationScore.RData")

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_Core_Detection_Workspace.RData")



# ##### From the cumulative plot of species, top 17 species (50 studies) can be considered as Major core species, while 67 species (24 studies) can be considered as minor core
# #top 17 species 
# saliva_CoreAssociationScore$species[1:17]
# #  [1] "Prevotella_salivae"          "Atopobium_parvulum"         
# #  [3] "Prevotella_pallens"          "Fusobacterium_periodonticum"
# #  [5] "Megasphaera_micronuciformis" "Prevotella_melaninogenica"  
# #  [7] "Porphyromonas_catoniae"      "Veillonella_dispar"         
# #  [9] "Prevotella_histicola"        "Haemophilus_parainfluenzae" 
# # [11] "Fusobacterium_nucleatum"     "Streptococcus_thermophilus" 
# # [13] "Actinomyces_graevenitzii"    "Campylobacter_concisus"     
# # [15] "Actinomyces_odontolyticus"   "Rothia_mucilaginosa"        
# # [17] "Prevotella_nanceiensis" 

# # top 67 species
# saliva_CoreAssociationScore$species[1:67]
# saliva_CoreKeyStoneDf3 <- saliva_CoreKeyStoneDf2[saliva_CoreKeyStoneDf2$studies_detected >= 24,c("studies_detected","species")]
# saliva_CoreKeyStoneDf3$Index <- rev(seq_len(nrow(saliva_CoreKeyStoneDf3)))

# library(ggplot2)
# library(ggrepel)


# pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_MinorCore_representation_in_total_no.studies.pdf",
#   height = 7, width = 18)

# ggplot(saliva_CoreKeyStoneDf3[saliva_CoreKeyStoneDf3$studies_detected >= 24, ],
#   aes(x = Index, y = studies_detected)
# ) +
#   geom_point(
#     color = "#2C7FB8",   # deep blue (clean & elegant)
#     size = 3.5,
#     alpha = 0.9
#   ) +
#   geom_text_repel(
#     aes(
#       label = rownames(saliva_CoreKeyStoneDf3)
#     ),
#     max.overlaps = 20,
#     size = 5,
#     box.padding = 0.15,
#     min.segment.length = 1
#   ) +
#   theme_bw(base_size = 14)

# dev.off()



# save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_Core_Detection_Workspace.RData")



################## Plot the Study vs core species (1/0) as a heatmap. This will show the distribution of core species across studies.

saliva_CoreKeyStoneDf3 <- saliva_CoreKeyStoneDf2[,1:99]
saliva_CoreKeyStoneDf3 <- saliva_CoreKeyStoneDf3[rowSums(saliva_CoreKeyStoneDf3)>19,]
saliva_CoreKeyStoneDf3 <- saliva_CoreKeyStoneDf3[,colSums(saliva_CoreKeyStoneDf3)>0]

library(pheatmap)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_heatmap_StudyWise_CoreSpecies.pdf", width = 40, height = 20)
pheatmap(saliva_CoreKeyStoneDf3,
         color = c("white", "#2C7FB8"),
         fontsize_row = 8,
         fontsize_col = 8,
         cellheight = 10,
         cellwidth = 16,
         cluster_rows = F,
         cluster_cols = T,
         border_color = "black",
         treeheight_row = 0,
         treeheight_col = 0
)
dev.off()

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_Core_Detection_Workspace.RData")

########## Get the carpet
ph <- pheatmap(
  saliva_CoreKeyStoneDf3,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  silent = TRUE
)

ordered_cols <- colnames(saliva_CoreKeyStoneDf3)[ph$tree_col$order]
saliva_CoreKeyStoneDf4 <- saliva_CoreKeyStoneDf3[, ordered_cols]

write.csv(saliva_CoreKeyStoneDf4, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_heatmap_carpet_StudyWise_CoreSpecies.csv")

################### Get the line plot for the core species distribution across studies. This will show the number of studies in which the species is detected as core species.
collective_CoreSpecies_distribution <- data.frame(species = rownames(saliva_CoreKeyStoneDf3), number_of_studies = rowSums(saliva_CoreKeyStoneDf3))

all(rownames(collective_CoreSpecies_distribution) == rownames(saliva_CoreKeyStoneDf3))
# its true. Now plot the line plot.

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_lineplot_CoreSpecies.pdf", width = 40, height = 10)
library(ggplot2)
ggplot(
  collective_CoreSpecies_distribution,
  aes(
    x = factor(species, levels = species),
    y = number_of_studies,
    group = 1
  )
) +
  geom_line() +
  geom_point(size = 8) +
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


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Saliva_Core_Detection_Workspace.RData")
