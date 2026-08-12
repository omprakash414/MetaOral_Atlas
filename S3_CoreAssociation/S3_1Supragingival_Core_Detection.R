

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


rm(list = ls(pattern = "saliva"))
rm(list = ls(pattern = "subgingival"))
rm(list = ls(pattern = "tongue"))
rm(list = ls(pattern = "buccal"))

########### Load all the fucntions required for the core detection.
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")

study_list <- unique(SpDf_supragingival_Control_associated$study_name)

########## Supragingival (0.80)
all(rownames(MetadataDf_supragingival_control) == rownames(SpDf_supragingival_Control_associated))
SpDf_supragingival_Control_associated$study_name <- MetadataDf_supragingival_control$study_name

supragingival_Core_Influencers <- keystoneInfluence(SpDf_supragingival_Control_associated, supragingival_AssociatedSpecies, 0.80)

## Summarising r2 df and pr df (pval df)
supragingival_summary_r2PrDf <- Summarise_r2pval_CoreKeystone(supragingival_Core_Influencers$ENV_fit_summary, rank_scale)
supragingival_r2_rankedDf <- supragingival_summary_r2PrDf$r2df_ranked
supragingival_r2_unrankedDf <- supragingival_summary_r2PrDf$r2df_unranked
supragingival_pval_Df <- supragingival_summary_r2PrDf$pvalue_df



########### Calculate the prevalence of the species in each of the study (same as done in S0 but this time selected species)
supragingival_prevalDf <- compute_detection(SpDf_supragingival_Control_associated,supragingival_AssociatedSpecies,"study_name",study_list)

########### Determining the threshold of r2 value to further define the core microbes and plot it.

supragingival_r2_threshold_df <- compute_EnvFit_r2_threshold(supragingival_prevalDf,supragingival_r2_rankedDf,supragingival_pval_Df,study_list,0.80)

## Take the mean of accuracies for each threshold across all the studies. and then plot it to see which threshold has max mean which is max accuracy which is r2 threshold to use
supragingival_r2_threshold_represented  <- data.frame(Threshold = as.numeric(rownames(supragingival_r2_threshold_df)), Accuracy  = rowMeans(supragingival_r2_threshold_df))

supragingival_r2_threshold_represented[which.max(supragingival_r2_threshold_represented[, 2]), 1]
# its 0.75

## Plot it
library(ggplot2)
supragingival_r2_threshold_represented$Threshold_f <- factor(supragingival_r2_threshold_represented$Threshold)
pdf(
  "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1supragingival_r2_threshold_detection.pdf",
  width = 6,
  height = 4
)

ggplot(
  supragingival_r2_threshold_represented,
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
supragingival_CoreKeyStoneDf <- data.frame(apply(supragingival_r2_rankedDf,2,function(x)(ifelse(x>=0.75,1,0))) * apply(supragingival_prevalDf,2,function(x)(ifelse(x>=0.80,1,0))))

write.csv(supragingival_CoreKeyStoneDf, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1supragingival_StudySpecific_CoreSpecies.csv")

############ Get the supragingival CoreKeyStones and their study wise distribution

supragingival_Core_Representation <- summarize_core_keystone_detection(supragingival_CoreKeyStoneDf)

supragingival_CoreKeyStoneDf2 <- supragingival_Core_Representation$core_keystone_by_species
supragingival_CoreKeyStoneDf_cumulative <- supragingival_Core_Representation$detection_summary


library(ggplot2)
## Plot these cumulative pattern
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1supragingival_cumulative_CoreKeyStone_detection_pattern.pdf",width = 15, height = 8)

# Create the dot plot
ggplot(supragingival_CoreKeyStoneDf_cumulative, aes(x = as.factor(studies_detected), y = cumulative_n_species)) +
  geom_bar(stat = "identity", fill = "mediumpurple2") +
  labs(x = "N Studies",
       y = "Cumulative CoreKeyStone (species)") +
  scale_y_continuous(breaks = seq(0, max(supragingival_CoreKeyStoneDf_cumulative$cumulative_n_species), by = 10)) +
  geom_text(aes(y = cumulative_n_species + 10, label = cumulative_n_species), vjust = 0.3, angle = 90) +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

dev.off()



############ Calcualte the Core Association Score
supragingival_CoreKeyStoneDf2$CoreAssociationScore <- rank_scale((supragingival_CoreKeyStoneDf2$studies_detected)/length(study_list))

supragingival_CoreAssociationScore <- supragingival_CoreKeyStoneDf2[,tail(colnames(supragingival_CoreKeyStoneDf2), 2)]
supragingival_CoreAssociationScore <- supragingival_CoreAssociationScore[order(supragingival_CoreAssociationScore[, 2], decreasing = TRUE), ]

save(supragingival_CoreAssociationScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1supragingival_CoreAssociationScore.RData")




save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Supragingival_Core_Detection_Workspace.RData")





################## Plot the Study vs core species (1/0) as a heatmap. This will show the distribution of core species across studies.

supragingival_CoreKeyStoneDf3 <- supragingival_CoreKeyStoneDf2[,1:23]
supragingival_CoreKeyStoneDf3 <- supragingival_CoreKeyStoneDf3[,colSums(supragingival_CoreKeyStoneDf3)>0]
supragingival_CoreKeyStoneDf3 <- supragingival_CoreKeyStoneDf3[rowSums(supragingival_CoreKeyStoneDf3)>0,]


library(pheatmap)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Supragingival_heatmap_StudyWise_CoreSpecies.pdf", width = 40, height = 20)
pheatmap(t(supragingival_CoreKeyStoneDf3),
         color = c("white", "#2C7FB8"),
         fontsize_row = 8,
         fontsize_col = 8,
         cellheight = 20,
         cellwidth = 16,
         cluster_rows = T,
         cluster_cols = F,
         border_color = "black",
         treeheight_row = 0,
         treeheight_col = 0
)
dev.off()


########## Get the carpet
ph <- pheatmap(
  t(supragingival_CoreKeyStoneDf3),
  cluster_rows = T,
  cluster_cols = F,
  silent = TRUE
)

ordered_cols <- rownames(t(supragingival_CoreKeyStoneDf3))[ph$tree_row$order]
supragingival_CoreKeyStoneDf4 <- supragingival_CoreKeyStoneDf3[,ordered_cols ]

write.csv(supragingival_CoreKeyStoneDf4, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Supragingival_heatmap_carpet_StudyWise_CoreSpecies.csv")

################### Get the line plot for the core species distribution across studies. This will show the number of studies in which the species is detected as core species.
collective_CoreSpecies_distribution <- data.frame(species = rownames(supragingival_CoreKeyStoneDf3), number_of_studies = rowSums(supragingival_CoreKeyStoneDf3))

all(rownames(collective_CoreSpecies_distribution) == rownames(supragingival_CoreKeyStoneDf3))
# its true. Now plot the line plot.

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Supragingival_lineplot_CoreSpecies.pdf", width = 40, height = 10)
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

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1Supragingival_Core_Detection_Workspace.RData")
