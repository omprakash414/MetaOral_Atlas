

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


rm(list = ls(pattern = "subgingival"))
rm(list = ls(pattern = "subgingival"))
rm(list = ls(pattern = "saliva"))
rm(list = ls(pattern = "buccal"))
########### Load all the fucntions required for the core detection.
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")


study_list <- unique(SpDf_tongue_tonsil_Control_associated$study_name)

########## tonguetonsil (0.80)
all(rownames(MetadataDf_tongue_tonsil_control) == rownames(SpDf_tongue_tonsil_Control_associated))
SpDf_tongue_tonsil_Control_associated$study_name <- MetadataDf_tongue_tonsil_control$study_name

tongue_tonsil_Core_Influencers <- keystoneInfluence(SpDf_tongue_tonsil_Control_associated, tonguetonsil_AssociatedSpecies, 0.80)


## Summarising r2 df and pr df (pval df)
tongue_tonsil_summary_r2PrDf <- Summarise_r2pval_CoreKeystone(tongue_tonsil_Core_Influencers$ENV_fit_summary, rank_scale)
tongue_tonsil_r2_rankedDf <- tongue_tonsil_summary_r2PrDf$r2df_ranked
tongue_tonsil_r2_unrankedDf <- tongue_tonsil_summary_r2PrDf$r2df_unranked
tongue_tonsil_pval_Df <- tongue_tonsil_summary_r2PrDf$pvalue_df


########### Calculate the prevalence of the species in each of the study (same as done in S0 but this time selected species)
tongue_tonsil_prevalDf <- compute_detection(SpDf_tongue_tonsil_Control_associated,tonguetonsil_AssociatedSpecies,"study_name",study_list)

########### Determining the threshold of r2 value to further define the core microbes and plot it.

tongue_tonsil_r2_threshold_df <- compute_EnvFit_r2_threshold(tongue_tonsil_prevalDf,tongue_tonsil_r2_rankedDf,tongue_tonsil_pval_Df,study_list,0.80)

## Take the mean of accuracies for each threshold across all the studies. and then plot it to see which threshold has max mean which is max accuracy which is r2 threshold to use
tongue_tonsil_r2_threshold_represented  <- data.frame(Threshold = as.numeric(rownames(tongue_tonsil_r2_threshold_df)), Accuracy  = rowMeans(tongue_tonsil_r2_threshold_df))

tongue_tonsil_r2_threshold_represented[which.max(tongue_tonsil_r2_threshold_represented[, 2]), 1]
# its 0.75

## Plot it
library(ggplot2)
tongue_tonsil_r2_threshold_represented$Threshold_f <- factor(tongue_tonsil_r2_threshold_represented$Threshold)
pdf(
  "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1tongue_tonsil_r2_threshold_detection.pdf",
  width = 6,
  height = 4
)

ggplot(
  tongue_tonsil_r2_threshold_represented,
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



########### Now calculate the Core Association Score using prevalance threshold as 0.80 and r2 threshold as 0.75
tongue_tonsil_CoreKeyStoneDf <- data.frame(apply(tongue_tonsil_r2_rankedDf,2,function(x)(ifelse(x>=0.75,1,0))) * apply(tongue_tonsil_prevalDf,2,function(x)(ifelse(x>=0.80,1,0))))

write.csv(tongue_tonsil_CoreKeyStoneDf, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1tonguetonsil_StudySpecific_CoreSpecies.csv")

############ Get the tongue_tonsil CoreKeyStones and their study wise distribution

tongue_tonsil_Core_Representation <- summarize_core_keystone_detection(tongue_tonsil_CoreKeyStoneDf)

tongue_tonsil_CoreKeyStoneDf2 <- tongue_tonsil_Core_Representation$core_keystone_by_species
tongue_tonsil_CoreKeyStoneDf_cumulative <- tongue_tonsil_Core_Representation$detection_summary


library(ggplot2)
## Plot these cumulative pattern
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1tongue_tonsil_cumulative_CoreKeyStone_detection_pattern.pdf",width = 15, height = 8)

# Create the dot plot
ggplot(tongue_tonsil_CoreKeyStoneDf_cumulative, aes(x = as.factor(studies_detected), y = cumulative_n_species)) +
  geom_bar(stat = "identity", fill = "mediumpurple2") +
  labs(x = "N Studies",
       y = "Cumulative CoreKeyStone (species)") +
  scale_y_continuous(breaks = seq(0, max(tongue_tonsil_CoreKeyStoneDf_cumulative$cumulative_n_species), by = 10)) +
  geom_text(aes(y = cumulative_n_species + 10, label = cumulative_n_species), vjust = 0.3, angle = 90) +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

dev.off()


############ Calcualte the Core Association Score
tongue_tonsil_CoreKeyStoneDf2$CoreAssociationScore <- rank_scale((tongue_tonsil_CoreKeyStoneDf2$studies_detected)/length(study_list))

tongue_tonsil_CoreAssociationScore <- tongue_tonsil_CoreKeyStoneDf2[,tail(colnames(tongue_tonsil_CoreKeyStoneDf2), 2)]
tongue_tonsil_CoreAssociationScore <- tongue_tonsil_CoreAssociationScore[order(tongue_tonsil_CoreAssociationScore[, 2], decreasing = TRUE), ]

save(tongue_tonsil_CoreAssociationScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1tongue_tonsil_CoreAssociationScore.RData")


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1TongueTonsil_Core_Detection_Workspace.RData")







################## Plot the Study vs core species (1/0) as a heatmap. This will show the distribution of core species across studies.

tongue_tonsil_CoreKeyStoneDf3 <- tongue_tonsil_CoreKeyStoneDf2[,1:18]
tongue_tonsil_CoreKeyStoneDf3 <- tongue_tonsil_CoreKeyStoneDf3[,colSums(tongue_tonsil_CoreKeyStoneDf3)>0]
tongue_tonsil_CoreKeyStoneDf3 <- tongue_tonsil_CoreKeyStoneDf3[rowSums(tongue_tonsil_CoreKeyStoneDf3)>0,]


library(pheatmap)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1TongueTonsil_heatmap_StudyWise_CoreSpecies.pdf", width = 40, height = 20)
pheatmap(t(tongue_tonsil_CoreKeyStoneDf3),
         color = c("white", "#2C7FB8"),
         fontsize_row = 8,
         fontsize_col = 8,
         cellheight = 20,
         cellwidth = 18,
         cluster_rows = T,
         cluster_cols = F,
         border_color = "black",
         treeheight_row = 0,
         treeheight_col = 0
)
dev.off()


########## Get the carpet
ph <- pheatmap(
  t(tongue_tonsil_CoreKeyStoneDf3),
  cluster_rows = T,
  cluster_cols = F,
  silent = TRUE
)

ordered_cols <- rownames(t(tongue_tonsil_CoreKeyStoneDf3))[ph$tree_row$order]
tongue_tonsil_CoreKeyStoneDf4 <- tongue_tonsil_CoreKeyStoneDf3[,ordered_cols ]

write.csv(tongue_tonsil_CoreKeyStoneDf4, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1TongueTonsil_heatmap_carpet_StudyWise_CoreSpecies.csv")

################### Get the line plot for the core species distribution across studies. This will show the number of studies in which the species is detected as core species.
collective_CoreSpecies_distribution <- data.frame(species = rownames(tongue_tonsil_CoreKeyStoneDf3), number_of_studies = rowSums(tongue_tonsil_CoreKeyStoneDf3))

all(rownames(collective_CoreSpecies_distribution) == rownames(tongue_tonsil_CoreKeyStoneDf3))
# its true. Now plot the line plot.

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1TongueTonsil_lineplot_CoreSpecies.pdf", width = 40, height = 5)
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


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1TongueTonsil_Core_Detection_Workspace.RData")


