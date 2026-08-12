

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
rm(list = ls(pattern = "supragingival"))
rm(list = ls(pattern = "tongue"))
rm(list = ls(pattern = "saliva"))
########### Load all the fucntions required for the core detection.
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")

study_list <- unique(SpDf_buccal_palate_Control_associated$study_name)

########## buccal (0.75)
all(rownames(MetadataDf_buccal_palate_control) == rownames(SpDf_buccal_palate_Control_associated))
SpDf_buccal_palate_Control_associated$study_name <- MetadataDf_buccal_palate_control$study_name

buccal_palate_Core_Influencers <- keystoneInfluence(SpDf_buccal_palate_Control_associated, buccal_AssociatedSpecies, 0.75)


## Summarising r2 df and pr df (pval df)
buccal_palate_summary_r2PrDf <- Summarise_r2pval_CoreKeystone(buccal_palate_Core_Influencers$ENV_fit_summary, rank_scale)
buccal_palate_r2_rankedDf <- buccal_palate_summary_r2PrDf$r2df_ranked
buccal_palate_r2_unrankedDf <- buccal_palate_summary_r2PrDf$r2df_unranked
buccal_palate_pval_Df <- buccal_palate_summary_r2PrDf$pvalue_df


########### Calculate the prevalence of the species in each of the study (same as done in S0 but this time selected species)
buccal_palate_prevalDf <- compute_detection(SpDf_buccal_palate_Control_associated,buccal_AssociatedSpecies,"study_name",study_list)

########### Determining the threshold of r2 value to further define the core microbes and plot it.

buccal_palate_r2_threshold_df <- compute_EnvFit_r2_threshold(buccal_palate_prevalDf,buccal_palate_r2_rankedDf,buccal_palate_pval_Df,study_list,0.75)

## Take the mean of accuracies for each threshold across all the studies. and then plot it to see which threshold has max mean which is max accuracy which is r2 threshold to use
buccal_palate_r2_threshold_represented  <- data.frame(Threshold = as.numeric(rownames(buccal_palate_r2_threshold_df)), Accuracy  = rowMeans(buccal_palate_r2_threshold_df))

buccal_palate_r2_threshold_represented[which.max(buccal_palate_r2_threshold_represented[, 2]), 1]
# its 0.75

## Plot it
library(ggplot2)
buccal_palate_r2_threshold_represented$Threshold_f <- factor(buccal_palate_r2_threshold_represented$Threshold)
pdf(
  "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1buccal_palate_r2_threshold_detection.pdf",
  width = 6,
  height = 4
)

ggplot(
  buccal_palate_r2_threshold_represented,
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



########### Now calculate the Core Association Score using prevalance threshold as 0.75 and r2 threshold as 0.75
buccal_palate_CoreKeyStoneDf <- data.frame(apply(buccal_palate_r2_rankedDf,2,function(x)(ifelse(x>=0.75,1,0))) * apply(buccal_palate_prevalDf,2,function(x)(ifelse(x>=0.75,1,0))))

write.csv(buccal_palate_CoreKeyStoneDf, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1buccalpalate_StudySpecific_CoreSpecies.csv")

############ Get the buccal_palate CoreKeyStones and their study wise distribution

buccal_palate_Core_Representation <- summarize_core_keystone_detection(buccal_palate_CoreKeyStoneDf)

buccal_palate_CoreKeyStoneDf2 <- buccal_palate_Core_Representation$core_keystone_by_species
buccal_palate_CoreKeyStoneDf_cumulative <- buccal_palate_Core_Representation$detection_summary


library(ggplot2)
## Plot these cumulative pattern
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1buccal_palate_cumulative_CoreKeyStone_detection_pattern.pdf",width = 15, height = 8)

# Create the dot plot
ggplot(buccal_palate_CoreKeyStoneDf_cumulative, aes(x = as.factor(studies_detected), y = cumulative_n_species)) +
  geom_bar(stat = "identity", fill = "mediumpurple2") +
  labs(x = "N Studies",
       y = "Cumulative CoreKeyStone (species)") +
  scale_y_continuous(breaks = seq(0, max(buccal_palate_CoreKeyStoneDf_cumulative$cumulative_n_species), by = 10)) +
  geom_text(aes(y = cumulative_n_species + 10, label = cumulative_n_species), vjust = 0.3, angle = 90) +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

dev.off()



############ Calcualte the Core Association Score
buccal_palate_CoreKeyStoneDf2$CoreAssociationScore <- rank_scale((buccal_palate_CoreKeyStoneDf2$studies_detected)/length(study_list))

buccal_palate_CoreAssociationScore <- buccal_palate_CoreKeyStoneDf2[,tail(colnames(buccal_palate_CoreKeyStoneDf2), 2)]
buccal_palate_CoreAssociationScore <- buccal_palate_CoreAssociationScore[order(buccal_palate_CoreAssociationScore[, 2], decreasing = TRUE), ]

save(buccal_palate_CoreAssociationScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1buccal_palate_CoreAssociationScore.RData")
write.csv(buccal_palate_CoreAssociationScore, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1buccal_palate_CoreAssociationScore.csv")
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1BuccalPalate_Core_Detection_Workspace.RData")






################## Plot the Study vs core species (1/0) as a heatmap. This will show the distribution of core species across studies.

buccal_palate_CoreKeyStoneDf3 <- buccal_palate_CoreKeyStoneDf2[,1:10]
buccal_palate_CoreKeyStoneDf3 <- buccal_palate_CoreKeyStoneDf3[,colSums(buccal_palate_CoreKeyStoneDf3)>0]
buccal_palate_CoreKeyStoneDf3 <- buccal_palate_CoreKeyStoneDf3[rowSums(buccal_palate_CoreKeyStoneDf3)>0,]


library(pheatmap)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1BuccalPalate_heatmap_StudyWise_CoreSpecies.pdf", width = 40, height = 20)
pheatmap(t(buccal_palate_CoreKeyStoneDf3),
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
  t(buccal_palate_CoreKeyStoneDf3),
  cluster_rows = T,
  cluster_cols = F,
  silent = TRUE
)

ordered_cols <- rownames(t(buccal_palate_CoreKeyStoneDf3))[ph$tree_row$order]
buccal_palate_CoreKeyStoneDf4 <- buccal_palate_CoreKeyStoneDf3[,ordered_cols ]

write.csv(buccal_palate_CoreKeyStoneDf4, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1BuccalPalate_heatmap_carpet_StudyWise_CoreSpecies.csv")

################### Get the line plot for the core species distribution across studies. This will show the number of studies in which the species is detected as core species.
collective_CoreSpecies_distribution <- data.frame(species = rownames(buccal_palate_CoreKeyStoneDf3), number_of_studies = rowSums(buccal_palate_CoreKeyStoneDf3))

all(rownames(collective_CoreSpecies_distribution) == rownames(buccal_palate_CoreKeyStoneDf3))
# its true. Now plot the line plot.

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1BuccalPalate_lineplot_CoreSpecies.pdf", width = 40, height = 10)
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


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1BuccalPalate_Core_Detection_Workspace.RData")
