#### This script is to determine thresholds for selecting oral species associated with the different body_site_category
#### It computes detection patterns and mean abundances of all species in each study separately
#### Then it evaluates a grid of thresholds to see how many species would be selected and how well they represent the oral microbiome
#### Finally, it selects a threshold combination based on the results

library(dplyr)
library(pheatmap)
library(ggplot2)

############ Load the code library functions
source("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/code_library_Metaoral.R")

############ Load the data
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_BuccalPalate_ControlCohort.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Saliva_ControlCohort.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Supragingival_ControlCohort.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_TongueTonsil_ControlCohort.RData")
load("/storage/omprakash/MetaOral_Analysis/S0_samples_distribution/S0_Subgingival_ControlCohort.RData")


# Remove the object which is not needed here.
rm(list = ls(pattern = "^SpDf_nonorm"))
gc()


############ first see if the rownames of metdata and species df are in same order for all subsite's dataframes and then do add the study_name column to species df
all(rownames(SpDf_saliva_Control) == rownames(MetadataDf_saliva_control)) # TRUE 
SpDf_saliva_Control$study_name <- MetadataDf_saliva_control$study_name

all(rownames(SpDf_supragingival_Control) == rownames(MetadataDf_supragingival_control)) # TRUE
SpDf_supragingival_Control$study_name <- MetadataDf_supragingival_control$study_name

all(rownames(SpDf_tongue_tonsil_Control) == rownames(MetadataDf_tongue_tonsil_control)) # TRUE
SpDf_tongue_tonsil_Control$study_name <- MetadataDf_tongue_tonsil_control$study_name

all(rownames(SpDf_subgingival_Control) == rownames(MetadataDf_subgingival_control)) # TRUE
SpDf_subgingival_Control$study_name <- MetadataDf_subgingival_control$study_name

all(rownames(SpDf_buccal_palate_Control) == rownames(MetadataDf_buccal_palate_control)) # TRUE
SpDf_buccal_palate_Control$study_name <- MetadataDf_buccal_palate_control$study_name


###################################################################################


############ compute the prevalence / detection pattern for each study
## First extract the species in each of the subsite
species_cols_saliva <- setdiff(colnames(SpDf_saliva_Control), "study_name")
list_sample_groups_saliva <- unique(SpDf_saliva_Control$study_name)

species_cols_supragingival <- setdiff(colnames(SpDf_supragingival_Control), "study_name")
list_sample_groups_supragingival <- unique(SpDf_supragingival_Control$study_name)

species_cols_tonguetonsil <- setdiff(colnames(SpDf_tongue_tonsil_Control), "study_name")
list_sample_groups_tonguetonsil <- unique(SpDf_tongue_tonsil_Control$study_name)

species_cols_subgingival <- setdiff(colnames(SpDf_subgingival_Control), "study_name")
list_sample_groups_subgingival <- unique(SpDf_subgingival_Control$study_name)

species_cols_buccal <- setdiff(colnames(SpDf_buccal_palate_Control), "study_name")
list_sample_groups_buccal <- unique(SpDf_buccal_palate_Control$study_name)

## Run the compute detection function for each subsite
AllSpeciesDetectionPattern_saliva <- compute_detection(SpDf_saliva_Control,species_cols_saliva,"study_name",list_sample_groups_saliva)
gc()
AllSpeciesDetectionPattern_supragingival <- compute_detection(SpDf_supragingival_Control,species_cols_supragingival,"study_name",list_sample_groups_supragingival)
AllSpeciesDetectionPattern_tonguetonsil <- compute_detection(SpDf_tongue_tonsil_Control,species_cols_tonguetonsil,"study_name",list_sample_groups_tonguetonsil)
AllSpeciesDetectionPattern_subgingival <- compute_detection(SpDf_subgingival_Control,species_cols_subgingival,"study_name",list_sample_groups_subgingival)
AllSpeciesDetectionPattern_buccal <- compute_detection(SpDf_buccal_palate_Control,species_cols_buccal,"study_name",list_sample_groups_buccal)
gc()



###################################################################################


############ compute Mean Abundance of each species in each study separately
AllSpeciesMeanAbundance_saliva <- compute_mean_abundance(SpDf_saliva_Control,species_cols_saliva,"study_name")
gc()
AllSpeciesMeanAbundance_supragingival <- compute_mean_abundance(SpDf_supragingival_Control,species_cols_supragingival,"study_name")
AllSpeciesMeanAbundance_tonguetonsil <- compute_mean_abundance(SpDf_tongue_tonsil_Control,species_cols_tonguetonsil,"study_name")
AllSpeciesMeanAbundance_subgingival <- compute_mean_abundance(SpDf_subgingival_Control,species_cols_subgingival,"study_name")
AllSpeciesMeanAbundance_buccal <- compute_mean_abundance(SpDf_buccal_palate_Control,species_cols_buccal,"study_name")
gc()



###################################################################################


############ Evaluate thresholds for each subsite
threshold_results_saliva <- evaluate_threshold_grid(AllSpeciesDetectionPattern_saliva ,AllSpeciesMeanAbundance_saliva,list_sample_groups_saliva)
df_numb_species_saliva <- threshold_results_saliva$df_numb_species
df_representation_90_plus_saliva <- threshold_results_saliva$df_representation_90_plus
df_representation_70_minus_saliva <- threshold_results_saliva$df_representation_70_minus
df_associated_identification_saliva <- threshold_results_saliva$df_associated_identification

threshold_results_supragingival <- evaluate_threshold_grid(AllSpeciesDetectionPattern_supragingival ,AllSpeciesMeanAbundance_supragingival,list_sample_groups_supragingival)
df_numb_species_supragingival <- threshold_results_supragingival$df_numb_species
df_representation_90_plus_supragingival <- threshold_results_supragingival$df_representation_90_plus
df_representation_70_minus_supragingival <- threshold_results_supragingival$df_representation_70_minus
df_associated_identification_supragingival <- threshold_results_supragingival$df_associated_identification

threshold_results_tonguetonsil <- evaluate_threshold_grid(AllSpeciesDetectionPattern_tonguetonsil ,AllSpeciesMeanAbundance_tonguetonsil,list_sample_groups_tonguetonsil)
df_numb_species_tonguetonsil <- threshold_results_tonguetonsil$df_numb_species
df_representation_90_plus_tonguetonsil <- threshold_results_tonguetonsil$df_representation_90_plus
df_representation_70_minus_tonguetonsil <- threshold_results_tonguetonsil$df_representation_70_minus
df_associated_identification_tonguetonsil <- threshold_results_tonguetonsil$df_associated_identification

threshold_results_subgingival <- evaluate_threshold_grid(AllSpeciesDetectionPattern_subgingival ,AllSpeciesMeanAbundance_subgingival,list_sample_groups_subgingival)
df_numb_species_subgingival <- threshold_results_subgingival$df_numb_species
df_representation_90_plus_subgingival <- threshold_results_subgingival$df_representation_90_plus
df_representation_70_minus_subgingival <- threshold_results_subgingival$df_representation_70_minus
df_associated_identification_subgingival <- threshold_results_subgingival$df_associated_identification

threshold_results_buccal <- evaluate_threshold_grid(AllSpeciesDetectionPattern_buccal ,AllSpeciesMeanAbundance_buccal,list_sample_groups_buccal)
df_numb_species_buccal <- threshold_results_buccal$df_numb_species
df_representation_90_plus_buccal <- threshold_results_buccal$df_representation_90_plus
df_representation_70_minus_buccal <- threshold_results_buccal$df_representation_70_minus
df_associated_identification_buccal <- threshold_results_buccal$df_associated_identification



###################################################################################



############ Plotting results saliva

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_saliva_df_associated_identification.pdf",width = 8, height = 7)
ggplot(df_associated_identification_saliva) +
  geom_point(
    aes(x = number_of_species, y = representation_90_plus), color = "#1B9E77", size = 2.8, alpha = 0.9) +
  geom_point(
    aes(x = number_of_species, y = representation_70_minus), color = "#D95F02", size = 2.8, alpha = 0.9) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Number of taxa (N)", y = "Proportion of microbiomes",
    title = "Cumulative representation of selected taxa across threshold combinations") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank())
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_saliva_df_representation_90_plus.pdf", width = 8, height = 6)
pheatmap(df_representation_90_plus_saliva,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5"))(100),
  display_numbers = df_numb_species_saliva,
  number_color = "black",
  fontsize_number = 8,
  fontsize_row = 10,
  fontsize_col = 10,
  legend = TRUE,
  main = "Proportion of microbiomes with ≥90% cumulative abundance")
dev.off()


############ Plotting results supragingival
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_supragingival_df_associated_identification.pdf",width = 8, height = 7)
ggplot(df_associated_identification_supragingival) +
  geom_point(
    aes(x = number_of_species, y = representation_90_plus), color = "#1B9E77", size = 2.8, alpha = 0.9) +
  geom_point(
    aes(x = number_of_species, y = representation_70_minus), color = "#D95F02", size = 2.8, alpha = 0.9) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Number of taxa (N)", y = "Proportion of microbiomes",
    title = "Cumulative representation of selected taxa across threshold combinations") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank())
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_supragingival_df_representation_90_plus.pdf", width = 8, height = 7)
pheatmap(df_representation_90_plus_supragingival,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5"))(100),
  display_numbers = df_numb_species_supragingival,
  number_color = "black",
  fontsize_number = 8,
  fontsize_row = 10,
  fontsize_col = 10,
  legend = TRUE,
  main = "Proportion of microbiomes with ≥90% cumulative abundance")
dev.off()


################ Plotting results tonguetonsil
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_tonguetonsil_df_associated_identification.pdf",width = 8, height = 7)
ggplot(df_associated_identification_tonguetonsil) +
  geom_point(
    aes(x = number_of_species, y = representation_90_plus), color = "#1B9E77", size = 2.8, alpha = 0.9) +
  geom_point(
    aes(x = number_of_species, y = representation_70_minus), color = "#D95F02", size = 2.8, alpha = 0.9) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Number of taxa (N)", y = "Proportion of microbiomes",
    title = "Cumulative representation of selected taxa across threshold combinations") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank())
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_tonguetonsil_df_representation_90_plus.pdf", width = 8, height = 7)
pheatmap(df_representation_90_plus_tonguetonsil,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5"))(100),
  display_numbers = df_numb_species_tonguetonsil,
  number_color = "black",
  fontsize_number = 8,
  fontsize_row = 10,
  fontsize_col = 10,
  legend = TRUE,
  main = "Proportion of microbiomes with ≥90% cumulative abundance")
dev.off()


############ Plotting results subgingival
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_subgingival_df_associated_identification.pdf",width = 8, height = 7)
ggplot(df_associated_identification_subgingival) +
  geom_point(
    aes(x = number_of_species, y = representation_90_plus), color = "#1B9E77", size = 2.8, alpha = 0.9) +
  geom_point(
    aes(x = number_of_species, y = representation_70_minus), color = "#D95F02", size = 2.8, alpha = 0.9) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Number of taxa (N)", y = "Proportion of microbiomes",
    title = "Cumulative representation of selected taxa across threshold combinations") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank())
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_subgingival_df_representation_90_plus.pdf", width = 8, height = 7)
pheatmap(df_representation_90_plus_subgingival,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5"))(100),
  display_numbers = df_numb_species_subgingival,
  number_color = "black",
  fontsize_number = 8,
  fontsize_row = 10,
  fontsize_col = 10,
  legend = TRUE,
  main = "Proportion of microbiomes with ≥90% cumulative abundance")
dev.off()


############ Plotting results buccal
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_buccal_df_associated_identification.pdf",width = 8, height = 7)
ggplot(df_associated_identification_buccal) +
  geom_point(
    aes(x = number_of_species, y = representation_90_plus), color = "#1B9E77", size = 2.8, alpha = 0.9) +
  geom_point(
    aes(x = number_of_species, y = representation_70_minus), color = "#D95F02", size = 2.8, alpha = 0.9) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Number of taxa (N)", y = "Proportion of microbiomes",
    title = "Cumulative representation of selected taxa across threshold combinations") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank())
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_buccal_df_representation_90_plus.pdf", width = 8, height = 7)
pheatmap(df_representation_90_plus_buccal,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#F7FBFF", "#a5cee7", "#67b0da", "#2378c2"))(50),
  display_numbers = df_numb_species_buccal,
  number_color = "black",
  fontsize_number = 8,
  fontsize_row = 10,
  fontsize_col = 10,
  legend = TRUE,
  main = "Proportion of microbiomes with ≥90% cumulative abundance")
dev.off()


###################################################################################

##### Based on the 70%minus and 90%plus select the threshold for species in every subsite (manually select the species)
saliva_AssociatedSpecies <- names(which(apply(AllSpeciesDetectionPattern_saliva,1,function(x)(length(x[x>=0.05])))/ncol(AllSpeciesDetectionPattern_saliva)>=0.15))[1:499]
supragingival_AssociatedSpecies <- names(which(apply(AllSpeciesDetectionPattern_supragingival,1,function(x)(length(x[x>=0.05])))/ncol(AllSpeciesDetectionPattern_supragingival)>=0.25))[1:301]
tonguetonsil_AssociatedSpecies <- names(which(apply(AllSpeciesDetectionPattern_tonguetonsil,1,function(x)(length(x[x>=0.05])))/ncol(AllSpeciesDetectionPattern_tonguetonsil)>=0.30))[1:266]
subgingival_AssociatedSpecies <- names(which(apply(AllSpeciesDetectionPattern_subgingival,1,function(x)(length(x[x>=0.05])))/ncol(AllSpeciesDetectionPattern_subgingival)>=0.45))[1:196] # 0.81 is the 90_plus and 0.09 is 70_minus
buccal_AssociatedSpecies <- names(which(apply(AllSpeciesDetectionPattern_buccal,1,function(x)(length(x[x>=0.05])))/ncol(AllSpeciesDetectionPattern_buccal)>=0.55))[1:166]

###################################################################################
## Now filter species profiles to only associated species for each body_site_category
SpDf_saliva_Control_associated <- SpDf_saliva_Control[,saliva_AssociatedSpecies]
SpDf_supragingival_Control_associated <- SpDf_supragingival_Control[,supragingival_AssociatedSpecies]
SpDf_tongue_tonsil_Control_associated <- SpDf_tongue_tonsil_Control[,tonguetonsil_AssociatedSpecies]
SpDf_subgingival_Control_associated <- SpDf_subgingival_Control[,subgingival_AssociatedSpecies]
SpDf_buccal_palate_Control_associated <- SpDf_buccal_palate_Control[,buccal_AssociatedSpecies]  

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Threshold_Species_Detection_Workspace.RData")
# ############ Save the associated species for each body_site_category
save(saliva_AssociatedSpecies,
     supragingival_AssociatedSpecies,
     tonguetonsil_AssociatedSpecies,
     subgingival_AssociatedSpecies,
     buccal_AssociatedSpecies,
     file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_SubsiteWise.RData")

write.csv(saliva_AssociatedSpecies, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_saliva.csv", row.names = FALSE)
write.csv(supragingival_AssociatedSpecies, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_supragingival.csv", row.names = FALSE)
write.csv(tonguetonsil_AssociatedSpecies, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_tonguetonsil.csv", row.names = FALSE)
write.csv(subgingival_AssociatedSpecies, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_subgingival.csv", row.names = FALSE)
write.csv(buccal_AssociatedSpecies, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Selected_Species_buccal.csv", row.names = FALSE)

save(SpDf_saliva_Control_associated,
     SpDf_supragingival_Control_associated,
     SpDf_tongue_tonsil_Control_associated,
     SpDf_subgingival_Control_associated,
     SpDf_buccal_palate_Control_associated,
     file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S1_Threshold_Species_Detection/S1_Species_Profiles_Associated_Species_SubsiteWise.RData")








