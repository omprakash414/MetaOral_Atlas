


########### Load all the scores computed from three different Analysis 
## Core Association Score
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S3_CoreAssociation/S3_1saliva_CoreAssociationScore.RData")
## Health Association Score
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S4_HealthAssociation_modified/S4_1Saliva_HealthAssociationScore.RData")


## Stability Association Score
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S5_StabilityAssociation/S5_2Saliva_StabilityAssociationScore.RData")


Saliva_DiseaseAnalysis_HealthScore <- Saliva_DiseaseAnalysis_HealthScore[rownames(saliva_CoreAssociationScore), , drop = FALSE]
Saliva_StabilityScore <- Saliva_StabilityScore[rownames(saliva_CoreAssociationScore), , drop = FALSE]

saliva_CoreAssociationScore$species <- NULL

########### Combine all the three dfs
Combined_Saliva_Scores <- cbind(saliva_CoreAssociationScore, Saliva_DiseaseAnalysis_HealthScore, Saliva_StabilityScore)
Combined_Saliva_Scores$species <- NULL
colnames(Combined_Saliva_Scores) <- c("CoreAssociationScore", "HealthAssociationScore", "StabilityAssociationScore")

########### Caclulate Combined Score
library(LaplacesDemon)
library(DescTools)

Combined_Saliva_Scores$HAC_Score <- apply(Combined_Saliva_Scores[,1:3],1,mean) * (1-apply(Combined_Saliva_Scores[,1:3],1,Gini))

rank_scale=function(x){
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

Combined_Saliva_Scores$HAC_Score_RankScaled <- rank_scale(Combined_Saliva_Scores$HAC_Score)
# Order based on the value in decreasing order
Combined_Saliva_Scores <- Combined_Saliva_Scores[order(Combined_Saliva_Scores$HAC_Score_RankScaled, decreasing = TRUE),]

save(Combined_Saliva_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedScores.RData")
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedAssociation_Workspace.RData")


#### create a temp df that has Combined_Saliva_Scores with extra columns of mean score of all three scores and 1-GINI of three scores.
temp_Combined_Saliva_Scores <- Combined_Saliva_Scores
temp_Combined_Saliva_Scores$mean_scores <- rowMeans(temp_Combined_Saliva_Scores[,c("CoreAssociationScore","HealthAssociationScore","StabilityAssociationScore")])
temp_Combined_Saliva_Scores$RewardGini <- 1-apply(temp_Combined_Saliva_Scores[,1:3],1,Gini)
temp_Combined_Saliva_Scores$HAC_Score_confirm <- temp_Combined_Saliva_Scores$mean_scores * temp_Combined_Saliva_Scores$RewardGini

write.csv(temp_Combined_Saliva_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Temp_Scores_Supplementary.csv")


########### Compare the Stability Score with Core Association Score

Combined_Saliva_Scores2 <- Combined_Saliva_Scores

Combined_Saliva_Scores2$quadrant <- with(
    Combined_Saliva_Scores2,
    ifelse(StabilityAssociationScore >= 0.7 & CoreAssociationScore >= 0.7, "Q1",
         ifelse(StabilityAssociationScore < 0.7 & CoreAssociationScore >= 0.7, "Q2",
                ifelse(StabilityAssociationScore < 0.7 & CoreAssociationScore < 0.7, "Q3", "Q4")))
)




library(ggplot2)
library(ggrepel)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CoreAssociation_StabilityAssociation.pdf", 
    width = 36, height = 14)

ggplot(Combined_Saliva_Scores2, aes(x = StabilityAssociationScore, y = CoreAssociationScore)) +
  geom_point(color = "black", size = 2) +  # keep dots black
  geom_text_repel(aes(label = rownames(Combined_Saliva_Scores2), color = quadrant),
                  size = 6, max.overlaps = 60) +  # color only labels
  geom_vline(xintercept = 0.7, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "solid", color = "black") +
  theme_bw() +
  xlab("Stability Association") +
  ylab("Core Association") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  scale_color_manual(values = c("Q1" = "purple", "Q2" = "red", "Q3" = "brown", "Q4" = "blue"))

dev.off()


Combined_Saliva_Scores2$quadrant2 <- with(
    Combined_Saliva_Scores2,
    ifelse(StabilityAssociationScore >= 0.7 & HealthAssociationScore >= 0.7, "Q1",
         ifelse(StabilityAssociationScore < 0.7 & HealthAssociationScore >= 0.7, "Q2",
                ifelse(StabilityAssociationScore < 0.7 & HealthAssociationScore < 0.7, "Q3", "Q4")))
)


########### Compare the Stability Score with Health Association Score
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_HealthAssociation_StabilityAssociation.pdf", 
    width = 28, height = 14)
ggplot(Combined_Saliva_Scores2, aes(x = StabilityAssociationScore, y = HealthAssociationScore)) +
  geom_point(color = "black", size = 2) +  # keep dots black
  geom_text_repel(aes(label = rownames(Combined_Saliva_Scores2), color = quadrant2),
                  size = 6, max.overlaps = 60) +  # color only labels
  geom_vline(xintercept = 0.7, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "solid", color = "black") +
  theme_bw() +
  xlab("Stability Association") +
  ylab("Health Association") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  scale_color_manual(values = c("Q1" = "purple", "Q2" = "red", "Q3" = "brown", "Q4" = "blue")) 
dev.off()


########### Compare the Health Score with Core Association Score
Combined_Saliva_Scores2$quadrant3 <- with(
    Combined_Saliva_Scores2,
    ifelse(HealthAssociationScore >= 0.7 & CoreAssociationScore >= 0.7, "Q1",
         ifelse(HealthAssociationScore < 0.7 & CoreAssociationScore >= 0.7, "Q2",
                ifelse(HealthAssociationScore < 0.7 & CoreAssociationScore < 0.7, "Q3", "Q4")))
)   

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CoreAssociation_HealthAssociation.pdf", 
    width = 28, height = 14)
ggplot(Combined_Saliva_Scores2, aes(x = HealthAssociationScore, y = CoreAssociationScore)) +
  geom_point(color = "black", size = 2) +  # keep dots black
  geom_text_repel(aes(label = rownames(Combined_Saliva_Scores2), color = quadrant3),
                  size = 6, max.overlaps = 60) +  # color only labels
  geom_vline(xintercept = 0.7, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "solid", color = "black") +
  theme_bw() +
  xlab("Health Association") +
  ylab("Core Association") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  scale_color_manual(values = c("Q1" = "purple", "Q2" = "red", "Q3" = "brown", "Q4" = "blue"))
dev.off()

write.csv(Combined_Saliva_Scores2, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1HealthAssociation_CoreAssociation_Scores.csv")
## Further, remove species with zero core association score and plot the scatter plot again. See line number 425 onwards

########### Extract the species in different quadrants

## Stability vs Core
rownames(Combined_Saliva_Scores2[Combined_Saliva_Scores2$StabilityAssociationScore >= 0.7 & Combined_Saliva_Scores2$CoreAssociationScore >= 0.7, ])
#  [1] "Fusobacterium_periodonticum"   "Lachnoanaerobaculum_umeaense" 
#  [3] "Campylobacter_concisus"        "Porphyromonas_catoniae"       
#  [5] "Actinomyces_graevenitzii"      "Prevotella_pallens"           
#  [7] "Prevotella_melaninogenica"     "Alloprevotella_rava"          
#  [9] "Eubacterium_sulci"             "Veillonella_rogosae"          
# [11] "Prevotella_veroralis"          "Catonella_morbi"              
# [13] "Prevotella_shahii"             "Eubacterium_yurii"            
# [15] "Haemophilus_parainfluenzae"    "Cardiobacterium_hominis"      
# [17] "Neisseria_elongata"            "Capnocytophaga_gingivalis"    
# [19] "Campylobacter_showae"          "Leptotrichia_buccalis"        
# [21] "Leptotrichia_hofstadii"        "Lachnoanaerobaculum_saburreum"
# [23] "Capnocytophaga_leadbetteri"    "Cardiobacterium_valvarum"     
# [25] "Actinomyces_odontolyticus"     "Haemophilus_pittmaniae"       
# [27] "Leptotrichia_wadei"            "Oribacterium_sinus"           
# [29] "Propionibacterium_propionicum" "Leptotrichia_goodfellowii"    
# [31] "Prevotella_oulorum"            "Prevotella_salivae"           
# [33] "Prevotella_loescheii"          "Leptotrichia_hongkongensis"   
# [35] "Gemella_sanguinis"             "Porphyromonas_gingivalis"     
# [37] "Leptotrichia_trevisanii"       "Stomatobaculum_longum"        
# [39] "Atopobium_parvulum"            "Capnocytophaga_sputigena"     
# [41] "Prevotella_maculosa"           "Solobacterium_moorei"         
# [43] "Veillonella_parvula"           "Tannerella_forsythia"         
# [45] "Alloprevotella_tannerae"       "Prevotella_denticola"         
# [47] "Prevotella_nanceiensis"        "Megasphaera_micronuciformis"  
# [49] "Capnocytophaga_granulosa"      "Prevotella_histicola"         
# [51] "Streptococcus_infantis"        "Peptostreptococcus_stomatis"  
# [53] "Selenomonas_infelix"           "Abiotrophia_defectiva"        
# [55] "Centipeda_periodontii"         "Lachnoanaerobaculum_orale"    
# [57] "Selenomonas_sputigena"         "Veillonella_atypica"          
# [59] "Prevotella_intermedia"

## Stability vs Health
rownames(Combined_Saliva_Scores2[Combined_Saliva_Scores2$StabilityAssociationScore >= 0.7 & Combined_Saliva_Scores2$HealthAssociationScore >= 0.7, ])
#  [1] "Fusobacterium_periodonticum"             
#  [2] "Lachnoanaerobaculum_umeaense"            
#  [3] "Campylobacter_concisus"                  
#  [4] "Porphyromonas_catoniae"                  
#  [5] "Actinomyces_graevenitzii"                
#  [6] "Prevotella_pallens"                      
#  [7] "Prevotella_melaninogenica"               
#  [8] "Alloprevotella_rava"                     
#  [9] "Eubacterium_sulci"                       
# [10] "Veillonella_rogosae"                     
# [11] "Prevotella_veroralis"                    
# [12] "Catonella_morbi"                         
# [13] "Prevotella_shahii"                       
# [14] "Eubacterium_yurii"                       
# [15] "Haemophilus_parainfluenzae"              
# [16] "Cardiobacterium_hominis"                 
# [17] "Neisseria_elongata"                      
# [18] "Capnocytophaga_gingivalis"               
# [19] "Campylobacter_showae"                    
# [20] "Leptotrichia_buccalis"                   
# [21] "Leptotrichia_hofstadii"                  
# [22] "Lachnoanaerobaculum_saburreum"           
# [23] "Streptococcus_lactarius"                 
# [24] "Prevotella_aurantiaca"                   
# [25] "Capnocytophaga_leadbetteri"              
# [26] "Cardiobacterium_valvarum"                
# [27] "Actinomyces_odontolyticus"               
# [28] "Haemophilus_pittmaniae"                  
# [29] "Propionibacterium_propionicum"           
# [30] "Neisseria_meningitidis"                  
# [31] "Streptococcus_mutans"                    
# [32] "Orientia_tsutsugamushi"                  
# [33] "Neisseria_oralis"                        
# [34] "Streptobacillus_moniliformis"            
# [35] "Veillonella_sp_T11011_6"                 
# [36] "Capnocytophaga_haemolytica"              
# [37] "Lachnospiraceae_bacterium_oral_taxon_096"
# [38] "Klebsiella_pneumoniae"                   
# [39] "Prevotella_scopos"                       
# [40] "Moryella_indoligenes"                    
# [41] "Simonsiella_muelleri"                    
# [42] "Treponema_refringens"                    
# [43] "Streptococcus_sp_A12"                    
# [44] "Mobiluncus_mulieris" 

## Stability vs Health vs Core
rownames(Combined_Saliva_Scores2[Combined_Saliva_Scores2$StabilityAssociationScore >= 0.7 & Combined_Saliva_Scores2$HealthAssociationScore >= 0.7 & Combined_Saliva_Scores2$CoreAssociationScore >= 0.7, ])
#  [1] "Fusobacterium_periodonticum"   "Lachnoanaerobaculum_umeaense" 
#  [3] "Campylobacter_concisus"        "Porphyromonas_catoniae"       
#  [5] "Actinomyces_graevenitzii"      "Prevotella_pallens"           
#  [7] "Prevotella_melaninogenica"     "Alloprevotella_rava"          
#  [9] "Eubacterium_sulci"             "Veillonella_rogosae"          
# [11] "Prevotella_veroralis"          "Catonella_morbi"              
# [13] "Prevotella_shahii"             "Eubacterium_yurii"            
# [15] "Haemophilus_parainfluenzae"    "Cardiobacterium_hominis"      
# [17] "Neisseria_elongata"            "Capnocytophaga_gingivalis"    
# [19] "Campylobacter_showae"          "Leptotrichia_buccalis"        
# [21] "Leptotrichia_hofstadii"        "Lachnoanaerobaculum_saburreum"
# [23] "Capnocytophaga_leadbetteri"    "Cardiobacterium_valvarum"     
# [25] "Actinomyces_odontolyticus"     "Haemophilus_pittmaniae"       
# [27] "Propionibacterium_propionicum"

## Health vs Core
rownames(Combined_Saliva_Scores2[Combined_Saliva_Scores2$HealthAssociationScore >= 0.7 & Combined_Saliva_Scores2$CoreAssociationScore >= 0.7, ])
#  [1] "Fusobacterium_periodonticum"   "Lachnoanaerobaculum_umeaense" 
#  [3] "Campylobacter_concisus"        "Porphyromonas_catoniae"       
#  [5] "Actinomyces_graevenitzii"      "Prevotella_pallens"           
#  [7] "Prevotella_melaninogenica"     "Alloprevotella_rava"          
#  [9] "Eubacterium_sulci"             "Veillonella_rogosae"          
# [11] "Prevotella_veroralis"          "Catonella_morbi"              
# [13] "Prevotella_shahii"             "Eubacterium_yurii"            
# [15] "Haemophilus_parainfluenzae"    "Cardiobacterium_hominis"      
# [17] "Neisseria_elongata"            "Capnocytophaga_gingivalis"    
# [19] "Campylobacter_showae"          "Leptotrichia_buccalis"        
# [21] "Leptotrichia_hofstadii"        "Lachnoanaerobaculum_saburreum"
# [23] "Capnocytophaga_leadbetteri"    "Cardiobacterium_valvarum"     
# [25] "Actinomyces_odontolyticus"     "Haemophilus_pittmaniae"       
# [27] "Propionibacterium_propionicum" "Neisseria_flavescens"         
# [29] "Corynebacterium_durum"         "Haemophilus_sputorum"         
# [31] "Treponema_vincentii"           "Aggregatibacter_segnis"       
# [33] "Kingella_oralis"               "Rothia_aeria"                 
# [35] "Neisseria_mucosa"              "Prevotella_nigrescens"        
# [37] "Streptococcus_sanguinis"       "Campylobacter_gracilis"       
# [39] "Treponema_socranskii"          "Streptococcus_thermophilus"

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedAssociation_Workspace.RData")



####################
####################
##### heatmap that have all the species with HAC_Score_RankScaled > 0.80 score
filtered_df <- Combined_Saliva_Scores2[
  Combined_Saliva_Scores2$HAC_Score_RankScaled >= 0.80,
]

library(ggplot2)
library(reshape2)

#---------------------------
# 1. Select required columns
#---------------------------
heat_df <- filtered_df[, c(
  "CoreAssociationScore",
  "HealthAssociationScore",
  "StabilityAssociationScore"
)]
heat_df$Taxa <- rownames(heat_df)

#-----------------------------------
# 2. Identify taxa with all >= 0.70
#-----------------------------------
taxa_flag <- apply(
  heat_df[, c(
    "CoreAssociationScore",
    "HealthAssociationScore",
    "StabilityAssociationScore"
  )],
  1,
  function(x) all(x >= 0.70)
)

taxa_color_df <- data.frame(
  Taxa = heat_df$Taxa,
  label_color = ifelse(taxa_flag, "red", "black")
)

write.csv(heat_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_Carpet_Heatmap_CombinedScores_HAC70_AllScores70.csv")
#---------------------------
# 3. Convert to long format
#---------------------------
heat_long <- melt(
  heat_df,
  id.vars = "Taxa",
  variable.name = "ScoreType",
  value.name = "Score"
)

heat_long$star <- ifelse(heat_long$Score >= 0.70, "*", "")

heat_long$Taxa <- factor(heat_long$Taxa, levels = heat_df$Taxa)
heat_long$ScoreType <- factor(
  heat_long$ScoreType,
  levels = c(
    "CoreAssociationScore",
    "HealthAssociationScore",
    "StabilityAssociationScore"
  )
)

#---------------------------
# 4. Plot heatmap
#---------------------------

pdf(
  "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_Heatmap_CombinedScores_HAC70_AllScores70.pdf",
  width = 14,
  height = 5  # slightly shorter height
)

ggplot(heat_long, aes(x = Taxa, y = ScoreType, fill = Score)) +
  geom_tile(color = "black", linewidth = 0.4) +
  geom_text(aes(label = star), size = 5) +

  geom_text(
    data = taxa_color_df,
    aes(
      x = Taxa,
      y = 0.5,
      label = Taxa,
      color = label_color
    ),
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_color_identity() +

  scale_fill_gradient(
    low  = "#e6eff7",
    high = "#3992d1",
    name = "Score"
  ) +

  scale_y_discrete(expand = expansion(add = c(1.8, 0))) +

  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 30, 10)
  ) +
  coord_cartesian(clip = "off")

dev.off()



save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedAssociation_Workspace.RData")





##########################
########################## Get the species that are present in the top 30 percentile for each of the scores.
### get 70th percentile cutoffs
cs_cutoff <- quantile(Combined_Saliva_Scores2$CoreAssociationScore, 0.65, na.rm = TRUE)
hs_cutoff <- quantile(Combined_Saliva_Scores2$HealthAssociationScore, 0.65, na.rm = TRUE)
ss_cutoff <- quantile(Combined_Saliva_Scores2$StabilityAssociationScore, 0.65, na.rm = TRUE)


### Filter species
selected_species <- rownames(
  Combined_Saliva_Scores2[
    Combined_Saliva_Scores2$CoreAssociationScore >= cs_cutoff &
    Combined_Saliva_Scores2$HealthAssociationScore >= hs_cutoff &
    Combined_Saliva_Scores2$StabilityAssociationScore >= ss_cutoff,
  ]
)

selected_species
length(selected_species)

### Plot it as well as save this in a csv file
top_35_percentile_df <- Combined_Saliva_Scores2[selected_species, ]

## get the heatmap for this top 35 percentile species for each of the scores. Values should be printed in the cells. While color should be based on the values it have i.e give shade of one color i.e from white to color1
library(pheatmap)
breaks = seq(0.40, 1.00, length.out = 101)
mat <- as.matrix(top_35_percentile_df[, 1:3])
breaks_percentile <- quantile(
  mat,
  probs = seq(0, 1, length.out = 101),
  na.rm = TRUE
)
breaks_percentile <- unique(breaks_percentile)
my_colors <- colorRampPalette(c("white","orange2"))(length(breaks_percentile) - 1)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Top35Percentile_Species_for_eachScore_Heatmap.pdf", width = 12, height = 30)
pheatmap(
  top_35_percentile_df[,1:3],
  color = my_colors,
  breaks = breaks_percentile,

  display_numbers = round(top_35_percentile_df[,1:3], 2),
  number_color = "black",
  fontsize_number = 18,

  fontsize_row = 18,
  fontsize_col = 18,
  cellheight = 30,
  cellwidth = 50,
  cluster_rows = T,
  cluster_cols = T,
  border_color = "black",
  treeheight_row = 0,
  treeheight_col = 0
)
dev.off()

########## Get the carpet
ph <- pheatmap(
  top_35_percentile_df[,1:3],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  silent = TRUE
)

ordered_cols <- colnames(top_35_percentile_df[,1:3])[ph$tree_col$order]
ordered_rows <- rownames(top_35_percentile_df)[ph$tree_row$order]
top_35_percentile_df <- top_35_percentile_df[ordered_rows, ordered_cols]

write.csv(top_35_percentile_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Top35Percentile_Species_for_eachScore.csv")


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedAssociation_Workspace.RData")






######### Plot the reduced species with zero core association score removed and plot the scatter plot again.
## Core vs health
Combined_Saliva_Scores2_filtered <- Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore > 0, ]
library(ggplot2)
library(ggrepel)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CoreAssociation_HealthAssociation_without_ZeroCS.pdf", 
    width = 28, height = 14)
ggplot(Combined_Saliva_Scores2_filtered, aes(x = HealthAssociationScore, y = CoreAssociationScore)) +
  geom_point(color = "black", size = 3) +  # keep dots black
  geom_text_repel(aes(label = rownames(Combined_Saliva_Scores2_filtered), color = quadrant3),
                  size = 7, max.overlaps = 60) +  # color only labels
  geom_vline(xintercept = 0.7, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "solid", color = "black") +
  theme_bw() +
  xlab("Health Association") +
  ylab("Core Association") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  scale_color_manual(values = c("Q1" = "purple", "Q2" = "red", "Q3" = "brown", "Q4" = "blue"))
dev.off()

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedAssociation_Workspace.RData")

## species with CS and HS >=0.70
CS_HS <- rownames(Combined_Saliva_Scores2_filtered[Combined_Saliva_Scores2_filtered$CoreAssociationScore >= 0.70 & Combined_Saliva_Scores2_filtered$HealthAssociationScore >= 0.70,])

########### Plot the scatter plot taking only species for Health vs stability association score. and then core vs stability score which qualify score > 0.50
## for Health vs stability
Combined_Saliva_Scores2_SSHS_filtered_50 <- Combined_Saliva_Scores2[Combined_Saliva_Scores2$HealthAssociationScore >= 0.60 & Combined_Saliva_Scores2$StabilityAssociationScore >= 0.60, ]

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_HealthAssociation_StabilityAssociation_Atleast60_score.pdf", 
    width = 20, height = 20)
ggplot(Combined_Saliva_Scores2_SSHS_filtered_50, aes(x = StabilityAssociationScore, y = HealthAssociationScore)) +
  geom_point(color = "black", size = 5) +  # keep dots black
  geom_text_repel(aes(label = rownames(Combined_Saliva_Scores2_SSHS_filtered_50), color = quadrant2),
                  size = 9, max.overlaps = 60) +  # color only labels
  geom_vline(xintercept = 0.7, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "solid", color = "black") +
  theme_bw() +
  xlab("Stability Association") +
  ylab("Health Association") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  scale_color_manual(values = c("Q1" = "purple", "Q2" = "red", "Q3" = "brown", "Q4" = "blue")) 
dev.off()

## species with HS and SS >=0.70
HS_SS <- rownames(Combined_Saliva_Scores2_SSHS_filtered_50[Combined_Saliva_Scores2_SSHS_filtered_50$HealthAssociationScore >= 0.70 & Combined_Saliva_Scores2_SSHS_filtered_50$StabilityAssociationScore >= 0.70,])


## for stability vs core association scores
Combined_Saliva_Scores2_SSCS_filtered_50 <- Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore >= 0.60 & Combined_Saliva_Scores2$StabilityAssociationScore >= 0.60, ]

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CoreAssociation_StabilityAssociation_Atleast_60Score.pdf", 
    width = 20, height = 20)

ggplot(Combined_Saliva_Scores2_SSCS_filtered_50, aes(x = StabilityAssociationScore, y = CoreAssociationScore)) +
  geom_point(color = "black", size = 5) +  # keep dots black
  geom_text_repel(aes(label = rownames(Combined_Saliva_Scores2_SSCS_filtered_50), color = quadrant),
                  size = 9, max.overlaps = 60) +  # color only labels
  geom_vline(xintercept = 0.7, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "solid", color = "black") +
  theme_bw() +
  xlab("Stability Association") +
  ylab("Core Association") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  scale_color_manual(values = c("Q1" = "purple", "Q2" = "red", "Q3" = "brown", "Q4" = "blue"))

dev.off()


# species with CS and SS >=0.70
CS_SS <- rownames(Combined_Saliva_Scores2_SSCS_filtered_50[Combined_Saliva_Scores2_SSCS_filtered_50$CoreAssociationScore >= 0.70 & Combined_Saliva_Scores2_SSCS_filtered_50$StabilityAssociationScore >=0.70,])


## common species in all the three combinations.
intersect(intersect(CS_HS, HS_SS), CS_SS)
#  [1] "Fusobacterium_periodonticum"   "Lachnoanaerobaculum_umeaense" 
#  [3] "Campylobacter_concisus"        "Porphyromonas_catoniae"       
#  [5] "Actinomyces_graevenitzii"      "Prevotella_pallens"           
#  [7] "Prevotella_melaninogenica"     "Alloprevotella_rava"          
#  [9] "Eubacterium_sulci"             "Veillonella_rogosae"          
# [11] "Prevotella_veroralis"          "Catonella_morbi"              
# [13] "Prevotella_shahii"             "Eubacterium_yurii"            
# [15] "Haemophilus_parainfluenzae"    "Cardiobacterium_hominis"      
# [17] "Neisseria_elongata"            "Capnocytophaga_gingivalis"    
# [19] "Campylobacter_showae"          "Leptotrichia_buccalis"        
# [21] "Leptotrichia_hofstadii"        "Lachnoanaerobaculum_saburreum"
# [23] "Capnocytophaga_leadbetteri"    "Cardiobacterium_valvarum"     
# [25] "Actinomyces_odontolyticus"     "Haemophilus_pittmaniae"       
# [27] "Propionibacterium_propionicum"

# confirm if these 27 are same as species with all three scores >=0.70
nrow(Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore >= 0.70 & Combined_Saliva_Scores2$HealthAssociationScore >= 0.70 & Combined_Saliva_Scores2$StabilityAssociationScore >= 0.70, ])
#confirmed, 27 species are same as the above list.



########## To get the phylogeny of species with sHACK score >=0.90, first extract the species
species_sHACK_90 <- rownames(Combined_Saliva_Scores[Combined_Saliva_Scores$HAC_Score_RankScaled >= 0.90, ])

species_clean <- gsub("_", " ", species_sHACK_90)

writeLines(species_clean,"/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_sHACK_species90.txt")

# for this list of speceis I have run the script to get the phylogeny of these species.
# Import that in other script and do the needful.

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedAssociation_Workspace.RData")



############# Calculate the HAC score for saliva. (HAC includes only Health scores and Core association scores. stability score is not included in the HAC score)
Saliva_DiseaseAnalysis_HealthScore2 <- Saliva_DiseaseAnalysis_HealthScore[rownames(saliva_CoreAssociationScore), , drop = FALSE]

all(rownames(Saliva_DiseaseAnalysis_HealthScore2) == rownames(saliva_CoreAssociationScore)) #TRUE

########### Combine all the three dfs
saliva_HAC_Score <- cbind(saliva_CoreAssociationScore, Saliva_DiseaseAnalysis_HealthScore2)
colnames(saliva_HAC_Score) <- c("CoreAssociationScore", "HealthAssociationScore")

########### Caclulate Combined Score
library(LaplacesDemon)
library(DescTools)

saliva_HAC_Score$HAC_Score <- apply(saliva_HAC_Score[,1:2],1,mean) * (1-apply(saliva_HAC_Score[,1:2],1,Gini))

saliva_HAC_Score$HAC_Score_RankScaled <- rank_scale(saliva_HAC_Score$HAC_Score)
# Order based on the value in decreasing order
saliva_HAC_Score <- saliva_HAC_Score[order(saliva_HAC_Score$HAC_Score_RankScaled, decreasing = TRUE),]

save(saliva_HAC_Score, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_HAC_Score_noHACK.RData")

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_1Saliva_CombinedAssociation_Workspace.RData")











######## Correlation of each of the quadrant in health and core association plot.
##### Quadrant 3
cor.test(Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore <0.7 & Combined_Saliva_Scores2$HealthAssociationScore<0.7,]$CoreAssociationScore,
          Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore <0.7 & Combined_Saliva_Scores2$HealthAssociationScore<0.7,]$HealthAssociationScore)

#         Pearson's product-moment correlation

# data:  Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore < 0.7 & Combined_Saliva_Scores2$HealthAssociationScore < 0.7, ]$CoreAssociationScore and Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore < 0.7 & Combined_Saliva_Scores2$HealthAssociationScore < 0.7, ]$HealthAssociationScore
# t = 0.76878, df = 314, p-value = 0.4426
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.06731039  0.15294551
# sample estimates:
#        cor 
# 0.04334422


##### Quadrant 2
cor.test(Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore >=0.7 & Combined_Saliva_Scores2$HealthAssociationScore<0.7,]$CoreAssociationScore,
          Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore >=0.7 & Combined_Saliva_Scores2$HealthAssociationScore<0.7,]$HealthAssociationScore)

#         Pearson's product-moment correlation

# data:  Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore >= 0.7 & Combined_Saliva_Scores2$HealthAssociationScore < 0.7, ]$CoreAssociationScore and Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore >= 0.7 & Combined_Saliva_Scores2$HealthAssociationScore < 0.7, ]$HealthAssociationScore
# t = 0.20075, df = 68, p-value = 0.8415
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2118482  0.2578368
# sample estimates:
#        cor 
# 0.02433726 


##### Quadrant 4
cor.test(Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore <0.7 & Combined_Saliva_Scores2$HealthAssociationScore>=0.7,]$CoreAssociationScore,
          Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore <0.7 & Combined_Saliva_Scores2$HealthAssociationScore>=0.7,]$HealthAssociationScore)

#         Pearson's product-moment correlation

# data:  Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore < 0.7 & Combined_Saliva_Scores2$HealthAssociationScore >= 0.7, ]$CoreAssociationScore and Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore < 0.7 & Combined_Saliva_Scores2$HealthAssociationScore >= 0.7, ]$HealthAssociationScore
# t = 2.2415, df = 71, p-value = 0.02812
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.02870773 0.45994113
# sample estimates:
#       cor 
# 0.2570771


##### Quadrant 1
cor.test(Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore >=0.7 & Combined_Saliva_Scores2$HealthAssociationScore>=0.7,]$CoreAssociationScore,
          Combined_Saliva_Scores2[Combined_Saliva_Scores2$CoreAssociationScore >=0.7 & Combined_Saliva_Scores2$HealthAssociationScore>=0.7,]$HealthAssociationScore)


Correlations of HS and CS in each quadrant are as follows:
             Cor        Pvalue
Quadrant 1   0.04        0.81
Quadrant 2   0.02        0.84
Quadrant 3   0.04        0.44
Quadrant 4   0.26        0.03