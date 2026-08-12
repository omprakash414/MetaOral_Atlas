
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4_AllScores_AllSubsites.RData")
load( "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Supragingival_data_for_prediction.RData")


######## Now start prediction using RF
library(caret)
library(randomForest)
library(dplyr)

run_rf_feature_selection <- function(
    species_df, metadata_df, outcome_col, 
    do_feature_selection = FALSE, rfe_sizes = c(5, 10, 20, 30)) {
  
  # -----------------------------
  # 1. Align metadata
  metadata_df <- metadata_df[match(rownames(species_df), rownames(metadata_df)), ]
  outcome <- metadata_df[[outcome_col]]
  full_df <- data.frame(Outcome = outcome, species_df)
  
  # Determine mode
  if (!is.numeric(outcome)) {
    full_df$Outcome <- as.factor(full_df$Outcome)
    mode_type <- "classification"
	print("Classification")
    }
  else {
    mode_type <- "regression"
	print("Regression")
  }
  
  # -----------------------------
  # FEATURE SELECTION MODE
  selected_features <- NULL
  rfe_res <- NULL
  near_zero_variance_removed <- NULL
  
  if (do_feature_selection) 
    {# -----------------------------
    # 2. Near Zero Variance (ONLY HERE)
    
    nzv <- nearZeroVar(full_df[, -1])
    near_zero_variance_removed <- colnames(full_df)[nzv + 1]
    
    species_filtered <- full_df[, -1][, -nzv, drop = FALSE]
    filtered_df <- data.frame(Outcome = full_df$Outcome, species_filtered)
    
    # -----------------------------
    # 3. RFE Feature Selection
    set.seed(123)
    ctrl <- rfeControl(functions = rfFuncs,method = "cv",number = 5)
    
    rfe_res <- rfe(x = filtered_df[, -1], y = filtered_df$Outcome, sizes = rfe_sizes, rfeControl = ctrl)
    
    selected_features <- predictors(rfe_res)
    
    # Keep only selected features
    rf_input <- data.frame(Outcome = filtered_df$Outcome, filtered_df[, selected_features, drop = FALSE])
  }
  else
    {
    # NO FEATURE SELECTION; NO NZV FILTERING
    rf_input <- full_df
    }
  
  # -----------------------------
  # 4. Random Forest model
  set.seed(123)
  # rf_model <- randomForest(Outcome ~ ., data = rf_input,importance = TRUE, ntree = 1000)
  x_rf <- as.matrix(rf_input[, -1, drop = FALSE])
  y_rf <- rf_input$Outcome

  rf_model <- randomForest(
    x = x_rf,
    y = y_rf,
    importance = TRUE,
    ntree = 1000
  )
  # -----------------------------
  # 5. Performance
  if (mode_type == "regression") 
    {
    performance <- cor.test(rf_model$predicted, rf_input$Outcome)
    }
  else {
    performance <- confusionMatrix(rf_model$predicted, rf_input$Outcome)
  }
  
  # -----------------------------
  # 6. Importance
  
  importance_df <- as.data.frame(importance(rf_model))
  importance_df$Species <- rownames(importance_df)
  
  # -----------------------------
  # 7. Return results
  return(list(mode = mode_type, near_zero_variance_removed = near_zero_variance_removed, rfe_results = rfe_res, 
              selected_features = selected_features, rf_model = rf_model, importance_table = importance_df, performance = performance))
}

supragingival_results1 <- run_rf_feature_selection(
  species_df = supragingival_MeanFunctionalProfile_groups,
  metadata_df = Combined_Supragingival_Scores,
  outcome_col = "HAC_Score")




############ Sort the importance table from lowest to highest %IncMSE
supragingival_combined_importance_table1 <- supragingival_results1$importance_table
supragingival_combined_importance_table1 <- supragingival_combined_importance_table1[order(supragingival_combined_importance_table1$`%IncMSE`), ,drop = FALSE]

## Numbers of top features to test
feature_numbers <- seq(from = 250,to = 20,by = -5)

## Empty data frame to store results
supragingival_rf_feature_number_results <- data.frame(Number_of_features = integer(),Correlation = numeric(),P_value = numeric(),CI_lower = numeric(),CI_upper = numeric(),stringsAsFactors = FALSE)

## Run models
for (n_features in feature_numbers) {

  message("Running RF with top ", n_features, " features")

  ## Select top features based on initial RF importance
  selected_features <- tail(rownames(supragingival_combined_importance_table1),n_features)

  ## Set seed so that results are reproducible
  set.seed(123)

  ## Fit RF model
  rf_result <- run_rf_feature_selection(
    species_df = supragingival_MeanFunctionalProfile_groups[,selected_features,drop = FALSE],
    metadata_df = Combined_Supragingival_Scores,
    outcome_col = "HAC_Score"
  )

  ## Extract cor.test output
  performance_test <- rf_result$performance

  ## Store results
  supragingival_rf_feature_number_results <- rbind(
    supragingival_rf_feature_number_results,
    data.frame(
      Number_of_features = n_features,
      Correlation = unname(performance_test$estimate),
      P_value = performance_test$p.value,
      CI_lower = performance_test$conf.int[1],
      CI_upper = performance_test$conf.int[2]
    )
  )
}

# Convert the numbers into two decimal places for better visualization in the barplots
supragingival_rf_feature_number_results$Correlation <- round(supragingival_rf_feature_number_results$Correlation,digits = 2)
supragingival_rf_feature_number_results$CI_lower <- round(supragingival_rf_feature_number_results$CI_lower,digits = 2)
supragingival_rf_feature_number_results$CI_upper <- round(supragingival_rf_feature_number_results$CI_upper,digits = 2)



### Now plot the Correlation barplots for each feature set, then CI lower barplot and then CI upper barplot.
library(ggplot2)
supragingival_rf_feature_number_results$Number_of_features <- factor(supragingival_rf_feature_number_results$Number_of_features,levels = supragingival_rf_feature_number_results$Number_of_features)
write.csv(supragingival_rf_feature_number_results, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3Supragingival_RF_Correlation_FeatureSelection.csv", row.names = TRUE)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3Supragingival_RF_Correlation_Barplot_FeatureSelection.pdf",width = 15, height = 6)

ggplot(supragingival_rf_feature_number_results,
       aes(x = Number_of_features, y = Correlation)) +
  geom_col(fill = "#a164c9") +
  geom_text(aes(label = sprintf("%.2f", Correlation)),
            vjust = 0, size = 3.5, angle = 90, hjust = -0.2) +
  coord_cartesian(ylim = c(0, 0.85)) +
  labs(
    x = "Number of top features",
    y = "Correlation",
    title = "Correlation across different numbers of top features"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3supragingival_RF_CI_Lower_Barplot_FeatureSelection.pdf",
    width = 15, height = 6)

ggplot(supragingival_rf_feature_number_results,
       aes(x = Number_of_features, y = CI_lower)) +
  geom_col(fill = "#59A14F") +
  geom_text(aes(label = sprintf("%.2f", CI_lower)),
            vjust = -0.5, size = 3.5, angle = 90, hjust = -0.2) +
  coord_cartesian(ylim = c(0, 0.85)) +
  labs(
    x = "Number of top features",
    y = "CI lower",
    title = "Lower confidence interval across different numbers of top features"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

dev.off()



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3supragingival_RF_CI_Upper_Barplot_FeatureSelection.pdf",
    width = 15, height = 6)

ggplot(supragingival_rf_feature_number_results,
       aes(x = Number_of_features, y = CI_upper)) +
  geom_col(fill = "#E15759") +
  geom_text(aes(label = sprintf("%.2f", CI_upper)),
            vjust = -0.5, size = 3.5, angle = 90, hjust = -0.2) +
  coord_cartesian(ylim = c(0, 0.85)) +
  labs(
    x = "Number of top features",
    y = "CI upper",
    title = "Upper confidence interval across different numbers of top features"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

dev.off()



### Now run the model for 25 of this features
supragingival_results2 <- run_rf_feature_selection(
  species_df = supragingival_MeanFunctionalProfile_groups[, tail(rownames(supragingival_combined_importance_table1),25)],
  metadata_df = Combined_Supragingival_Scores,
  outcome_col = "HAC_Score")

supragingival_plot_df <- data.frame(Predicted = supragingival_results2$rf_model$predicted,Actual = supragingival_results2$rf_model$y)
write.csv(supragingival_plot_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3Supragingival_Predicted_vs_Actual_25Features.csv", row.names = TRUE)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3Supragingival_Pred_vs_Actual_25Features.pdf",width = 3.5, height = 4)
ggplot(supragingival_plot_df, aes(x = Predicted, y = Actual)) +
  geom_point(color = "#64c977", size = 3, alpha = 0.8) + 
  geom_smooth(method = "lm", color = "#a164c9", fill = "#ae85ca", alpha = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Predicted HAC Score",
    y = "Actual HAC Score",
    title = "Predicted vs Actual HAC Score (25Features)"
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.x = element_text(color = "#117A65", size = 12),
    axis.title.y = element_text(color = "#117A65", size = 12),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = "top"
  )

dev.off()




########## Now see the direction and significance of the top 25 features using the correlaion test
top_25_selected_features <- tail(rownames(supragingival_combined_importance_table1),25)
supragingival_MeanFunctionalProfile_groups_filt <- supragingival_MeanFunctionalProfile_groups[, top_25_selected_features]
HAC_scores <- Combined_Supragingival_Scores[rownames(supragingival_MeanFunctionalProfile_groups_filt), "HAC_Score"]

supragingival_correlation_results <- data.frame(
  Feature = character(),
  Correlation_Coefficient = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
for(feature in top_25_selected_features) {
  feature_values <- supragingival_MeanFunctionalProfile_groups_filt[, feature]
  correlation_test <- cor.test(feature_values, HAC_scores)
    supragingival_correlation_results <- rbind(supragingival_correlation_results, data.frame(
        Feature = feature,
        Correlation_Coefficient = correlation_test$estimate,
        P_value = correlation_test$p.value
    ))

}

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3supragingival_HAC_prediction_Workspace.RData")


#### Now plot this correlation results as volcano plot (manuscript ready figure): Estimates on x axis, log(pvalue) on y axis
supragingival_correlation_results$qval <- p.adjust(supragingival_correlation_results$P_value, method = "fdr")
supragingival_correlation_results$Significance <- ifelse(supragingival_correlation_results$qval <= 0.05,"Significant","Non-significant")

library(ggplot2)
library(ggrepel)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3Supragingival_Correlation_VolcanoPlot_25Features.pdf",width = 9, height = 6)
ggplot(supragingival_correlation_results, aes(x = Correlation_Coefficient, y = -log10(qval), color = Significance)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(data = supragingival_correlation_results[supragingival_correlation_results$qval <= 0.05, , drop = FALSE],aes(label = Feature),size = 3.5,max.overlaps = 30,show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05),linewidth = 1,linetype = "solid",color = "black") +
  geom_vline(xintercept = 0,linewidth = 1,linetype = "solid",color = "black") +
  scale_color_manual(values = c("Significant" = "#8b33ff","Non-significant" = "grey70")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Correlation Coefficient",
    y = "-log10(Q-value)",
    title = "Correlation of Top 25 Features with sHACK Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.x = element_text(color = "#117A65", size = 12),
    axis.title.y = element_text(color = "#117A65", size = 12),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    legend.position = "top"
      )
dev.off()





###### Now get the actual annotations for these groups.
Supragingival_Top25_Annotations <- annotation_df[annotation_df$Groups %in% top_25_selected_features, ]
rownames(Supragingival_Top25_Annotations) <- Supragingival_Top25_Annotations$Groups

rownames(supragingival_correlation_results) <- supragingival_correlation_results$Feature
supragingival_correlation_results <- supragingival_correlation_results[rownames(Supragingival_Top25_Annotations),]

supragingival_combined_results <- cbind(supragingival_correlation_results, Supragingival_Top25_Annotations)
supragingival_combined_results$Feature <- NULL

Supragingival_MeanFun_Groups = supragingival_MeanFunctionalProfile_groups[,tail(rownames(supragingival_combined_importance_table1), 25)]

# Add manual annotation for the feature "kegg_mod__M00178"
supragingival_combined_results$Annotations[supragingival_combined_results$Feature_ids == "kegg_mod__M00178"] <- "Encodes the structural components of the bacterial 70S ribosome (comprising 30S small and 50S large subunits)"
supragingival_combined_results$Annotations[supragingival_combined_results$Feature_ids == "kegg_mod__M00357"] <- "Methanogenesis"
supragingival_combined_results$Annotations[supragingival_combined_results$Feature_ids == "kegg_mod__M00579"] <- "Phosphate acetyltransferase-acetate kinase pathway"
supragingival_combined_results$Annotations[supragingival_combined_results$Feature_ids == "kegg_rclass__RC00318"] <- paste0("N-ribosyl hydrolase reaction class involved in cleavage of ","nucleosides and related ribosylated metabolites")
supragingival_combined_results$Annotations[supragingival_combined_results$Feature_ids == "kegg_rclass__RC00543"] <- paste0("Sugar-acid dehydration reaction class involved in the conversion of ","hydroxy sugar acids to deoxy-keto sugar acids")
supragingival_combined_results$Annotations[supragingival_combined_results$Feature_ids == "kegg_rclass__RC02799"] <- paste0("Ammonia-amino acid interconversion reaction class involving ","glutamate, aspartate, and related nitrogen metabolism")

## Now add the Manually > functional Theme
library(dplyr)

supragingival_combined_results <- supragingival_combined_results %>%
  mutate(
    Functional_Theme = case_when(

      grepl("Ribosomal protein|ribosomal protein|70S ribosome",Annotations,ignore.case = TRUE) ~ "Ribosome structure and translation",

      grepl("tRNA synthetase|elongation factor EF-G|Elongation factor G",Annotations,ignore.case = TRUE) ~ "Translation and protein synthesis",

      grepl("aminotransferase|amino acid interconversion|glutamate|aspartate.*nitrogen",Annotations,ignore.case = TRUE) ~ "Amino acid and nitrogen metabolism",

      grepl("Ammonia channel|Ammonium Transporter|ammonium transporter",Annotations,ignore.case = TRUE) ~ "Ion transport and nitrogen homeostasis",

      grepl("SemiSWEET|Sugar transporter|Sugar efflux transporter",Annotations,ignore.case = TRUE) ~ "Sugar transport and exchange",

      grepl("Sugar-acid dehydration",Annotations,ignore.case = TRUE) ~ "Carbohydrate metabolism",

      grepl("Phosphate acetyltransferase-acetate kinase",Annotations,ignore.case = TRUE) ~ "Acetate production and energy metabolism",

      grepl("Methanogenesis",Annotations,ignore.case = TRUE) ~ "Methane metabolism",

      grepl("phosphopantothenoylcysteine decarboxylase|phosphopantothenate---cysteine ligase",Annotations,ignore.case = TRUE) ~ "Coenzyme A biosynthesis",

      grepl("Enoyl-\\(Acyl carrier protein\\) reductase",Annotations,ignore.case = TRUE) ~ "Fatty acid biosynthesis",

      grepl("N-ribosyl hydrolase|nucleosides|ribosylated metabolites",Annotations,ignore.case = TRUE) ~ "Nucleotide and nucleoside metabolism",

      grepl("filamentous hemagglutinin",Annotations,ignore.case = TRUE) ~ "Cell adhesion and host interaction",

      grepl("Cyclopropanoid cyclopropyl hydrolase",Annotations,ignore.case = TRUE) ~ "Lipid modification and metabolism",

      grepl("GyrI-like",Annotations,ignore.case = TRUE) ~ "Small-molecule binding and stress-associated functions",

      grepl("DUF1706|DUF4304|unknown function|Uncharacterized conserved protein", Annotations,ignore.case = TRUE) ~ "Domains of unknown function",

      is.na(Annotations) | trimws(Annotations) == "" ~ "Unannotated or insufficiently characterized",

      TRUE ~ "Other functional category")
  )




write.csv(supragingival_combined_results, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3Supragingival_Correlation_Results_Top25Features.csv", row.names = TRUE)
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3supragingival_HAC_prediction_Workspace.RData")


### Now plot the volcano plot with functional themes as colours and significant features
library(ggplot2)
library(ggrepel)

## Distinct colours for functional themes
supra_theme_colors <- c(
  "Ribosome structure and translation"                     = "#E41A1C",
  "Translation and protein synthesis"                      = "#377EB8",
  "Amino acid and nitrogen metabolism"                     = "#4DAF4A",
  "Ion transport and nitrogen homeostasis"                 = "#984EA3",
  "Lipid modification and metabolism"                      = "#FF7F00",
  "Sugar transport and exchange"                           = "#F781BF",
  "Domains of unknown function"                            = "#666666",
  "Coenzyme A biosynthesis"                                = "#A65628",
  "Cell adhesion and host interaction"                     = "#00A6A6",
  "Methane metabolism"                                     = "#6A3D9A",
  "Acetate production and energy metabolism"               = "#1B9E77",
  "Nucleotide and nucleoside metabolism"                   = "#D95F02",
  "Carbohydrate metabolism"                                = "#7570B3",
  "Fatty acid biosynthesis"                                = "#66A61E",
  "Small-molecule binding and stress-associated functions" = "#E7298A"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3Supragingival_Correlation_VolcanoPlot_25Features_Summarised.pdf",width = 14,height = 8.2)
ggplot(
  supragingival_combined_results,
  aes(x = Correlation_Coefficient,y = -log10(qval),color = Functional_Theme)) +
  geom_point(size = 5, alpha = 0.85) +
  geom_text_repel(data = supragingival_combined_results,
    aes(label = Functional_Theme),size = 5,max.overlaps = Inf,box.padding = 0.5,point.padding = 0.3,min.segment.length = 0,show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05),linewidth = 0.8,linetype = "dashed",color = "black") +
  geom_vline(xintercept = 0,linewidth = 0.8,linetype = "solid",color = "black") +
  scale_color_manual(values = supra_theme_colors,na.value = "grey50") +
  guides(color = guide_legend(title = "Functional theme",nrow = 3,byrow = TRUE,override.aes = list(size = 5, alpha = 1))) +
  labs(x = "Correlation coefficient",y = expression(-log[10]("Q-value")),title = "Correlation of Top 25 Functional Features with HAC Score") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5,face = "bold",size = 12),
    axis.title.x = element_text(color = "#117A65",size = 12),
    axis.title.y = element_text(color = "#117A65",size = 12),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(face = "bold")
  )
dev.off()


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_3supragingival_HAC_prediction_Workspace.RData")
