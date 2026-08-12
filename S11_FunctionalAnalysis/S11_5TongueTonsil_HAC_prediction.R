
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4_AllScores_AllSubsites.RData")
load( "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2TongueTonsil_data_for_prediction.RData")


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

tonguetonsil_results1 <- run_rf_feature_selection(
  species_df = tongue_MeanFunctionalProfile_groups,
  metadata_df = Combined_TongueTonsil_Scores,
  outcome_col = "HAC_Score")




############ Sort the importance table from lowest to highest %IncMSE
tonguetonsil_combined_importance_table1 <- tonguetonsil_results1$importance_table
tonguetonsil_combined_importance_table1 <- tonguetonsil_combined_importance_table1[order(tonguetonsil_combined_importance_table1$`%IncMSE`), ,drop = FALSE]

## Numbers of top features to test
feature_numbers <- seq(from = 250,to = 20,by = -5)

## Empty data frame to store results
tonguetonsil_rf_feature_number_results <- data.frame(Number_of_features = integer(),Correlation = numeric(),P_value = numeric(),CI_lower = numeric(),CI_upper = numeric(),stringsAsFactors = FALSE)

## Run models
for (n_features in feature_numbers) {

  message("Running RF with top ", n_features, " features")

  ## Select top features based on initial RF importance
  selected_features <- tail(rownames(tonguetonsil_combined_importance_table1),n_features)

  ## Set seed so that results are reproducible
  set.seed(123)

  ## Fit RF model
  rf_result <- run_rf_feature_selection(
    species_df = tongue_MeanFunctionalProfile_groups[,selected_features,drop = FALSE],
    metadata_df = Combined_TongueTonsil_Scores,
    outcome_col = "HAC_Score"
  )

  ## Extract cor.test output
  performance_test <- rf_result$performance

  ## Store results
  tonguetonsil_rf_feature_number_results <- rbind(
    tonguetonsil_rf_feature_number_results,
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
tonguetonsil_rf_feature_number_results$Correlation <- round(tonguetonsil_rf_feature_number_results$Correlation,digits = 2)
tonguetonsil_rf_feature_number_results$CI_lower <- round(tonguetonsil_rf_feature_number_results$CI_lower,digits = 2)
tonguetonsil_rf_feature_number_results$CI_upper <- round(tonguetonsil_rf_feature_number_results$CI_upper,digits = 2)



### Now plot the Correlation barplots for each feature set, then CI lower barplot and then CI upper barplot.
library(ggplot2)
tonguetonsil_rf_feature_number_results$Number_of_features <- factor(tonguetonsil_rf_feature_number_results$Number_of_features,levels = tonguetonsil_rf_feature_number_results$Number_of_features)
write.csv(tonguetonsil_rf_feature_number_results, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_RF_Correlation_FeatureSelection.csv", row.names = TRUE)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_RF_Correlation_Barplot_FeatureSelection.pdf",width = 15, height = 6)
ggplot(tonguetonsil_rf_feature_number_results,
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


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_RF_CI_Lower_Barplot_FeatureSelection.pdf", width = 15, height = 6)
ggplot(tonguetonsil_rf_feature_number_results,
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


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_RF_CI_Upper_Barplot_FeatureSelection.pdf", width = 15, height = 6)
ggplot(tonguetonsil_rf_feature_number_results,
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


### Now run the model for 20 of this features
tonguetonsil_results2 <- run_rf_feature_selection(
  species_df = tongue_MeanFunctionalProfile_groups[, tail(rownames(tonguetonsil_combined_importance_table1),20)],
  metadata_df = Combined_TongueTonsil_Scores,
  outcome_col = "HAC_Score")

tonguetonsil_plot_df <- data.frame(Predicted = tonguetonsil_results2$rf_model$predicted,Actual = tonguetonsil_results2$rf_model$y)
write.csv(tonguetonsil_plot_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_Predicted_vs_Actual_20Features.csv", row.names = TRUE)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_Pred_vs_Actual_20Features.pdf",width = 3.5, height = 4)
ggplot(tonguetonsil_plot_df, aes(x = Predicted, y = Actual)) +
  geom_point(color = "#64c977", size = 3, alpha = 0.8) + 
  geom_smooth(method = "lm", color = "#a164c9", fill = "#ae85ca", alpha = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Predicted HAC Score",
    y = "Actual HAC Score",
    title = "Predicted vs Actual HAC Score (20Features)"
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


########## Now see the direction and significance of the top 20 features using the correlaion test
top_20_selected_features <- tail(rownames(tonguetonsil_combined_importance_table1),20)
tongue_MeanFunctionalProfile_groups_filt <- tongue_MeanFunctionalProfile_groups[, top_20_selected_features]
HAC_scores <- Combined_TongueTonsil_Scores[rownames(tongue_MeanFunctionalProfile_groups_filt), "HAC_Score"]

tonguetonsil_correlation_results <- data.frame(
  Feature = character(),
  Correlation_Coefficient = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
for(feature in top_20_selected_features) {
  feature_values <- tongue_MeanFunctionalProfile_groups_filt[, feature]
  correlation_test <- cor.test(feature_values, HAC_scores)
    tonguetonsil_correlation_results <- rbind(tonguetonsil_correlation_results, data.frame(
        Feature = feature,
        Correlation_Coefficient = correlation_test$estimate,
        P_value = correlation_test$p.value
    ))

}

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_HAC_prediction_Workspace.RData")


## Now plot this correlation results as volcano plot (manuscript ready figure): Estimates on x axis, log(pvalue) on y axis
tonguetonsil_correlation_results$qval <- p.adjust(tonguetonsil_correlation_results$P_value, method = "fdr")
tonguetonsil_correlation_results$Significance <- ifelse(tonguetonsil_correlation_results$qval <= 0.05,"Significant","Non-significant")

library(ggplot2)
library(ggrepel)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_Correlation_VolcanoPlot_20Features.pdf",width = 9, height = 6)
ggplot(tonguetonsil_correlation_results, aes(x = Correlation_Coefficient, y = -log10(qval), color = Significance)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(data = tonguetonsil_correlation_results[tonguetonsil_correlation_results$qval <= 0.05, , drop = FALSE],aes(label = Feature),size = 3.5,max.overlaps = 30,show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05),linewidth = 1,linetype = "solid",color = "black") +
  geom_vline(xintercept = 0,linewidth = 1,linetype = "solid",color = "black") +
  scale_color_manual(values = c("Significant" = "#8b33ff","Non-significant" = "grey70")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Correlation Coefficient",
    y = "-log10(Q-value)",
    title = "Correlation of Top 20 Features with HAC Score"
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
TongueTonsil_Top20_Annotations <- annotation_df[annotation_df$Groups %in% top_20_selected_features, ]
rownames(TongueTonsil_Top20_Annotations) <- TongueTonsil_Top20_Annotations$Groups

rownames(tonguetonsil_correlation_results) <- tonguetonsil_correlation_results$Feature
tonguetonsil_correlation_results <- tonguetonsil_correlation_results[rownames(TongueTonsil_Top20_Annotations),]

tonguetonsil_combined_results <- cbind(tonguetonsil_correlation_results, TongueTonsil_Top20_Annotations)
tonguetonsil_combined_results$Feature <- NULL

TongueTonsil_MeanFun_Groups = tongue_MeanFunctionalProfile_groups[,tail(rownames(tonguetonsil_combined_importance_table1), 20)]


###### Now add the Functional Theme column
library(dplyr)

tonguetonsil_combined_results <- tonguetonsil_combined_results %>%
  mutate(
    Functional_Theme = case_when(
      grepl("Ribosomal protein|ribosomal protein",Annotations,ignore.case = TRUE) ~ "Ribosome structure and translation",
      grepl("Uracil phosphoribosyltransferase",Annotations,ignore.case = TRUE) ~ "Pyrimidine salvage and nucleotide metabolism",
      grepl("SecG|preprotein translocase",Annotations,ignore.case = TRUE) ~ "Protein secretion and membrane translocation",
      grepl("ArsR|arsenate|arsenite|antimonite",Annotations,ignore.case = TRUE) ~ "Metal resistance and transcriptional regulation",
      grepl("DNA primase|helicase|MutH|mismatch repair",Annotations,ignore.case = TRUE) ~ "DNA replication and repair",
      grepl("IS66|transposase",Annotations,ignore.case = TRUE) ~ "Mobile genetic elements and DNA recombination",
      grepl("CHAT domain",Annotations,ignore.case = TRUE) ~ "Proteolysis and peptidase activity",
      grepl("Ankyrin repeat|Leucine rich repeat",Annotations,ignore.case = TRUE) ~ "Protein interaction and repeat domains",
      grepl("AAA domain|Predicted ATPase",Annotations,ignore.case = TRUE) ~ "ATPase and ATP-dependent functions",
      grepl("YgiM|DUF1202|Uncharacterized conserved protein",Annotations,ignore.case = TRUE) ~ "Domains of unknown function",
      is.na(Annotations) | trimws(Annotations) == "" ~ "Unannotated or insufficiently characterized",
      TRUE ~ "Other functional category")
  )


write.csv(tonguetonsil_combined_results, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_Correlation_Results_Top20Features.csv", row.names = TRUE)
save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_HAC_prediction_Workspace.RData")


### Now plot the volcano plot with functional themes as colours and significant features
library(ggplot2)
library(ggrepel)

## Distinct colours for functional themes
tongue_theme_colors <- c(
  "Pyrimidine salvage and nucleotide metabolism"    = "#D95F02",
  "Ribosome structure and translation"              = "#E41A1C",
  "Domains of unknown function"                     = "#666666",
  "ATPase and ATP-dependent functions"              = "#A6CEE3",
  "Protein interaction and repeat domains"          = "#A65628",
  "Protein secretion and membrane translocation"    = "#4DAF4A",
  "Metal resistance and transcriptional regulation" = "#FF7F00",
  "DNA replication and repair"                      = "#8C510A",
  "Proteolysis and peptidase activity"              = "#377EB8",
  "Mobile genetic elements and DNA recombination"   = "#6A3D9A"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_Correlation_VolcanoPlot_20Features_Summarised.pdf",width = 14,height = 8.2)
ggplot(
  tonguetonsil_combined_results,
  aes(x = Correlation_Coefficient,y = -log10(qval),color = Functional_Theme)) +
  geom_point(size = 5, alpha = 0.85) +
  geom_text_repel(data = tonguetonsil_combined_results,
    aes(label = Functional_Theme),size = 5,max.overlaps = Inf,box.padding = 0.5,point.padding = 0.3,min.segment.length = 0,show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05),linewidth = 0.8,linetype = "dashed",color = "black") +
  geom_vline(xintercept = 0,linewidth = 0.8,linetype = "solid",color = "black") +
  scale_color_manual(values = tongue_theme_colors,na.value = "grey50") +
  guides(color = guide_legend(title = "Functional theme",nrow = 3,byrow = TRUE,override.aes = list(size = 5, alpha = 1))) +
  labs(x = "Correlation coefficient",y = expression(-log[10]("Q-value")),title = "Correlation of Top 20 Functional Features with HAC Score") +
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


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_5TongueTonsil_HAC_prediction_Workspace.RData")

