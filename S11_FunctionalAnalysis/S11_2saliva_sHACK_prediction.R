### This script is to predict the HAC score and sHACK score using the functional features. Before that extract the functional profile for only the overlapping species.

load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_1MeanFunctional_ConservationProfile_OverallGroups.RData")
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S6_CombinedAssociationScores/S6_4_AllScores_AllSubsites.RData")


## extract the matched species for each subsite
saliva_matched_taxa <- intersect(rownames(Combined_Saliva_Scores), rownames(MeanFunctionalProfile_groups))
# 366 species off 499

supragingival_matched_taxa <- intersect(rownames(Combined_Supragingival_Scores), rownames(MeanFunctionalProfile_groups))
# 220 species off 301

subgingival_matched_taxa <- intersect(rownames(Combined_Subgingival_Scores), rownames(MeanFunctionalProfile_groups))
# 140 species off 196

tongue_matched_taxa <- intersect(rownames(Combined_TongueTonsil_Scores), rownames(MeanFunctionalProfile_groups))
# 190 species off 266

overlapped_species <- unique(c(saliva_matched_taxa,supragingival_matched_taxa,subgingival_matched_taxa,tongue_matched_taxa))
write.csv(overlapped_species, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Overlapped_species_All.csv")

######## Now extract these matched species from the functional profile for each subsite and start prediction
saliva_MeanFunctionalProfile_groups <- MeanFunctionalProfile_groups[saliva_matched_taxa, ]
saliva_MeanFunctionalProfile_groups <- saliva_MeanFunctionalProfile_groups[,colSums(saliva_MeanFunctionalProfile_groups)>0]

supragingival_MeanFunctionalProfile_groups <- MeanFunctionalProfile_groups[supragingival_matched_taxa, ]
supragingival_MeanFunctionalProfile_groups <- supragingival_MeanFunctionalProfile_groups[,colSums(supragingival_MeanFunctionalProfile_groups)>0]

subgingival_MeanFunctionalProfile_groups <- MeanFunctionalProfile_groups[subgingival_matched_taxa, ]
subgingival_MeanFunctionalProfile_groups <- subgingival_MeanFunctionalProfile_groups[,colSums(subgingival_MeanFunctionalProfile_groups)>0]

tongue_MeanFunctionalProfile_groups <- MeanFunctionalProfile_groups[tongue_matched_taxa, ]
tongue_MeanFunctionalProfile_groups <- tongue_MeanFunctionalProfile_groups[,colSums(tongue_MeanFunctionalProfile_groups)>0]

write.csv(saliva_MeanFunctionalProfile_groups, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_MeanFunctionalProfile.csv")
write.csv(supragingival_MeanFunctionalProfile_groups, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Supragingival_MeanFunctionalProfile.csv")
write.csv(subgingival_MeanFunctionalProfile_groups, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Subgingival_MeanFunctionalProfile.csv")
write.csv(tongue_MeanFunctionalProfile_groups, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2TongueTonsil_MeanFunctionalProfile.csv")

save(supragingival_matched_taxa,supragingival_MeanFunctionalProfile_groups, annotation_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Supragingival_data_for_prediction.RData")
save(subgingival_matched_taxa,subgingival_MeanFunctionalProfile_groups, annotation_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Subgingival_data_for_prediction.RData")
save(tongue_matched_taxa,tongue_MeanFunctionalProfile_groups, annotation_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2TongueTonsil_data_for_prediction.RData")


rm(MeanFunctionalProfile_groups)


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2saliva_sHACK_prediction_Workspace.RData")



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



saliva_results1 <- run_rf_feature_selection(
  species_df = saliva_MeanFunctionalProfile_groups,
  metadata_df = Combined_Saliva_Scores,
  outcome_col = "HACK_Score")


# saliva_combined_importance_table2 <- saliva_results2$importance_table
# saliva_combined_importance_table2 <- saliva_combined_importance_table2[order(saliva_combined_importance_table2$`%IncMSE`), ]
# features_selected2 <- tail(rownames(saliva_combined_importance_table2),250)

# saliva_combined_importance_table3 <- saliva_results3$importance_table
# saliva_combined_importance_table3 <- saliva_combined_importance_table3[order(saliva_combined_importance_table3$`%IncMSE`), ]
# features_selected3 <- tail(rownames(saliva_combined_importance_table3),250)

# # Find common features
# common_features <- intersect(intersect(features_selected1, features_selected2), features_selected3)

# saliva_results2 <- run_rf_feature_selection(
#   species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),250)],
#   metadata_df = Combined_Saliva_Scores,
#   outcome_col = "HACK_Score")

# saliva_results2$performance
# #         Pearson's product-moment correlation

# # data:  rf_model$predicted and rf_input$Outcome
# # t = 19.021, df = 364, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.6506197 0.7539802
# # sample estimates:
# #       cor 
# # 0.7060404 

# saliva_results3 <- run_rf_feature_selection(
#   species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),200)],
#   metadata_df = Combined_Saliva_Scores,
#   outcome_col = "HACK_Score")

# saliva_results3$performance
# #         Pearson's product-moment correlation

# # data:  rf_model$predicted and rf_input$Outcome
# # t = 19.085, df = 364, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.6519770 0.7549956
# # sample estimates:
# #       cor 
# # 0.7072206


# saliva_results4 <- run_rf_feature_selection(
#   species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),150)],
#   metadata_df = Combined_Saliva_Scores,
#   outcome_col = "HACK_Score")

# saliva_results4$performance
# #         Pearson's product-moment correlation

# # data:  rf_model$predicted and rf_input$Outcome
# # t = 19.528, df = 364, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.6612681 0.7619327
# # sample estimates:
# #      cor 
# # 0.715291

# saliva_results5 <- run_rf_feature_selection(
#   species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),100)],
#   metadata_df = Combined_Saliva_Scores,
#   outcome_col = "HACK_Score")

# saliva_results5$performance
# #         Pearson's product-moment correlation

# # data:  rf_model$predicted and rf_input$Outcome
# # t = 19.554, df = 364, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.6617977 0.7623275
# # sample estimates:
# #       cor 
# # 0.7157506

# saliva_results6 <- run_rf_feature_selection(
#   species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),75)],
#   metadata_df = Combined_Saliva_Scores,
#   outcome_col = "HACK_Score")

# saliva_results6$performance
# #         Pearson's product-moment correlation

# # data:  rf_model$predicted and rf_input$Outcome
# # t = 19.461, df = 364, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.6598681 0.7608889
# # sample estimates:
# #       cor 
# # 0.7140759

# saliva_results7 <- run_rf_feature_selection(
#   species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),50)],
#   metadata_df = Combined_Saliva_Scores,
#   outcome_col = "HACK_Score")

# saliva_results7$performance
# #         Pearson's product-moment correlation

# # data:  rf_model$predicted and rf_input$Outcome
# # t = 19.066, df = 364, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.6515647 0.7546872
# # sample estimates:
# #       cor 
# # 0.7068621

# saliva_results8 <- run_rf_feature_selection(
#   species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),25)],
#   metadata_df = Combined_Saliva_Scores,
#   outcome_col = "HACK_Score")

# saliva_results8$performance
# #         Pearson's product-moment correlation

# # data:  rf_model$predicted and rf_input$Outcome
# # t = 17.299, df = 364, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.6112954 0.7243467
# # sample estimates:
# #       cor 
# # 0.6717133


# saliva_results9 <- run_rf_feature_selection(
#   species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),70)],
#   metadata_df = Combined_Saliva_Scores,
#   outcome_col = "HACK_Score")

# saliva_results9$performance
# #         Pearson's product-moment correlation

# # data:  rf_model$predicted and rf_input$Outcome
# # t = 19.463, df = 364, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.6599218 0.7609289
# # sample estimates:
# #       cor 
# # 0.7141224

############ Sort the importance table from lowest to highest %IncMSE
saliva_combined_importance_table1 <- saliva_results1$importance_table
saliva_combined_importance_table1 <- saliva_combined_importance_table1[order(saliva_combined_importance_table1$`%IncMSE`), ,drop = FALSE]

## Numbers of top features to test
feature_numbers <- seq(from = 250,to = 20,by = -5)

## Empty data frame to store results
saliva_rf_feature_number_results <- data.frame(Number_of_features = integer(),Correlation = numeric(),P_value = numeric(),CI_lower = numeric(),CI_upper = numeric(),stringsAsFactors = FALSE)

## Run models
for (n_features in feature_numbers) {

  message("Running RF with top ", n_features, " features")

  ## Select top features based on initial RF importance
  selected_features <- tail(rownames(saliva_combined_importance_table1),n_features)

  ## Set seed so that results are reproducible
  set.seed(123)

  ## Fit RF model
  rf_result <- run_rf_feature_selection(
    species_df = saliva_MeanFunctionalProfile_groups[,selected_features,drop = FALSE],
    metadata_df = Combined_Saliva_Scores,
    outcome_col = "HACK_Score"
  )

  ## Extract cor.test output
  performance_test <- rf_result$performance

  ## Store results
  saliva_rf_feature_number_results <- rbind(
    saliva_rf_feature_number_results,
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
saliva_rf_feature_number_results$Correlation <- round(saliva_rf_feature_number_results$Correlation,digits = 2)
saliva_rf_feature_number_results$CI_lower <- round(saliva_rf_feature_number_results$CI_lower,digits = 2)
saliva_rf_feature_number_results$CI_upper <- round(saliva_rf_feature_number_results$CI_upper,digits = 2)




### Now plot the Correlation barplots for each feature set, then CI lower barplot and then CI upper barplot.
library(ggplot2)
saliva_rf_feature_number_results$Number_of_features <- factor(saliva_rf_feature_number_results$Number_of_features,levels = saliva_rf_feature_number_results$Number_of_features)
write.csv(saliva_rf_feature_number_results, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_RF_Correlation_FeatureSelection.csv",row.names = FALSE)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_RF_Correlation_Barplot_FeatureSelection.pdf",width = 15, height = 6)

ggplot(saliva_rf_feature_number_results,
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


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_RF_CI_Lower_Barplot_FeatureSelection.pdf",
    width = 15, height = 6)

ggplot(saliva_rf_feature_number_results,
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



pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_RF_CI_Upper_Barplot_FeatureSelection.pdf",
    width = 15, height = 6)

ggplot(saliva_rf_feature_number_results,
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


### Now run the model for 40 of this features
saliva_results2 <- run_rf_feature_selection(
  species_df = saliva_MeanFunctionalProfile_groups[, tail(rownames(saliva_combined_importance_table1),40)],
  metadata_df = Combined_Saliva_Scores,
  outcome_col = "HACK_Score")

saliva_plot_df <- data.frame(Predicted = saliva_results2$rf_model$predicted,Actual = saliva_results2$rf_model$y)
write.csv(saliva_plot_df, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_Predicted_vs_Actual_40Features.csv",row.names = TRUE)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_Pred_vs_Actual_40Features.pdf",width = 3.5, height = 4)
ggplot(saliva_plot_df, aes(x = Predicted, y = Actual)) +
  geom_point(color = "#64c977", size = 3, alpha = 0.8) + 
  geom_smooth(method = "lm", color = "#a164c9", fill = "#ae85ca", alpha = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Predicted sHACK Score",
    y = "Actual sHACK Score",
    title = "Predicted vs Actual sHACK Score (40Features)"
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



########## Now see the direction and significance of the top 40 features using the correlaion test
top_40_selected_features <- tail(rownames(saliva_combined_importance_table1),40)
saliva_MeanFunctionalProfile_groups_filt <- saliva_MeanFunctionalProfile_groups[, top_40_selected_features]
sHACK_scores <- Combined_Saliva_Scores[rownames(saliva_MeanFunctionalProfile_groups_filt), "HACK_Score"]

saliva_correlation_results <- data.frame(
  Feature = character(),
  Correlation_Coefficient = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
for(feature in top_40_selected_features) {
  feature_values <- saliva_MeanFunctionalProfile_groups_filt[, feature]
  correlation_test <- cor.test(feature_values, sHACK_scores)
    saliva_correlation_results <- rbind(saliva_correlation_results, data.frame(
        Feature = feature,
        Correlation_Coefficient = correlation_test$estimate,
        P_value = correlation_test$p.value
    ))

}

save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2saliva_sHACK_prediction_Workspace.RData")


## Now plot this correlation results as volcano plot (manuscript ready figure): Estimates on x axis, log(qvalue) on y axis
saliva_correlation_results$qvalue <- p.adjust(saliva_correlation_results$P_value, method = "fdr")
saliva_correlation_results$Significance <- ifelse(saliva_correlation_results$qvalue <= 0.05,"Significant","Non-significant")

library(ggplot2)
library(ggrepel)
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_Correlation_VolcanoPlot_40Features.pdf",width = 9, height = 6)
ggplot(saliva_correlation_results, aes(x = Correlation_Coefficient, y = -log10(qvalue), color = Significance)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(data = saliva_correlation_results[saliva_correlation_results$qvalue <= 0.05, , drop = FALSE],aes(label = Feature),size = 3.5,max.overlaps = 30,show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05),linewidth = 1,linetype = "solid",color = "black") +
  geom_vline(xintercept = 0,linewidth = 1,linetype = "solid",color = "black") +
  scale_color_manual(values = c("Significant" = "#8b33ff","Non-significant" = "grey70")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Correlation Coefficient",
    y = "-log10(Q-value)",
    title = "Correlation of Top 40 Features with sHACK Score"
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
Saliva_Top40_Annotations <- annotation_df[annotation_df$Groups %in% top_40_selected_features, ]
rownames(Saliva_Top40_Annotations) <- Saliva_Top40_Annotations$Groups

rownames(saliva_correlation_results) <- saliva_correlation_results$Feature
saliva_correlation_results <- saliva_correlation_results[rownames(Saliva_Top40_Annotations),]

saliva_combined_results <- cbind(saliva_correlation_results, Saliva_Top40_Annotations)
saliva_combined_results$Feature <- NULL

Saliva_MeanFun_Groups = saliva_MeanFunctionalProfile_groups[,tail(rownames(saliva_combined_importance_table1), 40)]


##################### Lets group these functional groups into their respective categories or functional themes.
library(dplyr)

saliva_combined_results <- saliva_combined_results %>%
  mutate(
    Functional_Theme = case_when(

      ## Phosphate transport and regulation
      grepl("phosphate|PhoU|PhoR|PhoB",Annotations,ignore.case = TRUE) ~ "Phosphate acquisition and regulation",

      ## Magnesium and ammonium transport
      grepl("magnesium|Mg2\\+|ammonium transporter",Annotations,ignore.case = TRUE) ~ "Ion transport and homeostasis",

      ## Peptidases, amidases and hydrolases
      grepl("peptidase|hydrolase|amidohydrolase|amidase|acylase",Annotations,ignore.case = TRUE) ~ "Proteolysis and hydrolysis",

      ## Sugar transport
      grepl("sugar transporter|sugar efflux|SemiSWEET",Annotations,ignore.case = TRUE) ~ "Sugar transport and exchange",

      ## Tellurite-related functions
      grepl("tellurite|TehB",Annotations,ignore.case = TRUE) ~ "Stress resistance and methylation",

      ## DNA replication-associated function
      grepl("DNA polymerase processivity factor",Annotations,ignore.case = TRUE) ~ "DNA repair and replication",

      ## Membrane, cell wall and cell-surface functions
      grepl("flippase|LysM",Annotations,ignore.case = TRUE) ~ "Cell envelope and surface-associated functions",

      ## Signal transduction
      grepl("histidine kinase|WalK",Annotations,ignore.case = TRUE) ~ "Signal transduction and transcriptional regulation",

      ## Protein interaction
      grepl("ankyrin repeat",Annotations,ignore.case = TRUE) ~ "Protein interaction domains",

      ## Mobile genetic elements
      grepl("transposase",Annotations,ignore.case = TRUE) ~ "Mobile genetic elements",

      ## Oxidoreductase functions
      grepl("aldo/keto reductase",Annotations,ignore.case = TRUE) ~ "Oxidoreductase activity",

      ## ATPase-associated domain
      grepl("AAA domain",Annotations,ignore.case = TRUE) ~ "ATPase and ATP-dependent functions",

      ## Translation
      grepl("tRNA synthetases",Annotations,ignore.case = TRUE) ~ "Translation and tRNA metabolism",

      ## Cell-envelope stress
      grepl("PspC",Annotations,ignore.case = TRUE) ~ "Cell-envelope stress response",

      ## Domains of unknown function
      grepl("domain of unknown function|DUF",Annotations,ignore.case = TRUE) ~ "Domains of unknown function",

      ## Missing annotations
      is.na(Annotations) |trimws(Annotations) == "" |grepl("No data found",Annotations,ignore.case = TRUE) ~ "Unannotated or insufficiently characterized",

      ## Backup category
      TRUE ~ "Other functional category")
  )

### for this Unannotated or insufficiently characterized   find through manual search into databases and write its annotations
#            Total_features          Feature_ids   Annotations
# Group11386              1        ec__1.11.1.15 No data found
# Group11782              1         ec__3.6.3.27 No data found
# Group21479              1     kegg_mod__M00222              
# Group21520              1     kegg_mod__M00434              
# Group22060              1 kegg_rclass__RC00006              
# Group22122              1 kegg_rclass__RC00105  


saliva_combined_results$Annotations[saliva_combined_results$Feature_ids == "ec__1.11.1.15"] <- "Peroxiredoxin/thiol peroxidase; peroxide detoxification"
saliva_combined_results$Functional_Theme[saliva_combined_results$Feature_ids == "ec__1.11.1.15"] <- "Oxidative stress response and peroxide detoxification"

saliva_combined_results$Annotations[saliva_combined_results$Feature_ids == "ec__3.6.3.27"] <- paste0("ABC-type phosphate transporter; phosphate ABC transporter; ","phosphate-transporting ATPase (ambiguous)")
saliva_combined_results$Functional_Theme[saliva_combined_results$Feature_ids == "ec__3.6.3.27"] <- "Phosphate acquisition and regulation"

saliva_combined_results$Annotations[saliva_combined_results$Feature_ids == "kegg_mod__M00222"] <- "Phosphate transport system"
saliva_combined_results$Functional_Theme[saliva_combined_results$Feature_ids == "kegg_mod__M00222"] <- "Phosphate acquisition and regulation"

saliva_combined_results$Annotations[saliva_combined_results$Feature_ids == "kegg_mod__M00434"] <- "PhoR-PhoB (phosphate starvation response) two-component regulatory system"
saliva_combined_results$Functional_Theme[saliva_combined_results$Feature_ids == "kegg_mod__M00434"] <- "Phosphate acquisition and regulation"

saliva_combined_results$Annotations[saliva_combined_results$Feature_ids == "kegg_rclass__RC00006"] <- paste0("Amino acid-2-oxocarboxylic acid interconversion reaction class, ","including reductive amination and conversion of ","2-oxocarboxylic acids to CoA derivatives")
saliva_combined_results$Functional_Theme[saliva_combined_results$Feature_ids == "kegg_rclass__RC00006"] <- "Amino acid and 2-oxocarboxylic acid metabolism"

saliva_combined_results$Annotations[saliva_combined_results$Feature_ids == "kegg_rclass__RC00105"] <- paste0("Malic enzyme-associated oxidative decarboxylation of ","malate to pyruvate and CO2")
saliva_combined_results$Functional_Theme[saliva_combined_results$Feature_ids == "kegg_rclass__RC00105"] <- "Pyruvate and carbon metabolism"



save(saliva_combined_results, saliva_results2, Saliva_MeanFun_Groups, Combined_Saliva_Scores, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_Correlation_Results_Top40Features.RData")
write.csv(saliva_combined_results, file = "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_Correlation_Results_Top40Features.csv", row.names = TRUE)



### Now plot the volcano plot with functional themes as colours and significant features
library(ggplot2)
library(ggrepel)

## Distinct colours for functional themes
theme_colors <- c(
  "Cell envelope and surface-associated functions"        = "#E41A1C",
  "Phosphate acquisition and regulation"                  = "#377EB8",
  "Proteolysis and hydrolysis"                            = "#4DAF4A",
  "Signal transduction and transcriptional regulation"   = "#984EA3",
  "Stress resistance and methylation"                    = "#FF7F00",
  "Protein interaction domains"                          = "#A65628",
  "Sugar transport and exchange"                         = "#F781BF",
  "Oxidative stress response and peroxide detoxification"= "#00A6A6",
  "Ion transport and homeostasis"                        = "#FFD700",
  "Mobile genetic elements"                              = "#6A3D9A",
  "Amino acid and 2-oxocarboxylic acid metabolism"       = "#1B9E77",
  "Pyruvate and carbon metabolism"                       = "#D95F02",
  "Oxidoreductase activity"                              = "#7570B3",
  "Cell-envelope stress response"                        = "#66A61E",
  "Translation and tRNA metabolism"                      = "#E7298A",
  "DNA repair and replication"                           = "#8C510A",
  "ATPase and ATP-dependent functions"                   = "#A6CEE3",
  "Domains of unknown function"                          = "#666666"
)

pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2Saliva_Correlation_VolcanoPlot_40Features_Summarised.pdf",width = 14,height = 8.2)
ggplot(
  saliva_combined_results,
  aes(x = Correlation_Coefficient,y = -log10(qvalue),color = Functional_Theme)) +
  geom_point(size = 5, alpha = 0.85) +
  geom_text_repel(data = saliva_combined_results[saliva_combined_results$qvalue <= 0.05, ,drop = FALSE],
    aes(label = Functional_Theme),size = 5,max.overlaps = Inf,box.padding = 0.5,point.padding = 0.3,min.segment.length = 0,show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05),linewidth = 0.8,linetype = "dashed",color = "black") +
  geom_vline(xintercept = 0,linewidth = 0.8,linetype = "solid",color = "black") +
  scale_color_manual(values = theme_colors,na.value = "grey50") +
  guides(color = guide_legend(title = "Functional theme",nrow = 3,byrow = TRUE,override.aes = list(size = 5, alpha = 1))) +
  labs(x = "Correlation coefficient",y = expression(-log[10]("Q-value")),title = "Correlation of Top 40 Functional Features with sHACK Score") +
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



save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2saliva_sHACK_prediction_Workspace.RData")





# ############### Here Correlation cannot tell whether the the correlation is due to all species or only high HAC species or low species i.e it do not consider the absence of the feature in in any species.so we will do logistic regression to see the direction of the association and then we will do correlation (of only present function for a species) to see significance.
# saliva_two_part_results <- data.frame()

# for (feature in top_40_selected_features) {

#   feature_values <- saliva_MeanFunctionalProfile_groups_filt[, feature]

#   complete_cases <- complete.cases(feature_values,sHACK_scores)
#   feature_values_use <- feature_values[complete_cases]
#   sHACK_scores_use <- sHACK_scores[complete_cases]

#   # Part 1: Presence/absence

#   feature_presence <- as.integer(feature_values_use > 0)
#   detection_estimate <- NA
#   detection_pvalue <- NA

#   if (length(unique(feature_presence)) > 1) {

#     detection_model <- glm(feature_presence ~ sHACK_scores_use,family = binomial())

#     detection_summary <- summary(detection_model)$coefficients
#     detection_estimate <- detection_summary["sHACK_scores_use","Estimate"]

#     detection_pvalue <- detection_summary["sHACK_scores_use","Pr(>|z|)"]
#   }

#   # Part 2: Abundance among detected

#   detected <- feature_values_use > 0

#   abundance_rho <- NA
#   abundance_pvalue <- NA

#   if (
#     sum(detected) >= 3 &&
#     length(unique(feature_values_use[detected])) > 1
#   ) {

#     abundance_test <- cor.test(feature_values_use[detected],sHACK_scores_use[detected],method = "pearson",exact = FALSE)

#     abundance_rho <- unname(abundance_test$estimate)
#     abundance_pvalue <- unname(abundance_test$p.value)
#   }

#   saliva_two_part_results <- rbind(
#     saliva_two_part_results,
#     data.frame(
#       Feature = feature,
#       Number_Detected = sum(detected),
#       Detection_Frequency = mean(detected),
#       Detection_Estimate = detection_estimate,
#       Detection_P_value = detection_pvalue,
#       Abundance_Rho = abundance_rho,
#       Abundance_P_value = abundance_pvalue
#     )
#   )
# }

# # This saliva two part results will not be used. Its just for our understanding.


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2saliva_sHACK_prediction_Workspace.RData")



###########################
###################### Extended Functional Analysis
load("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S8_AssociationConsistency_PhyloTrees/S8_1Saliva_PCAdf.RData")

high_sHACK_sp <- rownames(df_pca_saliva[df_pca_saliva$sHACK >= 0.90,])
lower_sHACK_sp <- rownames(df_pca_saliva[df_pca_saliva$sHACK < 0.90,])
module2_high_sHACK_sp <- rownames(df_pca_saliva[df_pca_saliva$sHACK >= 0.90 & df_pca_saliva$cluster == 2,])
module2_sp <- rownames(df_pca_saliva[df_pca_saliva$cluster == 2,])
module3_sp <- rownames(df_pca_saliva[df_pca_saliva$cluster == 3,])



############################################################
## Analysis 1:
## Compare high-sHACK vs lower-sHACK taxa for top 40 features
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)

## Output directory
out_dir <- "/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis"

## Saliva_MeanFun_Groups should already contain the selected top 40 features
## rows = species, columns = functional features
functional_df <- Saliva_MeanFun_Groups

## Keep only species present in both functional data and df_pca_saliva
common_species <- intersect(rownames(functional_df), rownames(df_pca_saliva))

functional_df <- functional_df[common_species, , drop = FALSE]
saliva_meta_for_func <- df_pca_saliva[common_species, , drop = FALSE]

## Convert feature values to numeric safely
functional_df <- as.data.frame(functional_df)
functional_df[] <- lapply(functional_df, function(x) as.numeric(as.character(x)))

## Define high-sHACK and lower-sHACK groups
saliva_meta_for_func$sHACK_group <- ifelse(saliva_meta_for_func$sHACK >= 0.90,"High_sHACK","Lower_sHACK")
saliva_meta_for_func$sHACK_group <- factor(saliva_meta_for_func$sHACK_group,levels = c("High_sHACK", "Lower_sHACK"))

table(saliva_meta_for_func$sHACK_group) # 43 are high HACK and 323 are lower HACK

############################################################
## Feature-wise Wilcoxon test
############################################################

high_vs_lower_results <- data.frame()

for(feature in colnames(functional_df)) {
  
  test_df <- data.frame(
    Species = rownames(functional_df),
    Feature = feature,
    Feature_Value = functional_df[, feature],
    sHACK = saliva_meta_for_func[rownames(functional_df), "sHACK"],
    sHACK_group = saliva_meta_for_func[rownames(functional_df), "sHACK_group"]
  )
  
  test_df <- test_df[complete.cases(test_df), ]
  
  ## Run only if both groups are present
  if(length(unique(test_df$sHACK_group)) == 2) {
    
    wilcox_res <- wilcox.test(Feature_Value ~ sHACK_group,data = test_df,exact = FALSE)
    
    high_values <- test_df$Feature_Value[test_df$sHACK_group == "High_sHACK"]
    lower_values <- test_df$Feature_Value[test_df$sHACK_group == "Lower_sHACK"]
    
    high_vs_lower_results <- rbind(
      high_vs_lower_results,
      data.frame(
        Feature = feature,
        N_High_sHACK = length(high_values),
        N_Lower_sHACK = length(lower_values),
        Mean_High_sHACK = mean(high_values, na.rm = TRUE),
        Mean_Lower_sHACK = mean(lower_values, na.rm = TRUE),
        Median_High_sHACK = median(high_values, na.rm = TRUE),
        Median_Lower_sHACK = median(lower_values, na.rm = TRUE),
        Mean_Difference_High_minus_Lower =
          mean(high_values, na.rm = TRUE) - mean(lower_values, na.rm = TRUE),
        Median_Difference_High_minus_Lower =
          median(high_values, na.rm = TRUE) - median(lower_values, na.rm = TRUE),
        Wilcox_P_value = wilcox_res$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}

## FDR correction
high_vs_lower_results$Q_value <- p.adjust(high_vs_lower_results$Wilcox_P_value,method = "fdr")
## Direction
high_vs_lower_results$Direction <- ifelse(high_vs_lower_results$Mean_Difference_High_minus_Lower > 0,"Higher in High_sHACK","Higher in Lower_sHACK")

## Add annotations and functional themes if saliva_combined_results exists
if(exists("saliva_combined_results")) {
  
  annotation_add <- saliva_combined_results %>%
    rownames_to_column("Feature") %>%
    dplyr::select(
      Feature,
      Feature_ids,
      Annotations,
      Functional_Theme,
      Correlation_Coefficient,
      qvalue,
      Significance
    )
  
  high_vs_lower_results <- high_vs_lower_results %>%
    left_join(annotation_add, by = "Feature")
}

## Sort by FDR
high_vs_lower_results <- high_vs_lower_results %>% arrange(Q_value)

## Save result
write.csv(high_vs_lower_results,file = file.path(out_dir, "S11_2Saliva_Top40_High_sHACK_vs_Lower_sHACK_WilcoxResults.csv"),row.names = FALSE)

head(high_vs_lower_results)



############################################################
## Boxplot of all top 40 features
############################################################

plot_df <- functional_df %>%
  rownames_to_column("Species") %>%
  pivot_longer(cols = -Species,names_to = "Feature",values_to = "Feature_Value") %>%
  left_join(saliva_meta_for_func %>% rownames_to_column("Species") %>% dplyr::select(Species, sHACK, sHACK_group), by = "Species")

## Add functional theme for cleaner facet labels if available
if(exists("saliva_combined_results")) {
  
  theme_df <- saliva_combined_results %>%
    rownames_to_column("Feature") %>%
    dplyr::select(Feature, Functional_Theme, Annotations)
  
  plot_df <- plot_df %>%
    left_join(theme_df, by = "Feature")
  
  plot_df$Feature_Label <- paste0(plot_df$Functional_Theme, "\n", plot_df$Feature)
  
} else {
  plot_df$Feature_Label <- plot_df$Feature
}

pdf(file.path(out_dir, "S11_2Saliva_Top40_High_sHACK_vs_Lower_sHACK_Boxplots_AllFeatures.pdf"),width = 8,height = 14)
ggplot(plot_df, aes(x = sHACK_group, y = Feature_Value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.45) +
  facet_wrap(~ Feature_Label, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 12) +
  labs(
    x = "",
    y = "Genome-derived functional feature value",
    title = "Top 40 salivary functional features in high-sHACK versus lower-sHACK taxa"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

dev.off()




##############  Plot the box plot for each of the species set as below:
pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_1FeaturePlots/S11_2Saliva_Group24412_Boxplot.pdf",width = 4,height = 5)
boxplot(
  functional_df[intersect(high_sHACK_sp, rownames(functional_df)), 1],
  functional_df[intersect(lower_sHACK_sp, rownames(functional_df)), 1],
  functional_df[intersect(module2_high_sHACK_sp, rownames(functional_df)), 1],
  functional_df[intersect(module2_sp, rownames(functional_df)), 1],
  functional_df[intersect(module3_sp, rownames(functional_df)), 1],
  names = c("High sHACK","Low sHACK","High sHACK\nModule 2","Module 2","Module 3"),ylab = "Group24412",xlab = "",main = "Translation and tRNA metabolism\ntRNA synthetases",las = 2, outline = F, col = c("#2F6B9A", "#B85BCB", "#2E8B57", "#76C56E", "#D99A4E"),border = "black",cex.axis = 0.8)
dev.off()


pdf("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_1FeaturePlots/S11_2Saliva_Group22122_Boxplot.pdf",width = 4,height = 5)
boxplot(
  functional_df[intersect(high_sHACK_sp, rownames(functional_df)), 2],
  functional_df[intersect(lower_sHACK_sp, rownames(functional_df)), 2],
  functional_df[intersect(module2_high_sHACK_sp, rownames(functional_df)), 2],
  functional_df[intersect(module2_sp, rownames(functional_df)), 2],
  functional_df[intersect(module3_sp, rownames(functional_df)), 2],
  names = c("High sHACK","Low sHACK","High sHACK\nModule 2","Module 2","Module 3"),ylab = "Group22122",xlab = "",main = "",las = 2, outline = F, col = c("#2F6B9A", "#B85BCB", "#2E8B57", "#76C56E", "#D99A4E"),border = "black",cex.axis = 0.8)
dev.off()



plot_functional_boxplot <- function(col_num, functional_theme) {
  feature_name <- names(functional_df)[col_num]
  pdf(paste0("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_1FeaturePlots/","S11_2Saliva_", feature_name,"_Boxplot.pdf"),width = 4,height = 5)
  boxplot(
    functional_df[intersect(high_sHACK_sp, rownames(functional_df)), col_num],
    functional_df[intersect(lower_sHACK_sp, rownames(functional_df)), col_num],
    functional_df[intersect(module2_high_sHACK_sp, rownames(functional_df)), col_num],
    functional_df[intersect(module2_sp, rownames(functional_df)), col_num],
    functional_df[intersect(module3_sp, rownames(functional_df)), col_num],
    names = c("High sHACK","Low sHACK","High sHACK\nModule 2","Module 2","Module 3"),ylab = feature_name,xlab = "",main = functional_theme,las = 2,outline = FALSE,col = c("#2F6B9A","#B85BCB","#2E8B57","#76C56E","#D99A4E"),border = "black",cex.axis = 0.8,cex.main = 0.9)
  dev.off()
}


## Column 3 - Group25097
plot_functional_boxplot(3,"ATPase and ATP-dependent functions")

## Column 4 - Group15114
plot_functional_boxplot(4,"Ion transport and homeostasis")

## Column 5 - Group14604
plot_functional_boxplot(5,"Ion transport and homeostasis")

## Column 6 - Group25967
plot_functional_boxplot(6,"Domains of unknown function")

## Column 7 - Group24545
plot_functional_boxplot(7,"Ion transport and homeostasis")

## Column 8 - Group3202
plot_functional_boxplot(8,"Cell envelope and surface-associated functions")

## Column 9 - Group24636
plot_functional_boxplot(9,"DNA repair and replication")

## Column 10 - Group24071
plot_functional_boxplot(10,"Cell-envelope stress response")

## Column 11 - Group23105
plot_functional_boxplot(11,"Ion transport and homeostasis")

## Column 12 - Group11386
plot_functional_boxplot(12,"Oxidative stress response and peroxide detoxification")

## Column 13 - Group15757
plot_functional_boxplot(13,"Mobile genetic elements")

## Column 14 - Group11746
plot_functional_boxplot(14,"Proteolysis and hydrolysis")

## Column 15 - Group5189
plot_functional_boxplot(15,"Stress resistance and methylation")

## Column 16 - Group5593
plot_functional_boxplot(16,"Protein interaction domains")

## Column 17 - Group22954
plot_functional_boxplot(17,"Oxidoreductase activity")

## Column 18 - Group25899
plot_functional_boxplot(18,"Stress resistance and methylation")

## Column 19 - Group22060
plot_functional_boxplot(19,"Amino acid and 2-oxocarboxylic acid metabolism")

## Column 20 - Group3330
plot_functional_boxplot(20,"Phosphate acquisition and regulation")

## Column 21 - Group25837
plot_functional_boxplot(21,"Domains of unknown function")

## Column 22 - Group3655
plot_functional_boxplot(22,"Cell envelope and surface-associated functions")

## Column 23 - Group4110
plot_functional_boxplot(23,"Signal transduction and transcriptional regulation")

## Column 24 - Group16294
plot_functional_boxplot(24,"Phosphate acquisition and regulation")

## Column 25 - Group23969
plot_functional_boxplot(25,"Proteolysis and hydrolysis")

## Column 26 - Group11975
plot_functional_boxplot(26,"Stress resistance and methylation")

## Column 27 - Group3941
plot_functional_boxplot(27,"Proteolysis and hydrolysis")

## Column 28 - Group3335
plot_functional_boxplot(28,"Phosphate acquisition and regulation")

## Column 29 - Group14851
plot_functional_boxplot(29,"Phosphate acquisition and regulation")

## Column 30 - Group3561
plot_functional_boxplot(30,"Phosphate acquisition and regulation")

## Column 31 - Group14454
plot_functional_boxplot(31,"Phosphate acquisition and regulation")

## Column 32 - Group21479
plot_functional_boxplot(32,"Phosphate acquisition and regulation")

## Column 33 - Group21520
plot_functional_boxplot(33,"Phosphate acquisition and regulation")

## Column 34 - Group27096
plot_functional_boxplot(34,"Sugar transport and exchange")

## Column 35 - Group14455
plot_functional_boxplot(35,"Phosphate acquisition and regulation")

## Column 36 - Group14453
plot_functional_boxplot(36,"Phosphate acquisition and regulation")

## Column 37 - Group11782
plot_functional_boxplot(37,"Phosphate acquisition and regulation")

## Column 38 - Group14452
plot_functional_boxplot(38,"Phosphate acquisition and regulation")

## Column 39 - Group6023
plot_functional_boxplot(39,"Sugar transport and exchange")

## Column 40 - Group3413
plot_functional_boxplot(40,"Phosphate acquisition and regulation")


save.image("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis/S11_2saliva_sHACK_prediction_Workspace.RData")
