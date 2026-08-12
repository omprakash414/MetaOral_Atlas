############ Load required libraries for all functions below

library(reshape2)
library(vegan)
library(dplyr)
library(ade4)
library(pcaPP)
library(caret)
library(randomForest)
library(pROC)
library(MLmetrics)
# library(tidyverse) (done in local)
library(ggplot2)
library(ggrepel)
library(gplots)
library(adegraphics)
# library(writexl)  (done in local)
library(compositions)
library(sfsmisc)
library(MASS)
library(metafor)
library(reshape2)
library(readxl)
library(robustbase)
library(doParallel)
library(foreach)
library(igraph)
library(tidyr)
library(psych)
library(effsize)  # only if the cohen.d() part is used in meta_lm

############ Rank Scale Function 
rank_scale=function(x){
  # x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}




############  Rank Scale Function with NA handling
rank_scale1=function(x)
{
  x <- ifelse(is.na(x),NA,rank(x))
  min_x <- min(x[!is.na(x)])
  max_x <- max(x[!is.na(x)])
  y <- ifelse(is.na(x),NA,(x-min_x)/(max_x-min_x))
  return(y);
}

############ abundance detection computation function (Prevalance computation)
compute_detection <- function(data, species_cols, group_col, group_levels) {
  detection_matrix <- matrix(0,nrow = length(species_cols),ncol = length(group_levels),dimnames = list(species_cols, group_levels))
  
  for (i in seq_along(species_cols)) {
    sp <- species_cols[i]
    
    for (j in seq_along(group_levels)) {
      grp <- group_levels[j]
      grp_values <- data[data[[group_col]] == grp, sp]
      
      detection_matrix[i, j] <- sum(grp_values != 0) / length(grp_values)
    }
  }
  
  as.data.frame(detection_matrix)
}


############ Mean Abundance computation function
compute_mean_abundance <- function(data, species_cols, group_col) {
  groups <- unique(data[[group_col]])
  mean_abundance <- matrix(NA,nrow = length(species_cols),ncol = length(groups),dimnames = list(species_cols, groups))
  
  for (g in groups) {
    subset_data <- data[data[[group_col]] == g, species_cols, drop = FALSE]
    mean_abundance[, g] <- colMeans(subset_data)
  }
  
  return(mean_abundance)
}


############ Threshold evaluation function 
evaluate_threshold_grid <- function(
  AllSpeciesDetectionPattern,
  AllSpeciesMeanAbundance,
  list_sample_groups,
  study_thresholds     = seq(0.05, 0.95, by = 0.05),
  detection_thresholds = seq(0.05, 0.95, by = 0.05)
) {

  ## initialize result matrices
  df_numb_species <- matrix(
    NA,
    nrow = length(study_thresholds),
    ncol = length(detection_thresholds),
    dimnames = list(
      study_thresholds,
      detection_thresholds
    )
  )

  df_representation_90_plus  <- df_numb_species
  df_representation_70_minus <- df_numb_species

  ## precompute detection fraction
  detection_fraction <- apply(
    AllSpeciesDetectionPattern,
    1,
    function(x) length(x[x >= detection_thresholds[1]])
  )

  ## main loop
  for (j in seq_along(study_thresholds)) {

    study_perc <- study_thresholds[j]
    message("Study threshold: ", study_perc)

    for (i in seq_along(detection_thresholds)) {

      detect <- detection_thresholds[i]

      detection_fraction <- apply(
        AllSpeciesDetectionPattern,
        1,
        function(x) length(x[x >= detect])
      ) / ncol(AllSpeciesDetectionPattern)

      temp_species <- names(which(detection_fraction >= study_perc))
      df_numb_species[j, i] <- length(temp_species)

      if (length(temp_species) > 1) {

        df_representation_90_plus[j, i] <-
          sum(colSums(AllSpeciesMeanAbundance[temp_species, , drop = FALSE]) >= 0.90) /
          length(list_sample_groups)

        df_representation_70_minus[j, i] <-
          sum(colSums(AllSpeciesMeanAbundance[temp_species, , drop = FALSE]) < 0.70) /
          length(list_sample_groups)

      } else if (length(temp_species) == 1) {

        df_representation_90_plus[j, i] <-
          sum(AllSpeciesMeanAbundance[temp_species, ] >= 0.90) /
          length(list_sample_groups)

        df_representation_70_minus[j, i] <-
          sum(AllSpeciesMeanAbundance[temp_species, ] < 0.70) /
          length(list_sample_groups)

      } else {
        df_representation_90_plus[j, i]  <- 0
        df_representation_70_minus[j, i] <- 0
      }
    }
  }

  ## Create associated identification dataframe
  df_associated_identification <- data.frame("number_of_species"=as.numeric(df_numb_species),"representation_90_plus"=as.numeric(df_representation_90_plus),"representation_70_minus"=as.numeric(df_representation_70_minus))


  return(list(
    df_numb_species               = df_numb_species,
    df_representation_90_plus     = df_representation_90_plus,
    df_representation_70_minus    = df_representation_70_minus,
    df_associated_identification  = df_associated_identification
  ))
}



############ Jaccard similarity and representation computation function
compute_detection_jaccard_representation <- function(
  AllSpeciesDetectionPattern_saliva,
  saliva_AssociatedSpecies,
  SpDf_saliva_associated,
  list_sample_groups_saliva,
  detection_thresholds = seq(0.05, 0.95, by = 0.05)
) {

  ## -----------------------------
  ## Detection counts per threshold
  ## -----------------------------
  detection_with_diff_thresholds <- matrix(
    NA,
    length(detection_thresholds),
    ncol(AllSpeciesDetectionPattern_saliva)
  )
  rownames(detection_with_diff_thresholds) <- detection_thresholds
  colnames(detection_with_diff_thresholds) <- colnames(AllSpeciesDetectionPattern_saliva)

  for (i in 1:length(detection_thresholds)) {
    thres <- detection_thresholds[i]
    detection_with_diff_thresholds[i,] <-
      apply(
        AllSpeciesDetectionPattern_saliva[saliva_AssociatedSpecies, ],
        2,
        function(x) length(x[x >= thres])
      )
  }

  ## -----------------------------
  ## Jaccard & representation matrices
  ## -----------------------------
  df_jaccard <- matrix(
    NA,
    length(detection_thresholds),
    length(list_sample_groups_saliva)
  )
  rownames(df_jaccard) <- detection_thresholds
  colnames(df_jaccard) <- list_sample_groups_saliva

  df_representation <- matrix(
    NA,
    length(detection_thresholds),
    length(list_sample_groups_saliva)
  )
  rownames(df_representation) <- detection_thresholds
  colnames(df_representation) <- list_sample_groups_saliva

  ## -----------------------------
  ## Main computation loop
  ## -----------------------------
  for (i in 1:length(detection_thresholds)) {
    threshold <- detection_thresholds[i]
    print(threshold)

    for (j in 1:length(list_sample_groups_saliva)) {
      stdy_name <- list_sample_groups_saliva[j]

      temp_core <- rownames(AllSpeciesDetectionPattern_saliva)[
        which(AllSpeciesDetectionPattern_saliva[, stdy_name] >= threshold)
      ]

      if (length(temp_core) > 1) {

        temp_sp_profile <-
          SpDf_saliva_associated[
            SpDf_saliva_associated$study_name == stdy_name,
            temp_core
          ]

        temp_sp_profile <-
          temp_sp_profile[
            rowSums(temp_sp_profile) > 0,
            colSums(temp_sp_profile) > 0
          ]

        temp_jaccard <-
          as.matrix(vegdist(temp_sp_profile, method = "jaccard", binary = TRUE))

        diag(temp_jaccard) <- NA

        df_jaccard[i, j] <-
          median(apply(temp_jaccard, 1, function(x) x[!is.na(x)]))

        df_representation[i, j] <-
          median(rowSums(temp_sp_profile))

      } else {

        ## NOTE:
        ## Original logic preserved.
        ## temp_sp_profile is undefined here in the original code
        ## and may error if this branch is reached.

        df_jaccard[i, j] <- 0
        df_representation[i, j] <- median(sum(temp_sp_profile))
      }
    }
  }

  ## -----------------------------
  ## Melted pattern dataframe
  ## -----------------------------
  df_patterns <- melt(1 - df_jaccard)
  colnames(df_patterns) <- c("Threshold", "Study", "Jaccard_Similarity")

  df_patterns$Representation <- melt(df_representation)$value
  colnames(df_patterns)[4] <- "Representation"

  df_patterns <- df_patterns[df_patterns$Representation <= 1, ]

  ## -----------------------------
  ## Return everything
  ## -----------------------------
  return(list(
    df_jaccard = df_jaccard,
    df_representation = df_representation,
    df_patterns = df_patterns
  ))
}


############ Prevalent species computation function for a single dataframe
compute_prevalent_single_data <- function(data,threshold){		
  detection_percentage <- colSums(apply(data,2,function(x)(ifelse(x>0,1,0))))/nrow(data)
  highly_detected <- names(which(detection_percentage>=threshold))
  return(highly_detected)
}




############ Keystone influence (Core Microbes) computation function
keystoneInfluence <- function(species_profile, selected_species, prevalence_threshold) {
  
  ## Create an empty list to store the dfs, one df for one study which will be generated after the analysis
  output_all_stdy_dfs <- list()
  
  ## Create an empty list to store the species which have abundance above 70%
  all_dataset_prevalent_species <- list()
  
  ## Extract the unique study_names from the species profile
  unique_study_names <- unique(species_profile$study_name)
  
  countn <- 1
  ## Add the for loop on study_names
  for (stdy in unique_study_names) {
    tryCatch({
      print(paste0("############## ",countn," ",stdy, " #################"))
      countn <- countn + 1
      
      # separate species_profile for one study
      temp_species_profile <- species_profile[which(species_profile$study_name == stdy),]
      
      # remove study_name column from the species_profile
      temp_species_profile <- temp_species_profile[, !(colnames(temp_species_profile) %in% c("study_name", "body_site_category"))]
      
      # Replace NA value with 0 if any present in species_profile of current study
      temp_species_profile[is.na(temp_species_profile)] <- 0
      
      ## Storing the species which are detected in at least in the given threshold of samples
      all_dataset_prevalent_species[[stdy]] <- compute_prevalent_single_data(temp_species_profile, prevalence_threshold)
      
      # Normalize the species_profile of current study
      temp_species_profile_norm <- temp_species_profile / rowSums(temp_species_profile)
      
      # Replace NA value with 0 if any present in species_profile of current study
      temp_species_profile_norm[is.na(temp_species_profile_norm)] <- 0
      
      # remove rows whose row sums is 0
      temp_species_profile_norm <- temp_species_profile_norm[rowSums(temp_species_profile_norm) > 0, ]
      
      # Filter only species that are present in selected_species 
      inputData <- temp_species_profile_norm[, selected_species]
      
      ## Run the ENV-Fit on each species 
      # count <- 1
      
      # Create an empty dataframe to store the p-value and r value for each species for current study
      summary_envfit_onestudy <- as.data.frame(matrix(NA, nrow = length(selected_species), ncol = 2))
      colnames(summary_envfit_onestudy) <- c("r2", "p-value")
      rownames(summary_envfit_onestudy) <- selected_species
      
      ## Loop over the species in selected_species to find the ENV-Fit results.
      for (species in selected_species) {
        tryCatch({
          # printed <- paste0("#####################", count, " ", stdy, " ", species, "###########################")
          # count <- count + 1
          # print(printed)
          
          # removing selected species (1R)
          colIndex <- which(colnames(inputData) == species)
          newData <- inputData[, -colIndex]
          
          # removing empty rows
          newData <- newData[which(rowSums(newData) != 0), ]
          #print(dim(newData))
          
          # normalizing the data (2R)
          newData <- newData / rowSums(newData)
          
          # checking if the data got renormalized
          cat(length(which(rowSums(newData) == 0)))
          
          #cat("Creating distance matrix\n")
          distanceMatrix <- vegdist(newData, method = "bray")
          
          #print("distance matrix done")
          #cat("Creating dudi.pco\n")
          
          pco <- dudi.pco(distanceMatrix, scannf = FALSE)
          
          #print("pco is created")
          pcoPointsDf <- pco$li
          
          #cat("Generating model\n") (3R)
          # since, empty rows are removed the rownames(newData) are specifically mentioned to make sure the sample sequence is same while running the model
          model <- envfit(pcoPointsDf ~ inputData[rownames(newData), species])
          
          check <- as.data.frame(t(c(model$vector$r, model$vector$pvals)))
          colnames(check) <- c("r2", "p-value")
          
          summary_envfit_onestudy[species, "r2"] <- check[1, "r2"]
          summary_envfit_onestudy[species, "p-value"] <- check[1, "p-value"]
          
          #cat(species, "done\n")
        }, error = function(e) {
          cat("Error in processing species:", species, "in study:", stdy, "\n")
          print(e)
        })
      }
      cat(" done", "\n")
      output_all_stdy_dfs[[stdy]] <- summary_envfit_onestudy
      gc()
    }, error = function(e) {
      cat("Error in processing study:", stdy, "\n")
      print(e)
    })
  }
  final_output_list <- list(ENV_fit_summary = output_all_stdy_dfs, all_dataset_prevalent_species = all_dataset_prevalent_species)
  return(final_output_list)
}







##############  Summarise r2Df and Pval dfs from CoreKeyStone
Summarise_r2pval_CoreKeystone <- function(CoreInfluencers_EnvFit_output, rank_scale_fun = rank_scale) {
  
  ## Species list (assumed identical across studies)
  species_names <- rownames(CoreInfluencers_EnvFit_output[[1]])
  
  ## Initialize r2 and p-value dataframes
  r2rd_dataframe <- data.frame(species = species_names,stringsAsFactors = FALSE)
  rownames(r2rd_dataframe) <- species_names
  
  prd_dataframe <- data.frame(species = species_names,stringsAsFactors = FALSE)
  rownames(prd_dataframe) <- species_names
  
  ## Loop over studies
  for (study in names(CoreInfluencers_EnvFit_output)) {
    
    study_df <- CoreInfluencers_EnvFit_output[[study]]
    colnames(study_df) <- c(paste0(study, "_r2"),paste0(study, "_pval"))
    
    study_df$species <- rownames(study_df)
    
    ## Merge r2 values
    r2rd_dataframe <- left_join(r2rd_dataframe, study_df[, c("species", paste0(study, "_r2"))], by = "species")
    
    ## Merge p-values
    prd_dataframe <- left_join(prd_dataframe, study_df[, c("species", paste0(study, "_pval"))], by = "species")
  }
  
  ## Clean column names
  colnames(r2rd_dataframe) <- gsub("_r2", "", colnames(r2rd_dataframe))
  colnames(prd_dataframe)  <- gsub("_pval", "", colnames(prd_dataframe))
  
  ## Restore rownames and drop species column
  rownames(r2rd_dataframe) <- r2rd_dataframe$species
  rownames(prd_dataframe)  <- prd_dataframe$species
  
  r2rd_dataframe$species <- NULL
  prd_dataframe$species  <- NULL
  
  ## Rank-scale r2 dataframe
  r2rd_dataframe_ranked <- as.data.frame(apply(r2rd_dataframe, 2, rank_scale_fun))
  
  ## Return
  return(list(
    r2df_ranked = r2rd_dataframe_ranked,
    r2df_unranked = r2rd_dataframe,
    pvalue_df  = prd_dataframe
  ))
}


############ Significant overlap and shortlist computation function
compute_EnvFit_r2_threshold <- function(prevalDf,r2_rankedDf,pval_Df,study_list,preval_cutoff,pval_cutoff = 0.05,rank_threshold = seq(0.05, 0.95, by = 0.05)) {
  
  ## Initialize accuracy matrix
  df_sig_overlap <- matrix(0,nrow = length(rank_threshold),ncol = length(study_list))
  
  rownames(df_sig_overlap) <- rank_threshold
  colnames(df_sig_overlap) <- study_list
  
  ## Loop over R² thresholds
  for (i in seq_along(rank_threshold)) {
    
    thres <- rank_threshold[i]
    
    ## Predicted core-associated taxa:
    ## prevalence ≥ cutoff AND R² ≥ threshold
    temp_core_influence <-(apply(prevalDf, 2, function(x) ifelse(x >= preval_cutoff, 1, 0))) * (apply(r2_rankedDf, 2, function(x) ifelse(x >= thres, 1, 0)))
    
    ## True significant taxa (ground truth)
    temp_sig <- apply(pval_Df, 2, function(x) ifelse(x <= pval_cutoff, 1, 0))
    
    ## Compute accuracy per study
    for (j in seq_along(study_list)) {
      
      temp_r2_vec <- as.numeric(temp_core_influence[, j])
      temp_pr_vec <- as.numeric(temp_sig[, j])
      
      df_sig_overlap[i, j] <-
        (sum(temp_pr_vec == 1 & temp_r2_vec == 1) +   # True Positives
           sum(temp_pr_vec == 0 & temp_r2_vec == 0)) /  # True Negatives
        length(temp_r2_vec)
    }
  }
  
  return(df_sig_overlap)
}


############### Summarising the Core Key Stones, like Major and Minor core. 
summarize_core_keystone_detection <- function(core_keystone_df) {

  ## Ensure input is a data frame
  core_keystone_df <- as.data.frame(core_keystone_df)

  ## Number of studies in which each species is detected
  core_keystone_df$studies_detected <- rowSums(core_keystone_df)

  ## Order species by decreasing detection
  core_keystone_df <- core_keystone_df[
    order(core_keystone_df$studies_detected, decreasing = TRUE),
  ]

  ## Store species names explicitly
  core_keystone_df$species <- rownames(core_keystone_df)

  ## Summary:
  ## how many species are detected in how many studies
  summary_df <- core_keystone_df %>%
    group_by(studies_detected) %>%
    summarise(
      n_species = n(),
      species_list = paste(species, collapse = ","),
      .groups = "drop"
    )

  ## Order by decreasing number of studies detected
  summary_df <- summary_df[order(-summary_df$studies_detected), ]

  ## Cumulative number of core/keystone species
  summary_df$cumulative_n_species <- cumsum(summary_df$n_species)

  ## Reorder back by increasing studies_detected (often nicer for plotting)
  summary_df <- summary_df[order(summary_df$studies_detected), ]

  return(list(
    core_keystone_by_species = core_keystone_df,
    detection_summary = summary_df
  ))
}



############## Multiclass OOF 
RF_CV_Multiclass <- function(df,response_col,n_folds = 5,ntree = 1000,seed = 123) {

  set.seed(seed)
  ## Prepare data
  ## ---------------------------
  df[[response_col]] <- as.factor(df[[response_col]])

  formula_rf <- as.formula(paste(response_col, "~ ."))

  ## Train control
  ## ---------------------------
  ctrl <- trainControl(
    method = "cv",
    number = n_folds,
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = multiClassSummary
  )

  ## Train Random Forest
  ## ---------------------------
  rf_model <- train(
    formula_rf,
    data = df,
    method = "rf",
    ntree = ntree,
    trControl = ctrl
  )

  ## Extract out-of-fold predictions
  ## ---------------------------
  pred_df <- rf_model$pred

  ## Ensure factor levels match
  pred_df$obs <- factor(pred_df$obs, levels = levels(df[[response_col]]))

  ## Compute one-vs-all ROC & AUC
  ## ---------------------------
  classes <- levels(pred_df$obs)

  roc_list <- list()
  auc_df <- data.frame(Class = classes, AUC = NA_real_)

  for (cl in classes) {
    actual <- ifelse(pred_df$obs == cl, 1, 0)
    prob   <- pred_df[[cl]]

    roc_obj <- roc(response = actual, predictor = prob, quiet = TRUE)

    roc_list[[cl]] <- roc_obj
    auc_df$AUC[auc_df$Class == cl] <- as.numeric(auc(roc_obj))
  }

  ## Return EVERYTHING useful
  ## ---------------------------
  return(list(
    model = rf_model,                # trained caret model
    predictions = pred_df,           # out-of-fold predictions
    roc_objects = roc_list,           # pROC ROC objects
    auc_table = auc_df,               # AUC per class
    train_control = ctrl,             # caret control object
    response_levels = classes,        # class labels
    call = match.call()               # reproducibility
  ))
}




############ Compute distances
Compute_Distances <- function(species_norm, species_raw = NULL, pseudocount = 1e-6) {

  species_norm <- as.matrix(species_norm)

  if (is.null(species_raw)) {
    species_raw <- species_norm
  } else {
    species_raw <- as.matrix(species_raw)
  }

  ## Ensure same samples and features
  common_samples <- intersect(rownames(species_norm), rownames(species_raw))
  common_features <- intersect(colnames(species_norm), colnames(species_raw))

  species_norm <- species_norm[common_samples, common_features, drop = FALSE]
  species_raw  <- species_raw[common_samples, common_features, drop = FALSE]

  ## Remove zero-sum samples
  keep <- rowSums(species_norm) > 0
  species_norm <- species_norm[keep, , drop = FALSE]
  species_raw  <- species_raw[keep, , drop = FALSE]

  distances <- list()

  ## 1. Bray–Curtis (normalized)
  ## ==================================================
  distances$bray <- vegdist(species_norm, method = "bray")

  ## 2. Kendall-based distance (normalized)
  ##    Distance = (1 - Kendall correlation) / 2
  ## ==================================================
  kendall_cor <- cor.fk(t(species_norm))
  distances$kendall <- (1 - kendall_cor) / 2

  ## 3. Aitchison distance (CLR + Euclidean)
  ##    Uses non-normalized if provided
  ## ==================================================
  species_clr <- species_raw
  species_clr[species_clr == 0] <- pseudocount

  clr_mat <- clr(species_clr)
  distances$aitchison <- dist(clr_mat, method = "euclidean")

  ## -----------------------------
  ## Return distances only
  ## -----------------------------
  return(distances)
}




############# Compute PCoA for three distance measures.
PCoA_multi_distance <- function(species_norm,metadata_df,metadata_var,species_raw = NULL,n_axes = 10,pseudocount = 1e-6) {

  ## Checks
  ## -----------------------------
  stopifnot(metadata_var %in% colnames(metadata_df))

  ## 1. Compute distances (REUSE your function)
  ## -----------------------------
  if (is.null(species_raw)) {
  distances <- Compute_Distances(species_norm = species_norm)
} else {
  distances <- Compute_Distances(species_norm = species_norm,species_raw  = species_raw,pseudocount  = pseudocount)
}
  ## Extract sample order from Bray (reference)
  sample_ids <- labels(distances$bray)

  metadata_sub <- metadata_df[sample_ids, , drop = FALSE]

  ## 2. Run PCoA for each distance
  ## -----------------------------
  pcoa_bray <- dudi.pco(distances$bray,scannf = FALSE,nf = n_axes)

  pcoa_kendall <- dudi.pco(as.dist(distances$kendall),scannf = FALSE,nf = n_axes)

  pcoa_aitchison <- dudi.pco(distances$aitchison,scannf = FALSE,nf = n_axes)

  ## 3. Extract PC scores + metadata
  ## -----------------------------
  bray_scores <- as.data.frame(pcoa_bray$li)
  bray_scores[[metadata_var]] <- metadata_sub[[metadata_var]]

  kendall_scores <- as.data.frame(pcoa_kendall$li)
  kendall_scores[[metadata_var]] <- metadata_sub[[metadata_var]]

  aitchison_scores <- as.data.frame(pcoa_aitchison$li)
  aitchison_scores[[metadata_var]] <- metadata_sub[[metadata_var]]

  ## -----------------------------
  return(list(
    distances = list(bray = distances$bray,kendall = distances$kendall,aitchison = distances$aitchison),
    pcoa_scores = list(bray = bray_scores,kendall = kendall_scores,aitchison = aitchison_scores)
  ))
}






########################## Disease Analysis Specific functions: (From HACK Paper)

#############   Make the summary of data available for disease analysis
make_summary_single_subsite <- function(metadata_df,
                                        min_control = 10,
                                        min_disease = 10) {
  
  # Summarise counts per study
  summary_data <- metadata_df %>%
    group_by(study_name) %>%
    summarise(
      control_count = sum(study_condition == "Control"),
      disease_count = sum(study_condition == "Diseased"),
      diseases      = paste(unique(disease), collapse = ", "),
      .groups = "drop"
    )
  
  # Apply thresholds and 80–20 balance rule
  summary_data <- summary_data %>%
    filter(control_count >= min_control,
                  disease_count >= min_disease) %>%
    mutate(
      total = control_count + disease_count,
      prop_control = control_count / total,
      `80_20` = ifelse(prop_control < 0.2 | prop_control > 0.8, 0, 1)
    )
  
  # Filtered version (balanced studies only)
  summary_data_filtered <- summary_data %>%
    filter(`80_20` == 1)
  
  return(list(
    summary_raw = summary_data,
    summary_filtered = summary_data_filtered
  ))
}


library(dplyr)

#############   Make the sample list for one subsite (need output from make_summary_single_subsite)
make_sample_lists_single_subsite <- function(metadata_df,
                                             summary_info_subsite) {
  
  # Keep only studies that passed 80–20 rule
  valid_studies <- unique(summary_info_subsite$study_name)
  
  # Filter metadata for this subsite + valid studies
  subsite_metadata <- metadata_df %>%
    filter(
      study_name %in% valid_studies
    )
  
  # -------------------------------
  # CONTROL samples: study → samples
  # -------------------------------
  controlSampleList <- subsite_metadata %>%
    filter(study_condition == "Control") %>%
    group_by(study_name) %>%
    summarise(samples = list(sample_id), .groups = "drop")
  
  controlSampleList <- setNames(
    controlSampleList$samples,
    controlSampleList$study_name
  )
  
  # -------------------------------------
  # DISEASE samples: disease → study → samples
  # -------------------------------------
  diseaseSampleList <- list()
  
  disease_meta <- subsite_metadata %>%
    filter(study_condition == "Diseased")
  
  for (dis in unique(disease_meta$disease)) {
    diseaseSampleList[[dis]] <- list()
    
    for (study in unique(disease_meta$study_name[disease_meta$disease == dis])) {
      
      samples <- disease_meta %>%
        filter(disease == dis, study_name == study) %>%
        pull(sample_id)
      
      diseaseSampleList[[dis]][[study]] <- samples
    }
  }
  
  return(list(
    controlSampleList = controlSampleList,
    diseaseSampleList = diseaseSampleList
  ))
}





############# batch Wilcox Test function 
mannWhitney_batch = function(x,y){
  print("in man funciton")
  p_array <- NULL;
  medianDirection <- NULL;
  meanDirection <- NULL;
  
  z <- intersect(rownames(x),rownames(y));
  # print(rownames(x))
  for(i in 1:length(z))
  {
    printed= paste0(i," out of ", length(z))
    # print(printed)
    printed= paste0(i," ",z[i])
    # print(printed)
    p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
    
    medianDirection[i] <- ifelse(median(as.numeric(x[z[i],]),na.rm=TRUE) > median(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(median(as.numeric(x[z[i],]),na.rm=TRUE) < median(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    # print(medianDirection[i])
    meanDirection[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    # print(meanDirection[i])
    i <- i + 1;
  }
  out <- as.data.frame(cbind(medianDirection, meanDirection, p_array,p.adjust(p_array, method="fdr")));
  colnames(out)[4]= "q_array"
  
  
  rownames(out) <- z;
  # NA values will be either in p-value column or q-value column, if we take them as 1 it will still be insignificant and won't be considered.
  out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
  # print(head(out))
  return(t(out));
}


############# 

targetIntermediateOutput= function(controlSampleList, diseaseSampleList, AllCombinedSpProfile, uniqueDiseases, species_to_check, saveLocationString, fileName, normalizationRequired)
{
  # making sure that all the diseases from uniqueDiseases are present in diseaseSampleList
  if(length(setdiff(uniqueDiseases,names(diseaseSampleList)))>0)
  {
    print("The sample ids are not available for all the disesases you took please make sure diseaseSampleList contains all the uniqueDisease.")
    return()
  }
  
  # making sure that for all the studies taken for the disease profile their control abundance profile counterpart is available.
  
  allStudies= c()
  for(dis in names(diseaseSampleList))
  {
    allStudies= union(allStudies, names(diseaseSampleList[[dis]]))
  }
  
  if(length(setdiff(allStudies,names(controlSampleList)))>0)
  {
    print("all the studies that are taken for disease abundance profile, their control abundance profile counter part is not available")
    return()
  }
  
  
  # the list will contain the output for all the selected diseases.
  intermediateOutput= list() 
  
  cnt=1
  for(dis in uniqueDiseases)
  {
    printed= paste0(cnt," out of ", length(uniqueDiseases))
    print(printed)
    cnt= cnt+1
    
    # for the given disease taking all the disease and control samples from all the studies in which the disease was present.
    diseaseSpecificStudies= names(diseaseSampleList[[dis]])
    
    diseaseSpecificControlSamples= c()
    diseaseSpecificDiseaseSamples= c()
    
    for(study in diseaseSpecificStudies)
    {
      diseaseSpecificControlSamples= c(diseaseSpecificControlSamples, controlSampleList[[study]])
    }
    
    for(study in diseaseSpecificStudies)
    {
      diseaseSpecificDiseaseSamples= c(diseaseSpecificDiseaseSamples, diseaseSampleList[[dis]][[study]])
    }
    
    # removing empty rows from the control and disease species profile after taking 196 species.
    # for our disease analysis there are no empty rows, all the samples were considered.
    
    # global_allControlSamples<<- union(global_allControlSamples, diseaseSpecificControlSamples)
    # global_allDiseaseSamples<<- union(global_allDiseaseSamples, diseaseSpecificDiseaseSamples)
    
    controlSpProfile= data.frame()
    diseaseSpProfile= data.frame()
    
    print("working on specific species vector")
    
    # checking if there any species that are absent in the AllCombinedSpProfile, if yes then they have to be added with 0 values
    absentSpecies= setdiff(species_to_check,colnames(AllCombinedSpProfile))
    
    if(length(absentSpecies)>0)
    {
      print("some of the species_to_check species are absent in the species profile matrix")
      # print(absentSpecies)
      for(species in absentSpecies)
      {
        AllCombinedSpProfile[[species]]= with(AllCombinedSpProfile,0)
        printed= paste0("absent species is ", species)
        print(printed)
        
        # making sure that the species got added with the 0 values.
        print(unique(AllCombinedSpProfile[,species]))
      }
      
    }
    
    controlSpProfile= AllCombinedSpProfile[diseaseSpecificControlSamples, species_to_check]
    controlSpProfile= controlSpProfile[rowSums(controlSpProfile)>0,]
    # print(head(controlSpProfile))
    # normalizing the sp profile
    if(normalizationRequired==TRUE)
    {
      print("normalizing the control data")
      controlSpProfile= controlSpProfile/rowSums(controlSpProfile)  
      print(unique(rowSums(controlSpProfile)))
    }
    
    # transposing the data to make it compatible with the mannWhitney_batch funciton.
    controlSpProfile= as.data.frame(t(controlSpProfile))
    
    diseaseSpProfile= AllCombinedSpProfile[diseaseSpecificDiseaseSamples, species_to_check]
    diseaseSpProfile= diseaseSpProfile[rowSums(diseaseSpProfile)>0,]
    # normalizing the sp profile
    if(normalizationRequired==TRUE)
    {
      print("normalizing the disease data")
      diseaseSpProfile= diseaseSpProfile/rowSums(diseaseSpProfile)  
      print(unique(rowSums(diseaseSpProfile)))
    }
    
    # transposing the data to make it compatible with the mannWhitney_batch funciton.
    diseaseSpProfile= as.data.frame(t(diseaseSpProfile))
    
    printed= paste0("disease df dimensions =  ", dim(diseaseSpProfile))
    print(printed)
    
    printed= paste0("control df dimensions =  ", dim(controlSpProfile))
    print(printed)
    intermediateOutput[[dis]]= mannWhitney_batch(diseaseSpProfile,controlSpProfile)
    
  }
  
  save(intermediateOutput, file= paste0(saveLocationString,"\\",fileName))
  return(intermediateOutput)
}


###############

associationOutput= function(intermediateOutput, uniqueDiseases, outputFolder, rowDendogramed, colDendogramed)
{
  # creating the empty dataframe that contains information about the selected species(speceies_to_check) and the total number of diseases. It will contain association
  # value for all the species for each disease.
  meanBasedOutput= as.data.frame(matrix(nrow= nrow(intermediateOutput[[1]]), ncol= length(intermediateOutput)))
  rownames(meanBasedOutput)= rownames(intermediateOutput[[1]])
  colnames(meanBasedOutput)= names(intermediateOutput)
  
  cnt=1
  for(output in colnames(meanBasedOutput))
  {
    
    # getting the categorical directions as per q_value, p_value and meanDirection
    vector= c()
    for(rows in 1:nrow(meanBasedOutput))
    {
      if(intermediateOutput[[output]][rows,"meanDirection"]==0)
      {
        vector= c(vector,0)
      }
      else if(intermediateOutput[[output]][rows,"q_array"]<= 0.15)
      {
        vector= c(vector, intermediateOutput[[output]][rows,"meanDirection"]*3)
      }
      else if(intermediateOutput[[output]][rows,"p_array"]<= 0.05)
      {
        vector= c(vector, intermediateOutput[[output]][rows,"meanDirection"]*2)
      }
      else
      {
        vector= c(vector, intermediateOutput[[output]][rows,"meanDirection"]*1)
      }
      
    }
    meanBasedOutput[,output]= vector
    
  }
  
  print("getting the heatmap carpet")
  pdf(file = paste0(outputFolder, "\\heatmap.pdf"), width = 20, height = 15)
  heatmap2 <- heatmap.2(as.matrix(t(meanBasedOutput)) , density= "none", trace= "none", Rowv=rowDendogramed, Colv= colDendogramed)
  # print(names(heatmap2))
  # env$heatmap= heatmap2
  
  heatmapValues= heatmap2$carpet
  save(heatmapValues, file= paste0(outputFolder,"\\heatmapCarpet.RData"))
  write.table(rownames(heatmapValues), paste0(outputFolder,"\\AssociationHeatmap_colnames.txt"))
  write.table(colnames(heatmapValues), paste0(outputFolder,"\\AssociationHeatmap_rownames.txt"))
  
  dev.off()
  
  # getting the association count
  associationCount= data.frame("a")
  # checkMeanOutput <<- meanBasedOutput
  # it will get associations data for all 196 species instead of filtered one 
  # tempMeanOuptut<<- meanBasedOutput
  for(species in rownames(meanBasedOutput))
  {
    # print(species)
    currDf= as.data.frame(t(data.frame(table(unlist(meanBasedOutput[species,])))))
    columnNames= currDf[1,]
    currDf= as.data.frame(currDf[-1,])
    colnames(currDf)= columnNames
    currDf$species= species
    
    associationCount= merge(associationCount,currDf, all.x= TRUE, all.y= TRUE)
  }
  
  index = which(colnames(associationCount)=="X.a.")
  print(colnames(associationCount))
  associationCount= associationCount[,-index]
  
  species= associationCount$species
  
  index = which(colnames(associationCount)=="species")
  associationCount= associationCount[,-index]
  
  associationCount= as.data.frame(apply(associationCount,2,as.numeric))
  rownames(associationCount)= species
  
  associationCount[is.na(associationCount)]=0
  print("checkig if all the columns are present or not")
  for(columns in c("-3","-2","2","3","-1","0","1"))
  {
    if(!columns %in% colnames(associationCount))
    {
      associationCount[[columns]]= with(associationCount,0)  
      print(columns)
    }
    
    # making sure that the species got added with the 0 values.
    print(unique(associationCount[,columns]))
  }
  # tempAssociation <<- associationCount
  
  associationCount$totalSignificant= rowSums(associationCount[,c("-3","3")])
  associationCount$totalNegativeSignificant =
    rowSums(associationCount[, "-3", drop = FALSE])
  
  associationCount$totalPositiveSignificant =
    rowSums(associationCount[, "3", drop = FALSE])
  
  associationCount$totalNonSignificant= rowSums(associationCount[,c("-1","0","1")])
  associationCount= associationCount[,c("-3","-2","-1","0","1","2","3","totalSignificant","totalNegativeSignificant",
                                        "totalPositiveSignificant","totalNonSignificant")]
  
  associationCount$total= length(uniqueDiseases)
  
  associationCount$scaledDifference= ((associationCount$totalNegativeSignificant - associationCount$totalPositiveSignificant)/associationCount$total) * (1-(min(associationCount$totalNegativeSignificant,associationCount$totalPositiveSignificant) + 0.0001 /max(associationCount$totalNegativeSignificant,associationCount$totalPositiveSignificant)+ 0.0001))
  
  rank_scale=function(x)
  {
    # x <- rank(x);
    y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
    y <- ifelse(is.nan(y),0,y)
    return(y);
  }
  
  
  associationCount$ranked_scaledDifference= rank_scale(associationCount$scaledDifference)
  
  write_xlsx(associationCount, paste0(outputFolder,"\\associationCount.xlsx"))
  
  # for heatmap2
  heatmap2Df= associationCount
  heatmap2Df= heatmap2Df[rownames(heatmapValues),]
  
  heatmap2Df$species= rownames(heatmap2Df)
  heatmap2Df$species <- factor(heatmap2Df$species, levels = unique(heatmap2Df$species))
  
  # getting the names of only those species that will be labeled, if significant or non significant 
  # contribute to >=70% then it will be labeled
  heatmap2Df$marker= NA
  #View(heatmap2Df)
  counter=1
  
  for(species in rownames(heatmap2Df))
  {
    printed= paste0(counter," ", species)
    print(printed)
    counter= counter+1
    
    if(heatmap2Df[species,"totalSignificant"]==0)
    {
      heatmap2Df[species,"marker"]= "no"
    }
    
    else if(heatmap2Df[species,"totalSignificant"]/heatmap2Df[species,"total"]<0.30)
    {
      heatmap2Df[species,"marker"]= "no"
    }
    
    else if(heatmap2Df[species,"totalNegativeSignificant"]/heatmap2Df[species,"totalSignificant"]>=0.70)
    {
      heatmap2Df[species,"marker"]= "yes"
    }
    else if(heatmap2Df[species,"totalPositiveSignificant"]/heatmap2Df[species,"totalSignificant"]>=0.70)
    {
      heatmap2Df[species,"marker"]= "yes"
    }
    else
    {
      heatmap2Df[species,"marker"]= "no"
    }
    
  }
  
  heatmap2Df$diseaseAssociation= NA
  index= which(heatmap2Df$marker=="yes" & heatmap2Df$totalNegativeSignificant> heatmap2Df$totalPositiveSignificant)
  
  heatmap2Df[index, "diseaseAssociation"]= "Negative Association With Disease"
  
  index= which(heatmap2Df$marker=="yes" & heatmap2Df$totalNegativeSignificant < heatmap2Df$totalPositiveSignificant)
  heatmap2Df[index, "diseaseAssociation"]= "Positive Association With Disease"
  
  # saving the selected major taxa that are have yes in marker column and have more negatively
  # associated studies
  
  
  write_xlsx(heatmap2Df,paste0(outputFolder,"\\heatmap2Df.xlsx"))
  
  # melting the data
  df_long= heatmap2Df[,c("totalNegativeSignificant", "totalPositiveSignificant","species")] %>%
    pivot_longer(!species, names_to = "variable", values_to = "value")
  
  df_long2= heatmap2Df[,c("totalNegativeSignificant", "totalPositiveSignificant","marker")] %>%
    pivot_longer(!marker, names_to = "variable", values_to = "value")
  
  df_long$marker= df_long2$marker
  pdf(file = paste0(outputFolder, "/linePlot.pdf"), width = 20, height = 10)
  plot= ggplot(df_long,
               aes(x= species,
                   y= value, group=variable))+
    geom_line(aes(color= variable),alpha = 0.5,size=0.5)+
    geom_point() +
    scale_color_manual(values=c("green", "red"))+
    xlab("Species")+
    ylab("sutdy count")+
    theme(axis.text = element_text(color= "black"))+
    theme_bw()+
    theme(axis.text.x= element_text(angle= 90, vjust = 0.5, size=6))
  
  dev.off()
  return() 
}

#################

diseaseAnalysisIterationPipeline= function(iteration, allDiseases, controlSampleList, diseaseSampleList, MainOutputFolder, AllCombinedSpProfile, species_to_check)
{
  # it will store the information about what diseases are used in each iteration.
  iterationSpecificDiseaseInfo= list()
  
  for(repeation in 1:iteration)
  {
    printed= paste0("######################################################################################  ", repeation,
                    "  ######################################################################################")
    print(printed)
    
    # randomly picking 65% of diseases out of all 26 disesaes. 17 in this case.
    totalDiseases = length(allDiseases)
    totalSelected_diseases = round(totalDiseases * 0.65, 0)
    selectedIndex = sample(seq_len(totalDiseases),
                       totalSelected_diseases,
                       replace = FALSE)
					   
    uniqueDiseases= allDiseases[selectedIndex]
    iterationSpecificDiseaseInfo[[repeation]]= uniqueDiseases
    
    # getting all the disease samples for the selected diseases.
    AllDiseaseSamples= c()
    iterationSpecificStudies= c()
    
    for(dis in uniqueDiseases)
    {
      for(study in names(diseaseSampleList[[dis]]))
      {
        AllDiseaseSamples= union(AllDiseaseSamples, diseaseSampleList[[dis]][[study]]) 
      }
      iterationSpecificStudies= union(iterationSpecificStudies, names(diseaseSampleList[[dis]]))
    }
    
    # getting all the control samples for the selected diseases.
    AllControlSamples= c()
    for(study in iterationSpecificStudies)
    {
      AllControlSamples= c(AllControlSamples, controlSampleList[[study]])
    }
    
    # taking subset for control and disease samples as per the current iteration.
    subControlSampleList= controlSampleList[iterationSpecificStudies]
    subDiseaseSampleList= diseaseSampleList[uniqueDiseases]
    
    # removing the empty rows after filtering 196 species
    tempSpProfile= AllCombinedSpProfile[AllControlSamples, species_to_check]
    tempSpProfile= tempSpProfile[rowSums(tempSpProfile)>0,]
    AllControlSamples= rownames(tempSpProfile)
    
    
    tempSpProfile= AllCombinedSpProfile[AllDiseaseSamples,species_to_check]
    tempSpProfile= tempSpProfile[rowSums(tempSpProfile)>0,]
    AllDiseaseSamples= rownames(tempSpProfile)
    
    
    
    dir.create(paste0(MainOutputFolder,"\\",repeation))
    outputFolder= paste0(MainOutputFolder,"\\",repeation)
    
    # there will be difference in the count of the samples in the list and the sample vectors as we removed some of the samples from control and disease after filtering 196 species.
    save(subControlSampleList,
         subDiseaseSampleList,
         AllControlSamples,
         AllDiseaseSamples,
         file= paste0(outputFolder,"\\SampleInfo.RData"))
    
    
    intermediateOutput= targetIntermediateOutput(subControlSampleList, subDiseaseSampleList, AllCombinedSpProfile, uniqueDiseases, species_to_check, outputFolder, "mannWhitneyIntermediateOutput.RData", TRUE)
    
    associationOutput(intermediateOutput, uniqueDiseases, outputFolder, TRUE, TRUE)
  }
  
  # stroing the results of all iterations in the list.
  allIterationFolders= list.files(MainOutputFolder)
  allIterationFolders <- allIterationFolders[grepl("^[0-9]+$", allIterationFolders)]
  
  iterationOutputDf_direction= as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_direction)= species_to_check
  colnames(iterationOutputDf_direction)= allIterationFolders
  
  iterationOutputDf_divisionScore= as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_divisionScore)= species_to_check
  colnames(iterationOutputDf_divisionScore)= allIterationFolders
  
  iterationOutputDf_subtractionScore= as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_subtractionScore)= species_to_check
  colnames(iterationOutputDf_subtractionScore)= allIterationFolders
  
  iterationOutputDf_ranked_scaledDifference = as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_ranked_scaledDifference)= species_to_check
  colnames(iterationOutputDf_ranked_scaledDifference)= allIterationFolders
  
  iterationOutputDf_scaledDifference = as.data.frame(matrix(nrow= length(species_to_check), ncol= iteration))
  rownames(iterationOutputDf_scaledDifference)= species_to_check
  colnames(iterationOutputDf_scaledDifference)= allIterationFolders
  
  
  for(folder in allIterationFolders)
  {
    heatmapInfo= read_excel(paste0(MainOutputFolder,"\\",folder,"\\heatmap2Df.xlsx"))
    heatmapInfo= as.data.frame(heatmapInfo)
    rownames(heatmapInfo)= heatmapInfo$species
    
    iterationOutputDf_ranked_scaledDifference[species_to_check,folder]= heatmapInfo[species_to_check,"ranked_scaledDifference"]
    iterationOutputDf_scaledDifference[species_to_check,folder]= heatmapInfo[species_to_check,"scaledDifference"]
  }
  
  
  outputList= list()
  outputList[["iterationSpecificDiseaseNames"]]= iterationSpecificDiseaseInfo
  outputList[["iterationOutputDf_ranked_scaledDifference"]]= iterationOutputDf_ranked_scaledDifference
  outputList[["iterationOutputDf_scaledDifference"]]= iterationOutputDf_scaledDifference
  
  return(outputList) 
  dev.off()
}

################ Function to calculate the disease association score (modified version of the score used in GutHACK paper)
healthAssociation_iterations <- function(
    n_iterations = 10, study_fraction = 0.65,
    matched_studies,discovery_controls,discovery_diseased,metadata_df,species_profile,species_vector){
  
  all_iteration_scores <- matrix(NA,nrow = length(species_vector),ncol = n_iterations)
  rownames(all_iteration_scores) <- species_vector
  colnames(all_iteration_scores) <- paste0("Iter_",1:n_iterations)
  
  selected_studies_per_iteration <- list()
  df_directions_all <- list()  
  df_comparison_last <- list()
  for(iter in 1:n_iterations)
  {
    cat("\nIteration:",iter,"\n")
    ##############################
    ## randomly select 65% studies
    ##############################
    selected_studies <- sample(matched_studies,round(length(matched_studies) * study_fraction),replace = FALSE)
    selected_studies_per_iteration[[iter]] <- selected_studies
    
    ##############################
    ## comparison matrix
    ##############################
    
    df_comparison <- as.data.frame(matrix(NA,nrow = length(species_vector),ncol = length(selected_studies)))
    rownames(df_comparison) <- species_vector
    colnames(df_comparison) <- selected_studies
    
    ##############################
    ## study loop
    ##############################
    
    for(i in seq_along(selected_studies))
    {
      study <- selected_studies[i]
      print(paste0(i,". ",study))
      
      all_study_samples <-rownames(metadata_df[metadata_df$study_name == study,])
      
      study_controls <- intersect(discovery_controls,all_study_samples)
      
      study_diseased <-intersect(discovery_diseased,all_study_samples)
      
      if(length(study_controls) < 2 ||length(study_diseased) < 2){
        next
      }
      
      ##############################
      ## species loop
      ##############################
      
      for(species_name in species_vector)
      {
        control_values <- as.numeric(species_profile[study_controls,species_name])
        disease_values <-as.numeric(species_profile[study_diseased,species_name])
        
        temp_cohen_d <- cohen.d(control_values,disease_values)$estimate
        temp_wilcox <- wilcox.test(control_values,disease_values)$p.value
        df_comparison[species_name,study] <- sign(temp_cohen_d) * ifelse(temp_wilcox <= 0.05,3,ifelse(temp_wilcox <= 0.10,2,1))
      }
      
    }
    ##############################
    ## score calculation
    ##############################
    
    df_directions <- data.frame(
      control_increased =apply(df_comparison,1,function(x)length(x[!is.na(x) & x > 2])),
      control_depleted = apply(df_comparison,1,function(x)length(x[!is.na(x) & x < -2])))
    
    n_studies_used <- length(selected_studies)
    df_directions$score <- apply(df_directions,1,function(x) {((x[1] - x[2]) / n_studies_used) * (1 - ((min(x[1],x[2]) + 0.00001) / (max(x[1],x[2]) + 0.00001)))})
    all_iteration_scores[rownames(df_directions),iter] <- df_directions$score
    
    df_directions_all[[iter]] <- df_directions
    df_comparison_last[[iter]] <- df_comparison
  }
  
  all_iteration_scores <- as.data.frame(all_iteration_scores)
  all_iteration_scores_ranked <-as.data.frame(apply(all_iteration_scores,2,rank_scale))
  rownames(all_iteration_scores_ranked) <- rownames(all_iteration_scores)
  
  mean_scores <- data.frame(HealthAssociationScore = rowMeans(all_iteration_scores_ranked,na.rm = TRUE))
  
  return(
    list(iteration_scores = all_iteration_scores, RankedIteration_scores = all_iteration_scores_ranked, HealthAssociationScore = mean_scores,selected_studies = selected_studies_per_iteration,df_directions_all = df_directions_all, df_comparison_last = df_comparison_last))
  
}


################ Function to calculate the time-point based distance i.e followup distances.
Longitudinal_Microbiome_Distance <- function(metadata, rel_abund_matrix) {
  metadata$Jaccard_dist    <- NA
  metadata$Aitchison_dist  <- NA
  metadata$Kendall_dist    <- NA
  metadata$BrayCurtis_dist <- NA
  
  for (i in 1:nrow(metadata)) {
    samp <- metadata$Sample_ID[i]
    fu   <- metadata$Follow_up[i]
    
    if (is.na(fu) || fu %in% c("-", "_")) {
      next  # last timepoint → NA distances
    }
    
    df <- rel_abund_matrix[c(samp, fu), , drop = FALSE]
    
    # Jaccard (presence/absence on rel abund)
    jac_dist <- as.matrix(vegdist(df, method = "jaccard", binary = TRUE))
    jac <- round(jac_dist[1, 2], 3)
    
    # Euclidean ("Aitchison" in gut pipeline, but no clr)
    ait_dist <- as.matrix(vegdist(df, method = "euclidean"))
    ait <- round(ait_dist[1, 2], 3)
    
    # Kendall distance (1 - Kendall tau)
    kend_dist <- as.matrix(as.dist(1 - cor(t(df), method = "kendall")))
    kend <- round(kend_dist[1, 2], 3)
    
    # Bray–Curtis
    bray_dist <- as.matrix(vegdist(df, method = "bray"))
    bray <- round(bray_dist[1, 2], 3)
    
    metadata$Jaccard_dist[i]    <- jac
    metadata$Aitchison_dist[i]  <- ait
    metadata$Kendall_dist[i]    <- kend
    metadata$BrayCurtis_dist[i] <- bray
  }
  
  rownames(metadata) <- metadata$Sample_ID
  return(metadata)
}



################  Meta-analysis using robust linear models for multiple groups
compute_meta_lm_group <- function(data,feature_list,metadata_var,grouping_var,grouping_list)
{
  return_out <- as.data.frame(matrix(NA,length(feature_list),10))
  rownames(return_out) <- feature_list
  colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","consistency","qval","dir")
  for(i in 1:length(feature_list))
  {
    species_name <- feature_list[i]
    tryCatch(               
      expr = {      
        print(species_name)
        temp_res <- compute_meta_lm(data,species_name,metadata_var,grouping_var,grouping_list)
        return_out[i,"beta"] <- temp_res$model$beta
        return_out[i,"pval"] <- temp_res$model$pval
        return_out[i,"ci.ub"] <- temp_res$model$ci.ub
        return_out[i,"ci.lb"] <- temp_res$model$ci.lb
        return_out[i,"tau2"] <- temp_res$model$tau2
        return_out[i,"QE"] <- temp_res$model$QE
        return_out[i,"QEp"] <- temp_res$model$QEp
        return_out[i,"consistency"] <- length(which(sign(temp_res$df_studies$di) == sign(as.numeric(temp_res$model$beta))))/nrow(temp_res$df_studies)		
      },
      error = function(e){ 
        return_out[i,"beta"] <- 0
        return_out[i,"pval"] <- 1
        return_out[i,"ci.ub"] <- 0
        return_out[i,"ci.lb"] <- 0
        return_out[i,"tau2"] <- 1
        return_out[i,"QE"] <- 1
        return_out[i,"QEp"] <- 1
        return_out[i,"consistency"] <- 0
        print(e)
        print("Error observed. Moving to next")
      },
      finally = {            
        print("finally Executed")
      }
    )
    
  }
  return_out$qval <- p.adjust(return_out$pval,method="fdr")
  #return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.05,2*sign(return_out$beta),sign(return_out$beta)))
  return_out$dir <- ifelse(return_out$pval <= 0.05,2*sign(return_out$beta),sign(return_out$beta))  ## Updated 20.03.2026
  
  return(return_out)
}





###############  Meta-analysis using robust linear models for one group
compute_meta_lm <- function(data,var1,var2,grouping_variable,grouping_list)
{
  temp_meta <- data.frame(matrix(0,length(grouping_list),6))
  colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
  for(i in 1:length(grouping_list))
  {
    group <- grouping_list[i]
    temp_meta[i,1] <- group
    vec1 <- data[data[,grouping_variable]==group,var1]
    vec2 <- data[data[,grouping_variable]==group,var2]
    
    #print(paste0(group,",",length(vec1[abs(vec1)>0]),",",length(vec2[abs(vec2)>0])))
    if((length(vec1[abs(vec1)>0]) > 0)&&(length(vec2[abs(vec2)>0]) > 0))
    {
      #print(data[data[,grouping_variable]==group,c(var1,var2)])
      print(group)
      tryCatch(               
        expr = {    
          f <- as.formula(paste0(var1,"~",var2))
          temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
          summary_temp_rlm <- summary(temp_rlm)
          #print("Enter")
          print(summary_temp_rlm)
          temp_meta[i,2] <- summary_temp_rlm$coefficients[2,3]
          temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
          temp_meta[i,4] <- 1
          temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
          temp_meta[i,6] <- sign(temp_meta[i,2])
          
        },
        error = function(e){ 
          temp_meta[i,2] <- 0
          temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
          temp_meta[i,4] <- 1
          temp_meta[i,5] <- 1
          temp_meta[i,6] <- 1
          print(e)
          print("Error observed. Moving to next")
        },
        finally = {            
          print("finally Executed")
        }
      )
      
      
    }
    else
    {
      #temp_meta[i,2] <- 0
      #temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
      #temp_meta[i,4] <- 1
      #temp_meta[i,5] <- 1
      #temp_meta[i,6] <- 1
    }
  }
  temp_meta <- temp_meta[!is.na(temp_meta[,"ti"])&(temp_meta[,"ti"] != 0),]
  grouping_list <- temp_meta$dataset
  print(temp_meta)
  temp_meta <- mutate(temp_meta,study_id=grouping_list)
  rownames(temp_meta) <- grouping_list
  
  #temp_meta <- temp_meta %>% select(study_id, ri:ni)
  temp_meta <- escalc(measure="ZPCOR",mi=mi,ni=ni,ti=ti,data=temp_meta)
  res <- rma(yi, vi, data=temp_meta)
  res$ids <- rownames(temp_meta)
  res$slabs <- rownames(temp_meta)
  return_list <- list("df_studies"=temp_meta,"model"=res)
  return(return_list)
}




###############  Meta-analysis using robust linear models for multiple groups
Meta_lm_Iterrative <- function(data,species,dist_var,study_var,n_studies,n_iter = 10,seed = 1) {
  set.seed(seed)
  studies <- unique(data[[study_var]])

  iter_results <- vector("list", n_iter)

  for (i in seq_len(n_iter)) {
    selected_studies <- sample(studies, n_studies, replace = FALSE)

    iter_results[[i]] <- compute_meta_lm_group(
      data = data,
      feature_list = species,
      metadata_var = dist_var,
      grouping_var = study_var,
      grouping_list = selected_studies
    )
  }

  names(iter_results) <- paste0("iter", seq_len(n_iter))
  return(iter_results)
}



################  Compute Meta Stability Score
Compute_Meta_Stability <- function(
  iter_list,
  species,
  use_qval = TRUE
) {
  # Collect beta & qval matrices
  beta_mat <- sapply(iter_list, function(x) x$beta)
  qval_mat <- sapply(iter_list, function(x) if (use_qval) x$qval else x$pval)
  cons_mat <- sapply(iter_list, function(x) x$consistency)

  rownames(beta_mat) <- species
  rownames(qval_mat) <- species
  rownames(cons_mat) <- species

  # Direction preserved (same sign flips as your code)
  beta_mat <- -1 * beta_mat
  qval_mat <- -log10(qval_mat)

  # Stability score (same formula)
  stability_mat <- mapply(
    function(b, q, c) rank_scale1(b * q * c),
    as.data.frame(beta_mat),
    as.data.frame(qval_mat),
    as.data.frame(cons_mat),
    SIMPLIFY = FALSE
  )
  stability_mat <- as.data.frame(stability_mat)
  rownames(stability_mat) <- species

  # Summary statistics
  list(
    beta_summary = data.frame(
      median_beta = apply(beta_mat, 1, median, na.rm = TRUE),
      iqr_beta    = apply(beta_mat, 1, IQR, na.rm = TRUE)
    ),
    qval_summary = data.frame(
      median_qval = apply(qval_mat, 1, median, na.rm = TRUE),
      iqr_qval    = apply(qval_mat, 1, IQR, na.rm = TRUE)
    ),
    stability = data.frame(
      mean = apply(stability_mat, 1, mean, na.rm = TRUE),
      IQR  = apply(stability_mat, 1, IQR, na.rm = TRUE)
    ),
    stability_mat = stability_mat
  )
}




########### Meta-analysis using robust linear models for two variables
compute_meta_lm2 <- function(data, var1, var2, grouping_variable, grouping_list) {
  subset_data <- data[data[[grouping_variable]] %in% grouping_list, ]
  split_data <- split(subset_data, subset_data[[grouping_variable]])
  
  n <- length(grouping_list)
  temp_meta <- data.frame(
    dataset = grouping_list,
    ti = numeric(n),
    ni = integer(n),
    mi = rep(1, n),
    pi = numeric(n),
    di = numeric(n),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(grouping_list)) {
    group <- grouping_list[i]
    group_data <- split_data[[group]]
    vec1 <- group_data[[var1]]
    vec2 <- group_data[[var2]]
    
    if (sum(abs(vec1) > 0) > 0 && sum(abs(vec2) > 0) > 0) {
      f <- as.formula(paste0(var1, " ~ ", var2))
      temp_rlm <- rlm(f, data = group_data)
      s <- summary(temp_rlm)
      temp_meta$ti[i] <- s$coefficients[2, 3]
      temp_meta$ni[i] <- nrow(group_data)
      temp_meta$pi[i] <- f.robftest(temp_rlm, var = var2)$p.value
      temp_meta$di[i] <- sign(temp_meta$ti[i])
    } else {
      temp_meta$ti[i] <- 0
      temp_meta$ni[i] <- nrow(group_data)
      temp_meta$pi[i] <- 1
      temp_meta$di[i] <- 1
    }
  }
  
  rownames(temp_meta) <- grouping_list
  temp_meta <- escalc(measure = "ZPCOR", mi = mi, ni = ni, ti = ti, data = temp_meta)
  
  res <- rma(yi, vi, data = temp_meta)  # Use faster method
  res$ids <- rownames(temp_meta)
  res$slabs <- rownames(temp_meta)
  
  list(df_studies = temp_meta, model = res)
}


############# Meta-analysis using robust linear models for two variables over multiple species
Rem_Network2 <- function(data, species_list, group_name, study_list, feature_list) {
  #data <- subset(data, select = -c(disease, body_site_category, sample_id))
  # Subset and prepare data
  species_data <- data[, species_list]
  species_data$group <- data[, group_name]
  
  # Initialize result matrices
  est_matrix <- as.data.frame(matrix(0, length(species_list), length(feature_list)))
  rownames(est_matrix) <- species_list
  colnames(est_matrix) <- feature_list
  
  pval_matrix <- as.data.frame(matrix(1, length(species_list), length(feature_list)))
  rownames(pval_matrix) <- species_list
  colnames(pval_matrix) <- feature_list
  
  consistency_matrix <- as.data.frame(matrix(0, length(species_list), length(feature_list)))
  rownames(consistency_matrix) <- species_list
  colnames(consistency_matrix) <- feature_list
  
  # Generate index grid
  index_grid <- expand.grid(i = 1:length(species_list), j = 1:length(feature_list))
  
  # Setup parallel backend
  # cores <- parallel::detectCores() - 1
  cores <- max(1, floor(parallel::detectCores() * 2/3))
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Export variables to all workers
  clusterExport(cl, c("compute_meta_lm2","group_name", "study_list", "index_grid"), envir = environment())
  
  # Run parallel computation
  results <- foreach(idx = 1:nrow(index_grid), .combine = rbind, 
                     .packages = c("MASS", "metafor", "robustbase", "dplyr","sfsmisc")) %dopar% {
                       i <- index_grid[idx, "i"]
                       j <- index_grid[idx, "j"]
                       species1 <- species_list[i]
                       species2 <- feature_list[j]
                       
                       est <- NA
                       pval <- 1
                       consistency <- 0
                       
                       message(paste("Running:", species1, "vs", species2))  # This will now show
                       
                       tryCatch({
                         if (species1 != species2) {
                           temp_rem <- compute_meta_lm2(data, species1, species2, group_name, study_list)
                           if (!is.null(temp_rem$model$beta)) {
                             est <- as.numeric(temp_rem$model$beta)
                             pval <- temp_rem$model$pval
                             consistency <- length(which(sign(temp_rem$df_studies$di) == sign(est))) / 
                               length(temp_rem$df_studies$di)
                           }
                         }
                       }, error = function(e) {
                         message(paste("ERROR:", species1, species2, "->", conditionMessage(e)))
                       })
                       
                       data.frame(i = i, j = j, est = est, pval = pval, consistency = consistency)
                     }
  
  
  # Stop the cluster
  stopCluster(cl)
  
  # Fill result matrices
  for (k in 1:nrow(results)) {
    i <- results$i[k]
    j <- results$j[k]
    species1 <- species_list[i]
    species2 <- feature_list[j]
    
    est_matrix[species1, species2] <- results$est[k]
    pval_matrix[species1, species2] <- results$pval[k]
    consistency_matrix[species1, species2] <- results$consistency[k]
  }
  
  # FDR correction
  qval_matrix <- apply(pval_matrix, 2, function(x) p.adjust(x, method = "fdr"))
  
  # Direction matrix
  dir_matrix <- as.data.frame(matrix(0, length(species_list), length(feature_list)))
  rownames(dir_matrix) <- species_list
  colnames(dir_matrix) <- feature_list
  
  for (i in 1:length(species_list)) {
    for (j in 1:length(feature_list)) {
      if ((qval_matrix[i, j] <= 0.05) && (consistency_matrix[i, j] >= 0.70)) {
        dir_matrix[i, j] <- sign(est_matrix[i, j])
      }
    }
  }
  
  return(list(
    est = est_matrix,
    pval = pval_matrix,
    consistency = consistency_matrix,
    qval = qval_matrix,
    dir = dir_matrix
  ))
}




########## Melt the Adjecency matrix to get the edge list
Melt_Adjacency_Matrix <- function(mat) {
  # Ensure input is a matrix
  mat <- as.matrix(mat)

  # Get indices of upper triangular (excluding diagonal)
  ut_idx <- upper.tri(mat, diag = FALSE)

  # Build edge list
  data.frame(
    Node1 = rownames(mat)[row(mat)[ut_idx]],
    Node2 = colnames(mat)[col(mat)[ut_idx]],
    Weight = mat[ut_idx],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}







########### embeddings
compute_association_consistency <- function(data,species_list,study_list,total_studies)
{
  sp_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
  rownames(sp_matrix) <- species_list
  colnames(sp_matrix) <- species_list
  
  sn_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
  rownames(sn_matrix) <- species_list
  colnames(sn_matrix) <- species_list
  
  final_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
  rownames(final_matrix) <- species_list
  colnames(final_matrix) <- species_list
  
  for(i in 1:total_studies)
  {
    study <- study_list[i]
    print(study)
    temp_corr <- corr.test(data[data$study_name == study,species_list],method="spearman",use="pairwise.complete",adjust="fdr")
    r_matrix <- apply(temp_corr$r,2,function(x)(ifelse(is.na(x),0,x)))
    p_matrix <- apply(temp_corr$p,2,function(x)(ifelse(is.na(x),1,x)))
    dir_matrix <- apply(r_matrix,2,sign) * apply(p_matrix,2,function(x)(ifelse(x<=0.1,1,0)))
    sp_matrix <- sp_matrix + apply(dir_matrix,2,function(x)(ifelse(x==1,1,0)))
    sn_matrix <- sn_matrix + apply(dir_matrix,2,function(x)(ifelse(x== -1,1,0)))
  }
  
  final_matrix <- ((sp_matrix - sn_matrix)/total_studies)*(1-(pmin(sp_matrix,sn_matrix)+1)/(pmax(sp_matrix,sn_matrix)+1))
  
  return_list <- list("final_matrix"=final_matrix,"sp_matrix"=sp_matrix,"sn_matrix"=sn_matrix)
}


## same above function but in parallel processing
compute_association_consistency_parallel <- function(
  data,
  species_list,
  study_list,
  total_studies
) {

  library(foreach)
  library(doParallel)
  library(psych)

  # ---- setup parallel backend ----
  cores <- max(1, floor(parallel::detectCores() * 0.9))
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  clusterExport(
    cl,
    varlist = c("data", "species_list"),
    envir = environment()
  )

  # ---- parallel loop (NO combine) ----
  results <- foreach(
    i = 1:total_studies,
    .packages = "psych"
  ) %dopar% {

    study <- study_list[i]

    temp_corr <- corr.test(
      data[data$study_name == study, species_list],
      method = "spearman",
      use = "pairwise.complete",
      adjust = "fdr"
    )

    r_matrix <- temp_corr$r
    p_matrix <- temp_corr$p

    r_matrix[is.na(r_matrix)] <- 0
    p_matrix[is.na(p_matrix)] <- 1

    dir_matrix <- sign(r_matrix) * (p_matrix <= 0.1)

    list(
      sp = (dir_matrix == 1),
      sn = (dir_matrix == -1)
    )
  }

  stopCluster(cl)

  # ---- initialize accumulation matrices ----
  sp_matrix <- matrix(0, length(species_list), length(species_list))
  sn_matrix <- matrix(0, length(species_list), length(species_list))

  # ---- sum across studies ----
  for (res in results) {
    sp_matrix <- sp_matrix + res$sp
    sn_matrix <- sn_matrix + res$sn
  }

  rownames(sp_matrix) <- species_list
  colnames(sp_matrix) <- species_list
  rownames(sn_matrix) <- species_list
  colnames(sn_matrix) <- species_list

  # ---- final consistency score ----
  final_matrix <- ((sp_matrix - sn_matrix) / total_studies) *
    (1 - (pmin(sp_matrix, sn_matrix) + 1) /
         (pmax(sp_matrix, sn_matrix) + 1))

  return_list <- list(
  final_matrix = as.data.frame(final_matrix),
  sp_matrix = as.data.frame(sp_matrix),
  sn_matrix = as.data.frame(sn_matrix))

return(return_list)
}



########### Add single variable from metadata df to species df
adding_variable_to_SpDf <- function(species_profile, metadata_df, variable = "study_name") {
  
  ## 1. Check variable existence
  if (!variable %in% colnames(metadata_df)) {
    stop(paste("Variable", variable, "not found in metadata_df"))
  }
  
  ## 2. Check rownames existence
  if (is.null(rownames(species_profile)) || is.null(rownames(metadata_df))) {
    stop("Both species_profile and metadata_df must have rownames")
  }
  
  ## 3. Take common samples
  common_samples <- intersect(rownames(species_profile), rownames(metadata_df))
  
  if (length(common_samples) == 0) {
    stop("No common rownames between species_profile and metadata_df")
  }
  
  ## 4. Subset BOTH data frames
  species_profile <- species_profile[common_samples, , drop = FALSE]
  metadata_df     <- metadata_df[common_samples, , drop = FALSE]
  
  ## 5. Ensure identical order
  metadata_df <- metadata_df[rownames(species_profile), , drop = FALSE]
  
  ## 6. Add variable
  species_profile[[variable]] <- metadata_df[[variable]]
  
  ## 7. Post-condition check
  stopifnot(identical(rownames(species_profile), rownames(metadata_df)))
  
  return(species_profile)
}


