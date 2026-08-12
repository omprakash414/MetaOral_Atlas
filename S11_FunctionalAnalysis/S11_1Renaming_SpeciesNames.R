# This script is to update the species names in the functional profiles based on the mapping provided in the mismatched species name dataframes. It will also log the changes made to a text file for reference.

setwd("/storage/omprakash/MetaOral_Analysis/MetaOral_Data_Analysis/S11_FunctionalAnalysis")

library(dplyr)

## import the mismatched species name dataframes
load("Mismatched_names_NCBI_taxonomy.RData")

## import the functional profile dfs
load("MeanFunction_ConservationProfile_OverallGroups.RData")

update_species_names <- function(df, df_name) {
  
  species_names <- names(df)
  
  # ---------------- Helper Functions ----------------
  contains_numbers <- function(x) grepl("[0-9]", x)
  
  apply_mapping <- function(names_vec, mapping_df, old_col, new_col, mark = "") {
    # Standardize new names
    mapping_df[[new_col]] <- gsub(" ", "_", mapping_df[[new_col]])
    
    # Filter out names with numbers
    valid_idx <- !contains_numbers(mapping_df[[old_col]]) &
      !contains_numbers(mapping_df[[new_col]])
    mapping_df <- mapping_df[valid_idx, ]
    
    # Create named vector: new_name → old_name
    map_vec <- setNames(mapping_df[[old_col]], mapping_df[[new_col]])
    
    # Identify which species names to replace
    replace_idx <- names_vec %in% names(map_vec) & !contains_numbers(names_vec)
    
    # Apply mapping
    names_vec[replace_idx] <- map_vec[names_vec[replace_idx]]
    
    # Optionally mark names
    if (mark != "") {
      names_vec[replace_idx] <- paste0(names_vec[replace_idx], mark)
    }
    
    return(names_vec)
  }
  
  # ---------------- Step 1: 14 → 24 ----------------
  species_names <- apply_mapping(species_names, mismatched_14_24, old_col = "names_14", new_col = "names_24")
  
  # ---------------- Step 2: 18 → 24 ----------------
  species_names <- apply_mapping(species_names, mismatched_18_24, old_col = "names_18", new_col = "names_24", mark = "%%%")
  
  # ---------------- Step 3: 14 → 18 ----------------
  species_names <- apply_mapping(species_names, mismatched_14_18, old_col = "names_14", new_col = "names_18")
  
  # ---------------- Final Formatting Cleanup ----------------
  species_names <- gsub(" ", "_", species_names)
  species_names <- gsub("\\[|\\]", "", species_names)
  
  # ---------------- Log Name Changes ----------------
  changes <- data.frame(Original = names(df),Updated  = species_names,stringsAsFactors = FALSE)
  changes <- changes[changes$Original != changes$Updated, ]
  
  write.table(
    changes,
    file = paste0("SpeciesRenamed_", df_name, "_log.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("\n***********************\n")
  cat("Changes for", df_name, ":\n")
  print(changes)
  cat("\n***********************\n")
  
  # ---------------- Update Dataframe and Return ----------------
  names(df) <- species_names
  names(df) <- gsub("%", "", names(df))  # Remove any marking symbols
  
  return(df)
}



################# Update species names in both dataframes
#################
species_conservation_profile_t <- data.frame(t(species_conservation_profile))
species_conservation_profile_t <- update_species_names(species_conservation_profile_t, "species_conservation_profile")

MeanFunctionalProfile_groups_t <- data.frame(t(MeanFunctionalProfile_groups))
MeanFunctionalProfile_groups_t <- update_species_names(MeanFunctionalProfile_groups_t, "MeanFunctionalProfile_groups")

rm(MeanFunctionalProfile_groups,species_conservation_profile)

MeanFunctionalProfile_groups <- data.frame(t(MeanFunctionalProfile_groups_t))
species_conservation_profile <- data.frame(t(species_conservation_profile_t))

rm(species_conservation_profile_t,MeanFunctionalProfile_groups_t)
gc()  # Clean up memory


# Changes for MeanFunctionalProfile_groups :
#                Original                   Updated
# 4700 Sarcina_ventriculi Clostridium_ventriculi%%%


save(MeanFunctionalProfile_groups, species_conservation_profile, annotation_df, file = "S11_1MeanFunctional_ConservationProfile_OverallGroups.RData")
save.image("S11_1Renaming_SpeciesNames_workspace.RData")

