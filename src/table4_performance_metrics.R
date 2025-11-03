# .rs.restartR(clean = TRUE)
rm(list=ls())

library(here)
library(mgcv)
library(pROC)
library(dplyr)

################################## Setup ##################################

cat("========================================\n")
cat("Presence Model Cross-Validation Summary\n")
cat("========================================\n\n")

# Load combined data
load(here("output", "all_combined_data.rda"))

cat("Reading model outputs...\n")

# Get all results files
results_files <- list.files(here('output/output_maps'), 
                            pattern = '^results_light_.*\\.rds$', 
                            full.names = TRUE)

# Read all results and unwrap rasters
all_results <- lapply(results_files, function(f) {
  result <- readRDS(f)
  result$raster <- terra::unwrap(result$raster)
  result$raster_raw <- terra::unwrap(result$raster_raw)
  return(result)
})
names(all_results) <- gsub('results_light_|\\.rds', '', basename(results_files))

cat("Loaded", length(all_results), "species models\n\n")

# Save information for exporting new objects later
existing_objects <- ls(envir = .GlobalEnv)

################################## Helper Functions ##################################

# Calculate classification metrics
calculate_metrics <- function(observed, predicted_prob, threshold) {
  predicted_binary <- ifelse(predicted_prob >= threshold, 1, 0)
  
  # Confusion matrix components
  tp <- sum(observed == 1 & predicted_binary == 1)
  tn <- sum(observed == 0 & predicted_binary == 0)
  fp <- sum(observed == 0 & predicted_binary == 1)
  fn <- sum(observed == 1 & predicted_binary == 0)
  
  # Metrics
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  
  # AUC
  roc_obj <- roc(observed, predicted_prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  return(c(AUC = auc_val, 
           Sensitivity = sensitivity, 
           Specificity = specificity, 
           Accuracy = accuracy))
}

################################## Cross-Validation ##################################

cat("Running 5-fold cross-validation...\n\n")

# Initialize results storage
cv_results_list <- list()
overall_results_list <- list()

for(species in names(all_results)) {
  cat("Processing", species, "...\n")
  
  # Get threshold from results
  species_result <- all_results[[species]]
  optimal_threshold <- species_result$threshold
  
  # Load the full model file to get presence_model and model_data
  model_file <- here("output", "output_GAMs", paste0(species, "_models.rds"))
  if(!file.exists(model_file)) {
    cat("  Warning: Model file not found, skipping\n")
    next
  }
  
  species_models <- readRDS(model_file)
  presence_model <- species_models$presence_model
  model_data <- species_models$model_data
  
  # Prepare data
  all_vars <- all.vars(formula(presence_model))
  complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  
  # Create binary presence variable
  complete_data$presence <- ifelse(complete_data$cover > 0, 1, 0)
  
  cat("  Data: N =", nrow(complete_data), 
      "| Presence =", sum(complete_data$presence == 1),
      "| Absence =", sum(complete_data$presence == 0), "\n")
  
  # ========== OVERALL METRICS (full dataset) ==========
  overall_pred <- predict(presence_model, newdata = complete_data, type = "response")
  overall_metrics <- calculate_metrics(complete_data$presence, overall_pred, optimal_threshold)
  overall_results_list[[species]] <- overall_metrics
  
  # ========== CROSS-VALIDATION ==========
  n <- nrow(complete_data)
  fold_ids <- sample(rep(1:5, length.out = n))
  
  cv_metrics <- matrix(NA, nrow = 5, ncol = 4)
  colnames(cv_metrics) <- c("AUC", "Sensitivity", "Specificity", "Accuracy")
  
  for(fold in 1:5) {
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    # Fit model on training data
    presence_fit <- gam(formula(presence_model), 
                        data = complete_data[train_idx, ], 
                        family = binomial())
    
    # Predict on test data
    test_pred <- predict(presence_fit, 
                         newdata = complete_data[test_idx, ], 
                         type = "response")
    test_obs <- complete_data$presence[test_idx]
    
    # Calculate metrics using the optimal threshold from full model
    cv_metrics[fold, ] <- calculate_metrics(test_obs, test_pred, optimal_threshold)
  }
  
  # Average across folds
  cv_results_list[[species]] <- colMeans(cv_metrics)
  
  cat("  CV complete\n")
}

cat("\n✓ Cross-validation complete for all species\n\n")

################################## Format Results ##################################

# Create CV summary table
cv_summary <- do.call(rbind, cv_results_list)
cv_summary_df <- as.data.frame(cv_summary)
cv_summary_df$Species <- rownames(cv_summary_df)
cv_summary_df <- cv_summary_df[, c("Species", "AUC", "Sensitivity", "Specificity", "Accuracy")]

# Create overall summary table
overall_summary <- do.call(rbind, overall_results_list)
overall_summary_df <- as.data.frame(overall_summary)
overall_summary_df$Species <- rownames(overall_summary_df)
overall_summary_df <- overall_summary_df[, c("Species", "AUC", "Sensitivity", "Specificity", "Accuracy")]

# Define display order and names
species_codes <- c("agaricia", "colpophyllia", "dendrogyra", "dichocoenia", 
                   "diploria", "eusmilia", "madracis", "meandrina", 
                   "montastraea", "mycetophyllia", "orbicella", "porites", 
                   "pseudodiploria", "siderastrea", "stephanocoenia")

display_names <- c("Agaricia", "C. natans", "D. cylindrus", "D. stokesii", 
                   "D. labyrinthiformis", "E. fastigiata", "Madracis", 
                   "Meandrina", "M. cavernosa", "Mycetophyllia", "Orbicella", 
                   "Porites", "Pseudodiploria", "Siderastrea", "S. intersepta")

sus_groups <- c("LS", "HS", "HS", "HS", "HS", "HS", "LS", "HS", 
                "MS", "HS", "MS", "LS", "HS", "LS", "LS")

# Create mapping
species_mapping <- data.frame(
  code = species_codes,
  display = display_names,
  sus = sus_groups,
  stringsAsFactors = FALSE
)

# Reorder and format CV results
cv_summary_df <- cv_summary_df %>%
  mutate(Species_lower = tolower(Species)) %>%
  left_join(species_mapping, by = c("Species_lower" = "code")) %>%
  select(Taxon = display, Sus = sus, AUC, Sensitivity, Specificity, Accuracy) %>%
  arrange(match(Taxon, display_names))

# Reorder and format overall results
overall_summary_df <- overall_summary_df %>%
  mutate(Species_lower = tolower(Species)) %>%
  left_join(species_mapping, by = c("Species_lower" = "code")) %>%
  select(Taxon = display, Sus = sus, AUC, Sensitivity, Specificity, Accuracy) %>%
  arrange(match(Taxon, display_names))

# Calculate group means for CV
ls_indices <- which(cv_summary_df$Sus == "LS")
ms_indices <- which(cv_summary_df$Sus == "MS")
hs_indices <- which(cv_summary_df$Sus == "HS")

cv_ls_mean <- data.frame(
  Taxon = "Mean", Sus = "",
  AUC = mean(cv_summary_df$AUC[ls_indices]),
  Sensitivity = mean(cv_summary_df$Sensitivity[ls_indices]),
  Specificity = mean(cv_summary_df$Specificity[ls_indices]),
  Accuracy = mean(cv_summary_df$Accuracy[ls_indices])
)

cv_ms_mean <- data.frame(
  Taxon = "Mean", Sus = "",
  AUC = mean(cv_summary_df$AUC[ms_indices]),
  Sensitivity = mean(cv_summary_df$Sensitivity[ms_indices]),
  Specificity = mean(cv_summary_df$Specificity[ms_indices]),
  Accuracy = mean(cv_summary_df$Accuracy[ms_indices])
)

cv_hs_mean <- data.frame(
  Taxon = "Mean", Sus = "",
  AUC = mean(cv_summary_df$AUC[hs_indices]),
  Sensitivity = mean(cv_summary_df$Sensitivity[hs_indices]),
  Specificity = mean(cv_summary_df$Specificity[hs_indices]),
  Accuracy = mean(cv_summary_df$Accuracy[hs_indices])
)

cv_grand_mean <- data.frame(
  Taxon = "Grand Mean", Sus = "",
  AUC = mean(cv_summary_df$AUC),
  Sensitivity = mean(cv_summary_df$Sensitivity),
  Specificity = mean(cv_summary_df$Specificity),
  Accuracy = mean(cv_summary_df$Accuracy)
)

# Calculate group means for overall
overall_ls_mean <- data.frame(
  Taxon = "Mean", Sus = "",
  AUC = mean(overall_summary_df$AUC[ls_indices]),
  Sensitivity = mean(overall_summary_df$Sensitivity[ls_indices]),
  Specificity = mean(overall_summary_df$Specificity[ls_indices]),
  Accuracy = mean(overall_summary_df$Accuracy[ls_indices])
)

overall_ms_mean <- data.frame(
  Taxon = "Mean", Sus = "",
  AUC = mean(overall_summary_df$AUC[ms_indices]),
  Sensitivity = mean(overall_summary_df$Sensitivity[ms_indices]),
  Specificity = mean(overall_summary_df$Specificity[ms_indices]),
  Accuracy = mean(overall_summary_df$Accuracy[ms_indices])
)

overall_hs_mean <- data.frame(
  Taxon = "Mean", Sus = "",
  AUC = mean(overall_summary_df$AUC[hs_indices]),
  Sensitivity = mean(overall_summary_df$Sensitivity[hs_indices]),
  Specificity = mean(overall_summary_df$Specificity[hs_indices]),
  Accuracy = mean(overall_summary_df$Accuracy[hs_indices])
)

overall_grand_mean <- data.frame(
  Taxon = "Grand Mean", Sus = "",
  AUC = mean(overall_summary_df$AUC),
  Sensitivity = mean(overall_summary_df$Sensitivity),
  Specificity = mean(overall_summary_df$Specificity),
  Accuracy = mean(overall_summary_df$Accuracy)
)

# Combine with group means
cv_with_means <- rbind(
  cv_summary_df[1:5, ],
  cv_ls_mean,
  cv_summary_df[6:7, ],
  cv_ms_mean,
  cv_summary_df[8:15, ],
  cv_hs_mean,
  cv_grand_mean
)

overall_with_means <- rbind(
  overall_summary_df[1:5, ],
  overall_ls_mean,
  overall_summary_df[6:7, ],
  overall_ms_mean,
  overall_summary_df[8:15, ],
  overall_hs_mean,
  overall_grand_mean
)

# Format for display
cv_table <- cv_with_means
cv_table[, 3:6] <- lapply(cv_table[, 3:6], function(x) sprintf("%.2f", x))

overall_table <- overall_with_means
overall_table[, 3:6] <- lapply(overall_table[, 3:6], function(x) sprintf("%.2f", x))

################################## Print Results ##################################

cat("\n========================================\n")
cat("OUT-OF-SAMPLE (5-Fold CV) METRICS\n")
cat("========================================\n\n")
print(cv_table, row.names = FALSE)

cat("\n\n========================================\n")
cat("OVERALL (Full Dataset) METRICS\n")
cat("========================================\n\n")
print(overall_table, row.names = FALSE)

cat("\n✓ Analysis complete\n")