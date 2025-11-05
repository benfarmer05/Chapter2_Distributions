  
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
  threshold_list <- list()  # NEW: Store thresholds
  
  for(species in names(all_results)) {
    cat("Processing", species, "...\n")
    
    # Get threshold from results
    species_result <- all_results[[species]]
    optimal_threshold <- species_result$threshold
    threshold_list[[species]] <- optimal_threshold  # NEW: Store threshold
    
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
  
  # NEW: Create threshold summary table
  threshold_summary <- data.frame(
    Species = names(threshold_list),
    Threshold = unlist(threshold_list),
    stringsAsFactors = FALSE
  )
  
  # Define display order and names (GROUPED BY SUSCEPTIBILITY)
  species_codes <- c("agaricia", "madracis", "porites", "siderastrea", "stephanocoenia",  # LS
                     "montastraea", "orbicella",  # MS
                     "colpophyllia", "dendrogyra", "dichocoenia", "diploria",  # HS
                     "eusmilia", "meandrina", "mycetophyllia", "pseudodiploria")
  
  display_names <- c("Agaricia", "Madracis", "Porites", "Siderastrea", "S. intersepta",  # LS
                     "M. cavernosa", "Orbicella",  # MS
                     "C. natans", "D. cylindrus", "D. stokesii", "D. labyrinthiformis",  # HS
                     "E. fastigiata", "Meandrina", "Mycetophyllia", "Pseudodiploria")
  
  sus_groups <- c("LS", "LS", "LS", "LS", "LS",  # LS group
                  "MS", "MS",  # MS group
                  "HS", "HS", "HS", "HS", "HS", "HS", "HS", "HS")  # HS group
  
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
  
  # NEW: Reorder and format threshold results
  threshold_summary_df <- threshold_summary %>%
    mutate(Species_lower = tolower(Species)) %>%
    left_join(species_mapping, by = c("Species_lower" = "code")) %>%
    select(Taxon = display, Threshold) %>%
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
  
  # NEW: Add placeholder rows for threshold means
  threshold_placeholder <- data.frame(Taxon = "Mean", Threshold = NA)
  threshold_grand_placeholder <- data.frame(Taxon = "Grand Mean", Threshold = NA)
  
  # NEW: Combine thresholds with means
  threshold_with_means <- rbind(
    threshold_summary_df[1:5, ],
    threshold_placeholder,
    threshold_summary_df[6:7, ],
    threshold_placeholder,
    threshold_summary_df[8:15, ],
    threshold_placeholder,
    threshold_grand_placeholder
  )
  
  # Combine with group means (LS rows 1-5, MS rows 6-7, HS rows 8-15)
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
  
  # NEW: Create combined table with Threshold column
  combined_results <- cv_summary_df %>%
    left_join(threshold_summary_df, by = "Taxon") %>%  # Add threshold
    left_join(overall_summary_df, by = "Taxon", suffix = c("_cv", "_all")) %>%
    select(
      Taxon,
      Threshold,  # NEW: Threshold column
      AUC_All = AUC_all,
      AUC_Test = AUC_cv,
      Sensitivity_All = Sensitivity_all,
      Sensitivity_Test = Sensitivity_cv,
      Specificity_All = Specificity_all,
      Specificity_Test = Specificity_cv,
      Accuracy_All = Accuracy_all,
      Accuracy_Test = Accuracy_cv
    )
  
  # ========== VALIDATION CHECKS ==========
  cat("\n========================================\n")
  cat("DATA VALIDATION CHECKS\n")
  cat("========================================\n\n")
  
  # Check 1: Verify no NAs in final table (excluding Threshold which may have NAs for means)
  na_check <- sapply(combined_results[, -c(1, 2)], function(x) sum(is.na(x)))
  if(any(na_check > 0)) {
    cat("⚠ WARNING: NA values detected in results!\n")
    print(na_check[na_check > 0])
  } else {
    cat("✓ Check 1 PASSED: No NA values in results\n")
  }
  
  # Check 2: Verify all species from results are in final table
  species_in_results <- names(all_results)
  species_in_table <- tolower(gsub("\\.", " ", combined_results$Taxon))
  species_in_table <- sapply(strsplit(species_in_table, " "), function(x) x[1])
  
  # Create mapping for comparison
  name_map <- setNames(species_codes, display_names)
  species_codes_in_table <- sapply(combined_results$Taxon, function(x) name_map[x])
  
  missing_from_table <- setdiff(species_in_results, species_codes_in_table)
  if(length(missing_from_table) > 0) {
    cat("⚠ WARNING: Species in results but missing from table:\n")
    print(missing_from_table)
  } else {
    cat("✓ Check 2 PASSED: All species from results are in final table\n")
  }
  
  # Check 3: Verify correct data mapping by checking a few species
  cat("\n✓ Check 3: Verifying data integrity for sample species...\n")
  sample_species <- c("agaricia", "orbicella", "pseudodiploria")
  
  for(sp in sample_species) {
    if(sp %in% names(cv_results_list) && sp %in% names(overall_results_list)) {
      # Get display name
      display_name <- species_mapping$display[species_mapping$code == sp]
      
      # Get values from original lists
      cv_auc_original <- cv_results_list[[sp]]["AUC"]
      overall_auc_original <- overall_results_list[[sp]]["AUC"]
      
      # Get values from final table
      table_row <- combined_results[combined_results$Taxon == display_name, ]
      cv_auc_table <- as.numeric(table_row$AUC_Test)
      overall_auc_table <- as.numeric(table_row$AUC_All)
      
      # Compare
      cv_match <- abs(cv_auc_original - cv_auc_table) < 0.001
      overall_match <- abs(overall_auc_original - overall_auc_table) < 0.001
      
      if(cv_match && overall_match) {
        cat(sprintf("  ✓ %s: CV AUC = %.3f, Overall AUC = %.3f (CORRECT)\n", 
                    display_name, cv_auc_table, overall_auc_table))
      } else {
        cat(sprintf("  ✗ %s: MISMATCH DETECTED!\n", display_name))
        cat(sprintf("    Original CV: %.3f vs Table: %.3f\n", cv_auc_original, cv_auc_table))
        cat(sprintf("    Original Overall: %.3f vs Table: %.3f\n", overall_auc_original, overall_auc_table))
      }
    }
  }
  
  # Check 4: Verify species order matches expected grouping
  cat("\n✓ Check 4: Verifying species order (grouped by susceptibility)...\n")
  expected_order <- display_names
  actual_order <- combined_results$Taxon
  if(all(expected_order == actual_order)) {
    cat("  ✓ Species order is correct (LS → MS → HS)\n")
  } else {
    cat("  ✗ WARNING: Species order does not match expected!\n")
    cat("  Expected:", paste(expected_order, collapse=", "), "\n")
    cat("  Actual:", paste(actual_order, collapse=", "), "\n")
  }
  
  # Check 5: Verify AUC values are in valid range [0, 1]
  cat("\n✓ Check 5: Verifying AUC values are in valid range...\n")
  auc_values <- c(combined_results$AUC_All, combined_results$AUC_Test)
  if(all(auc_values >= 0 & auc_values <= 1)) {
    cat("  ✓ All AUC values are in valid range [0, 1]\n")
  } else {
    cat("  ✗ WARNING: Some AUC values are out of range!\n")
  }
  
  # NEW: Check 6: Verify thresholds are in valid range [0, 1]
  cat("\n✓ Check 6: Verifying threshold values are in valid range...\n")
  threshold_values <- combined_results$Threshold[!is.na(combined_results$Threshold)]
  if(all(threshold_values >= 0 & threshold_values <= 1)) {
    cat("  ✓ All threshold values are in valid range [0, 1]\n")
  } else {
    cat("  ✗ WARNING: Some threshold values are out of range!\n")
  }
  
  cat("\n========================================\n")
  cat("VALIDATION COMPLETE\n")
  cat("========================================\n")
  
  # Format for display (2 decimal places for metrics, 3 for threshold)
  display_table <- combined_results
  display_table$Threshold <- ifelse(is.na(display_table$Threshold), 
                                    "", 
                                    sprintf("%.3f", display_table$Threshold))
  display_table[, 3:10] <- lapply(display_table[, 3:10], function(x) sprintf("%.2f", x))
  
  # Print as one block for proper alignment
  {
    cat("\n========================================\n")
    cat("MODEL PERFORMANCE METRICS\n")
    cat("========================================\n\n")
    cat(sprintf("%-20s %9s %17s %17s %17s %17s\n",
                "Taxon", "Threshold", "AUC", "Sensitivity", "Specificity", "Accuracy"))
    cat(sprintf("%-20s %9s %8s %8s %8s %8s %8s %8s %8s %8s\n",
                "", "", "All", "Test", "All", "Test", "All", "Test", "All", "Test"))
    cat(paste(rep("-", 110), collapse = ""), "\n")
    
    for(i in 1:nrow(display_table)) {
      cat(sprintf("%-20s %9s %8s %8s %8s %8s %8s %8s %8s %8s\n",
                  display_table$Taxon[i],
                  display_table$Threshold[i],
                  display_table$AUC_All[i],
                  display_table$AUC_Test[i],
                  display_table$Sensitivity_All[i],
                  display_table$Sensitivity_Test[i],
                  display_table$Specificity_All[i],
                  display_table$Specificity_Test[i],
                  display_table$Accuracy_All[i],
                  display_table$Accuracy_Test[i]))
    }
    
    cat("\n✓ Analysis complete\n")
    cat("\nNote: 'All' = full dataset performance; 'Test' = 5-fold CV average\n")
    cat("      Threshold = optimal classification probability threshold\n")
  }
  
  # Save table to CSV (with appropriate decimal places)
  csv_output <- combined_results
  csv_output$Threshold <- round(csv_output$Threshold, 3)
  csv_output[, 3:10] <- lapply(csv_output[, 3:10], function(x) round(x, 2))
  output_csv <- here("output", "output_figures_tables", "cv_performance_table.csv")
  write.csv(csv_output, output_csv, row.names = FALSE)
  cat("\n✓ Table saved to:", output_csv, "\n")