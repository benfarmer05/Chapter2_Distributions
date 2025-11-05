  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(mgcv)
  library(gratia)
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(gam.hp)
  library(foreach)
  library(doParallel)
  
  ################################## SETTINGS ##################################
  
  # Toggle for including partial residuals in draw plots
  INCLUDE_RESIDUALS <- FALSE  # Set to FALSE to show only smooth functions
  
  # Percentile limits for y-axis (only used if INCLUDE_RESIDUALS = TRUE)
  # Set to c(0, 1) for no trimming, or e.g. c(0.02, 0.98) to trim 2% from each tail
  Y_PERCENTILES <- c(0.02, 0.98)  # Lower and upper percentiles
  
  # Variable importance method
  # Options: "chisq" (fast, uses chi-squared from summary) or "gam.hp" (slow but rigorous)
  IMPORTANCE_METHOD <- "chisq"  # Set to "chisq" or "gam.hp"
  
  # Parallel processing settings
  USE_PARALLEL <- TRUE  # Set to TRUE to process species in parallel
  N_CORES <- 4  # Number of cores to use (adjust based on your system)
  
  ################################## Setup ##################################
  
  cat("========================================\n")
  cat("GAM Predictor Analysis\n")
  cat("Variable importance method:", IMPORTANCE_METHOD, "\n")
  cat("========================================\n\n")
  
  cat("Plot settings:\n")
  cat("  Include residuals:", INCLUDE_RESIDUALS, "\n")
  if(INCLUDE_RESIDUALS) {
    cat("  Y-axis percentile limits:", Y_PERCENTILES[1], "to", Y_PERCENTILES[2], "\n")
  }
  cat("  Importance method:", IMPORTANCE_METHOD, "\n")
  cat("  Parallel processing:", USE_PARALLEL, "\n")
  if(USE_PARALLEL) {
    cat("  Number of cores:", N_CORES, "\n")
  }
  cat("\n")
  
  # Create output directory
  output_dir <- here("output", "output_figures_tables", "predictorsummaries")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Species mapping
  species_codes <- c("agaricia", "madracis", "porites", "siderastrea", "stephanocoenia",
                     "montastraea", "orbicella",
                     "colpophyllia", "dendrogyra", "dichocoenia", "diploria",
                     "eusmilia", "meandrina", "mycetophyllia", "pseudodiploria")
  
  display_names <- c("Agaricia", "Madracis", "Porites", "Siderastrea", "S. intersepta",
                     "M. cavernosa", "Orbicella",
                     "C. natans", "D. cylindrus", "D. stokesii", "D. labyrinthiformis",
                     "E. fastigiata", "Meandrina", "Mycetophyllia", "Pseudodiploria")
  
  code_4letter <- c("AGAR", "MADR", "PORI", "SIDE", "SINT",
                    "MCAV", "ORBI",
                    "CNAT", "DCYL", "DSTO", "DLAB",
                    "EFAS", "MEAN", "MYCE", "PSEU")
  
  species_mapping <- data.frame(code = species_codes, 
                                display = display_names,
                                code_4letter = code_4letter,
                                stringsAsFactors = FALSE)
  
  ################################## Helper Functions ##################################
  
  # Clean predictor names from GAM output
  clean_predictor_name <- function(name) {
    # Remove s() wrapper
    name <- gsub("^s\\(([^)]+)\\)$", "\\1", name)
    
    # Map full predictor names to shortened codes
    predictor_map <- c(
      "depth_bathy" = "BAT",
      "mean_dir" = "DIR",
      "mean_SST" = "TPM",
      "slope" = "SLP",
      "TPI" = "TPI",
      "mean_kd490" = "490",
      "mean_PAR" = "PAR",
      "max_BOV" = "BOV",
      "dist_to_deep" = "DEP",
      "complexity" = "COM",
      "range_SST" = "TPR",
      "planform_curv" = "PCV",
      "mean_spm" = "SPM",
      "SAPA" = "SAP",
      "mean_chla" = "CHL",
      "max_Hsig" = "HSX",
      "mean_Hsig" = "HSM",
      "VRM" = "VRM",
      "aspect" = "ASP",
      "dist_to_land" = "LND"
    )
    
    # Return mapped name if exists, otherwise return original
    if(name %in% names(predictor_map)) {
      return(predictor_map[name])
    } else {
      return(name)
    }
  }
  
  # Extract predictor importance from GAM
  get_predictor_importance <- function(model, method = "chisq") {
    
    if(method == "gam.hp") {
      # Use gam.hp for hierarchical partitioning (slow but rigorous)
      tryCatch({
        hp_result <- gam.hp(model, type = "dev")
        
        # Extract the hierarchical.partitioning table
        hp_table <- hp_result$hierarchical.partitioning
        
        # Get predictor names and their Individual contributions
        pred_names <- rownames(hp_table)
        ind_values <- hp_table[, "Individual"]
        
        # Clean predictor names
        names(ind_values) <- sapply(pred_names, clean_predictor_name)
        
        # Convert to percentage
        rel_importance <- ind_values * 100
        
        return(rel_importance)
      }, error = function(e) {
        warning(paste("gam.hp failed:", e$message))
        return(numeric(0))
      })
      
    } else if(method == "chisq") {
      # Use chi-squared from summary (fast)
      s <- summary(model)
      
      # Extract chi-squared values for smooth terms
      if(!is.null(s$s.table) && nrow(s$s.table) > 0) {
        chi_sq <- s$s.table[, "Chi.sq"]
        names(chi_sq) <- sapply(rownames(s$s.table), clean_predictor_name)
        
        # Calculate relative importance (as percentage)
        total_chi <- sum(chi_sq, na.rm = TRUE)
        if(total_chi > 0) {
          rel_importance <- (chi_sq / total_chi) * 100
          return(rel_importance)
        }
      }
      
      # Return empty if no smooth terms
      return(numeric(0))
      
    } else {
      stop("Invalid importance method. Use 'chisq' or 'gam.hp'")
    }
  }
  
  # Format importance table with highlighting
  format_importance_table <- function(importance_df) {
    numeric_cols <- names(importance_df)[names(importance_df) != "Species"]
    
    # For each species, identify which predictor is #1 and its value
    top_predictors <- character(nrow(importance_df))
    max_vals <- numeric(nrow(importance_df))
    
    for(i in 1:nrow(importance_df)) {
      x <- as.numeric(importance_df[i, numeric_cols])
      if(all(is.na(x))) {
        top_predictors[i] <- NA
        max_vals[i] <- NA
      } else {
        max_idx <- which.max(x)
        top_predictors[i] <- numeric_cols[max_idx]
        max_vals[i] <- x[max_idx]
      }
    }
    
    # Count how many species have each predictor as #1
    predictor_counts <- table(top_predictors[!is.na(top_predictors)])
    
    # For each predictor that's #1 for at least one species,
    # find the maximum value among its "resident taxa"
    predictor_max_vals <- sapply(names(predictor_counts), function(pred) {
      species_with_this_top <- which(top_predictors == pred)
      max(importance_df[species_with_this_top, pred], na.rm = TRUE)
    })
    
    # Order predictors by: 
    # 1) Count (descending) - how many species have it as #1
    # 2) Max value among residents (descending) - tiebreaker
    # 3) Alphabetically - final tiebreaker
    predictors_with_counts <- data.frame(
      predictor = names(predictor_counts),
      count = as.numeric(predictor_counts),
      max_val = predictor_max_vals,
      stringsAsFactors = FALSE
    )
    predictors_with_counts <- predictors_with_counts[order(-predictors_with_counts$count,
                                                           -predictors_with_counts$max_val,
                                                           predictors_with_counts$predictor), ]
    
    # Predictors that were never #1
    never_top <- setdiff(numeric_cols, names(predictor_counts))
    
    if(length(never_top) > 0) {
      # For never-#1 predictors, find max value across all species
      never_top_max_vals <- sapply(never_top, function(pred) {
        max(importance_df[, pred], na.rm = TRUE)
      })
      
      # Order by max value (descending), then alphabetically
      never_top_order <- order(-never_top_max_vals, never_top)
      never_top <- never_top[never_top_order]
    }
    
    # Final column order
    col_order <- c(predictors_with_counts$predictor, never_top)
    
    # Create a sort key for rows:
    # 1) Which predictor is #1 (ordered by our col_order)
    # 2) The value of that predictor (descending)
    predictor_rank <- match(top_predictors, col_order)
    
    # Sort rows by: predictor rank (ascending), then max value (descending)
    row_order <- order(predictor_rank, -max_vals, na.last = TRUE)
    
    # Apply row sorting
    importance_df <- importance_df[row_order, ]
    
    # Create formatted version with no decimals (rounded to integer)
    formatted <- importance_df
    formatted[numeric_cols] <- lapply(formatted[numeric_cols], function(x) {
      ifelse(is.na(x), "", sprintf("%.0f", round(x)))
    })
    
    # Reorder columns
    formatted <- formatted[, c("Species", col_order)]
    importance_df <- importance_df[, c("Species", col_order)]
    
    return(list(formatted = formatted, importance = importance_df))
  }
  
  ################################## Process Models ##################################
  
  cat("Processing models...\n\n")
  
  # Function to process a single species
  process_species <- function(i, species_codes, display_names, code_4letter, 
                              output_dir, INCLUDE_RESIDUALS, Y_PERCENTILES, IMPORTANCE_METHOD) {
    species <- species_codes[i]
    display_name <- display_names[i]
    code_4 <- code_4letter[i]
    
    cat("Processing", display_name, "...\n")
    
    model_file <- here("output", "output_GAMs", paste0(species, "_models.rds"))
    if(!file.exists(model_file)) {
      cat("  Warning: Model file not found, skipping\n")
      return(list(
        species = display_name,
        code_4 = code_4,
        skipped = TRUE,
        presence_imp = NULL,
        abundance_imp = NULL
      ))
    }
    
    models <- readRDS(model_file)
    
    # === PRESENCE MODEL ===
    presence_model <- models$presence_model
    
    # Draw plot
    pdf(file.path(output_dir, paste0(species, "_presence_draws.pdf")), width = 12, height = 8)
    
    if(INCLUDE_RESIDUALS) {
      pres <- residuals(presence_model, type = "working")
      y_lower <- quantile(pres, Y_PERCENTILES[1], na.rm = TRUE)
      y_upper <- quantile(pres, Y_PERCENTILES[2], na.rm = TRUE)
      print(draw(presence_model, residuals = TRUE, rug = FALSE) & 
              ylim(y_lower, y_upper))
    } else {
      print(draw(presence_model, residuals = FALSE))
    }
    
    dev.off()
    
    # Predictor importance
    presence_imp <- get_predictor_importance(presence_model, method = IMPORTANCE_METHOD)
    
    # === ABUNDANCE MODEL ===
    abundance_model <- models$abundance_model
    
    # Draw plot
    pdf(file.path(output_dir, paste0(species, "_abundance_draws.pdf")), width = 12, height = 8)
    
    if(INCLUDE_RESIDUALS) {
      pres <- residuals(abundance_model, type = "working")
      y_lower <- quantile(pres, Y_PERCENTILES[1], na.rm = TRUE)
      y_upper <- quantile(pres, Y_PERCENTILES[2], na.rm = TRUE)
      print(draw(abundance_model, residuals = TRUE, rug = FALSE) & 
              ylim(y_lower, y_upper))
    } else {
      print(draw(abundance_model, residuals = FALSE))
    }
    
    dev.off()
    
    # Predictor importance
    abundance_imp <- get_predictor_importance(abundance_model, method = IMPORTANCE_METHOD)
    
    cat("  Complete\n")
    
    return(list(
      species = display_name,
      code_4 = code_4,
      skipped = FALSE,
      presence_imp = presence_imp,
      abundance_imp = abundance_imp
    ))
  }
  
  # Process models (parallel or sequential)
  if(USE_PARALLEL) {
    # Set up parallel backend
    cl <- makeCluster(N_CORES)
    registerDoParallel(cl)
    
    cat("Running in parallel on", N_CORES, "cores...\n\n")
    
    # Export necessary objects and functions to workers
    clusterExport(cl, c("here", "output_dir", "INCLUDE_RESIDUALS", "Y_PERCENTILES",
                        "get_predictor_importance", "clean_predictor_name", "IMPORTANCE_METHOD"))
    clusterEvalQ(cl, {
      library(here)
      library(mgcv)
      library(gratia)
      library(gam.hp)
    })
    
    # Process species in parallel
    results <- foreach(i = seq_along(species_codes), 
                       .packages = c("here", "mgcv", "gratia", "gam.hp")) %dopar% {
                         process_species(i, species_codes, display_names, code_4letter,
                                         output_dir, INCLUDE_RESIDUALS, Y_PERCENTILES, IMPORTANCE_METHOD)
                       }
    
    # Stop cluster
    stopCluster(cl)
    
  } else {
    # Sequential processing
    cat("Running sequentially...\n\n")
    results <- lapply(seq_along(species_codes), function(i) {
      process_species(i, species_codes, display_names, code_4letter,
                      output_dir, INCLUDE_RESIDUALS, Y_PERCENTILES, IMPORTANCE_METHOD)
    })
  }
  
  # Extract results
  presence_importance_list <- list()
  abundance_importance_list <- list()
  skipped_species <- character(0)
  
  for(result in results) {
    if(result$skipped) {
      skipped_species <- c(skipped_species, result$species)
    } else {
      presence_importance_list[[result$code_4]] <- result$presence_imp
      abundance_importance_list[[result$code_4]] <- result$abundance_imp
    }
  }
  
  cat("\n✓ All models processed\n")
  
  # Report skipped species
  if(length(skipped_species) > 0) {
    cat("\n⚠ Skipped species (model files not found):\n")
    for(sp in skipped_species) {
      cat("  -", sp, "\n")
    }
  }
  cat("\n")
  
  ################################## Create Importance Tables ##################################
  
  cat("Creating predictor importance tables...\n")
  
  # Get all unique predictors
  all_presence_preds <- unique(unlist(lapply(presence_importance_list, names)))
  all_abundance_preds <- unique(unlist(lapply(abundance_importance_list, names)))
  
  # Create presence importance matrix
  presence_matrix <- matrix(NA, 
                            nrow = length(presence_importance_list),
                            ncol = length(all_presence_preds),
                            dimnames = list(names(presence_importance_list), all_presence_preds))
  
  for(sp in names(presence_importance_list)) {
    imp <- presence_importance_list[[sp]]
    if(length(imp) > 0) {
      presence_matrix[sp, names(imp)] <- imp
    }
  }
  
  presence_df <- as.data.frame(presence_matrix)
  presence_df$Species <- rownames(presence_df)
  presence_df <- presence_df[, c("Species", setdiff(names(presence_df), "Species"))]
  
  # Create abundance importance matrix
  abundance_matrix <- matrix(NA,
                             nrow = length(abundance_importance_list),
                             ncol = length(all_abundance_preds),
                             dimnames = list(names(abundance_importance_list), all_abundance_preds))
  
  for(sp in names(abundance_importance_list)) {
    imp <- abundance_importance_list[[sp]]
    if(length(imp) > 0) {
      abundance_matrix[sp, names(imp)] <- imp
    }
  }
  
  abundance_df <- as.data.frame(abundance_matrix)
  abundance_df$Species <- rownames(abundance_df)
  abundance_df <- abundance_df[, c("Species", setdiff(names(abundance_df), "Species"))]
  
  # Format and sort tables
  presence_result <- format_importance_table(presence_df)
  abundance_result <- format_importance_table(abundance_df)
  
  # Save formatted tables (1 decimal place)
  write.csv(presence_result$formatted, 
            file.path(output_dir, "presence_predictor_importance.csv"),
            row.names = FALSE)
  
  write.csv(abundance_result$formatted,
            file.path(output_dir, "abundance_predictor_importance.csv"),
            row.names = FALSE)
  
  # Save numeric values rounded to integers for Excel
  presence_numeric <- presence_result$importance
  abundance_numeric <- abundance_result$importance
  
  # Round to integers
  numeric_cols_pres <- names(presence_numeric)[names(presence_numeric) != "Species"]
  numeric_cols_abund <- names(abundance_numeric)[names(abundance_numeric) != "Species"]
  
  presence_numeric[numeric_cols_pres] <- lapply(presence_numeric[numeric_cols_pres], 
                                                function(x) round(x, 0))
  abundance_numeric[numeric_cols_abund] <- lapply(abundance_numeric[numeric_cols_abund], 
                                                  function(x) round(x, 0))
  
  write_xlsx(list(Presence = presence_numeric,
                  Abundance = abundance_numeric),
             file.path(output_dir, "predictor_importance_tables.xlsx"))
  
  cat("✓ Tables saved\n")
  
  ################################## Print Summary ##################################
  
  cat("\n========================================\n")
  cat("PRESENCE MODEL - TOP 3 PREDICTORS\n")
  cat("========================================\n\n")
  
  for(sp in rownames(presence_result$importance)) {
    vals <- presence_result$importance[sp, -1]
    vals <- vals[!is.na(vals)]
    if(length(vals) == 0) {
      cat(sprintf("%-6s: No predictors\n", sp))
      next
    }
    top3 <- sort(vals, decreasing = TRUE)[1:min(3, length(vals))]
    cat(sprintf("%-6s: %s\n", sp, 
                paste(paste0(names(top3), " (", sprintf("%.1f", top3), "%)"), 
                      collapse = ", ")))
  }
  
  cat("\n========================================\n")
  cat("ABUNDANCE MODEL - TOP 3 PREDICTORS\n")
  cat("========================================\n\n")
  
  for(sp in rownames(abundance_result$importance)) {
    vals <- abundance_result$importance[sp, -1]
    vals <- vals[!is.na(vals)]
    if(length(vals) == 0) {
      cat(sprintf("%-6s: No predictors\n", sp))
      next
    }
    top3 <- sort(vals, decreasing = TRUE)[1:min(3, length(vals))]
    cat(sprintf("%-6s: %s\n", sp,
                paste(paste0(names(top3), " (", sprintf("%.1f", top3), "%)"),
                      collapse = ", ")))
  }
  
  cat("\n✓ Analysis complete\n")
  cat("\nOutput saved to:", output_dir, "\n")
  cat("  - Individual species draw plots (PDF)\n")
  cat("  - Predictor importance tables (CSV and XLSX)\n")
  if(IMPORTANCE_METHOD == "gam.hp") {
    cat("\nVariable importance calculated using gam.hp hierarchical partitioning\n")
    cat("Citation: Lai et al. (2024) Plant Diversity 46(4):542-546\n")
    cat("DOI: 10.1016/j.pld.2024.06.002\n")
  } else if(IMPORTANCE_METHOD == "chisq") {
    cat("\nVariable importance calculated using relative chi-squared contributions\n")
    cat("Citation: Wood, S.N. (2017) Generalized Additive Models: An Introduction with R (2nd ed.)\n")
  }