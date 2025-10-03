library(scam)


################################## try curves for P/A ##################################

#STEP 1: Bin data to 100m grid

# Convert your data to a SpatVector with appropriate CRS
orbicella_vect <- vect(orbicella_model_data, 
                       geom = c("lon", "lat"), 
                       crs = "EPSG:4326")

# Project to a metric CRS (adjust UTM zone for your study area)
orbicella_vect <- project(orbicella_vect, "EPSG:32617")

# Create a 100m resolution raster template covering your data extent
ext_orbicella <- ext(orbicella_vect)
r_template <- rast(ext_orbicella, resolution = 100, crs = "EPSG:32617")

# Assign cell IDs to each point
cell_ids <- cells(r_template, orbicella_vect)[, "cell"]
orbicella_df <- as.data.frame(orbicella_vect, geom = "XY")
orbicella_df$cell_id <- cell_ids

# Aggregate by grid cell: ANY presence = 1, all absences = 0
# For abundance, take mean of cover_prop where present
orbicella_binned <- orbicella_df %>%
  group_by(cell_id) %>%
  summarise(
    present = as.integer(max(present, na.rm = TRUE)),
    cover_prop = mean(cover_prop[cover_prop > 0], na.rm = TRUE),  # Mean of presences only
    depth_bathy = mean(depth_bathy, na.rm = TRUE),
    aspect = mean(aspect, na.rm = TRUE),
    slope = mean(slope, na.rm = TRUE),
    complexity = mean(complexity, na.rm = TRUE),
    SAPA = mean(SAPA, na.rm = TRUE),
    dir_at_max_hsig = mean(dir_at_max_hsig, na.rm = TRUE),
    range_SST = mean(range_SST, na.rm = TRUE),
    mean_Hsig = mean(mean_Hsig, na.rm = TRUE),
    mean_SST = mean(mean_SST, na.rm = TRUE),
    mean_PAR = mean(mean_PAR, na.rm = TRUE),
    mean_chla = mean(mean_chla, na.rm = TRUE),
    mean_kd490 = mean(mean_kd490, na.rm = TRUE),
    dist_to_land = mean(dist_to_land, na.rm = TRUE),
    planform_curv = mean(planform_curv, na.rm = TRUE),
    lon = mean(x, na.rm = TRUE),
    lat = mean(y, na.rm = TRUE),
    .groups = "drop"
  )

# Check the aggregation
cat("Original data rows:", nrow(orbicella_model_data), "\n")
cat("Binned data rows:", nrow(orbicella_binned), "\n")

# 
# PRESENCE/ABSENCE MODEL WITH CALIBRATION
#

cat("\n\n===== PRESENCE/ABSENCE MODEL =====\n\n")

# Fit GAM with binned data
# orbicella_gam_binned <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
#                               s(slope) + s(complexity) +
#                               s(dir_at_max_hsig, bs = 'cc') +
#                               s(range_SST) + s(lon, lat),
#                             data = orbicella_binned,
#                             family = binomial())
orbicella_gam_binned <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                      s(planform_curv) + s(SAPA) +
                                      s(dir_at_max_hsig, bs = 'cc') +
                                      s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
                                      s(range_SST) +
                                      s(dist_to_land),
                                    data = orbicella_binned,
                                    # weights = weights_vec,
                                    # select = TRUE,
                                    family = binomial())

summary(orbicella_gam_binned)

# Extract the data actually used in the model fit
model_data_used_pa <- orbicella_gam_binned$model

# Get fitted values
fitted.response.pa <- predict(orbicella_gam_binned, type = "response")
fitted.link.pa <- predict(orbicella_gam_binned, type = "link")

# Create calibration plot
pdata.pa <- data.frame(fits = fitted.response.pa, 
                       observed = model_data_used_pa$present)

calibration_plot_pa <- ggplot(pdata.pa, aes(fits, jitter(observed, amount = 0.05))) +
  scale_x_continuous("GAM fitted values", limits = c(0, 1)) +
  scale_y_continuous("observed occurrences", limits = c(0, 1)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, col = "black") +
  geom_smooth(col = "grey60") +
  theme_bw() +
  ggtitle("P/A Calibration - Raw Predictions")

print(calibration_plot_pa)

# Calibration statistics
cali.gam.pa <- glm(model_data_used_pa$present ~ fitted.link.pa, 
                   family = binomial)

cat("\n=== P/A Calibration Statistics ===\n")
cat("Ideal: Intercept ≈ 0, Slope ≈ 1\n\n")
print(summary(cali.gam.pa))

# Fit calibration model
cal.scam.pa <- scam(model_data_used_pa$present ~ s(fitted.response.pa, bs = "mpi"), 
                    family = binomial)

# Get calibrated predictions
corr.fitted.response.pa <- predict(cal.scam.pa, type = "response")

# Plot calibrated predictions
pdata.corr.pa <- data.frame(corr = corr.fitted.response.pa,
                            observed = model_data_used_pa$present)

calibration_plot_corrected_pa <- ggplot(pdata.corr.pa, aes(corr, jitter(observed, amount = 0.05))) +
  scale_x_continuous("GAM: calibrated values", limits = c(0, 1)) +
  scale_y_continuous("observed occurrences", limits = c(0, 1)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, col = "black") +
  geom_smooth(col = "grey60") +
  theme_bw() +
  ggtitle("P/A Calibration - Calibrated Predictions")

print(calibration_plot_corrected_pa)












################################## abundance - curve types ##################################

#UNIFIED CALIBRATION CORRECTION FRAMEWORK
library(mgcv)
library(scam)
library(quantreg)

calibrate_gam_predictions <- function(
    fitted_model,
    method = c("isotonic", "spline", "quantile", "sigmoid"),
    shrinkage = 0.5,      # For sigmoid method: pull toward calibration (0-1)
    L_low = 0.8,          # For sigmoid method: correction below flexion point
    L_high = 1.2,         # For sigmoid method: correction above flexion point
    k = 50,               # For sigmoid method: transition steepness
    tau_seq = seq(0.1, 0.9, by = 0.1),  # For quantile method
    plot = TRUE
) {
  
  # Extract predictions and observations
  pred <- predict(fitted_model, type = "response")
  obs <- fitted_model$model[[1]]  # More robust than assuming name
  
  # Original R²
  r2_original <- cor(obs, pred)^2
  
  # Initialize corrected predictions
  pred_corrected <- pred
  
  #METHOD SELECTION
  method <- match.arg(method)
  
  if (method == "isotonic") {
    # Isotonic regression: monotonic, non-parametric calibration
    # This is often the most robust for calibration
    iso_model <- isoreg(pred, obs)
    
    # Function to apply correction
    correction_fn <- function(new_pred) {
      # Interpolate using isotonic regression
      corrected <- approx(iso_model$x, iso_model$yf, xout = new_pred, 
                          rule = 2, ties = mean)$y
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
  } else if (method == "spline") {
    # Smooth calibration curve with shape constraints
    # Use cyclic cubic regression spline for smooth calibration
    cal_data <- data.frame(pred = pred, obs = obs)
    cal_model <- gam(obs ~ s(pred, bs = "cr", k = 10), 
                     data = cal_data, method = "REML")
    
    correction_fn <- function(new_pred) {
      corrected <- predict(cal_model, newdata = data.frame(pred = new_pred), 
                           type = "response")
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
  } else if (method == "quantile") {
    # Improved quantile approach: use conditional median as primary calibration
    # Weight toward median but allow for heteroscedastic correction
    
    # Fit quantiles
    qr_model <- rq(obs ~ pred, tau = tau_seq)
    
    # Get median (50th percentile) if available, otherwise central quantile
    median_tau <- which.min(abs(tau_seq - 0.5))
    qr_median_model <- rq(obs ~ pred, tau = tau_seq[median_tau])
    
    correction_fn <- function(new_pred) {
      # Primary correction from median
      corrected <- predict(qr_median_model, newdata = data.frame(pred = new_pred))
      
      # Optional: blend with mean of quantiles for stability
      qr_all <- predict(qr_model, newdata = data.frame(pred = new_pred))
      qr_mean <- rowMeans(qr_all, na.rm = TRUE)
      
      # Weighted average favoring median
      corrected <- 0.7 * corrected + 0.3 * qr_mean
      
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
  } else if (method == "sigmoid") {
    # Your flexible sigmoid approach - improved parameterization
    
    # Fit calibration line
    bestfitline_model <- lm(obs ~ pred)
    a <- coef(bestfitline_model)[1]
    b <- coef(bestfitline_model)[2]
    
    # Flexion point (where calibration crosses identity)
    x0 <- a / (1 - b)
    
    correction_fn <- function(new_pred) {
      # Target from calibration line
      pred_target <- a + b * new_pred
      
      # Shrink toward target
      pred_shrunk <- new_pred + shrinkage * (pred_target - new_pred)
      
      # Multiplicative correction
      sigmoid_val <- 1 / (1 + exp(-k * (pred_shrunk - x0)))
      correction_factor <- L_low + (L_high - L_low) * sigmoid_val
      
      corrected <- pred_shrunk * correction_factor
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
    cat("Sigmoid parameters:\n")
    cat("  Flexion point (x0):", round(x0, 4), "\n")
    cat("  Shrinkage:", shrinkage, "\n")
    cat("  L_low:", L_low, "| L_high:", L_high, "\n\n")
  }
  
  # Calculate corrected R²
  r2_corrected <- cor(obs, pred_corrected)^2
  
  # Calculate calibration metrics
  slope_original <- coef(lm(obs ~ pred))[2]
  slope_corrected <- coef(lm(obs ~ pred_corrected))[2]
  
  #plotting
  if (plot) {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    lim <- c(0, max(c(obs, pred, pred_corrected)))
    
    # Before
    plot(pred, obs, xlim = lim, ylim = lim, pch = 16, col = rgb(0, 0, 0, 0.3),
         xlab = "Predicted", ylab = "Observed", main = "Original")
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(obs ~ pred), col = "black", lwd = 1, lty = 2)
    text(lim[1], lim[2] * 0.95, 
         paste0("R² = ", round(r2_original, 3), "\nSlope = ", round(slope_original, 3)),
         adj = 0, cex = 0.9)
    
    # After
    plot(pred_corrected, obs, xlim = lim, ylim = lim, pch = 16, 
         col = rgb(0, 0.4, 0.8, 0.3),
         xlab = "Predicted (Corrected)", ylab = "Observed",
         main = paste("Calibrated:", method))
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(obs ~ pred_corrected), col = "darkblue", lwd = 1, lty = 2)
    text(lim[1], lim[2] * 0.95,
         paste0("R² = ", round(r2_corrected, 3), "\nSlope = ", round(slope_corrected, 3)),
         adj = 0, cex = 0.9)
    
    par(mfrow = c(1, 1))
  }
  
  #results
  cat("\n=== CALIBRATION SUMMARY ===\n")
  cat("Method:", method, "\n")
  cat("Original  - R²:", round(r2_original, 3), "| Slope:", round(slope_original, 3), "\n")
  cat("Corrected - R²:", round(r2_corrected, 3), "| Slope:", round(slope_corrected, 3), "\n\n")
  
  # Return correction function and diagnostics
  list(
    method = method,
    correction_function = correction_fn,
    r2_original = r2_original,
    r2_corrected = r2_corrected,
    slope_original = slope_original,
    slope_corrected = slope_corrected,
    pred_corrected = pred_corrected,
    observations = obs,
    predictions_original = pred
  )
}


################################## optimize curve types ##################################
#OPTIMIZED CALIBRATION WITH BIAS CORRECTION
library(mgcv)
library(scam)
library(quantreg)

calibrate_gam_predictions <- function(
    fitted_model,
    method = c("isotonic", "spline", "quantile", "sigmoid", "sigmoid_optimized"),
    optimize_params = FALSE,  # NEW: auto-optimize parameters
    shrinkage = 0.5,
    L_low = 0.8,
    L_high = 1.2,
    k = 50,
    tau_seq = seq(0.1, 0.9, by = 0.1),
    plot = TRUE
) {
  
  # Extract predictions and observations
  pred <- predict(fitted_model, type = "response")
  obs <- fitted_model$model[[1]]
  
  # Original R²
  r2_original <- cor(obs, pred)^2
  
  # Initialize
  pred_corrected <- pred
  method <- match.arg(method)
  
  #OBJECTIVE FUNCTION FOR OPTIMIZATION
  calc_calibration_metrics <- function(pred_cal, obs_cal) {
    # Multiple metrics to assess calibration quality
    
    # 1. Calibration slope (want close to 1)
    slope <- coef(lm(obs_cal ~ pred_cal))[2]
    slope_penalty <- (slope - 1)^2
    
    # 2. Calibration intercept (want close to 0)
    intercept <- coef(lm(obs_cal ~ pred_cal))[1]
    intercept_penalty <- intercept^2
    
    # 3. Binned calibration error (check bias across range)
    n_bins <- 10
    pred_bins <- cut(pred_cal, breaks = quantile(pred_cal, probs = seq(0, 1, length.out = n_bins + 1)),
                     include.lowest = TRUE, labels = FALSE)
    
    bin_bias <- sapply(1:n_bins, function(i) {
      if (sum(pred_bins == i, na.rm = TRUE) > 0) {
        mean(obs_cal[pred_bins == i], na.rm = TRUE) - mean(pred_cal[pred_bins == i], na.rm = TRUE)
      } else {
        0
      }
    })
    
    # Mean absolute binned bias (want close to 0)
    mab <- mean(abs(bin_bias))
    
    # 4. R² (want high, but less important than calibration)
    r2 <- cor(obs_cal, pred_cal)^2
    
    # Combined loss (weighted to prioritize calibration over R²)
    loss <- 10 * slope_penalty + 5 * intercept_penalty + 20 * mab - 0.5 * r2
    
    return(list(
      loss = loss,
      slope = slope,
      intercept = intercept,
      mab = mab,
      r2 = r2,
      bin_bias = bin_bias
    ))
  }
  
  #SIGMOID WITH OPTIMIZATION
  if (method == "sigmoid_optimized" || (method == "sigmoid" && optimize_params)) {
    
    cat("Optimizing sigmoid parameters...\n")
    
    # Fit calibration line
    bestfitline_model <- lm(obs ~ pred)
    a <- coef(bestfitline_model)[1]
    b <- coef(bestfitline_model)[2]
    x0_init <- a / (1 - b)
    
    # Objective function for optimization
    sigmoid_objective <- function(params) {
      shrinkage_opt <- plogis(params[1])  # Constrain to [0, 1]
      L_low_opt <- exp(params[2])         # Constrain to positive
      L_high_opt <- exp(params[3])        # Constrain to positive
      k_opt <- exp(params[4])             # Constrain to positive
      x0_opt <- plogis(params[5]) * max(pred)  # Constrain to data range
      
      # Apply correction
      pred_target <- a + b * pred
      pred_shrunk <- pred + shrinkage_opt * (pred_target - pred)
      sigmoid_val <- 1 / (1 + exp(-k_opt * (pred_shrunk - x0_opt)))
      correction_factor <- L_low_opt + (L_high_opt - L_low_opt) * sigmoid_val
      pred_test <- pred_shrunk * correction_factor
      pred_test <- pmax(0.001, pmin(0.999, pred_test))
      
      # Calculate loss
      metrics <- calc_calibration_metrics(pred_test, obs)
      return(metrics$loss)
    }
    
    # Initial parameters (transformed)
    init_params <- c(
      qlogis(shrinkage),
      log(L_low),
      log(L_high),
      log(k),
      qlogis(x0_init / max(pred))
    )
    
    # Optimize
    opt_result <- optim(
      par = init_params,
      fn = sigmoid_objective,
      method = "Nelder-Mead",
      control = list(maxit = 1000)
    )
    
    # Extract optimized parameters
    shrinkage <- plogis(opt_result$par[1])
    L_low <- exp(opt_result$par[2])
    L_high <- exp(opt_result$par[3])
    k <- exp(opt_result$par[4])
    x0 <- plogis(opt_result$par[5]) * max(pred)
    
    cat("Optimized parameters:\n")
    cat("  Flexion point (x0):", round(x0, 4), "\n")
    cat("  Shrinkage:", round(shrinkage, 4), "\n")
    cat("  L_low:", round(L_low, 4), "| L_high:", round(L_high, 4), "\n")
    cat("  k:", round(k, 2), "\n\n")
    
    # Apply with optimized parameters
    correction_fn <- function(new_pred) {
      pred_target <- a + b * new_pred
      pred_shrunk <- new_pred + shrinkage * (pred_target - new_pred)
      sigmoid_val <- 1 / (1 + exp(-k * (pred_shrunk - x0)))
      correction_factor <- L_low + (L_high - L_low) * sigmoid_val
      corrected <- pred_shrunk * correction_factor
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
  } else if (method == "sigmoid") {
    # Original sigmoid without optimization
    bestfitline_model <- lm(obs ~ pred)
    a <- coef(bestfitline_model)[1]
    b <- coef(bestfitline_model)[2]
    x0 <- a / (1 - b)
    
    correction_fn <- function(new_pred) {
      pred_target <- a + b * new_pred
      pred_shrunk <- new_pred + shrinkage * (pred_target - new_pred)
      sigmoid_val <- 1 / (1 + exp(-k * (pred_shrunk - x0)))
      correction_factor <- L_low + (L_high - L_low) * sigmoid_val
      corrected <- pred_shrunk * correction_factor
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
  } else if (method == "isotonic") {
    iso_model <- isoreg(pred, obs)
    
    correction_fn <- function(new_pred) {
      corrected <- approx(iso_model$x, iso_model$yf, xout = new_pred, 
                          rule = 2, ties = mean)$y
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
  } else if (method == "spline") {
    cal_data <- data.frame(pred = pred, obs = obs)
    cal_model <- gam(obs ~ s(pred, bs = "cr", k = 10), 
                     data = cal_data, method = "REML")
    
    correction_fn <- function(new_pred) {
      corrected <- predict(cal_model, newdata = data.frame(pred = new_pred), 
                           type = "response")
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
  } else if (method == "quantile") {
    qr_model <- rq(obs ~ pred, tau = tau_seq)
    median_tau <- which.min(abs(tau_seq - 0.5))
    qr_median_model <- rq(obs ~ pred, tau = tau_seq[median_tau])
    
    correction_fn <- function(new_pred) {
      corrected <- predict(qr_median_model, newdata = data.frame(pred = new_pred))
      qr_all <- predict(qr_model, newdata = data.frame(pred = new_pred))
      qr_mean <- rowMeans(qr_all, na.rm = TRUE)
      corrected <- 0.7 * corrected + 0.3 * qr_mean
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
  }
  
  # Calculate metrics
  r2_corrected <- cor(obs, pred_corrected)^2
  slope_original <- coef(lm(obs ~ pred))[2]
  slope_corrected <- coef(lm(obs ~ pred_corrected))[2]
  intercept_corrected <- coef(lm(obs ~ pred_corrected))[1]
  
  # Binned bias analysis
  metrics_corrected <- calc_calibration_metrics(pred_corrected, obs)
  
  #ENHANCED PLOTTING WITH BIAS DIAGNOSTICS
  if (plot) {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    lim <- c(0, max(c(obs, pred, pred_corrected)))
    
    # Plot 1: Original
    plot(pred, obs, xlim = lim, ylim = lim, pch = 16, col = rgb(0, 0, 0, 0.3),
         xlab = "Predicted", ylab = "Observed", main = "Original")
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(obs ~ pred), col = "black", lwd = 1, lty = 2)
    text(lim[1], lim[2] * 0.95, 
         paste0("R² = ", round(r2_original, 3), "\nSlope = ", round(slope_original, 3)),
         adj = 0, cex = 0.9)
    
    # Plot 2: Corrected
    plot(pred_corrected, obs, xlim = lim, ylim = lim, pch = 16, 
         col = rgb(0, 0.4, 0.8, 0.3),
         xlab = "Predicted (Corrected)", ylab = "Observed",
         main = paste("Calibrated:", method))
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(obs ~ pred_corrected), col = "darkblue", lwd = 1, lty = 2)
    text(lim[1], lim[2] * 0.95,
         paste0("R² = ", round(r2_corrected, 3), 
                "\nSlope = ", round(slope_corrected, 3),
                "\nIntercept = ", round(intercept_corrected, 3)),
         adj = 0, cex = 0.8)
    
    # Plot 3: Binned bias - Original
    n_bins <- 10
    pred_bins_orig <- cut(pred, breaks = quantile(pred, probs = seq(0, 1, length.out = n_bins + 1)),
                          include.lowest = TRUE, labels = FALSE)
    bin_centers_orig <- sapply(1:n_bins, function(i) mean(pred[pred_bins_orig == i], na.rm = TRUE))
    bin_bias_orig <- sapply(1:n_bins, function(i) {
      mean(obs[pred_bins_orig == i], na.rm = TRUE) - mean(pred[pred_bins_orig == i], na.rm = TRUE)
    })
    
    plot(bin_centers_orig, bin_bias_orig, type = "b", pch = 19, col = "black",
         xlab = "Predicted (binned)", ylab = "Bias (Obs - Pred)",
         main = "Original: Bias by Prediction Range", ylim = range(c(bin_bias_orig, 0)) * 1.2)
    abline(h = 0, col = "red", lwd = 2, lty = 2)
    grid()
    
    # Plot 4: Binned bias - Corrected
    pred_bins_corr <- cut(pred_corrected, breaks = quantile(pred_corrected, probs = seq(0, 1, length.out = n_bins + 1)),
                          include.lowest = TRUE, labels = FALSE)
    bin_centers_corr <- sapply(1:n_bins, function(i) mean(pred_corrected[pred_bins_corr == i], na.rm = TRUE))
    bin_bias_corr <- metrics_corrected$bin_bias
    
    plot(bin_centers_corr, bin_bias_corr, type = "b", pch = 19, col = "darkblue",
         xlab = "Predicted (binned)", ylab = "Bias (Obs - Pred)",
         main = "Corrected: Bias by Prediction Range", ylim = range(c(bin_bias_corr, 0)) * 1.2)
    abline(h = 0, col = "red", lwd = 2, lty = 2)
    grid()
    text(min(bin_centers_corr), max(bin_bias_corr) * 0.9,
         paste0("Mean Abs Bias = ", round(metrics_corrected$mab, 4)),
         adj = 0, cex = 0.9)
    
    par(mfrow = c(1, 1))
  }
  
  #RESULTS
  cat("\n=== CALIBRATION SUMMARY ===\n")
  cat("Method:", method, "\n")
  cat("Original  - R²:", round(r2_original, 3), "| Slope:", round(slope_original, 3), "\n")
  cat("Corrected - R²:", round(r2_corrected, 3), "| Slope:", round(slope_corrected, 3),
      "| Intercept:", round(intercept_corrected, 3), "\n")
  cat("Mean Absolute Bias (binned):", round(metrics_corrected$mab, 4), "\n\n")
  
  # Return
  list(
    method = method,
    correction_function = correction_fn,
    r2_original = r2_original,
    r2_corrected = r2_corrected,
    slope_original = slope_original,
    slope_corrected = slope_corrected,
    intercept_corrected = intercept_corrected,
    mean_abs_bias = metrics_corrected$mab,
    pred_corrected = pred_corrected,
    observations = obs,
    predictions_original = pred
  )
}

#USAGE

# Automatically optimize sigmoid parameters
result_optimized <- calibrate_gam_predictions(
  orbicella_gam_abundance_beta,
  method = "sigmoid_optimized",
  plot = TRUE
)

# Or optimize with custom starting values
result_custom <- calibrate_gam_predictions(
  orbicella_gam_abundance_beta,
  method = "sigmoid",
  optimize_params = TRUE,
  shrinkage = 0.7,
  L_low = 0.6,
  L_high = 1.5,
  k = 100,
  plot = TRUE
)

# Apply to new data
# new_pred <- predict(orbicella_gam_abundance_beta, newdata = new_data, type = "response")
# new_pred_corrected <- result_optimized$correction_function(new_pred)

#usage examples

# Try isotonic regression (usually best for calibration)
result_iso <- calibrate_gam_predictions(
  orbicella_gam_abundance_beta,
  method = "isotonic",
  plot = TRUE
)

# Try sigmoid approach with your parameters
result_sig <- calibrate_gam_predictions(
  orbicella_gam_abundance_beta,
  method = "quantile",
  shrinkage = 0.5,
  L_low = 0.8,
  L_high = 1.2,
  k = 50,
  plot = TRUE
)

#APPLY TO NEW PREDICTIONS
# For out-of-sample predictions:
# new_pred <- predict(orbicella_gam_abundance_beta, newdata = new_data, type = "response")
# new_pred_corrected <- result_iso$correction_function(new_pred)






################################## Waldock 2022 rescaling ##################################

# Simple linear rescaling calibration (from the 2021 abundance SDM review)
calibrate_rescale <- function(fitted_model, plot = TRUE) {
  
  pred <- predict(fitted_model, type = "response")
  obs <- fitted_model$model[[1]]
  
  # Linear rescaling formula from the paper
  min_obs <- min(obs)
  max_obs <- max(obs)
  min_pred <- min(pred)
  max_pred <- max(pred)
  
  pred_corrected <- (pred - min_pred) / (max_pred - min_pred) * (max_obs - min_obs) + min_obs
  
  # Metrics
  r2_orig <- cor(obs, pred)^2
  r2_corr <- cor(obs, pred_corrected)^2
  slope_orig <- coef(lm(obs ~ pred))[2]
  slope_corr <- coef(lm(obs ~ pred_corrected))[2]
  mae_orig <- mean(abs(obs - pred))
  mae_corr <- mean(abs(obs - pred_corrected))
  
  if (plot) {
    # Set common axis limits
    lim <- c(0, max(c(obs, pred, pred_corrected)))
    
    par(mfrow = c(1, 2))
    
    # Original
    plot(pred, obs, pch = 16, col = rgb(0.3, 0.3, 0.3, 0.3),
         xlim = lim, ylim = lim,
         main = "Original", xlab = "Predicted", ylab = "Observed")
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(obs ~ pred), col = "black", lwd = 2, lty = 2)
    text(lim[1] + 0.01, lim[2] * 0.98, 
         sprintf("R²=%.3f\nSlope=%.3f", r2_orig, slope_orig),
         adj = c(0, 1))
    
    # Rescaled
    plot(pred_corrected, obs, pch = 16, col = rgb(0, 0.4, 0.7, 0.3),
         xlim = lim, ylim = lim,
         main = "Rescaled", xlab = "Predicted (Rescaled)", ylab = "Observed")
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(obs ~ pred_corrected), col = "darkblue", lwd = 2, lty = 2)
    text(lim[1] + 0.01, lim[2] * 0.98,
         sprintf("R²=%.3f\nSlope=%.3f", r2_corr, slope_corr),
         adj = c(0, 1), col = "darkblue")
    
    par(mfrow = c(1, 1))
  }
  
  cat(sprintf("\nORIGINAL: R²=%.3f  Slope=%.3f  MAE=%.4f\n", 
              r2_orig, slope_orig, mae_orig))
  cat(sprintf("RESCALED: R²=%.3f  Slope=%.3f  MAE=%.4f\n\n", 
              r2_corr, slope_corr, mae_corr))
  
  list(
    pred_corrected = pred_corrected,
    min_obs = min_obs,
    max_obs = max_obs,
    min_pred = min_pred,
    max_pred = max_pred
  )
}

# Usage
result_rescale <- calibrate_rescale(orbicella_gam_abundance_beta, plot = TRUE)

# Apply to new data
calibrate_new <- function(new_pred, calib_result) {
  (new_pred - calib_result$min_pred) / 
    (calib_result$max_pred - calib_result$min_pred) * 
    (calib_result$max_obs - calib_result$min_obs) + 
    calib_result$min_obs
}
################################## Dan / gam/lm method ##################################

library(mgcv)
library(ggplot2)
library(gridExtra)

#SIMPLE REGRESSION CALIBRATION
calibrate_regression <- function(fitted_model, 
                                 method = c("lm", "gam", "loess"),
                                 intercept = FALSE,
                                 plot = TRUE) {
  
  method <- match.arg(method)
  
  pred <- predict(fitted_model, type = "response")
  obs <- fitted_model$model[[1]]
  
  # Fit calibration model: obs ~ f(pred)
  if (method == "lm") {
    if (intercept) {
      calib_model <- lm(obs ~ pred)
    } else {
      calib_model <- lm(obs ~ pred + 0)
    }
  } else if (method == "gam") {
    if (intercept) {
      calib_model <- gam(obs ~ s(pred, k = 10))
    } else {
      calib_model <- gam(obs ~ s(pred, k = 10) + 0)
    }
  } else if (method == "loess") {
    calib_model <- loess(obs ~ pred, span = 0.3)
  }
  
  # Create correction function
  correction_fn <- function(new_pred) {
    predict(calib_model, newdata = data.frame(pred = new_pred))
  }
  
  pred_corrected <- correction_fn(pred)
  
  # Metrics
  r2_orig <- cor(obs, pred)^2
  r2_corr <- cor(obs, pred_corrected)^2
  slope_orig <- coef(lm(obs ~ pred))[2]
  slope_corr <- coef(lm(obs ~ pred_corrected))[2]
  
  # Plot
  if (plot) {
    df_orig <- data.frame(pred = pred, obs = obs)
    df_corr <- data.frame(pred = pred_corrected, obs = obs)
    df_resid <- data.frame(pred = rep(pred, 2),
                           residual = c(obs - pred, obs - pred_corrected),
                           type = rep(c("Original", "Calibrated"), each = length(pred)))
    
    lim <- c(0, max(c(obs, pred_corrected)))
    
    # Plot 1: Original
    p1 <- ggplot(df_orig, aes(x = pred, y = obs)) +
      geom_point(alpha = 0.3, size = 2) +
      geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "black", 
                  linetype = "dashed", linewidth = 1) +
      coord_fixed(xlim = lim, ylim = lim) +
      labs(title = "Original", x = "Predicted", y = "Observed") +
      annotate("text", x = lim[1], y = lim[2] * 0.98, 
               label = sprintf("R²=%.3f\nSlope=%.3f", r2_orig, slope_orig),
               hjust = 0, vjust = 1, size = 3.5) +
      theme_bw(base_size = 11) +
      theme(aspect.ratio = 1)
    
    # Plot 2: Corrected
    p2 <- ggplot(df_corr, aes(x = pred, y = obs)) +
      geom_point(alpha = 0.3, size = 2, color = "#0066B3") +
      geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "darkblue", 
                  linetype = "dashed", linewidth = 1) +
      coord_fixed(xlim = lim, ylim = lim) +
      labs(title = paste0("Calibrated (", method, ")"), 
           x = "Predicted (Corrected)", y = "Observed") +
      annotate("text", x = lim[1], y = lim[2] * 0.98,
               label = sprintf("R²=%.3f\nSlope=%.3f", r2_corr, slope_corr),
               hjust = 0, vjust = 1, size = 3.5, color = "darkblue") +
      theme_bw(base_size = 11) +
      theme(aspect.ratio = 1)
    
    # Plot 3: Correction curve
    pred_seq <- seq(min(pred), max(pred), length.out = 200)
    df_curve <- data.frame(x = pred_seq, y = correction_fn(pred_seq))
    
    p3 <- ggplot(df_curve, aes(x = x, y = y)) +
      geom_line(color = "darkblue", linewidth = 1.5) +
      geom_abline(slope = 1, intercept = 0, color = "red", 
                  linetype = "dashed", linewidth = 1) +
      labs(title = "Correction Function", 
           x = "Original Prediction", y = "Corrected Prediction") +
      theme_bw(base_size = 11)
    
    # Plot 4: Residuals
    p4 <- ggplot(df_resid, aes(x = pred, y = residual, color = type)) +
      geom_point(alpha = 0.3, size = 1.5) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
      scale_color_manual(values = c("Original" = "gray40", "Calibrated" = "#0066B3")) +
      labs(title = "Residuals", x = "Prediction", y = "Residual", color = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position = c(0.85, 0.9),
            legend.background = element_rect(fill = "white", color = NA))
    
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  }
  
  cat(sprintf("\nORIGINAL:   R²=%.3f  Slope=%.3f\n", r2_orig, slope_orig))
  cat(sprintf("CALIBRATED: R²=%.3f  Slope=%.3f\n\n", r2_corr, slope_corr))
  
  list(
    correction_function = correction_fn,
    calibration_model = calib_model,
    pred_corrected = pred_corrected
  )
}

#USAGE

# GAM calibration - NO intercept (forced through origin)
result_gam <- calibrate_regression(orbicella_gam_abundance_beta, 
                                   method = "gam", 
                                   intercept = FALSE,
                                   plot = TRUE)

# GAM calibration - WITH intercept
result_gam_int <- calibrate_regression(orbicella_gam_abundance_beta, 
                                       method = "gam", 
                                       intercept = TRUE,
                                       plot = TRUE)

# Linear calibration - NO intercept
result_lm <- calibrate_regression(orbicella_gam_abundance_beta, 
                                  method = "lm",
                                  intercept = FALSE,
                                  plot = TRUE)

# Linear calibration - WITH intercept
result_lm <- calibrate_regression(orbicella_gam_abundance_beta, 
                                  method = "lm",
                                  intercept = TRUE,
                                  plot = TRUE)

# Linear calibration - NO intercept
result_loess <- calibrate_regression(orbicella_gam_abundance_beta, 
                                     method = "loess",
                                     plot = TRUE)

# Apply to new data
# new_calibrated <- result_gam$correction_function(new_predictions)

################################## quantile method ##################################
library(mgcv)

#SIMPLE QUANTILE MAPPING CALIBRATION
calibrate_quantiles <- function(fitted_model, n_quantiles = 100, plot = TRUE) {
  
  pred <- predict(fitted_model, type = "response")
  obs <- fitted_model$model[[1]]
  
  # Map prediction quantiles to observation quantiles
  probs <- seq(0, 1, length.out = n_quantiles)
  pred_q <- quantile(pred, probs)
  obs_q <- quantile(obs, probs)
  
  # Create correction function
  correction_fn <- approxfun(pred_q, obs_q, method = "linear", rule = 2)
  pred_corrected <- correction_fn(pred)
  
  # Metrics
  r2_orig <- cor(obs, pred)^2
  r2_corr <- cor(obs, pred_corrected)^2
  slope_orig <- coef(lm(obs ~ pred))[2]
  slope_corr <- coef(lm(obs ~ pred_corrected))[2]
  mae_orig <- mean(abs(obs - pred))
  mae_corr <- mean(abs(obs - pred_corrected))
  
  # Plot
  if (plot) {
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
    lim <- c(0, max(c(obs, pred_corrected)))
    
    # Original
    plot(pred, obs, xlim = lim, ylim = lim, pch = 16, 
         col = rgb(0.3, 0.3, 0.3, 0.3), main = "Original",
         xlab = "Predicted", ylab = "Observed")
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(obs ~ pred), col = "black", lwd = 2, lty = 2)
    text(0.01, lim[2] * 0.98, 
         sprintf("R²=%.3f\nSlope=%.3f\nMAE=%.4f", r2_orig, slope_orig, mae_orig),
         adj = c(0, 1), cex = 0.8)
    
    # Corrected
    plot(pred_corrected, obs, xlim = lim, ylim = lim, pch = 16,
         col = rgb(0, 0.4, 0.7, 0.3), main = "Calibrated",
         xlab = "Predicted (Calibrated)", ylab = "Observed")
    abline(0, 1, col = "red", lwd = 2)
    abline(lm(obs ~ pred_corrected), col = "darkblue", lwd = 2, lty = 2)
    text(0.01, lim[2] * 0.98,
         sprintf("R²=%.3f\nSlope=%.3f\nMAE=%.4f", r2_corr, slope_corr, mae_corr),
         adj = c(0, 1), cex = 0.8, col = "darkblue")
    
    # Correction curve
    pred_seq <- seq(0, max(pred) * 1.1, length.out = 500)
    plot(pred_seq, correction_fn(pred_seq), type = "l", lwd = 3, 
         col = "darkblue", main = "Correction Function",
         xlab = "Original Prediction", ylab = "Corrected Prediction")
    abline(0, 1, col = "red", lwd = 2, lty = 2)
    points(pred_q, obs_q, pch = 19, cex = 0.5, col = rgb(1, 0.5, 0, 0.6))
    grid(col = "gray85")
    
    # Residuals
    plot(pred, obs - pred, pch = 16, col = rgb(0.3, 0.3, 0.3, 0.3),
         main = "Residuals", xlab = "Prediction", ylab = "Obs - Pred",
         ylim = range(c(obs - pred, obs - pred_corrected)))
    points(pred, obs - pred_corrected, pch = 16, col = rgb(0, 0.4, 0.7, 0.3))
    abline(h = 0, col = "red", lwd = 2, lty = 2)
    legend("topright", c("Original", "Calibrated"),
           col = c(rgb(0.3, 0.3, 0.3, 0.6), rgb(0, 0.4, 0.7, 0.6)),
           pch = 16, bty = "n")
    
    par(mfrow = c(1, 1))
  }
  
  # Summary
  cat(sprintf("\nORIGINAL:   R²=%.3f  Slope=%.3f  MAE=%.4f\n", 
              r2_orig, slope_orig, mae_orig))
  cat(sprintf("CALIBRATED: R²=%.3f  Slope=%.3f  MAE=%.4f\n\n", 
              r2_corr, slope_corr, mae_corr))
  
  # Create comparison dataframe
  comparison_df <- data.frame(
    observed = obs,
    predicted_original = pred,
    predicted_corrected = pred_corrected,
    residual_original = obs - pred,
    residual_corrected = obs - pred_corrected
  )
  
  list(
    correction_function = correction_fn,
    pred_corrected = pred_corrected,
    observations = obs,
    predictions_original = pred,
    comparison_df = comparison_df
  )
}

#USAGE
# result <- calibrate_quantiles(orbicella_gam_abundance_beta, plot = TRUE)
result <- calibrate_quantiles(orbicella_gam_abundance_beta, 
                              n_quantiles = 100, 
                              plot = TRUE)

# View the comparison dataframe
test = result$comparison_df


# Create quantile comparison plot
pred <- result$predictions_original
obs <- result$observations
n_quantiles = 100

probs <- seq(0, 1, length.out = 100)
pred_q <- quantile(pred, probs)
obs_q <- quantile(obs, probs)

# Plot the quantiles
par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))

# Plot 1: Quantile values
plot(probs * n_quantiles, pred_q, type = "l", lwd = 2.5, col = "gray40",
     xlab = "Percentile", ylab = "Value",
     main = "Quantile Distributions", ylim = c(0, max(obs_q)))
lines(probs * n_quantiles, obs_q, lwd = 2.5, col = "red")
legend("topleft", legend = c("Predicted", "Observed"),
       col = c("gray40", "red"), lwd = 2.5, bty = "n")
grid()

# Plot 2: Direct comparison (the mapping pairs)
plot(pred_q, obs_q, type = "p", pch = 19, cex = 0.8, col = "darkblue",
     xlab = "Predicted Quantiles", ylab = "Observed Quantiles",
     main = "Quantile Pairs (What Gets Mapped)")
abline(0, 1, col = "red", lwd = 2, lty = 2)
grid()

par(mfrow = c(1, 1))
# Apply to new data:
# new_calibrated <- result$correction_function(new_predictions)