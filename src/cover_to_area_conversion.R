  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(broom)
  library(plotly)
  
  ################################## setup ##################################
  
  #import workspace from upstream script
  load(here("output/process_species_data.RData"))
  
  ################################## main ##################################
  
  # Recreate the shared PSUs vector
  demo_psus <- combined_demo_data_averaged_psu %>% 
    distinct(PSU) %>% 
    pull(PSU)
  
  benthic_psus <- combined_benthic_data_averaged_psu %>% 
    distinct(PSU) %>% 
    pull(PSU)
  
  # Find PSUs present in both datasets
  shared_psus <- intersect(demo_psus, benthic_psus)
  
  # Display results for verification
  cat("PSUs in demo data:", length(demo_psus), "\n")
  cat("PSUs in benthic data:", length(benthic_psus), "\n")
  cat("PSUs in both datasets:", length(shared_psus), "\n")
  
  
  # More comprehensive version with additional metadata
  shared_psu_data <- data.frame(PSU = shared_psus) %>%
    # Join with demo data
    left_join(
      combined_demo_data_averaged_psu %>% 
        select(PSU, dataset, lat, lon, depth, SA_density), 
      by = "PSU"
    ) %>%
    # Join with benthic data
    left_join(
      combined_benthic_data_averaged_psu %>% 
        select(PSU, dataset, cover), 
      by = "PSU",
      suffix = c("_demo", "_benthic")
    ) %>%
    # Clean up and organize columns
    select(PSU, 
           demo_dataset = dataset_demo,
           benthic_dataset = dataset_benthic,
           lat, lon, depth,
           SA_density, 
           cover) %>%
    arrange(PSU)
  
  print(shared_psu_data)
  
  # Summary statistics
  cat("Summary of shared PSU data:\n")
  cat("Number of shared PSUs:", nrow(shared_psu_data), "\n")
  cat("SA_density range:", round(range(shared_psu_data$SA_density, na.rm = TRUE), 4), "\n")
  cat("Cover range:", round(range(shared_psu_data$cover, na.rm = TRUE), 2), "\n")
  
  
  # Create the scatter plot
  ggplot(shared_psu_data, aes(x = SA_density, y = cover)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3) +
    labs(
      title = "Relationship between SA Density and Cover",
      subtitle = paste("Based on", nrow(shared_psu_data), "shared PSUs"),
      x = "SA Density",
      y = "Cover (%)",
      caption = "Blue line shows linear trend with 95% confidence interval"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  
  
  # Calculate correlation
  correlation <- cor(shared_psu_data$SA_density, shared_psu_data$cover, use = "complete.obs")
  
  # Enhanced plot with correlation info
  ggplot(shared_psu_data, aes(x = SA_density, y = cover)) +
    geom_point(alpha = 0.7, size = 2, color = "darkblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
    labs(
      title = "Relationship between SA Density and Cover",
      subtitle = paste("n =", nrow(shared_psu_data), "PSUs | r =", round(correlation, 3)),
      x = "SA Density",
      y = "Cover (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  
  # Print correlation info
  cat("Correlation between SA_density and cover:", round(correlation, 4), "\n")
  
  
  
  
  
  # Create a linear model
  sa_cover_model <- lm(SA_density ~ cover, data = shared_psu_data)
  
  # Check model summary
  summary(sa_cover_model)
  
  # Function to predict SA_density from cover
  predict_sa_density <- function(cover_values) {
    predict(sa_cover_model, newdata = data.frame(cover = cover_values))
  }
  
  # Example usage
  predicted_sa <- predict_sa_density(c(10, 20, 30))
  print(predicted_sa)
  
  
  # Calculate mean ratio
  mean_ratio <- mean(shared_psu_data$SA_density / shared_psu_data$cover, na.rm = TRUE)
  
  # Function using ratio
  predict_sa_density_ratio <- function(cover_values) {
    cover_values * mean_ratio
  }
  
  cat("Mean SA_density/Cover ratio:", round(mean_ratio, 4), "\n")
  
  
  
  
  
  # Add predictions to your shared data for comparison
  shared_psu_data <- shared_psu_data %>%
    mutate(
      predicted_sa_lm = predict(sa_cover_model),
      predicted_sa_ratio = cover * mean_ratio,
      residual_lm = SA_density - predicted_sa_lm
    )
  
  # Check model performance
  cat("Linear Model R-squared:", round(summary(sa_cover_model)$r.squared, 3), "\n")
  cat("RMSE:", round(sqrt(mean(shared_psu_data$residual_lm^2, na.rm = TRUE)), 4), "\n")
  
  # Plot to compare methods
  ggplot(shared_psu_data, aes(x = SA_density)) +
    geom_point(aes(y = predicted_sa_lm), color = "blue", alpha = 0.6) +
    geom_point(aes(y = predicted_sa_ratio), color = "red", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Predicted vs Actual SA Density",
         subtitle = "Blue = Linear Model, Red = Ratio Method",
         x = "Actual SA Density", y = "Predicted SA Density")
  
  
  
  
  
  # Adjusted depth bins for your data range (0 to ~30m)
  shared_psu_data <- shared_psu_data %>%
    mutate(
      depth_bin = cut(depth, 
                      breaks = c(0, 7.5, 15, 22.5, 30.1),  # Adjusted for 0-30m range
                      labels = c("0-7.5m", "7.5-15m", "15-22.5m", "22.5-30m"),
                      include.lowest = TRUE)
    )
  
  # Calculate correlation by depth bin
  depth_correlations <- shared_psu_data %>%
    group_by(depth_bin) %>%
    summarise(
      n = n(),
      correlation = cor(SA_density, cover, use = "complete.obs"),
      mean_depth = mean(depth, na.rm = TRUE),
      mean_sa_density = mean(SA_density, na.rm = TRUE),
      mean_cover = mean(cover, na.rm = TRUE),
      .groups = 'drop'
    )
  
  print(depth_correlations)
  
  # Modern approach using nest() and map()
    
  depth_models <- shared_psu_data %>%
    filter(!is.na(depth_bin)) %>%
    group_by(depth_bin) %>%
    nest() %>%
    mutate(
      model = map(data, ~ lm(SA_density ~ cover, data = .x)),
      model_summary = map(model, summary),
      model_glance = map(model, glance),
      r_squared = map_dbl(model_glance, "r.squared"),
      slope = map_dbl(model, ~ coef(.x)[2]),
      intercept = map_dbl(model, ~ coef(.x)[1]),
      n_obs = map_dbl(data, nrow)
    )
  
  # View the results
  depth_models %>% 
    select(depth_bin, n_obs, r_squared, slope, intercept) %>%
    print()
  
  # Recalculate quartiles based on actual data range
  depth_quartiles <- quantile(shared_psu_data$depth, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  print("Actual depth quartiles:")
  print(round(depth_quartiles, 2))
  
  shared_psu_data <- shared_psu_data %>%
    mutate(
      depth_quartile = cut(depth, 
                           breaks = depth_quartiles,
                           labels = c("Q1 (Shallow)", "Q2", "Q3", "Q4 (Deep)"),
                           include.lowest = TRUE)
    )
  
  quartile_correlations <- shared_psu_data %>%
    group_by(depth_quartile) %>%
    summarise(
      n = n(),
      depth_range = paste(round(min(depth, na.rm = TRUE), 1), "-", 
                          round(max(depth, na.rm = TRUE), 1), "m"),
      correlation = cor(SA_density, cover, use = "complete.obs"),
      .groups = 'drop'
    )
  
  print(quartile_correlations)
  
  # Updated plots with better depth range
  # Plot 1: Correlation vs depth
  ggplot(depth_correlations, aes(x = mean_depth, y = correlation)) +
    geom_point(aes(size = n), alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = "Correlation between SA Density and Cover by Depth",
         x = "Mean Depth (m)", 
         y = "Correlation Coefficient",
         size = "Sample Size") +
    xlim(0, 32) +  # Adjusted x-axis limit
    theme_minimal()
  
  # Plot 2: Faceted scatter plots by depth bin
  ggplot(shared_psu_data, aes(x = SA_density, y = cover)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    facet_wrap(~depth_bin, scales = "free") +
    labs(title = "SA Density vs Cover Relationship by Depth",
         x = "SA Density", 
         y = "Cover (%)") +
    theme_minimal()
  
  # Plot 3: 3D visualization of depth effect
    
  plot_3d <- plot_ly(shared_psu_data, 
                     x = ~SA_density, 
                     y = ~cover, 
                     z = ~depth,
                     color = ~depth,
                     type = "scatter3d",
                     mode = "markers") %>%
    layout(title = "3D View: SA Density, Cover, and Depth (0-30m)",
           scene = list(xaxis = list(title = "SA Density"),
                        yaxis = list(title = "Cover (%)"),
                        zaxis = list(title = "Depth (m)", range = c(0, 32))))
  
  plot_3d
  
  # Rest of the statistical testing remains the same
  depth_interaction_model <- lm(SA_density ~ cover * depth, data = shared_psu_data)
  simple_model <- lm(SA_density ~ cover + depth, data = shared_psu_data)
  no_depth_model <- lm(SA_density ~ cover, data = shared_psu_data)
  
  # Compare models
  anova(no_depth_model, simple_model, depth_interaction_model)
  
  # Summary of interaction model
  summary(depth_interaction_model)
  
  # Calculate R-squared improvement
  cat("R-squared without depth:", round(summary(no_depth_model)$r.squared, 3), "\n")
  cat("R-squared with depth:", round(summary(simple_model)$r.squared, 3), "\n")
  cat("R-squared with interaction:", round(summary(depth_interaction_model)$r.squared, 3), "\n")
  
  