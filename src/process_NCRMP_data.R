  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())

  library(here)
  library(tidyverse)
  
  source(here("src/functions.R"))
  
  ################################## setup ##################################
  
  # load_spat_objects(directory = 'output/output_import_merge_rasters/') #call function
  # load(here('output', 'output_import_merge_rasters/import_merge_rasters_workspace.RData')) #load workspace from upstream script
  
  ################################## load NCRMP data ##################################
  
  load(here('data/USVI_2013_benthic_cover.rda'))
  load(here('data/USVI_2015_benthic_cover.rda'))
  load(here('data/USVI_2017_benthic_cover.rda'))
  load(here('data/USVI_2013_coral_demographics.rda'))
  load(here('data/USVI_2015_coral_demographics.rda'))
  load(here('data/USVI_2017_coral_demographics.rda'))
  
  ################################## summarize cover & SA ##################################
  
  USVI_benthic_cover = bind_rows(
    USVI_2013_benthic_cover,
    USVI_2015_benthic_cover,
    USVI_2017_benthic_cover
  )
  
  USVI_demo = bind_rows(
    USVI_2013_coral_demographics,
    USVI_2015_coral_demographics,
    USVI_2017_coral_demographics
  )
  
  #refactor variables and remove unnecessary ones
  USVI_benthic_cover = USVI_benthic_cover %>%
    mutate(PRIMARY_SAMPLE_UNIT = as.factor(PRIMARY_SAMPLE_UNIT)) %>%
    select(-REGION, -STATION_NR, -RUGOSITY_CD, -WTD_RUG, -MAPGRID_NR, -HABITAT_CD, -STRAT, -SUB_REGION_NAME, -SUB_REGION_NR,
           -ZONE_NAME, -ZONE_NR, -MPA_NAME, -MPA_NR, -ADMIN, -PROT, -DEPTH_STRAT,
           -METERS_COMPLETED, -MIN_DEPTH, -MAX_DEPTH)
  USVI_demo = USVI_demo %>%
    mutate(PRIMARY_SAMPLE_UNIT = as.factor(PRIMARY_SAMPLE_UNIT)) %>%
    select(-REGION, -STATION_NR, -RUGOSITY_CD, -WTD_RUG, -MAPGRID_NR, -HABITAT_CD, -STRAT, -SUB_REGION_NAME, -SUB_REGION_NR,
           -ZONE_NAME, -ZONE_NR, -MPA_NAME, -MPA_NR, -ADMIN, -PROT, -DEPTH_STRAT,
           -METERS_COMPLETED, -MIN_DEPTH, -MAX_DEPTH, -N, -JUV, -BLEACH_CONDITION, -DISEASE)
  
  #refactor date
  USVI_benthic_cover = USVI_benthic_cover %>%
    mutate(
      date = as.POSIXct(
        paste(sprintf("%02d.%02d.%02d", MONTH, DAY, YEAR %% 100)), 
        format = "%m.%d.%y"
      )
    ) %>%
    select(-YEAR, -MONTH, -DAY) %>%
    select(date, everything())
  USVI_demo = USVI_demo %>%
    mutate(
      date = as.POSIXct(
        paste(sprintf("%02d.%02d.%02d", MONTH, DAY, YEAR %% 100)), 
        format = "%m.%d.%y"
      )
    ) %>%
    select(-YEAR, -MONTH, -DAY) %>%
    select(date, everything())
  
  #calculate cover
  USVI_benthic_cover <- USVI_benthic_cover %>%
    mutate(cover = HARDBOTTOM_P + SOFTBOTTOM_P + RUBBLE_P) %>%
    select(-HARDBOTTOM_P, -SOFTBOTTOM_P, -RUBBLE_P)
  
  #calculate susceptible tissue surface area (SA)
  # NOTE - could consider trimming down predicted SA for OANN since it has a lot of dead space in reality
  #
  # remove corals that had recently suffered complete mortality upon survey, but keep NAs
  USVI_demo = USVI_demo %>%
    filter(is.na(OLD_MORT + RECENT_MORT) | OLD_MORT + RECENT_MORT < 100)
  # set '0' cm heights as an arbitrarily small amount (cm) to allow the hemi-ellipsoid estimation to function correctly. coral recruits 
  #   have very little tangible height, and were recorded underwater as 0 cm height. there is also
  #   an instance of '-1' height but it is an acroporid and gets filtered out downstream anyways
  USVI_demo = USVI_demo %>%
    mutate(HEIGHT = ifelse(HEIGHT == 0, 0.01, HEIGHT))
  #
  # hemi-ellipsoid estimation (Knud Thomsen approximation; see Xu 2009 but also Holstein 2015).
  #   p is a dimensionless constant; all else in square cm
  USVI_demo = USVI_demo %>%
    mutate(
      a = HEIGHT,
      b = PERP_DIAMETER / 2,
      c = MAX_DIAMETER / 2,
      p = 1.6075,
      SA_colony = 2 * pi * (((a * b)^p + (a * c)^p + (b * c)^p) / 3)^(1 / p),
      SA_colony = SA_colony / 10000,  # Convert from square cm to square meters
      SA = SA_colony*(1-(OLD_MORT + RECENT_MORT)/100)
    ) %>%
    select(-a, -b, -c, -p)
  #
  # remove unnecessary variables
  USVI_demo = USVI_demo %>%
    select(-MAX_DIAMETER, -PERP_DIAMETER, -HEIGHT, -OLD_MORT, -RECENT_MORT, -SA_colony)
  
  #rename variables
  USVI_benthic_cover <- USVI_benthic_cover %>%
    rename_with(tolower) %>% #convert to lowercase
    rename_with(~ tolower(gsub("_", "", .))) %>% #remove underscores
    rename(
      PSU = primarysampleunit,
      lat = latdegrees,
      lon = londegrees,
      code = covercatcd,
      spp = covercatname
    )
  USVI_demo <- USVI_demo %>%
    rename_with(tolower) %>% #convert to lowercase
    rename_with(~ tolower(gsub("_", "", .))) %>% #remove underscores
    rename(
      PSU = primarysampleunit,
      lat = latdegrees,
      lon = londegrees,
      code = speciescd,
      spp = speciesname
    )
  
  #preserve absence data
  # NOTE - this is a little tricky with demo data since possible species codes
  #         expanded with each sampling year. should be fine since we are binning species coarsely
  #         by susceptibility group, though
  all_PSU_info_cover <- USVI_benthic_cover %>%
    distinct(PSU, lat, lon, date)
  all_PSU_info_demo <- USVI_demo %>%
    distinct(PSU, lat, lon, date, spp)
  
  #filter out non-scleractinian cover
  # NOTE - removing 'OTH SPE.' - even if this might include corals, we don't know which species
  #         - also removing 'OTH CORA' - again because we do not know which coral this was
  levels(USVI_benthic_cover$code)
  USVI_benthic_cover = USVI_benthic_cover %>%
    filter(!code %in% c('BAR SUB.', 'CLI SPE.', 'CYA SPE.', 'DIC SPE.', 'GOR ENCR', 'GOR GORG', 'HAL SPE.',
                              'LOB SPE.', 'MAC CALC', 'MAC FLES', 'MAG SPE.', 'MIL SPE.', 'PAL SPE.', 'PEY SPE.',
                              'RHO CRUS', 'SPO OTHE', 'TUR FREE', 'TUR SEDI', 'RAM SPE.',
                              'OTH SPE.', 'OTH CORA')) %>%
    mutate(code = droplevels(code),
           spp = droplevels(spp)) #drop factor levels which no longer are associated with any data
  #filter out unidentified scleractinians
  levels(USVI_demo$code)
  USVI_demo = USVI_demo %>%
    filter(!code %in% c('SCL SPE.')) %>%
    mutate(code = droplevels(code),
           spp = droplevels(spp)) #drop factor levels which no longer are associated with any data
  
  #break corals into susceptibility groups
  # Questionable corals:
  # - PCLI
  # - Scolymia
  # - Agaricia
  # - Madracis
  # - Helioceris
  # - Manicina
  # - Siderastrea radians
  # - Favia fragum
  # - Isophyllia
  # - Solenastrea
  # - Tubastrea coccinea
  # - Oculina
  # - Porites
  # - Mycetophyllia
  levels(USVI_benthic_cover$code)
  USVI_benthic_cover = USVI_benthic_cover %>%
    mutate(
      susc = case_when(
        code %in% c('AGA AGAR', 'AGA FRAG', 'AGA GRAH', 'AGA HUMI', 'AGA LAMA', 'AGA SPE.', 'MAD AURE',
                   'MAD DECA', 'MAD SPE.', 'POR ASTR', 'POR DIVA', 'POR FURC', 'POR PORI', 'POR SPE.',
                   'SID RADI', 'SID SIDE', 'SID SPE.', 'STE INTE', 'AGA TENU', 'POR COLO') ~ 'low',
        code %in% c('MON CAVE', 'ORB ANNU', 'ORB FAVE', 'ORB FRAN', 'ORB SPE.', 'SOL BOUR', 'SOL SPE.',
                   'ORB ANCX') ~ 'moderate',
        code %in% c('COL NATA', 'DEN CYLI', 'DIC STOK', 'DIP LABY', 'EUS FAST', 'MEA MEAN', 'MYC ALIC', 'MYC FERO', 
                   'PSE CLIV', 'PSE SPE.', 'PSE STRI', 'MEA JACK', 'MYC REES', 'MEA SPE.') ~ 'high',
        code %in% c('ACR CERV', 'ACR PALM') ~ 'Unaffected',
        code %in% c('SCO CUBE', 'SCO SPE.', 'HEL CUCU', 'MAN AREO', 'FAV FRAG', 'ISO RIGI', 'ISO SINU',
                   'TUB COCC') ~ 'Unknown'
      )
    )
  #
  # Find codes in demo but not in benthic_cover
  levels(USVI_demo$code)
  benthic_levels <- levels(USVI_benthic_cover$code)
  demo_levels <- levels(USVI_demo$code)
  in_demo_only <- setdiff(demo_levels, benthic_levels)
  cat("Coral codes in USVI_demo but not in USVI_benthic_cover:\n")
  print(in_demo_only)
  #
  USVI_demo = USVI_demo %>%
    mutate(
      susc = case_when(
        code %in% c('AGA AGAR', 'AGA FRAG', 'AGA GRAH', 'AGA HUMI', 'AGA LAMA', 'AGA SPE.', 'MAD AURE',
                    'MAD DECA', 'MAD SPE.', 'POR ASTR', 'POR DIVA', 'POR FURC', 'POR PORI', 'POR SPE.',
                    'SID RADI', 'SID SIDE', 'SID SPE.', 'STE INTE', 'AGA TENU', 'POR COLO',
                    'MAD PHAR', 'POR BRAN', 'MAD FORM') ~ 'low',
        code %in% c('MON CAVE', 'ORB ANNU', 'ORB FAVE', 'ORB FRAN', 'ORB SPE.', 'SOL BOUR', 'SOL SPE.',
                    'ORB ANCX') ~ 'moderate',
        code %in% c('COL NATA', 'DEN CYLI', 'DIC STOK', 'DIP LABY', 'EUS FAST', 'MEA MEAN', 'MYC ALIC', 'MYC FERO', 
                    'PSE CLIV', 'PSE SPE.', 'PSE STRI', 'MEA JACK', 'MYC REES', 'MEA SPE.',
                    'MYC SPE.', 'MEA DANA', 'MYC DANA', 'MYC LAMA') ~ 'high',
        code %in% c('ACR CERV', 'ACR PALM', 'ACR PROL') ~ 'Unaffected',
        code %in% c('SCO CUBE', 'SCO SPE.', 'HEL CUCU', 'MAN AREO', 'FAV FRAG', 'ISO RIGI', 'ISO SINU',
                    'TUB COCC', 'OCU DIFF', 'MUS ANGU', 'SCO LACE', 'ISO SPE.', 'OCU SPE.') ~ 'Unknown'
      )
    )
  
  #summaries
  cover_by_susc <- USVI_benthic_cover %>%
    group_by(PSU, susc) %>%
    summarise(total_cover = sum(cover, na.rm = TRUE), .groups = 'drop')
  #
  total_cover_by_PSU <- USVI_benthic_cover %>%
    group_by(PSU) %>%
    summarise(total_cover = sum(cover, na.rm = TRUE), .groups = 'drop')
  #
  SA_by_susc <- USVI_demo %>%
    group_by(PSU, susc) %>%
    summarise(total_SA = sum(sa, na.rm = TRUE), .groups = 'drop')
  #
  total_SA_by_PSU <- USVI_demo %>%
    group_by(PSU) %>%
    summarise(total_SA = sum(sa, na.rm = TRUE), .groups = 'drop')
  
  
  # STOPPING POINT - 2 June 2025. am about to produce plots of SA, then should figure out 'absence' for
  #   demo surveys, then will collate the demo & LPI surveys together, then append bathymetry & derived
  #   values to dataframe, then produce a simple test GAM. finally will need to figure out prediction
  #   onto a uniform grid which will "talk" to the habitat grid plugged into the CMS
  
  
  # Updated summaries that include zero-cover PSU's
  cover_by_susc_complete <- USVI_benthic_cover %>%
    group_by(PSU, susc) %>%
    summarise(total_cover = sum(cover, na.rm = TRUE), .groups = 'drop') %>%
    complete(PSU = all_PSU_info_cover$PSU, 
             susc = c('low', 'moderate', 'high', 'Unaffected', 'Unknown'), 
             fill = list(total_cover = 0)) %>%
    left_join(all_PSU_info_cover, by = "PSU")  # Add back location info
  
  total_cover_by_PSU_complete <- USVI_benthic_cover %>%
    group_by(PSU) %>%
    summarise(total_cover = sum(cover, na.rm = TRUE), .groups = 'drop') %>%
    right_join(all_PSU_info_cover, by = "PSU") %>%
    mutate(total_cover = replace_na(total_cover, 0))
  
  hist(total_cover_by_PSU_complete$total_cover)
  
  ################################## plot total cover ##################################
  
  # Alternative version with density overlay
  ggplot(total_cover_by_PSU_complete, aes(x = total_cover)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_density(color = "red", size = 1) +
    scale_x_continuous(
      breaks = seq(0, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE), by = 5),
      minor_breaks = seq(0, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE), by = 1)
    ) +
    labs(
      title = "Distribution of Total Coral Cover by PSU",
      subtitle = "With density curve overlay",
      x = "Total Coral Cover (%)",
      y = "Density"
    ) +
    theme_minimal()
  
  # Bar plot version for ultimate clarity at low values
  low_cover_data <- total_cover_by_PSU_complete %>%
    filter(total_cover <= 5) %>%
    mutate(cover_rounded = round(total_cover * 4) / 4)  # Round to nearest 0.25%
  
  ggplot(low_cover_data, aes(x = factor(cover_rounded))) +
    geom_bar(fill = "darkgreen", alpha = 0.7, color = "white") +
    labs(
      title = "Distribution of Low Coral Cover Sites (≤5%)",
      subtitle = "Each bar represents 0.25% intervals - perfect for distinguishing 0% vs 1% sites",
      x = "Total Coral Cover (%) - Rounded to nearest 0.25%",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10)
    )
  library(ggplot2)
  library(dplyr)
  
  # Basic detailed histogram
  ggplot(total_cover_by_PSU_complete, aes(x = total_cover)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
    scale_x_continuous(
      breaks = seq(0, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE), by = 5),
      minor_breaks = seq(0, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE), by = 1)
    ) +
    labs(
      title = "Distribution of Total Coral Cover by PSU",
      subtitle = "USVI Benthic Cover Data (2013, 2015, 2017)",
      x = "Total Coral Cover (%)",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Alternative version with density overlay
  ggplot(total_cover_by_PSU_complete, aes(x = total_cover)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_density(color = "red", size = 1) +
    scale_x_continuous(
      breaks = seq(0, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE), by = 5),
      minor_breaks = seq(0, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE), by = 1)
    ) +
    labs(
      title = "Distribution of Total Coral Cover by PSU",
      subtitle = "With density curve overlay",
      x = "Total Coral Cover (%)",
      y = "Density"
    ) +
    theme_minimal()
  
  # Version with fine bin width to distinguish low values clearly
  ggplot(total_cover_by_PSU_complete, aes(x = total_cover)) +
    geom_histogram(binwidth = 0.5, fill = "coral", alpha = 0.7, color = "white", size = 0.3) +
    scale_x_continuous(
      breaks = seq(0, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE), by = 1),
      minor_breaks = seq(0, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE), by = 0.5),
      limits = c(-0.25, max(total_cover_by_PSU_complete$total_cover, na.rm = TRUE) + 1)
    ) +
    labs(
      title = "Distribution of Total Coral Cover by PSU",
      subtitle = "Bin width = 0.5% to clearly distinguish zero vs low cover sites",
      x = "Total Coral Cover (%)",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor.x = element_line(color = "grey90", size = 0.3),
      panel.grid.major.x = element_line(color = "grey70", size = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )
  
  # Zoomed-in version focusing on low values (0-10%)
  ggplot(total_cover_by_PSU_complete, aes(x = total_cover)) +
    geom_histogram(binwidth = 0.25, fill = "steelblue", alpha = 0.7, color = "white", size = 0.3) +
    scale_x_continuous(
      breaks = seq(0, 10, by = 0.5),
      minor_breaks = seq(0, 10, by = 0.25),
      limits = c(-0.125, 10)
    ) +
    labs(
      title = "Distribution of Coral Cover: Focus on Low Values (0-10%)",
      subtitle = "Bin width = 0.25% to clearly separate zero, near-zero, and low cover",
      x = "Total Coral Cover (%)",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor.x = element_line(color = "grey90", size = 0.3),
      panel.grid.major.x = element_line(color = "grey60", size = 0.5),
      axis.text.x = element_text(size = 10)
    )
  
  # Bar plot version for ultimate clarity at low values
  low_cover_data <- total_cover_by_PSU_complete %>%
    filter(total_cover <= 5) %>%
    mutate(cover_rounded = round(total_cover * 4) / 4)  # Round to nearest 0.25%
  
  ggplot(low_cover_data, aes(x = factor(cover_rounded))) +
    geom_bar(fill = "darkgreen", alpha = 0.7, color = "white") +
    labs(
      title = "Distribution of Low Coral Cover Sites (≤5%)",
      subtitle = "Each bar represents 0.25% intervals - perfect for distinguishing 0% vs 1% sites",
      x = "Total Coral Cover (%) - Rounded to nearest 0.25%",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10)
    )
  
  # Summary statistics to help interpret the histogram
  summary_stats <- total_cover_by_PSU_complete %>%
    summarise(
      mean_cover_all = mean(total_cover, na.rm = TRUE),
      mean_cover_nonzero = mean(total_cover[total_cover > 0], na.rm = TRUE),
      median_cover = median(total_cover, na.rm = TRUE),
      median_cover_nonzero = median(total_cover[total_cover > 0], na.rm = TRUE),
      sd_cover_nonzero = sd(total_cover[total_cover > 0], na.rm = TRUE),
      min_cover = min(total_cover, na.rm = TRUE),
      max_cover = max(total_cover, na.rm = TRUE),
      zero_cover_sites = sum(total_cover == 0, na.rm = TRUE),
      nonzero_cover_sites = sum(total_cover > 0, na.rm = TRUE),
      total_sites = n(),
      percent_zero_sites = round(100 * zero_cover_sites / total_sites, 1)
    )
  
  print(summary_stats)
  
  # Additional breakdown by cover categories
  cover_categories <- total_cover_by_PSU_complete %>%
    mutate(
      cover_category = case_when(
        total_cover == 0 ~ "Zero cover (0%)",
        total_cover > 0 & total_cover <= 1 ~ "Very low (0-1%)",
        total_cover > 1 & total_cover <= 5 ~ "Low (1-5%)",
        total_cover > 5 & total_cover <= 10 ~ "Moderate (5-10%)",
        total_cover > 10 & total_cover <= 25 ~ "High (10-25%)",
        total_cover > 25 ~ "Very high (>25%)"
      )
    ) %>%
    count(cover_category) %>%
    mutate(percentage = round(100 * n / sum(n), 1))
  
  print("Cover category breakdown:")
  print(cover_categories)
  
  ################################## plot by susceptibility ##################################
  
  # Summary statistics by susceptibility group
  susc_summary <- cover_by_susc_complete %>%
    mutate(susc = factor(susc, levels = c("low", "moderate", "high", "Unaffected", "Unknown"))) %>%
    group_by(susc) %>%
    summarise(
      mean_cover_all = mean(total_cover, na.rm = TRUE),
      mean_cover_nonzero = mean(total_cover[total_cover > 0], na.rm = TRUE),
      median_cover_nonzero = median(total_cover[total_cover > 0], na.rm = TRUE),
      sd_cover_nonzero = sd(total_cover[total_cover > 0], na.rm = TRUE),
      min_cover = min(total_cover, na.rm = TRUE),
      max_cover = max(total_cover, na.rm = TRUE),
      zero_cover_sites = sum(total_cover == 0, na.rm = TRUE),
      nonzero_cover_sites = sum(total_cover > 0, na.rm = TRUE),
      total_sites = n(),
      percent_zero_sites = round(100 * zero_cover_sites / total_sites, 1),
      .groups = 'drop'
    )
  
  print("Summary statistics by susceptibility group:")
  print(susc_summary)
  
  # Histogram by susceptibility group
  ggplot(cover_by_susc_complete %>% 
           mutate(susc = factor(susc, levels = c("low", "moderate", "high", "Unaffected", "Unknown"))), 
         aes(x = total_cover, fill = susc)) +
    geom_histogram(binwidth = 0.5, alpha = 0.7, color = "white", size = 0.2) +
    facet_wrap(~susc, scales = "free_y", ncol = 2) +
    scale_x_continuous(
      breaks = seq(0, max(cover_by_susc_complete$total_cover, na.rm = TRUE), by = 2),
      minor_breaks = seq(0, max(cover_by_susc_complete$total_cover, na.rm = TRUE), by = 1)
    ) +
    scale_fill_viridis_d(name = "Susceptibility") +
    labs(
      title = "Distribution of Coral Cover by Susceptibility Group",
      subtitle = "Bin width = 0.5% for each susceptibility category",
      x = "Coral Cover (%) by Susceptibility Group",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor.x = element_line(color = "grey90", size = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      strip.text = element_text(size = 10, face = "bold")
    )
  
  # Box plot comparing susceptibility groups
  ggplot(cover_by_susc_complete %>% 
           filter(total_cover > 0) %>%
           mutate(susc = factor(susc, levels = c("low", "moderate", "high", "Unaffected", "Unknown"))), 
         aes(x = susc, y = total_cover, fill = susc)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
    scale_fill_viridis_d(name = "Susceptibility") +
    labs(
      title = "Coral Cover by Susceptibility Group (Non-zero sites only)",
      subtitle = "Box plots with individual data points overlaid",
      x = "Susceptibility Group",
      y = "Coral Cover (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Stacked bar chart showing proportion of zero vs non-zero sites by susceptibility
  susc_zero_nonzero <- cover_by_susc_complete %>%
    mutate(
      susc = factor(susc, levels = c("low", "moderate", "high", "Unaffected", "Unknown")),
      has_cover = ifelse(total_cover > 0, "Has coral cover", "Zero coral cover")
    ) %>%
    count(susc, has_cover) %>%
    group_by(susc) %>%
    mutate(
      total_sites = sum(n),
      percentage = round(100 * n / total_sites, 1)
    )
  
  ggplot(susc_zero_nonzero, aes(x = susc, y = n, fill = has_cover)) +
    geom_col(position = "stack", alpha = 0.8) +
    geom_text(aes(label = paste0(n, "\n(", percentage, "%)")), 
              position = position_stack(vjust = 0.5), 
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("Zero coral cover" = "grey60", "Has coral cover" = "coral"),
      name = "Cover Status"
    ) +
    labs(
      title = "Sites with vs without Coral Cover by Susceptibility Group",
      subtitle = "Numbers and percentages shown for each category",
      x = "Susceptibility Group",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )