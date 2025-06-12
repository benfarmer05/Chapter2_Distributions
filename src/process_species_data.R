  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(readxl)
  
  ################################## setup ##################################
  
  # NOTE - in general, should I be including conditions (bleached, diseased, etc.) as "coral cover"?
  #          - also will need to figure out how compatible each survey type is. some of these are
  #             long-term monitoring rather than random, which could mean they artificially inflate
  #             (or reduce the quality of output of) coral cover due to bias toward rich sites
  #      - also, for LTR sites we need to consider that they have repeated measures. so, we may want
  #           to average across observations?
  
  # NOTE - why are there no corals marked properly as JUV in NCRMP data?
  
  # NOTE - what is the meaning of intercept methods in the TCRMP data?
  
  # NOTE - need to add in PR NCRMP!
  
  # NOTE - 3 June 2025 - other datasets possibly to include:
  #
  # - NODICE demo
  # - SESAP demo
  # - DeepLion demo
  # - VINPS / CSUN
  # - PR FEMA (2018)
  # - PR DNER (1999 to 2020)
  # - NCCOS pre-NCRMP in VI/PR (2001 to 2012)
  # - EPA for VI (2006 to 2009)
  #
  # - CSUN Pete Edmunds data is great, but a lot of it is pre-NCRMP. do we want to mix years like that?
  #     - good bets might be 'Population Projections' or 'Landscape-scale' from:
  #         https://coralreefs.csun.edu/data/data-catalogue/
  #     - but in particular, species-level data 1987-present is located in first hit
  #         (Virgin Islands National Park: Coral Reef: Population Dynamics: Scleractinian corals)
  #         on that page
  #     - consider though that the above does not include lat/lons. they should roughly match with
  #         those found in:
  #         'Coral cover at six sites on the south coast of St. John, USVI from 1992 to 2019', found in
  #         https://www.bco-dmo.org/project/2272, though (from https://coralreefs.csun.edu/data/)
  
  load(here('data/USVI_2013_benthic_cover.rda'))
  load(here('data/USVI_2015_benthic_cover.rda'))
  load(here('data/USVI_2017_benthic_cover.rda'))
  load(here('data/USVI_2013_coral_demographics.rda'))
  load(here('data/USVI_2015_coral_demographics.rda'))
  load(here('data/USVI_2017_coral_demographics.rda'))
  DCRMP = read_xlsx(here("data/DCRMP_Master_BenthicCover_Apr2022.xlsx"), sheet = 'BenthicData')
  DCRMP_metadata = read_xlsx(here("data/DCRMP_Master_BenthicCover_Apr2022.xlsx"), sheet = 'SiteMetadata')
  DCRMP_codes = read_xlsx(here("data/DCRMP_Master_BenthicCover_Apr2022.xlsx"), sheet = 'BenthicCodes')
  SESAP = read_xls(here("data/SESAP_Benthic MasterFINAL_toJMP_Zs_150122.xls"))
  NODICE_45m = read_xlsx(here("data/NODICE_BenthicMaster_45m.xlsx"))
  NODICE_60m = read_xlsx(here("data/NODICE_BenthicMaster_60m.xlsx"))
  NODICE_100m = read_xlsx(here("data/NODICE_Benthic Master_100m.xlsx"))
  DeepLion = read_xlsx(here("data/DeepLion_Master_BenthicCover_Jun2019.xlsx"), sheet = 'BenthicCover')
  DeepLion_metadata = read_xlsx(here("data/DeepLion_Master_BenthicCover_Jun2019.xlsx"), sheet = 'SiteMetadata')
  DeepLion_codes = read_xlsx(here("data/DeepLion_Master_BenthicCover_Jun2019.xlsx"), sheet = 'BenthicCodes')
  CSUN_random = read.csv(here("data/CSUN_USVI_random_20241026.csv"))
  CSUN_yawzi_tektite = read.csv(here("data/CSUN_USVI_yawzi_tektite_20241025.csv"))
  TCRMP_benthic = read_xlsx(here("data/TCRMP_Master_Benthic_Cover_Feb2022.xlsx"), sheet = 'BenthicData')
  TCRMP_benthic_metadata = read_xlsx(here("data/TCRMP_Master_Benthic_Cover_Feb2022.xlsx"), sheet = 'SiteMetadata')
  TCRMP_benthic_codes = read_xlsx(here("data/TCRMP_Master_Benthic_Cover_Feb2022.xlsx"), sheet = 'BenthicCodes')
  TCRMP_health_metadata = read_xlsx(here("data/TCRMP_Master_Coral_Health_Jan2022.xlsx"), sheet = 'SiteMetadata')
  TCRMP_health_codes = read_xlsx(here("data/TCRMP_Master_Coral_Health_Jan2022.xlsx"), sheet = 'HealthCodes')
  
  #special case for TCRMP health data; both dates and the width variable are funny
  # NOTE - since we don't have a sample day variable, just setting each 'day' to first of the month of
  #         that year. coarse, but shouldn't really matter for our purposes
  # NOTE - there is at least one site where the sample year does not match the sample date. for now,
  #         just relying on the sample year variable
  # NOTE - filtering by 2013 and later to match NCRMP time domain, but can choose another cut-off
  col_types <- rep("guess", 69)  # All columns default to "guess"
  col_types[12] <- "numeric"     # Force Width (column 12) to be numeric
  col_types[68] <- "numeric"     # Force OldMort (column 68) to be numeric
  col_types[69] <- "numeric"     # Force RecMort (column 69) to be numeric
  TCRMP_health <- read_xlsx(
    here("data/TCRMP_Master_Coral_Health_Jan2022.xlsx"), 
    sheet = 'HealthDataRaw',
    col_types = col_types,
    range = cell_cols("A:BQ")  # A through BQ = 69 columns
  ) %>%
    select(Location:Height, OldMort, RecMort, SampleYear, SampleMonth) %>%  # Include the year/month columns
    mutate(
      # Create date as first day of the month using SampleYear and SampleMonth
      date = as.POSIXct(paste(SampleYear, SampleMonth, "01", sep = "-"), format = "%Y-%m-%d")
    ) %>%
    # # Remove the original year and month columns since we now have date
    # select(-SampleYear, -SampleMonth) %>%
    # Move date to the beginning
    relocate(date, .before = everything()) %>%
    # Filter for 2013 and onward
    filter(year(date) >= 2013)
  
  # NOTE / VERY IMPORTANT
  # - if I choose to go back and look at any data prior to 2013, would need to make that everything I do far downstream with wrangling
  #     species still holds up (i.e., not throwing out any important taxa-level data)
  
  ################################## wrangle CSUN ##################################
  
  CSUN_random <- CSUN_random %>%
    mutate(
      across(-c(site, quadrat, year), as.numeric),
      site = as.factor(site)
    )
  
  ################################## wrangle DeepLion ##################################
  
  # Identify metadata columns (first 7 columns based on your preview)
  metadata_cols <- c("SiteCode", "FilmDate", "MPA", "SampleZone", "NoPts", "AnalysisBy", "AnalysisDate")
  
  # Transform to long format
  DeepLion_long <- DeepLion %>%
    pivot_longer(
      cols = -all_of(metadata_cols),
      names_to = "species_code", 
      values_to = "cover"
    ) %>%
    # Convert cover to numeric and filter out NA values
    mutate(cover = as.numeric(cover)) %>%
    filter(!is.na(cover))
  
  DeepLion_long <- DeepLion_long %>%
    rename(
      location_id = SiteCode,
      date = FilmDate,
      analysis_by = AnalysisBy,
      analysis_date = AnalysisDate,
      data_points = NoPts,
      species = species_code
    ) %>%
    select(date, location_id, MPA, SampleZone, analysis_by, analysis_date, 
           data_points, species, cover) %>%
    arrange(location_id, species)
  
  #refactor locations and drop unnecessary variables
  DeepLion_long <- DeepLion_long %>%
    mutate(location_id = as.factor(location_id)) %>%
    select(-MPA, -SampleZone, -analysis_by, -analysis_date, -data_points)
  
  #bring in lat/lons
  DeepLion_long <- DeepLion_long %>%
    left_join(DeepLion_metadata %>% select(SiteCode, Latitude, Longitude), 
              by = c("location_id" = "SiteCode")) %>%
    rename(lat = Latitude, lon = Longitude)
  
  #refactor species names and drop non-coral species
  DeepLion_long <- DeepLion_long %>%
    left_join(DeepLion_codes %>% select(Code, Meaning, Category), 
              by = c("species" = "Code")) %>%
    filter(Category == "Coral") %>%
    rename(species_name = Meaning) %>%
    select(-Category, -species)
  
  # Rename columns to lowercase to match NODICE format
  DeepLion_long <- DeepLion_long %>%
    rename_with(tolower) %>%
    rename_with(~ tolower(gsub(":", "", .))) %>% #remove underscores
    rename(
      PSU = location_id,
      spp = species_name
    )
  
  ################################## wrangle NODICE ##################################
  
  # Remove the extra 'Depth (f)' row from NODICE_60m to match format of other datasets
  depth_f_row <- which(NODICE_60m[,1] == "Depth (f)")
  if(length(depth_f_row) > 0) {
    NODICE_60m <- NODICE_60m[-depth_f_row, ]
  }
  
  #remove lat/lons with no data (assuming they were not surveyed)
  date_filming_row <- 5  # NOTE - hard-coded; using NA in 'date of filming' for reference
  cols_to_keep <- as.vector(!is.na(NODICE_100m[date_filming_row, ]))
  cols_to_keep[1] <- TRUE  # Always keep the first column (species names)
  NODICE_100m <- NODICE_100m[, cols_to_keep]
  
  # Function to transform each dataset
  transform_nodice_dataset <- function(nodice_data) {
    
    metadata_rows <- c("Latitude", "Longitude", "Transect Number:", "Tape Number:", 
                       "Date of Filming:", "Transect Type:", "Depth (m)", 
                       "# of data points:", "Analysis by:", "Date of Analysis:", 
                       "Transect Length:")
    
    species_start_row <- 13 #this should be dynamic ideally; reference if there are bugs
    metadata_df <- nodice_data[1:(species_start_row-1), -1]  # Exclude first column for metadata
    species_df <- nodice_data[species_start_row:nrow(nodice_data), ]
    location_ids <- colnames(nodice_data)[-1]  # Exclude first column
    
    all_metadata <- data.frame()
    for(location_id in location_ids) {
      col_data <- metadata_df[[location_id]]
      
      # Extract key metadata - adjust row indices based on your actual data structure
      location_metadata <- data.frame(
        location_id = location_id,
        lat = as.numeric(col_data[1]),  # Latitude row
        lon = as.numeric(col_data[2]),  # Longitude row
        date = col_data[5],             # Date of Filming row
        depth_m = as.numeric(col_data[7]), # Depth row
        transect_type = col_data[6],    # Transect Type row
        analysis_date = col_data[10],   # Date of Analysis row
        data_points = as.numeric(col_data[8]) # Number of data points
      )
      
      all_metadata <- rbind(all_metadata, location_metadata)
    }
    
    # Transform species data to long format
    # Use the first column as species names, pivot the rest
    species_long <- species_df %>%
      # Create species column from first column, then select only location columns
      mutate(species = species_df[,1]) %>%
      select(-1) %>%  # Remove the first column since we saved it as 'species'
      mutate(species = species) %>%  # Add species back as a column
      pivot_longer(cols = -species, 
                   names_to = "location_id", 
                   values_to = "cover") %>%
      # Convert cover to numeric (keep all values, including 0)
      mutate(cover = as.numeric(cover)) %>%
      # Only remove rows where cover is NA (missing data)
      filter(!is.na(cover))
    
    # Join metadata with species data
    nodice_transformed <- species_long %>%
      left_join(all_metadata, by = "location_id") %>%
      # Reorder columns to match USVI format
      select(date, location_id, lat, lon, depth_m, transect_type, 
             analysis_date, data_points, species, cover) %>%
      # Sort by location and species
      arrange(location_id, species)
    
    return(nodice_transformed)
  }
  
  # Transform each dataset
  NODICE_45m_transformed <- transform_nodice_dataset(NODICE_45m)
  NODICE_60m_transformed <- transform_nodice_dataset(NODICE_60m)
  NODICE_100m_transformed <- transform_nodice_dataset(NODICE_100m)
  
  # Combine all three datasets
  NODICE_combined <- bind_rows(NODICE_45m_transformed, NODICE_60m_transformed, NODICE_100m_transformed)
  
  # Continue with the rest of the processing on the combined dataset
  #keep only corals
  NODICE_combined$species <- NODICE_combined$species$`Location:` #drop weird nested tibble
  NODICE_combined <- NODICE_combined %>% filter(grepl("- coral", species))
  
  #refactor variables and remove unnecessary ones
  NODICE_long = NODICE_combined %>%
    mutate(location_id = as.factor(location_id)) %>%
    select(-date, -transect_type, -analysis_date, -data_points) 
  
  #rename variables
  NODICE_long <- NODICE_long %>%
    rename_with(tolower) %>% #convert to lowercase
    rename_with(~ tolower(gsub(":", "", .))) %>% #remove underscores
    rename(PSU = location_id, depth = depth_m)
  
  #refactor coral names
  NODICE_long <- NODICE_long %>%
    mutate(
      spp = gsub("\\s*\\([^)]+\\)\\s*-.*$", "", species),
      spp = trimws(spp)
    ) %>%
    select(-species)
  
  ################################## wrangle SESAP ##################################
  
  #convert from wide to long format
  SESAP_long <- SESAP %>%
    pivot_longer(
      cols = -c(`Shelf Strata`, `Seaward Strata`, `Location:`, Latitude, Longitude, 
                `Transect Number:`, `Tape Number:`, `Date of Filming:`, `Transect Type:`, 
                `Depth (m)`, `# of data points:`, `Analysis by:`, `Date of Analysis:`, 
                `Transect Length`),
      names_to = "species_full_name", 
      values_to = "cover"
    ) %>%
    filter(grepl("- coral$", species_full_name, ignore.case = TRUE)) %>% #keep only corals
    mutate(
      spp = gsub("\\s*\\([^)]+\\)\\s*-.*$", "", species_full_name),
      spp = trimws(spp) #clean up species labels
    )
  
  #refactor variables and remove unnecessary ones
  SESAP_long = SESAP_long %>%
    mutate(`Location:` = as.factor(`Location:`)) %>%
    select(-`Seaward Strata`, -`Shelf Strata`, -`Transect Number:`, -`Tape Number:`, -`Transect Type:`,
           -`# of data points:`, -`Analysis by:`, -species_full_name, -`Date of Filming:`,
           -`Date of Analysis:`) 
  
  #rename variables
  SESAP_long <- SESAP_long %>%
    rename_with(tolower) %>% #convert to lowercase
    rename_with(~ tolower(gsub(":", "", .))) %>% #remove underscores
    rename(
      PSU = location,
      lat = latitude,
      lon = longitude,
      meterscompleted = `transect length`,
      depth = `depth (m)`
    )
  
  #convert depth to positive format
  SESAP_long = SESAP_long %>%
    mutate(depth = abs(depth))
  
  ################################## wrangle DCRMP ##################################
  
  # Identify metadata columns (first 6 columns based on your preview)
  metadata_cols <- c("SiteCode", "SampleYear", "SurveyDate", "NoPts", "AnalysisBy", "AnalysisDate")
  
  # Transform to long format
  DCRMP_long <- DCRMP %>%
    pivot_longer(
      cols = -all_of(metadata_cols),
      names_to = "species_code",
      values_to = "cover"
    ) %>%
    # Convert cover to numeric and filter out NA values
    mutate(cover = as.numeric(cover)) %>%
    filter(!is.na(cover))
  
  DCRMP_long <- DCRMP_long %>%
    rename(
      location_id = SiteCode,
      year = SampleYear,
      date = SurveyDate,
      analysis_by = AnalysisBy,
      analysis_date = AnalysisDate,
      data_points = NoPts,
      species = species_code
    ) %>%
    select(date, location_id, analysis_by, analysis_date,
           data_points, species, cover) %>%
    arrange(location_id, species)
  
  #refactor locations and drop unnecessary variables
  DCRMP_long <- DCRMP_long %>%
    # Remove the as.factor() conversion - keep as numeric
    select(-analysis_by, -analysis_date, -data_points)
  
  #bring in lat/lons
  # NOTE - this step introduces a many-to-many warning and makes the dataframe slightly longer.
  #         something to keep in mind
  DCRMP_long <- DCRMP_long %>%
    left_join(DCRMP_metadata %>% select(SiteCode, Latitude, Longitude, DepthM),
              by = c("location_id" = "SiteCode")) %>%
    rename(lat = Latitude, lon = Longitude, depth = DepthM)
  
  # Convert to factor after the join if you still want it as a factor
  DCRMP_long <- DCRMP_long %>%
    mutate(location_id = as.factor(location_id))
  
  #refactor species names and drop non-coral species
  DCRMP_long <- DCRMP_long %>%
    left_join(DCRMP_codes %>% select(Code, Meaning, Category),
              by = c("species" = "Code")) %>%
    filter(Category == "Coral") %>%
    rename(species_name = Meaning) %>%
    select(-Category, -species)

  # Rename columns to lowercase to match NODICE format
  DCRMP_long <- DCRMP_long %>%
    rename_with(tolower) %>%
    rename_with(~ tolower(gsub(":", "", .))) %>% #remove underscores
    rename(
      PSU = location_id,
      spp = species_name
    )
  
  # Check for duplicates in DCRMP_long data
  duplicates_dcrmp <- DCRMP_long %>%
    group_by(date, PSU, lat, lon, spp) %>%
    summarise(
      count = n(), 
      cover_values = paste(cover, collapse = ", "),
      unique_covers = n_distinct(cover),
      .groups = "drop"
    ) %>%
    filter(count > 1) %>%
    arrange(desc(count))
  
  cat("=== DCRMP DUPLICATE ANALYSIS ===\n")
  cat("Number of species-sampling event combinations with multiple records:", nrow(duplicates_dcrmp), "\n")
  
  # Show the duplicates
  print(duplicates_dcrmp)
  
  # Identify true duplicates vs multiple observations
  true_duplicates_dcrmp <- duplicates_dcrmp %>%
    filter(unique_covers == 1)
  
  multiple_observations_dcrmp <- duplicates_dcrmp %>%
    filter(unique_covers > 1)
  
  cat("\nTrue duplicates (same cover value):", nrow(true_duplicates_dcrmp), "\n")
  cat("Multiple observations (different cover values):", nrow(multiple_observations_dcrmp), "\n")
  
  # Remove duplicates - keep only the first occurrence of each combination
  DCRMP_long_cleaned <- DCRMP_long %>%
    distinct(date, PSU, lat, lon, spp, .keep_all = TRUE)
  
  cat("\n=== CLEANING RESULTS ===\n")
  cat("Original DCRMP_long rows:", nrow(DCRMP_long), "\n")
  cat("Cleaned DCRMP_long rows:", nrow(DCRMP_long_cleaned), "\n")
  cat("Rows removed:", nrow(DCRMP_long) - nrow(DCRMP_long_cleaned), "\n")
  
  # Verify no duplicates remain
  remaining_duplicates <- DCRMP_long_cleaned %>%
    group_by(date, PSU, lat, lon, spp) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(count > 1)
  
  cat("Duplicates remaining after cleaning:", nrow(remaining_duplicates), "\n")
  
  # Show first few rows of cleaned data
  cat("\n=== FIRST 10 ROWS OF CLEANED DCRMP DATA ===\n")
  print(head(DCRMP_long_cleaned, 10))
  
  DCRMP_long = DCRMP_long_cleaned
  
  ################################## wrangle TCRMP ##################################
  
  # NOTE / STOPPING POINT - 6 June 2025 - will need to fix dates for TCRMP benthic too!!
  
  # benthic
  #
  # Identify metadata columns (first 7 columns based on your preview)
  metadata_cols <- c("SampleYear", "SampleMonth", "Period", "Location", "FilmDate", "NoPts",
                     "AnalysisBy", "AnalysisDate", "Transect")
  
  # Transform to long format
  TCRMP_benthic_long <- TCRMP_benthic %>%
    pivot_longer(
      cols = -all_of(metadata_cols),
      names_to = "species_code", 
      values_to = "cover"
    ) %>%
    # Convert cover to numeric and filter out NA values
    mutate(cover = as.numeric(cover)) %>%
    filter(!is.na(cover))
  
  TCRMP_benthic_long <- TCRMP_benthic_long %>%
    rename(
      year = SampleYear,
      month = SampleMonth,
      location_id = Location,
      date = FilmDate,
      analysis_by = AnalysisBy,
      analysis_date = AnalysisDate,
      data_points = NoPts,
      species = species_code,
      transect = Transect
    ) %>%
    filter(year(date) >= 2013) %>% # NOTE - mismatch of sample year with date in Ginsburg Fringe
    select(year, month, date, transect, location_id, analysis_by, analysis_date, data_points, species, cover) %>%
    arrange(location_id, species)
  
  #refactor locations and drop unnecessary variables
  TCRMP_benthic_long <- TCRMP_benthic_long %>%
    mutate(location_id = as.factor(location_id)) %>%
    select(-year, -month, -analysis_by, -analysis_date, -data_points)
  
  #bring in lat/lons
  TCRMP_benthic_long <- TCRMP_benthic_long %>%
    left_join(TCRMP_benthic_metadata %>% select(Location, Latitude, Longitude, Depth), 
              by = c("location_id" = "Location")) %>%
    rename(lat = Latitude, lon = Longitude)
  
  #refactor species names and drop non-coral species
  TCRMP_benthic_long <- TCRMP_benthic_long %>%
    left_join(TCRMP_benthic_codes %>% select(Code, Meaning, Category), 
              by = c("species" = "Code")) %>%
    filter(Category == "Coral") %>%
    rename(species_name = Meaning) %>%
    select(-Category, -species)
  
  # Rename columns to lowercase to match NODICE format
  TCRMP_benthic_long <- TCRMP_benthic_long %>%
    rename_with(tolower) %>%
    rename_with(~ tolower(gsub(":", "", .))) %>% #remove underscores
    rename(
      PSU = location_id,
      spp = species_name
    )
  
  #demographic (health)
  #
  #shift columns to lowercase
  TCRMP_health_long <- TCRMP_health %>%
    rename_with(tolower) %>%
    rename_with(~ tolower(gsub(":", "", .))) %>% #remove underscores
    rename(
      PSU = location
    )
  
  #refactor locations and drop unnecessary variables
  TCRMP_health_long <- TCRMP_health_long %>%
    mutate(PSU = as.factor(PSU)) %>%
    select(-sampleyear, -samplemonth, -period, -sampletype, -recorder)
  
  #bring in lat/lons
  TCRMP_health_long <- TCRMP_health_long %>%
    left_join(TCRMP_health_metadata %>% select(Location, Latitude, Longitude, Depth), 
              by = c("PSU" = "Location")) %>%
    rename(lat = Latitude, lon = Longitude, depth = Depth)
  
  #refactor species names
  TCRMP_health_long <- TCRMP_health_long %>%
    left_join(TCRMP_health_codes %>% select(Code, Meaning, Category), 
              by = c("spp" = "Code")) %>%
    filter(Category == "Coral") %>%
    rename(code = spp) %>%
    rename(spp = Meaning) %>%
    select(-Category)
  
  #drop species code
  TCRMP_health_long_withcode = TCRMP_health_long
  TCRMP_health_long = TCRMP_health_long %>%
    select(-code)
  
  #### Check for duplicates in TCRMP_benthic_long data
  duplicates_TCRMP_benthic <- TCRMP_benthic_long %>%
    group_by(date, transect, PSU, lat, lon, spp) %>% # NOTE - can take out 'transect' to see how data is structured by transect
    summarise(
      count = n(), 
      cover_values = paste(cover, collapse = ", "),
      unique_covers = n_distinct(cover),
      any_non_zero = any(cover != 0),  # Check if THIS species has any non-zero values
      .groups = "drop"
    ) %>%
    filter(count > 1) %>%
    arrange(desc(count))
  #
  cat("=== TCRMP_benthic DUPLICATE ANALYSIS ===\n")
  cat("Number of species-sampling event combinations with multiple records:", nrow(duplicates_TCRMP_benthic), "\n")
  #
  # Show the duplicates
  print(duplicates_TCRMP_benthic)
  #
  # True duplicates = same cover value AND at least one non-zero value for this species
  true_duplicates_TCRMP_benthic <- duplicates_TCRMP_benthic %>%
    filter(unique_covers == 1 & any_non_zero)
  #
  # Multiple observations = different values OR all values are zero
  multiple_observations_TCRMP_benthic <- duplicates_TCRMP_benthic %>%
    filter(unique_covers > 1 | !any_non_zero)
  #
  cat("\nTrue duplicates (same cover value with non-zero values):", nrow(true_duplicates_TCRMP_benthic), "\n")
  cat("Multiple observations (different values or all zeros):", nrow(multiple_observations_TCRMP_benthic), "\n")
  
  ################################## wrangle NCRMP ##################################
  
  #benthic
  #
  NCRMP_benthic = bind_rows(
    USVI_2013_benthic_cover,
    USVI_2015_benthic_cover,
    USVI_2017_benthic_cover
  )
  
  #refactor variables and remove unnecessary ones
  NCRMP_benthic_long = NCRMP_benthic %>%
    mutate(PRIMARY_SAMPLE_UNIT = as.factor(PRIMARY_SAMPLE_UNIT)) %>%
    select(-REGION, -STATION_NR, -RUGOSITY_CD, -WTD_RUG, -MAPGRID_NR, -HABITAT_CD, -STRAT, -SUB_REGION_NAME, -SUB_REGION_NR,
           -ZONE_NAME, -ZONE_NR, -MPA_NAME, -MPA_NR, -ADMIN, -PROT, -DEPTH_STRAT)
  
  #calculate average depth (this will come from bathymetry later anyways, but interesting to compare)
  NCRMP_benthic_long = NCRMP_benthic_long %>%
    mutate(depth = (MIN_DEPTH + MAX_DEPTH)/2) %>%
    select(-MIN_DEPTH, -MAX_DEPTH)
  
  #refactor date
  NCRMP_benthic_long = NCRMP_benthic_long %>%
    mutate(
      date = as.POSIXct(
        paste(sprintf("%02d.%02d.%02d", MONTH, DAY, YEAR %% 100)), 
        format = "%m.%d.%y"
      )
    ) %>%
    select(-YEAR, -MONTH, -DAY) %>%
    select(date, everything())
  
  #calculate cover
  NCRMP_benthic_long <- NCRMP_benthic_long %>%
    mutate(cover = HARDBOTTOM_P + SOFTBOTTOM_P + RUBBLE_P) %>%
    select(-HARDBOTTOM_P, -SOFTBOTTOM_P, -RUBBLE_P)
  
  #rename variables
  NCRMP_benthic_long <- NCRMP_benthic_long %>%
    rename_with(tolower) %>% #convert to lowercase
    rename_with(~ tolower(gsub("_", "", .))) %>% #remove underscores
    rename(
      PSU = primarysampleunit,
      lat = latdegrees,
      lon = londegrees,
      code = covercatcd,
      spp = covercatname
    )
  
  #filter out non-scleractinian cover
  # NOTE - removing 'OTH SPE.' - even if this might include corals, we don't know which species
  #         - also removing 'OTH CORA' - again because we do not know which coral this was
  levels(NCRMP_benthic_long$code)
  NCRMP_benthic_long = NCRMP_benthic_long %>%
    filter(!code %in% c('BAR SUB.', 'CLI SPE.', 'CYA SPE.', 'DIC SPE.', 'GOR ENCR', 'GOR GORG', 'HAL SPE.',
                        'LOB SPE.', 'MAC CALC', 'MAC FLES', 'MAG SPE.', 'MIL SPE.', 'PAL SPE.', 'PEY SPE.',
                        'RHO CRUS', 'SPO OTHE', 'TUR FREE', 'TUR SEDI', 'RAM SPE.',
                        'OTH SPE.', 'OTH CORA')) %>%
    mutate(code = droplevels(code),
           spp = droplevels(spp)) #drop factor levels which no longer are associated with any data
  
  #drop the code entirely
  NCRMP_benthic_long_withcode = NCRMP_benthic_long
  NCRMP_benthic_long = NCRMP_benthic_long %>%
    select(-code)
  
  #demo
  #
  NCRMP_demo = bind_rows(
    USVI_2013_coral_demographics,
    USVI_2015_coral_demographics,
    USVI_2017_coral_demographics
  )
  
  #refactor variables and remove unnecessary ones
  NCRMP_demo_long = NCRMP_demo %>%
    mutate(PRIMARY_SAMPLE_UNIT = as.factor(PRIMARY_SAMPLE_UNIT)) %>%
    select(-REGION, -STATION_NR, -RUGOSITY_CD, -WTD_RUG, -MAPGRID_NR, -HABITAT_CD, -STRAT, -SUB_REGION_NAME, -SUB_REGION_NR,
           -ZONE_NAME, -ZONE_NR, -MPA_NAME, -MPA_NR, -ADMIN, -PROT, -DEPTH_STRAT, -N, -JUV,
           -BLEACH_CONDITION, -DISEASE)
  
  #calculate average depth (this will come from bathymetry later anyways, but interesting to compare)
  NCRMP_demo_long = NCRMP_demo_long %>%
    mutate(depth = (MIN_DEPTH + MAX_DEPTH)/2) %>%
    select(-MIN_DEPTH, -MAX_DEPTH)
  
  #refactor further
  NCRMP_demo_long <- NCRMP_demo_long %>%
    rename_with(tolower) %>% #convert to lowercase
    rename_with(~ tolower(gsub("_", "", .))) %>% #remove underscores
    rename(
      PSU = primarysampleunit,
      lat = latdegrees,
      lon = londegrees,
      code = speciescd,
      spp = speciesname,
      length = maxdiameter,
      width = perpdiameter
    )
  
  # Convert year, month, day to POSIXct date
  NCRMP_demo_long <- NCRMP_demo_long %>%
    mutate(
      # Create date string in format "m.d.yy" (matching your example format)
      date_string = paste0(month, ".", day, ".", str_sub(as.character(year), 3, 4)),
      # Convert to POSIXct
      date = as.POSIXct(date_string, format = "%m.%d.%y")
    ) %>%
    # Remove the original year, month, day columns and the temporary date_string
    select(-year, -month, -day, -date_string) %>%
    # Move date to leftmost position
    relocate(date, .before = everything())
  
  #drop species code
  NCRMP_demo_long_withcode = NCRMP_demo_long
  NCRMP_demo_long = NCRMP_demo_long %>%
    select(-code)
  
  # # NOTE - the below was from when I was grouping by species with just NCRMP data. now, doing so *after*
  # #         collating all datasets
  # levels(NCRMP_benthic_long$code)
  # NCRMP_benthic_long = NCRMP_benthic_long %>%
  #   mutate(
  #     susc = case_when(
  #       code %in% c('AGA AGAR', 'AGA FRAG', 'AGA GRAH', 'AGA HUMI', 'AGA LAMA', 'AGA SPE.', 'MAD AURE',
  #                  'MAD DECA', 'MAD SPE.', 'POR ASTR', 'POR DIVA', 'POR FURC', 'POR PORI', 'POR SPE.',
  #                  'SID RADI', 'SID SIDE', 'SID SPE.', 'STE INTE', 'AGA TENU', 'POR COLO') ~ 'low',
  #       code %in% c('MON CAVE', 'ORB ANNU', 'ORB FAVE', 'ORB FRAN', 'ORB SPE.', 'SOL BOUR', 'SOL SPE.',
  #                  'ORB ANCX') ~ 'moderate',
  #       code %in% c('COL NATA', 'DEN CYLI', 'DIC STOK', 'DIP LABY', 'EUS FAST', 'MEA MEAN', 'MYC ALIC', 'MYC FERO',
  #                  'PSE CLIV', 'PSE SPE.', 'PSE STRI', 'MEA JACK', 'MYC REES', 'MEA SPE.') ~ 'high',
  #       code %in% c('ACR CERV', 'ACR PALM') ~ 'Unaffected',
  #       code %in% c('SCO CUBE', 'SCO SPE.', 'HEL CUCU', 'MAN AREO', 'FAV FRAG', 'ISO RIGI', 'ISO SINU',
  #                  'TUB COCC') ~ 'Unknown'
  #     )
  #   )
  # # Find codes in demo but not in benthic_cover
  # levels(NCRMP_demo_long$code)
  # benthic_levels <- levels(NCRMP_benthic_long$code)
  # demo_levels <- levels(NCRMP_demo_long$code)
  # in_demo_only <- setdiff(demo_levels, benthic_levels)
  # cat("Coral codes in NCRMP_demo_long but not in NCRMP_benthic_long:\n")
  # print(in_demo_only)
  # #
  # NCRMP_demo_long = NCRMP_demo_long %>%
  #   mutate(
  #     susc = case_when(
  #       code %in% c('AGA AGAR', 'AGA FRAG', 'AGA GRAH', 'AGA HUMI', 'AGA LAMA', 'AGA SPE.', 'MAD AURE',
  #                   'MAD DECA', 'MAD SPE.', 'POR ASTR', 'POR DIVA', 'POR FURC', 'POR PORI', 'POR SPE.',
  #                   'SID RADI', 'SID SIDE', 'SID SPE.', 'STE INTE', 'AGA TENU', 'POR COLO',
  #                   'MAD PHAR', 'POR BRAN', 'MAD FORM') ~ 'low',
  #       code %in% c('MON CAVE', 'ORB ANNU', 'ORB FAVE', 'ORB FRAN', 'ORB SPE.', 'SOL BOUR', 'SOL SPE.',
  #                   'ORB ANCX') ~ 'moderate',
  #       code %in% c('COL NATA', 'DEN CYLI', 'DIC STOK', 'DIP LABY', 'EUS FAST', 'MEA MEAN', 'MYC ALIC', 'MYC FERO', 
  #                   'PSE CLIV', 'PSE SPE.', 'PSE STRI', 'MEA JACK', 'MYC REES', 'MEA SPE.',
  #                   'MYC SPE.', 'MEA DANA', 'MYC DANA', 'MYC LAMA') ~ 'high',
  #       code %in% c('ACR CERV', 'ACR PALM', 'ACR PROL') ~ 'Unaffected',
  #       code %in% c('SCO CUBE', 'SCO SPE.', 'HEL CUCU', 'MAN AREO', 'FAV FRAG', 'ISO RIGI', 'ISO SINU',
  #                   'TUB COCC', 'OCU DIFF', 'MUS ANGU', 'SCO LACE', 'ISO SPE.', 'OCU SPE.') ~ 'Unknown'
  #     )
  #   )
  
  
  ################################## collate benthic data ##################################
  
  # NOTE - eventually may need to consider meters completed for coral cover surveys, and/or #/points
  #   used for CPCe analysis
  
  # Function to standardize datasets
  standardize_benthic_dataset <- function(df, dataset_name) {
    # Add source dataset column
    df$dataset <- dataset_name
    
    # Standardize column names based on what's available
    if ("depth_m" %in% names(df)) {
      df <- df %>% rename(depth = depth_m)
    }
    
    # Ensure all datasets have the same core columns
    # If date column doesn't exist, set to NA
    if (!"date" %in% names(df)) {
      df$date <- NA
    }
    
    # If depth column doesn't exist, set to NA
    if (!"depth" %in% names(df)) {
      df$depth <- NA
    }
    
    # Select and reorder columns consistently
    df <- df %>%
      select(dataset, date, PSU, lat, lon, cover, spp, depth) %>%
      # Move dataset to leftmost position
      relocate(dataset, .before = everything())
    
    return(df)
  }
  
  # Standardize each dataset
  DCRMP_std <- standardize_benthic_dataset(DCRMP_long, "DCRMP")
  DeepLion_std <- standardize_benthic_dataset(DeepLion_long, "DeepLion")
  NODICE_std <- standardize_benthic_dataset(NODICE_long, "NODICE")
  SESAP_std <- standardize_benthic_dataset(SESAP_long, "SESAP")
  TCRMP_std <- standardize_benthic_dataset(TCRMP_benthic_long, "TCRMP_benthic")
  NCRMP_std = standardize_benthic_dataset(NCRMP_benthic_long, "NCRMP_benthic")
  
  # Combine all datasets
  combined_benthic_data <- bind_rows(
    DCRMP_std,
    DeepLion_std,
    NODICE_std,
    SESAP_std,
    TCRMP_std,
    NCRMP_std
  )
  
  # Convert spp to factor
  combined_benthic_data$spp <- as.factor(combined_benthic_data$spp)
  
  # Show all species levels
  print("\nAll species levels:")
  spp_levels <- levels(combined_benthic_data$spp)
  print(spp_levels)
  print(paste("Total number of species:", length(spp_levels)))
  
  ### check duplicates - should just be TCRMP
  duplicates <- combined_benthic_data %>%
    group_by(dataset, PSU, date, lat, lon, spp) %>%
    summarise(count = n(), cover_values = paste(cover, collapse = ", "), .groups = "drop") %>%
    filter(count > 1)
  #
  # See how many sampling events you actually have
  n_sampling_events <- combined_benthic_data %>%
    select(dataset, PSU, date, lat, lon) %>%
    distinct() %>%
    nrow()
  #
  # See how many unique species
  n_species <- length(unique(combined_benthic_data$spp))
  #
  cat("Unique sampling events:", n_sampling_events)
  cat("Unique species:", n_species) 
  cat("Expected complete grid size:", n_sampling_events * n_species)
  
  ### Check for repeat (different dates) visits to same PSUs across datasets
  # NOTE - all looks fine. a few repeat visits in the DeepLion/DCRMP for whatever reason
  cat("\n\n=== MULTIPLE PSU VISITS CHECK ===\n")
  #
  # For datasets with dates: DeepLion_long, DCRMP_long, NCRMP_benthic_long
  datasets_with_dates <- c("DeepLion", "DCRMP", "NCRMP_benthic")
  #
  for(dataset_name in datasets_with_dates) {
    dataset_data <- combined_benthic_data %>%
      filter(dataset == dataset_name)
    
    if(nrow(dataset_data) > 0) {
      multiple_visits <- dataset_data %>%
        select(PSU, date) %>%
        distinct() %>%
        group_by(PSU) %>%
        summarise(
          visit_count = n(),
          dates = paste(unique(date), collapse = ", "),
          .groups = "drop"
        ) %>%
        filter(visit_count > 1) %>%
        arrange(desc(visit_count))
      
      cat(paste0("\n", dataset_name, " - PSUs visited multiple times: ", nrow(multiple_visits), "\n"))
      if(nrow(multiple_visits) > 0) {
        print(multiple_visits)
      }
    } else {
      cat(paste0("\n", dataset_name, " - No data found\n"))
    }
  }
  #
  # For datasets without dates: SESAP_long, NODICE_long (check for >61 entries per PSU)
  datasets_without_dates <- c("SESAP", "NODICE")
  #
  for(dataset_name in datasets_without_dates) {
    dataset_data <- combined_benthic_data %>%
      filter(dataset == dataset_name)
    
    if(nrow(dataset_data) > 0) {
      psu_counts <- dataset_data %>%
        group_by(PSU) %>%
        summarise(
          entry_count = n(),
          species_count = n_distinct(spp),
          .groups = "drop"
        ) %>%
        filter(entry_count > 61) %>%
        arrange(desc(entry_count))
      
      cat(paste0("\n", dataset_name, " - PSUs with >61 entries: ", nrow(psu_counts), "\n"))
      if(nrow(psu_counts) > 0) {
        print(psu_counts)
      }
    } else {
      cat(paste0("\n", dataset_name, " - No data found\n"))
    }
  }
  #
  cat("\n=== END MULTIPLE VISITS CHECK ===\n")  
  
  ################################## fill in benthic absences ##################################
  
  ### Check for incomplete absences by sampling event (fewer than 61 species per sampling event)
  # NOTE - everything looks good! only NCRMP had true missing absences in the first place
  # Check for incomplete sampling events (fewer than 61 species per sampling event)
  cat("=== CHECKING FOR INCOMPLETE SAMPLING EVENTS ===\n")
  #
  # For datasets with dates (most datasets)
  datasets_with_dates <- combined_benthic_data %>%
    filter(!dataset %in% c("SESAP", "NODICE")) %>%
    group_by(dataset, PSU, date, lat, lon) %>%
    summarise(
      species_count = n_distinct(spp),
      total_records = n(),
      .groups = "drop"
    ) %>%
    filter(species_count < 61) %>%
    arrange(dataset, species_count)
  
  cat("Sampling events with <61 species (datasets with dates):", nrow(datasets_with_dates), "\n")
  if(nrow(datasets_with_dates) > 0) {
    print(datasets_with_dates)
  }
  #
  # For SESAP and NODICE (no dates, use PSU/lat/lon)
  sesap_nodice <- combined_benthic_data %>%
    filter(dataset %in% c("SESAP", "NODICE")) %>%
    group_by(dataset, PSU, lat, lon) %>%
    summarise(
      species_count = n_distinct(spp),
      total_records = n(),
      .groups = "drop"
    ) %>%
    filter(species_count < 61) %>%
    arrange(dataset, species_count)
  #
  cat("\nSampling events with <61 species (SESAP/NODICE):", nrow(sesap_nodice), "\n")
  if(nrow(sesap_nodice) > 0) {
    print(sesap_nodice)
  }
  
  ### Fill in absences (mainly matters for NCRMP but also doing it for all datasets
  #     since there are different factor levels by dataset)
  # Step 1: Get all unique species from your dataset
  all_species <- unique(combined_benthic_data$spp)
  cat("Total unique species found:", length(all_species), "\n")
  cat("First 10 species:\n")
  print(head(all_species, 10))
  
  # Step 2: Create sampling event identifiers
  # This combines all the grouping variables that define a unique sampling event
  sampling_events <- combined_benthic_data %>%
    select(dataset, PSU, date, lat, lon) %>%
    distinct()
  
  cat("\nTotal unique sampling events:", nrow(sampling_events), "\n")
  
  # Step 3: Create complete grid of all possible combinations
  complete_grid <- sampling_events %>%
    crossing(spp = all_species)
  
  cat("Expected total rows after completion:", nrow(complete_grid), "\n")
  
  # Step 4: Join with original data and fill missing values
  completed_data <- complete_grid %>%
    left_join(combined_benthic_data, 
              by = c("dataset", "PSU", "date", "lat", "lon", "spp")) %>%
    mutate(
      # Create indicator for inferred absence BEFORE filling missing values
      # If cover is NA, it means this record wasn't in the original data
      inferred_absence = ifelse(is.na(cover), "Y", "N"),
      # Fill missing cover values with 0
      cover = ifelse(is.na(cover), 0, cover)
      # Note: depth remains as is - NA for inferred records, original values for real records
    ) %>%
    arrange(dataset, PSU, date, lat, lon, spp)
  
  # Step 5: Summary of changes
  original_rows <- nrow(combined_benthic_data)
  final_rows <- nrow(completed_data)
  inferred_records <- sum(completed_data$inferred_absence == "Y")
  
  cat("\n=== SUMMARY ===\n")
  cat("Original rows:", original_rows, "\n")
  cat("Final rows:", final_rows, "\n")
  cat("Records added (inferred absences):", inferred_records, "\n")
  cat("Percentage of data that is inferred:", 
      round(100 * inferred_records / final_rows, 1), "%\n")
  
  combined_benthic_data = completed_data
  
  ################################## collate demo data ##################################
  
  # Function to standardize demographic/health datasets
  standardize_demo_dataset <- function(df, dataset_name) {
    # Add source dataset column
    df$dataset <- dataset_name
    
    # # Standardize column names based on what's available
    # # Handle different depth/distance column names
    # if ("meterscompleted" %in% names(df)) {
    #   df <- df %>% rename(depth = meterscompleted)
    # }
    
    # Handle different mortality column names
    if ("recmort" %in% names(df)) {
      df <- df %>% rename(recentmort = recmort)
    }
    
    # Ensure all datasets have the same core columns
    # Add missing columns with NA if they don't exist
    required_cols <- c("date", "PSU", "lat", "lon", "depth", "spp", 
                       "length", "width", "height", "oldmort", "recentmort")
    
    for (col in required_cols) {
      if (!col %in% names(df)) {
        df[[col]] <- NA
      }
    }
    
    # Handle any additional columns that might be useful
    # Keep sampledate if it exists (from TCRMP)
    if ("sampledate" %in% names(df)) {
      df <- df %>% select(dataset, all_of(required_cols), sampledate, everything())
    } else {
      df <- df %>% select(dataset, all_of(required_cols), everything())
    }
    
    # Move dataset to leftmost position
    df <- df %>% relocate(dataset, .before = everything())
    
    return(df)
  }
  
  # Standardize each dataset
  TCRMP_health_std <- standardize_demo_dataset(TCRMP_health_long, "TCRMP_health")
  NCRMP_demo_std <- standardize_demo_dataset(NCRMP_demo_long, "NCRMP_demo")
  
  # Combine the datasets
  combined_demo_data <- bind_rows(
    TCRMP_health_std,
    NCRMP_demo_std
  )
  
  # Convert spp to factor
  combined_demo_data$spp <- as.factor(combined_demo_data$spp)
  
  # Show all species levels
  print("\nAll species levels:")
  spp_levels <- levels(combined_demo_data$spp)
  print(spp_levels)
  print(paste("Total number of species:", length(spp_levels)))
  
  
  # STOPPING POINT 6 June 2025- almost at the point of summarizing everything together. will still need
  #   to go back and get things like meters covered maybe...or just confirm. and also get any
  #   datasets I'm still missing
  
  ################################## fill in demo absences ##################################
  
  # NOTE / STOPPING POINT - 9 June 2025: should check for duplication in demo at the END of
  #           the introduction of absence data, to make the dataframe(s) easy to parse
  
  # NOTE - do I actually even need to do this before summarizing? Maybe?
  
  
  
  ################################## summarize cover & SA ##################################
  
  #calculate susceptible tissue surface area (SA)
  # NOTE - could consider trimming down predicted SA for OANN since it has a lot of dead space in reality
  #         - but, it also has more or equal SA potentially in lobes than it would as an ellipsoid?
  #
  # calculate total mortality, then remove corals that had suffered 100% mortality upon survey
  #   NOTE - NAs treated as 0% mortality here
  combined_demo_data_grouped = combined_demo_data %>%
    mutate(totalmort = coalesce(oldmort, 0) + coalesce(recentmort, 0)) %>%
    filter(totalmort < 100)
  #
  # set '0' cm heights as an arbitrarily small amount (cm) to allow the hemi-ellipsoid estimation to function correctly. coral recruits 
  #   have very little tangible height, and were recorded underwater as 0 cm height. there is also
  #   an instance of '-1' height but it is an acroporid and gets filtered out downstream anyways
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    mutate(height = ifelse(height <= 0, 0.01, height))
  #
  # hemi-ellipsoid estimation (Knud Thomsen approximation; see Xu 2009 but also Holstein 2015).
  #   p is a dimensionless constant; all else in square cm
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    mutate(
      a = height,
      b = width / 2,
      c = length / 2,
      p = 1.6075,
      SA_colony = 2 * pi * (((a * b)^p + (a * c)^p + (b * c)^p) / 3)^(1 / p),
      SA_colony = SA_colony / 10000,  # Convert from square cm to square meters
      SA = SA_colony*(1-(totalmort)/100)
    ) %>%
    select(-a, -b, -c, -p)
  #
  # remove unnecessary variables
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    select(-length, -width, -height, -oldmort, -recentmort, -SA_colony)
  
  #filter out corals simply marked as 'juvenile' or 'Coral' - since we don't know what species they were
  combined_benthic_data_grouped = combined_benthic_data
  combined_benthic_data_grouped = combined_benthic_data_grouped %>%
    filter(!spp %in% c('Juvenile coral spp.', 'Coral spp.', 'Coral juvenile', 'Hard Coral, unknown spp.')) %>%
    filter(!grepl('Millepora', spp)) %>% #also drop hydrozoans
    mutate(spp = droplevels(spp))
  
  #filter 'Montastraea spp' since none were marked as present anyways, and this taxonomy is ambiguous
  combined_benthic_data_grouped = combined_benthic_data_grouped %>%
    filter(!spp %in% c('Montastraea species', 'Montastraea spp.')) %>%
    mutate(spp = droplevels(spp))
  
  #break corals into susceptibility groups
  # NOTE - made some assumptions here; please see 'questionable' groupings below. a main assumption was
  #         that corals in Mussidae are all considered highly susceptible, even if they are so small/rare
  #         that we can't confirm this from the field
  #
  # NOTE - assuming Oculina diffusa (1 occurrence) and Tubastraea coccinea (also 1 colony) are 'Unaffected'
  #
  # NOTE - Questionable corals:
  # - Pseudodiploria clivosa (see Spadafore 2021 for some more reference on this)
  # - Scolymia
  # - Agaricia
  # - Madracis
  # - Helioseris
  # - Manicina
  # - Siderastrea radians
  # - Favia fragum
  # - Isophyllia
  # - Solenastrea (see Spadafore 2021 for some more reference on this)
  # - Tubastrea coccinea
  # - Oculina
  # - Porites
  # - Mycetophyllia
  # - Mussa
  #
  # NOTE - how do we deal with situations like 'Montastraea spp' absences ? pre 2012ish or so, this included
  #         orbicellids, which we now know should not be the case. important distinction...may need to
  #         simply toss that species level entirely. it is all absences, at least
  #
  levels(combined_benthic_data_grouped$spp)
  combined_benthic_data_grouped = combined_benthic_data_grouped %>%
    mutate(
      susc = case_when(
        spp %in% c('Agaricia agaricites', 'Agaricia fragilis', 'Agaricia grahamae', 'Agaricia humilis',
                   'Agaricia lamarcki', 'Agaricia species', 'Agaricia spp', 'Agaricia spp.',
                   'Agaricia tenuifolia', 'Agaricia undata', 'Branching Porites spp.', 'Madracis auretenra',
                    'Madracis decactis', 'Madracis formosa', 'Madracis mirabilis', 'Madracis pharensis',
                   'Madracis spp', 'Madracis spp.', 'Porites astreoides', 'Porites branching species',
                   'Porites branneri', 'Porites colonensis', 'Porites divaricata', 'Porites furcata',
                   'Porites porites', 'Porites spp', 'Siderastrea radians', 'Siderastrea siderea',
                   'Siderastrea species', 'Siderastrea spp', 'Siderastrea spp.',
                   'Stephanocoenia intercepta', 'Stephanocoenia intersepta',
                   'Helioceris cucullata', 'Helioseris cucullata', 'Leptoseris cucullata') ~ 'low',
        spp %in% c('Montastraea annularis', 'Montastraea annularis complex', 'Montastraea cavernosa',
                   'Montastraea faveolata', 'Montastraea franksi', 'Montastraea species',
                   'Montastraea spp.', 'Orbicella annularis', 'Orbicella annularis species complex',
                   'Orbicella faveolata', 'Orbicella franksi', 'Orbicella franksii',
                   'Orbicella species complex', 'Orbicella spp', 'Solenastrea bournoni',
                   'Solenastrea hyades', 'Solenastrea spp') ~ 'moderate',
        spp %in% c('Colpophyllia natans', 'Dendrogyra cylindrus', 'Dichocoenia stokesii',
                   'Diploria labyrinthiformis', 'Eusmilia fastigiata', 'Meandrina danae',
                   'Meandrina jacksoni', 'Meandrina meandrites', 'Meandrina spp', 'Mycetophyllia aliciae',
                   'Mycetophyllia danaana', 'Mycetophyllia daniana', 'Mycetophyllia ferox',
                   'Mycetophyllia lamarckiana', 'Mycetophyllia reesi', 'Mycetophyllia species',
                   'Mycetophyllia spp.', 'Pseudodiploria clivosa', 'Pseudodiploria spp', 'Diploria strigosa',
                   'Diploria clivosa', 'Pseudodiploria strigosa', 'Isophyllastrea rigida',
                   'Isopyhyllastrea rigida', 'Isophyllia sinuosa', 'Manicina areolata',
                   'Mussa angulosa', 'Scolymia cubensis', 'Scolymia lacera', 'Scolymia species', 'Scolymia spp',
                   'Scolymia spp.', 'Manicina areolata', 'Favia fragum') ~ 'high',
        spp %in% c('Acropora cervicornis', 'Acropora palmata', 'Acropora prolifera',
                   'Oculina diffusa', 'Tubastraea coccinea') ~ 'Unaffected'
      )
    ) %>%
    mutate(susc = as.factor(susc))
  
  # Update species names to current taxonomy and correct spelling errors,and also collapse species to
  #   genus where appropriate
  #     NOTE - the largest effect of collapsing here is on Agaricia, because there is so much A. undata
  #             at mesophotic depth. should carefully consider effect this has on distribution modeling
  #          - but also a large effect on Porites, since there is now no distinction between branching and
  #             P. astreoides
  #          - effect on Madracis may be important because of certain species being prevalent deep
  #     NOTE - am treating orbicella complex as its own species! a tough choice, but might make the most sense
  #             given available data. importantly, am including the few 'Orbicella spp' entries from NCRMP
  #             in this category, and also merging ALL of O. franksi & O. faveolata into it, since these
  #             are commonly extremely hard to distinguish or may be mostly franksi anyways
  combined_benthic_data_grouped = combined_benthic_data_grouped %>%
    mutate(
      spp = case_when(
        spp == 'Isophyllastrea rigida' ~ 'Isophyllia rigida',
        spp == 'Montastraea annularis' ~ 'Orbicella annularis',
        spp == 'Montastraea annularis' ~ 'Orbicella annularis',
        spp == 'Montastraea annularis complex' ~ 'Orbicella annularis',
        spp == 'Orbicella annularis species complex' ~ 'Orbicella annularis',
        spp == 'Montastraea faveolata' ~ 'Orbicella comp.',
        spp == 'Montastraea franksi' ~ 'Orbicella comp.',
        spp == 'Orbicella species complex' ~ 'Orbicella comp.',
        spp == 'Orbicella spp' ~ 'Orbicella comp.',
        spp == 'Orbicella faveolata' ~ 'Orbicella comp.',
        spp == 'Orbicella franksi' ~ 'Orbicella comp.',
        spp == 'Orbicella franksii' ~ 'Orbicella comp.',
        spp == 'Diploria clivosa' ~ 'Pseudodiploria clivosa',
        spp == 'Diploria strigosa' ~ 'Pseudodiploria strigosa',
        spp == 'Helioceris cucullata' ~ 'Helioseris cucullata',
        spp == 'Leptoseris cucullata' ~ 'Helioseris cucullata',
        spp == 'Isopyhyllastrea rigida' ~ 'Isophyllia rigida',
        spp == 'Stephanocoenia intercepta' ~ 'Stephanocoenia intersepta',
        TRUE ~ spp  # keep all other species names unchanged
      )
    ) %>%
    mutate(spp = factor(spp)) %>%  # convert back to factor
    mutate(spp = droplevels(spp))  # drop unused factor levels
  levels(combined_benthic_data_grouped$spp)
  
  # collapse species to genus where appropriate
  #     NOTE - the largest effect of collapsing here is on Agaricia, because there is so much A. undata
  #             at mesophotic depth. should carefully consider effect this has on distribution modeling
  #          - but also a large effect on Porites, since there is now no distinction between branching and
  #             P. astreoides
  #          - effect on Madracis may be important because of certain species being prevalent deep
  combined_benthic_data_grouped = combined_benthic_data_grouped %>%
    mutate(
      spp = case_when(
        grepl('Agaricia', spp) ~ 'Agaricia spp',
        grepl('Porites', spp) ~ 'Porites spp',
        grepl('Madracis', spp) ~ 'Madracis spp',
        grepl('Meandrina', spp) ~ 'Meandrina spp',
        grepl('Mycetophyllia', spp) ~ 'Mycetophyllia spp',
        grepl('Pseudodiploria', spp) ~ 'Pseudodiploria spp',
        grepl('Scolymia', spp) ~ 'Scolymia spp',
        grepl('Siderastrea', spp) ~ 'Siderastrea spp',
        grepl('Solenastrea', spp) ~ 'Solenastrea spp',
        TRUE ~ spp  # keep all other species names unchanged
      )
    ) %>%
    mutate(spp = factor(spp)) %>%  # convert back to factor
    mutate(spp = droplevels(spp))  # drop unused factor levels
  levels(combined_benthic_data_grouped$spp)
  
  #drop the unaffected species
  #   - Acropora: like 35 occurrences total
  #   - Oculina & Tubastraea: 1 occurrence each
  combined_benthic_data_grouped = combined_benthic_data_grouped %>%
    filter(!susc %in% c('Unaffected')) %>%
    mutate(susc = droplevels(susc))
  
  #quick read-out of species by susceptibility group
  combined_benthic_data_grouped %>%
    distinct(spp, susc) %>%
    arrange(susc, spp) %>%
    group_by(susc) %>%
    summarise(species = paste(spp, collapse = "\n  - ")) %>%
    mutate(output = paste0(susc, ":\n  - ", species)) %>%
    pull(output) %>%
    cat(sep = "\n\n")
  
  #filter out unidentified scleractinians from demo data
  # NOTE - these were just 2 very small colonies in 2013. but of course, also there were plenty of "absences"
  #         of unidentifiable coral...just not sure that has real ecological meaning
  levels(combined_demo_data_grouped$spp)
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    filter(!spp %in% c('Scleractinia spp')) %>%
    filter(!grepl('Millepora', spp)) %>% #also drop hydrozoans
    filter(!grepl('Isophyllia spp', spp)) %>% #drop since no occurences and makes groupings ambiguous
    filter(!grepl('Oculina spp', spp)) %>% #drop since no occurences and makes groupings ambiguous
    mutate(spp = droplevels(spp)) #drop factor levels which no longer are associated with any data
  
  #break corals into susceptibility groups
  # NOTE - same as grouping for coral cover (benthic data)
  levels(combined_demo_data_grouped$spp)
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    mutate(
      susc = case_when(
        spp %in% c('Agaricia agaricites', 'Agaricia fragilis', 'Agaricia grahamae', 'Agaricia humilis',
                   'Agaricia lamarcki', 'Agaricia species', 'Agaricia spp', 'Agaricia spp.',
                   'Agaricia tenuifolia', 'Agaricia undata', 'Branching Porites spp.', 'Madracis auretenra',
                   'Madracis decactis', 'Madracis formosa', 'Madracis mirabilis', 'Madracis pharensis',
                   'Madracis spp', 'Madracis spp.', 'Porites astreoides', 'Porites branching species',
                   'Porites branneri', 'Porites colonensis', 'Porites divaricata', 'Porites furcata',
                   'Porites porites', 'Porites spp', 'Siderastrea radians', 'Siderastrea siderea',
                   'Siderastrea species', 'Siderastrea spp', 'Siderastrea spp.',
                   'Stephanocoenia intercepta', 'Stephanocoenia intersepta',
                   'Helioceris cucullata', 'Helioseris cucullata', 'Leptoseris cucullata') ~ 'low',
        spp %in% c('Montastraea annularis', 'Montastraea annularis complex', 'Montastraea cavernosa',
                   'Montastraea faveolata', 'Montastraea franksi', 'Montastraea species',
                   'Montastraea spp.', 'Orbicella annularis', 'Orbicella annularis species complex',
                   'Orbicella faveolata', 'Orbicella franksi', 'Orbicella franksii',
                   'Orbicella species complex', 'Orbicella spp', 'Solenastrea bournoni',
                   'Solenastrea hyades', 'Solenastrea spp') ~ 'moderate',
        spp %in% c('Colpophyllia natans', 'Dendrogyra cylindrus', 'Dichocoenia stokesii',
                   'Diploria labyrinthiformis', 'Eusmilia fastigiata', 'Meandrina danae',
                   'Meandrina jacksoni', 'Meandrina meandrites', 'Meandrina spp', 'Mycetophyllia aliciae',
                   'Mycetophyllia danaana', 'Mycetophyllia daniana', 'Mycetophyllia ferox',
                   'Mycetophyllia lamarckiana', 'Mycetophyllia reesi', 'Mycetophyllia species',
                   'Mycetophyllia spp.', 'Mycetophyllia spp', 'Pseudodiploria clivosa', 'Pseudodiploria spp', 'Diploria strigosa',
                   'Diploria clivosa', 'Pseudodiploria strigosa', 'Isophyllastrea rigida',
                   'Isopyhyllastrea rigida', 'Isophyllia sinuosa', 'Manicina areolata',
                   'Mussa angulosa', 'Scolymia cubensis', 'Scolymia lacera', 'Scolymia species', 'Scolymia spp',
                   'Scolymia spp.', 'Manicina areolata', 'Favia fragum') ~ 'high',
        spp %in% c('Acropora cervicornis', 'Acropora palmata', 'Acropora prolifera',
                   'Oculina diffusa', 'Tubastraea coccinea') ~ 'Unaffected'
      )
    ) %>%
    mutate(susc = as.factor(susc))
  
  # Update species names to current taxonomy and correct spelling errors, and also collapse species to
  #   genus where appropriate
  #     NOTE - the largest effect of collapsing here is on Agaricia, because there is so much A. undata
  #             at mesophotic depth. should carefully consider effect this has on distribution modeling
  #          - but also a large effect on Porites, since there is now no distinction between branching and
  #             P. astreoides
  #          - effect on Madracis may be important because of certain species being prevalent deep
  #     NOTE - am treating orbicella complex as its own species! a tough choice, but might make the most sense
  #             given available data. importantly, am including the few 'Orbicella spp' entries from NCRMP
  #             in this category, and also merging ALL of O. franksi & O. faveolata into it, since these
  #             are commonly extremely hard to distinguish or may be mostly franksi anyways
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    mutate(
      spp = case_when(
        spp == 'Isophyllastrea rigida' ~ 'Isophyllia rigida',
        spp == 'Montastraea annularis' ~ 'Orbicella annularis',
        spp == 'Montastraea annularis' ~ 'Orbicella annularis',
        spp == 'Montastraea annularis complex' ~ 'Orbicella annularis',
        spp == 'Orbicella annularis species complex' ~ 'Orbicella annularis',
        spp == 'Montastraea faveolata' ~ 'Orbicella comp.',
        spp == 'Montastraea franksi' ~ 'Orbicella comp.',
        spp == 'Orbicella species complex' ~ 'Orbicella comp.',
        spp == 'Orbicella spp' ~ 'Orbicella comp.',
        spp == 'Orbicella faveolata' ~ 'Orbicella comp.',
        spp == 'Orbicella franksi' ~ 'Orbicella comp.',
        spp == 'Orbicella franksii' ~ 'Orbicella comp.',
        spp == 'Diploria clivosa' ~ 'Pseudodiploria clivosa',
        spp == 'Diploria strigosa' ~ 'Pseudodiploria strigosa',
        spp == 'Helioceris cucullata' ~ 'Helioseris cucullata',
        spp == 'Leptoseris cucullata' ~ 'Helioseris cucullata',
        spp == 'Isopyhyllastrea rigida' ~ 'Isophyllia rigida',
        spp == 'Stephanocoenia intercepta' ~ 'Stephanocoenia intersepta',
        TRUE ~ spp  # keep all other species names unchanged
      )
    ) %>%
    mutate(spp = factor(spp)) %>%  # convert back to factor
    mutate(spp = droplevels(spp))  # drop unused factor levels
  levels(combined_demo_data_grouped$spp)
  
  # collapse species to genus where appropriate
  #     NOTE - the largest effect of collapsing here is on Agaricia, because there is so much A. undata
  #             at mesophotic depth. should carefully consider effect this has on distribution modeling
  #          - but also a large effect on Porites, since there is now no distinction between branching and
  #             P. astreoides
  #          - effect on Madracis may be important because of certain species being prevalent deep
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    mutate(
      spp = case_when(
        grepl('Agaricia', spp) ~ 'Agaricia spp',
        grepl('Porites', spp) ~ 'Porites spp',
        grepl('Madracis', spp) ~ 'Madracis spp',
        grepl('Meandrina', spp) ~ 'Meandrina spp',
        grepl('Mycetophyllia', spp) ~ 'Mycetophyllia spp',
        grepl('Pseudodiploria', spp) ~ 'Pseudodiploria spp',
        grepl('Scolymia', spp) ~ 'Scolymia spp',
        grepl('Siderastrea', spp) ~ 'Siderastrea spp',
        grepl('Solenastrea', spp) ~ 'Solenastrea spp',
        TRUE ~ spp  # keep all other species names unchanged
      )
    ) %>%
    mutate(spp = factor(spp)) %>%  # convert back to factor
    mutate(spp = droplevels(spp))  # drop unused factor levels
  levels(combined_demo_data_grouped$spp)
  
  #drop the unaffected species
  #   - Acropora: like 35 occurrences total
  #   - Oculina & Tubastraea: 1 occurrence each
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    filter(!susc %in% c('Unaffected')) %>%
    mutate(susc = droplevels(susc))
  
  #quick read-out of species by susceptibility group
  combined_demo_data_grouped %>%
    distinct(spp, susc) %>%
    arrange(susc, spp) %>%
    group_by(susc) %>%
    summarise(species = paste(spp, collapse = "\n  - ")) %>%
    mutate(output = paste0(susc, ":\n  - ", species)) %>%
    pull(output) %>%
    cat(sep = "\n\n")

  #filter only 2013 to November 2018 pre-SCTLD
  # NOTE - consider if this is the filter we want!!
  combined_benthic_data_grouped = combined_benthic_data_grouped %>%
    filter(year(date) >= 2013 & date <= as.Date("2018-11-30"))
  combined_demo_data_grouped = combined_demo_data_grouped %>%
    filter(year(date) >= 2013 & date <= as.Date("2018-11-30"))
  
  # NOTE / STOPPING POINT - 12 June 2025
  # - there is something very strange in the data. the number of observations differs wildly between species within 
  #     a PSU, when it really should be the exact same number. something funky going on upstream
  # - this also really messes with the averages, since there are too many or a random number of zeros. need to fix!
  # - one consideration would simply be to "fill in" the zeros (absences) only AFTER finalizing the species factor
  #     levels. this way, it is easy to confirm all data are in the proper format for "continuous/absence"
  
  #summaries
  # cover_by_susc <- combined_benthic_data_grouped %>%
  #   group_by(PSU, susc) %>%
  #   summarise(total_cover = sum(cover, na.rm = TRUE), .groups = 'drop')
  # #
  # total_cover_by_PSU <- combined_benthic_data_grouped %>%
  #   group_by(PSU) %>%
  #   summarise(total_cover = sum(cover, na.rm = TRUE), .groups = 'drop')
  # Average cover for each species within each PSU across dates, with sample size
  avg_cover_by_species_PSU <- combined_benthic_data_grouped %>%
    group_by(PSU, spp) %>%
    summarise(
      avg_cover = mean(cover, na.rm = TRUE),
      N = n(),  # Total number of transect observations
      .groups = 'drop'
    )
  
  # Then sum across species to get total average cover per PSU
  total_avg_cover_by_PSU <- avg_cover_by_species_PSU %>%
    group_by(PSU) %>%
    summarise(
      total_avg_cover = sum(avg_cover, na.rm = TRUE),
      .groups = 'drop'
    )
  #
  SA_by_susc <- combined_demo_data_grouped %>%
    group_by(PSU, susc) %>%
    summarise(total_SA = sum(sa, na.rm = TRUE), .groups = 'drop')
  #
  total_SA_by_PSU <- combined_demo_data_grouped %>%
    group_by(PSU) %>%
    summarise(total_SA = sum(sa, na.rm = TRUE), .groups = 'drop')
  
  
  # STOPPING POINT - 2 June 2025. am about to produce plots of SA, then should figure out 'absence' for
  #   demo surveys, then will collate the demo & LPI surveys together, then append bathymetry & derived
  #   values to dataframe, then produce a simple test GAM. finally will need to figure out prediction
  #   onto a uniform grid which will "talk" to the habitat grid plugged into the CMS
  
  #preserve absence data
  all_PSU_info_cover <- NCRMP_benthic_long %>%
    distinct(PSU, lat, lon, date)
  
  # Updated summaries that include zero-cover PSU's
  cover_by_susc_complete <- NCRMP_benthic_long %>%
    group_by(PSU, susc) %>%
    summarise(total_cover = sum(cover, na.rm = TRUE), .groups = 'drop') %>%
    complete(PSU = all_PSU_info_cover$PSU, 
             susc = c('low', 'moderate', 'high', 'Unaffected', 'Unknown'), 
             fill = list(total_cover = 0)) %>%
    left_join(all_PSU_info_cover, by = "PSU")  # Add back location info
  
  total_cover_by_PSU_complete <- NCRMP_benthic_long %>%
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