  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(tidyverse)
  library(readxl)
  library(terra)
  
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
  load(here('data/PRICO_2014_benthic_cover.rda'))
  load(here('data/PRICO_2016_benthic_cover.rda'))
  load(here('data/PRICO_2014_coral_demographics.rda'))
  load(here('data/PRICO_2016_coral_demographics.rda'))
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
  # Find the Transect column position and force it to text
  # You'll need to check which column number Transect is - let's assume it's column 3
  col_types[9] <- "text"         # Force Transect to be text to handle mixed values
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
  
  #there were some issues with the read-in date format - just setting it to 10/23/2015 manually
  #   - some filming dates were actually 10/22/2015, but doesn't matter for this analysis
  NODICE_long = NODICE_long %>%
    mutate(date = as.POSIXct("2015-10-23"))
  
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
  # NOTE - can remove transect number here, since it was always '1'
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
  
  #as with NODICE, there were some issues with the read-in date format - just setting it to 2013-08-22 manually
  #   - some filming dates were actually out to at least 2013-10-04, but doesn't matter for this analysis
  SESAP_long = SESAP_long %>%
    mutate(date = as.POSIXct("2013-08-22"))
  
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
  #filter out coral recruitment surveys - these only looked at very small juvenile corals,
  #   and are not suitable for averaging with regular surveys
  TCRMP_health = TCRMP_health %>%
    filter(Period != 'Juvenile')
  
  # #drop all instances where there is no transect information. this problem is specific to Black Point, 
  # #   December 2017 where for whatever reason there was a very large amount of row entries with no transect
  # #   number noted. very odd; makes it harder to clearly process data downstream so filtering out here
  # # NOTE - consider if this should have been done!
  # TCRMP_health = TCRMP_health %>%
  #   filter(!is.na(Transect))
  
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
    USVI_2013_benthic_cover %>% mutate(region_id = "VI"),
    USVI_2015_benthic_cover %>% mutate(region_id = "VI"),
    USVI_2017_benthic_cover %>% mutate(region_id = "VI"),
    PRICO_2014_benthic_cover %>% mutate(region_id = "PR"),
    PRICO_2016_benthic_cover %>% mutate(region_id = "PR")
  )
  
  #refactor variables and remove unnecessary ones
  NCRMP_benthic_long = NCRMP_benthic %>%
    mutate(
      PRIMARY_SAMPLE_UNIT = paste(PRIMARY_SAMPLE_UNIT, region_id, sep = "_"),
      PRIMARY_SAMPLE_UNIT = as.factor(PRIMARY_SAMPLE_UNIT)
    ) %>%
    select(-REGION, -STATION_NR, -RUGOSITY_CD, -WTD_RUG, -MAPGRID_NR, -HABITAT_CD, -STRAT, -SUB_REGION_NAME, -SUB_REGION_NR,
           -ZONE_NAME, -ZONE_NR, -MPA_NAME, -MPA_NR, -ADMIN, -PROT, -DEPTH_STRAT, -region_id)

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
  
  #clip to extent of Chapter 3 hydrodynamic modeling
  # NOTE - could instead choose to use ALL of Puerto Rico data - but not sure this makes sense for our work
  #
  #read the CMS-formatted hydrodynamic domain extent
  hydro_extent = vect(here("output/", "hydro_domain_extent.shp"))
  #
  # Create a spatial points object from the NCRMP data
  ncrmp_points_proj <- vect(NCRMP_benthic_long, 
                       geom = c("lon", "lat"), 
                       crs = "EPSG:4269")  # NAD83
  #
  # Convert hydro_extent from 0-360° to -180-180° format before projecting
  # First convert back to geographic coordinates if it's projected
  hydro_extent_geo <- project(hydro_extent, "EPSG:4269")
  #
  # Fix the longitude coordinates in the hydro extent
  hydro_coords <- crds(hydro_extent_geo)
  hydro_coords[, 1] <- ifelse(hydro_coords[, 1] > 180, 
                              hydro_coords[, 1] - 360, 
                              hydro_coords[, 1])
  #
  # Create a new hydro extent with corrected coordinates
  hydro_extent_corrected <- vect(hydro_coords, type = "polygons", crs = "EPSG:4269")
  #
  # Find which points are within the hydro extent
  points_within <- relate(ncrmp_points_proj, hydro_extent_corrected, "within")
  #
  # Filter the original dataframe to keep only rows within the extent
  NCRMP_benthic_long <- NCRMP_benthic_long[points_within, ]
  #
  # Plot the results
  plot(hydro_extent_corrected, main = "NCRMP Points Within Hydro Extent", 
       col = "lightblue", border = "blue", lwd = 2)
  #
  # Add all original points in gray
  plot(ncrmp_points_proj, add = TRUE, col = "gray", pch = 16, cex = 0.5)
  #
  # Add points within extent in red
  ncrmp_within_proj <- ncrmp_points_proj[points_within]
  plot(ncrmp_within_proj, add = TRUE, col = "red", pch = 16, cex = 0.8)
  #
  # Add legend
  legend("topright", 
         legend = c("Hydro Extent", "All NCRMP Points", "Points Within Extent"),
         col = c("lightblue", "gray", "red"),
         pch = c(15, 16, 16),
         cex = 0.8)  
  test <- vect(NCRMP_benthic_long,
                            geom = c("lon", "lat"),
                            crs = "EPSG:4269")  # NAD83
  plot(test)
  
  #drop species code
  NCRMP_benthic_long_withcode = NCRMP_benthic_long
  NCRMP_benthic_long = NCRMP_benthic_long %>%
    select(-code)
  
  #demo
  #
  NCRMP_demo = bind_rows(
    USVI_2013_coral_demographics %>% mutate(region_id = 'VI'),
    USVI_2015_coral_demographics %>% mutate(region_id = 'VI'),
    USVI_2017_coral_demographics %>% mutate(region_id = 'VI'),
    PRICO_2014_coral_demographics %>% mutate(region_id = 'PR'),
    PRICO_2016_coral_demographics %>% mutate(region_id = 'PR')
  )
  
  #refactor variables and remove unnecessary ones
  NCRMP_demo_long = NCRMP_demo %>%
    mutate(
      PRIMARY_SAMPLE_UNIT = paste(PRIMARY_SAMPLE_UNIT, region_id, sep = "_"),
      PRIMARY_SAMPLE_UNIT = as.factor(PRIMARY_SAMPLE_UNIT)
    ) %>%
    select(-REGION, -STATION_NR, -RUGOSITY_CD, -WTD_RUG, -MAPGRID_NR, -HABITAT_CD, -STRAT, -SUB_REGION_NAME, -SUB_REGION_NR,
           -ZONE_NAME, -ZONE_NR, -MPA_NAME, -MPA_NR, -ADMIN, -PROT, -DEPTH_STRAT, -N, -JUV,
           -BLEACH_CONDITION, -DISEASE, -region_id)
  
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
  
  #clip to extent of Chapter 3 hydrodynamic modeling
  #
  # Create a spatial points object from the NCRMP data
  ncrmp_points_proj <- vect(NCRMP_demo_long, 
                            geom = c("lon", "lat"), 
                            crs = "EPSG:4269")  # NAD83
  #
  # Convert hydro_extent from 0-360° to -180-180° format before projecting
  # First convert back to geographic coordinates if it's projected
  hydro_extent_geo <- project(hydro_extent, "EPSG:4269")
  #
  # Find which points are within the hydro extent
  points_within <- relate(ncrmp_points_proj, hydro_extent_corrected, "within")
  #
  # Filter the original dataframe to keep only rows within the extent
  NCRMP_demo_long <- NCRMP_demo_long[points_within, ]
  #
  # Plot the results
  plot(hydro_extent_corrected, main = "NCRMP Points Within Hydro Extent", 
       col = "lightblue", border = "blue", lwd = 2)
  #
  # Add all original points in gray
  plot(ncrmp_points_proj, add = TRUE, col = "gray", pch = 16, cex = 0.5)
  #
  # Add points within extent in red
  ncrmp_within_proj <- ncrmp_points_proj[points_within]
  plot(ncrmp_within_proj, add = TRUE, col = "red", pch = 16, cex = 0.8)
  #
  # Add legend
  legend("topright", 
         legend = c("Hydro Extent", "All NCRMP Points", "Points Within Extent"),
         col = c("lightblue", "gray", "red"),
         pch = c(15, 16, 16),
         cex = 0.8)  
  test <- vect(NCRMP_demo_long,
                            geom = c("lon", "lat"),
                            crs = "EPSG:4269")  # NAD83
  plot(test)
  
  #drop species code
  NCRMP_demo_long_withcode = NCRMP_demo_long
  NCRMP_demo_long = NCRMP_demo_long %>%
    select(-code)
  
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
    
    # If transect column doesn't exist, set to NA
    if (!"transect" %in% names(df)) {
      df$transect <- NA
    }
    
    # Select and reorder columns consistently
    df <- df %>%
      select(dataset, date, PSU, transect, lat, lon, cover, spp, depth) %>%
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
  
  ################################## check incomplete benthic absences ##################################
  
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
    required_cols <- c("date", "PSU", "transect", "lat", "lon", "depth", "spp", 
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
  
  ################################## define susceptibility ##################################
  
  ## BENTHIC
  #
  #filter out corals simply marked as 'juvenile' or 'Coral' - since we don't know what species they were
  combined_benthic_data_trimmed = combined_benthic_data
  combined_benthic_data_trimmed = combined_benthic_data_trimmed %>%
    filter(!spp %in% c('Juvenile coral spp.', 'Coral spp.', 'Coral juvenile', 'Hard Coral, unknown spp.',
                       'Scleractinia spp')) %>%
    filter(!grepl('Millepora', spp)) %>% #also drop hydrozoans
    mutate(spp = droplevels(spp))
  
  #filter only 2013 to November 2018 pre-SCTLD
  # NOTE - consider if this is the filter we want!!
  combined_benthic_data_trimmed = combined_benthic_data_trimmed %>%
    filter(year(date) >= 2013 & date <= as.Date("2018-11-30"))
  
  #filter 'Montastraea spp' since none were marked as present anyways, and this taxonomy is ambiguous
  combined_benthic_data_trimmed = combined_benthic_data_trimmed %>%
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
  # - Cladocora
  # - Porites
  # - Mycetophyllia
  # - Mussa
  #
  # NOTE - how do we deal with situations like 'Montastraea spp' absences ? pre 2012ish or so, this included
  #         orbicellids, which we now know should not be the case. important distinction...may need to
  #         simply toss that species level entirely. it is all absences, at least
  #
  levels(combined_benthic_data_trimmed$spp)
  combined_benthic_data_trimmed = combined_benthic_data_trimmed %>%
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
                   'Oculina diffusa', 'Tubastraea coccinea', 'Cladocora abruscula') ~ 'Unaffected'
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
  combined_benthic_data_trimmed = combined_benthic_data_trimmed %>%
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
  levels(combined_benthic_data_trimmed$spp)
  
  # collapse species to genus where appropriate
  #     NOTE - the largest effect of collapsing here is on Agaricia, because there is so much A. undata
  #             at mesophotic depth. should carefully consider effect this has on distribution modeling
  #          - but also a large effect on Porites, since there is now no distinction between branching and
  #             P. astreoides
  #          - effect on Madracis may be important because of certain species being prevalent deep
  combined_benthic_data_trimmed = combined_benthic_data_trimmed %>%
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
  levels(combined_benthic_data_trimmed$spp)
  
  #drop the unaffected species
  #   - Acropora: like 35 occurrences total
  #   - Oculina & Tubastraea: 1 occurrence each
  combined_benthic_data_trimmed = combined_benthic_data_trimmed %>%
    filter(!susc %in% c('Unaffected')) %>%
    mutate(susc = droplevels(susc)) %>%
    mutate(spp = droplevels(spp))
  
  #quick read-out of species by susceptibility group
  combined_benthic_data_trimmed %>%
    distinct(spp, susc) %>%
    arrange(susc, spp) %>%
    group_by(susc) %>%
    summarise(species = paste(spp, collapse = "\n  - ")) %>%
    mutate(output = paste0(susc, ":\n  - ", species)) %>%
    pull(output) %>%
    cat(sep = "\n\n")
  
  ## DEMO
  #
  # NOTE / UPDATE
  #   - the more I think about it, the more I think the TCRMP intercept data can't be used at all
  #   - the problem is, all the small corals will be totally lost with this approach, and there is no
  #       way to properly compare with belt transect methodology
  #   - that said, I should still be able to use the 10, 50, and/or 100 cm belt transect data!
  #       - a consideration here though is that belts were only used in situations with already-low density.
  #         so, may bias a bit towards low surface area estimates. still seems like a better solution
  #
  
  #filter only 2013 to November 2018 pre-SCTLD
  combined_demo_data = combined_demo_data %>%
    filter(year(date) >= 2013 & date <= as.Date("2018-11-30"))
  
  #drop all TCRMP health data collected using transect-intercept methods
  combined_demo_data <- combined_demo_data %>%
    filter(method != "intercept" | is.na(method))
  
  #calculate 2D area of survey (disregards rugosity)
  # NOTE - 13 June 2025
  #   - This is where I will go back and fill in better information once I have it. NCRMP is fine,
  #       it is just TCRMP where there is ambiguity
  combined_demo_data <- combined_demo_data %>%
    mutate(
      surveyarea = case_when(
        # If meterscompleted has a value, use it (NCRMP surveys that were all 1 meter wide)
        !is.na(meterscompleted) ~ meterscompleted,

        # # If method is 'intercept', surveyarea is 10 (TCRMP surveys; always went to 10 m length)
        # method == "intercept" ~ 10,

        # If method is '10 cm belt', area is 10m * 0.1m = 1 square meter
        method == "10 cm belt" ~ (10 * 0.1),

        # If method is '50 cm belt', area is 10m * 0.5m = 5 square meters
        method == "50 cm belt" ~ (10 * 0.5),

        # If method is '100 cm belt', area is 10m * 1.0m = 10 square meters
        method == "100 cm belt" ~ (10 * 1),

        # Default case (if none of the above conditions are met)
        TRUE ~ NA_real_
      )
    )
  
  #drop 10-cm belt data. it looks like it inflates surface area density within a transect too much,
  #   unfortunately
  combined_demo_data = combined_demo_data %>%
    filter(method != '10 cm belt' | is.na(method))
  
  #calculate susceptible tissue surface area (SA)
  # NOTE - could consider trimming down predicted SA for OANN since it has a lot of dead space in reality
  #         - but, it also has more or equal SA potentially in lobes than it would as an ellipsoid?
  #
  # calculate total mortality, then remove corals that had suffered 100% mortality upon survey
  #   NOTE - NAs treated as 0% mortality here
  combined_demo_data_trimmed = combined_demo_data %>%
    mutate(totalmort = coalesce(oldmort, 0) + coalesce(recentmort, 0)) %>%
    filter(totalmort < 100)
  #
  # set '0' cm heights as an arbitrarily small amount (cm) to allow the hemi-ellipsoid estimation to function correctly. coral recruits 
  #   have very little tangible height, and were recorded underwater as 0 cm height. there is also
  #   an instance of '-1' height but it is an acroporid and gets filtered out downstream anyways
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
    mutate(height = ifelse(height <= 0, 0.01, height))
  #
  # hemi-ellipsoid estimation (Knud Thomsen approximation; see Xu 2009 but also Holstein 2015).
  #   p is a dimensionless constant; all else in square cm
  #     NOTE / IMPORTANT - at least at some point during TCRMP, colonies <10 cm were not assessed for mortality
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
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
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
    select(-length, -width, -height, -oldmort, -recentmort, -SA_colony)
  
  #filter out unidentified scleractinians from demo data
  # NOTE - these were just 2 very small colonies in 2013. but of course, also there were plenty of "absences"
  #         of unidentifiable coral...just not sure that has real ecological meaning
  levels(combined_demo_data_trimmed$spp)
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
    filter(!spp %in% c('Scleractinia spp', 'Other coral')) %>%
    filter(!grepl('Millepora', spp)) %>% #also drop hydrozoans
    filter(!grepl('Isophyllia spp|Oculina spp', spp)) %>% #drop since no occurrences and makes groupings ambiguous
    mutate(spp = droplevels(spp)) #drop factor levels which no longer are associated with any data
  
  #break corals into susceptibility groups
  # NOTE - same as grouping for coral cover (benthic data)
  levels(combined_demo_data_trimmed$spp)
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
    mutate(
      susc = case_when(
        grepl("Madracis", spp) | #handle strange string for Madracis
        spp %in% c('Agaricia agaricites', 'Agaricia fragilis', 'Agaricia grahamae', 'Agaricia humilis',
                   'Agaricia lamarcki', 'Agaricia species', 'Agaricia spp', 'Agaricia spp.',
                   'Agaricia tenuifolia', 'Agaricia undata', 'Undaria spp', 'Branching Porites spp.',
                   'Madracis auretenra', 'Madracis decactis', 'Madracis formosa', 'Madracis mirabilis',
                   'Madracis pharensis', 'Madracis spp', 'Madracis spp.', 'Porites astreoides',
                   'Porites branching species', 'Porites branneri', 'Porites colonensis', 'Porites divaricata',
                   'Porites furcata', 'Porites porites', 'Porites spp', 'Siderastrea radians',
                   'Siderastrea siderea', 'Siderastrea species', 'Siderastrea spp', 'Siderastrea spp.',
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
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
    mutate(
      spp = case_when(
        spp == 'Undaria spp' ~ 'Agaricia spp',
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
  levels(combined_demo_data_trimmed$spp)
  
  # collapse species to genus where appropriate
  #     NOTE - the largest effect of collapsing here is on Agaricia, because there is so much A. undata
  #             at mesophotic depth. should carefully consider effect this has on distribution modeling
  #          - but also a large effect on Porites, since there is now no distinction between branching and
  #             P. astreoides
  #          - effect on Madracis may be important because of certain species being prevalent deep
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
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
  levels(combined_demo_data_trimmed$spp)
  
  #set transect as a factor
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
    mutate(transect = as.factor(transect))
  levels(combined_demo_data_trimmed$transect)
  
  #drop the unaffected species
  #   - Acropora: like 35 occurrences total
  #   - Oculina & Tubastraea: 1 occurrence each
  combined_demo_data_trimmed = combined_demo_data_trimmed %>%
    filter(!susc %in% c('Unaffected')) %>%
    mutate(susc = droplevels(susc)) %>% # drop unused factor levels
    mutate(spp = droplevels(spp))
  
  #quick read-out of species by susceptibility group
  combined_demo_data_trimmed %>%
    distinct(spp, susc) %>%
    arrange(susc, spp) %>%
    group_by(susc) %>%
    summarise(species = paste(spp, collapse = "\n  - ")) %>%
    mutate(output = paste0(susc, ":\n  - ", species)) %>%
    pull(output) %>%
    cat(sep = "\n\n")
  
  
  ################################## absences: benthic ##################################
  
  # Get all unique species levels (should be 24)
  all_species_levels <- levels(combined_benthic_data_trimmed$spp)
  
  # Create master susceptibility lookup BEFORE splitting datasets
  master_susc_lookup <- combined_benthic_data_trimmed %>% 
    select(spp, susc) %>% 
    distinct()
  
  # Function to complete species data for each dataset
  complete_species_data <- function(data) {
    # Check if transect column exists and has non-NA values
    has_transect <- "transect" %in% names(data) && any(!is.na(data$transect))
    
    if (has_transect) {
      # Has transect level (like TCRMP)
      grouping_vars <- c("dataset", "date", "PSU", "transect")
      
      # Get sampling units with metadata
      sampling_units <- data %>%
        group_by(dataset, date, PSU, transect) %>%
        summarise(
          lat = first(lat),
          lon = first(lon),
          depth = first(depth),
          .groups = "drop"
        )
      
      # Get actual cover data
      cover_data <- data %>%
        group_by(dataset, date, PSU, transect, spp) %>%
        summarise(actual_cover = sum(cover, na.rm = TRUE), .groups = "drop")
      
      join_vars <- c("dataset", "date", "PSU", "transect", "spp")
      
    } else {
      # No transect level (other datasets)
      grouping_vars <- c("dataset", "date", "PSU")
      
      # Get sampling units with metadata
      sampling_units <- data %>%
        group_by(dataset, date, PSU) %>%
        summarise(
          lat = first(lat),
          lon = first(lon),
          depth = first(depth),
          .groups = "drop"
        ) %>%
        mutate(transect = NA)  # Add transect column as NA
      
      # Get actual cover data
      cover_data <- data %>%
        group_by(dataset, date, PSU, spp) %>%
        summarise(actual_cover = sum(cover, na.rm = TRUE), .groups = "drop")
      
      join_vars <- c("dataset", "date", "PSU", "spp")
    }
    
    # Create complete grid of all combinations
    complete_grid <- sampling_units %>%
      crossing(spp = factor(all_species_levels, levels = all_species_levels)) %>%
      # Use the MASTER susceptibility lookup instead of dataset-specific one
      left_join(master_susc_lookup, by = "spp") %>%
      # Left join with original data to get actual cover values
      left_join(cover_data, by = join_vars) %>%
      # Fill in missing values
      mutate(
        cover = ifelse(is.na(actual_cover), 0, actual_cover),
        inferred_absence = case_when(
          is.na(actual_cover) ~ "Y",  # Species not recorded
          actual_cover == 0 ~ "N",   # Species recorded as 0
          TRUE ~ "N"                 # Species recorded with positive cover
        )
      ) %>%
      select(-actual_cover) %>%
      # Reorder columns to match original
      select(dataset, date, PSU, transect, lat, lon, cover, spp, depth, susc, inferred_absence)
    
    return(complete_grid)
  }
  
  # Apply to each dataset separately
  datasets <- split(combined_benthic_data_trimmed, combined_benthic_data_trimmed$dataset)
  combined_benthic_data_complete <- map_dfr(datasets, complete_species_data)
  
  # Check the results
  cat("Original data dimensions:", dim(combined_benthic_data_trimmed), "\n")
  cat("Complete data dimensions:", dim(combined_benthic_data_complete), "\n")
  cat("Number of inferred absences:", sum(combined_benthic_data_complete$inferred_absence == "Y"), "\n")
  
  # Verify all species are present for each sampling unit
  sampling_unit_check <- combined_benthic_data_complete %>%
    group_by(dataset, date, PSU, transect) %>%
    summarise(
      n_species = n_distinct(spp),
      .groups = "drop"
    )
  
  cat("Species count per sampling unit (should all be 24):\n")
  print(table(sampling_unit_check$n_species))
  
  # Additional check: verify susceptibility data is complete
  susc_check <- combined_benthic_data_complete %>%
    group_by(spp) %>%
    summarise(
      n_na_susc = sum(is.na(susc)),
      unique_susc = n_distinct(susc, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\nSusceptibility data check:\n")
  print(susc_check)
  if(any(susc_check$n_na_susc > 0)) {
    cat("WARNING: Some species have missing susceptibility data!\n")
  } else {
    cat("All species have susceptibility data - SUCCESS!\n")
  }  
  
  combined_benthic_data_summed = combined_benthic_data_complete
  
  ################################## spp-grouping ##################################
  
  # Average coral cover across transects and dates for each PSU
  combined_benthic_data_averaged <- combined_benthic_data_summed %>%
    group_by(dataset, PSU, spp) %>%
    summarise(
      # Average cover across all transects/dates within each PSU
      cover = mean(cover, na.rm = TRUE),
      # Keep other variables (taking first value since they should be consistent within PSU)
      lat = first(lat),
      lon = first(lon),
      depth = first(depth),
      susc = first(susc),
      # Set date to NA if we're averaging across multiple samples, otherwise keep original date
      date = if_else(n() > 1, as.Date(NA), first(date)),
      # Set transect to NA for everything since we're now at PSU level
      transect = NA_real_,
      .groups = 'drop'
    ) %>%
    # Reorder columns to match original structure
    select(dataset, date, PSU, transect, lat, lon, cover, spp, depth, susc)  
  
  # test = combined_benthic_data_summed %>%
  #   filter(spp == 'Meandrina spp',
  #          PSU == 'Black Point') %>%
  #   summarise(avg = mean(cover, na.rm = TRUE))
  # test
  
  #summation by PSU/susc, which accounts for repeat observations but also 
  #   multiple transects in the case of TCRMP. also handles lack of 'date' information
  #   for NODICE & SESAP
  #
  # Sum cover values with dataset-specific grouping by susceptibility group
  combined_benthic_data_summed_susc = combined_benthic_data_summed %>%
    group_by(dataset, PSU, date,
             transect = if_else(dataset == "TCRMP_benthic", transect, NA_real_),
             susc  # Group by susceptibility instead of species
    ) %>%
    summarise(
      # Keep other variables (taking first value since they should be consistent within groups)
      lat = first(lat),
      lon = first(lon),
      cover = sum(cover, na.rm = TRUE),  # Sum all species cover within susceptibility group
      depth = first(depth),
      # Count number of species contributing to this susceptibility group
      species_count = n_distinct(spp),
      .groups = 'drop'
    ) %>%
    # Reorder columns to match structure
    select(dataset, date, PSU, transect, lat, lon, cover, susc, depth, species_count)  
  
  # Average coral cover across transects and dates for each PSU by susceptibility group
  combined_benthic_data_averaged_susc <- combined_benthic_data_summed_susc %>%
    group_by(dataset, PSU, susc) %>%
    summarise(
      # Average cover across all transects/dates within each PSU
      cover = mean(cover, na.rm = TRUE),
      # Keep other variables (taking first value since they should be consistent within PSU)
      lat = first(lat),
      lon = first(lon),
      depth = first(depth),
      # Average species count across transects
      species_count = mean(species_count, na.rm = TRUE),
      # Set date to NA if we're averaging across multiple samples, otherwise keep original date
      date = if_else(n() > 1, as.Date(NA), first(date)),
      # Set transect to NA for everything since we're now at PSU level
      transect = NA_real_,
      .groups = 'drop'
    ) %>%
    # Reorder columns to match structure
    select(dataset, date, PSU, transect, lat, lon, cover, susc, depth, species_count)
  
  ## DEMO
  #
  #summation by PSU/spp, which accounts for repeat observations but also 
  #   multiple transects in the case of TCRMP
  #
  # NOTE - there are some outliers here, mainly on the 10-cm belt transects, because a few very large colonies
  #         were considered part of the transect even though it is very narrow. will likely need to cut
  #         those dates out, or simply all 10-cm belt transects for consistency
  #
  # Sum SA values with dataset-specific grouping
  combined_demo_data_summed <- combined_demo_data_trimmed %>%
    group_by(
      dataset, PSU, date,
      # # Keep original date column for grouping, but handle NAs properly
      # date = if_else(dataset %in% c("SESAP", "NODICE"), as.Date(NA), date),
      # Keep original transect, but set to NA for non-TCRMP datasets
      # transect = if_else(dataset == "TCRMP_health", transect, NA_real_),
      transect = if_else(dataset == "TCRMP_health", transect, NA_character_), # Changed to NA_character_
      spp
    ) %>%
    summarise(
      # Keep other variables (taking first value since they should be consistent within groups)
      lat = first(lat),
      lon = first(lon),
      SA = sum(SA, na.rm = TRUE),
      depth = first(depth),
      susc = first(susc),
      method = first(method),
      meterscompleted = first(meterscompleted),
      surveyarea = first(surveyarea),
      .groups = 'drop'
    ) %>%
    # Reorder columns to match original structure
    select(dataset, date, PSU, transect, lat, lon, depth, spp, method, meterscompleted, surveyarea,
           SA, susc)  
  
  #calculate transect- and species-level surface area which is scaled to 2D transect area
  combined_demo_data_summed = combined_demo_data_summed %>%
    mutate(SA_density = SA / surveyarea)
  
  # Average coral SA across transects and dates for each PSU
  combined_demo_data_averaged <- combined_demo_data_summed %>%
    group_by(dataset, PSU, spp) %>%
    summarise(
      # Average SA across all transects/dates within each PSU
      SA = mean(SA, na.rm = TRUE),
      SA_density = mean(SA_density, na.rm = TRUE),
      # Keep other variables (taking first value since they should be consistent within PSU)
      lat = first(lat),
      lon = first(lon),
      depth = first(depth),
      susc = first(susc),
      # Set date to NA if we're averaging across multiple samples, otherwise keep original date
      date = if_else(n() > 1, as.Date(NA), first(date)),
      # Set transect to NA for everything since we're now at PSU level
      transect = NA_real_,
      .groups = 'drop'
    ) %>%
    # Reorder columns to match original structure
    select(dataset, date, PSU, transect, lat, lon, SA, spp, depth, susc, SA_density)  
  
  #fill in absences
  #
  # Get all unique species across all datasets
  all_species <- combined_demo_data_averaged %>%
    distinct(spp) %>%
    pull(spp)
  
  # Get all unique combinations of dataset/PSU (with their associated metadata)
  all_locations <- combined_demo_data_averaged %>%
    select(dataset, PSU, lat, lon, depth) %>%
    distinct()
  
  # Create complete grid of all location/species combinations
  complete_grid <- expand_grid(
    location_data = all_locations,
    spp = all_species
  ) %>%
    # Unnest the location data
    unnest(location_data)
  
  # Join with existing data to fill in absences
  combined_demo_data_averaged <- complete_grid %>%
    left_join(combined_demo_data_averaged, 
              by = c("dataset", "PSU", "spp", "lat", "lon", "depth")) %>%
    # Fill in missing values
    mutate(
      # Mark inferred absences
      inferred_absence = if_else(is.na(SA), "Y", "N"),
      # Fill in zeros for SA and SA_density where missing
      SA = if_else(is.na(SA), 0, SA),
      SA_density = if_else(is.na(SA_density), 0, SA_density),
      # For inferred absences, we need to assign susc values
      # You may want to modify this logic based on your species data
      susc = if_else(is.na(susc), "unknown", susc),
      # Keep date and transect as NA for inferred absences
      date = if_else(inferred_absence == "Y", as.Date(NA), date),
      transect = if_else(inferred_absence == "Y", NA_real_, transect)
    ) %>%
    # Reorder columns to match your structure plus new column
    select(dataset, date, PSU, transect, lat, lon, SA, spp, depth, susc, SA_density, inferred_absence)

  #summation by PSU/susc, which accounts for repeat observations but also 
  #   multiple transects in the case of TCRMP
  #
  # Sum SA values with dataset-specific grouping by susceptibility group
  combined_demo_data_summed_susc <- combined_demo_data_trimmed %>%
    group_by(
      dataset, PSU, date,
      # # Keep original date column for grouping, but handle NAs properly
      # date = if_else(dataset %in% c("SESAP", "NODICE"), as.Date(NA), date),
      # Keep original transect, but set to NA for non-TCRMP datasets
      transect = if_else(dataset == "TCRMP_health", transect, NA_character_),
      susc  # Group by susceptibility instead of species
    ) %>%
    summarise(
      # Keep other variables (taking first value since they should be consistent within groups)
      lat = first(lat),
      lon = first(lon),
      SA = sum(SA, na.rm = TRUE),  # Sum all species within susceptibility group
      depth = first(depth),
      method = first(method),
      meterscompleted = first(meterscompleted),
      surveyarea = first(surveyarea),
      # Count number of species contributing to this susceptibility group
      species_count = n_distinct(spp),
      .groups = 'drop'
    ) %>%
    # Reorder columns to match structure
    select(dataset, date, PSU, transect, lat, lon, depth, susc, method, meterscompleted, surveyarea,
           SA, species_count)  
  
  #calculate transect- and susceptibility-level surface area which is scaled to 2D transect area
  combined_demo_data_summed_susc = combined_demo_data_summed_susc %>%
    mutate(SA_density = SA / surveyarea)
  
  # Average coral SA across transects and dates for each PSU by susceptibility group
  combined_demo_data_averaged_susc <- combined_demo_data_summed_susc %>%
    group_by(dataset, PSU, susc) %>%
    summarise(
      # Average SA across all transects/dates within each PSU
      SA = mean(SA, na.rm = TRUE),
      SA_density = mean(SA_density, na.rm = TRUE),
      # Keep other variables (taking first value since they should be consistent within PSU)
      lat = first(lat),
      lon = first(lon),
      depth = first(depth),
      # Average species count across transects
      species_count = mean(species_count, na.rm = TRUE),
      # Set date to NA if we're averaging across multiple samples, otherwise keep original date
      date = if_else(n() > 1, as.Date(NA), first(date)),
      # Set transect to NA for everything since we're now at PSU level
      transect = NA_character_,
      .groups = 'drop'
    ) %>%
    # Reorder columns to match structure
    select(dataset, date, PSU, transect, lat, lon, SA, susc, depth, SA_density, species_count)
  
  #fill in absences
  #
  # Get all unique species across all datasets
  all_susc <- combined_demo_data_averaged_susc %>%
    distinct(susc) %>%
    pull(susc)
  
  # Get all unique combinations of dataset/PSU (with their associated metadata)
  all_locations <- combined_demo_data_averaged_susc %>%
    select(dataset, PSU, lat, lon, depth) %>%
    distinct()
  
  # Create complete grid of all location/species combinations
  complete_grid <- expand_grid(
    location_data = all_locations,
    susc = all_susc
  ) %>%
    # Unnest the location data
    unnest(location_data)
  
  # Join with existing data to fill in absences
  combined_demo_data_averaged_susc <- complete_grid %>%
    left_join(combined_demo_data_averaged_susc, 
              by = c("dataset", "PSU", "susc", "lat", "lon", "depth")) %>%
    # Fill in missing values
    mutate(
      # Mark inferred absences
      inferred_absence = if_else(is.na(SA), "Y", "N"),
      # Fill in zeros for SA and SA_density where missing
      SA = if_else(is.na(SA), 0, SA),
      SA_density = if_else(is.na(SA_density), 0, SA_density),
      # Keep date and transect as NA for inferred absences
      date = if_else(inferred_absence == "Y", as.Date(NA), date),
      transect = if_else(inferred_absence == "Y", NA_character_, transect)  # Changed to NA_character_
      # Note: removed the susc line since it's already correctly populated from the grid
    ) %>%
    # Reorder columns to match your structure plus new column
    select(dataset, date, PSU, transect, lat, lon, SA, depth, susc, SA_density, inferred_absence)  
  
  ################################## site-level grouping ##################################
  
  # NOTE - there are zero NCRMP PSU's with absolute zero coral coverage. I am wondering if some were
  #         in fact recorded but filtered before hosting on Github? if so, would be very useful to retrieve
  #         missing absence data
  
  # BENTHIC
  #
  #summation by PSU only, pooling all species together
  #   accounts for repeat observations and multiple transects in the case of TCRMP
  #
  # Sum cover values across all species with dataset-specific grouping
  combined_benthic_data_summed_psu <- combined_benthic_data_summed %>%
    group_by(
      dataset, PSU, date,
      # # Keep original date column for grouping, but handle NAs properly
      # date = if_else(dataset %in% c("SESAP", "NODICE"), as.Date(NA), date),
      # Keep original transect, but set to NA for non-TCRMP datasets
      transect = if_else(dataset == "TCRMP_benthic", transect, NA_real_)
      # No species or susceptibility grouping - pooling all species
    ) %>%
    summarise(
      # Keep other variables (taking first value since they should be consistent within groups)
      lat = first(lat),
      lon = first(lon),
      cover = sum(cover, na.rm = TRUE),  # Sum all species cover together
      depth = first(depth),
      # Count total number of species and susceptibility groups
      total_species_count = n_distinct(spp),
      high_susc_count = sum(susc == "high", na.rm = TRUE),
      moderate_susc_count = sum(susc == "moderate", na.rm = TRUE),
      low_susc_count = sum(susc == "low", na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Reorder columns
    select(dataset, date, PSU, transect, lat, lon, depth, cover, 
           total_species_count, high_susc_count, moderate_susc_count, low_susc_count)  
  
  # Average total coral cover across transects and dates for each PSU
  combined_benthic_data_averaged_psu <- combined_benthic_data_summed_psu %>%
    group_by(dataset, PSU) %>%
    summarise(
      # Average cover across all transects/dates within each PSU
      cover = mean(cover, na.rm = TRUE),
      # Keep other variables (taking first value since they should be consistent within PSU)
      lat = first(lat),
      lon = first(lon),
      depth = first(depth),
      # Average species counts across transects
      total_species_count = mean(total_species_count, na.rm = TRUE),
      high_susc_count = mean(high_susc_count, na.rm = TRUE),
      moderate_susc_count = mean(moderate_susc_count, na.rm = TRUE),
      low_susc_count = mean(low_susc_count, na.rm = TRUE),
      # Set date to NA if we're averaging across multiple samples, otherwise keep original date
      date = if_else(n() > 1, as.Date(NA), first(date)),
      # Set transect to NA for everything since we're now at PSU level
      transect = NA_real_,
      .groups = 'drop'
    ) %>%
    # Reorder columns
    select(dataset, date, PSU, transect, lat, lon, cover, depth,
           total_species_count, high_susc_count, moderate_susc_count, low_susc_count)
  
  # DEMO
  #
  #summation by PSU only, pooling all species together
  #   accounts for repeat observations and multiple transects in the case of TCRMP
  #
  # Sum SA values across all species with dataset-specific grouping
  combined_demo_data_summed_psu <- combined_demo_data_trimmed %>%
    group_by(
      dataset, PSU, date,
      # # Keep original date column for grouping, but handle NAs properly
      # date = if_else(dataset %in% c("SESAP", "NODICE"), as.Date(NA), date),
      # Keep original transect, but set to NA for non-TCRMP datasets
      transect = if_else(dataset == "TCRMP_health", transect, NA_character_)
      # No species or susceptibility grouping - pooling all species
    ) %>%
    summarise(
      # Keep other variables (taking first value since they should be consistent within groups)
      lat = first(lat),
      lon = first(lon),
      SA = sum(SA, na.rm = TRUE),  # Sum all species together
      depth = first(depth),
      method = first(method),
      meterscompleted = first(meterscompleted),
      surveyarea = first(surveyarea),
      # Count total number of species and susceptibility groups
      total_species_count = n_distinct(spp),
      high_susc_count = sum(susc == "high", na.rm = TRUE),
      moderate_susc_count = sum(susc == "moderate", na.rm = TRUE),
      low_susc_count = sum(susc == "low", na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Reorder columns
    select(dataset, date, PSU, transect, lat, lon, depth, method, meterscompleted, surveyarea,
           SA, total_species_count, high_susc_count, moderate_susc_count, low_susc_count)  
  
  #calculate transect-level total surface area scaled to 2D transect area
  combined_demo_data_summed_psu = combined_demo_data_summed_psu %>%
    mutate(SA_density = SA / surveyarea)
  
  # Average total coral SA across transects and dates for each PSU
  combined_demo_data_averaged_psu <- combined_demo_data_summed_psu %>%
    group_by(dataset, PSU) %>%
    summarise(
      # Average SA across all transects/dates within each PSU
      SA = mean(SA, na.rm = TRUE),
      SA_density = mean(SA_density, na.rm = TRUE),
      # Keep other variables (taking first value since they should be consistent within PSU)
      lat = first(lat),
      lon = first(lon),
      depth = first(depth),
      # Average species counts across transects
      total_species_count = mean(total_species_count, na.rm = TRUE),
      high_susc_count = mean(high_susc_count, na.rm = TRUE),
      moderate_susc_count = mean(moderate_susc_count, na.rm = TRUE),
      low_susc_count = mean(low_susc_count, na.rm = TRUE),
      # Set date to NA if we're averaging across multiple samples, otherwise keep original date
      date = if_else(n() > 1, as.Date(NA), first(date)),
      # Set transect to NA for everything since we're now at PSU level
      transect = NA_character_,
      .groups = 'drop'
    ) %>%
    # Reorder columns
    select(dataset, date, PSU, transect, lat, lon, SA, depth, SA_density, 
           total_species_count, high_susc_count, moderate_susc_count, low_susc_count)
  
  
  # NOTE - from what I can tell, absences are automatically filtered down to the site-level for both
  #         benthic and demo data
  
  
  ################################## plot total cover ##################################
  
  
  # STOPPING POINT - 2 June 2025. am about to produce plots, then append bathymetry & derived
  #   values to dataframes, then produce simple test GAMs. finally will need to figure out prediction
  #   onto a uniform grid which will "talk" to the habitat grid plugged into the CMS
  
  
  hist(combined_benthic_data_averaged_psu$cover)
  hist(combined_demo_data_averaged_psu$SA_density)
  
  # Bar plot version for ultimate clarity at low values
  low_cover_data <- combined_benthic_data_averaged_psu %>%
    filter(cover <= 5) %>%
    mutate(cover_rounded = round(cover * 4) / 4)  # Round to nearest 0.25%

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
  
  # Basic detailed histogram
  ggplot(combined_benthic_data_averaged_psu, aes(x = cover)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
    scale_x_continuous(
      breaks = seq(0, max(combined_benthic_data_averaged_psu$cover, na.rm = TRUE), by = 5),
      minor_breaks = seq(0, max(combined_benthic_data_averaged_psu$cover, na.rm = TRUE), by = 1)
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
  ggplot(combined_benthic_data_averaged_psu, aes(x = cover)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_density(color = "red", size = 1) +
    scale_x_continuous(
      breaks = seq(0, max(combined_benthic_data_averaged_psu$cover, na.rm = TRUE), by = 5),
      minor_breaks = seq(0, max(combined_benthic_data_averaged_psu$cover, na.rm = TRUE), by = 1)
    ) +
    labs(
      title = "Distribution of Total Coral Cover by PSU",
      subtitle = "With density curve overlay",
      x = "Total Coral Cover (%)",
      y = "Density"
    ) +
    theme_minimal()
  
  # Version with fine bin width to distinguish low values clearly
  ggplot(combined_benthic_data_averaged_psu, aes(x = cover)) +
    geom_histogram(binwidth = 0.5, fill = "coral", alpha = 0.7, color = "white", size = 0.3) +
    scale_x_continuous(
      breaks = seq(0, max(combined_benthic_data_averaged_psu$cover, na.rm = TRUE), by = 1),
      minor_breaks = seq(0, max(combined_benthic_data_averaged_psu$cover, na.rm = TRUE), by = 0.5),
      limits = c(-0.25, max(combined_benthic_data_averaged_psu$cover, na.rm = TRUE) + 1)
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
  ggplot(combined_benthic_data_averaged_psu, aes(x = cover)) +
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
  low_cover_data <- combined_benthic_data_averaged_psu %>%
    filter(cover <= 5) %>%
    mutate(cover_rounded = round(cover * 4) / 4)  # Round to nearest 0.25%

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
  summary_stats <- combined_benthic_data_averaged_psu %>%
    summarise(
      mean_cover_all = mean(cover, na.rm = TRUE),
      mean_cover_nonzero = mean(cover[cover > 0], na.rm = TRUE),
      median_cover = median(cover, na.rm = TRUE),
      median_cover_nonzero = median(cover[cover > 0], na.rm = TRUE),
      sd_cover_nonzero = sd(cover[cover > 0], na.rm = TRUE),
      min_cover = min(cover, na.rm = TRUE),
      max_cover = max(cover, na.rm = TRUE),
      zero_cover_sites = sum(cover == 0, na.rm = TRUE),
      nonzero_cover_sites = sum(cover > 0, na.rm = TRUE),
      total_sites = n(),
      percent_zero_sites = round(100 * zero_cover_sites / total_sites, 1)
    )

  print(summary_stats)

  # Additional breakdown by cover categories
  cover_categories <- combined_benthic_data_averaged_psu %>%
    mutate(
      cover_category = case_when(
        cover == 0 ~ "Zero cover (0%)",
        cover > 0 & cover <= 1 ~ "Very low (0-1%)",
        cover > 1 & cover <= 5 ~ "Low (1-5%)",
        cover > 5 & cover <= 10 ~ "Moderate (5-10%)",
        cover > 10 & cover <= 25 ~ "High (10-25%)",
        cover > 25 ~ "Very high (>25%)"
      )
    ) %>%
    count(cover_category) %>%
    mutate(percentage = round(100 * n / sum(n), 1))

  print("Cover category breakdown:")
  print(cover_categories)

  
  
  
  
  
  # Stacked bar chart showing proportion of zero vs non-zero sites by susceptibility
  susc_zero_nonzero_psu <- combined_benthic_data_averaged_psu %>%
    mutate(
      # Create high/moderate/low categories based on the count columns
      # You may need to adjust this logic based on how you want to categorize PSUs
      dominant_susc = case_when(
        high_susc_count >= moderate_susc_count & high_susc_count >= low_susc_count ~ "high",
        moderate_susc_count >= low_susc_count ~ "moderate",
        TRUE ~ "low"
      ),
      dominant_susc = factor(dominant_susc, levels = c("low", "moderate", "high")),
      has_cover = ifelse(cover > 0, "Has coral cover", "Zero coral cover")
    ) %>%
    count(dominant_susc, has_cover) %>%
    group_by(dominant_susc) %>%
    mutate(
      total_sites = sum(n),
      percentage = round(100 * n / total_sites, 1)
    )
  
  ggplot(susc_zero_nonzero_psu, aes(x = dominant_susc, y = n, fill = has_cover)) +
    geom_col(position = "stack", alpha = 0.8) +
    geom_text(aes(label = paste0(n, "\n(", percentage, "%)")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("Zero coral cover" = "grey60", "Has coral cover" = "coral"),
      name = "Cover Status"
    ) +
    labs(
      title = "PSUs with vs without Coral Cover by Dominant Susceptibility Group",
      subtitle = "Numbers and percentages shown for each category",
      x = "Dominant Susceptibility Group",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  
  # Alternative: Group by dataset instead of susceptibility
  dataset_zero_nonzero_psu <- combined_benthic_data_averaged_psu %>%
    mutate(
      has_cover = ifelse(cover > 0, "Has coral cover", "Zero coral cover")
    ) %>%
    count(dataset, has_cover) %>%
    group_by(dataset) %>%
    mutate(
      total_sites = sum(n),
      percentage = round(100 * n / total_sites, 1)
    )
  
  ggplot(dataset_zero_nonzero_psu, aes(x = dataset, y = n, fill = has_cover)) +
    geom_col(position = "stack", alpha = 0.8) +
    geom_text(aes(label = paste0(n, "\n(", percentage, "%)")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("Zero coral cover" = "grey60", "Has coral cover" = "coral"),
      name = "Cover Status"
    ) +
    labs(
      title = "PSUs with vs without Coral Cover by Dataset",
      subtitle = "Numbers and percentages shown for each category",
      x = "Dataset",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  
  
  # Alternative: Group by dataset instead of susceptibility
  dataset_zero_nonzero_psu_demo <- combined_demo_data_averaged_psu %>%
    mutate(
      has_corals = ifelse(SA_density > 0, "Has corals", "No corals")
    ) %>%
    count(dataset, has_corals) %>%
    group_by(dataset) %>%
    mutate(
      total_sites = sum(n),
      percentage = round(100 * n / total_sites, 1)
    )
  
  ggplot(dataset_zero_nonzero_psu_demo, aes(x = dataset, y = n, fill = has_corals)) +
    geom_col(position = "stack", alpha = 0.8) +
    geom_text(aes(label = paste0(n, "\n(", percentage, "%)")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("Zero corals" = "grey60", "Has corals" = "coral"),
      name = "Demo/Health Status"
    ) +
    labs(
      title = "Demo/Health survey PSUs with vs without corals by Dataset",
      subtitle = "Numbers and percentages shown for each category",
      x = "Dataset",
      y = "Number of PSUs"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ################################## plot by susceptibility ##################################

  # Summary statistics by susceptibility group
  susc_summary <- combined_benthic_data_averaged_susc %>%
    mutate(susc = factor(susc, levels = c("low", "moderate", "high"))) %>%
    group_by(susc) %>%
    summarise(
      mean_cover_all = mean(cover, na.rm = TRUE),
      mean_cover_nonzero = mean(cover[cover > 0], na.rm = TRUE),
      median_cover_nonzero = median(cover[cover > 0], na.rm = TRUE),
      sd_cover_nonzero = sd(cover[cover > 0], na.rm = TRUE),
      min_cover = min(cover, na.rm = TRUE),
      max_cover = max(cover, na.rm = TRUE),
      zero_cover_sites = sum(cover == 0, na.rm = TRUE),
      nonzero_cover_sites = sum(cover > 0, na.rm = TRUE),
      total_sites = n(),
      percent_zero_sites = round(100 * zero_cover_sites / total_sites, 1),
      .groups = 'drop'
    )

  print("Summary statistics by susceptibility group:")
  print(susc_summary)

  # Histogram by susceptibility group
  ggplot(combined_benthic_data_averaged_susc %>%
           mutate(susc = factor(susc, levels = c("low", "moderate", "high", "Unaffected", "Unknown"))),
         aes(x = cover, fill = susc)) +
    geom_histogram(binwidth = 0.5, alpha = 0.7, color = "white", size = 0.2) +
    facet_wrap(~susc, scales = "free_y", ncol = 2) +
    scale_x_continuous(
      breaks = seq(0, max(combined_benthic_data_averaged_susc$cover, na.rm = TRUE), by = 2),
      minor_breaks = seq(0, max(combined_benthic_data_averaged_susc$cover, na.rm = TRUE), by = 1)
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
  ggplot(combined_benthic_data_averaged_susc %>%
           filter(cover > 0) %>%
           mutate(susc = factor(susc, levels = c("low", "moderate", "high", "Unaffected", "Unknown"))),
         aes(x = susc, y = cover, fill = susc)) +
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
  susc_zero_nonzero <- combined_benthic_data_averaged_susc %>%
    mutate(
      susc = factor(susc, levels = c("low", "moderate", "high", "Unaffected", "Unknown")),
      has_cover = ifelse(cover > 0, "Has coral cover", "Zero coral cover")
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
  
  
  ################################## Save output ##################################
  
  # # # Save specific combined datasets in an .rda file for stats downstream
  # save(combined_benthic_data_averaged, 
  #      combined_benthic_data_averaged_susc, 
  #      combined_benthic_data_averaged_psu,
  #      combined_demo_data_averaged, 
  #      combined_demo_data_averaged_susc, 
  #      combined_demo_data_averaged_psu,
  #      file = here("output", "all_combined_data.rda"))
  
  # #pass workspace to downstream script
  # save.image(file = here("output", "process_species_data.RData"))