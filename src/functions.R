  
  
  # Function to save SpatRaster and SpatVector objects
  save_spat_objects <- function(output_dir = "output") {
    # Get all objects in the workspace
    all_objects <- ls(envir = .GlobalEnv)
    
    # Create output directory if it doesn't exist
    dir.create(here(output_dir), showWarnings = FALSE)
    
    # Loop through each object
    for (obj_name in all_objects) {
      obj <- get(obj_name)
      
      # Check if the object is a SpatRaster or SpatVector
      # if (inherits(obj, "SpatRaster") || inherits(obj, "SpatVector")) {
      if (inherits(obj, c("SpatRaster", "SpatVector", "SpatExtent"))) {
        # Create a filename based on the object name
        file_name <- paste0(obj_name, ".rds")
        file_path <- here(output_dir, file_name)
        
        # Save the object to an RDS file
        saveRDS(obj, file_path)
        
        message("Saved ", obj_name, " to ", file_path)
      }
    }
  }
  
  # Function to load SpatRaster and SpatVector objects
  # NOTE - 21 May 2025, updated to skip files with, e.g., no matching .tif instead of killing the process. also includes support for
  #   SpatExtent files
  load_spat_objects <- function(directory = here("output")) {
    
    # List all .rds files in the directory
    rds_files <- list.files(directory, pattern = "\\.rds$", full.names = TRUE)
    
    # Loop through each file and load the object
    for (file in rds_files) {
      tryCatch({
        # Load the RDS file
        obj <- readRDS(file)
        
        # Check if the object is a spatial type we care about
        if (inherits(obj, c("SpatRaster", "SpatVector", "SpatExtent"))) {
          # Determine the name of the object from the file name
          obj_name <- tools::file_path_sans_ext(basename(file))
          
          # Assign the loaded object to the global environment
          assign(obj_name, obj, envir = .GlobalEnv)
          
          # Print status message
          message(paste("Loaded:", obj_name, "of class", class(obj)))
        } else {
          warning(paste("File", basename(file), "does not contain a recognized terra object."))
        }
      }, error = function(e) {
        warning(paste("Skipping", basename(file), "due to error:", conditionMessage(e)))
      })
    }
  }
  