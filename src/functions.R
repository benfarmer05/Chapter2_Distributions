  
  
  # Function to save SpatRaster and SpatVector objects
  save_spat_objects <- function(output_dir = "output") {
    all_objects <- ls(envir = .GlobalEnv)
    dir.create(here(output_dir), showWarnings = FALSE)
    
    for (obj_name in all_objects) {
      obj <- get(obj_name)
      
      if (inherits(obj, c("SpatRaster", "SpatVector", "SpatExtent"))) {
        # Check if SpatRaster has values
        if (inherits(obj, "SpatRaster")) {
          if (!terra::hasValues(obj)) {
            warning("Object ", obj_name, " is a SpatRaster with no values")
            next  # Skip saving empty rasters
          }
        }
        
        file_name <- paste0(obj_name, ".rds")
        file_path <- here(output_dir, file_name)
        
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
  
  
  # Save only *new* objects
  save_new_objects <- function(output_dir = "output", existing_objects) {
    current_objects <- ls(envir = .GlobalEnv)
    new_objects <- setdiff(current_objects, existing_objects)
    
    if (length(new_objects) == 0) {
      message("No new objects to save")
      return()
    }
    
    dir.create(here(output_dir), showWarnings = FALSE)
    
    for (obj_name in new_objects) {
      obj <- get(obj_name)
      
      if (inherits(obj, c("SpatRaster", "SpatVector", "SpatExtent"))) {
        # Check if SpatRaster has values
        if (inherits(obj, "SpatRaster") && !terra::hasValues(obj)) {
          warning("Object ", obj_name, " is a SpatRaster with no values")
          next
        }
        
        file_name <- paste0(obj_name, ".rds")
        file_path <- here(output_dir, file_name)
        saveRDS(obj, file_path)
        message("Saved NEW spatial object: ", obj_name)
      }
    }
    
    # Save non-spatial new objects separately
    new_non_spatial <- new_objects[!sapply(new_objects, function(x) inherits(get(x), c("SpatRaster", "SpatVector", "SpatExtent")))]
    
    if (length(new_non_spatial) > 0) {
      save(list = new_non_spatial, file = here(output_dir, "create_habitat_grid_workspace.RData"))
      message("Saved ", length(new_non_spatial), " new non-spatial objects to create_habitat_grid_workspace.RData")
    }
  }
  
