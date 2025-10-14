
  save_spat_objects <- function(output_dir = "output") {
    all_objects <- ls(envir = .GlobalEnv)
    dir.create(here(output_dir), showWarnings = FALSE)

    spatial_metadata <- list()

    for (obj_name in all_objects) {
      obj <- get(obj_name)

      if (inherits(obj, c("SpatRaster", "SpatVector", "SpatExtent"))) {
        message("Processing: ", obj_name)

        if (inherits(obj, "SpatRaster")) {
          if (!terra::hasValues(obj)) {
            warning("Object ", obj_name, " is a SpatRaster with no values")
            next
          }
        }

        # Debug CRS for SpatRaster and SpatVector
        if (inherits(obj, c("SpatRaster", "SpatVector"))) {
          current_crs <- terra::crs(obj)

          message("  CRS length: ", nchar(current_crs))
          message("  CRS starts with: ", substr(current_crs, 1, 50), "...")
          message("  is.na(current_crs): ", is.na(current_crs))
          message("  current_crs == '': ", current_crs == "")

          # Store CRS info
          spatial_metadata[[obj_name]] <- list(
            crs = current_crs,
            class = class(obj)[1]
          )
          message("  ✅ Added to metadata")
        }

        # Save the object
        file_name <- paste0(obj_name, ".rds")
        file_path <- here(output_dir, file_name)
        saveRDS(obj, file_path)
        message("  Saved to ", file_path)
      }
    }

    message("\n=== FINAL CHECK ===")
    message("spatial_metadata length: ", length(spatial_metadata))
    message("spatial_metadata names: ", paste(names(spatial_metadata), collapse = ", "))

    # Always try to save metadata
    if (length(spatial_metadata) > 0) {
      metadata_file <- here(output_dir, "spatial_metadata.rds")
      saveRDS(spatial_metadata, metadata_file)
      message("Saved metadata to: ", metadata_file)
    } else {
      message("No metadata to save!")
    }
  }

  load_spat_objects <- function(directory = here("output")) {
    rds_files <- list.files(directory, pattern = "\\.rds$", full.names = TRUE)
    
    message("=== Starting load_spat_objects ===")
    message("Directory: ", directory)
    message("Found ", length(rds_files), " .rds files")
    
    # Load spatial metadata if available
    metadata_file <- here(directory, "spatial_metadata.rds")
    spatial_metadata <- if (file.exists(metadata_file)) {
      message("Loading spatial metadata from: ", metadata_file)
      readRDS(metadata_file)
    } else {
      message("No spatial metadata file found")
      NULL
    }
    
    crs_restored_count <- 0
    unpacked_count <- 0
    
    for (file in rds_files) {
      message("\n--- Processing file: ", basename(file), " ---")
      
      if (basename(file) == "spatial_metadata.rds") {
        message("Skipping metadata file")
        next
      }
      
      tryCatch({
        message("Reading RDS file...")
        obj <- readRDS(file)
        
        message("Object class: ", class(obj)[1])
        
        # NEW: Check if object is a Packed spatial object (S4 class) and unwrap it
        if (inherits(obj, c("PackedSpatRaster", "PackedSpatVector", "PackedSpatExtent"))) {
          message("🔍 Detected ", class(obj)[1], " - attempting to unwrap...")
          obj <- terra::unwrap(obj)
          message("✅ Successfully unwrapped!")
          message("Unwrapped class: ", class(obj)[1])
          unpacked_count <- unpacked_count + 1
        }
        
        message("Checking if object inherits spatial classes...")
        message("Is SpatRaster? ", inherits(obj, "SpatRaster"))
        message("Is SpatVector? ", inherits(obj, "SpatVector"))
        message("Is SpatExtent? ", inherits(obj, "SpatExtent"))
        
        if (inherits(obj, c("SpatRaster", "SpatVector", "SpatExtent"))) {
          obj_name <- tools::file_path_sans_ext(basename(file))
          message("Object name will be: ", obj_name)
          
          # CRS restoration for spatial objects
          if (inherits(obj, c("SpatRaster", "SpatVector"))) {
            current_crs <- terra::crs(obj)
            message("Current CRS: ", substr(current_crs, 1, 50), "...")
            expected_metadata <- spatial_metadata[[obj_name]]
            
            # Restore CRS if corrupted and we have backup
            if (!is.null(expected_metadata)) {
              message("Found metadata for this object")
              if (current_crs != expected_metadata$crs) {
                message("Restoring CRS...")
                terra::crs(obj) <- expected_metadata$crs
                crs_restored_count <- crs_restored_count + 1
              } else {
                message("CRS matches metadata - no restoration needed")
              }
            } else {
              message("No metadata found for this object")
            }
          }
          
          message("Assigning to global environment...")
          assign(obj_name, obj, envir = .GlobalEnv)
          message("✅ Loaded: ", obj_name)
        } else {
          message("⚠️ Object does not inherit from spatial classes - skipping")
        }
      }, error = function(e) {
        message("❌ ERROR: ", conditionMessage(e))
        warning("Skipping ", basename(file), " due to error: ", conditionMessage(e))
      })
    }
    
    message("\n=== Summary ===")
    if (unpacked_count > 0) {
      message("📦 Unpacked ", unpacked_count, " packed spatial objects")
    }
    if (crs_restored_count > 0) {
      message("✅ Restored CRS for ", crs_restored_count, " objects")
    }
    message("=== Done ===")
  }
  
  save_new_objects <- function(output_dir = "output", existing_objects) {
    current_objects <- ls(envir = .GlobalEnv)
    new_objects <- setdiff(current_objects, existing_objects)
    
    if (length(new_objects) == 0) {
      message("No new objects to save")
      return()
    }
    
    dir.create(here(output_dir), showWarnings = FALSE)
    
    spatial_metadata <- list()
    
    for (obj_name in new_objects) {
      obj <- get(obj_name)
      
      if (inherits(obj, c("SpatRaster", "SpatVector", "SpatExtent"))) {
        message("Processing: ", obj_name)
        
        if (inherits(obj, "SpatRaster")) {
          if (!terra::hasValues(obj)) {
            warning("Object ", obj_name, " is a SpatRaster with no values")
            next
          }
        }
        
        # Debug CRS for SpatRaster and SpatVector
        if (inherits(obj, c("SpatRaster", "SpatVector"))) {
          current_crs <- terra::crs(obj)
          
          message("  CRS length: ", nchar(current_crs))
          message("  CRS starts with: ", substr(current_crs, 1, 50), "...")
          message("  is.na(current_crs): ", is.na(current_crs))
          message("  current_crs == '': ", current_crs == "")
          
          # Store CRS info
          spatial_metadata[[obj_name]] <- list(
            crs = current_crs,
            class = class(obj)[1]
          )
          message("  ✅ Added to metadata")
        }
        
        # Save the object
        file_name <- paste0(obj_name, ".rds")
        file_path <- here(output_dir, file_name)
        saveRDS(obj, file_path)
        message("  Saved NEW spatial object: ", obj_name, " to ", file_path)
      }
    }
    
    message("\n=== FINAL CHECK ===")
    message("spatial_metadata length: ", length(spatial_metadata))
    message("spatial_metadata names: ", paste(names(spatial_metadata), collapse = ", "))
    
    # Always try to save metadata
    if (length(spatial_metadata) > 0) {
      metadata_file <- here(output_dir, "spatial_metadata.rds")
      saveRDS(spatial_metadata, metadata_file)
      message("Saved metadata to: ", metadata_file)
    } else {
      message("No metadata to save!")
    }
    
    # Save non-spatial new objects separately
    new_non_spatial <- new_objects[!sapply(new_objects, function(x) inherits(get(x), c("SpatRaster", "SpatVector", "SpatExtent")))]
    
    if (length(new_non_spatial) > 0) {
      # Extract the directory name to create an appropriate workspace filename
      dir_name <- basename(output_dir)
      workspace_filename <- paste0(dir_name, "_workspace.RData")
      workspace_path <- here(output_dir, workspace_filename)
      
      save(list = new_non_spatial, file = workspace_path)
      message("Saved ", length(new_non_spatial), " new non-spatial objects to ", workspace_filename)
      message("Non-spatial objects saved: ", paste(new_non_spatial, collapse = ", "))
    } else {
      message("No new non-spatial objects to save")
    }
  }
  