
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

    # Load spatial metadata if available
    metadata_file <- here(directory, "spatial_metadata.rds")
    spatial_metadata <- if (file.exists(metadata_file)) {
      readRDS(metadata_file)
    } else {
      NULL
    }

    crs_restored_count <- 0

    for (file in rds_files) {
      if (basename(file) == "spatial_metadata.rds") next

      tryCatch({
        obj <- readRDS(file)

        if (inherits(obj, c("SpatRaster", "SpatVector", "SpatExtent"))) {
          obj_name <- tools::file_path_sans_ext(basename(file))

          # CRS restoration for spatial objects
          if (inherits(obj, c("SpatRaster", "SpatVector"))) {
            current_crs <- terra::crs(obj)
            expected_metadata <- spatial_metadata[[obj_name]]

            # Restore CRS if corrupted and we have backup
            if (!is.null(expected_metadata) && current_crs != expected_metadata$crs) {
              terra::crs(obj) <- expected_metadata$crs
              crs_restored_count <- crs_restored_count + 1
            }
          }

          assign(obj_name, obj, envir = .GlobalEnv)
          message("Loaded: ", obj_name)
        }
      }, error = function(e) {
        warning("Skipping ", basename(file), " due to error: ", conditionMessage(e))
      })
    }

    if (crs_restored_count > 0) {
      message("✅ Restored CRS for ", crs_restored_count, " objects")
    }
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
      save(list = new_non_spatial, file = here(output_dir, "create_habitat_grid_workspace.RData"))
      message("Saved ", length(new_non_spatial), " new non-spatial objects to create_habitat_grid_workspace.RData")
    }
  }
