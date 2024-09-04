

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
    if (inherits(obj, "SpatRaster") || inherits(obj, "SpatVector")) {
      # Create a filename based on the object name
      file_name <- paste0(obj_name, ".rds")
      file_path <- here(output_dir, file_name)
      
      # Save the object to an RDS file
      saveRDS(obj, file_path)
      
      message("Saved ", obj_name, " to ", file_path)
    }
  }
}
