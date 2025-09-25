  
  # .rs.restartR(clean = TRUE)
  
  library(sf)
  library(here)
  library(terra) 
  library(tidyterra)
  library(ggplot2)
  library(tmap)
  library(rayshader) #this requires installation of XQuartz on MacOS, and possibly OpenGL if it isn't installed
  library(scico)
  library(RColorBrewer)
  library(progress)
  
  source(here("src/functions.R"))
  
  ################################## setup ##################################
  
  load_spat_objects(directory = 'output/output_import_merge_rasters_higher-res/') #call function
  load(here('output', 'output_import_merge_rasters_higher-res/import_merge_rasters_workspace.RData')) #load workspace from upstream script
  
  existing_objects <- ls(envir = .GlobalEnv)
  
  ################################## Create mask <50 m ##################################
  
  #read the CMS-formatted hydrodynamic domain extent
  hydro_extent = vect(here("output/", "hydro_domain_extent.shp"))
  hydro_extent <- project(hydro_extent, projected_crs)
  
  # #redacted - not using that merge. may need to apply hydro extent differently now ?
  # seamask <- app(bathy_merged3_crm_reefdepth, fun = function(x) {
  #   ifelse(x < 0 & x > -50, 0, 1)
  # })
  
  # NOTE - 4 JULY 2025
  #   - return to this if we want to limit depth range for computational reasons (at least for CMS version)
  # seamask <- app(bathy_final, fun = function(x) {
  #   ifelse(x < 0, 0, 1)
  # })
  seamask <- app(bathy_final, fun = function(x) {
    ifelse(x < 0 & x > -60, 0, 1)
  })
  
  seamask <- crop(seamask, hydro_extent)
  seamask <- mask(seamask, hydro_extent)
  
  unique(values(seamask)) #0 is reefy depths, 1 is deep ocean, white is land
  
  # plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # plot_extents = ext(260000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(305000, 330000, 2020000, 2035000) #for investigating STJ
  # plot_extents = ext(317000, 350000, 2030000, 2050000) #for investigating Tortola
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # plot_extents = ext(300000, 340000, 1940000, 1980000) #for investigating St Croix
  plot_extents = ext(320000, 322000, 1962000, 1964000) #for investigating Altona Lagoon, STX
  # plot_extents = ext(279000, 310000, 2010000, 2050000) #for investigating St Thomas
  # plot_extents = ext(240000, 275000, 2000000, 2040000) #for investigating Culebra
  # plot_extents = ext(120000, 220000, 2020000, 2060000) #for investigating northern PR
  
  plot(seamask, 
       # ext = plot_extents, #e_crm,
       legend=TRUE)
  
  # Create a color palette for bathymetry - blues from light to dark
  bathy_colors <- colorRampPalette(c("lightcyan", "cyan", "deepskyblue", "royalblue", "navy"))(100)
  
  # plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # plot_extents = ext(300000, 340000, 1940000, 1980000) #for investigating St Croix
  # plot_extents = ext(290000, 330000, 2040000, 2080000) #for investigating St Thomas
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  
  # Plot the bathymetry first
  bathy_final_clamp <- clamp(bathy_final, lower = -60, upper = 0)
  plot(bathy_final_clamp,  #bathy_merged3_crm_reefdepth
       main="Bathymetry with Sea Mask and Hydrological Extent", 
       col=bathy_colors,
       # ext = plot_extents, #e_crm,
       legend=TRUE)
  
  ################################## Create habitat grid ##################################
  
  # NOTE - somewhere here near grid creation, could choose to also "snap" to existing operational 650-m grid 
  #         (and/or the actual NCRMP grid). Could even get creative and produce a second snap to the PR NCRMP
  #         grid and then make an intersection area (with polygon slivers)at the MCD where PR & USVI NCRMP grids
  #         meet. but that might just not be super necessary
  
  # Create a 650 x 650 m grid that spans the extent of seamask
  # First, get the extent of the seamask
  ext_habitat <- ext(seamask)
  # ext_habitat <- ext(hydro_extent_proj)
  
  # Define the grid resolution (650 m)
  grid_res <- 650
  
  # Calculate the number of cells needed
  ncol_grid <- ceiling((ext_habitat$xmax - ext_habitat$xmin) / grid_res)
  nrow_grid <- ceiling((ext_habitat$ymax - ext_habitat$ymin) / grid_res)
  
  # Create a new raster with the desired resolution
  grid_rast <- rast(ext_habitat, ncol=ncol_grid, nrow=nrow_grid) #, res = grid_res
  
  # Create a new raster with the desired resolution and SAME projection
  grid_rast <- rast(
    xmin = ext_habitat$xmin,  # Start exactly at western edge
    ymin = ext_habitat$ymin,
    xmax = ext_habitat$xmin + (ncol_grid * grid_res),  # Calculate max extent based on resolution
    ymax = ext_habitat$ymax,  # Align with northern edge
    resolution = grid_res,
    crs = projected_crs  # Use the same CRS as the seamask
  )
  # Assign cell numbers (or 1) to create a grid
  values(grid_rast) <- 1:ncell(grid_rast)  # or simply values(grid_rast) <- 1
  
  # Make sure the grid and seamask have the same CRS
  #   NOTE - will need to decide on a universal CRS that is compatible with global coordinates in USCROMS hydro
  crs(grid_rast) <- crs(projected_crs)
  
  # Convert grid raster to vector (polygons)
  grid_poly <- as.polygons(grid_rast)
  
  #assign unique IDs (NOTE - SUPER IMPORTANT)
  names(grid_poly)  # See existing field names
  names(grid_poly)[names(grid_poly) == "lyr.1"] <- "unique_ID" #rename
  any(duplicated(grid_poly$unique_ID))  # Should return FALSE if IDs are unique (which is good!)
  
  # Now clip the grid to only include the reefy depths present in seamask
  # First, identify reefy depths in the seamask (assuming these are specific values)
  # For example, if reefy depths are between -60 and 0 meters:
  reefy_mask <- seamask
  reefy_mask[!(seamask >= -60 & seamask <= 0)] <- NA  # Adjust these values as needed
  
  #eliminate artifacts introduced in the far north (if needed)
  y_threshold <- 2090000  # Adjust this value as needed
  reefy_mask[terra::yFromCell(reefy_mask, 1:terra::ncell(reefy_mask)) > y_threshold] <- NA
  
  # Convert the reefy areas to polygons
  reefy_poly <- as.polygons(reefy_mask, dissolve=TRUE)
  
  # # NOTE - USE THE BELOW ONLY IF YOU WANT A COOKIE-CUTTER GRID
  # # Clip the grid to include only cells that intersect with reefy areas
  # # Using intersect to keep only grid cells that overlap with reefy areas
  # grid_reefy <- intersect(grid_poly, reefy_poly)
  
  # Alternatively, you could select grid cells based on the center point
  grid_reefy <- grid_poly[reefy_poly]
  
  #check #/polygon squares
  nrow(grid_reefy)
  
  # Plot to check results
  plot(grid_poly, main="Full Grid (650m)", border = adjustcolor("lightgray", alpha.f = 0.2))
  plot(reefy_poly, col="lightblue", main="Reefy Areas", add = TRUE) # add=TRUE,
  plot(grid_reefy, col= NA, border = adjustcolor("red", alpha.f = 0.7), add=TRUE)
  
  # Create a standalone plot of just the final clipped grid
  plot(grid_reefy, col="lightgreen", border="red", main="650m Grid Clipped to Reefy Depths")
  
  
  
  
  
  # STOPPING POINT - 7 JULY 2025
  # 1.) don't think I need 650-m "SDM" version at all. just work with native 50-m grid
  # 2.) but here, should be clipping to hydro extent for "CMS" version
  # 3.) may need to consider depth range again - extending to 50 m already introduces some problems. 90 m
  #       is even worse
  
  
  
  # compare with April 2025 operational 650-m resolution grid from QGIS
  polys_apr2025_operational <- vect(here("output", "polys_apr2025_operational.shp"))
  
  #check #/polygon squares
  nrow(polys_apr2025_operational)
  
  # Make sure the CRS matches other spatial objects
  if (crs(polys_apr2025_operational) != crs(projected_crs)) {
    polys_apr2025_operational <- project(polys_apr2025_operational, crs(projected_crs))
    cat("Projections were different - reprojected reef polygons to match grid_reefy\n")
  }
  
  # Convert terra SpatVector objects to sf for plotting
  reefy_poly_sf <- st_as_sf(reefy_poly)
  grid_reefy_sf <- st_as_sf(grid_reefy)
  polys_apr2025_operational_sf <- st_as_sf(polys_apr2025_operational)
  
  ggplot() +
    geom_sf(data = reefy_poly_sf, fill = "lightblue", color = NA, alpha = 0.8) +
    geom_sf(data = grid_reefy_sf, fill = "lightgreen", color = "red", size = 0.1, alpha = 0.1) +
    geom_sf(data = polys_apr2025_operational_sf, alpha = 0.2, fill = rgb(1, 0.5, 0, 0.5), color = "darkred", size = 0.4) +
    scale_fill_identity() +
    scale_color_identity() +
    labs(
      title = "Reef and Grid Overlay",
      caption = "Source: Spatial data in terra SpatVector format"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  ################################## Shave off points completely on land ##################################
  
  landmask <- vect(here("output", "landmask_dissolved.shp"))
  landmask <- project(landmask, projected_crs)
  
  # Convert terra SpatVector objects to sf for plotting
  landmask_sf <- st_as_sf(landmask)
  # hydro_extent_proj_sf = st_as_sf(hydro_extent_proj)
  
  # plot_extents <- list(xmin = 280000, xmax = 310000, ymin = 2000000, ymax = 2040000)  # south of STT
  # plot_extents <- list(xmin = 270000, xmax = 290000, ymin = 2000000, ymax = 2040000)  # MCD
  plot_extents <- list(xmin = 300000, xmax = 340000, ymin = 2000000, ymax = 2050000)  # STJ
  # plot_extents <- list(xmin = 220000, xmax = 260000, ymin = 2000000, ymax = 2010000)  # Vieques
  # plot_extents <- list(xmin = 341000, xmax = 379000, ymin = 2057000, ymax = 2078000)  # Anegada
  # plot_extents <- list(xmin = 300000, xmax = 340000, ymin = 1940000, ymax = 1980000)  # St Croix
  # plot_extents <- list(xmin = 280000, xmax = 320000, ymin = 2000000, ymax = 2040000)  # St Thomas
  # plot_extents <- list(xmin = 240000, xmax = 280000, ymin = 2000000, ymax = 2040000)  # Mona Island
  
  ggplot() +
    # geom_spatraster(data = seamask_binary) +
    geom_sf(data = reefy_poly_sf, fill = "purple", color = NA, alpha = 0.8) +
    geom_sf(data = landmask_sf, fill = "lightblue", color = NA, alpha = 0.8) +
    geom_sf(data = grid_reefy_sf, fill = "lightgreen", color = "red", size = 0.1, alpha = 0.1) +
    # geom_sf(data = hydro_extent_proj_sf, fill = "lightpink", color = "red", size = 0.1, alpha = 0.1) +
    # geom_sf(data = polys_apr2025_operational_sf, alpha = 0.2, fill = rgb(1, 0.5, 0, 0.5), color = "darkred", size = 0.4) +
    scale_fill_identity() +
    scale_color_identity() +
    labs(
      title = "Reef and Grid Overlay",
      caption = "Source: Spatial data in terra SpatVector format"
    ) +
    theme_minimal() +
    coord_sf(xlim = c(plot_extents$xmin, plot_extents$xmax),
             ylim = c(plot_extents$ymin, plot_extents$ymax)) +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  # Keep only grid cells NOT completely within the landmask
  inside_land <- relate(grid_reefy, landmask, "within")  # returns logical vector
  grid_reefy_no_land <- grid_reefy[!inside_land]
  plot(grid_reefy)
  plot(grid_reefy_no_land)
  
  nrow(grid_reefy_no_land)
  
  # Convert terra SpatVector objects to sf for plotting
  grid_reefy_no_land_sf <- st_as_sf(grid_reefy_no_land)
  
  ggplot() +
    geom_sf(data = reefy_poly_sf, fill = "purple", color = NA, alpha = 0.8) +
    geom_sf(data = landmask_sf, fill = "lightblue", color = NA, alpha = 0.8) +
    geom_sf(data = grid_reefy_no_land_sf, fill = "lightgreen", color = "red", size = 0.1, alpha = 0.1) +
    # geom_sf(data = polys_apr2025_operational_sf, alpha = 0.2, fill = rgb(1, 0.5, 0, 0.5), color = "darkred", size = 0.4) +
    scale_fill_identity() +
    scale_color_identity() +
    labs(
      title = "Reef and Grid Overlay",
      caption = "Source: Spatial data in terra SpatVector format"
    ) +
    # coord_sf(xlim = x_limits, ylim = y_limits, expand = FALSE) +   # zoom in here!
    theme_minimal() +
    coord_sf(xlim = c(plot_extents$xmin, plot_extents$xmax),
             ylim = c(plot_extents$ymin, plot_extents$ymax)) +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  ################################## Place habitat centroids ##################################
  
  # inspect attribute tables
  names(grid_poly)  # See existing field names
  head(values(grid_poly))
  head(values(grid_reefy))
  head(values(grid_reefy_no_land))
  
  #inspect polygon counts
  nrow(grid_poly)
  nrow(grid_reefy)
  nrow(grid_reefy_no_land)
  
  #create centroids
  centroids <- centroids(grid_reefy_no_land)
  
  # Convert both to sf for ggplot
  centroids_sf <- st_as_sf(centroids)
  
  # plot_extents <- list(xmin = 280000, xmax = 310000, ymin = 2000000, ymax = 2040000)  # south of STT
  # plot_extents <- list(xmin = 270000, xmax = 290000, ymin = 2000000, ymax = 2040000)  # MCD
  # plot_extents <- list(xmin = 300000, xmax = 340000, ymin = 2000000, ymax = 2050000)  # STJ
  # plot_extents <- list(xmin = 220000, xmax = 260000, ymin = 2000000, ymax = 2010000)  # Vieques
  # plot_extents <- list(xmin = 341000, xmax = 379000, ymin = 2057000, ymax = 2078000)  # Anegada
  # plot_extents <- list(xmin = 300000, xmax = 340000, ymin = 1940000, ymax = 1980000)  # St Croix
  plot_extents <- list(xmin = 280000, xmax = 320000, ymin = 2000000, ymax = 2040000)  # St Thomas
  # plot_extents <- list(xmin = 240000, xmax = 280000, ymin = 2000000, ymax = 2040000)  # Mona Island
  
  #plot centroids
  ggplot() +
    geom_sf(data = landmask_sf, fill = "lightblue", color = NA, alpha = 0.8) +
    geom_sf(data = grid_reefy_no_land_sf, fill = NA, color = "gray30", size = 0.3) +
    geom_sf(data = centroids_sf, color = "red", size = 0.3, alpha = 0.5) +  # small dots
    theme_minimal() +
    # coord_sf(xlim = c(300000, 320000), ylim = c(1960000, 1970000)) +
    coord_sf(xlim = c(plot_extents$xmin, plot_extents$xmax), 
             ylim = c(plot_extents$ymin, plot_extents$ymax)) +
    labs(title = "Habitat Grid with Centroids")  
  
  ################################## Push centroids away from land ##################################
  
  # Step 1: Identify grid cells that touch land
  touches_land <- relate(grid_reefy_no_land, landmask, relation = "intersects") #|> lengths() > 0
  
  # Step 2: Split grid into touching and non-touching
  grid_touching <- grid_reefy_no_land[touches_land]
  grid_notouching <- grid_reefy_no_land[!touches_land]
  plot(grid_reefy)
  plot(grid_reefy_no_land)
  plot(grid_notouching)
  
  # Create a function to determine if points are on land or water
  point_is_on_land <- function(point, landmask) {
    return(relate(point, landmask, "intersects"))
  }
  
  # Modified farthest_from_land function
  farthest_from_land <- function(polygon, land, n = 500) {
    # Generate many random points within the polygon
    candidates <- spatSample(polygon, size = n, method = "random")
    
    # Determine which points are on land
    on_land <- point_is_on_land(candidates, land)
    
    if (any(on_land)) {
      # If some points are on land, find points that are on water
      water_points <- candidates[!on_land]
      
      if (length(water_points) > 0) {
        # If we have water points, use those
        # Get the farthest water point from land
        dists <- distance(water_points, land)
        farthest_point <- water_points[which.max(dists), ]
      } else {
        # All points are on land, push towards nearest coast
        # Create buffer around landmask and find nearest point on the buffer
        buffer <- buffer(land, width=-0.001)  # Negative buffer moves inward
        boundary <- erase(land, buffer)  # This gives us the coast
        
        # Find closest point to the coast
        dist_to_coast <- distance(candidates, boundary)
        closest_point <- candidates[which.min(dist_to_coast), ]
        farthest_point <- closest_point
      }
    } else {
      # All points are in water, find farthest from land
      dists <- distance(candidates, land)
      farthest_point <- candidates[which.max(dists), ]
    }
    
    return(farthest_point)
  }
  
  # Step 4: Apply function to all touching grid squares
  pb <- progress_bar$new( #initialize progress bar
    format = "Processing [:bar] :percent in :elapsed",
    total = nrow(grid_touching), clear = FALSE, width = 60
  )
  
  #apply function with progress
  centroids_touching <- lapply(1:nrow(grid_touching), function(i) {
    pb$tick()
    farthest_from_land(grid_touching[i], landmask)
  }) |> do.call(what = rbind)
  
  # # Add IDs back ?
  
  # Step 5: For non-touching polygons, use regular centroid
  centroids_nontouching <- centroids(grid_notouching)
  centroids_nontouching$grid_reefy_no_land <- grid_notouching$grid_reefy_no_land
  
  # Step 6: Combine both
  centroids_all <- rbind(centroids_touching, centroids_nontouching)
  
  # convert to sf and plot
  centroids_all_sf <- st_as_sf(centroids_all)
  
  # plot_extents <- list(xmin = 280000, xmax = 310000, ymin = 2000000, ymax = 2040000)  # south of STT
  # plot_extents <- list(xmin = 270000, xmax = 290000, ymin = 2000000, ymax = 2040000)  # MCD
  # plot_extents <- list(xmin = 300000, xmax = 340000, ymin = 2000000, ymax = 2050000)  # STJ
  # plot_extents <- list(xmin = 220000, xmax = 260000, ymin = 2000000, ymax = 2010000)  # Vieques
  # plot_extents <- list(xmin = 341000, xmax = 379000, ymin = 2057000, ymax = 2078000)  # Anegada
  plot_extents <- list(xmin = 300000, xmax = 340000, ymin = 1940000, ymax = 1980000)  # St Croix
  # plot_extents <- list(xmin = 280000, xmax = 320000, ymin = 2000000, ymax = 2040000)  # St Thomas
  # plot_extents <- list(xmin = 240000, xmax = 280000, ymin = 2000000, ymax = 2040000)  # Mona Island
  
  ggplot() +
    geom_sf(data = landmask_sf, fill = "tan", color = NA, alpha = 0.5) +
    geom_sf(data = grid_reefy_no_land_sf, fill = NA, color = "gray70", alpha = 0.5) +
    # geom_sf(data = centroids_all_sf, color = "red", size = 1.5, alpha = 0.5) +
    geom_sf(data = centroids_all_sf, color = "red", size = 0.3, alpha = 0.5) +
    theme_minimal() +
    coord_sf(xlim = c(plot_extents$xmin, plot_extents$xmax), 
             ylim = c(plot_extents$ymin, plot_extents$ymax)) +
    labs(title = "Repositioned points")
  
  ggplot() +
    geom_sf(data = landmask_sf, fill = "tan", color = NA, alpha = 0.5) +
    geom_sf(data = grid_reefy_no_land_sf, fill = NA, color = "gray70", alpha = 0.5) +
    # geom_sf(data = centroids_sf, color = "red", size = 1.5, alpha = 0.5) +
    geom_sf(data = centroids_sf, color = "red", size = 0.3, alpha = 0.5) +
    theme_minimal() +
    coord_sf(xlim = c(plot_extents$xmin, plot_extents$xmax), 
             ylim = c(plot_extents$ymin, plot_extents$ymax)) +
    labs(title = "Original point locations")
  
  #visualize repositioning with arrows if desired
  # First, create a dataframe with both original and new positions
  # Original centroids (before repositioning)
  original_centroids_touching <- centroids(grid_touching)
  original_centroids_nontouching <- centroids(grid_notouching)
  original_centroids_all <- rbind(original_centroids_touching, original_centroids_nontouching)

  # Convert both to sf for plotting
  original_centroids_sf <- st_as_sf(original_centroids_all)
  repositioned_centroids_sf <- st_as_sf(centroids_all)

  # Create a dataframe for the arrows
  arrow_data <- data.frame(
    x_start = st_coordinates(original_centroids_sf)[,1],
    y_start = st_coordinates(original_centroids_sf)[,2],
    x_end = st_coordinates(repositioned_centroids_sf)[,1],
    y_end = st_coordinates(repositioned_centroids_sf)[,2]
  )
  
  # Only show arrows where there was actual movement (distance > small threshold)
  movement_distance <- sqrt((arrow_data$x_end - arrow_data$x_start)^2 +
                              (arrow_data$y_end - arrow_data$y_start)^2)
  arrow_data <- arrow_data[movement_distance > 10, ]  # Adjust threshold as needed
  
  # Create the zoomed plot with arrows
  ggplot() +
    geom_sf(data = landmask_sf, fill = "tan", color = NA, alpha = 0.5) +
    geom_sf(data = grid_reefy_no_land_sf, fill = NA, color = "gray70", alpha = 0.5) +
    # Add arrows showing movement
    geom_segment(data = arrow_data,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.1, "cm")),
                 color = "blue",
                 alpha = 0.7,
                 size = 0.3) +
    # Add repositioned centroids
    geom_sf(data = repositioned_centroids_sf, color = "red", size = 0.3, alpha = 0.5) +
    # Add original centroids (optional - might make plot too busy)
    geom_sf(data = original_centroids_sf, color = "green", size = 0.2, alpha = 0.3) +
    coord_sf(xlim = c(plot_extents$xmin, plot_extents$xmax), 
             ylim = c(plot_extents$ymin, plot_extents$ymax)) +
    theme_minimal()
  
  
  ################################## output new 650_m csv for MATLAB / CMS ##################################
  
  # Transform centroids to geographic coordinates
  centroids_all_geographic <- project(centroids_all, geographic_crs)
  centroids_ordered <- centroids_all_geographic[match(grid_reefy_no_land$unique_ID, centroids_all_geographic$unique_ID), ]
  
  # Visual check setup
  use_random_seed <- TRUE
  if (use_random_seed) set.seed(as.numeric(Sys.time())) else set.seed(123)
  
  test_ids <- sample(grid_reefy_no_land$unique_ID, 100)
  test_grid <- project(grid_reefy_no_land[grid_reefy_no_land$unique_ID %in% test_ids, ], geographic_crs)
  test_centroids <- centroids_ordered[match(test_ids, centroids_ordered$unique_ID), ]
  
  # Plot extents toggle
  plot_extents <- list(xmin = -65.3, xmax = -64.9, ymin = 18.3, ymax = 18.4)  # south of STT
  # plot_extents <- list(xmin = -65.1, xmax = -64.8, ymin = 18.3, ymax = 18.4)  # MCD
  # plot_extents <- list(xmin = -64.8, xmax = -64.3, ymin = 18.3, ymax = 18.5)  # STJ
  # plot_extents <- list(xmin = -65.4, xmax = -65.0, ymin = 18.1, ymax = 18.2)  # Vieques
  # plot_extents <- list(xmin = -64.5, xmax = -64.2, ymin = 18.7, ymax = 18.8)  # Anegada
  # plot_extents <- list(xmin = -64.9, xmax = -64.5, ymin = 17.7, ymax = 17.9)  # St Croix
  # plot_extents <- list(xmin = -65.1, xmax = -64.8, ymin = 18.3, ymax = 18.4)  # St Thomas
  # plot_extents <- list(xmin = -65.3, xmax = -65.0, ymin = 18.3, ymax = 18.4)  # Mona Island
  # plot_extents <- NULL  # Use full extent
  
  # Create plot
  p <- ggplot() +
    geom_sf(data = st_as_sf(project(grid_reefy_no_land, geographic_crs)), fill = NA, color = "lightgray", size = 0.1, alpha = 0.8) +
    geom_sf(data = st_as_sf(centroids_ordered), color = "pink", size = 0.1, alpha = 0.8) +
    geom_sf(data = st_as_sf(test_grid), fill = NA, color = "blue", size = 1) +
    geom_sf(data = st_as_sf(test_centroids), color = "red", size = 2) +
    geom_sf_text(data = st_as_sf(test_grid), aes(label = unique_ID), size = 3, color = "blue", fontface = "bold") +
    geom_sf_text(data = st_as_sf(test_centroids), aes(label = unique_ID), size = 3, color = "red", fontface = "bold", nudge_y = 0.01) +
    theme_minimal() +
    labs(title = "Visual check: Do red points (with IDs) fall within their blue polygons (with IDs)?")
  
  # if (!is.null(plot_extents)) p <- p + coord_sf(xlim = c(plot_extents$xmin, plot_extents$xmax), ylim = c(plot_extents$ymin, plot_extents$ymax))
  print(p)
  
  # Create CSV
  coords_matrix <- geom(centroids_ordered)[, c("x", "y")]
  repositioned_data <- data.frame(
    unique_ID = grid_reefy_no_land$unique_ID,
    x = coords_matrix[, 1],
    y = coords_matrix[, 2]
  )
  
  write.csv(repositioned_data, here("output", "points_650_none-on-land.csv"), row.names = FALSE)
  cat("Created CSV with", nrow(repositioned_data), "points\n")
  head(repositioned_data)
  
  ################################## Save objects/workspace ##################################
  
  # #updated way to handle saving of new objects
  # save_new_objects("output/output_create_habitat_grid", existing_objects)
  