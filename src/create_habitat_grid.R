  
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
  
  # # Print GDAL, PROJ, and GEOS versions linked with sf
  # sf::sf_extSoftVersion()
  # 
  # # Print GDAL version linked with terra
  # terra::gdal()
  
  ################################## setup ##################################

  load_spat_objects(directory = 'output/output_import_merge_rasters/') #call function
  load(here('output', 'output_import_merge_rasters/import_merge_rasters_workspace.RData')) #load workspace from upstream script

  # # all various things I tried to get things working. not currently required on my M1 Macbook, newly updated
  # install.packages("sf", dependencies = TRUE) #this may be required after updating to Sequoia (I reinstalled R via homebrew due to issues with Sequoia and it broke some more things)
  # tools::package_dependencies(c("sf", "here", "terra", "tidyterra", "ggplot2", "tmap", "rayshader", "scico", "RColorBrewer"), recursive = TRUE) #check all dependencies for below packages
  # install.packages(c("sf", "here", "terra", "tidyterra", "ggplot2", "tmap", "rayshader", "scico", "RColorBrewer"), dependencies = TRUE) #more dependency installs, if required
  # Sys.getenv("LD_LIBRARY_PATH")
  # Sys.setenv(LD_LIBRARY_PATH = "/opt/homebrew/lib:/opt/homebrew/include:/opt/homebrew/share") #this may also be required, followed by possibly installing sf and terra from source (below)
  # Sys.setenv(PROJ_LIB = "/opt/homebrew/share/proj") #this was temporarily needed after updating M1 Macbook to Sequoia OS
  # Sys.setenv(GDAL_DATA = "/opt/homebrew/share/gdal") #this was temporarily needed after updating M1 Macbook to Sequoia OS

  # # What I had to do to get R spatial stuff working after updating to Sequoia was:
  # #   1.) Tried uninstalling everything macports, which was probably dumb because now QGIS has to be entirely reinstalled (and may
  # #       may not work now anyways through macports yet because of delays in dependency updates after Sequoia - we shall see)
  # #        - Update Oct 2024: looks like QGIS is working fine after installing the launcher and also macports version of it
  # #   2.) Tried installing R through homebrew, this was a cluster and most things installed poorly or not at all
  # #   3.) Tried installing R again through regular Mac ARM, and then installing spatial packages from source as below. this seems to work
  # # Install from source, forcing R to use the latest versions
  # install.packages("sf", type = "source") #this might not work yet for Sequoia? also maybe doesn't matter if the code runs I guess. for now
  # install.packages("terra", type = "source")
  # #check current R version
  # R.version.string

  ################################## Test simple crm-usvi ##################################
  
  # STOPPING POINT - 8 April 2025
  #   - okay new plan!! this might work!!
  #       - PR comes from crm_USVI which I pulled off NOAA servers in 2024
  #       - STT/STJ/STX come from crm_USVI which Dan pulled off NOAA servers in (?? 2019?)
  #       - Chose this because these are the best combination of resolution, no major artifacts, and extent
  #       - Will just have to figure out how to properly merge them together. this will involve clipping the 2024 crm at the PR/STT edge
  #
  #   - okay, end of day update and it's actually even more complicated than I thought
  #       - this dataviewer (https://www.ncei.noaa.gov/maps/bathymetry/?layers=nos_hydro&minx=-65.5684&maxx=-64.4215&miny=17.8471&maxy=18.906)
  #           shows CUDEM vs CRM, 3, 10, 30, and 90 m resolutions.
  #       - 30 m and 90 m res CRM are very solid, and also comprehensive (other than Anegada). I would simply use them for the entire
  #           domain, except for some slight problems:
  #         - Dan's crm (where/when did it come from ???) actually has much better resolution, and less artifacts, in the MCD & north drop
  #         - BUT, the recent 30 m crm has a bit better res south of STJ - sort of? it's a toss up there
  #         - but DEFINITELY, the recent 30 m crm has far better resolution in Coral Bay & Haulover (STJ). also east end bay Tortola. but mosaicing artifacts introduced at east end
  #           - question is - is that resolution at Coral Bay / Haulover / East End Bay worth it?? would mean stitching crm-2024 for PR to MCD, with crm-2019 for MCD to Coral Bay, with crm-2024 again for Coral Bay to Anegada
  #           - is that 3-part stitch even possible? maybe after merging and comprehensive resampling?
  #           - try it! if that seems impossible though, simply stitch crm-2024 for PR with crm-2019 for USVI
  #           - and worst case, even that doesn't work well, so I just straight up use crm-2024 for everything (artifacts and all). maybe just with Anegada tacked on with crm-2019
  #           - remember that, absolute worst case, I can always throw out "training" data at NCRMP / DCRMP points that exist over problematic sections of bathymetry
  
  ################################## Create mask <50 m ##################################
  
  seamask <- app(bathy_merged3_crm_reefdepth, fun = function(x) {
    ifelse(x < 0 & x > -50, 0, 1)
  })
  
  unique(values(seamask)) #0 is reefy depths, 1 is deep ocean, white is land

  plot(seamask, 
       main="Merge #2",
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  # ### TESTING ###
  # # # Plot the raster first
  # # plot(seamask, main="Sea Mask with Hydrological Extent Overlay",
  # #      col=c("white", "lightblue"), legend=FALSE)
  # # 
  # # # Add the polygon on top with borders only, no fill
  # # plot(hydro_extent_proj, add=TRUE, border="darkblue", lwd=2)
  # 
  # # reverse
  # plot(hydro_extent_proj, main = "Sea Mask with Hydrological Extent Overlay",
  #      col=c("white", "lightblue"), legend = FALSE)
  # # plot(seamask, add = TRUE, border = 'darkblue', lwd = 2)
  # # plot(bathy_PR_East_cropped, add = TRUE, border = 'darkblue', lwd = 2)
  # plot(bathy_PR_East_clipped, add = TRUE, border = 'darkblue', lwd = 2)
  # 
  # # Add a legend
  # legend("topright", 
  #        legend=c("Sea (value=1)", "Land (value=0)", "Hydro Extent"),
  #        fill=c("lightblue", "white", NA),
  #        border=c(NA, NA, "darkblue"),
  #        lwd=c(NA, NA, 2))
  
  # Create a color palette for bathymetry - blues from light to dark
  bathy_colors <- colorRampPalette(c("lightcyan", "cyan", "deepskyblue", "royalblue", "navy"))(100)
  
  # Plot the bathymetry first
  plot(bathy_merged3_crm_reefdepth, 
       main="Bathymetry with Sea Mask and Hydrological Extent", 
       col=bathy_colors,
       legend=TRUE)
  
  # # Add the seamask with transparency
  # plot(seamask, 
  #      col=c(NA, rgb(0, 0, 1, 0.3)), # NA for 0 values, transparent blue for 1 values
  #      legend=FALSE, 
  #      add=TRUE)
  
  # Add the polygon outline on top
  plot(hydro_extent_proj, 
       add=TRUE, 
       border="darkred", 
       lwd=2)
  
  # # Add the landmask
  # plot(landmask, 
  #      add=TRUE, 
  #      border="darkgreen", 
  #      lwd=2)
  
  # Add a legend
  legend("topright", 
         legend=c("Bathymetry", "Sea Mask", "Hydro Extent"),
         fill=c("deepskyblue", rgb(0, 0, 1, 0.3), NA),
         border=c(NA, NA, "darkred"),
         lwd=c(NA, NA, 2))
  
  ################################## Create habitat grid ##################################
  
  # NOTE - somewhere here near grid creation, could choose to also "snap" to existing operational 650-m grid 
  #         (and/or the actual NCRMP grid). Could even get creative and produce a second snap to the PR NCRMP
  #         grid and then make an intersection area (with polygon slivers)at the MCD where PR & USVI NCRMP grids
  #         meet. but that might just not be super necessary
  
  # Create a 650 x 650 m grid that spans the extent of seamask
  # First, get the extent of the seamask
  # ext_habitat <- ext(seamask)
  ext_habitat <- ext(hydro_extent_proj)
  
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
  # For example, if reefy depths are between -50 and 0 meters:
  reefy_mask <- seamask
  reefy_mask[!(seamask >= -50 & seamask <= 0)] <- NA  # Adjust these values as needed
  
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
    # geom_sf(data = polys_apr2025_operational_sf, alpha = 0.2, fill = rgb(1, 0.5, 0, 0.5), color = "darkred", size = 0.4) +
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
  
  # Ensure CRS matches
  #  NOTE - come back to this if projection issues are suspected
  if (crs(landmask) != crs(projected_crs)) {
    landmask <- project(landmask, crs(projected_crs))
    cat("Reprojected landmask to match CRS of grid_reefy\n")
  }
  
  # Convert terra SpatVector objects to sf for plotting
  landmask_sf <- st_as_sf(landmask)
  hydro_extent_proj_sf = st_as_sf(hydro_extent_proj)
  
  ggplot() +
    # geom_spatraster(data = seamask_binary) +
    geom_sf(data = reefy_poly_sf, fill = "purple", color = NA, alpha = 0.8) +
    geom_sf(data = landmask_sf, fill = "lightblue", color = NA, alpha = 0.8) +
    geom_sf(data = grid_reefy_sf, fill = "lightgreen", color = "red", size = 0.1, alpha = 0.1) +
    geom_sf(data = hydro_extent_proj_sf, fill = "lightpink", color = "red", size = 0.1, alpha = 0.1) +
    # geom_sf(data = polys_apr2025_operational_sf, alpha = 0.2, fill = rgb(1, 0.5, 0, 0.5), color = "darkred", size = 0.4) +
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
  
  # Keep only grid cells NOT completely within the landmask
  inside_land <- relate(grid_reefy, landmask, "within")  # returns logical vector
  grid_reefy_no_land <- grid_reefy[!inside_land]
  plot(grid_reefy)
  plot(grid_reefy_no_land)
  
  nrow(grid_reefy_no_land)
  
  # Convert terra SpatVector objects to sf for plotting
  grid_reefy_no_land_sf <- st_as_sf(grid_reefy_no_land)
  
  ggplot() +
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
  
  #plot centroids
  ggplot() +
    geom_sf(data = landmask_sf, fill = "lightblue", color = NA, alpha = 0.8) +
    geom_sf(data = grid_reefy_no_land_sf, fill = NA, color = "gray30", size = 0.3) +
    geom_sf(data = centroids_sf, color = "red", size = 0.3, alpha = 0.5) +  # small dots
    theme_minimal() +
    coord_sf(xlim = c(300000, 320000), ylim = c(1960000, 1970000)) +
    labs(title = "Habitat Grid with Centroids")  
  
  ################################## Push centroids away from land ##################################
  
  # NOTE - need to update "push" function so that points that started on land don't sometimes get
  #         pushed INTO land. they need to get pushed out to sea
  
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
  
  # # Add IDs back
  # centroids_touching$grid_reefy_no_land <- grid_touching$grid_reefy_no_land
  
  # Step 5: For non-touching polygons, use regular centroid
  centroids_nontouching <- centroids(grid_notouching)
  centroids_nontouching$grid_reefy_no_land <- grid_notouching$grid_reefy_no_land
  
  # Step 6: Combine both
  centroids_all <- rbind(centroids_touching, centroids_nontouching)
  
  # convert to sf and plot
  centroids_all_sf <- st_as_sf(centroids_all)
  
  ggplot() +
    geom_sf(data = landmask_sf, fill = "tan", color = NA, alpha = 0.5) +
    geom_sf(data = grid_reefy_no_land_sf, fill = NA, color = "gray70", alpha = 0.5) +
    # geom_sf(data = centroids_all_sf, color = "red", size = 1.5, alpha = 0.5) +
    geom_sf(data = centroids_all_sf, color = "red", size = 0.3, alpha = 0.5) +
    theme_minimal() +
    # coord_sf(xlim = c(300000, 320000), ylim = c(1960000, 1970000)) +
    labs(title = "Repositioned Centroids Away from Land")
  
  ggplot() +
    geom_sf(data = landmask_sf, fill = "tan", color = NA, alpha = 0.5) +
    geom_sf(data = grid_reefy_no_land_sf, fill = NA, color = "gray70", alpha = 0.5) +
    # geom_sf(data = centroids_sf, color = "red", size = 1.5, alpha = 0.5) +
    geom_sf(data = centroids_sf, color = "red", size = 0.3, alpha = 0.5) +
    theme_minimal() +
    # coord_sf(xlim = c(300000, 320000), ylim = c(1960000, 1970000)) +
    labs(title = "Repositioned Centroids Away from Land")
  
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
    theme_minimal() +
    labs(title = "Centroid Movement Away from Land - Vieques Focus") +
    # Set limits to focus on Vieques
    coord_sf(xlim = c(300000, 320000), ylim = c(1960000, 1970000))
  
  ################################## Save objects/workspace ##################################
  
  # #save terra objects and then workspace for use in downstream scripts
  # save_spat_objects() #call from functions.R
  # save.image(file = here("output", "create_habitat_grid_workspace.RData"))
  