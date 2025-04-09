# Chain sequence for running scripts:
#   - import_merge_rasters.R
#   - calculate_bathy_rasters.R


#           Consider:

# Main script (e.g., main_script.R)

# Source the first script
source("script1.R")

# Source the second script
source("script2.R")

# Source the third script
source("script3.R")

# Continue with additional scripts as needed



- from calculate_bathy_rasters.R:

  # STOPPING POINT - 3 SEPT 2024
  #   - yes okay the more I think about it, why not just use the native NCRMP grid(s) directly to calculate bathymetry, slope, aspect etc.
  #       within their grid squares? this might be really computationally expensive, I don't know. but I'll have to figure out something
  #   - and anyways, next steps will require
  #       - test how bad the mosaicing / artifact patterns matter in Blondeau's NOAA product (produce slope/complexity and plot)
  #       - compare with NOAA CRM exports. bring in mine and Holstein's to R, clip to 0-50 m, merge them (may need to clip PR against USVI
  #         first). then produce slope/complexity and plot. compare plot with that of Blondeau's product
  #       - lastly, and actually I'll do this first, is porting data from Allen Coral Atlas and seeing how that bathy looks! and turbidity
  
    # STOPPING POINT - 3 SEPT 2024
  # - above, I was finishing up plotting the derived bathymetric products and seeing the effect of the weird patterns in the Blondeau data.
  #   so far, signs point toward me needing to make a better merged bathymetry because of the issues present. should be doable
  # - ... and below was some old scripting from 2022 that I am adopting for the here and now
  
#   4 SEP 2024
  #   Okay, now I'm starting to understand the computational limits I'm working against. Merging the bathymetry or sample frame rasters
  #     across regions is kind of out of the question, at least as long as I am interested in Puerto Rico being in the same analysis as the
  #     VI. So, what makes the most sense is deriving all bathymetry products at 2 m resolution, within each section of the bathymetry (if
  #     that is even possible ?). Assuming that is possible though, then we have the problem of all the issues in the Blondeau files.
  #     ideally, I would just use a bathymetry data source that doesn't have the weird patterns. maybe, no matter the source of the bathy,
  #     I'll need to break it up into chunks, derive the products necessary slope, aspect, etc.) and then aggregate their mean & SD within
  #     the relevant sample frames? pretty messy. phew
  #
  #   - big problem: what is in fact coastline is for some reason underwater in the Blondeau output. I could mask it using other bathy
  #     and set that section to NA. but I would want to be confident that there aren't major issues right next to the "new" correct coastline
  #   - another thought, just in general, is clipping all rasters to 50 m depth and seeing how much that helps with processing times. not
  #       only clipping but also using 'ifel' as above
  #   - last thought...would be great to simply extract the areas of highest resolution from all the rasters I have to work with. but is this
  #       feasible? would be easiest if I had proper access to the underlying tiles (of different resolutions) making up the products I've
  #       been given or found. yeehaw
  #
  # ### Port in Allen Coral Atlas bathymetry & turbidity data. maybe try this using CLI in the future - https://openoceans.xyz/projects/pycoral/
  # # okay upon looking at the bathy - it's very grainy. I guess that makes sense for satellite imaging, but either way cannot use.
  # # note that it is in centimeters, and >depth = positive values.
  # #   - the turbidity raster also doesn't seem like it has great resolution either. very focused on super nearshore reefs. could be useful
  # #     for PR though ? not sure
  # bathy_ACA = rast(here("data/Allen_Coral_Atlas/PR_USVI-20240830214543/Bathymetry---composite-depth/bathymetry_0.tif"))
  # turbidity_ACA = rast(here("data/Allen_Coral_Atlas/PR_USVI-20240830214543/Turbidity-2019/turbidity-annual_0.tif"))
  # bathy_ACA <- clamp(bathy_ACA, lower = 0, upper = 50, values = TRUE)
  # tm_shape(bathy_ACA) +
  #   tm_raster(style = 'cont')
  # tm_shape(turbidity_ACA) +
  #   tm_raster(style = 'cont')
  # 
  # bathy_ACA_Culebra = rast(here("data/Allen_Coral_Atlas/Culebrita-20240830215241/Bathymetry---composite-depth/bathymetry_0.tif"))
  # turbidity_ACA_Culebra = rast(here("data/Allen_Coral_Atlas/Culebrita-20240830215241/Turbidity-2019/turbidity-annual_0.tif"))
  # bathy_ACA_Culebra <- clamp(bathy_ACA_Culebra, lower = 0, upper = 50, values = TRUE)
  # tm_shape(bathy_ACA_Culebra) +
  #   tm_raster(style = 'cont')
  # turbidity_ACA_Culebra <- clamp(turbidity_ACA_Culebra, lower = 0, upper = 100, values = TRUE)
  # tm_shape(turbidity_ACA_Culebra) +
  #   tm_raster(style = 'cont')


