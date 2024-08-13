Documentation of what I did to integrate the shapefiles from online

TNC:
https://sites.google.com/view/caribbean-marine-maps
https://www.nature.org/en-us/newsroom/caribbean-marine-maps/

-	Brought in shapefiles for PR, DR, BVI, USVI. Extracted the first 5 raster values (Reef Crest, Fore Reef, Back Reef, Coral/Algae, and Spur and Groove Reef). According to Supplementary Material from Schill et al. (2021) in Remote Sensing. Some coral are in the Hardbottom Algae areas, as well as in Seagrass and Boulders and Rocks, however including these layers would have vastly overestimated shallow coral coverage since corals are generally quite sparse in these areas.
-	Methods:
o	Extract only raster values I want, to a new raster (https://gis.stackexchange.com/questions/49800/extracting-raster-values-and-creating-new-raster-in-qgis)
o	Polygonize that raster (now “reefs” are 1s and everything else is 0s)
o	Extract just the 1s to a new polygon (Extract by Attributes) (https://gis.stackexchange.com/questions/93237/eliminate-0-values-from-raster)

UNEP-WCMC, WorldFish Centre, WRI, TNC:
https://data.unep-wcmc.org/pdfs/1/WCMC_008_Global_Distribution_of_Coral_Reefs.pdf?1617121809
https://data.unep-wcmc.org/datasets/1


Insular shelf south of St Thomas and St John:
https://coastalscience.noaa.gov/project/habitat-map-insular-shelf-south-of-st-thomas-st-john/
https://www.fisheries.noaa.gov/inport/item/50381


Habitat drafts from Dan:
These appear to be a merged file from:
https://coastalscience.noaa.gov/project/habitat-mapping-us-virgin-islands/
https://products.coastalscience.noaa.gov/collections/benthic/e95usvi_pr/
https://products.coastalscience.noaa.gov/collections/benthic/e93stcroix/
https://coastalscience.noaa.gov/project/benthic-habitat-mapping-buck-island-reef-national-monument-st-croix/
https://data.doi.gov/dataset/digital-geologic-map-of-st-croix-and-buck-island-reef-national-monument-u-s-virgin-islands
-	and Hind Bank MCD map (unclear where this is housed online)

SDM from Katharine Egan:
Shared through email with me bfarme8@lsu.edu
NOAA BVI:
https://www.fisheries.noaa.gov/inport/item/40578
https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0049955

PR:
https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0217139
https://coastalscience.noaa.gov/news/2019-puerto-rico-coral-reef-monitoring-data-now-available/


DCRMP:
Want to use this?


STT+STX sampleframe extent:
272699.9688,348350.0313,1950849.75,2037350.0000 [EPSG:26920]


Lo-res net:
270099.9688,376375.1563,2005498.5672,2065950.0000 [EPSG:26920]

Net over BVI:
270100.0000,448700.0000,2005500.0000,2137950.0000 [EPSG:26920]

270100.0000,378100.0000,2005500.0000,2066700.0000 [EPSG:26920]

270100.0000,378100.0000,2065950.0000,2127150.0000 [EPSG:26920]

376700.0000,401900.0000,2005500.0000,2065950.0000 [EPSG:26920]




Possible future sources:

-	NOAA NCOSS: PR/USVI Prioritization Mapping Inventory
o	https://noaa.maps.arcgis.com/home/webmap/viewer.html?webmap=3849f7c4a0ca4fe8a54a4bc853da86b5
-	Check out: Coral Health Survey Basic mapping (SCTLD) (2019-2021)
o	This includes the DR!
o	https://us-shallow-coral-reef-mapping-priorities-noaa.hub.arcgis.com/pages/results

