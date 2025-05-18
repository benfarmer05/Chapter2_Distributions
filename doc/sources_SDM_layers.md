-	CARICOOS Curents
o	https://www.caricoos.org/currents/forecast/FVCOM/VIHR/Currents
-	CARICOOS Waves
o	https://www.caricoos.org/waves/forecast/SWAN/PRVI/hsig
-	Also just CARICOOS in general

- NOAA SWAN
- Go to here: https://www.caricoos.org/data-download
	- then go to ERDDAP, search 'wave', and look at the SWAN models. only 1-km resolution though

-	Copernicus data using esa snap
o	https://step.esa.int/main/download/snap-download/



-	NASA SeaDAS Ocean Color
o	https://seadas.gsfc.nasa.gov/


-	Mary Stoll (https://uwpcc.ocean.washington.edu/person/Mary_Margaret_Stoll) did a NOAA PMEL seminar on species distribution models!
o	CMEMS global ocean 3d BBP / CHL-a
	https://data.marine.copernicus.eu/product/MULTIOBS_GLO_BIO_BGC_3D_REP_015_010/description
o	BGC Argo
	https://www.aoml.noaa.gov/biogeochemical-argo-program/



Consider re: SST & chl-a:
- Archer et al. 2025, Nature, show amazing graphic (750-m resolution) of chl-a and SST. If I can wrangle these for USVI/PR, would be very nice. Then I can get current direction/speed if I want, and 
- The source was VIIRS satellite data, some resources for that:
	- https://coralreefwatch.noaa.gov/product/oc/ (Daily 750m VIIRS Satellite Ocean Color Monitoring).
		- problem is this one cuts off Tortola
	- https://lpdaac.usgs.gov/data/get-started-data/collection-overview/missions/s-npp-nasa-viirs-overview/
	- https://oceandata.sci.gsfc.nasa.gov/ (NASA Ocean Color database; has VIIRS in it)
		- Can check out NASA Earthdata, AppEEARS, etc. for visualizing and downloading data for specific areal extents. But the resulting data might be super hard to work with in R / MATLAB :(