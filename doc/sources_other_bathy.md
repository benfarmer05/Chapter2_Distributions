-	Allen Coral Atlas 10m bath (check this one out for BVI especially!!):
o	https://allencoralatlas.org/methods/
-	If I ever want to reach out to NOAA about most recent bathymetry, Jen Kraus is probably the contact!







Correspondence with Tyler Smith on 20 June 2024:

Hi Tyler,
 
I wanted to return to the conversation about bathymetry, as I noticed some projection skew between the NCRMP sample frame grids and the CRM bathymetry I have been using. Dan shared with me ‘crm_usvi.tif’ years ago, and he said I should reach out to you to figure out what the source of the data was. From what I can tell online (https://www.ncei.noaa.gov/products/seafloor-mapping), it came from the Coastal Relief Model binned around the USVI (which I found using this grid extract tool: https://www.ncei.noaa.gov/maps/grid-extract/).
 
Does that sound right? At 1 arc-second resolution, I think that would be about 30 m, which looks right to me in many spots. Let me know if you’d like me to re-share Dan’s “VI_Shapes” folder on Box to take a closer look. What I may do in any case, is merge a CRM that spans the VI to PR, and from there merge in higher-resolution LIDAR sets (like this one: https://www.ngdc.noaa.gov/nos/H12001-H14000/H12271.html) where they are available.
 
One issue that may be related, or not, to what I’m seeing with skew is that the bathymetry sets tend to be in NAD83 CRS, while NCRMP is in WGS84.
 
Best regards, 




General dashboard:
https://www.ncei.noaa.gov/products/seafloor-mapping

1 arc-second CRM:
https://www.ncei.noaa.gov/maps/grid-extract/
-	North: 19.032
-	South: 16.988
-	East: -64.246
-	West: -68.022

Another dashboard (NCEI Bathymetric Data Viewer):
https://www.ncei.noaa.gov/maps/bathymetry/?layers=nos_hydro&minx=-65.5684&maxx=-64.4215&miny=17.8471&maxy=18.906
-	Above coordinates work as well for selecting Grid Extract here. Access to a LOT of bathymetry, up to 1 m resolution in coastal areas. If interested in domain all the way to the Dominican Republic:
o	North: 19.032
o	South: 16.988
o	East: -64.246
o	West: -68.022

US Bathymetry Coverage and Gap Analysis:
https://iocm.noaa.gov/seabed-2030-bathymetry.html
