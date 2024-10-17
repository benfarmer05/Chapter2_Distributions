-	Allen Coral Atlas 10m bath (check this one out for BVI especially!!):
o	https://allencoralatlas.org/methods/
-	If I ever want to reach out to NOAA about most recent bathymetry, Jen Kraus is probably the contact!

- Blue Habitats is intriguing, but maybe not needed for my applications: https://bluehabitats.org/





Correspondence with Tyler Smith on 20 June 2024:

Hi Tyler,
 
I wanted to return to the conversation about bathymetry, as I noticed some projection skew between the NCRMP sample frame grids and the CRM bathymetry I have been using. Dan shared with me ‘crm_usvi.tif’ years ago, and he said I should reach out to you to figure out what the source of the data was. From what I can tell online (https://www.ncei.noaa.gov/products/seafloor-mapping), it came from the Coastal Relief Model binned around the USVI (which I found using this grid extract tool: https://www.ncei.noaa.gov/maps/grid-extract/).
 
Does that sound right? At 1 arc-second resolution, I think that would be about 30 m, which looks right to me in many spots. Let me know if you’d like me to re-share Dan’s “VI_Shapes” folder on Box to take a closer look. What I may do in any case, is merge a CRM that spans the VI to PR, and from there merge in higher-resolution LIDAR sets (like this one: https://www.ngdc.noaa.gov/nos/H12001-H14000/H12271.html) where they are available.
 
One issue that may be related, or not, to what I’m seeing with skew is that the bathymetry sets tend to be in NAD83 CRS, while NCRMP is in WGS84.
 
Best regards, 




# General dashboard:
https://www.ncei.noaa.gov/products/seafloor-mapping

# 1 arc-second CRM:
https://www.ncei.noaa.gov/maps/grid-extract/
-	North: 19.032
-	South: 16.988
-	East: -64.246
-	West: -68.022

# Another dashboard (NCEI Bathymetric Data Viewer):
https://www.ncei.noaa.gov/maps/bathymetry/?layers=nos_hydro&minx=-65.5684&maxx=-64.4215&miny=17.8471&maxy=18.906
-	Above coordinates work as well for selecting Grid Extract here. Access to a LOT of bathymetry, up to 1 m resolution in coastal areas. If interested in domain all the way to the Dominican Republic:
o	North: 19.032
o	South: 16.988
o	East: -64.246
o	West: -68.022

# US Bathymetry Coverage and Gap Analysis:
https://iocm.noaa.gov/seabed-2030-bathymetry.html

# Another another dashboard (info from Laughlin Siceloff on September 5, 2024 <laughlin.siceloff@noaa.gov to bfarme8@lsu.edu>
#	- This one may be useful for filling in any gaps, or at least understanding where tiles in Jiangang's mosaic likely came from. It does not always match up exactly from what I can tell, but it's pretty close to what Jeremiah shared with me. Will take some real effort
#		to put together a final mosaic that is suitable for SDM work at high resolution
Hi Laughlin,
 
Yes, I think the latest of Jiangang’s work are what we received – it is a pretty large ESRI geodatabase with 4 tiles for PR and then one for STTSTJ and STX. I think Jeremiah just clipped it to 90 m depth, and we’re planning on limiting output to 50 m anyways. The stitching patterns and a few spots where land is “underwater” around Lameshur are spots where we may run into trouble with the SDMs, but the landmask at least should be easy to solve.
 
And thank you! I hadn’t found nowCOAST or Blue Topo yet. Would be great to hear more about the deep coral modeling later on. Our model will hopefully be able to project scleractinian cover to mesophotic depths, and maybe will be a supplemental additional resource.
 
And yes please, would be great to hear what you all are developing in the fall with mosaic products. I’ll check out the individual tiles / GitHub repo.
 
Best regards,
 
Ben Farmer (he/they)
PhD Candidate, NSF Graduate Research Fellow
Seascape Ecology Lab
Louisiana State University
(859) 475-3903
signature_2650033374
 
From: Laughlin Siceloff - NOAA Affiliate <laughlin.siceloff@noaa.gov>
Date: Thursday, September 5, 2024 at 16:04
To: Benjamin Farmer <bfarme8@lsu.edu>
Cc: Daniel Holstein <dholstein1@lsu.edu>
Subject: Re: Bathymetry

Hey guys,
 
I'm guessing you received Jiangang's mosaics from Jeremiah?  My understanding is that they've been a work in progress and he's made a few improvements but maybe there are more to make.  I haven't done any analysis with them although we use them for basic NCRMP needs.  He's incorporating all the latest NCCOS multibeam and Lidar that we have, so the raw ingredients are there at least.    
 
We are likely going to use OCM's Blue Topo mosaic tiles for our Caribbean deep coral modeling and derivate analyses.  We haven't started that process yet so I can't report back, but should happen this fall.  If you're not familiar with this product, check it out here: https://nowcoast.noaa.gov/
 
If you play with the symbology you can see the individual tiles, which are varying resolution but are continually updated with new bathy.  It's only set up on this web tool to DL individual tiles of interest but they have github coding to download in bulk.  NCCOS will be making various mosaic products with these tiles in the fall and I can keep you posted on that progress.    
 
Thanks,
L
 
On Thu, Sep 5, 2024 at 12:23 PM Benjamin Farmer <bfarme8@lsu.edu> wrote:
Hi Laughlin,
 
Was great diving with you at Flower Gardens! In the time since we talked about bathymetry a bit on the cruise, I received a geodatabase from Jeremiah Blondeau that has been really helpful. It just has a few mosaic/stitching patterns in it that have a downstream effect on derived bathymetry variables (slope, slope of slope, aspect, etc.).
 
I know you had mentioned working on a similar SDM approach as our lab is, and I am wondering if you have run into any issues like the above. I already reached out to Jeremiah, but attached are examples that I sent him as well, from west and south of St Thomas. If you might be able to share the bathymetry you’re working with (assuming it is different than what I’ve got), it would be great to compare. I am also happy to share anything I am working on if it would helpful!
 
Best regards,
 
Ben Farmer (he/they)
PhD Candidate, NSF Graduate Research Fellow
Seascape Ecology Lab
Louisiana State University
(859) 475-3903
signature_3914545081

 
--
Laughlin Siceloff
CSS-Inc., under contract to NOAA
Biogeography Branch
National Centers for Coastal Ocean Science
NOAA National Ocean Service
Cell: 919-225-7148