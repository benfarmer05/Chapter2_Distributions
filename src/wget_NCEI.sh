
#useful for reading in all files from somewhere like the NCRMP NCEI page

#wget -r -np -nH --cut-dirs=8 -R "index.html*" -P /Volumes/dissertation_Box_backup/sample_frames/NCEI_inport/PR2021 https://www.ncei.noaa.gov/data/oceans/archive/arc0209/0272313/1.1/data/0-data/NCRMP_PR_2021_Benthics/Data_Sets/Sample_Frames/



#!/bin/bash

# Base URL
BASE_URL="https://www.fisheries.noaa.gov/inport/item/70663"

# Destination directory
DEST_DIR="/Volumes/dissertation_Box_backup/sample_frames/NCEI_inport/PR2021"

# Run wget to download all files
wget -r -np -nH --cut-dirs=8 -R "index.html*" -P $DEST_DIR $BASE_URL
