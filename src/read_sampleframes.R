
#useful for reading in .rda files from somewhere like the NCRMP GitHub

library(sp)
library(here)

# List all .rda files in the current directory
rda_files <- list.files(pattern = "\\.rda$")

# Load each .rda file
for (file in rda_files) {
  load(file)
}

# Optional: Print the names of loaded objects
print(ls())
