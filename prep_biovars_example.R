## Prep PRISM climate data into biovars

#load libraries
library(tidyverse)
library(dismo)
library(raster)
library(rgdal)
library(here)

#load the function
source(here("scripts/biovars_addtl.R"))

###### LOAD CLIMATE DATA

# these represent 30-year climate normals, downloaded from PRISM. These variables span CONUS at 4km resolution. 
# you can use any ppt, tmin, tmax variables that you are interested in.

# raster::brick() is necessary for the biovars command. 

ppt_normal <- brick(here("data/climate/prism/normal4km/pptnormal_brick.tif"))
tmin_normal <- brick(here("data/climate/prism/normal4km/tminnormal_brick.tif"))
tmax_normal <- brick(here("data/climate/prism/normal4km/tmaxnormal_brick.tif"))

############ CALC ORIGINAL BIOCLIMATIC VARIABLES ##################

# use biovars function to calculation the 19 bioclimatic variables
# for more info on biovars, see: http://www.worldclim.org/bioclim

prism_biovars <- biovars(ppt_normal, tmin_normal, tmax_normal)

#inspect - make sure they look okay
str(prism_biovars)
plot(prism_biovars)

#export
writeRaster(prism_biovars, file =here("data/climate/prism/biovars/prism_biovars_US.tif"), options = "INTERLEAVE=BAND", overwrite=TRUE)

#make sure file can be read back in
#prism_biovars <- brick(paste0(here("data/climate/prism/biovars/prism_biovars_US_"), res, ".tif"))


############ CALC ADDITIONAL BIOCLIMATIC VARIABLES ##################

# use adapted function to calculation an additional bioclimatic variables
# additional bioclim variables include: 
# BIO20: Mean Temperature of Second-Wettest Quarter
# BIO21: Mean Temperature of Second-Driest Quarter
# BIO22: Mean Temperature of Second-Warmest Quarter
# BIO23: Mean Temperature of Second-Coldest Quarter 
# BIO24: Mean Precipitation of Second-Wettest Quarter
# BIO25: Mean Precipitation of Second-Driest Quarter
# BIO26: Mean Precipitation of Second-Warmest Quarter
# BIO27: Mean Precipitation of Second-Coldest Quarter 

#calculate
prism_biovars_addtl <- biovars_addtl(ppt_normal, tmin_normal, tmax_normal)

# #inspect
# prism_biovars_addtl
# plot(prism_biovars_addtl)

#export
writeRaster(prism_biovars_addtl, file =here("data/climate/prism/biovars/prism_biovars_addtl_US.tif"), options = "INTERLEAVE=BAND", overwrite=TRUE)

# make sure file can be read back in 
prism_biovars_addtl <- brick(here("data/climate/prism/biovars/prism_biovars_addtl_US_.tif"))
