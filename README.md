# Calculating additional bioclimatic variables for species distribution modeling

##Bioclim variables : what are we working with? 

There is a suite of [19 bioclimatic variables](<https://www.worldclim.org/bioclim>) available from WorldClim, and also calculable in R with the biovars() function in the [dismo](<https://www.rdocumentation.org/packages/dismo/versions/1.1-4/topics/biovars>) package. 

These bioclimatic variables are commonly used in species distribution models, and include the following indices: 

BIO1 = Annual Mean Temperature
BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
BIO3 = Isothermality (BIO2/BIO7) (* 100)
BIO4 = Temperature Seasonality (standard deviation *100)
BIO5 = Max Temperature of Warmest Month
BIO6 = Min Temperature of Coldest Month
BIO7 = Temperature Annual Range (BIO5-BIO6)
BIO8 = Mean Temperature of Wettest Quarter
BIO9 = Mean Temperature of Driest Quarter
BIO10 = Mean Temperature of Warmest Quarter
BIO11 = Mean Temperature of Coldest Quarter
BIO12 = Annual Precipitation
BIO13 = Precipitation of Wettest Month
BIO14 = Precipitation of Driest Month
BIO15 = Precipitation Seasonality (Coefficient of Variation)
BIO16 = Precipitation of Wettest Quarter
BIO17 = Precipitation of Driest Quarter
BIO18 = Precipitation of Warmest Quarter
BIO19 = Precipitation of Coldest Quarter

However, I was developing a model of the distribution of a winter annual grass - which we know to be dependent on fall precipitation and spring temperature. Fall precipitation and spring temperature are unlikely to be captured by the above bioclimatic variables, which represent maxima and minima of temperature and precipitation. 

As a result, I developed a script to calculate additional bioclimatic variables that capture climate other than quarterly minima and maxima.

### What are the additional bioclimatic variables? 

To the list of 19 existing bioclimatic variables above, I added calculations of temperature and precipitation dynamics for the second-wettest / second-hottest quarter, etc. The additional variables are as follows: 

BIO20: Mean Temperature of Second-Wettest Quarter

BIO21: Mean Temperature of Second-Driest Quarter

BIO22: Mean Temperature of Second-Warmest Quarter

BIO23: Mean Temperature of Second-Coldest Quarter 

BIO24: Mean Precipitation of Second-Wettest Quarter

BIO25: Mean Precipitation of Second-Driest Quarter

BIO26: Mean Precipitation of Second-Warmest Quarter

BIO27: Mean Precipitation of Second-Coldest Quarter 



## Implementation

Download biovars_addtl.R to your location of choice. 

Load the libraries and the function. 

```R#load libraries
# load libraries
library(dismo)
library(raster)
library(rgdal)
library(here)

#load the function
source(here("scripts/biovars_addtl.R"))
```



Load your climate data of choice. 

```R
###### LOAD CLIMATE DATA

# these represent 30-year climate normals, downloaded from PRISM. These variables span CONUS at 4km resolution. 
# you can use any ppt, tmin, tmax variables that you are interested in.

# raster::brick() is necessary for the biovars command. 

ppt_normal <- brick(here("data/climate/prism/normal4km/pptnormal_brick.tif"))
tmin_normal <- brick(here("data/climate/prism/normal4km/tminnormal_brick.tif"))
tmax_normal <- brick(here("data/climate/prism/normal4km/tmaxnormal_brick.tif"))

```

biovars_addtl is implemented identically to dismo::biovars. 

```R
############ CALC ORIGINAL BIOCLIMATIC VARIABLES ##################

# use biovars function to calculation the 19 bioclimatic variables
# for more info on biovars, see: http://www.worldclim.org/bioclim

prism_biovars <- biovars(ppt_normal, tmin_normal, tmax_normal)

############ CALC ADDITIONAL BIOCLIMATIC VARIABLES ##################

# use adapted function to calculation an additional bioclimatic variables

#calculate
prism_biovars_addtl <- biovars_addtl(ppt_normal, tmin_normal, tmax_normal)

```



Crop, export, manipulate, and use them from here as you see fit!



## How it works

Working from the original [biovars](<https://github.com/cran/dismo/blob/master/R/biovars.R>) command, I manipulated the function to identify the second-highest (or second-lowest) value for a non-overlapping quarter. I.e., I want the months contined in the second warmest quarter to all be distinct from the months contained in the warmest quarter, etc.



The original biovars function looks at a vector of 12 values, with one value representing each month. biovars calculates rolling, 3-month windows. To identify, for example, the mean temperature of the wettest quarter, biovars starts with inputs of minimum temperature, maximum temperature, and precip. 

```R
tmin
tmax
prec
```

The window function calculates a rolling, 3-month mean of a variable of interest. 

```R
 window <- function(x)  { 
              lng <- length(x)
              x <- c(x,  x[1:3])
              m <- matrix(ncol=3, nrow=lng)
              for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
              apply(m, MARGIN=1, FUN=sum)
            }
```

The window function is used to create input temp and precip quarterly values— still vectors of 12 for each variable. 

```R
tavg <- (tmin + tmax) / 2
wet <- t(apply(prec, 1, window))
tmp <- t(apply(tavg, 1, window)) / 3
```

To identify, for example, the Mean Temperature of the Wettest Quarter, biovars identifies the quarterly value with the maximum precip, and extracts the temperature.

```R
# P8. Mean Temperature of Wettest Quarter 
              wetqrt <- cbind(1:nrow(p), as.integer(apply(wet, 1, which.max)))
              p[,8] <- tmp[wetqrt]
```



But, since I'm interested in the Mean Temperature of the **Second-Wettest Quarter** (that doesn't overlap with the Wettest Quarter), I need to a) identify the indices that overlap with the wettest quarter and b) calculate the second-highest value, after those indices are omitted. 

I do that using the following functions: 

```R
					#function that identifies which vectors to omit
            omit_ind <- function(x) {
              if(x > 2 & x <= 12) {
                c(x-2, x-1, x, x + 1, x+2)
              }
              else if(x == 2) {
                c(x-1, x, x + 1, x+2, 10 + x)
              }
              else if(x == 1) {
                c(x, x + 1, x+2, 10 + x, 11+x)
              }
            }
            
            #function that identifies second-highest or lowest value
            minmax2 <- function(x, minmax) {
              
              #na handling, if it makes it through
              if(length(which.max(x))== 0){
                NA}
              
              else if(length(which.max(x)) != 0){
                if(minmax == "max") {
                m1 <- as.numeric(which.max(x))
                omit <- omit_ind(m1)
                max(x[-c(omit)])
              }
                else if(minmax == "min") {
                m1 <- as.numeric(which.min(x))
                omit <- omit_ind(m1)
                min(x[-c(omit)])
              }
              }
            }
```

So, where the original biovars function identifies the minimum or the maximum, I identify the second-highest non-overlapping value. 

E.g., 

```R
# P20. Mean Temperature of Second-Wettest Quarter 
              wetqrt2 <- cbind(1:nrow(p), as.integer(apply(wet, 1, function(x) which(x==minmax2(x, "max"))[1])))
              p[,(20-19)] <- tmp[wetqrt2]
```



## Authors and acknowledgements

* **Lila Leatherman** 

Built directly on top of dismo::biovars . 



####dismo::biovars authorship: 

Author: Robert Hijmans

November 2009

License GPL3

based on:

MkBCvars.AML 

Author Robert Hijmans

January 2006  

Museum of Vertebrate Zoology, UC Berkeley

## License

TBD

