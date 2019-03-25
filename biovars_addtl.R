# Author: Lila Leatherman
# Jan 2019
# License TBS
#
# based on:
# Author: Robert Hijmans
# November 2009
# License GPL3
#
# based on:
# MkBCvars.AML 
# Author Robert Hijmans
# January 2006  
# Museum of Vertebrate Zoology, UC Berkeley
#
# Version 1.0
#
# function to create additional BIOCLIM variables from 
# monthly T-min, T-max, and Precipitation data
#
# standard bioclim variables include: 
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (P2/P7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (P5-P6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter
# 
# additional bioclim variables include: 
# BIO20: Mean Temperature of Second-Wettest Quarter
# BIO21: Mean Temperature of Second-Driest Quarter
# BIO22: Mean Temperature of Second-Warmest Quarter
# BIO23: Mean Temperature of Second-Coldest Quarter 
# BIO24: Mean Precipitation of Second-Wettest Quarter
# BIO25: Mean Precipitation of Second-Driest Quarter
# BIO26: Mean Precipitation of Second-Warmest Quarter
# BIO27: Mean Precipitation of Second-Coldest Quarter 
# 
# The aim of these additional variables is to represent spring or fall climate, 
# which are important for the phenology and distribution of annual grasses. 
# 
# The original summary Bioclimatic variables are after:
#   Nix, 1986. A biogeographic analysis of Australian elapid snakes. In: R. Longmore (ed.).
#      Atlas of elapid snakes of Australia. Australian Flora and Fauna Series 7.
#      Australian Government Publishing Service, Canberra.
#
# and Expanded following the ANUCLIM manual
# 

## I have retained the naming of the additional biovars to align with the orignal ones. 
## However, the raster stack will only contain 8 layers.

#if (!isGeneric("biovars_addtl")) {
  setGeneric("biovars_addtl", function(prec, tmin, tmax, ...)
    standardGeneric("biovars_addtl"))
#}	


setMethod('biovars_addtl', signature(prec='vector', tmin='vector', tmax='vector'), 
          function(prec, tmin, tmax) {
            biovars_addtl(t(as.matrix(prec)), t(as.matrix(tmin)), t(as.matrix(tmax)))
          }
)


setMethod('biovars_addtl', signature(prec='Raster', tmin='Raster', tmax='Raster'), 
          function(prec, tmin, tmax, filename='', progress='', ...) {
            
            if (nlayers(prec) != 12) stop('nlayers(prec) is not 12')
            if (nlayers(tmin) != 12) stop('nlayers(tmin) is not 12')
            if (nlayers(tmax) != 12) stop('nlayers(tmax) is not 12')
            
            
            # temporary fix to avoid warning
            compareRaster(prec, tmin, tmax)
            
            out <- brick(prec, values=FALSE)
            out@data@nlayers <- as.integer(8)
            names(out) <- paste('bio', 20:27, sep="")
            
            filename <- trim(filename)
            if (!canProcessInMemory(out, 8)) {
              if (filename == '') { 
                filename <- rasterTmpFile()
              }
            } 
            if (filename == "") {
              v <- matrix(nrow=ncell(out), ncol=8)
            } else {
              out <- writeStart(out, filename, ...)
            }	
            
            tr <- blockSize(out, n=nlayers(out)+36)
            pb <- pbCreate(tr$n, ...)	
            for (i in 1:tr$n) {
              prc <- getValues(prec, tr$row[i], tr$nrows[i])
              tmn <- getValues(tmin, tr$row[i], tr$nrows[i])
              tmx <- getValues(tmax, tr$row[i], tr$nrows[i])
              p <- biovars_addtl(prc, tmn, tmx)
              if (filename != "") {
                out <- writeValues(out, p, tr$row[i])
              } else {
                start <- (tr$row[i]-1) * out@ncols + 1
                end <- (tr$row[i]+tr$nrows[i]-1) * out@ncols
                v[start:end,] <- p
              }
            }
            
            if (filename == "") {
              out <- setValues(out, v)
            } else {
              out <- writeStop(out)
            }	
            return(out)
          }
)




setMethod('biovars_addtl', signature(prec='matrix', tmin='matrix', tmax='matrix'), 
          function(prec, tmin, tmax) {
            
            if (nrow(prec) != nrow(tmin) | nrow(tmin) != nrow(tmax) ) {
              stop('prec, tmin and tmax should have same length')
            }
            
            if (ncol(prec) != ncol(tmin) | ncol(tmin) != ncol(tmax) ) {
              stop('prec, tmin and tmax should have same number of variables (columns)')
            }
            
            # can't have missing values in a row
            nas <- apply(prec, 1, function(x){ any(is.na(x)) } )
            nas <- nas | apply(tmin, 1, function(x){ any(is.na(x)) } )
            nas <- nas | apply(tmax, 1, function(x){ any(is.na(x)) } )
            p <- matrix(nrow=nrow(prec), ncol=8)
            colnames(p) = paste('bio', 20:27, sep='')
            if (all(nas)) { return(p) }
            
            prec[nas,] <- NA
            tmin[nas,] <- NA
            tmax[nas,] <- NA
            
            window <- function(x)  { 
              lng <- length(x)
              x <- c(x,  x[1:3])
              m <- matrix(ncol=3, nrow=lng)
              for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
              apply(m, MARGIN=1, FUN=sum)
            }
            
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
            
            tavg <- (tmin + tmax) / 2
            
            # # P1. Annual Mean Temperature 
            # p[,1] <- apply(tavg,1,mean)
            # # P2. Mean Diurnal Range(Mean(period max-min)) 
            # p[,2] <- apply(tmax-tmin, 1, mean)
            # # P4. Temperature Seasonality (standard deviation) 
            # p[,4] <- 100 * apply(tavg, 1, sd)
            # # P5. Max Temperature of Warmest Period 
            # p[,5] <- apply(tmax,1, max)
            # # P6. Min Temperature of Coldest Period 
            # p[,6] <- apply(tmin, 1, min)
            # # P7. Temperature Annual Range (P5-P6) 
            # p[,7] <- p[,5] - p[,6]
            # # P3. Isothermality (P2 / P7) 
            # p[,3] <- 100 * p[,2] / p[,7]
            # # P12. Annual Precipitation 
            # p[,12] <- apply(prec, 1, sum)
            # # P13. Precipitation of Wettest Period 
            # p[,13] <-  apply(prec, 1, max)
            # # P14. Precipitation of Driest Period 
            # p[,14] <-  apply(prec, 1, min)
            # # P15. Precipitation Seasonality(Coefficient of Variation) 
            # # the "1 +" is to avoid strange CVs for areas where mean rainfaill is < 1)
            # p[,15] <- apply(prec+1, 1, cv)
            
            # precip by quarter (3 months)		
            wet <- t(apply(prec, 1, window))
            # # P16. Precipitation of Wettest Quarter 
            # p[,16] <- apply(wet, 1, max)
            # # P17. Precipitation of Driest Quarter 
            # p[,17] <- apply(wet, 1, min)
            
            # P24. Precipitation of Second-Wettest Quarter
            p[,(24-19)] <- apply(wet, 1, function(x) minmax2(x, "max"))
            # P25. Precipitation of Second-Driest Quarter
            p[,(25-19)] <- apply(wet, 1, function(x) minmax2(x, "min"))
            
            tmp <- t(apply(tavg, 1, window)) / 3
            
            if (all(is.na(wet))) {
              p[,(20-19)] <- NA		
              p[,(21-19)] <- NA		
            } else {
              # # P8. Mean Temperature of Wettest Quarter 
              # wetqrt <- cbind(1:nrow(p), as.integer(apply(wet, 1, which.max)))
              # p[,8] <- tmp[wetqrt]
              # # P9. Mean Temperature of Driest Quarter 
              # dryqrt <- cbind(1:nrow(p), as.integer(apply(wet, 1, which.min)))
              # p[,9] <- tmp[dryqrt]
              
              # P20. Mean Temperature of Second-Wettest Quarter 
              wetqrt2 <- cbind(1:nrow(p), as.integer(apply(wet, 1, function(x) which(x==minmax2(x, "max"))[1])))
              p[,(20-19)] <- tmp[wetqrt2]
              # P21. Mean Temperature of Second-Driest Quarter 
              dryqrt2 <- cbind(1:nrow(p), as.integer(apply(wet, 1, function(x) which(x==minmax2(x, "min"))[1])))
              p[,(21-19)] <- tmp[dryqrt2]
            }
            # # P10 Mean Temperature of Warmest Quarter 
            # p[,10] <- apply(tmp, 1, max)
            # 
            # # P11 Mean Temperature of Coldest Quarter
            # p[,11] <- apply(tmp, 1, min) 
            
            # P22. Mean Temperature of Second-Warmest Quarter 
            p[,(22-19)] <- apply(tmp, 1, function(x) minmax2(x, "max"))
            
            # P23 Mean Temperature of Second-Coldest Quarter
            p[,(23-19)] <- apply(tmp, 1, function(x) minmax2(x, "min")) 
            
            if (all(is.na(tmp))) {
              p[,(26-19)] <- NA		
              p[,(27-19)] <- NA
            } else {
              # # P18. Precipitation of Warmest Quarter 
              # hot <- cbind(1:nrow(p), as.integer(apply(tmp, 1, which.max)))
              # p[,18] <- wet[hot]
              # # P19. Precipitation of Coldest Quarter 
              # cold <- cbind(1:nrow(p), as.integer(apply(tmp, 1, which.min)))
              # p[,19] <- wet[cold]
              
              # P26. Precipitation of Second-Warmest Quarter 
              hot2 <- cbind(1:nrow(p), as.integer(apply(tmp, 1, function(x) which(x==minmax2(x, "max"))[1])))
              p[,(26-19)] <- wet[hot2]
              # P27. Precipitation of Second-Coldest Quarter 
              cold2 <- cbind(1:nrow(p), as.integer(apply(tmp, 1, function(x) which(x==minmax2(x, "min"))[1])))
              p[,(27-19)] <- wet[cold2]
            }
            
            return(p)	
          }
)

