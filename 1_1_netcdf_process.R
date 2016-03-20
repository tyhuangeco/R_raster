
library(RNetCDF)
library(raster)
library(zoo) 

ncreshape <- function(filename) { 
  
  nc    <- open.nc(filename) 
  
  var_na <- strsplit(sub(".*./", "", filename), "_")[[1]][1]  # pick up the variable name from nc file i.e. gpp & pr 
  var    <- var.get.nc(nc, var_na)
  
  lat   <- rev(var.get.nc(nc, "lat"))                         # make lat from positive to negative 
  lon   <- var.get.nc(nc, "lon") 
  
  ref   <- att.get.nc(nc, "time", "units")             # get the starting date point of nc file i.e.: "days since 1861-01-01 00:00:00"
  tag   <- var.get.nc(nc, "time")                      # get the date based on ref point i.e.:  "15.5   45.0   74.5"  
  cal   <- as.yearmon((utcal.nc(ref, tag, type="s")))  # convert ref + tag format into things like "Jan 1861" (package use: RNetCDF, zoo)
  
  Nz    <- as.Date(utcal.nc(ref, tag, type="s"))       # convert ref + tag format into "1861-01-16"
  Nx    <- dim(var)[1]
  Ny    <- dim(var)[2]
  var_n <- aperm(var[, Ny:1, ], c(2, 1, 3))            # reverse data along latitude dimension, and then transpose
                                                       # the reason that we don't rotate longtitude here is because
                                                       # we left it to rotate() function which also can change the 
                                                       # longitude axses in raster. well, it is not a very good reasoon, but 
                                                       # raster package is weird too ! 
  
  layers <- brick(var_n)
  extent(layers) <- c(min(lon), max(lon), 
                      min(lat), max(lat))
  projection(layers) <- "+proj=utm +zone=48 +datum=WGS84"
  names(layers)      <- cal
  layers             <- setZ(layers, Nz, name='time')
  
  return(layers)
  
}

#CHAPTER 1 GPP data process######################################################

#1 extract GPP array from NetCDF data and build brick objects, extract time vector
#and build up Date objects for further application################################

files <- list.files(pattern=c("gpp_Lmon", ".nc"), full.name=T)

m    <- lapply(files, ncreshape)
n    <- stack()
Zval <- as.Date(vector())

for (i in seq_along(m)) {
  n    <- addLayer(n, m[[i]]) 
  Zval <- c(Zval, getZ(m[[i]]))
}

gpp <- brick(n)
gpp <- setZ(gpp, Zval) 
#gpp <- rotate(gpp)

#2 identify common time scope of CMIP5 data output and FLUXnet data###############

spi_nc  <- brick("spi3_6_12_1deg_cru_ts_3_21_1949_2012.nc")
flx     <- brick("GPP_GL.nc")

t_flx   <- as.yearmon(getZ(flx))
t_gpp   <- as.yearmon(getZ(gpp))
t_spi   <- as.yearmon(getZ(spi_nc))
t_com   <- as.yearmon(Reduce(intersect, list(t_gpp, t_flx, t_spi)))
idx_com <- match(t_com, t_gpp)

gpp  <- gpp[[idx_com]]
Zval <- Zval[idx_com]

gpp <- setZ(gpp, Zval) 

#3 some CMIP5 GPP dataset use 0 instead of NA to represent ocean grid cell, here provide
#a way to convert ocean cell into NA value##############################################

carbon  <- getValues(gpp) 
carbon1 <- getValues(gpp[[2:13]])
na_ind  <- which(apply(carbon1, 1, sum) == 0)
carbon[na_ind, ] <- NA

gpp <- setValues(gpp, carbon) 

#CHAPTER 2 Precipitation data process###############################################

#1 extract precipitation array from NetCDF data and build brick objects, extract time vector
#and build up Date objects for further application##########################################

files <- list.files(pattern=c("pr_Amon", ".nc"), full.name=T)

m    <- lapply(files, ncreshape)
n    <- stack()
Zval <- as.Date(vector())

for (i in seq_along(m)) {
  n    <- addLayer(n, m[[i]]) 
  Zval <- c(Zval, getZ(m[[i]]))
}

pr <- brick(n)
pr <- setZ(pr, Zval)
#pr <- rotate(pr)

#2 identify common time scope of CMIP5 precipitation data############################

pr   <- pr[[idx_com]]
Zval <- Zval[idx_com]

pr <- setZ(pr, Zval) 

#3 select pr data in terreatrial region by using gpp data#############################

for (n in 1: dim(pr)[3]) { 
  
  pr[[n]] <- mask(pr[[n]], gpp[[n]])
  
}

#4 ###################################################################################



time_cal <- seq(as.Date("1982-01-01"), length=length(time)+1, by='month')
dayofmon <- as.vector(diff(time_cal))

for (n in 1: dim(pr)[3]) { 
  pr_dsb[[n]] <- pr_dsb[[n]] * (60*60*24*dayofmon[n])
}































