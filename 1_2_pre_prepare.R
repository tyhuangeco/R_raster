
library(raster)
library(zoo)
library(RNetCDF)

#0 functions to deal with NCnetcdf dataset####################################################

ncreshape <- function(filename) { 
  
  nc    <- open.nc(filename) 
  
  var_na <- strsplit(sub(".*./", "", filename), "_")[[1]][1]  # pick up the variable name from nc file i.e. gpp & pr 
  var    <- var.get.nc(nc, var_na)
  
  lat   <- rev(var.get.nc(nc, "lat"))                  # make lat from positive to negative 
  lon   <- var.get.nc(nc, "lon") 
  
  ref   <- att.get.nc(nc, "time", "units")             # get the starting date point of nc file i.e.: "days since 1861-01-01 00:00:00"
  tag   <- var.get.nc(nc, "time")                      # get the date based on ref point i.e.:  "15.5   45.0   74.5"  
  UTC   <- as.yearmon((utcal.nc(ref, tag, type="s")))  # convert ref + tag format into "Jan 1861" 
  Zval  <- as.Date(utcal.nc(ref, tag, type="s"))       # convert ref + tag format into "1861-01-16"
  
  Nx    <- dim(var)[1]
  Ny    <- dim(var)[2]
  var_n <- aperm(var[, Ny:1, ], c(2, 1, 3)) #reverse data along latitude dimension, and then transpose
  #the reason that we don't rotate longtitude here is because
  #we left it to rotate() function which also can change the 
  #longitude axses in raster. well, it is not a very good reasoon, but 
  #raster package is weird too ! 
  
  layers <- brick(var_n)
  extent(layers) <- c(min(lon), max(lon), 
                      min(lat), max(lat))
  projection(layers) <- "+proj=utm +zone=48 +datum=WGS84"
  names(layers)      <- UTC
  layers             <- setZ(layers, Zval, name='time')
  
  return(layers)
  
}


#1 make raw pr nc data in standard shape############################
files <- list.files(pattern=c("pr_Amon_Had", ".nc"), full.name=T)

m    <- lapply(files, ncreshape)
n    <- stack()
Zval <- as.Date(vector())

for (i in seq_along(m)) {
  n    <- addLayer(n, m[[i]]) 
  Zval <- c(Zval, getZ(m[[i]]))
}

pr_Had <- brick(n)
pr_Had <- setZ(pr_Had, Zval)

#2 identify time scoep of pr_Had output based on fluxnet data#######

spi_nc  <- brick("spi3_6_12_1deg_cru_ts_3_21_1949_2012.nc")
flx     <- brick("GPP_GL.nc")
Had     <- brick("GPP_Had.grd")
extent(Had) <- extent(pr_Had) 

t_flx    <- as.yearmon(getZ(flx))
t_pr_Had <- as.yearmon(getZ(pr_Had))
t_spi    <- as.yearmon(getZ(spi_nc))
t_com    <- as.yearmon(Reduce(intersect, list(t_pr_Had, t_flx, t_spi)))
idx_com  <- match(t_com, t_pr_Had)

pr_Had   <- pr_Had[[idx_com]]

#3 cut pr data by using gpp data####################################

for (n in 1: dim(pr_Had)[3]) { 
  
  pr_Had[[n]] <- mask(pr_Had[[n]], Had[[n]])
  
}

#4 convert the precipitation unit from kg/m^2/s to mm/month#########

begin    <- as.POSIXlt(getZ(pr_Had)[1])
st_point <- as.Date(begin) - begin[["mday"]] + 1

time_cal <- seq(st_point, length=length(names(pr_Had))+1, by='month')
dayofmon <- as.vector(diff(time_cal))


for (n in 1: dim(pr_Had)[3]) { 

   pr_Had[[n]] <- pr_Had[[n]] * (60*60*24*dayofmon[n])
  
}

writeRaster(pr_Had,  filename='pr_Had.grd') 













 
