install.packages("rmarkdown")
install.packages("knitr") 
library(rmarkdown)
library(knitr)
# load the ncdf4 package
library(ncdf4)
library(ncdf.tools)

# set path and filename
ncpath <- "C:/geog490/data/"
ncname <- "icec.mon.mean"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "icec" 

# open a netCDF file
ncin <- nc_open(ncfname)
print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
head(lon)
lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
head(lat)
print(c(nlon,nlat))

# get time
time <- ncvar_get(ncin,"time")
time
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt
tunits

Years <- seq(1850, (2019-0.008333), by = 0.08333)
head(Years)

# convert time -- split the time units string into fields
print(tunits)
tustr <- strsplit(tunits$value, " ")
ptime <- convertDateNcdf2R(time, unlist(tustr)[1], origin = 
                             as.POSIXct(unlist(tustr)[3], tz = "UTC"), time.format = "%Y-%m-%d")
head(time); tail(time)
head(ptime); tail(ptime) 

# get ice data
var_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(var_array)

# close the netCDF file
nc_close(ncin)


# set up for calculating area averages--define two functions

# radians function (convert degrees to radians)
radians <- function(x) {
  radians <- (2/360) * pi * x
}


grid_cell_area <- function(lat, lon) {
  nlat <- length(lat)
  nlon <- length(lon)
  radius <- 6371.0087714
  area <- rep(NA, nlat)
  dlon <- abs(lon[1] - lon[2]) # longitudes assumed to be equally spaced
  # check direction of latitudes
  if (lat[1] <= 0.0) {
    #print("increasing")
    dlat1 <- lat[1] - -90.0
    for (k in (1:nlat)) {
      if (k < nlat) {
        dlat2 <- (lat[k+1] - lat[k]) /2
      } else {
        dlat2 <- 90.0 - lat[nlat]
      }
      area[k] <- (radius ^ 2) * (radians(dlon) * 
                                   (sin(radians(lat[k] + dlat2)) - sin(radians(lat[k] - dlat1))) )
      dlat1 <- dlat2
    }
    
  } else {
    #print("decreasing")
    dlat1 <- 90.0 - lat[1]
    for (k in (1:nlat)) {
      if (k < nlat) {
        dlat2 <- (lat[k] - lat[k+1]) /2
      } else {
        dlat2 <- lat[nlat] - -90.0
      }
      area[k] <- (radius ^ 2) * (radians(dlon) * 
                                   (sin(radians(lat[k] + dlat1)) - sin(radians(lat[k] - dlat2))) )
      dlat1 <- dlat2
    }
  }
  return(area)
}

# get grid cell areas (in km^2)
area <- grid_cell_area(lat, lon)
area

# get area-weighted average of temperature in Northern Hemisphere
area_ave <- rep(0, nt)

# loop over times
for (n in (1:nt)) {
  # loop over longitudes and latitudes
  sum <- 0.0    # sum of weighted values
  wsum <- 0.0   # sum of weights (areas)
  for (j in (1:nlon)) {
    for (k in (1:nlat)) {
      if (!is.na(var_array[j,k,n])) {
        # grid point isn't missing (i.e. not ocean), so accumulate values
        sum <- sum + var_array[j,k,n] * area[k]
        wsum <- wsum + area[k]
      }
    }
  }
  # calculate weighted average
  area_ave[n] <- sum / wsum
}
print(area_ave)

# get total land area above 0 degC (in Northern Hemisphere)
Ice_Concentration_Totals <- rep(0, nt)

# loop over times
for (n in (1:nt)) {
  # loop over longitudes and latitudes
  areasum <- 0.0   # sum of areas
  for (j in (1:nlon)) {
    for (k in (1:nlat)) {
      if (!is.na(var_array[j,k,n])) {
        # grid point isn't missing (i.e. not ocean), so accumulate areas if above 0 degC
        Ice_Concentration_Totals[n] <- Ice_Concentration_Totals[n] + area[k]*var_array[j,k,n]
      }
    }
  }
}
print(Ice_Concentration_Totals)

yrs <- 169
mnth <- 12 
ann_ave <- rep(0, yrs)
for (n in (1: yrs)) {
  for (m in (1:mnth)) {
    ann_ave[n] <- ann_ave[n] + Ice_Concentration_Totals[(n-1) * mnth + m]
  }
  ann_ave[n] = ann_ave[n]/mnth
}

yr <- seq(1850, 2018, by = 1)

plot(yr, ann_ave, type="l", main = "Annual Averages of Ice Concentration",
     xlab = "Year", ylab = "Average Ice Concentration")
