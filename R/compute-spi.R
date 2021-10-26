library(tidyverse)
library(sf)
library(raster)
library(ncdf4)

#define ancellary functions
wherenearest   <- function(val,matrix){
  dist  = abs(matrix-val)
  index = which.min(dist)
  return( index )
}

roi = st_read('/home/zhoylman/Desktop/RussianRiver/RussianRiver.shp') %>%
  st_transform(4326)

points = raster('/home/zhoylman/soil-moisture-validation-data/gridmet/precip_gridmet.nc')[[1]] %>%
  crop(., roi) %>%
  mask(., roi) %>%
  rasterToPoints() %>%
  as.data.frame() %>%
  as_tibble() %>%
  dplyr::select(x,y) 
  
nc_precip  = nc_open('/home/zhoylman/soil-moisture-validation-data/gridmet/precip_gridmet.nc')
lon.precip = ncvar_get(nc_precip,varid='lon')
lat.precip = ncvar_get(nc_precip,varid='lat')

precip_time = read_csv('/home/zhoylman/soil-moisture-validation-data/gridmet/precip_gridmet_time.csv')

#exctract topo vals from disk
data_precip = apply(points,MARGIN=1,FUN=function(x){  
  return(ncvar_get(nc_precip,varid='precip',start=c(wherenearest(x['x'] %>% as.numeric(),lon.precip),wherenearest(x['y'] %>% as.numeric(),lat.precip),1),count=c(1,1,-1))) 
}) %>%
  as_tibble() %>%
  mutate(time = precip_time$datetime) %>%
  pivot_longer(cols = -c(time))

nc_close(nc_precip)
