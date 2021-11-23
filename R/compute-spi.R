library(tidyverse)
library(sf)
library(raster)
library(ncdf4)
library(lubridate)
library(foreach)
library(doParallel)

#define ancellary functions
wherenearest   <- function(val,matrix){
  dist  = abs(matrix-val)
  index = which.min(dist)
  return( index )
}

#source all relevant drought functions
source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/master/processing/ancillary-functions/R/drought-functions.R')

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

sum = data_precip %>%
  group_by(time) %>%
  summarise(sum = sum(value))

write_csv(sum, file = '/home/zhoylman/temp/russian_precip_sum_mm_232px.csv')

#add in data length conditional (30 years)
compute_spi = function(index, data, timescale){
  data = data %>%
    mutate(yday = yday(time)) %>%
    #slice out data older than index
    slice(1:index)
  
  #compute indexes for time breaks
  first_date_breaks = which(data$yday == data$yday[index])
  second_date_breaks = first_date_breaks-(timescale-1)
  
  #if there are negative indexes remove last year (incomplete data range)
  #change this to remove all indexes from both vectors that are negative
  if(!all(second_date_breaks < 0)){
    pos_index = which(second_date_breaks > 0)
    first_date_breaks = first_date_breaks[c(pos_index)]
    second_date_breaks = second_date_breaks[c(pos_index)]
  }
  
  #create slice vectors and group by vectors
  for(j in 1:length(first_date_breaks)){
    if(j == 1){
      slice_vec = seq(second_date_breaks[j],first_date_breaks[j], by = 1)
      group_by_vec = rep(j,(first_date_breaks[j] - second_date_breaks[j]+1))
    }
    else{
      slice_vec = append(slice_vec, seq(second_date_breaks[j],first_date_breaks[j], by = 1))
      group_by_vec = append(group_by_vec, rep(j,(first_date_breaks[j] - second_date_breaks[j]+1)))
    }
  }
  
  processed = data %>%
    slice(slice_vec) %>%
    tibble::add_column(group_by_vec = group_by_vec)%>%
    #group by the group_by_vec
    group_by(group_by_vec)%>%
    #summation of data (and convert to mm)
    dplyr::summarise(sum = sum(value/10, na.rm = T)) %>%
    mutate(spi = gamma_fit_spi(sum, export_opts = 'SPI', return_latest = F))
  
  return(processed$spi[length(processed$spi)])
}

# data = data_precip %>%
#   filter(name == 'V1')
# index = (length(data$time)-365):length(data$time)
# timescale = 30

cells = unique(data_precip$name)

cl = makeCluster(30)
registerDoParallel(cl)

out = foreach(i = 1:length(cells), .packages = c('tidyverse', 'lubridate')) %dopar% {
  temp_data = data_precip %>% 
    filter(name == cells[i])
  
  index = which(temp_data$time == as.Date('2010-01-01')): length(temp_data$time)
  #index = which(temp_data$time == as.Date('2021-10-01')): length(temp_data$time)
  
  export = tibble(t_30 = apply(index %>% as.data.frame, MARGIN = 1, FUN = compute_spi, data = temp_data, timescale = 30),
         t_60 = apply(index %>% as.data.frame, MARGIN = 1, FUN = compute_spi, data = temp_data, timescale = 60),
         t_180 = apply(index %>% as.data.frame, MARGIN = 1, FUN = compute_spi, data = temp_data, timescale = 180),
         t_365 = apply(index %>% as.data.frame, MARGIN = 1, FUN = compute_spi, data = temp_data, timescale = 365),
         t_720 = apply(index %>% as.data.frame, MARGIN = 1, FUN = compute_spi, data = temp_data, timescale = 720)) %>%
    mutate(time = temp_data$time[index],
           id = cells[i]) %>%
    dplyr::select(id, time, t_30, t_60, t_180, t_365, t_720)
  export
}

stopCluster(cl)

final = bind_rows(out) %>%
  pivot_longer(cols = -c(id, time))

write_csv(final, file = '/home/zhoylman/gnss-drought-data/russian-spi-pixel-specific.csv')

summary = final %>%
  group_by(time, name) %>%
  summarise(mean_spi = mean(value, na.rm = T)) 

write_csv(summary %>% 
            pivot_wider(names_from = name, values_from = mean_spi) %>%
            dplyr::select(time, t_30, t_60, t_180, t_365, t_720), file = '/home/zhoylman/gnss-drought-data/russian-spi-summary.csv')

ggplot(data = summary %>% filter(name == 't_720' & time > as.Date('2013-01-01') & time < as.Date('2017-01-01')))+
  geom_line(aes(x = time, y = mean_spi, color = name))
 