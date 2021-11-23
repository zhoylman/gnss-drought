library(tidyverse)
#source all relevant drought functions
source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/master/processing/ancillary-functions/R/drought-functions.R')

gnss_data = list.files('/home/zhoylman/gnss-drought/data/gnss', full.names = T) %>%
  lapply(., read_csv) %>%
  lapply(., mutate, date = X1) %>%
  lapply(., select, -X1) %>%
  lapply(., select, c(date, vertical))

name_short = list.files('/home/zhoylman/gnss-drought/data/gnss')

for(i in 1:3){
  gnss_data[[i]] = gnss_data[[i]] %>%
    mutate(id = name_short[i])
}


final = gnss_data %>%
  bind_rows() %>%
  mutate(yday = lubridate::yday(date)) %>%
  select(id, date, yday, vertical) %>%
  group_by(yday, id) %>%
  mutate(generalizedLogistic_vertical = glo_fit_spei(vertical, export_opts = 'SPEI', return_latest = F),
         zScore_vertical = (vertical - mean(vertical, na.rm = T))/sd(vertical, na.rm = T))

#sanity check
test = final %>%
  filter(id == 'P190_UNR_offrem1.txt' & yday == 152)

ggplot()+
  geom_point(data = final, aes(x = vertical, generalizedLogistic_vertical, color = id), shape = 21)+
  facet_wrap(~id, ncol = 1)+
  labs(x = 'Raw Vertical', y = 'Drought Standardized Vertical (generalized logistic)')+
  theme_bw()

write_csv(final, '/home/zhoylman/gnss-drought/data/output/drought_standardized_gnss.csv')
