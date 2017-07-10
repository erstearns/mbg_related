##############################################
# Code author: Erin Stearns
# Code objective: Get polygon data in format for MBG - polygon resampling
# Date: 5.8.2017
#############################################


rm(list = ls()) 

# ----------------------------------------------- set up environment -------------------------------------------
my_lib <- ("C:/Users/stearns7/Documents/R/win-library/3.3")
.libPaths(my_lib)
pacman::p_load(fields, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, 
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, seegMBG, sp, rgeos, maptools, rgeos, rgdal)

save_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/mbg/polygon_resampling')

# ------------------------------------------------ load data ---------------------------------------------


clean <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/gahi_join/gahi_clean_fin_5.5.csv',stringsAsFactors = F)
polys <- clean[point==0,]
table(polys$adm0)

ghana <- polys[adm0 == 'Ghana',]

#pulling in sampling code
sampling <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/cleaning/geo_cleaning/LF_polygon_classification_GAHI_Feb_22_EAC.csv', stringsAsFactors = F)
sampl_1 <- sampling[,c('gahi_uid', 'Sampling'), with=F]

ghana_smpl <- merge(ghana, sampl_1, by.x = 'uid_from_source_site', by.y = 'gahi_uid')


ghana_mbg <- ghana_smpl[,c('gahi_join_uid', 'year_start', 'adm0', 'data_source', 'point', 'latitude', 'longitude', 'shapefile', 'location_code', 
                      'Sampling', 'pop_mf', 'np_mf', 'prev_mf'), with=F]

ghana_mbg[,year:=year_start]
ghana_mbg[year_start < 2000, year:= 2000]

ghana_mbg[,cluster_id:=seq(1:nrow(ghana_mbg))]

newnames <- c('lf_uid', 'start_year', 'country', 'source', 'point', 'latitude', 'longitude',
              'shapefile', 'location_code', 'sampling', 
              'N', 'lf_pos', 'lf_prev', 'year', 'cluster_id')

setnames(ghana_mbg, old = names(ghana_mbg), new = newnames)


write.csv(ghana_mbg, file = paste0(save_dir, '/ghana_polys.csv'))
