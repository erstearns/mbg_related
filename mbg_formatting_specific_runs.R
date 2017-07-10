##############################################
# Code author: Erin Stearns
# Code objective: Prepare data for MBG 
# Date: 4.3.2017
#############################################

rm(list = ls()) 

# ------------------------------------------- set up environment -------------------------------------------
pacman::p_load(fields, tidyverse, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, 
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, sp, rgeos, maptools)

file_dir_gahi <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/gahi_join')
file_dir_sr <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/sys_rev')
save_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/mbg')
needs_review_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/cleaning/needs_review')

#------------------------------------------ load data ----------------------------------------------------

gahi_mbg <- fread(file = paste0(save_dir, '/gahi_mbg_4.3.17.csv'), stringsAsFactors = F)
gahi_mbg <- gahi_mbg[,-c('V1')]


#data from liz, cleaned up & points pushed to land 4.4.17
ge <- fread(file = paste0(save_dir, '/sr_gahi_ptscorrect.csv'), stringsAsFactors = F)
ge <- ge[,-c('V1')]
# ---------------------------------------- pts only run, 2000-2015 -----------
ge_pts <- ge[point == 1,]

#changing cols
mbg_input <- ge[, c(1:4, 14:15, 10:12)]

#adjusting years

mbg_input1 <- mbg_input[,year:=start_year]
mbg_input1 <- mbg_input1[year <= 2000, year := 2000]

write.csv(mbg_input1, file = paste0(save_dir, '/ok_4.4almost4.5.csv'))

#----------------------------------------- points only run -----------------------------------------------

#specifically doing this for Ethiopia run
gahi_pts <- gahi_mbg[point == 1,]

#changing cols
mbg_input <- gahi_pts[, c(1:7, 10:12)]

mbg_input[, year := start_year]
mbg_input1 <- mbg_input[country == 'Ethiopia' & year <= 2000, year := 2000]

iso <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/spatial_data/iso_key.csv', stringsAsFactors = F)

iso <- iso[,-c(2,4)]
#latest
almostdata <- fread(file = paste0(save_dir, '/lf_pos_4ethrun.4.4.17.csv'), stringsAsFactors = F)

iso_data <- as.data.table(left_join(almostdata, iso, by ='country'))
iso_data <- as.data.table(left_join(mbg_input, iso, by ='country'))


iso_data2 <- iso_data[country == 'Democratic Republic of the Congo' | country == 'Congo', iso3 := 'COD']
iso_data3 <- iso_data2[country == 'Micronesia (Federated States of)', iso3 := 'FSM']
iso_data4 <- iso_data3[country == 'Republic of Korea', iso3 := 'KOR']
iso_data5 <- iso_data4[country == 'United Republic of Tanzania', iso3 := 'TZA']
iso_data6 <- iso_data5[country == 'Venezuela', iso3 := 'VEN']
iso_data7 <- iso_data6[country == 'Wallis and Futuna', iso3 := 'WLF']
iso_data8 <- iso_data7[country == "CoteDivoire", iso3 := 'CIV']

iso_data9 <- iso_data8[!is.na(iso3),]
iso_data10 <- iso_data9[,-c('country')]
iso_data10[,country:=iso3]
iso_data11 <- iso_data10[,-c('iso3')]
iso_data12 <- iso_data11[,c(1:2,10,3, 5:9)]
iso_data12[,year:=start_year]
iso_data12[start_year <= 2000, year:= 2000]
iso_data12[,weight:=1]

write.csv(iso_data12,file = paste0(save_dir, '/gahi_only_4.11.csv'))

write.csv(mbg_input, file = paste0(save_dir, '/lf_pos_4ethrun.4.4.17.csv'))
check <- fread(file = paste0(save_dir, '/lf_pos_4ethrun.4.4.17.csv'), stringsAsFactors = F)

# --------------------------------------- pulling in background poly data for sri lanka ---------------------

bgdata <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/spatial_data/bg_lka.csv', stringsAsFactors = F)

ge_join <- ge[,c(1:5, 14:15, 6:8, 10:12)]
  
data <- rbind(ge_join, bgdata)
  
#looking to see where i need to add fake data

#sri lanka - points in all necessary time periods
lka <- data[country == 'LKA' & point == 1,]
table(lka$start_year)

#ethiopia - no, need 2000 and 2005

#ghana - need 2015 data points
gha <- data[country == 'GHA' & point == 1,]
table(gha$start_year)

#togo - need 2015 data points
tgo <- data[country == 'TGO' & point == 1,]
table(tgo$start_year)

#benin - need 2005, 2010, 2015 data points
ben <- data[country == 'BEN' & point == 1,]
table(ben$start_year)

#creating binned year var
data1 <- data[, year := start_year]
data2 <- data1[year <= 2000, year := 2000]

write.csv(data2, file = paste0(save_dir, 'needsfakedata.csv'))


#load in fake data 

fake <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/cleaning/fake_data_to_get_to_run.csv', stringsAsFactors = F)

#bind to main

data_wfake <- rbind(data2, fake)

write.csv(data_wfake, file = paste0(save_dir, '/all_data_4.5.17.csv'))


#polygons for sri lanka
lka <- data_wfake[country=='LKA',]
lka_poly <- lka[point==0 & shapefile != 'NULL',]
write.csv(lka_poly, file = paste0(save_dir, '/lka_poly.csv'))



benin <- data_wfake[country=='BEN',]
benin_pts <- benin[,c('latitude', 'longitude')]
points(benin_pts)


#plotting
world <- readOGR('J:/WORK/11_geospatial/06_original shapefiles/GAUL_admin/admin0/g2015_2014_0', 'g2015_2014_0')
world <- fortify(world, region='GEO_ID')
lf_data_pts <- data_wfake[point==1, c('latitude', 'longitude')]

plot(world)

# ------------------------------------------------ Loading Liz's data as-is --------------------------------------------
gah_ext <- fread(file = ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/gahi_sr_fromLiz.csv'), stringsAsFactors = F)
gah_ext <- gah_ext[,-c('V1')]

#subset to points only
a <- gah_ext[point ==1,]

#creating MBG fields
b <- a[,lf_uid := UID_plot]
c <- b[,start_year := year_start]  
d <- c[,year:=year_start]
d1 <- d[year_start <= 2000, year:=2000]
e <- d1[,weight:=1]
f <- e[,country:=iso3]

#mbg format
g <- f[,c(24:25, 28, 22, 5:6, 11:13, 27, 26)]

#just want to see if any fields are null that shouldn't be - 5 rows with null lf_pos and prev
mbgcols <- names(g)
nas <- g[, lapply(.SD, function(x) length(x[is.na(x)])), .SDcols=mbgcols]

#going to drop the 5 rows with null lf_pos and prev
g1 <- g[!is.na(lf_pos),]

#save
write.csv(g1, file = paste0(save_dir, '/from_liz.csv'))
