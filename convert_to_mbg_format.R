##############################################
# Code author: Erin Stearns
# Code objective: Prepare data for MBG format
# Date: 4.3.2017
#############################################

rm(list = ls()) 

# ------------------------------------------- set up environment -------------------------------------------
#load packages
pacman::p_load(fields, tidyverse, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, 
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, sp, rgeos, maptools)
#useful directories
file_dir_gahi <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/gahi_join')
file_dir_sr <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/sys_rev')
save_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/mbg')
needs_review_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/cleaning/needs_review')

# ------------------------------------------ load data ----------------------------------------------------

#make sure when loading data, that have source appropriately assigned

to_convert <- fread(file = paste0(file_dir_gahi, '/rtg_4_3.csv'), stringsAsFactors = F)
to_convert <- to_convert[,-c('V1')]

to_convert[gahi_site == 1 & gahi_restr == 1, data_source := 'gahi_both']
to_convert[is.na(data_source) & gahi_restr == 1, data_source := 'gahi_restr']
to_convert[is.na(data_source) & gahi_site == 1, data_source := 'gahi_site']

# ------------------------------------------ get rows with necessary geographic information --------------
a <- copy(to_convert)

# rows without necessary geo information
no_geo_info <- a[is.na(latitude) & is.na(shapefile),] 
write.csv(no_geo_info, file = paste0(needs_review_dir, '/needs_lat_long_or_shapefile_4.3.csv'))

# rows with necessary geo information
a2 <- a[(!is.na(latitude) & !is.na(longitude)) | (!is.na(shapefile) & !is.na(location_code)),]

# ------------------------------------------ get rows with lf prevalence information --------------------
b <- copy(a2)

#create 2 vectors - 1 for prev vars, 1 for clinical vars
dxprevvars <- c('mf', 'mf_filtration', 'ict', 'elisa', 'br', 'chamber', 'knott', 'ifi', 'pcr', 'new_rapid_test')
clinicvars <- c('hydrocele', 'lymph', 'clinic')

#need to make a fix
b[, prev_knott := prev_knott_]
b <- b[,-c('prev_knott_')]

b[, prev_dx := '']
b[, clinic_dx := '']

#looping over vars to create indicator var detailing what tests a given row contains
for (d in (dxprevvars)) {
  dx_var = paste0('prev_', d)
  newvar = d
  b[(!is.na(get(dx_var)) & !grepl(newvar, prev_dx)), prev_dx := paste0(prev_dx, " ", newvar)]
}


for (c in (clinicvars)) {
  dx_var = paste0('pop_', c)
  newvar = c
  b[(!is.na(get(dx_var)) & !grepl(newvar, clinic_dx)), clinic_dx := paste0(clinic_dx, " ", newvar)]
}

#subset to rows with prevalence measured
b1 <- b[prev_dx != '',]

#drop all clinical vars
#vector of all clinical vars
drop_clinvars <- names(b1)[grepl(('hydrocele|lymph|clinic'), names(b1))]

b2 <- b1[,-drop_clinvars, with = F]

# ------------------------------------------ get single dx field --------------------------------
c <- copy(b2)

#dropping rows with knott, ifi, new_rapid_test, pcr

c1 <- c[prev_dx != ' knott',]
c2 <- c1[prev_dx != ' ifi',]
c3 <- c2[prev_dx != ' new_rapid_test',]
c4 <- c3[prev_dx != ' pcr',]
c5 <- c4[prev_dx != ' knott ifi',]

# create MBG fields
d <- copy(c5)

# for reference
prev_metrics <- table(d$prev_dx)
write.csv(prev_metrics, file = paste0(file_dir_gahi, '/dx_db.csv'))

d[,N:=0]
d[,lf_pos := 0]
d[,lf_prev := 0]

# convert all dx vars into integers
dxvar4convert <- names(d)[grepl(('mf|mf_filtration|ict|elisa|br|chamber|knott|ifi|pcr|new_rapid_test'), names(d))]

for (i in dxvar4convert) {
  d[, (i) := as.numeric(get(i))]
}

#populating for elisa
d[prev_dx == ' elisa', N := pop_elisa]
d[prev_dx == ' elisa', lf_pos := np_elisa]
d[prev_dx == ' elisa', lf_prev := prev_elisa]

#populating for mf_filtration
d1 <- copy(d)

d1[prev_dx == ' mf_filtration', N := pop_mf_filtration]
d1[prev_dx == ' mf_filtration', lf_pos := np_mf_filtration]
d1[prev_dx == ' mf_filtration', lf_prev := prev_mf_filtration]

#populating for ict - making this the default in the presence of other prev measures
d2 <- copy(d1)

d2[grepl('ict', d2$prev_dx), N := pop_ict]
d2[grepl('ict', d2$prev_dx), lf_pos := np_ict]
d2[grepl('ict', d2$prev_dx), lf_prev := prev_ict]

#populating for mf
d3 <- copy(d2)

d3[!grepl('ict', d2$prev_dx) & grepl('mf', d2$prev_dx) & prev_dx != ' mf_filtration', N := pop_mf]
d3[!grepl('ict', d2$prev_dx) & grepl('mf', d2$prev_dx) & prev_dx != ' mf_filtration', lf_pos := np_mf]
d3[!grepl('ict', d2$prev_dx) & grepl('mf', d2$prev_dx) & prev_dx != ' mf_filtration', lf_prev := prev_mf]

#168 rows without numerator or denominator information, need to be checked
d3_no <- d3[is.na(N),]
write.csv(d3_no, file = paste0(needs_review_dir, '/no_denom_num_4.3.17.csv'))

#only take rows with denom and num
d4 <- d3[!is.na(N) & !is.na(lf_pos),]

# --------------------------------------------------- bringing in sampling and cluster var ------------------------------------
e <- copy(d4)

sampling <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/cleaning/geo_cleaning/LF_polygon_classification_GAHI_Feb_22_EAC.csv', stringsAsFactors = F)
sampling <- sampling[!is.na(Sampling),]
sampling <- sampling[,c('gahi_uid', 'Sampling', 'n_clusters'), with = F]

e1 <- left_join(e, sampling, by = c('uid_from_source'= 'gahi_uid'))

# --------------------------------------------------- subsetting to only necessary fields -------------------------------------
f <- copy(e1)

mbg_fields <- c("gahi_join_uid", "year_start", "GAUL_ADM0_NAME", "data_source", "point", "latitude", "longitude", "shapefile", "location_code", 
                "N", "lf_pos", "lf_prev")

f1 <- as.data.table(f[, mbg_fields])

#rename
f1[, lf_uid := gahi_join_uid]
f1[, start_year := year_start]
f1[, country := GAUL_ADM0_NAME]
f1[, source := data_source]

f2 <- f1[, -c('gahi_join_uid', 'year_start', 'GAUL_ADM0_NAME', 'data_source'), with = F]

#reorder
neworder <- c('lf_uid', 'start_year', 'country', 'source', 'point', 'latitude', 'longitude', 'shapefile', 'location_code', 'N', 'lf_pos', 
              'lf_prev')

setcolorder(f2, neworder)

#save

write.csv(f2, file = paste0(save_dir, '/gahi_mbg_4.3.17.csv'))

gmbg <- fread(file = paste0(save_dir, '/gahi_mbg_4.3.17.csv'), stringsAsFactors = F)
# -------------------------------------------------- bringing in gahi and sr data from Liz ----------------------------------

fromliz <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/gahi_sr_fromLiz.csv', stringsAsFactors = F)

fromliz <- fromliz[,-c('V1', 'poly_reference', 'ADM0_name', 'poly_id_field', 'poly_id_field_name.x', 'poly_id_field_name.y', 'n', 'max_n', 'new_file', 'shapenames', 'poly_ref_cl')]

g <- copy(fromliz)

g[,lf_uid:=UID_plot]
g[,country:=iso3]
g[,start_year := year_start]

#join to sampling
eg1 <- as.data.table(left_join(g, sampling, by = c('uid_from_source' = 'gahi_uid')))

eg1[,sampling:=Sampling]

eg2 <- eg1[, -c('uid_from_source', 'year_start', 'iso3', 'UID_plot', 'Sampling')]
neword1 <- c('lf_uid', 'start_year', 'country', 'source', 'point', 'latitude', 'longitude', 'shapefile', 'location_code', 'sampling', 'n_clusters', 'N', 'lf_pos', 
             'lf_prev')

setcolorder(eg2, neword1)

#drop rows with null lf_pos & lf_prev

eg3 <- eg2[!is.na(lf_pos),]

pt <- eg3[point == 1,]

write.csv(eg3, file = paste0(save_dir, '/liz_gahi_extr_mbg_4.4.17.csv'))



#missing shapefile paths for important polygons, need to join to my datasets and see if that fills things out
shps <- as.data.table(f[,c('uid_from_source', 'ADM0_name', 'shapefile', 'location_code')])
shps[,shp_file := shapefile]
shps[,loc_code := location_code]
shps[,uid:=uid_from_source]

shps1 <- shps[,-c('shapefile', 'location_code', 'uid_from_source')]


#make a join to liz's data

wshps <- as.data.table(left_join(g, shps1, by = c('uid_from_source' = 'uid')))

# added little


# need to run nearestLand on points
# ----------------------------------- push to land ---------------------------------------

t <- copy(post_geomatch)

#making sure point var is up to date
table(t$point)
dblcheck <- t[!is.na(latitude) & !is.na(longitude),]

# data ---------------
t <- fread(file ='H:/ntds/LF/Dataliz_gahi_extr_mbg_4.4.17.csv')
# points
lf_pts <- t[point ==1 ,]
lf_pts1 <- lf_pts[,c('latitude','longitude')]
names(lf_pts1) <- c('y', 'x') 

#lf polygons for joining later
lf_poly <- t[point == 0,]

#raster grid
grid <- raster("J:/WORK/11_geospatial/01_covariates/02_Oxford/01_Global_Masks/Land_Sea_Masks/CoastGlobal_5k.tif")

#running nearest land -----------------
lf_corrected <- nearestLand(lf_pts1, grid, 100) 

# count how many were moved - 3144 points moved!
sum(!is.na(lf_corrected[, 1]))

#now joining to main data set

lf_rightpts <- cbind(lf_pts, lf_corrected)

# now creating an'adjusted' field
lf_rightpts$adjusted <- ifelse(!is.na(lf_rightpts$x) & !is.na(lf_rightpts$y), 1, 0)
lf_rightpts <- as.data.table(lf_rightpts)

# now creating a single lat and single long field

lf_updated <- lf_rightpts[!is.na(lf_rightpts$x), lat:=x]
lf_updated <- lf_rightpts[!is.na(lf_rightpts$y), long:=y]
lf_updated <- lf_rightpts[is.na(lf_rightpts$lat), lat:=latitude]
lf_updated <- lf_rightpts[is.na(lf_rightpts$long), long:=longitude]

lf_updated <- lf_updated[,-c("x", "y")]

#join to polygon data
lf_poly <- as.data.frame(lf_poly)
lf_updated <- as.data.frame(lf_updated)

lf_updated1 <- as.data.table(rbind.fill(lf_updated, lf_poly))


#drop original lat and long fields 
lf_updated2 <- lf_updated1[,-c("latitude", "longitude")]

#rename new lat and long to latitude and longitude

lf_updated3 <- lf_updated2[!is.na(lf_updated2$lat), latitude := lat]
lf_updated3 <- lf_updated2[is.na(lf_updated2$lat), latitude := NA]

lf_updated3 <- lf_updated2[!is.na(lf_updated2$long), longitude := long]
lf_updated3 <- lf_updated2[is.na(lf_updated2$long), longitude := NA]

# drop old lat and long
lf_updated4 <- lf_updated3[,-c("lat", "long")]

#save

write.csv(lf_updated4, file = paste0(save_dir, '/some_geo_cl_left_but_most_done.4.3.csv'))#ran on virtual machine, newest on h:drive

write.csv(lf_updated4, file = 'H:/ntds/LF/Data/sr_gahi_ptscorrect.csv')

up_to_date <- fread(file = paste0(save_dir, '/sr_gahi_ptscorrect.csv'), stringsAsFactors = F)
up_to_date <- up_to_date[,-c('V1')]






