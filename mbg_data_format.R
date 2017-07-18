##############################################
# Code author: Erin Stearns
# Code objective: Prepare data for MBG 
# Date: 5.22.17
#############################################

rm(list = ls()) 

# ------------------------------------------- set up environment -------------------------------------------
my_lib <- ("C:/Users/stearns7/Documents/R/win-library/3.3")
.libPaths(my_lib)
pacman::p_load(fields, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, ggmap, devtools, scale,
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, seegMBG, sp, rgeos, maptools, rgeos, rgdal)

latest_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/all/')
save_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/mbg/formatted/')
data_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/latest/')

#------------------------------------------ load data ----------------------------------------------------

#all <- fread(file = paste0(latest_dir, '/lf_data_5.22.csv'), stringsAsFactors = F) #9096 rows
#all <- all[,-c('V1')]

#data from liz, cleaned up & points pushed to land 4.4.17 - 7246
#liz <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/mbg_data/currently_in_use/lf_pos_4.18.csv', stringsAsFactors = F)
#liz <- liz[,-c('V1')]

#latest (6.1.17) -- liz recovered some rows missing year info
#all <- fread(file = paste0(latest_dir, '/lf_data_fin.csv'), stringsAsFactors = F) #9096 rows
#all <- all[,-c('V1')]

# latest clean lf dataset (as of 6/28/2017, last edited on 6/28/2017 at 14:17, assuming by Kevin)
#all <- fread(paste0(data_dir, 'lf_data_fin.csv'), stringsAsFactors = F)

# latest clean lf dataset (as of 7/18/2017, last edited on 7/12/2017 at 11:13(not by me)) -- 9128 rows
all <- fread(paste0(data_dir, 'lf_data_fin.csv'), stringsAsFactors = F)
#all <- all[,-c('V1')]


# ----------------------------------------- transform full dataset to MBG format -------------------------
a <- copy(all)

no_prev <- a[prev_dx=='',]
#write.csv(no_prev, file = paste0(data_dir, 'no_prev_7.18.csv'))

a1 <- a[prev_dx != '',] #9015 with latest (6/28)
a1 <- a1[,-c(68:76), with=FALSE]

# ---------------------------------------- getting single set of prev measures --------------------------
b <- copy(a1)

#creating mbg prev fields
b[,N:=-99]
b[,lf_pos:=-99]
b[,lf_prev:=-99]

#get vector of dx fields to loop over
popvars <- names(b)[grepl(("pop"), names(b))]
npvars <- names(b)[grepl(("np"), names(b))]
prevvars <- names(b)[grepl(("prev"), names(b))]
remove <- c('prev_dx','no_prev','lf_prev')
prevvars <- prevvars[!prevvars %in% remove]

loopvars <- c(popvars,npvars, prevvars)

b1 <- copy(b)

#creating indicator to tell loop what dx to use to fill mbg prev fields
b1[grepl('mf', b1$prev_dx) & !is.na(pop_mf) & pop_mf > 0, dx_to_use:='mf']
b1[grepl('ict', b1$prev_dx) & !is.na(pop_ict) & pop_ict > 0 & is.na(dx_to_use), dx_to_use:='ict']
b1[grepl('elisa', b1$prev_dx) & !is.na(pop_elisa) & pop_elisa > 0 & is.na(dx_to_use), dx_to_use:='elisa']
b1[grepl('br', b1$prev_dx) & !is.na(pop_br) & pop_br > 0 & is.na(dx_to_use), dx_to_use:='br']
b1[grepl('pcr', b1$prev_dx) & !is.na(pop_pcr) & pop_pcr > 0 & is.na(dx_to_use), dx_to_use:='pcr']

#rows missing denominators or numerators or report a dx as 'other' or 'density -- 404
set_aside <- b1[is.na(dx_to_use),]
#write.csv(set_aside, file=paste0(data_dir,'miss_denom_oth_dens_rows.csv'))

#rows with all necessary prev info
b2 <- b1[!is.na(dx_to_use),] #8611

b3 = sapply(sort(unique(b2$dx_to_use)), function(i) {
  b2[dx_to_use==i, N := get(paste0('pop_', i))]
  b2[dx_to_use==i, lf_pos := get(paste0('np_', i))]
  b2[dx_to_use==i, lf_prev := get(paste0('prev_', i))]
})

#save summary of what diagnostics used in this dataset
summ_dx <- table(b2$dx_to_use)
#write.csv(summ_dx, file=paste0(data_dir, 'summary/summ_dx_used_7.18.csv'))

b2[,dx:=dx_to_use]

# ---------------------------------------- adding necessary fields ---------------------------------------
c <- copy(b2)

c[,weight:=1]
c[,cluster_id:=1]
c1[,lf_uid:=rownames(c1)]
c1[,source:=data_source]

c1 <- c[,original_year:=year_start]

c1[,year:=year_start]
c1[year_start <=2000, year:=2000]

#drop rows with N=0, & year start =0
c2 <- c1[,no_year:=0]
c2[year_start==0, no_year:=1]
c2[is.na(year_start), no_year:=1]
c2[,no_denom:=0]
c2[N==0, no_denom:=1]

c3 <- c2[no_year ==0 & no_denom ==0,]

#place those in need of attention aside
c4 <- c2[no_year == 1 | no_denom == 1,]
#write.csv(c4, file=paste0(data_dir, '/data_needs_review/no_year_7.18.csv'))

# -------------------------------------- assigning iso3 codes ------------------------------------------
d <- copy(c3)

d[,country:=adm0]

#need iso3 codes for R shiny and data plots
iso3 <- fread(file='C:/Users/stearns7/OneDrive - UW Office 365/spatial_data/mbg_gaul_iso3_key_lf_adapted.csv', stringsAsFactors = F)

#removing accent from cote d'ivoire
#iso3 <- iso3[official_name=="Cote D'Ivoire", official_name:="Cote d'Ivoire"]
#iso3 <- iso3[official_name == 'United Republic of Tanzania', official_name:="Tanzania"]
#iso3 <- iso3[official_name == 'Venezuela (Bolivarian Republic of)', official_name:="Venezuela"]
#iso3 <- iso3[official_name == 'Viet Nam', official_name:="Vietnam"]
#iso3 <- iso3[official_name == 'Wallis and Futuna Islands', official_name:="Wallis and Futuna"]

#write.csv(iso3, file = 'C:/Users/stearns7/OneDrive - UW Office 365/spatial_data/mbg_gaul_iso3_key_lf_adapted.csv')

#match to data
d1 <- merge(d, iso3, by.x = 'adm0', by.y='official_name')

#make iso3 code country var
d2 <- d1[, country:=iso3]

# ------------------------------------ deal with sampling (temporary) ----------------------------------
d3 <- copy(d2)
d3[point==0,sampling:=99]
d3[point==0 & sampling_prod==0, sampling:=0]
d3[point==0 & sampling_prod==1, sampling:=1]
d3[point==0 & sampling==99, sampling:=1]


#save latest
write.csv(d3, file=paste0(data_dir, 'lf_data_pre_mbg_trim&poly_rev_drop_7.18.csv'))

#subset to polygon rows that have been reviewed by Liz
d4 <- d3[poly_rev != 1 | is.na(poly_rev),]

# ------------------------------------ drop rows with N < 20  ----------------------------------

d5 <- copy(d4)

#rows with N < 20
lt20 <- d5[N < 20,]
#write.csv(lt20, file=paste0(data_dir, 'data_needs_review/n_less_than_20_7.18.csv'))

#dropping the rows with N < 20
d6 <- d5[N >= 20,]

#Add lf_uid
row.names(d6) <- 1:nrow(d6)
d6$lf_uid <- 1:nrow(d6)

#setting source
d6[,source:=data_source]

# -------------------------------------- trimming fields -------------------------------------------------
e <- copy(d6)

e1 <- e[,c('lf_uid', 'source', 'original_year', 'country', 
           'point', 'latitude', 'longitude', 'shapefile', 'location_code', 'sampling', 'N', 
           'lf_pos', 'lf_prev', 'weight', 'cluster_id', 'year'), with=F]

write.csv(e1,file = paste0(save_dir,'lf_pos_all.csv'))

#points only
e2 <- e1[point==1,]
e2 <- e2[,-c('point', 'shapefile', 'location_code', 'sampling'), with=F]

write.csv(e2,file = paste0(save_dir,'lf_pos_pts_only.csv'))



