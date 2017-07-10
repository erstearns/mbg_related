##############################################
# Code author: Erin Stearns
# Code objective: Prepare data for MBG 
# Date: 5.22.17
#############################################

rm(list = ls()) 

# ------------------------------------------- set up environment -------------------------------------------
my_lib <- ("C:/Users/stearns7/Documents/R/win-library/3.3")
.libPaths(my_lib)
pacman::p_load(fields, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, 
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, seegMBG, sp, rgeos, maptools, rgeos, rgdal)

latest_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/all')
save_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/mbg/formatted')
data_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/latest/')

#------------------------------------------ load data ----------------------------------------------------

all <- fread(file = paste0(latest_dir, '/lf_data_5.22.csv'), stringsAsFactors = F) #9096 rows
#all <- all[,-c('V1')]

#data from liz, cleaned up & points pushed to land 4.4.17 - 7246
liz <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/mbg_data/currently_in_use/lf_pos_4.18.csv', stringsAsFactors = F)
#liz <- liz[,-c('V1')]

#latest (6.1.17) -- liz recovered some rows missing year info
all <- fread(file = paste0(latest_dir, '/lf_data_fin.csv'), stringsAsFactors = F) #9096 rows
#all <- all[,-c('V1')]

# latest clean lf dataset (as of 6/28/2017, last edited on 6/28/2017 at 14:17, assuming by Kevin)
all <- fread(paste0(data_dir, 'lf_data_fin.csv'), stringsAsFactors = F)

# ----------------------------------------- transform full dataset to MBG format -------------------------
a <- copy(all)

no_prev <- a[prev_dx=='',]
write.csv(no_prev, file = paste0(data_dir, 'no_prev.csv'))

a1 <- a[prev_dx != '',] #9015 with latest (6/28)
a1 <- a1[,-c(69:77), with=FALSE]

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
#write.csv(summ_dx, file=paste0(data_dir, 'summ_dx_used.csv'))

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
#write.csv(c4, file=paste0(data_dir, 'no_year.csv'))

# -------------------------------------- assigning iso3 codes ------------------------------------------
d <- copy(c3)

d[,country:=adm0]

#need iso3 codes for R shiny and data plots
iso3 <- fread(file='C:/Users/stearns7/OneDrive - UW Office 365/spatial_data/mbg_gaul_iso3_key.csv', stringsAsFactors = F)

#removing accent from cote d'ivoire
iso3 <- iso3[official_name=="Cote D'Ivoire", official_name:="Cote d'Ivoire"]
iso3 <- iso3[official_name == 'United Republic of Tanzania', official_name:="Tanzania"]
iso3 <- iso3[official_name == 'Venezuela (Bolivarian Republic of)', official_name:="Venezuela"]
iso3 <- iso3[official_name == 'Viet Nam', official_name:="Vietnam"]
iso3 <- iso3[official_name == 'Wallis and Futuna Islands', official_name:="Wallis and Futuna"]

write.csv(iso3, file = 'C:/Users/stearns7/OneDrive - UW Office 365/spatial_data/mbg_gaul_iso3_key_lf_adapted.csv')

#match to data
d1 <- merge(d, iso3, by.x = 'adm0', by.y='official_name')

# ------------------------------------ deal with sampling (temporary) ----------------------------------
d2 <- copy(d1)
d2[point==0,sampling:=99]
d2[point==0 & sampling_prod==0, sampling:=0]
d2[point==0 & sampling_prod==1, sampling:=1]
d2[point==0 & sampling==99, sampling:=1]

# -------------------------------------- trimming fields -------------------------------------------------
e <- copy(d2)

e1 <- e[,c('lf_uid', 'source', 'original_year', 'country', 
           'point', 'latitude', 'longitude', 'shapefile', 'location_code', 'sampling', 'N', 
           'lf_pos', 'lf_prev', 'weight', 'cluster_id', 'year'), with=F]

write.csv(e1,file = paste0(data_dir,'lf_pos_all.csv'))



#scrap



# ---------------------------------------- subset to prev only rows -- no sequela --------------------------
a <- copy(all)

no_prev <- a[prev_dx=='',]
#write.csv(no_prev, file = paste0(latest_dir, '/no_prev.csv'))

a1 <- a[prev_dx != '',] #8914 with latest (6/1)

# ---------------------------------------- getting single set of prev measures --------------------------
b <- copy(a1)

b[,N:=NA]
b[,lf_pos:=NA]
b[,lf_prev:=NA]

#mf rows
b_mf <- b[grep('mf', b$prev_dx),] #3198 rows
b_mf[,N:=pop_mf]
b_mf[,lf_pos:=np_mf]
b_mf[,lf_prev:=prev_mf]

#set aside rows that report prev only -- 22 rows
b_mf_prev_only <- b_mf[is.na(N),]
#write.csv(b_mf_prev_only, file = paste0(latest_dir, '/prev_only_mf.csv'))

b_mf_fin <- b_mf[!is.na(N),] #3176 rows

#ict rows
mf_uids <- b_mf[['uid_from_source']]
b_next <- b[!uid_from_source %in% mf_uids,]

b_ict <- b_next[grep('ict', b_next$prev_dx),] #5421 rows
b_ict[,N:=pop_ict]
b_ict[,lf_pos:=np_ict]
b_ict[,lf_prev := prev_ict]

b_ict_prev_only <- b_ict[is.na(N),] #35
#write.csv(b_ict_prev_only, file = paste0(latest_dir, '/prev_only_ict.csv'))

b_ict_fin <- b_ict[!is.na(N),] #5386

#elisa rows
ict_uids <- b_ict[['uid_from_source']]

b_next2 <- b_next[!uid_from_source %in% ict_uids,]

b_elisa <- b_next2[grep('elisa', b_next2$prev_dx),] #166 rows
b_elisa[,N:=pop_elisa]
b_elisa[,lf_pos:=np_elisa]
b_elisa[,lf_prev:=prev_elisa]

#br rows
elisa_uids <- b_elisa[['uid_from_source']]

b_next3 <- b_next2[!uid_from_source %in% elisa_uids,]

b_br <- b_next3[grep('br', b_next3$prev_dx),] #21 rows
b_br[,N:=pop_br]
b_br[,lf_pos:=np_br]
b_br[,lf_prev:=prev_br]

#pcr rows
br_uids <- b_br[['uid_from_source']]

b_next4 <- b_next3[!uid_from_source %in% br_uids,]

b_pcr <- b_next4[grep('pcr', b_next4$prev_dx),] #3 rows
b_pcr[,N:=pop_pcr]
b_pcr[,lf_pos:=np_pcr]
b_pcr[,lf_prev:=prev_pcr]

#other rows
pcr_uids <- b_pcr[['uid_from_source']]

b_next5 <- b_next4[!uid_from_source %in% pcr_uids,]

#all that remain are either 'other' or report density only
write.csv(b_next5, file = paste0(latest_dir, '/other_dens_only.csv'))

#bind all back together
mbg_list <- list(b_mf_fin, b_ict_fin, b_elisa, b_br, b_pcr)
b2 <- rbindlist(mbg_list, use.names = T) #8752 rows
#write.csv(b2, file = paste0(latest_dir, '/mbg_prevs_populated_6.2.csv'))

# ---------------------------------------- adding necessary fields ---------------------------------------
c <- copy(b2)

c[,weight:=1]
c[,cluster_id:=1]

c1 <- c[,original_year:=NA]
c1[,original_year:=year_start]

c1[,year:=year_start]
c1[year_start <=2000, year:=2000]

c1[,lf_uid:=rownames(c1)]

c1[,source:=data_source]

c1[,country:=adm0]
#write.csv(c1, file = paste0(save_dir, '/all_mbg_nec_fields_plus.6.2.csv'))

#c1 <- fread(file = paste0(save_dir, '/all_mbg_nec_fields_plus.5.22.csv'), stringsAsFactors = F)
#c1 <- c1[,-c('V1')]

#drop rows with N=0, & year start =0
c2 <- c1[,no_year:=0]
c2[year_start==0, no_year:=1]
c2[is.na(year_start), no_year:=1]
c2[,no_denom:=0]
c2[N==0, no_denom:=1]

c3 <- c2[no_year ==0 & no_denom ==0,]
# c3[adm0_code==90, adm0 := 'Gambia']
# c3[adm0_code == 257, adm0 := 'Tanzania']
# c3[adm0_code == 246, adm0 := 'Trinidad and Tobago']

#place those in need of attention aside
c4 <- c2[no_year == 1 | no_denom == 1,]
#write.csv(c4, file=paste0(latest_dir, '/no_denom_no_year.csv'))

#need iso3 codes for R shiny and data plots
iso3 <- fread(file='C:/Users/stearns7/OneDrive - UW Office 365/spatial_data/iso3.csv', stringsAsFactors = F)

iso3 <- iso3[,c('name', 'ISO3166-1-Alpha-3', 'GAUL'), with=F]
iso3[GAUL==202, name:='Republic of Korea']
iso3[`ISO3166-1-Alpha-3` == 'CIV', name:="Cote d'Ivoire"]
iso3[`ISO3166-1-Alpha-3` == 'COG', name:="Congo"]
iso3[`ISO3166-1-Alpha-3` == 'COD', name:="Democratic Republic of the Congo"]
iso3[`ISO3166-1-Alpha-3` == 'FSM', name:="Micronesia (Federated States of)"]
iso3[`ISO3166-1-Alpha-3` == 'TTO', name:="Trinidad and Tobago"]
iso3[`ISO3166-1-Alpha-3` == 'TZA', name:="Tanzania"]

#match to data
try1 <- merge(c3, iso3, by.x = 'adm0', by.y='name')

try2 <- try1[,-c('country')]
try2[,country:=`ISO3166-1-Alpha-3`]

# -------------------------------------- trimming fields -------------------------------------------------
d <- copy(try2)

d1 <- d[,c('lf_uid', 'source', 'original_year', 'country', 'point', 'latitude', 'longitude', 'shapefile', 'location_code', 'N', 
           'lf_pos', 'lf_prev', 'weight', 'cluster_id', 'year'), with=F]

write.csv(d1,file = paste0(save_dir,'/lf_data_pts&poly_6.1.csv'))


# ---------------------------------------- pts only run, 2000-2015 -----------
e <- copy(d1)

all_pts <- e[point==1,]

#trim
all_pts1 <- all_pts[,-c('point', 'shapefile', 'location_code'), with=F]

write.csv(all_pts1, file = paste0(save_dir, '/pts_only_6.1.csv'))

# --------------------------------------- prev loess curve --------------
pts <- fread(file = paste0(save_dir, '/pts_only_6.1.csv'), stringsAsFactors = F)
pts_ssa <- fread(file = 'C:/Users/stearns7/Desktop/input_data.csv', stringsAsFactors = F)
pts_poly <- fread(file = paste0(save_dir,'/lf_data_pts&poly_5.25.csv'), stringsAsFactors = F)


#lf prev & year -- all
ggplot(pts, aes(x=year, y=lf_prev)) +
  geom_point(shape=16) +
  geom_smooth()

# for SSA
ggplot(pts_ssa, aes(x=year, y=lf_prev)) +
  geom_point(shape=16) +
  geom_smooth()

# for pts and poly
ggplot(pts_poly, aes(x=year, y=lf_prev)) +
  geom_point(shape=16) +
  geom_smooth()

ssa <- pts_ssa[['country']]

ptspoly_ssa <- pts_poly[country %in% ssa,]

# for pts and poly in ssa
ggplot(ptspoly_ssa, aes(x=year, y=lf_prev)) +
  geom_point(shape=16) +
  geom_smooth()