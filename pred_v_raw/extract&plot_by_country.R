##############################################
# Code author: Erin Stearns
# Code objective: investigating relationship between raw and predicted LF prevalence
# Date: 6.21.17
#############################################

rm(list = ls()) 

# ------------------------------------------- set up environment -------------------------------------------
my_lib <- ("C:/Users/stearns7/Documents/R/win-library/3.3")
.libPaths(my_lib)
pacman::p_load(fields, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, ggmap, devtools, scale,
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, seegMBG, sp, rgeos, maptools, rgeos, rgdal)
source('C:/Users/stearns7/Documents/git_repos/gen_fxns/geo_fxns/extract_ras2pts.R')
source('C:/Users/stearns7/Documents/git_repos/gen_fxns/multiplot_fxn.R')

data_dir <- 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/investigations/prediction_v_raw/7_10_2017_run/'
shape_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/spatial_data/africa_0/')

oth_dir <- 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/output_rasters/africa/2017_06_06_14_35_57/'
data_dir2 <- 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/investigations/prediction_v_raw/6_6_2017_run/'

# ------------------------------------------- load data ----------------------------------------------------
#africa admin0 shapefile
afr <- shapefile(paste0(shape_dir,'africa_adm0.shp'))

# ---------- 7/10 run -------------
#input data
lf <- fread(file = paste0(data_dir, 'input_data.csv'), stringsAsFactors = F)
lf_pts <- lf[point==1,]

#load mbg output raster from 7/10 run
#pred <- brick(paste0(data_dir, 'lf_pos_all_dat_mean_raster_full.tif'))

# --------- 6/6 run --------------
#load data from 6/6 run
model2pts <- fread(file = paste0(ras_dir, 'raster_vals_to_original_data_pts.csv'), stringsAsFactors = F)
model2pts <- model2pts[RASTERVALU != -9999,]
model2pts[,lf_prev:=lf_pos/N]

# formatting 6/6 data
setnames(model2pts, 'RASTERVALU', 'old_pred_prev')
old <- model2pts[,old_diff:=old_pred_prev - lf_prev]

#load mbg output raster from 6/6 run
pred <- brick(paste0(oth_dir, 'lf_pos_mean_raster.tif'))

# ------------------------------------------- extracting ras to points ----------------------------------------------------
#set data set to use:
lf_pts <- copy(old)

#need to make period map
period_map <- make_period_map(modeling_periods = c(2000:2015))

#need to place raster in a list
pred_ras <- c('pred')

#copy point data over
df <- copy(lf_pts)

#extract the rasters by points
locations = SpatialPoints(coords = as.matrix(df[,.(longitude,latitude)]), proj4string = CRS(proj4string(pred)))

pred_values = as.data.frame(extract(pred, locations))

#combine the dataset with the covariate extractions
df = cbind(df, pred_values)

#reshape to fill out the columns-- these just get the tv cov names in a way that is needed for two different efforts
#make this less likey to cause a hiccup
pred_colnames = names(pred) #unlisted
pred_colist = lapply(pred_colnames, function(x) grep(x, names(df), value = T))

df1 <- copy(df)

#getting pred from corresponding data year

df1[year==2000,pred_prev:=lf_pos_mean_raster.1]
df1[year==2001, pred_prev:=lf_pos_mean_raster.2]
df1[year==2002, pred_prev:=lf_pos_mean_raster.3]
df1[year==2003, pred_prev:=lf_pos_mean_raster.4]
df1[year==2004, pred_prev:=lf_pos_mean_raster.5]
df1[year==2005, pred_prev:=lf_pos_mean_raster.6]
df1[year==2006, pred_prev:=lf_pos_mean_raster.7]
df1[year==2007, pred_prev:=lf_pos_mean_raster.8]
df1[year==2008, pred_prev:=lf_pos_mean_raster.9]
df1[year==2009, pred_prev:=lf_pos_mean_raster.10]
df1[year==2010, pred_prev:=lf_pos_mean_raster.11]
df1[year==2011, pred_prev:=lf_pos_mean_raster.12]
df1[year==2012, pred_prev:=lf_pos_mean_raster.13]
df1[year==2013, pred_prev:=lf_pos_mean_raster.14]
df1[year==2014, pred_prev:=lf_pos_mean_raster.15]
df1[year==2015, pred_prev:=lf_pos_mean_raster.16]


df2 <- df1[!is.na(pred_prev),]
df2 <- df2[,c(2:15,32), with=FALSE]

#validating observed prevalence field
df2[,lf_prev:=lf_pos/N]

#creating diff field
df2[,diff:=pred_prev-lf_prev]

#save
#write.csv(df2, file=paste0(data_dir2, 'preds_extracted2inputdata.csv'))

#run through plotting code

# ------------------------------------- creating time series plots by country -------------

plot_dir <- paste0(data_dir2,'country_plots/')

#joining to get actual country name to an iso3 key
iso3 <- fread(file='C:/Users/stearns7/OneDrive - UW Office 365/spatial_data/mbg_gaul_iso3_key_lf_adapted.csv',stringsAsFactors = F)
iso3_mod <- iso3[,c(2, 5), with=F]

df3 <- copy(df2)
df3 <- merge(df3, iso3_mod, by.x = 'country', by.y='iso3')

setnames(df3, 'short_name', 'name')
df3[,observed_prevalence:=lf_pos/N]
df3[,predicted_prevalence:=pred_prev]

#both

p.list = lapply(sort(unique(df3$name)), function(i){
  basic_obs <- ggplot(df3[name==i], aes(x=year, y=observed_prevalence)) +
    geom_point(shape=16) + 
    facet_wrap(~as.factor(name)) +
    labs(x = 'Year', y = 'Observed Prevalence') +
    scale_x_continuous(limits = range(df3$year)) + 
    scale_y_continuous(limits = range(df3$observed_prevalence))
  loes_obs <- basic_obs + geom_smooth()
  
  basic_pred <- ggplot(df3[name==i], aes(x=year, y=predicted_prevalence)) +
    geom_point(shape=16) + 
    facet_wrap(~as.factor(name)) +
    labs(x = 'Year', y = 'Predicted Prevalence') +
    scale_x_continuous(limits = range(df3$year)) + 
    scale_y_continuous(limits = range(df3$predicted_prevalence))
  loes_pred <- basic_pred + geom_smooth()
  
  basic_comp <- ggplot(df3[name==i], aes(x=observed_prevalence, y=predicted_prevalence)) +
    geom_point(shape=16) + 
    facet_wrap(~as.factor(name)) +
    labs(x = 'Observed Prevalence', y = 'Predicted Prevalence') +
    scale_x_continuous(limits = range(df3$observed_prevalence)) + 
    scale_y_continuous(limits = range(df3$predicted_prevalence))
  loes_comp <- basic_comp + geom_smooth()
  
  #ggsave(loes_obs, filename = paste0(ras_dir, 'by_country_plots/', i, 'loes_obs.png'))
  #ggsave(loes_pred, filename = paste0(ras_dir, 'by_country_plots/', i, 'loes_pred.png'))
  #ggsave(loes_comp, filename = paste0(ras_dir, 'by_country_plots/', i, 'loes_comp.png'))
  
  #all <- multiplot(loes_obs, loes_pred, loes_comp, cols=2)
  
  png(filename = paste0(plot_dir, i, '_obs_v_pred.png'))
  multiplot(loes_obs, loes_pred, loes_comp, cols=2)
  dev.off()
  
})

#add year bins

df3$year_bin <- ifelse(df3$year < 2004, '2000',
                        ifelse((df3$year > 2003 & df3$year < 2008), '2005',
                               ifelse(df3$year > 2007 & df3$year < 2012, '2010', '2015')))

# basic maps

x <- get_map(location = 'africa', zoom=3, color="bw")

ggmap(x) +
  geom_point(data=df3, aes(x=longitude, y=latitude, fill = diff), colour='black',pch=21) +
  scale_fill_gradient2(low=('blue'), mid="white", high=('red'), midpoint = 0) + 
  facet_wrap(~year_bin, ncol = 2)


# -------------------------- want to look at trends overall -----------------------------
