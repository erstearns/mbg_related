##############################################
# Code author: Erin Stearns
# Code objective: investigating relationship between raw and predicted LF prevalence
# Date: 6.21.17
#############################################

rm(list = ls()) 

# ------------------------------------------- set up environment -------------------------------------------
my_lib <- ("C:/Users/stearns7/Documents/R/win-library/3.3")
.libPaths(my_lib)
pacman::p_load(fields, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, 
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, seegMBG, sp, rgeos, maptools, rgeos, rgdal)

ras_dir <- 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/output_rasters/africa/2017_06_06_14_35_57/'

# ------------------------------------------- load data ----------------------------------------------------
model2pts <- fread(file = paste0(ras_dir, 'raster_vals_to_original_data_pts.csv'), stringsAsFactors = F)
model2pts <- model2pts[RASTERVALU != -9999,]
model2pts[,lf_prev_o:=lf_pos/N]

logmodel2pts <- fread(file = paste0(ras_dir, 'log_raster_vals_to_original_data_pts.csv'), stringsAsFactors = F)
logmodel2pts <- logmodel2pts[RASTERVALU != -9999,]

# -------------------------------------- quick viz for f3 mtg 6/21 -------------------------

basic <- ggplot(model2pts, aes(x=lf_prev_o, y=RASTERVALU)) +
  geom_point(shape=16) + 
  labs(x = 'Original LF Prevalence', y = 'Predicted LF Prevalence')
loes <- basic + geom_smooth()
by_year <- basic + geom_point(aes(colour = factor(year)), size=2)
ggsave(basic, filename = paste0(plot_dir, 'basic_', p, '.png'))
ggsave(loes, filename = paste0(plot_dir, 'loes_', p, '.png'))
ggsave(by_year, filename = paste0(plot_dir, 'by_year_', p, '.png'))

# -------------------------------------- creating a diff field -----------------------------
a <- copy(model2pts)

a[,diff:=(RASTERVALU-lf_prev_o)]

# ------------------------------------- creating time series plots by country -------------
b <- copy(a)

# for just raw and pred separately over time by country
b[,predicted_prevalence:=RASTERVALU]
b[,observed_prevalence:=lf_prev_o]

write.csv(b, file=paste0(ras_dir, 'updated_raw_pred.csv'))

#joining to get actual country name to an iso3 key
iso3 <- fread(file='C:/Users/stearns7/OneDrive - UW Office 365/spatial_data/iso3.csv',stringsAsFactors = F)
iso3_mod <- iso3[,c(1, 5), with=F]
iso3_mod[,iso3:=`ISO3166-1-Alpha-3`]
iso3_mod <- iso3_mod[,-2]

b2 <- merge(b, iso3_mod, by.x = 'country', by.y='iso3')
b2[country=='CIV', name := "Cote D'Ivoire"]

b2 <- b2[name != 'Fiji',]

plotvars <- names(b2[,c(17:18)])

print(summary(b$lf_prev_o))

plot_dir <- 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/output_rasters/africa/2017_06_06_14_35_57/country_plots/'

#both

p.list = lapply(sort(unique(b2$name)), function(i){
  basic_obs <- ggplot(b2[name==i], aes(x=year, y=observed_prevalence)) +
    geom_point(shape=16) + 
    facet_wrap(~as.factor(name)) +
    labs(x = 'Year', y = 'Observed Prevalence') +
    scale_x_continuous(limits = range(b2$year)) + 
    scale_y_continuous(limits = range(b2$observed_prevalence))
  loes_obs <- basic_obs + geom_smooth()
  
  basic_pred <- ggplot(b2[name==i], aes(x=year, y=predicted_prevalence)) +
    geom_point(shape=16) + 
    facet_wrap(~as.factor(name)) +
    labs(x = 'Year', y = 'Predicted Prevalence') +
    scale_x_continuous(limits = range(b2$year)) + 
    scale_y_continuous(limits = range(b2$predicted_prevalence))
  loes_pred <- basic_pred + geom_smooth()
  
  basic_comp <- ggplot(b2[name==i], aes(x=observed_prevalence, y=predicted_prevalence)) +
    geom_point(shape=16) + 
    facet_wrap(~as.factor(name)) +
    labs(x = 'Observed Prevalence', y = 'Predicted Prevalence') +
    scale_x_continuous(limits = range(b2$observed_prevalence)) + 
    scale_y_continuous(limits = range(b2$predicted_prevalence))
  loes_comp <- basic_comp + geom_smooth()
  
  #ggsave(loes_obs, filename = paste0(ras_dir, 'by_country_plots/', i, 'loes_obs.png'))
  #ggsave(loes_pred, filename = paste0(ras_dir, 'by_country_plots/', i, 'loes_pred.png'))
  #ggsave(loes_comp, filename = paste0(ras_dir, 'by_country_plots/', i, 'loes_comp.png'))
  
  #all <- multiplot(loes_obs, loes_pred, loes_comp, cols=2)
  
  png(filename = paste0(plot_dir, i, '_obs_v_pred.png'))
  multiplot(loes_obs, loes_pred, loes_comp, cols=2)
  dev.off()
  
})









#scrap
for(var in unique(b$country)){
  for (p in plotvars) {
    basic <- ggplot(b[country==var,], aes(x=b[country==var, year], y=b[country==var,p])) +
      geom_point(shape=16) + 
      labs(x = 'Year', y = (p))
    loes <- basic + geom_smooth()
    #by_year <- basic + geom_point(aes(colour = factor(year)), size=2)
    ggsave(basic, filename = paste0(ras_dir, 'by_country_plots/basic_',var, '_', p, '.png'))
    ggsave(loes, filename = paste0(ras_dir, 'by_country_plots/loes_', var, '_', p, '.png'))
    #ggsave(by_year, filename = paste0(ras_dir, 'by_year_', var, '_', p, '.png'))
  }
}

for (p in plotvars) {
  basic <- ggplot(b, aes(x=year, y=b[[p]])) +
    geom_point(shape=16) + 
    facet_grid(~as.factor(country)) +
    labs(x = 'Year', y = (p))
  loes <- basic + geom_smooth()
  #by_year <- basic + geom_point(aes(colour = factor(year)), size=2)
  ggsave(basic, filename = paste0(ras_dir, 'by_country_plots/basic_',var, '_', p, '.png'))
  ggsave(loes, filename = paste0(ras_dir, 'by_country_plots/loes_', var, '_', p, '.png'))
  #ggsave(by_year, filename = paste0(ras_dir, 'by_year_', var, '_', p, '.png'))
}

for (p in plotvars) {
  p.list = lapply(sort(unique(b2$name)), function(i){
    basic <- ggplot(b2[name==i], aes(x=year, y=b[[p]])) +
      geom_point(shape=16) + 
      facet_wrap(~as.factor(name)) +
      labs(x = 'Year', y = (p))
    scale_x_continuous(limits = range(b2$year))
    scale_y_continuous(limits = range(b2[[p]]))
    loes <- basic + geom_smooth()
    #by_year <- basic + geom_point(aes(colour = factor(year)), size=2)
    ggsave(basic, filename = paste0(ras_dir, 'by_country_plots/basic_',i, '_', p, '.png'))
    ggsave(loes, filename = paste0(ras_dir, 'by_country_plots/loes_', i, '_', p, '.png'))
    #ggsave(by_year, filename = paste0(ras_dir, 'by_year_', var, '_', p, '.png'))
  })
}





