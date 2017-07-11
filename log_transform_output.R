##############################################
# Code author: Erin Stearns
# Code objective: log transforming raster values of MBG output
# Date: 6.15.17
#############################################

rm(list = ls()) 

# ------------------------------------------- set up environment -------------------------------------------
my_lib <- ("C:/Users/stearns7/Documents/R/win-library/3.3")
.libPaths(my_lib)
pacman::p_load(fields, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, 
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, seegMBG, sp, rgeos, maptools, rgeos, rgdal)

ras_dir <- 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/output_rasters/africa/2017_06_06_14_35_57/'

# ------------------------------------------- load data ----------------------------------------------------
mbg_out <- paste0(ras_dir, 'lf_pos_mean_raster.tif')

m <- stack(mbg_out)

# ------------------------------------------- explore a bit ------------------------------------------------
#min vals
cellStats(m, stat='min', na.rm=T)

#max vals
cellStats(m, stat='max', na.rm=T)

# ------------------------------------------- log transform -------------------------------------------------
l <- log(m)
#min vals
cellStats(l, stat='min', na.rm=T)

#max vals
cellStats(l, stat='max', na.rm=T)

writeRaster(l, paste0(ras_dir, '/log_trans_lf_pos_mean_raster'), 'GTiff')













