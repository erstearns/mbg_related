pacman::p_load(rgeos, data.table, raster, rgdal, INLA, seegSDM, seegMBG, plyr, sp, foreign)

fixed_effects <- "africa_lf + evi + LST_day + total_pop + aridity + elevation + lf_vector + annual_precip + PET"

# Hard-code central directories
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
central_cov_dir <- paste0(root,'/WORK/11_geospatial/01_covariates/09_MBG_covariates/') # central folder per Lucas
central_tv_covs <- c('evi','lights_new','LST_day','total_pop','urban_rural', 'land_cover') #removed U5m relevant stuff - 'rates','malaria','fertility',
central_ntv_covs <- c('access', 'Gecon','irrigation','LF','LF_vector','reservoirs','aridity','elevation','annual_precip','PET')


selected_covs <- strsplit(fixed_effects," ")
selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]
num_covs <- length(selected_covs)
lcovs <- list()
for(i in 1:num_covs) {
  message(selected_covs[i])
  this_cov <- selected_covs[i]
  if(this_cov %in% central_ntv_covs) { # Add if it is from the temporally non-varying list.
    lcovs[[i]] <- get(this_cov)
  }
  if(this_cov %in% central_tv_covs) { # Add if it is from the temporally varying list.
    lcovs[[i]] <- get(this_cov)
  }
  print(this_cov)
  print(lcovs)
  print(get(this_cov))
  names(lcovs)[[i]] <- this_cov
}