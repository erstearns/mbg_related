###########################################################################################
# author: erin stearns
# objective: investigate covariate trends
# date: 10 July 2017
###########################################################################################

## Set repo location
repo <- '/share/code/geospatial/stearns7/mbg/'
indicator_group <- 'lf'
indicator <- 'lf_pos'

## Load libraries and miscellaneous MBG project functions.
setwd(repo)
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                      paste0(root,'temp/geospatial/geos_packages'),
                      ('/share/code/geospatial/stearns7/mbg_pkgs_conda'))
#    dependent on the machine where the user runs the code.
.libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library().
#    Necessary for seeg libraries.
source('mbg_central/mbg_functions.R')                   # Functions to run MBG model.
source('mbg_central/prep_functions.R')                  # Functions to setup MBG run
source('mbg_central/covariate_functions.R')             # Functions to prep and transform 5*5 covariates
source('mbg_central/misc_functions.R')                  # Other computational MBG-related functions.
source('mbg_central/post_estimation_functions.R')
source('mbg_central/gbd_functions.R')
source('mbg_central/shiny_functions.R')
source('mbg_central/holdout_functions.R')
source('mbg_central/categorical_variable_functions.R')
source('mbg_central/validation_functions.R')
source('mbg_central/seegMBG_transform_functions.R')     # Using Roy's edit for now that can take temporally varying covariates,
source('mbg_central/stacking_functions.R')

# ers functions
source('lf/check_loc_results.R')
source('lf/misc_fxns_lf.R')

package_list <- c('foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

install.packages("rasterVis")
install.packages("ggplot2")
install.packages('scales')
library(rasterVis)
library(ggplot2)
library(scales)
# --------------------------------------- loading and formatting everything  ----------------------------------

config <- load_config(repo = repo,
                      indicator_group = indicator_group,
                      indicator = indicator)

## bring in run date from parent script
run_date <- make_time_stamp(time_stamp)

## Create directory structure for this model run
create_dirs(indicator_group = indicator_group,
            indicator = indicator)

#years to use
year_list <- eval(parse(text=year_list))
## Define modeling space (right now, just in years)
period_map <- make_period_map(modeling_periods = c(min(year_list):max(year_list)))

#start loop here
countries <- c('bfa', 'cmr', 'gha')
for (c in countries){
      ## Load simple polygon template to model over
      gaul_list <- get_gaul_codes(c) # TO DO ------ CHANGE HERE TO COUNTRY/i if loop bfa, cmr, gha---------------------
      simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                                 buffer = 0.4)
      subset_shape   <- simple_polygon_list[[1]]
      simple_polygon <- simple_polygon_list[[2]]
      raster_list    <- build_simple_raster_pop(subset_shape)
      simple_raster  <- raster_list[['simple_raster']]
      pop_raster     <- raster_list[['pop_raster']]

      ## Load input data and make sure it is cropped to modeling area
      df <- load_input_data(indicator = indicator,
                            simple = simple_polygon,
                            removeyemen = TRUE,
                            pathaddin = pathaddin,
                            years = 'annual',
                            withtag = TRUE,
                            datatag = datatag,
                            use_share = FALSE)

      ## Add GAUL_CODE and region to df given the lat/longs. We need this to stratify in holdout function.
      df <- add_gauls_regions(df = df,
                             simple_raster = simple_raster)


      #make covariates conditional
      cov_layers = NULL

      #Load in covariates
      # Pull all covariate bricks/layers
      if(nchar(fixed_effects)> 0){
        selected_fixed_effects <- strsplit(fixed_effects," ")
        selected_fixed_effects <- selected_fixed_effects[[1]][selected_fixed_effects[[1]] != "+"]
        selected_measures <- strsplit(fixed_effects_measures," ")
        selected_measures <- selected_measures[[1]][selected_measures[[1]] != "+"]
        cov_layers <- load_and_crop_covariates_annual(covs = selected_fixed_effects,
                                                   measures = selected_measures,
                                                   simple_polygon = simple_polygon,
                                                   start_year  = min(year_list),
                                                   end_year    = max(year_list),
                                                   interval_mo = 12,
                                                  agebin=1)
      }

      #setting up
      all_cov_layers <- cov_layers
      all_fixed_effects = paste(names(all_cov_layers), collapse = " + ")
      the_covs = format_covariates(all_fixed_effects)

      #copy the dataset to avoid unintended namespace conflicts
      the_data = copy(df)

      #extract covariates to the points and subset data where its missing covariate values
      cs_covs = extract_covariates(the_data, all_cov_layers, return_only_results = T, centre_scale = T, period_var = 'year', period_map = period_map)

      the_data = cbind(the_data, cs_covs[[1]])
      covs_cs_df = cs_covs[[2]]

      the_data = na.omit(the_data, c(indicator, 'N', the_covs)) #this will drop rows with NA covariate values

      write.csv(the_data, file = paste0(root, 'temp/stearns7/lf/covariates/covs_extracted_to_pt_data/', run_date, '_covs2data.csv')) #add country bit here

      # --------------------------------------- running analysis ----------------------------------
      #separate static and time-varying covariates
      ntv_var_names <- c('access','irrigation', 'elevation', 'distrivers', 'growingseason')
      tv_vars <- all_cov_layers[!names(all_cov_layers) %in% ntv_var_names]
      ntv_vars <- all_cov_layers[names(all_cov_layers) %in% ntv_var_names]

      names(tv_vars) <- c('CRUTS_Aridity', 'CRUTS_DiurnalTemperatureRange', 'CRUTS_FrostDayFrequency', 'CRUTS_PET', 'CRUTS_Precipitation', 'CRUTS_AveDailyMinTemp',
                          'CRUTS_AveDailyMeanTemp', 'CRUTS_AveDailyMaxTemp', 'CRUTS_WetDayFrequency', 'NightTimeLights', 'EVI', 'Urbanicity', 'MODIS_AveLandSurfaceTemp',
                          'MODIS_DaytimeLandSurfaceTemp', 'MODIS_DiurnalDifferenceInLandSurfaceTemp', 'MODIS_NighttimeLandSurfaceTemp', 'NEX_NDVI', 'MODIS_TCB', 'MODIS_TCW', 'Worldpop',
                          'ITN_Coverage', 'MDA')

      years <- c(2000:2015)

      names(tv_vars[[1]]) <- as.character(c(paste0('Year:',2000:2015)))

      theme_set(theme_bw())
      gplot(tv_vars[['crutsard']]) + geom_tile(aes(fill = value)) +
        facet_wrap(~ variable) +
        scale_fill_gradient2(low = muted('blue'), mid = ("white"), high = muted('red'), midpoint = 0.5, na.value = "grey50") +
        coord_equal()


      #plot each covariate band
      for (i in 1:(length(tv_vars))){

        names(tv_vars[[i]]) <- as.character(c(paste0('Year:',2000:2015)))

        cov_plot <- gplot(tv_vars[[i]]) + geom_tile(aes(fill = value)) +
          facet_wrap(~ variable) +
          scale_fill_gradient(low = "white", high = "red", na.value = "grey50") +
          coord_equal() +
          ggtitle(paste0(names(tv_vars)[[i]], ' in ', c))

        ggsave(paste0(c, "_", (names(tv_vars)[[i]]),".png"), cov_plot, path = paste0(root, 'temp/stearns7/lf/covariates/'))

      }

      # plot time invariant covs
      for (i in 1:(length(ntv_vars))){
          ntv <- gplot(ntv_vars[[i]]) + geom_tile(aes(fill = value)) +
            scale_fill_gradient(low = "white", high = "red", na.value = "grey50") +
            coord_equal() +
            ggtitle(paste0(names(ntv_vars)[[i]], ' in ', c))

          ggsave(paste0(c, "_", (names(ntv_vars)[[i]]),".png"), ntv, path = paste0(root, 'temp/stearns7/lf/covariates/'))

      }
}
