###########################################################################################
###########################################################################################
## Run MBG model for Lymphatic filariasis prevalence
## Nick Graetz's LRI code adapted by Erin Stearns
## source('/share/code/geospatial/stearns7/mbg/lf/launch.R')
## qsub -e /share/temp/sgeoutput/stearns7/errors -o /share/temp/sgeoutput/stearns7/output -cwd -pe multi_slot 30 -P proj_geospatial -N lri_raw /share/code/geospatial/ngraetz/mbg/lri/r_shell.sh /share/code/geospatial/ngraetz/mbg/lri/launch.R   
###########################################################################################
###########################################################################################

## Set repo location and indicator group
repo <- '/share/code/geospatial/stearns7/mbg/'
indicator_group <- 'lf'
indicator <- 'lf_pos'
Regions=c('LKA')

## Load libraries and miscellaneous MBG project functions.
setwd(repo)
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
package_lib <- paste0(root,'/temp/geospatial/packages') # Library for all MBG versioned packages. Ensures that none of this code is
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
#   TODO: will need to send pull request to seegMBG of github
package_list <- c('foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

## Read config file and save all parameters in memory
config <- load_config(repo = repo,
                      indicator_group = indicator_group,
                      indicator = indicator)

## Create run date in correct format
run_date <- make_time_stamp(time_stamp)

## Create directory structure for this model run
create_dirs(indicator_group = indicator_group,
            indicator = indicator)

## Load simple polygon template to model over
countries <- c('LKA')
gaul_list <- get_gaul_codes(countries)
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
                      years = 'annual')

## Add GAUL_CODE and region to df given the lat/longs. We need this to stratify in holdout function.
df <- add_gauls_regions(df = df,
                        simple_raster = simple_raster)

df <- df[start_year <= 2000, year:=2000]

## Run function to create holdouts (returns list of data.frames with an additional "folds" column)
table(df$region, df$year)
# long_col = 'longitude'
# lat_col = 'latitude'
# n_folds = as.numeric(n_folds)
# stratum_qt <- make_folds(data = df, n_folds = n_folds, spat_strat = spat_strat,
#                          temp_strat = temp_strat, strat_cols = "region",
#                          ts = 20, mb = 10)

## Set strata as character vector of each strata (in my case, just stratifying by region whereas U5M stratifies by region/age)
strata <- Regions

## ~~~~~~~~~~~~~~~~~~~~~~~~  Parallel MBG  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~ Submit job by strata/holdout  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get gbd estimates for raking
# gbd <- load_gbd_data(gbd_type = "output",
#                      gbd_name = 322,
#                      gaul_list = get_gaul_codes('africa'),
#                      measure_id = 5,
#                      age_group_id = 1,
#                      metric_id = 3)

## This is where we begin to do things by region
if(model_type == 'full')  parallel_script <- 'parallel_model_full'
if(model_type == 'child') parallel_script <- 'parallel_model_child'
if(model_type == 'raw')   parallel_script <- 'parallel_model_raw'
if(model_type == 'null')  parallel_script <- 'parallel_model_null'
slots <- as.numeric(slots)
sharedir       <- sprintf('/share/geospatial/mbg/%s/%s',indicator_group,indicator)

loopvars <- NULL
for(r in strata){
  for(holdout in 0) {
    
    qsub <- make_qsub(code = parallel_script,
                      reg = r,
                      saveimage = TRUE,
                      test = TRUE,
                      holdout = holdout)
    
    system(qsub)
    
    loopvars <- rbind(loopvars, c(r,holdout))
    
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~ Post-Estimation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## check to make sure models are done before continuing
waitformodelstofinish()

## Save strata for Shiny to use in producing aggregated fit statistics
dir.create(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/fit_stats'))
save(strata, file = paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/fit_stats/strata.RData'))

# Post estimation to be done by strata (in my case, this is just region)
# To-Do: parallelize this loop
for(reg in strata){
  message(reg)
  
  load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/', indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))
  
  ## get aggregated estimates for all admin0. Aggregate to level you rake to
  simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = TRUE)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
  ## Pull 2000-2015 annual population brick using new covariates function
  pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',                
                                                       measures = 'total',          
                                                       simple_polygon = simple_polygon,
                                                       start_year  = 2000,
                                                       end_year    = 2015,
                                                       interval_mo = 12,
                                                       agebin=1)
  pop_raster_annual  <- pop_raster_annual[[1]]
  pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
  pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
  pop_raster_annual  <- mask(pop_raster_annual, simple_raster)
  
  ## Copy GBD values across years (because GBD2015 nonfatal didn't produce annual... could interpolate later)
  # new_gbd <- list()
  # for(this_year in c(2000:2015)) {
  #   if(this_year %in% 2000:2004) copied_data <- gbd[year == 2000, ]
  #   if(this_year %in% 2005:2009) copied_data <- gbd[year == 2005, ]
  #   if(this_year %in% 2010:2014) copied_data <- gbd[year == 2010, ]
  #   if(this_year %in% 2015:2015) copied_data <- gbd[year == 2015, ]
  #   copied_data <- copied_data[, year := this_year]
  #   new_gbd[[as.character(this_year)]] <- copied_data
  # }
  # new_gbd <- rbindlist(new_gbd)
  
  ## Create population weights using the annual brick and feed custom year argument to aggregation function
  pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                          simple_raster = simple_raster,
                                          pop_raster    = pop_raster_annual,
                                          gaul_list     = get_gaul_codes(reg))
  
  cond_sim_raw_adm0 <- make_condSim(pop_wts_object = pop_wts_adm0,
                                    gaul_list      = get_gaul_codes(reg),
                                    admin_level    = 0,
                                    cell_pred      = cell_pred,
                                    summarize      = TRUE,
                                    years          = c(2000:2015))
  
  
  ## Get raking factors
  ## Try logit rake function from Bobby
  # rf <- calc_logit_raking_factors(pop_wts_object = pop_wts_adm0,
  #                                 gaul_list      = get_gaul_codes(reg),
  #                                 admin_level    = 0,
  #                                 cell_pred      = cell_pred)
  ## Get raking factors (old method)
  # rf   <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
  #                             gaul_list   = gaul_list,
  #                             rake_to     = new_gbd)
  
  # raked_cell_pred <- rake_predictions(raking_factors = rf,
  #                                     pop_wts_object = pop_wts_adm0,
  #                                     cell_pred      = cell_pred,
  #                                     logit_rake     = FALSE)
  
  ## summarize raked predictions for each cell
  mean_raster <- make_cell_pred_summary( draw_level_cell_pred = cell_pred,
                                         mask                 = simple_raster,
                                         return_as_raster     = TRUE,
                                         summary_stat         = 'median')
  
  # raked_mean_raster <- make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
  #                                              mask                 = simple_raster,
  #                                              return_as_raster     = TRUE,
  #                                              summary_stat         = 'median')
  
  #assign(sprintf('%s_rf',reg),rf)
  assign(sprintf('%s_mean_raster',reg),mean_raster)
  #assign(sprintf('%s_raked_mean_raster',reg),raked_mean_raster)
  rm(mean_raster); rm(cell_pred); rm(pop_wts) #rm(raked_mean_raster); rm(raked_cell_pred); 
  
}


# combine regions raster (child only for now)
# rf <- do.call(rbind.fill, list(name_rf,
#                                essa_rf,
#                                sssa_rf,
#                                wssa_rf,
#                                cssa_rf))
m = do.call(raster::merge,list(LKA_mean_raster))
# m_raked = do.call(raster::merge,list(name_raked_mean_raster,
#                                      essa_raked_mean_raster,
#                                      sssa_raked_mean_raster,
#                                      wssa_raked_mean_raster,
#                                      cssa_raked_mean_raster))
#save_post_est(rf,'csv','rf')
save_post_est(m,'raster','mean_raster')
#save_post_est(m_raked,'raster','mean_raked_raster')



## ~~~~~~~~~~~~~~~~~~~~~~~~ Shiny Diagnostics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~http://mbg-viz.duckdns.org:3456/~~~~~~~~~~~~~~~~~~~~~~~

# Combine model image history from each reason (these Shiny functions assume a single model image)
all_fixed_effects <- paste(fixed_effects, stacked_fixed_effects, sep=" + ") #gbd_fixed_effects, mbg_fixed_effects, 
combine_region_image_history(indicator = indicator,
                             indicator_group = indicator_group,
                             run_date = run_date,
                             fixed_effects = all_fixed_effects)

# Load combined image history
load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/', run_date, '.RData'))

# Load all Africa polygon for plotting
simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes('LKA'), buffer = 0.4, subset_only = TRUE)
subset_shape     <- simple_polygon_list[[1]]
simple_polygon   <- simple_polygon_list[[2]]

# Plot Shiny stuff
shiny_data_and_preds(gaul_list = get_gaul_codes(countries),
                     run_date = run_date,
                     indicator = indicator,
                     
                     indicator_group = indicator_group,
                     pred_file = 'lf_pos_mean_raster.tif',
                     layer_name = 'lf_pos_mean_raster.') 

# shiny_raked(gaul_list = get_gaul_codes('africa'),
#             run_date = run_date,
#             indicator = indicator,
#             indicator_group = indicator_group,
#             pred_file = 'lf_pos_mean_raked_raster.tif',
#             layer_name = 'lf_pos_mean_raked_raster.') 

shiny_cov_layers(fixed_effects = all_fixed_effects,
                 gaul_list = gaul_list,
                 run_date = run_date,
                 indicator = indicator,
                 indicator_group = indicator_group)

shiny_raking_map()

shiny_data_scatter(df = df,
                   run_date = run_date,
                   indicator = indicator,
                   indicator_group = indicator_group,
                   year_var = 'start_year')


## ~~~~~~~~~~~~~~~~~~~~~~~~ Mark model best ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~ Also would need to run Shiny functions on run_date = "best" ~~~

# save_best_model(indicator_group = indicator_group,
#                 indicator = indicator,
#                 run_date = '2016_10_16_19_27_10',
#                 pred_file = 'has_lri_mean_raster.tif')

##############################################################################
##############################################################################
##############################################################################