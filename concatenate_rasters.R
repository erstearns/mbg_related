#########################################
# code author: erin stearns
# code intent: mangaged to get pieces of a whole to run, now need to make into 1 raster
# date: 4.18.17
########################################


# ---------------- load libraries and source code --------------------------------------------------
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
package_lib <- paste0(root,'/temp/geospatial/packages') # Library for all MBG versioned packages. Ensures that none of this code is
#    dependent on the machine where the user runs the code.
.libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library(). 
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
source('mbg_central/seegMBG_transform_functions.R')   
package_list <- c('foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

# --------------- set locations ------------------------------------------------------------------
ras.dir <- ('/home/j/temp/stearns7/lf/mbg/lf_mbg_model_runs/pcs_4.18/')


# --------------- read in each raster ------------------------------------------------------------

name_mean_raster <- raster(paste0(ras.dir, 'lf_pos_prediction_eb_bin0_name_0.grd'))
essa_mean_raster <- raster(paste0(ras.dir, 'lf_pos_prediction_eb_bin0_essa_0.grd'))
wssa_mean_raster <- raster(paste0(ras.dir, 'lf_pos_prediction_eb_bin0_wssa_0.grd'))
csssa_mean_raster <- raster(paste0(ras.dir, 'lf_pos_prediction_eb_bin0_csssa_0.grd'))

#merge

Regions=c('cssa','essa','name', 'wssa') #csssa doesn't look right

strata <- Regions

for(reg in strata){
  message(reg)

  load(paste0('/home/j/temp/stearns7/lf/mbg/lf_mbg_model_runs/pcs_4.18/lf_pos_cell_draws_eb_bin0_', reg, '_0.RData'))

  ## get aggregated estimates for all admin0. Aggregate to level you rake to
  simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = TRUE)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]

  mean_raster <- make_cell_pred_summary( draw_level_cell_pred = cell_pred,
                                         mask                 = simple_raster,
                                         return_as_raster     = TRUE,
                                         summary_stat         = 'median')

  assign(sprintf('%s_mean_raster',reg),mean_raster)
  rm(mean_raster); rm(cell_pred) #rm(pop_wts) #rm(raked_mean_raster); rm(raked_cell_pred);

}


#  ------ TO DO: CHANGE HERE IF CHANGING REGIONS TO MERGE

m = do.call(raster::merge,list(name_mean_raster,
                               essa_mean_raster,
                               wssa_mean_raster,
                               cssa_mean_raster))

#write to a raster
writeRaster(
  m,
  file = paste0(root, 'temp/stearns7/lf/mbg/lf_mbg_model_runs/pcs_4.18/africa_sans_sssa_mean_raster'),
  format='GTiff',
  overwrite = TRUE)
