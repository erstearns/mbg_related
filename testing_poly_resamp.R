#################################
# code author: erin stearns
# code intent: polygon resampling
#################################


## Set repo location and indicator group
repo <- '/share/code/geospatial/stearns7/mbg/'
indicator_group <- 'lf'
indicator <- 'lf_pos'


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
gaul_list <- c('94')
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                           buffer = 0.4)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]


# -------------------------- loading polygon data ------------------------
#data location
folder <- paste0(root, 'temp/stearns7/lf/mbg/resample_data')
csv <- paste0(folder, '/ghana_polys.csv')

# polygon data
ghana_polys <- fread(file = paste0(root, '/temp/stearns7/lf/mbg/resample_data/ghana_polys.csv'), stringsAsFactors = F)

all_poly_data <- copy(ghana_polys)

source('lf/poly_resamp.R')

all_poly_data <- resample_polygons_pop(data = all_poly_data,
                                   cores = 1,
                                   indic = 'lf_pos',
                                   density = 0.001, #this is where we can change things up to determine the number of points that get sprinkled
                                   use_1k_popraster = TRUE
)


write.csv(all_poly_data, file = paste0(root, '/temp/stearns7/lf/mbg/resample_data/ghana_polys_resamp0.001.csv'))


# --------------------- trying risk and pop resampling -------------------


#cano environmental suitability raster
path.to.risk.rast <- (paste0(root, 'WORK/11_geospatial/01_covariates/00_MBG_STANDARD/lf/mean/synoptic/lf_mean_synoptic.tif'))


# polygon data
ghana_polys <- fread(file = paste0(root, '/temp/stearns7/lf/mbg/resample_data/ghana_polys_risk.csv'), stringsAsFactors = F)

## we need to set up the raster dict so we can use the shapefiles of
## our choice for sampling point locations
## column 1 is entries in raster_col in df
## column 2 is path to those rasters
rd <- data.frame(tag = c(1,0), #swapped out"risk", "population"
                 raster.loc = c(path.to.risk.rast,
                                paste0(root,'/WORK/11_geospatial/01_covariates/09_MBG_covariates/WorldPop_total_global_stack.tif')))
rd$raster.loc <- as.character(rd$raster.loc)

all_poly_data <- copy(ghana_polys)

source('lf/poly_resamp.R')

all_poly_data <- resample_polygons(data = all_poly_data,
                                   cores = 1,
                                   indic = 'lf_pos',
                                   density = 0.1, #this is where we can change things up to determine the number of points that get sprinkled
                                   raster_dict = rd,
                                   raster_col = "sampling",
                                   use_1k_popraster = TRUE,
                                   regions = gaul_list
)                                   

write.csv(all_poly_data, file = paste0(root, '/temp/stearns7/lf/mbg/resample_data/ghana_polys_resamp_risk_0.1.csv'))






