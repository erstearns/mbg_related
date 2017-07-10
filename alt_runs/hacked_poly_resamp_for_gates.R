#hacked poly resamp to get going for gates

#################################
# code author: aaron osgood-zmmerman, modified by erin stearns
# code intent: double wrapper for lf data that needs to have polygons sampled by pop and or by risk
#################################


## Set repo location and indicator group ------------------------------
repo <- '/share/code/geospatial/stearns7/mbg/'
indicator_group <- 'lf'
indicator <- 'lf_pos'

## Load libraries and miscellaneous MBG project functions.
setwd(repo)
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")

#if on a geos node:
#package_lib <- paste0(root, '/temp/geospatial/geos_packages')

#if on cluster prod
package_lib <- paste0(root,'/temp/geospatial/packages') # Library for all MBG versioned packages. Ensures that none of this code is
##    dependent on the machine where the user runs the code.
.libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library().
##    Necessary for seeg libraries.
source('mbg_central/mbg_functions.R')                   # Functions to run MBG model.
source('mbg_central/prep_functions.R')                  # Functions to setup MBG run
source('mbg_central/covariate_functions.R')             # Functions to prep and transform 5*5 covariates
source('mbg_central/misc_functions.R')                  # Other computational MBG-related functions.
source('mbg_central/post_estimation_functions.R')
source('mbg_central/gbd_functions.R')
source('mbg_central/shiny_functions.R')
source('mbg_central/holdout_functions.R')
source('mbg_central/validation_functions.R')
source('mbg_central/categorical_variable_functions.R')
source('mbg_central/stacking_functions.R') ###
source('mbg_central/polygon_functions.R')
source('mbg_central/seegMBG_transform_functions.R')     # Using Roy's edit for now that can take temporally varying covariates,
#   TODO: will need to send pull request to seegMBG of github
package_list <- c('survey', 'foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr', 'foreach', 'doParallel')
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

#path extension to distinguish this as a preprocessing run
pathaddin = paste0('point_only_preproc_run') 

## Load simple polygon template to model over
gaul_list <- c('231')
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                           buffer = 0.4)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]


#######################
## DONE FITTING INLA ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################

# Grabbing 'risk' raster----------------
source('lf/poly_resamp.R')

#temp for gates run 
path.to.risk.rast <- (paste0(root,'temp/stearns7/lf/mbg/lka/sri_lanka.tif'))

#poly data
all_poly_data <- fread(file = paste0(root, 'temp/stearns7/lf/mbg/lka/lka_poly_4.7.csv'), stringsAsFactors = F)

## we need to set up the raster dict so we can use the shapefiles of
## our choice for sampling point locations
## column 1 is entries in raster_col in df
## column 2 is path to those rasters
rd <- data.frame(tag = c(1,0), #swapped out"risk", "population"
                 raster.loc = c(path.to.risk.rast,
                                paste0(root,'/WORK/11_geospatial/01_covariates/09_MBG_covariates/WorldPop_total_global_stack.tif')))
rd$raster.loc <- as.character(rd$raster.loc)

## resample here

all_poly_data <- resample_polygons(data = all_poly_data,
                                   cores = 1,
                                   indic = 'lf_pos',
                                   density = 10, #this is where we can change things up to determine the number of points that get sprinkled
                                   raster_dict = rd,
                                   raster_col = "sampling",
                                   use_1k_popraster = TRUE # - can I do this and will it overwrite the 'TRUE' that is the function default?
)

#all_poly_data <- all_poly_data[, lf_pos := lf_pos_count / N] #not sure what the point of this is

test <- all_poly_data

# get point data - 'point' data in all_poly_data does not have original values for lf_pos and N
all_point_data <- all_collapsed

#names(all_point_data)[names(all_point_data)=="lat"] <- "latitude" - these are already named correctly
#names(all_point_data)[names(all_point_data)=="long"] <- "longitude"
#all_point_data <- all_point_data[, lf_pos_count := lf_pos * N] # not sure what this is for
all_point_data <- all_point_data[, names(all_poly_data), with=FALSE]
#all_poly_data <- all_poly_data[, lf_pos := lf_pos * N]

all_collapsed_pp <- rbind(all_point_data, all_poly_data)

############################ already formatted prior to ers
## 5. Format and save
names(all_collapsed)[names(all_collapsed)=='start_year'] <- 'year'
all_collapsed <- all_collapsed[, c('survey_series', 'year','country','latitude','longitude','N','lf_pos','weight'), with = FALSE]
all_collapsed <- all_collapsed[, cluster_id := .GRP, by = c('year','latitude','longitude','country','survey_series')]
names(all_collapsed)[names(all_collapsed)=='survey_series'] <- 'source'
##########################################


