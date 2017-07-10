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

config <- load_config(repo = repo,
                      indicator_group = indicator_group,
                      indicator = indicator)

run_date <- make_time_stamp(time_stamp)

## Create directory structure for this model run
create_dirs(indicator_group = indicator_group,
            indicator = indicator)

pathaddin = paste0('point_only_preproc_run') #path extension to distinguish this as a preprocessing run

# temp: testing polygon code so reloading successful file - ## locate the path to the mean raster. we'll need that
output_dir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/2017_03_08_16_27_37')
path.to.risk.rast <- paste0(output_dir, '/', indicator,'_prediction_eb')

#######################
## DONE FITTING INLA ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################

all_poly_data <- fread(paste0(root,'/temp/stearns7/lf/mbg/lka/all_lka.csv')) #change for each country/source

## 4. Resample polygon to points via Roy's U5M method
all_poly_data <- all_poly_data[(shapefile!="COL_DHS_2005" & shapefile!="GTM_DHS_1987") | is.na(shapefile), ]

## we need to set up the raster dict so we can use the shapefiles of
## our choice for sampling point locations
## column 1 is entries in raster_col in df
## column 2 is path to those rasters
rd <- data.frame(tag = c(1,0), #swapped out"risk", "population"
                 raster.loc = c(path.to.risk.rast,
                                paste0(root,'/WORK/11_geospatial/01_covariates/09_MBG_covariates/WorldPop_total_global_stack.tif')))

## resample here
source('mbg_central/polygon_functions.R')
all_poly_data <- resample_polygons(data = all_poly_data,
                                       cores = 1,
                                       indic = 'lf_pos',
                                       density = 10,
                                       raster_dict = rd,
                                       raster_col = "sampling"
                                       #use_1k_popraster = FALSE - can I do this and will it overwrite the 'TRUE' that is the function default?
)


test <- all_poly_data

#
collapsed <- fread(paste0(root,'/temp/stearns7/lf/mbg/lka/all_lka.csv'), stringsAsFactors = F))

all_point_data <- collapsed[point == 1, ]
all_point_data <- all_point_data[, weight := 1]
all_point_data <- all_point_data[, pseudocluster := FALSE]
#names(all_point_data)[names(all_point_data)=="lat"] <- "latitude" - these are already named correctly
#names(all_point_data)[names(all_point_data)=="long"] <- "longitude"
#all_point_data <- all_point_data[, lf_pos_count := lf_pos * N] # not sure what this is for
all_point_data <- all_point_data[, names(all_poly_data), with=FALSE]
all_poly_data <- all_poly_data[, lf_pos := lf_pos * N]
all_collapsed <- rbind(all_point_data, all_poly_data)
names(all_collapsed)[names(all_collapsed)=='start_year'] <- 'year'
all_collapsed <- all_collapsed[, c('survey_series', 'year','country','latitude','longitude','N','lf_pos','weight'), with = FALSE]

## 5. Format and save
names(all_collapsed)[names(all_collapsed)=='start_year'] <- 'year'
all_collapsed <- all_collapsed[, c('survey_series', 'year','country','latitude','longitude','N','lf_pos','weight'), with = FALSE]
all_collapsed <- all_collapsed[, cluster_id := .GRP, by = c('year','latitude','longitude','country','survey_series')]
names(all_collapsed)[names(all_collapsed)=='survey_series'] <- 'source'
all_collapsed <- all_collapsed[order(cluster_id,latitude,longitude,lf_pos,year,source)]
## Replace year with period 1998-2002, 2003-2007, 2008-2012, 2013-2017
all_collapsed <- subset(all_collapsed, year >= 1998)
names(all_collapsed)[names(all_collapsed) == "year"] = "original_year"
all_collapsed <- all_collapsed[original_year >= 1998 & original_year <= 2002, year := 2000]
all_collapsed <- all_collapsed[original_year >= 2003 & original_year <= 2007, year := 2005]
all_collapsed <- all_collapsed[original_year >= 2008 & original_year <= 2012, year := 2010]
all_collapsed <- all_collapsed[original_year >= 2013 & original_year <= 2017, year := 2015]
all_collapsed <- all_collapsed[, latitude := as.numeric(latitude)]
all_collapsed <- all_collapsed[, longitude := as.numeric(longitude)]
all_collapsed <- all_collapsed[!is.na(latitude)]
all_collapsed <- all_collapsed[!is.na(longitude)]
all_collapsed <- all_collapsed[latitude>=-90 & latitude<=90]
all_collapsed <- all_collapsed[longitude>=-180 & longitude<=180]
all_collapsed <- all_collapsed[, lf_pos := round(lf_pos, 0)]
## In clusters where LRI > N (due to tiny samples and every child having LRI), cap at N
all_collapsed <- all_collapsed[lf_pos > N, lf_pos := N]
write.csv(all_collapsed, file = paste0(root,"/WORK/11_geospatial/10_mbg/input_data/lf_pos.csv"), row.names=FALSE)
