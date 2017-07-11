#############################################
# Author: Erin Stearns
# Objective: Pull out BRT output for examination
# Date: 6/20/2017
############################################
## Set repo location and indicator group
repo <- '/share/code/geospatial/stearns7/mbg/'
indicator_group <- 'lf'
indicator <- 'lf_pos'
Regions = c('wssa_lf', 'cssa_lf', 'essa_lf') #c('lf_endem_afr') # #created new region - essa_lf - , added zim to essa; 'essa','name','wssa','cssa'; trying new region of all lf endemic africa countries

## Load libraries and miscellaneous MBG project functions.
setwd(repo)
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                      paste0(root,'temp/geospatial/geos_packages'),
                      paste0(root,'temp/geospatial/packages'))
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

# ers functions
source('lf/check_loc_results.R')

package_list <- c('foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

rd <- 'XXXX' #TO-DO input run date of desired output to examine
ind <- indicator # indicator
ind_gp <- indicator_group # indicator_group
reg <- Regions
age = 0
holdout = 0
vallevel = "" ## validation level if used in qsub

load(paste0(root, 'temp/stearns7/lf/mbg/lf_mbg_model_runs/2017_06_06_14_35_57/child_model_list_', reg, '_', holdout))



for (r in reg){
  load(paste0('/share/geospatial/mbg/', ind_gp, '/', ind, '/output/', rd,
            sprintf('/child_model_list_%s_%i.RData', reg, holdout)))
  smo <- child_models
  gbm <- smo[['gbm']]
