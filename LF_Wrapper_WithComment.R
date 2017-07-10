###########################################################################################
###########################################################################################
## Run MBG model for LF prevalence in Africa
## Code author: Nick Graetz, modified for LF specific-model by Erin Stearns 
###########################################################################################
###########################################################################################

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
source('mbg_central/seegMBG_transform_functions.R')     # Using Roy's edit for now that can take temporally varying covariates,
#   TODO: will need to send pull request to seegMBG of github
source('mbg_central/stacking_functions.R')

package_list <- c('rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr')
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

## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load simple polygon template to model over
gaul_list <- get_gaul_codes('africa')
simple_polygon <- load_simple_polygon(gaul_list = gaul_list,
                                      buffer = 0.4)

## Load list of raster inputs (pop and simple)
raster_list <- build_simple_raster_pop(subset_shape)
simple_raster <- raster_list[['simple_raster']]
pop_raster <- raster_list[['pop_raster']]

## Load input data and make sure it is cropped to modeling area
df <- load_input_data(indicator = indicator,
                      simple = simple_polygon)

## Build spatial mesh over modeling area
mesh_s <- build_space_mesh(d = df,
                           simple = simple_polygon,
                           max_edge = mesh_s_max_edge,
                           mesh_offset = mesh_s_offset)

## Build temporal mesh (standard for now)
mesh_t <- build_time_mesh()

## Load covariates and crop
cov_layers <- load_and_crop_covariates(fixed_effects = fixed_effects, 
                                       simple_raster = simple_raster)

#ES - not using: gbd_cov_layers <- load_gbd_covariates(gbd_fixed_effects, c(2000,2005,2010,2015), gaul_list)

#ES - not using: mbg_cov_layers <- load_mbg_covariates(mbg_fixed_effects, simple_raster)

## Combine all the covariate layers we want to use and combine our fixed effects formula
#ES - Only interested in cov_layers: all_cov_layers <- c(cov_layers, gbd_cov_layers, mbg_cov_layers)
#ES - Only interested in cov_layers: all_fixed_effects <- paste(fixed_effects, gbd_fixed_effects, mbg_fixed_effects, sep=" + ")

all_cov_layers <- c(cov_layers) 
all_fixed_effects <- paste(fixed_effects, sep=" + ")


###STACK THE THINGS###

if(length(all_fixed_effects)==1){
  #split back into parts
  the_covs = unlist(tstrsplit(all_fixed_effects,"\\+"))
}

#copy the dataset to avoid unintended namespace conflicts
the_data = copy(df)

#shuffle the data into six folds
the_data = the_data[sample(nrow(the_data)),]
the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=10,labels=FALSE)] #make six folds

#extract covariates to the points and subset data where its missing covariate values
the_data = extract_covariates(the_data, all_cov_layers)
the_covs = trimws(the_covs)
the_data = na.omit(the_data, c(indicator, 'N', the_covs)) #this will drop rows with NA covariate values

#Fit a gam model
gam_results = fit_gam_child_model(df = the_data, #data frame
                                  model_name = 'gam', #identifier for the child model-- needs to be consistent through the steps
                                  fold_id_col = 'fold_id',
                                  covariates = all_fixed_effects, #rhs formula
                                  additional_terms = 'year', #column(s) in df that should be included in the fit. Ideally, there is a raster companion to the column. These columns are not splined
                                  weight_column = NULL, #column in the data frame specifying weights
                                  bam =F, #should mgcv::bam be called rather than gam?
                                  spline_args = list(bs = 'ts', k = 3), #spline arguments to be applied to the columns specified by the covariates argument
                                  auto_model_select =T, #should the function override bam calls in low data situations (boosts probability of convergence)
                                  indicator = indicator, #outcome variable column
                                  indicator_family = 'binomial', #family of the outcome. Binomial and Gaussian are implemented.
                                  cores = slots) #number of compute cores available

#Fit a GBM/BRT model
gbm_results = fit_gbm_child_model(df = the_data,
                                  model_name = 'gbm',
                                  fold_id_col = 'fold_id',
                                  covariates = all_fixed_effects,
                                  weight_column = NULL,
                                  tc = 4, #tree complexity
                                  lr = 0.005, #learning rate
                                  bf = 0.75, #bag fraction
                                  indicator = indicator,
                                  indicator_family = indicator_family,
                                  cores = slots)
#Fit a random forest
rf_results = fit_rf_child_model(df = the_data,
                                model_name = 'rf',
                                fold_id_col = 'fold_id',
                                covariates = all_fixed_effects,
                                additional_terms = NULL,
                                ntree = 1000, #number of trees in the forest
                                indicator = indicator,
                                indicator_family = indicator_family,
                                cores = slots)
#fit penalized regressions
lasso_results = fit_glmnet_child_model(df = the_data,
                                       model_name = 'lasso',
                                       fold_id_col = 'fold_id',
                                       covariates = all_fixed_effects,
                                       additional_terms = NULL,
                                       weight_column = NULL,
                                       alpha = 1, #The elasticnet mixing parameter, with 0≤α≤ 1. 0 is ridge, 1 is lasso
                                       indicator = indicator,
                                       indicator_family = 'binomial',
                                       cores = slots)

ridge_results = fit_glmnet_child_model(df = the_data,
                                       model_name = 'ridge',
                                       fold_id_col = 'fold_id',
                                       covariates = all_fixed_effects,
                                       additional_terms = NULL,
                                       weight_column = NULL,
                                       alpha = 0, #The elasticnet mixing parameter, with 0≤α≤ 1. 0 is ridge, 1 is lasso
                                       indicator = indicator,
                                       indicator_family = 'binomial',
                                       cores = slots)

enet_results = fit_glmnet_child_model(df = the_data,
                                      model_name = 'enet',
                                      fold_id_col = 'fold_id',
                                      covariates = all_fixed_effects,
                                      additional_terms = NULL,
                                      weight_column = NULL,
                                      alpha = .5, #The elasticnet mixing parameter, with 0≤α≤ 1. 0 is ridge, 1 is lasso
                                      indicator = indicator,
                                      indicator_family = 'binomial',
                                      cores = slots)


#combine the children models
the_data = cbind(the_data, gam_results[[1]], gbm_results[[1]], rf_results[[1]], lasso_results[[1]],ridge_results[[1]],enet_results[[1]])

#fit gam stacker
stacked_results = gam_stacker(df = the_data, #the dataset in data table
                              model_names=c('gam','gbm','rf','lasso','ridge','enet'), #models to be stacked
                              weight_column = NULL, #prefixes of the models to be stacked
                              bam = F, #whether or not bam should be used
                              spline_args = list(bs = 'ts', k = 3),
                              indicator = indicator, #the indicator of analysis
                              indicator_family = indicator_family,
                              cores = slots)




#return the stacked rasters
stacked_rasters = make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                     period = NULL, #period of analysis, NULL returns all periods
                                     child_models = list(gam_results[[2]], gbm_results[[2]], rf_results[[2]], lasso_results[[2]],ridge_results[[2]],enet_results[[2]]), #model objects fitted to full data
                                     stacker_model = stacked_results[[2]],
                                     indicator_family = 'binomial',
                                     return_children = F)

the_data =cbind(the_data, stacked_results[[1]]) #may be unnecessary and overwritten later in code

#copy things back over to df
df = copy(the_data)

#alter the cov list and fixed effects
all_fixed_effects = 'stacked_results'

#########################################################################



## Save all inputs for MBG model into correct location on /share
save_mbg_input(indicator = indicator, 
               indicator_group = indicator_group,
               df = df,
               simple_raster = simple_raster,
               mesh_s = mesh_s,
               mesh_t = mesh_t, 
               cov_list = stacked_rasters)
load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/', run_date, '.RData'))

## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Generate MBG formula for INLA call
mbg_formula <- build_mbg_formula(fixed_effects = all_fixed_effects)

## Create SPDE INLA stack
input_data <- build_mbg_data_stack(df = df,
                                   fixed_effects = all_fixed_effects,
                                   mesh_s = mesh_s,
                                   mesh_t = mesh_t)
stacked_input <- input_data[[1]]
spde <- input_data[[2]]
cs_df <- input_data[[3]]

## Generate other inputs necessary
outcome=df[[indicator]] # N+_i - event obs in cluster
N=df$N                  # N_i - total obs in cluster

## Fit MBG model
model_fit <- fit_mbg(indicator_family = indicator_family,
                     stack.obs = stacked_input,
                     spde = spde,
                     cov = outcome,
                     N = N,
                     int_prior_mn = intercept_prior,
                     f_mbg = mbg_formula,
                     run_date = run_date,
                     keep_inla_files = keep_inla_files,
                     cores = slots)

predict_mbg <- predict_mbg(res_fit = model_fit,
                           cs_df = cs_df,
                           mesh_s = mesh_s,
                           mesh_t = mesh_t,
                           cov_list = stacked_rasters,
                           samples = samples,
                           simple_raster = simple_raster,
                           transform = transform)
mean_ras <- predict_mbg[[1]]
sd_ras <- predict_mbg[[2]]
cell_pred <- predict_mbg[[3]]

## Save MBG outputs in standard outputs folder structure
save_mbg_preds(config = config,
               time_stamp = time_stamp,
               run_date = run_date,
               mean_ras = mean_ras,
               sd_ras = sd_ras,
               res_fit = model_fit,
               cell_pred = cell_pred,
               df = df)

## ~~~~~~~~~~~~~~~~~~~~~~~~ Post-Estimation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get aggregated estimates for all admin0. Aggregate to level you rake to
#pop_wts_adm0 <- make_population_weights(admin_level   = 0,
#                                        simple_raster = simple_raster,
#                                        pop_raster    = pop_raster,
#                                        gaul_list     = gaul_list)

#cond_sim_raw_adm0 <- make_condSim(pop_wts_object = pop_wts_adm0,
#                                  cell_pred      = cell_pred,
#                                  summarize      = TRUE,
#                                  admin_level    = 0)

#ES - Not relevant for LF
#gbd <- load_gbd_data(gbd_type = "output",
#                     gbd_name = 322,
#                     gaul_list = gaul_list,
#                     measure_id = 6,
#                     age_group_id = 1,
#                     metric_id = 3)

## Get raking factors
#rf   <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
#                            rake_to     = gbd)

#raked_cell_pred <- rake_predictions(raking_factors = rf,
#                                    pop_wts_object = pop_wts_adm0,
#                                    cell_pred      = cell_pred)

## summarize raked predictions for each cell
mean_raster <- make_cell_pred_summary( draw_level_cell_pred = cell_pred,
                                       mask                 = simple_raster,
                                       return_as_raster     = TRUE,
                                       summary_stat         = 'mean')

#raked_mean_raster <- make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
#                                             mask                 = simple_raster,
#                                             return_as_raster     = TRUE,
#                                             summary_stat         = 'mean')

## save important post estimation things
save_post_est(mean_raster,'raster','mean_raster')
#save_post_est(raked_mean_raster,'raster','mean_raked_raster')
#condsims = mget(grep('cond_sim_raked_adm*', ls(), value = TRUE))
#save_post_est(condsims,'rdata','mean_raked_values_all_adm')
#save_post_est(rf,'csv','raking_factors')

## Clean up temporary INLA dirs
cleanup_inla_scratch(run_date)

## ~~~~~~~~~~~~~~~~~~~~~~~~ Shiny Diagnostics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~http://mbg-viz.duckdns.org:3456/~~~~~~~~~~~~~~~~~~~~~~~

# Plot Shiny stuff
shiny_data_and_preds(gaul_list = get_gaul_codes('africa'),
                     run_date = run_date,
                     indicator = indicator,
                     indicator_group = indicator_group,
                     pred_file = 'lf_pos_mean_raster.tif',
                     layer_name = 'lf_pos_mean_raster.') 

#ES - Not relevant for LF
# shiny_raked(gaul_list = get_gaul_codes('africa'),
#             run_date = run_date,
#             indicator = indicator,
#             indicator_group = indicator_group,
#             pred_file = 'has_lri_mean_raked_raster.tif',
#             layer_name = 'has_lri_mean_raked_raster.') 

shiny_cov_layers(fixed_effects = all_fixed_effects,
                 gaul_list = gaul_list,
                 run_date = run_date,
                 indicator = indicator,
                 indicator_group = indicator_group)

#shiny_raking_map()

##############################################################################
##############################################################################
##############################################################################

