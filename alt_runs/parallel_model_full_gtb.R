reg=as.character(commandArgs()[3])
age=as.numeric(commandArgs()[4])
run_date=as.character(commandArgs()[5])
test=as.character(commandArgs()[6])
holdout=as.character(commandArgs()[7])
indicator=as.character(commandArgs()[8])
indicator_group=as.character(commandArgs()[9])

pathaddin = paste0('_bin',age,'_',reg,'_',holdout)
message(pathaddin)
message(run_date)
message(test)

make_country_fe <- 1

## Load image saved from make_qsub in launch_lri_by_region
load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))


#load an image of the main environment
load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))
slots = as.numeric(slots)

#reload functions
setwd(repo)
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
package_lib <- paste0(root,'/temp/geospatial/packages') # Library for all MBG versioned packages. Ensures that none of this code is
#   dependent on the machine where the user runs the code.
.libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library(). Necessary for seeg libraries.
source('mbg_central/mbg_functions.R')                   # Functions to run MBG model.
source('mbg_central/prep_functions.R')                  # Functions to setup MBG run
source('mbg_central/covariate_functions.R')             # Functions to prep and transform 5*5 covariates
source('mbg_central/misc_functions.R')                  # Other computational MBG-related functions.
source('mbg_central/post_estimation_functions.R')
source('mbg_central/gbd_functions.R')
source('mbg_central/shiny_functions.R')
source('mbg_central/stacking_functions.R')
source('mbg_central/categorical_variable_functions.R')
source('mbg_central/seegMBG_transform_functions.R')
source('mbg_central/validation_functions.R') 
source('/share/code/geospatial/fever/stackin/stacking_functions.R')
source('/share/code/geospatial/fever/mbg/fever/function_tests/categorical_variable_functions.R') #functions to allow categorical variables
#source(paste0(repo,'fever/function_tests/update_predict_mbg.R')) #an improved prediction function
#call stacking functions

#functions come from the saved imaged
package_list <- c('rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr','foreign')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cores_to_use <- round(slots*.5)

## Load simple polygon template to model over
gaul_list <- get_gaul_codes(reg)
suppressMessages(suppressWarnings(simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                                                             buffer = 0.4)))
subset_shape     <- simple_polygon_list[[1]]
simple_polygon   <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
raster_list <- suppressMessages(suppressWarnings(build_simple_raster_pop(subset_shape)))
simple_raster <- raster_list[['simple_raster']]
pop_raster <- raster_list[['pop_raster']]

## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
if(holdout!=0) {
  df <- as.data.table(stratum_qt[[paste0('region__',reg)]])
  oos_df <- df[fold == holdout, ]
  i <- length(unique(oos_df$year))
  periods <- data.frame(group = rep(1:i,5),years = rep(sort(unique(oos_df$year)),5))
  oos_df$period <- match(oos_df$year, periods$years) # add these to df
  df <- df[fold != holdout, ]
}
if(holdout==0) {
  df <- load_input_data(indicator = indicator,
                        simple = simple_polygon,
                        removeyemen = TRUE,
                        pathaddin = pathaddin,
                        years = 'annual')
}

## Define modeling space (right now, just in years)
period_map <- make_period_map(modeling_periods = c(2000:2015)) 

#make covariates conditional
cov_layers = NULL
gbd_cov_layers = NULL
mbg_cov_layers = NULL

# Pull all covariate bricks/layers
if(nchar(fixed_effects)> 0){
  selected_fixed_effects <- strsplit(fixed_effects," ")
  selected_fixed_effects <- selected_fixed_effects[[1]][selected_fixed_effects[[1]] != "+"]
  selected_measures <- strsplit(fixed_effects_measures," ")
  selected_measures <- selected_measures[[1]][selected_measures[[1]] != "+"]
  cov_layers <- load_and_crop_covariates_annual(covs = selected_fixed_effects,                
                                                measures = selected_measures,          
                                                simple_polygon = simple_polygon,
                                                start_year  = 2000,
                                                end_year    = 2015,
                                                interval_mo = 12,
                                                agebin=1)
}
if(nchar(gbd_fixed_effects)>0){
  selected_gbd_fixed_effects <- strsplit(gbd_fixed_effects," ")
  selected_gbd_fixed_effects <- selected_gbd_fixed_effects[[1]][selected_gbd_fixed_effects[[1]] != "+"]
  # Create list of gaul codes for region + any countries in the data that aren't in region (buffer zone)
  gbd_gaul_list <- unique(c(gaul_convert(unique(df[, country])), gaul_list))
  # Use a layer from the geospatial cov_layers list as a template (raster of simple_polygon) to rasterize GBD cov spdfs.
  gbd_cov_layers <- suppressMessages(suppressWarnings(load_gbd_covariates(gbd_fixed_effects = selected_gbd_fixed_effects, 
                                                                          year_ids = c(2000:2015), 
                                                                          gaul_list = gbd_gaul_list,
                                                                          template = cov_layers[[1]][[1]])))
}
if(nchar(mbg_fixed_effects)>0){
  mbg_cov_layers <- suppressMessages(suppressWarnings(load_mbg_covariates(mbg_fixed_effects = mbg_fixed_effects, 
                                                                          simple_polygon = simple_polygon)))
}

## Combine all the covariate layers we want to use and combine our fixed effects formula
#all_cov_layers <- c(cov_layers, gbd_cov_layers, mbg_cov_layers)

#update all cov layers with an indicator variable on country
all_cov_layers <- c(cov_layers, gbd_cov_layers)
source('mbg_central/stacking_functions.R')

if(make_country_fe){
  fe_gaul_list <- unique(c(gaul_convert(unique(df[, country])), gaul_list))  
  fe_template = cov_layers[[1]][[1]]
  suppressMessages(suppressWarnings(simple_polygon_list <- load_simple_polygon(gaul_list = fe_gaul_list,
                                                                               buffer = 0.4,
                                                                               subset_only = TRUE)))
  fe_subset_shape     <- simple_polygon_list[[1]]
  gaul_code <- rasterize(fe_subset_shape, fe_template, field = 'GAUL_CODE')
  gaul_code = setNames(gaul_code,'gaul_code')
  gaul_code = create_categorical_raster(gaul_code)
  
  #update covlayers and add country fixed effects to the
  all_cov_layers = update_cov_layers(all_cov_layers, gaul_code)
}

#sort out other categorical variables
# if(nchar(categorical_vars) > 0){
#   #convert categorical vars to a multi object vector
#   cat_vars = unlist(strsplit(categorical_vars, " ", fixed = T))
#   cat_ras = lapply(cat_vars, function(x) create_categorical_raster(all_cov_layers[[x]]))
#   all_cov_layers = update_cov_layers(unlist(cat_ras))
#   all_cov_layers = all_cov_layers[!names(all_cov_layers) %in% cat_vars] #remove the pre blocked out ones
# }

#add latitude and longitude as covariates
# lat_ras = setNames(init(all_cov_layers[[1]][[1]], 'y'), 'lat_ras')
# long_ras = setNames(init(all_cov_layers[[1]][[1]], 'x'), 'long_ras')
# space_interact = setNames((lat_ras + abs(cellStats(lat_ras, stat = 'min'))) * (long_ras + abs(cellStats(long_ras, stat = 'min'))), 'space_interact')
# 
# #update all_cov_layers
# all_cov_layers = c(unlist(all_cov_layers), lat_ras = lat_ras, long_ras = long_ras,space_interact = space_interact)

#regenerate all fixed effects equation from the cov layers
all_fixed_effects = paste(names(all_cov_layers), collapse = " + ")

###STACK THE THINGS###
the_covs = format_covariates(all_fixed_effects)

#copy the dataset to avoid unintended namespace conflicts
the_data = copy(df)

#shuffle the data into six folds
n_stack_folds <- 5
the_data = the_data[sample(nrow(the_data)),]
the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=as.numeric(n_stack_folds),labels=FALSE)] #make folds

#extract covariates to the points and subset data where its missing covariate values
cs_covs = extract_covariates(the_data, all_cov_layers, return_only_results = T, centre_scale = T, period_var = 'year', period_map = period_map)

the_data = cbind(the_data, cs_covs[[1]])
covs_cs_df = cs_covs[[2]]

the_data = na.omit(the_data, c(indicator, 'N', the_covs)) #this will drop rows with NA covariate values

#Fit a gam model
gam = fit_gam_child_model(df = the_data, #data frame
                          model_name = 'gam', #identifier for the child model-- needs to be consistent through the steps
                          fold_id_col = 'fold_id',
                          covariates = all_fixed_effects, #rhs formula
                          additional_terms = 'year', #column(s) in df that should be included in the fit. Ideally, there is a raster companion to the column. These columns are not splined
                          weight_column = 'weight', #column in the data frame specifying weights
                          bam =F, #should mgcv::bam be called rather than gam?
                          spline_args = list(bs = 'ts', k = 3), #spline arguments to be applied to the columns specified by the covariates argument
                          auto_model_select =T, #should the function override bam calls in low data situations (boosts probability of convergence)
                          indicator = indicator, #outcome variable column
                          indicator_family = indicator_family, #family of the outcome. Binomial and Gaussian are implemented.
                          cores = 10) #number of compute cores available

#Fit a GBM/BRT model
gbm = fit_gbm_child_model(df = the_data,
                          model_name = 'gbm',
                          fold_id_col = 'fold_id',
                          covariates = all_fixed_effects,
                          weight_column = 'weight',
                          tc = 2, #tree complexity, change back to 4 for real runs
                          lr = 0.005, #learning rate
                          bf = 0.75, #bag fraction
                          indicator = indicator,
                          indicator_family = indicator_family,
                          cores = cores_to_use)

#fit some nets
#lasso
lasso = fit_glmnet_child_model(df = the_data,
                               model_name = 'lasso',
                               covariates =all_fixed_effects,
                               fold_id_col = 'fold_id',
                               additional_terms = NULL,
                               indicator_family = indicator_family,
                               indicator = indicator,
                               cores = cores_to_use,
                               alpha = 0,
                               weight_column = 'weight')


#ridge
ridge = fit_glmnet_child_model(df = the_data,
                               model_name = 'ridge',
                               covariates = all_fixed_effects,
                               fold_id_col = 'fold_id',
                               additional_terms = NULL,
                               indicator_family = indicator_family,
                               indicator = indicator,
                               cores = cores_to_use,
                               alpha = 1,
                               weight_column = 'weight')

#enet
enet = fit_glmnet_child_model(df = the_data,
                              model_name = 'enet',
                              covariates = all_fixed_effects,
                              fold_id_col = 'fold_id',
                              additional_terms = NULL,
                              indicator_family = indicator_family,
                              indicator = indicator,
                              cores = cores_to_use,
                              alpha = .5,
                              weight_column = 'weight')


#combine the children models
the_data = cbind(the_data, do.call(cbind, lapply(list(gam,gbm,lasso,ridge,enet), function(x) x[[1]])))
child_model_names = c('gam', 'gbm','lasso','ridge','enet')
child_model_objs = setNames(lapply(list(gam,gbm,lasso,ridge,enet), function(x) x[[2]]), child_model_names)

#fit stacker
stacked_results = gam_stacker(the_data, #the dataset in data table format
                              model_names= child_model_names, #prefixes of the models to be stacked
                              indicator = indicator, #the indicator of analysis
                              indicator_family = indicator_family) #indicator family (e.g. binomial)

#return the stacked rasters
stacked_rasters = make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                     period = min(period_map[, period_id]):max(period_map[, period_id]), #period of analysis, NULL returns all periods
                                     child_models = child_model_objs, #model objects fitted to full data
                                     stacker_model = stacked_results[[2]],
                                     indicator_family = indicator_family,
                                     return_children = T,
                                     centre_scale_df = covs_cs_df)

# if(as.logical(gpr_stack)){
#   all_fixed_effects = paste0(names(stacked_rasters[2:length(stacked_rasters)]), collapse = ' + ')
# } else{
#   all_fixed_effects = 'stacked_results'
# }
#all_fixed_effects = 'stacked_results'
#all_fixed_effects = paste('gbm','lasso','ridge','enet', names(gaul_code), sep=" + ")
#eval(parse(text=stacked_fixed_effects))

#all_fixed_effects = paste(stacked_fixed_effects, paste(names(gaul_code), collapse = " + "), sep=" + ")
all_fixed_effects <- stacked_fixed_effects

#copy things back over to df
#df = copy(the_data)

#remove the covariate columns so that there are no name conflicts when they get added back in
#df = df[,paste0(the_covs) := rep(NULL, length(the_covs))]

df <- load_input_data(indicator = indicator,
                      simple = simple_polygon,
                      removeyemen = TRUE,
                      years = 'annual')

## Build spatial mesh over modeling area
mesh_s <- build_space_mesh(d = df,
                           simple = simple_polygon,
                           max_edge = mesh_s_max_edge,
                           mesh_offset = mesh_s_offset)

## Build temporal mesh (standard for now)
mesh_t <- build_time_mesh(periods=eval(parse(text=mesh_t_knots)))

#create a full raster list to carry though to the shiny/next steps
full_raster_list = c(unlist(stacked_rasters),unlist(all_cov_layers))
child_mod_ras = full_raster_list[child_model_names]

## Save all inputs for MBG model into correct location on /share
save_mbg_input(indicator = indicator, 
               indicator_group = indicator_group,
               df = df,
               simple_raster = simple_raster,
               mesh_s = mesh_s,
               mesh_t = mesh_t, 
               cov_list = full_raster_list,
               pathaddin = pathaddin,
               run_date = run_date,
               child_model_names = child_model_names,
               all_fixed_effects = all_fixed_effects,
               period_map = period_map) #specify by region

#reload data an prepare for MBG
load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/', run_date, pathaddin, '.RData'))

## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# if(as.logical(gpr_stack)){
#   #for stacking, overwrite the columns matching the model_names so that we can trick inla into being our stacker
#   df = df[,paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]
# }

## Generate MBG formula for INLA call
mbg_formula <- build_mbg_formula(fixed_effects = all_fixed_effects)

## Create SPDE INLA stack
input_data <- build_mbg_data_stack(df = df,
                                   fixed_effects = all_fixed_effects,
                                   mesh_s = mesh_s,
                                   mesh_t = mesh_t)

#combine all the inputs
stacked_input <- input_data[[1]]
spde <- input_data[[2]]
cs_df <- input_data[[3]]

#binomial data
## Generate other inputs necessary
outcome=df[[indicator]] # N+_i - event obs in cluster
N=df$N                  # N_i - total obs in cluster
weights=df$weight

#catch in case there is no weight column
if(is.null(weights)){
  weights = rep(1,nrow(df))
}


## Fit MBG model
source('mbg_central/mbg_functions.R')
model_fit <- fit_mbg(indicator_family = indicator_family,
                     stack.obs = stacked_input,
                     spde = spde,
                     cov = outcome,
                     N = N,
                     int_prior_mn = intercept_prior,
                     f_mbg = mbg_formula,
                     run_date = run_date,
                     keep_inla_files = keep_inla_files,
                     cores = cores_to_use,
                     wgts = weights)

#predict some MBG
source('mbg_central/mbg_functions.R')
predict_mbg <- predict_mbg(res_fit = model_fit,
                           cs_df = cs_df,
                           mesh_s = mesh_s,
                           mesh_t = mesh_t,
                           cov_list = full_raster_list,
                           samples = samples,
                           simple_raster = simple_raster,
                           transform = transform)

## Save MBG outputs in standard outputs folder structure
mean_ras <- predict_mbg[[1]]
sd_ras <- predict_mbg[[2]]
cell_pred <- predict_mbg[[3]]

##################################################################################################
############################### SAVE AND CROSS-VAL ##### #########################################
##################################################################################################
## Save MBG outputs in standard outputs folder structure
save_mbg_preds(config = config,
               time_stamp = time_stamp,
               run_date = run_date,
               mean_ras = mean_ras,
               sd_ras = sd_ras,
               res_fit = model_fit,
               cell_pred = cell_pred,
               df = df,
               pathaddin = pathaddin)

if(holdout != 0){
  message('Doing out of sample predictive validity over quadtree aggregates')
  
  aggval <- aggregate_validation_dev(holdoutlist         = stratum_qt,
                                     cell_draws_filename = paste0('%s_cell_draws_eb_bin%i_%s_%i.RData'),
                                     years               = c(2000,2005,2010,2015),
                                     indicator_group     = indicator_group,
                                     indicator           = indicator,
                                     run_date            = run_date,
                                     reg                 = reg,
                                     holdout             = holdout,
                                     iter                = 1)
  
  # get validation metrics
  draws <- fit_stats_ho_id(draws       = aggval,
                           draw_prefix = 'phat_',
                           observed    = 'p')
  draws$covered = NULL # not really coverage!
  
  save_post_est(draws,'csv', paste0('validation',pathaddin))
  save_post_est(aggval,'csv',paste0('aggval',    pathaddin))
  
}



# ## Calculate fit statistics unless holdout==0 which indicates the full data run.
# if(holdout!=0) {
#   job_fit_stats <- fit_stats(is_data = df,
#                              oos_data = oos_df,
#                              indicator = indicator,
#                              indicator_group = indicator_group,
#                              run_date = run_date,
#                              pathaddin = pathaddin)
# } else {
#   #build summaries
#   raster_list    <- build_simple_raster_pop(subset_shape)
#   simple_raster  <- raster_list[['simple_raster']]
#   pop_raster     <- raster_list[['pop_raster']]
# 
# 
#   ## summarize raked predictions for each cell
#   mean_raster <- make_cell_pred_summary( draw_level_cell_pred = cell_pred,
#                                            mask                 = simple_raster,
#                                            return_as_raster     = TRUE,
#                                            summary_stat         = 'mean')
#   #fix naming of the mean raster
#   assign(paste0(reg,'_mean_raster'), mean_raster)
#   rm(mean_raster)
# 
#   #estimate for national and subnationl units
#   for(ad in 0){
#     pop_wts <- make_population_weights(admin_level   = ad,
#                                        simple_raster = simple_raster,
#                                        pop_raster    = pop_raster,
#                                        gaul_list     = get_gaul_codes(reg))
# 
#     condsim <-    make_condSim(pop_wts_object = pop_wts,
#                                cell_pred      = cell_pred,
#                                gaul_list      = get_gaul_codes(reg),
#                                admin_level    = ad,
#                                summarize      = TRUE)
# 
#     condsim <- data.table(cbind(split_geo_names(as.matrix(condsim)),mean=unname(condsim) ))
#     assign(paste0('cond_sim_adm',ad,'_',reg),condsim)
#     rm(condsim)
# 
#   }
# 
# 
#   #save a summary object the includes the mean raster and conditional sims
#   print(paste0('Saving to: ', '/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/', indicator, '_result_summaries_', reg, '.RData'))
# 
#   save(list = c(paste0(reg,'_mean_raster'), #mean raster for the region
#        paste0('cond_sim_adm',0:2,'_',reg), 'full_raster_list'), #the conditional simulations for admin 0, 1 and 2
#        file = paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/', indicator, '_result_summaries_', reg, '.RData')
#     )
# }
