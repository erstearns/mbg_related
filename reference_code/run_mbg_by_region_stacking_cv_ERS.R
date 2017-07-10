
## GET SET UP FROM THE QSUB CALL

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

## Load image saved from make_qsub in launch_lri_by_region
load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))

## Load libraries and miscellaneous MBG project functions.
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
  source('mbg_central/seegMBG_transform_functions.R')     # Using Roy's edit for now that can take temporally varying covariates,
#   TODO: will need to send pull request to seegMBG of github
  package_list <- c('rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr')
  for(package in package_list) {
    library(package, lib.loc = package_lib, character.only=TRUE)
  }


## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
                          pathaddin = pathaddin)
  }


##################################################################################################
############################### STACKING #########################################################
##################################################################################################

#make covariates conditional
  cov_layers = NULL
  gbd_cov_layers = NULL
  mbg_cov_layers = NULL

# Pull all covariate bricks/layers
  if(nchar(fixed_effects)> 0){
    cov_layers <- suppressMessages(suppressWarnings(load_and_crop_covariates(fixed_effects = fixed_effects, 
                                                                             simple_raster = simple_raster)))
  }
  if(nchar(gbd_fixed_effects)>0){
    gbd_cov_layers <- suppressMessages(suppressWarnings(load_gbd_covariates(gbd_fixed_effects, c(2000,2005,2010,2015), gaul_list)))
  }
  if(nchar(mbg_fixed_effects)>0){
    mbg_cov_layers <- suppressMessages(suppressWarnings(load_mbg_covariates(mbg_fixed_effects, simple_raster)))
  }


## Combine all the covariate layers we want to use and combine our fixed effects formula
  all_cov_layers <- c(cov_layers, gbd_cov_layers, mbg_cov_layers)
  effects = c(fixed_effects, gbd_fixed_effects, mbg_fixed_effects)
  all_fixed_effects <- paste(effects[!is.na(effects)&effects!=""], collapse=" + ")

#update all cov layers with an indicator variable on country
  gaul_code = setNames(simple_raster,'gaul_code')
  gaul_code = create_categorical_raster(gaul_code)

#update covlayers and fixed effects
  all_cov_layers = update_cov_layers(all_cov_layers, gaul_code)
  cat_effects = paste(names(gaul_code), collapse = " + ")  

###STACK THE THINGS###

  if(length(all_fixed_effects)==1){
    #split back into parts
    the_covs = unlist(tstrsplit(all_fixed_effects,"\\+"))
  }

#copy the dataset to avoid unintended namespace conflicts
  the_data = copy(df)

#shuffle the data into five folds
  the_data = the_data[sample(nrow(the_data)),]
  the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=5,labels=FALSE)] #make five folds

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

#fit some nets
#lasso
  lasso = fit_glmnet_child_model(df = the_data,
                                 model_name = 'lasso',
                                 covariates = paste0(all_fixed_effects, ' + ', cat_effects),
                                 fold_id_col = 'fold_id',
                                 additional_terms = NULL,
                                 indicator_family = 'binomial',
                                 indicator = indicator, 
                                 cores = 20,
                                 alpha = 0)


#ridge
  ridge = fit_glmnet_child_model(df = the_data,
                                 model_name = 'ridge',
                                 covariates = paste0(all_fixed_effects, ' + ', cat_effects),
                                 fold_id_col = 'fold_id',
                                 additional_terms = NULL,
                                 indicator_family = 'binomial',
                                 indicator = indicator, 
                                 cores = 20,
                                 alpha = 1)

#enet
  enet = fit_glmnet_child_model(df = the_data,
                                model_name = 'enet',
                                covariates = paste0(all_fixed_effects, ' + ', cat_effects),
                                fold_id_col = 'fold_id',
                                additional_terms = NULL,
                                indicator_family = 'binomial',
                                indicator = indicator, 
                                cores = 1,
                                alpha = .5)

#combine the children models
  the_data = cbind(the_data, gam_results[[1]], gbm_results[[1]], rf_results[[1]], lasso[[1]], ridge[[1]], enet[[1]])

#fit stacker
  stacked_results = gam_stacker(the_data, #the dataset in data table format
                                model_names=eval(parse(text=child_stacker_model_names)), #prefixes of the models to be stacked
                                indicator = indicator, #the indicator of analysis
                                indicator_family = indicator_family) #indicator family (e.g. binomial)

#return the stacked rasters
  stacked_rasters = make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks 
                                       period = NULL, #period of analysis, NULL returns all periods
                                       child_models = eval(parse(text=child_stacker_models)), #model objects fitted to full data
                                       stacker_model = stacked_results[[2]],
                                       indicator_family = 'binomial',
                                       return_children = T)

  the_data =cbind(the_data, stacked_results[[1]])

#copy things back over to df
  df = copy(the_data)
  

##################################################################################################
############################### FIT/PREDICT WITH STACKER #########################################
##################################################################################################

  # Combine all fixed effects and cov layers
  effects = c(fixed_effects, gbd_fixed_effects, mbg_fixed_effects, stacked_fixed_effects)
  all_fixed_effects <- paste(effects[!is.na(effects)&effects!=""], collapse=" + ")
  all_cov_layers[[stacked_fixed_effects]] <- stacked_rasters[[stacked_fixed_effects]]
  final_cov_layers <- all_cov_layers
  
  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d = df,
                             simple = simple_polygon,
                             max_edge = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset)
  
  ## Build temporal mesh (standard for now)
  mesh_t <- build_time_mesh()
  
  ## Save all inputs for MBG model into correct location on /share
  save_mbg_input(indicator = indicator, 
                 indicator_group = indicator_group,
                 df = df,
                 simple_raster = simple_raster,
                 mesh_s = mesh_s,
                 mesh_t = mesh_t, 
                 cov_list = final_cov_layers,
                 pathaddin = stacker_pathaddin) #specify by region
  
  #reload data an prepare for MBG
  load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/', run_date, stacker_pathaddin, '.RData'))
  
  
  
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
  weights=df$weight
  
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
                       cores = slots,
                       wgts = weights)
  
  source('mbg_central/mbg_functions.R')
  predict_mbg <- predict_mbg(res_fit = model_fit,
                             cs_df = cs_df,
                             mesh_s = mesh_s,
                             mesh_t = mesh_t,
                             cov_list = final_cov_layers,
                             samples = samples,
                             simple_raster = simple_raster,
                             transform = transform)
  
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
                 pathaddin = stacker_pathaddin)
  
  ## Calculate fit statistics unless holdout==0 which indicates the full data run.
  if(holdout!=0) {
    job_fit_stats <- fit_stats(is_data = df,
                               oos_data = oos_df,
                               indicator = indicator,
                               indicator_group = indicator_group,
                               run_date = run_date,
                               pathaddin = stacker_pathaddin)
  }
  

## Remove temporary image
file.remove(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))


