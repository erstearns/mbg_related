###STACK THE THINGS###

if(length(all_fixed_effects)==1){
  #split back into parts
  the_covs = unlist(tstrsplit(all_fixed_effects,"\\+"))
}

#copy the dataset to avoid unintended namespace conflicts
the_data = copy(df)

#shuffle the data into six folds
the_data = the_data[sample(nrow(the_data)),] #Shuffles data
the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=10,labels=FALSE)] #make 10 folds - divides shuffled data into 11 - could bring back down to 6 once over 2000 data pts for AFR; should be dynamically coded - specify in config file

#extract covariates to the points and subset data where its missing covariate values
the_data = extract_covariates(the_data, all_cov_layers) #extracts covariate values at each data point; reconcile_timevarying default = true; integrating lags may be necessary to use historical pts
the_covs = trimws(the_covs)
the_data = na.omit(the_data, c(indicator, 'N', the_covs)) #this will drop rows with NA covariate values

#Fit a gam model
gam_results = fit_gam_child_model(df = the_data, #data frame
                                  model_name = 'gam', #identifier for the child model-- needs to be consistent through the steps; also pre-fix for some results returned
                                  fold_id_col = 'fold_id', #separate from fold id nick and aaron
                                  covariates = all_fixed_effects, #rhs formula
                                  additional_terms = 'year', #column(s) in df that should be included in the fit. Ideally, there is a raster companion to the column. These columns are not splined; will automatically add any binary vars; if need to specify others, talk to daniel
                                  weight_column = NULL, #column in the data frame specifying weights; character string ref that references a column in the datatable to be specified with weight_column = "weight' and put weight col in input dataset
                                  bam =F, #should mgcv::bam be called rather than gam? separate more scalable impl. than GAM
                                  spline_args = list(bs = 'ts', k = 3), #spline arguments to be applied to the columns specified by the covariates argument; all covs splined in same way
                                  auto_model_select =T, #should the function override bam calls in low data situations (boosts probability of convergence); irrel w/ bam = f
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
                                       alpha = 1, #The elasticnet mixing parameter, with 0≤α≤ 1. 0 is ridge, 1 is lasso; penalizes based on model type
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

#fit gam stacker  #GAM model being fit on cross-validated data, then makes predictions using results from full child models
stacked_results = gam_stacker(df = the_data, #the dataset in data table
                              model_names=c('gam','gbm','rf','lasso','ridge','enet'), #models to be stacked
                              weight_column = NULL, #prefixes of the models to be stacked
                              bam = F, #whether or not bam should be used
                              spline_args = list(bs = 'ts', k = 3),
                              indicator = indicator, #the indicator of analysis
                              indicator_family = indicator_family,
                              cores = slots)




#return the stacked rasters; if elastic net throws an error, cov layers may not line up; 
stacked_rasters = make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                     period = NULL, #period of analysis, NULL returns all periods
                                     child_models = list(gam_results[[2]], gbm_results[[2]], rf_results[[2]], lasso_results[[2]],ridge_results[[2]],enet_results[[2]]), #model objects fitted to full data
                                     stacker_model = stacked_results[[2]],
                                     indicator_family = 'binomial',
                                     return_children = F) #if set to T, will return rasters for each child model

the_data =cbind(the_data, stacked_results[[1]])

#copy things back over to df
df = copy(the_data)

#alter the cov list and fixed effects
all_fixed_effects = 'stacked_results'

#########################################################################