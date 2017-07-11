## Set repo location and indicator group
repo <- '/share/code/geospatial/stearns7/mbg/'
indicator_group <- 'lf'
indicator <- 'lf_pos'
Regions = c('wssa')  #created new region - essa_lf - , added zim to essa; 'essa','name','wssa','cssa'; trying new region of all lf endemic africa countries


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
gaul_list <- get_gaul_codes(Regions) # TO DO ------ CHANGE HERE -----------------------------------------------
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                           buffer = 0.4)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]

## Load input data and make sure it is cropped to modeling area
df_list <- load_input_data(indicator = indicator,
                           simple = simple_polygon,
                           removeyemen = TRUE,
                           years= 'annual',
                           update_run_date = TRUE)

df <- df_list[[1]]
run_date <- df_list[[2]]

## Add GAUL_CODE and region to df given the lat/longs. We need this to stratify in holdout function.
df <- add_gauls_regions(df = df,
                        simple_raster = simple_raster)
## Run function to create holdouts (returns list of data.frames with an additional "folds" column)
table(df$region, df$year)

## Set strata as character vector of each strata (in my case, just stratifying by region whereas U5M stratifies by region/age)
strata <- Regions

year_list <- eval(parse(text=year_list))

## Define modeling space (right now, just in years)
period_map <- make_period_map(modeling_periods = year_list) 

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
                                                start_year  = min(year_list),
                                                end_year    = max(year_list),
                                                interval_mo = 12,
                                                agebin=1)
}

#update all cov layers with an indicator variable on country
all_cov_layers <- c(cov_layers, gbd_cov_layers)
source('mbg_central/stacking_functions.R')

if(use_child_country_fes == TRUE | use_inla_country_fes == TRUE){
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

#regenerate all fixed effects equation from the cov layers
all_fixed_effects = paste(names(all_cov_layers), collapse = " + ")
if(use_child_country_fes==FALSE) all_fixed_effects = paste(names(all_cov_layers)[!grepl("gaul_code_", names(all_cov_layers))], collapse = " + ")

# Figure out which models we're going to use
child_model_names <- stacked_fixed_effects %>% 
  gsub(" ", "", .) %>%
  strsplit(., "+", fixed=T) %>%
  unlist

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
if ('gam' %in% child_model_names) {

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

}

#Fit a GBM/BRT model
if ('gbm' %in% child_model_names) {

  gbm = fit_gbm_child_model(df = the_data,
                            model_name = 'gbm',
                            fold_id_col = 'fold_id',
                            covariates = all_fixed_effects,
                            weight_column = 'weight',
                            tc = 3, #tree complexity, change back to 4 for real runs
                            lr = 0.005, #learning rate
                            bf = 0.75, #bag fraction
                            indicator = indicator,
                            indicator_family = indicator_family,
                            cores = 10)

}

#fit some nets
#lasso
if ('lasso' %in% child_model_names) {

  lasso = fit_glmnet_child_model(df = the_data,
                                 model_name = 'lasso',
                                 covariates =all_fixed_effects,
                                 fold_id_col = 'fold_id',
                                 additional_terms = NULL,
                                 indicator_family = indicator_family,
                                 indicator = indicator,
                                 cores = 10,
                                 alpha = 0,
                                 weight_column = 'weight')

}

#ridge
if ('ridge' %in% child_model_names) {

  ridge = fit_glmnet_child_model(df = the_data,
                                 model_name = 'ridge',
                                 covariates = all_fixed_effects,
                                 fold_id_col = 'fold_id',
                                 additional_terms = NULL,
                                 indicator_family = indicator_family,
                                 indicator = indicator,
                                 cores = 10,
                                 alpha = 1,
                                 weight_column = 'weight')

}

#enet
if ('enet' %in% child_model_names) {

  enet = fit_glmnet_child_model(df = the_data,
                                model_name = 'enet',
                                covariates = all_fixed_effects,
                                fold_id_col = 'fold_id',
                                additional_terms = NULL,
                                indicator_family = indicator_family,
                                indicator = indicator,
                                cores = 10,
                                alpha = .5,
                                weight_column = 'weight')
                 
}

#combine the children models
the_data = cbind(the_data, do.call(cbind, lapply(lapply(child_model_names, 'get'), function(x) x[[1]])))
child_model_objs = setNames(lapply(lapply(child_model_names, 'get'), function(x) x[[2]]), child_model_names)

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
if(use_inla_country_fes) all_fixed_effects = paste(stacked_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep=" + ")


#copy things back over to df
df = copy(the_data)

#remove the covariate columns so that there are no name conflicts when they get added back in
df = df[,paste0(the_covs) := rep(NULL, length(the_covs))]

## Double-check that gaul codes get dropped before extracting in save_mbg_input()
df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]

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



#for stacking, overwrite the columns matching the model_names so that we can trick inla into being our stacker
df = df[,paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]

## Generate MBG formula for INLA call
if(use_inla_nugget == TRUE) df$IID.ID <- 1:nrow(df)
mbg_formula <- build_mbg_formula(fixed_effects = all_fixed_effects,
                                 add_nugget = use_inla_nugget)



## Create SPDE INLA stack
input_data <- build_mbg_data_stack(df = df,
                                   fixed_effects = all_fixed_effects,
                                   mesh_s = mesh_s,
                                   mesh_t = mesh_t)

#combine all the inputs
stacked_input <- input_data[[1]]
spde <- input_data[[2]]
cs_df <- input_data[[3]]

# Add nugget if toggled
if(use_inla_nugget==TRUE){
  stacked_input$effects$data$IID.ID <- 1:stacked_input$effects$nrow
  stacked_input$effects$ncol <- c(stacked_input$effects$ncol, 1)
  names(stacked_input$effects$ncol)[length(names(stacked_input$effects$ncol))] <- "IID.ID"
  stacked_input$effects$names[["IID.ID"]] <- "IID.ID"
}

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



tic("Getting covariate weights")

x <- (fixed_effects %>% 
        gsub(" ", "", .) %>%
        strsplit(., "+", fixed=T))%>%
  unlist
cov_wts <- get.cov.wts(stacked.model.objects = child_model_objs,
                       inla.fit = model_fit,
                       raw.covs = x,
                       fit.data = the_data)
# Write to csv
output_dir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, "/")
output_file_covs <- paste0(output_dir, "cov.wts_", indicator, "_", run_date, ".csv")
write.csv(cov_wts,
          file      = output_file_covs,
          row.names = F)

toc(log=T)
#may need to change 'fit.data' bc same name for data that went into stacker as data cbinded after


## testing things here
smo = child_model_objs
rc = x
fit = model_fit
X = the_data


zed <- list(gbm, lasso, ridge, enet)


rownames(imp.mat) <- c(child_model_names, "INLA COMBINED")

x = seq(1:10)
y = seq(11:20)
z <- seq(21:30)

xx <- matrix(c(x,y), nrow = 10, ncol=10)
yy <- matrix(c(y,z), nrow = 10, ncol=10)
zz <- matrix(c(z,x), nrow = 10, ncol=10)

xyz <- list(xx,yy,zz)


t <- lapply(seq_along(xyz), function(x) names(xyz[x]))
t <- lapply(seq_along(xyz), function(i) names(xyz)[i])

t <- sapply(xyz, paste0, collapse="")

smo_names <- for (i in length(smo)) {
  name = smo[[i]][2]
  x[paste0(x,"", name)]
}
