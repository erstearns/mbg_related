## Load and crop covariates to the modeled area
load_and_crop_covariates <- function(fixed_effects, simple_raster, agebin=1) {
  
  # Hard-code central directories
  root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
  central_cov_dir <- paste0(root,'/WORK/11_geospatial/01_covariates/09_MBG_covariates/') # central folder per Lucas
  central_tv_covs <- c('evi','lights_new','LST_day','total_pop','rates','malaria','fertility')
  central_ntv_covs <- c('irrigation', 'access')
  central_ntv_ml_covs <- c('africa_lf1') #ES - not temporally variant, multi_layer raster: data held in more than 1 layer
  
  # Load all temporally-varying covariates
  evi             <- brick(paste0(central_cov_dir, 'EVI_stack.tif'))
  lights_new      <- brick(paste0(central_cov_dir, 'NTL_stack.tif'))
  LST_day         <- brick(paste0(central_cov_dir, 'LST_day_stack.tif'))
  total_pop       <- brick(paste0(central_cov_dir, 'WorldPop_allStages_stack.tif'))
  # Load all temporally-nonvarying covariates
  access          <- brick(paste0(central_cov_dir, 'synoptic_stack.tif'))$synoptic_stack.1
  irrigation      <- brick(paste0(central_cov_dir, 'synoptic_stack.tif'))$synoptic_stack.2
  
  #ES - Load all temporally non-varying, multi-layer raster covariates 
  africa_lf1      <- brick(paste0(central_cov_dir, 'africa_lf1.tif')) 
  
  
  # some u5m additions
  malaria         <- brick(paste0(central_cov_dir, 'malaria_infant_death_rate_stack.tif'))
  values(malaria)=log(as.matrix(malaria)+.01)
  fertility         <- brick(paste0(central_cov_dir, 'fertility_stack.tif'))
  u5m_dir <-paste0(root,'/temp/geospatial/U5M_africa/')
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_m0.Rdata'))
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_m1_11.Rdata'))
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_2q1.Rdata'))
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_2q3.Rdata'))
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_5q0.Rdata'))
  rates = get(paste0('rates_',agebin))
  names(rates)=paste0('rates.',1:4)
  
  # Add names to layers
  names(access) <- "access"
  names(irrigation) <- "irrigation"
  for(c in central_tv_covs){
    tmp=get(c)
    names(tmp)=rep(paste0(c,'.',1:4))
    assign(c,tmp)
  }
  
  # Construct list of covariates to GAM and use in model from fixed_effects parameter equation.
  selected_covs <- strsplit(fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]
  num_covs <- length(selected_covs)
  lcovs <- list()
  for(i in 1:num_covs) {
    this_cov <- selected_covs[i]
    if(this_cov %in% central_ntv_covs) { # Add if it is from the temporally non-varying list.
      lcovs[[i]] <- get(this_cov)
    }
    if(this_cov %in% central_tv_covs) { # Add if it is from the temporally varying list.
      lcovs[[i]] <- get(this_cov)
    }
    if(this_cov %in% central_ntv_ml_covs) { #ES - Add if it is from the temporally non-varying, multi-layer list.
      lcovs[[i]] <- get(this_cov) 
    }
    names(lcovs)[i] <- this_cov
  }
  
  # Make sure covariate layers line up with raster we are modeling over
  for(l in 1:length(lcovs)) {
    #message(names(lcovs)[l])
    lcovs[[l]]  <- crop(lcovs[[l]], extent(simple_raster))
    lcovs[[l]]  <- setExtent(lcovs[[l]], simple_raster)
    lcovs[[l]]  <- mask(lcovs[[l]], simple_raster)
  }
  
  return(lcovs)
  
}

# save a csv of contributions of different variables in brt
save_brt_contributions <- function(brt_mod_obj = trans_covs[[1]],
                                   a=age,r=reg,returnx=FALSE){
  x=data.table(rbind(cbind(brt_mod_obj[[1]]$contributions,year=2000,age=a,reg=r),
                     cbind(brt_mod_obj[[2]]$contributions,year=2005,age=a,reg=r),
                     cbind(brt_mod_obj[[3]]$contributions,year=2010,age=a,reg=r),
                     cbind(brt_mod_obj[[4]]$contributions,year=2015,age=a,reg=r)))
  x$var=gsub(pattern='\\.[0-9]',replacement='',x=x$var) # remove .# in varnames
  fwrite(x,paste0('/share/geospatial/mbg/', indicator_group, '/',
                  indicator, '/output/', run_date,'/',
                  'brt_contribs_bin',a,'_',r,'.csv'))
  if(returnx) return(x)
}



## Function to take input data and covariate layers and fit GAM models
#   Arguments:
#     df = data.table with variables "latitude", "longitude", "N", and specified indicator
#     lcovs = list of raster layers or bricks (for temporally varying). Currently assumes
#             four layers named like U5M periods.
#   Returns: Two-item list (bricks of time-varying and non-time-varying covariates)
#   Notes: This function currently won't work if you have anything other than *four* distinct periods in your data and covariates.
#           This will be made more flexible.
gam_covs <- function(df, indicator_family, lcovs) {
  coords <- df[, c('longitude', 'latitude'), with=FALSE]
  coords$lat=as.numeric(coords$latitude)
  coords$long=as.numeric(coords$longitude)
  coords <- coords[, c('long', 'lat'), with=FALSE]
  
  if(indicator_family=="binomial") response <- cbind(died = df[, get(indicator)], lived = df[, N] - df[, get(indicator)])
  if(indicator_family=="gaussian") response <- cbind(outcome = df[, get(indicator)])
  
  extra_data <- data.frame(year = df$year)
  
  # fit gam
  # This should take a few minutes
  system.time(trans <- gamTrans(coords=coords,
                                response=response,
                                covs=lcovs,
                                family = indicator_family,
                                extra_terms = ~ year,
                                extra_data = extra_data,
                                bam = TRUE,
                                predict = TRUE,
                                condition = NULL,
                                condition_covs = NULL,
                                s_args = list(bs = 'ts', k = 3),
                                samfrac = 0.1,
                                use.chol = TRUE))
  
  ### IF trans$trans IS RETURNED AS A LIST, IT IS SPLIT INTO TEMPORAL AND NON-TEMPORAL COVARIATES
  if(class(trans$trans)=='list') {
    temporal=TRUE
  } else {
    temporal=FALSE
  }
  
  message("CLAMP AND SAVE")
  if(!temporal){ ######################################################## - check into; 'clamp'?
    # THEY ALL PASS
    # use chi-squared stats to determine covariate usefulness
    keep <- which(summary(trans$model)$chi.sq > 0.1)
    trans_ras <- trans$trans[[keep]]
    
    # clamp covariates
    # find most extreme vaaues of transofmred covariates that were observed
    vals <- extract(trans_ras, coords[idx_fit, ])
    sry <- apply(vals, 2, range, na.rm = TRUE)
    
    # clamp the covariates to these values
    for (i in 1:nlayers(trans_ras)) {
      range <- sry[, colnames(sry) == names(trans_ras)[i]]
      trans_ras[[i]][trans_ras[[i]] < range[1]] <- range[1]
      trans_ras[[i]][trans_ras[[i]] > range[2]] <- range[2]
    }
    
    return(trans_ras)
    
  }
  
  # temporally varying covariates are present, save them all separately
  # non varying ones will be save in the same covs_transformed location as before
  if(temporal){
    
    # first clamp and save non temporally varying
    message('time invariant covariates')
    trans_ras=trans$trans$nT_vals_trans
    # clamp covariates
    # find most extreme vaaues of transofmred covariates that were observed
    vals <- extract(trans_ras, coords)
    sry <- apply(vals, 2, range, na.rm = TRUE)
    
    
    # clamp the covariates to these values
    for (i in 1:nlayers(trans_ras)) {
      range <- sry[, colnames(sry) == names(trans_ras)[i]]
      trans_ras[[i]][trans_ras[[i]] < range[1]] <- range[1]
      trans_ras[[i]][trans_ras[[i]] > range[2]] <- range[2]
    }
    
    # If you only specify one non-varying term, it simply gets name "layer" in the GAM function. Rather than fixing it in there
    # I'm just going to check if that's the case and rename it here.
    all_rasters <- list()
    if(length(names(trans_ras))==1) {
      for(cov in c('access','irrigation')) {
        if(cov %in% names(lcovs)) {
          names(trans_ras) <- cov
          all_rasters[[cov]] <- trans_ras
        }
      }
    }
    
    # Now, this same process for the individual temporally varying covariates
    count <- length(all_rasters) + 1
    for(n in names(trans$trans$T_vals_trans)){
      message(n)
      message(count)
      
      trans_ras=trans$trans$T_vals_trans[[n]]
      
      # clamp covariates
      # find most extreme vaaues of transformed covariates that were observed
      vals <- extract(trans_ras, coords)
      sry <- apply(vals, 2, range, na.rm = TRUE)
      
      
      # clamp the covariates to these values
      for (i in 1:nlayers(trans_ras)) {
        range <- sry[, colnames(sry) == names(trans_ras)[i]]
        trans_ras[[i]][trans_ras[[i]] < range[1]] <- range[1]
        trans_ras[[i]][trans_ras[[i]] > range[2]] <- range[2]
      }
      
      all_rasters[[n]] <- trans_ras
      count <- count + 1
      
    }
    
    return(all_rasters)
    
  }
  
}



#################################################################################
### Takes Covariates and Data and returns Boosted Regression Trees Outputs
## Inputs:
# df = data.table with variables "latitude", "longitude", "N", and specified indicator
# lcovs = list of raster layers or bricks (for temporally varying).
#         output of load_and_crop_covariates()
# years: analysis years. defaults to 2000, 2005, 2010, 2015
# weight: character of weight variable name in df, defaults to null
# tc: tree.complexity for BRT, defaults to 4
# lr: learning.rate for BRT, defaults to 0.005
# bf: bagging.fraction for BRT, defaults to 0.75
# return_only_raster: if TRUE, only returns raster of results, otherwise returns a list of BRT model object and prediction raster
## Outputs: see above
# Note: depends on gbm and dismo packages
#################################################################################
brt_covs     <- function(df,
                         indicator_family   = indicator_family,
                         lcovs              = cov_layers,
                         years              = c(2000,2005,2010,2015),
                         weight             = NULL,
                         tc                 = 4,
                         lr                 = 0.005,
                         bf                 = 0.75,
                         return_only_raster = TRUE) {
  
  
  
  require(gbm)
  require(dismo)
  #require(parallel)
  #require(foreach)
  #require(doMC)
  
  # take lcovs, see which are bricks (time varying) and which arent
  ntv <- tv <- c()
  for(i in 1:length(lcovs))
    if(inherits(lcovs[[i]],"RasterLayer")) ntv=c(ntv,i) else tv=c(tv,i)
  
  
  # make sure the tv covariates have the same number of years as the years argument
  for(i in tv){
    y=dim(lcovs[[i]])[3]
    if(y==length(years)){
      message(sprintf('The time-varying covariate `%s` has %i years of data. Matching to argument `years`: %s. With layer 1 as earliest year. If this is an issue, please rearrange your input RasterBrick for time-varying covariates. \n\n',names(lcovs)[i], y,paste(years,collapse=', ')))
    } else stop(sprintf('The time-varying covariate `%s` has %i years of data. Which does not match the number of years in argument `years`: %s.',names(lcovs)[i], y,paste(years,collapse=', ')))
  }
  
  # run BRT by years
  message(sprintf('Running BRT on %i years of data independently. Result will be a RasterBrick with %i layers.',length(years),length(years)))
  
  #registerDoMC(cores=length(years))
  #out <- foreach(i = 1:length(years),
  #        .packages=c('dismo', 'gbm', 'raster'),
  #        .export=c('df','years','indicator_family','weight','ntv','tv','lcovs')
  #        ) %dopar% {
  x<-list()
  for(i in 1:length(years)){
    # subset only to years we need
    d <- df[df$year==years[i],]
    d <- na.omit(d)
    
    # keep only the variables of interest
    coords           <- d[, c('longitude', 'latitude'), with=FALSE]
    coords$latitude  <- as.numeric(coords$latitude)
    coords$longitude <- as.numeric(coords$longitude)
    
    # BRT function we use has no binomial, only pois with offset
    if(indicator_family %in% c('binomial','poisson'))  {
      indicator_family = 'poisson'
      offset =log(d$N)
      message('WARNING: For Poisson to work, need to round decimals in the response')
      d[,eval(indicator):=round(d[,eval(indicator),with=FALSE],0)]
    } else offset = NULL
    
    d <- d[,c(indicator,weight),with=FALSE]
    
    # extract the values of the covariates
    c <- brick(lcovs[ntv])
    for(j in tv)
      c <- addLayer(c,lcovs[[j]][[i]])
    d <- cbind(d,extract(c,coords))
    
    
    # learning brt
    set.seed(123)
    # TODO: throw a try-catch so if some years work it at least will return that, if it doesnt it will try different things (like changing the learning rate. )
    mod <- try(
      gbm.step(data             = d,
               gbm.y            = 1,
               gbm.x            = names(c),
               offset           = offset,
               family           = indicator_family,
               weights          = weight,
               tree.complexity  = tc,
               learning.rate    = lr,
               bag.fraction     = bf),silent=TRUE)
    if(is.null(mod)){
      message('First BRT attempt failed. Lowering Learning Rate by 1/10')
      mod <- try(
        gbm.step(data             = d,
                 gbm.y            = 1,
                 gbm.x            = names(c),
                 offset           = offset,
                 family           = indicator_family,
                 weights          = weight,
                 tree.complexity  = tc,
                 learning.rate    = lr*.1,
                 bag.fraction     = bf))
    }
    if(is.null(mod)){
      message('Second BRT attempt failed. Lowering Original Learning rate by 1/1000 AGAIN')
      mod <- try(
        gbm.step(data             = d,
                 gbm.y            = 1,
                 gbm.x            = names(c),
                 offset           = offset,
                 family           = indicator_family,
                 weights          = weight,
                 tree.complexity  = tc,
                 learning.rate    = lr*.001,
                 bag.fraction     = bf))
    }
    if(is.null(mod)){
      message('Third BRT attempt failed. Slow learn plus low tree complexity')
      mod <- try(
        gbm.step(data             = d,
                 gbm.y            = 1,
                 gbm.x            = names(c),
                 offset           = offset,
                 family           = indicator_family,
                 weights          = weight,
                 tree.complexity  = 2,
                 learning.rate    = lr*.001,
                 bag.fraction     = bf))
    }
    
    if(is.null(mod)) stop('ALL BRT ATTEMPTS FAILED')
    
    # save prediction rasters
    p <- predict(c,mod,n.trees=mod$gbm.call$best.trees,type='response')
    
    # save the outputs
    x[[paste0('m_',i)]]=mod
    x[[paste0('brt.',i)]]=p
    #x
    
  } # closes years parallel loop
  
  
  #extract objects and prediction surfaces as two separate lists and save them
  for(i in 1:length(x))
    assign(names(x)[i],x[[i]])
  
  m=(mget(paste0('m_',1:length(years))))
  p=(mget(paste0('brt.',1:length(years))))
  
  
  
  if(return_only_raster) {
    return(brick(p))
  } else {
    return(list(m,brick(p)))
  }
  
  
}



load_gbd_covariates <- function(gbd_fixed_effects,
                                year_ids,
                                gaul_list) {
  
  message(paste0('pulling covs from GBD db for ', gbd_fixed_effects, '...'))
  
  gbd_names <- strsplit(gbd_fixed_effects," ")
  gbd_names <- gbd_names[[1]][gbd_names[[1]] != "+"]
  
  gbd_cov <- function(gbd_name) {
    
    message(paste0('querying covariate_name_short == ', gbd_name, '...'))
    message('make sure this is a valid covariate_name_short from GBD')
    gbd <- load_gbd_data(gbd_type = "covariate",
                         gbd_name = gbd_name,
                         gaul_list = gaul_list,
                         measure_id = 6,
                         age_group_id = 1,
                         metric_id = 3,
                         year_ids = year_ids)
    
    rasterize_gbd_cov <- function(year_option) {
      year_subset <- gbd[year==year_option,]
      year_shape <- merge(subset_shape, year_subset, by.x="GAUL_CODE", by.y="name")
      year_raster <- rasterize(year_shape, simple_raster, field="mean")
      return(year_raster)
    }
    
    brick <- lapply(year_ids, rasterize_gbd_cov)
    brick <- brick(brick)
    names(brick) <- gsub('layer', gbd_name, names(brick))
    
    message(paste0('returning rasterized ', gbd_name))
    return(brick)
    
  }
  
  all_gbd_layers <- lapply(gbd_names, gbd_cov)
  names(all_gbd_layers) <- gbd_names
  
  return(all_gbd_layers)
  
}

load_mbg_covariates <- function(mbg_fixed_effects, simple_raster) {
  
  # Hard-code MBG covariates that have best models available
  mbg_covs <- c('edu_0')
  
  # Load MBG covariates available
  edu_0 <- brick('/share/geospatial/mbg/education/edu_0/output/best/edu_0.tif')
  
  # Construct list of covariates to GAM and use in model from fixed_effects parameter equation.
  selected_covs <- strsplit(mbg_fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]
  num_covs <- length(selected_covs)
  lcovs <- list()
  for(i in 1:num_covs) {
    this_cov <- selected_covs[i]
    lcovs[[i]] <- get(this_cov)
    names(lcovs[[1]]) <- gsub('period_', paste0(this_cov,'.'), names(lcovs[[1]]))
    names(lcovs)[i] <- this_cov
  }
  
  # Make sure covariate layers line up with raster we are modeling over
  for(l in 1:length(lcovs)) {
    lcovs[[l]]  <- crop(lcovs[[l]], extent(simple_raster))
    lcovs[[l]]  <- setExtent(lcovs[[l]], simple_raster)
    lcovs[[l]]  <- mask(lcovs[[l]], simple_raster)
  }
  
  return(lcovs)
}
