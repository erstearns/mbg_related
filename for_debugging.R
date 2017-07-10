##########################
# erin stearns
# de-bugging mbg runs
# 24 Febe - 8 March 2017
##########################

rm(list = ls())

# setting the stage ----------------------------------------------------------------------- SET UP -------------------
#packages for testing things locally
pacman::p_load(fields, foreign, rgeos, data.table, raster, rgdal, INLA, seegSDM, seegMBG, plyr, dplyr, tidyr, magrittr)
setwd('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/R_WD')

main_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/main_in_use')
spatial_main <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/main_in_use/spatial_data')
priority8 <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/priority8')



# for data viz in 'collapse' script ------------------------------------------------------- DATA COVERAGE MAP --------------------
df <- coverage_data
var <- 'lf_prev'
out_file <- "/snfs1/temp/stearns7/lf/mbg/test.png"
title <- 'NTD'
year_min <- '1980'
year_max <- '2016'
year_var <- 'start_year'
region <- 'africa'
sum_by <- 'n'
cores <- 10
indicator <- 'lf_prev'
high_is_bad <- TRUE
return_maps <- TRUE
legend_title <- "Prevalence"

# for gbm model fit - sri lanka ----------------------------------------------------------- BRT ISSUE - sri lanka data ---------------------------
df <- the_data
model_name <- 'gbm'
fold_id_col <- 'fold_id'
covariates <- all_fixed_effects
weight_column <- 'weight'
tc <- 4 #tree complexity
lr <- 0.0005 #learning rate
bf <- 0.75 #bag fraction
indicator <- indicator
indicator_family <- indicator_family
cores <- slots
additional_terms = NULL
bam =F
spline_args = list(bs = 'ts', k = 3)
auto_model_select =T


#writing out the truncated sri lanka dataset after omitting rows with na cov vals -------- EXTRACT covs to data issues for Sri Lanka ----------------

write.csv(the_data, file = '/home/j/temp/stearns7/lf/mbg/sri_lanka_truncated_pts_3.2.csv')

#reading it in to compare to:
sl_pts <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/mbg_data/sri_lanka/lf_pos.csv')

#compare to:
sl_pts_trun <- fread(file = 'J:/temp/stearns7/lf/mbg/sri_lanka_truncated_pts_3.2.csv')
save(sl_pts_trun, file = paste0(priority8, 'weird_sri_trunc.rda'))
sl_pts_trun$in_trun <- 1

#join the 2 and see where in_trun is NA

sher <- join(sl_pts, sl_pts_trun)

# rows being lost

losing <- sher[is.na(in_trun),]
write.csv(losing, file = paste0(priority8, '/losing_lka_pts.csv'))

#covariate extraction -------------------------------------------------------------------- EXTRACT function ----------------------------------------

for(i in names(all_cov_layers)){
   #extract the rasters by points
  locations = SpatialPoints(coords = as.matrix(df[,.(longitude,latitude)]), proj4string = CRS(proj4string(all_cov_layers[[1]])))
  cov_values = as.data.frame(lapply(all_cov_layers, function(x) extract(x, locations)))
  df = cbind(df, cov_values)
  return(df)
}

# pop raster
pop <- raster("J:/temp/stearns7/lf/mbg/input/afripop_&_asiapop_1km.grd")

plot(pop)


# polygon resampling
# temp: testing polygon code so reloading successful file - ## locate the path to the mean raster. we'll need that
output_dir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/2017_03_02_12_57_01')
path.to.risk.rast <- paste0(output_dir, '/', indicator,'_prediction_eb')


all_poly_data <- fread(paste0(root,'/temp/stearns7/lf/mbg/lka/all_lka.csv')) #change for each country/source

data = all_poly_data
shp_path = '/home/j/WORK/11_geospatial/05_survey shapefile library/Shapefile directory'
ignore_warnings = TRUE
cores           = 1
indic           = 'lf_pos'
unique_vars     = NULL
density         = 0.001
perpixel        = TRUE
prob            = TRUE
use_1k_popraster= TRUE
raster_dict = rd
raster_col = "sampling"

# why is this shapefile loading NA values for GAUL code -----------------------------------------------------------------------------
test <- shapefile("J:/WORK/11_geospatial/05_survey shapefile library/Shapefile directory/SriLanka_LF_Gunawardena_2014.shp")

summary(test)

#or these -- 3. - cluster address:
shp_path <-  '/home/j/WORK/11_geospatial/05_survey shapefile library/Shapefile directory'
g2002 <- shapefile(paste0(shp_path, '/g2015_2002_2', '.shp'))
g2008 <- shapefile(paste0(shp_path, '/g2015_2008_2', '.shp'))
g2010 <- shapefile(paste0(shp_path, '/g2015_2010_2', '.shp'))
sri_custom <- shapefile(paste0(shp_path, '/SriLanka_LF_Gunawardena_2014', '.shp'))

#for local investigation
shp_path <-  'J:/WORK/11_geospatial/05_survey shapefile library/Shapefile directory'
g2002 <- shapefile(paste0(shp_path, '/g2015_2002_2', '.shp'))
g2008 <- shapefile(paste0(shp_path, '/g2015_2008_2', '.shp'))
g2010 <- shapefile(paste0(shp_path, '/g2015_2010_2', '.shp'))
sri_custom <- shapefile(paste0(shp_path, '/SriLanka_LF_Gunawardena_2014', '.shp'))

#poly <- fread(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/mbg_data/sri_lanka/poly_lka.csv', stringsAsFactors = F)
all <- fread(file = 'J:/temp/stearns7/lf/mbg/lka/all_lka.csv', stringsAsFactors = F)
poly <- all[point == 0,]

#checking that all poly ids exist in shapefiles
x <- g2002[g2002$GAUL_CODE == 25848, ] #exists

g2008x <- as.data.table(g2008)
g2008_polys <- poly[shapefile == 'g2015_2008_2', location_code]
y <- g2008x[GAUL_CODE %in% g2008_polys,] #all there

z <- g2010[g2010$GAUL_CODE == 25852,] #all there

# visualizing gtb data --------------------------------------------------------------------------
y <- lf_mbg_pts[country == 'Ghana'| country == 'Togo' | country == 'Benin',]

quilt.plot(y$latitude, y$longitude, y$N)

#losing gtb points same as sri lanka, want to investigate 3.8.17 ----------------------- EXTRACT covs to data issues - gtb -----------------

#original = y above
y <- fread(file = paste0(priority8, '/mbg_ptsonly.csv'), stringsAsFactors = F)
y1 <- y[country == 'Ghana'| country == 'Togo' | country == 'Benin',]

#compare to:
gtb_trun <- fread(file = 'J:/temp/stearns7/lf/mbg/gtb_truncated_pts_3.8.csv', stringsAsFactors = F)
gtb_trun$in_trun <- 1

trun_uid <- gtb_trun[,cluster_id]

y2 <- y1[!cluster_id %in% trun_uid,]
write.csv(y2, file = paste0(main_dir, '/weird/pts_losing_gtb_3.8.csv'))

#join the 2 and see where in_trun is NA
sher <- full_join(y1, gtb_trun)


# for gbm model fit - gtb -------------------------------------------------------------- BRT ISSUE - gbt data -------------------------------------------------
df <- the_data
model_name <- 'gbm'
fold_id_col <- 'fold_id'
covariates <- all_fixed_effects
weight_column <- 'weight'
tc <- 4 #tree complexity
lr <- 0.0005 #learning rate
bf <- 0.75 #bag fraction
indicator <- indicator
indicator_family <- indicator_family
cores <- slots
additional_terms = NULL
bam =F
spline_args = list(bs = 'ts', k = 3)
auto_model_select =T

#if run line by line, runs through (???) why?
library(mgcv)
#remove scoping surprises
df = copy(df)

#start by fitting the full gam
message('Fitting the Full GAM model')
full_model = fit_gam(df,
                     covariates = covariates,
                     additional_terms = additional_terms,
                     weight_column = weight_column,
                     bam = bam,
                     spline_args = spline_args,
                     auto_model_select =auto_model_select,
                     indicator = indicator,
                     indicator_family = indicator_family,
                     cores = cores)

#add a name to the game object
full_model$model_name = model_name

#fit the child/cv rfs
folds = unique(df[,get(fold_id_col)])

for(fff in folds){
  baby_model = fit_gam(df = df[get(fold_id_col) != fff,],
                       covariates=covariates,
                       additional_terms = additional_terms,
                       bam = bam,
                       weight_column = weight_column,
                       spline_args = spline_args,
                       auto_model_select = auto_model_select,
                       indicator = indicator,
                       indicator_family = indicator_family,
                       cores = cores)

  #fill in the data
  df[get(fold_id_col)==fff, paste0(model_name,'_cv_pred') := predict(baby_model, df[get(fold_id_col)==fff,],type = 'response')]
}

#predict using full model fit earlier
df[,paste0(model_name,'_full_pred') := predict(full_model,df,type = 'response')]

#return a subset of the columns. Full pred denotes the fit from the full model. CV pred is the OOS stuff
suffixes = c('_full_pred', '_cv_pred')
return_cols = paste0(model_name, suffixes)

gbm <- (setNames(list(df[,return_cols,with = F], full_model),c('dataset',paste0(model_name))))

# other countries
all <- fread(file = paste0(priority8, '/mbg_ptsonly.csv'), stringsAsFactors = F)

zam <- all[country == 'Zambia',]
quilt.plot(zam$longitude, zam$latitude, zam$N)

# polygon function not working --------------------------------------------------------- RESAMPLE POLYGONS - not working ----------------------


data <-  all_poly_data
cores <-  1
indic <-  'lf_pos'
raster_dict <-  rd
raster_col <-  "sampling"
shp_path <-  '/home/j/WORK/11_geospatial/05_survey shapefile library/Shapefile directory'
ignore_warnings <-  TRUE
unique_vars     <-  NULL
density         <-  10
perpixel        <-  TRUE
prob            <- TRUE
use_1k_popraster <- TRUE


# polygon function not working --------------------------------------------------------- RESAMPLE POLYGONS, getPoints - not working ----------------------

shape    = poly
raster   =
n = density
perpixel = perpixel
prob = prob


#exp
shape    = poly
raster   = popraster[[3]]
n        = density
perpixel = perpixel
prob     = prob
# same issue as above, testing locally
#
raster   = ifelse(length(raster_list[[which.rast]]) == 1, # <- issue here ers
                  raster_list[[1]],
                  raster_list[[3]])
# blargggghhh
n.rast <- data.frame(n = 1:n_chunk, rast = pars[, raster_col])

which.rast <- as.character(n.rast[2])
x <- as.numeric(n.rast[1])
poly_name <- paste(pars[x,c('shapefile', 'location_code')], collapse = '__')
poly <- polys[[poly_name]]

#changing args - removing dens as arg and just assigning n to be what we assigned in original function call

getPointsWrapper.w.dict <- function(x){ #this is the one that should be running ------------------------------------
  which.rast <- as.character(x[2])
  x <- as.numeric(x[1])
  poly_name <- paste(pars[x,c('shapefile', 'location_code')], collapse = '__')
  poly <- polys[[poly_name]]

  ## NOTE: simply resampling on 2010 population data for now for convenience. Could redo by year, but wouldnt make much of difference anyway.
  if (is.null(poly)) {
    ## make a dummy, empty points dataframe
    points <- data.frame(longitude = NA, latitude  = NA, weight    = NA)[0, ]
  } else {
    rast <-
      points <- try(
        getPoints(shape    = poly,
                  raster   = raster_list[[which.rast]],
                  n        = density,
                  perpixel = perpixel,
                  prob     = prob) )
    if(inherits(points,'try-error')){
      points <- data.frame(longitude = NA, latitude  = NA, weight    = NA)[0, ]
    } else {
      colnames(points) <- c('longitude', 'latitude', 'weight')
    }
  }
  return(points)
}

 #chunk stuff
 points <- chunk_points[[1]]
