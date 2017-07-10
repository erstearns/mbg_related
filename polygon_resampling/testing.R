if(run_points_only == 'TRUE'){
  message("You're config file said to run a points-only model, so we're doing just that!")
  model_script_pts <- 'launch'
  # name the points only indicator (to be specified in the config file)
  #indicator_pts <- points_only_indicator
  # set up remaining qsub args
  slots <- as.numeric(slots)
  sharedir       <- sprintf('/share/geospatial/mbg/%s/%s',indicator_group,indicator)
  fin_pts_file_name <- fin_pts_only_run_name
  njobs <- 1
  # qsub points only model
  make_qsub_par(ig = indicator_group,
                indic = config_name
                indic2pass = indicator,
                code = model_script_pts,
                geo_nodes = use_geos_nodes,
                log_location = 'sharedir',
                fin_file_name = fin_pts_file_name)
  # waiting for raster
  output_dir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date)
  Sys.sleep(600)
  check_loc_results(c(1:njobs), check_dir=output_dir, prefix="fin_",postfix=".txt")
} else {
  message("You're config file said not to run a points only model, so we're kicking it to polygon resampling! ")
  next
}


#testing polygon subsetting

d <- read.csv()
lf_pos_all <- copy(e1)
d <- copy(e1)


d <- d[point==0,]
d <- d[!is.na(shapefile),]
message(paste0('You have ',nrow(d), ' rows of data in your dataset.'))


region <- 'africa'

if(region=='africa'){
  africa_iso <- c('DZA','AGO','SHN','BEN','BWA',
                  'BFA','BDI','CMR','CPV','CAF',
                  'TCD','COM','COG','COD','DJI',
                  'EGY','GNQ','ERI','ETH','GAB',
                  'GMB','GHA','GIN','GNB','CIV',
                  'KEN','LSO','LBR','LBY','MDG',
                  'MWI','MLI','MRT','MUS','MYT',
                  'MAR','MOZ','NAM','NER','NGA',
                  'STP','REU','RWA','STP','SEN',
                  'SYC','SLE','SOM','ZAF','SSD',
                  'SHN','SDN','SWZ','TZA','TGO',
                  'TUN','UGA','COD','ZMB','TZA',
                  'ZWE')
  d1 <- d[country %in% africa_iso,]
  message(paste0('There are ', nrow(d1), ' polygons in ', region))
}else{
  message('Region needs to be defined')
}

# testing out how to delete certain files from a specific folder
output_dir <- 'C:/Users/stearns7/Documents/test1/'
del <- list.files(output_dir, pattern = "fin_")

unlink(paste0(output_dir,del))


# testing something
x <- read.csv(file = 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/polygon_resampling/testing_new_structure/lf_pos_all_dat.csv')
x1 <- x[[1]]

df <- df_list[[1]]
run_date <- df_list[[2]]


#testing something else
config <- fread('C:/Users/stearns7/Documents/git_repos/new_mbg/mbg/lf/config_lf_pos.csv', stringsAsFactors = F)
as.character(config$fixed_effects)
selected_fixed_effects <-
  strsplit(config$fixed_effects," ")[[1]][config$strsplit(fixed_effects," ")[[1]] != "+"]

