cs_covs = extract_covariates(the_data, all_cov_layers, return_only_results = T, centre_scale = T, period_var = 'year', period_map = period_map)


df, covariate_list, id_col = NULL, reconcile_timevarying = T, period_var = NULL, return_only_results = F, centre_scale = F, period_map = NULL

df <- the_data
covariate_list <- all_cov_layers
id_col = NULL
reconcile_timevarying = T
period_var = 'year'
return_only_results = T
centre_scale = T
period_map = period_map

#predict_mbg

res_fit <-  model_fit
cs_df <-  cs_df
mesh_s <-  mesh_s
mesh_t <-  mesh_t
cov_list <-  full_raster_list
samples <-  samples
simple_raster <-  simple_raster
transform <-  transform



set.seed(45)
require(data.table)
DT <- data.table(
  i_1 = c(1:5, NA), 
  i_2 = c(NA,6,7,8,9,10), 
  f_1 = factor(sample(c(letters[1:3], NA), 6, TRUE)), 
  f_2 = factor(c("z", "a", "x", "c", "x", "x"), ordered=TRUE), 
  c_1 = sample(c(letters[1:3], NA), 6, TRUE), 
  d_1 = as.Date(c(1:3,NA,4:5), origin="2013-09-01"), 
  d_2 = as.Date(6:1, origin="2012-01-01"))
# add a couple of list cols
DT[, l_1 := DT[, list(c=list(rep(i_1, sample(5,1)))), by = i_1]$c]
DT[, l_2 := DT[, list(c=list(rep(c_1, sample(5,1)))), by = i_1]$c]

melt(DT, id=1:2, measure="f_1")

melt(DT, id=1:2, measure="1_2")

x <- DT[,c(5:6)]

y <- DT[, names(x), with = FALSE]



#mbg model output

#directories
mbg_out_dir <- ('J:/temp/stearns7/lf/mbg/lf_mbg_model_runs')
write_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/maps/lf_mbg_model_runs/4_6_2017')

#files
mbg_output <- ('/2017_04_06_09_50_25/lf_pos_prediction_eb_bin0_essa_0.grd')
mbg_gtb <- ('/2017_04_06_18_27_20/lf_pos_prediction_eb.grd')
mbg_dn <- ('/2017_04_06_21_27_17/lf_pos_prediction_eb.grd')
west <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/maps/lf_mbg_model_runs/4_7_17_westafrica/lf_pos_prediction_eb.grd')
ssa_essa <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/maps/lf_mbg_model_runs/4_7_17_ssa/lf_pos_prediction_eb.grd')
ssa_wssa <- 
ssa_name <- 
ssa_cssa <- 
ssa_sssa <- 

#reading in 
essa <- raster(paste0(mbg_out_dir, mbg_output))
gtb <- raster(paste0(mbg_out_dir, mbg_gtb))
dn <- raster(paste0(mbg_out_dir, mbg_dn))
wssa <- raster(west)


# plotting
plot(gtb)
plot(wssa)


#making into brick
gtb_brick <- brick(gtb)

#writing to a tiff file
writeRaster(gtb_brick, filename = paste0(write_dir, '/gtb.tif'), format='GTiff', overwrite=T)

eth <- readOGR('J:/WORK/11_geospatial/05_survey shapefile library/Shapefile directory', '2007_ETH_Admin2_Altered')
eth <- fortify(eth, region = "GEO_ID")

plot(essa)
plot(eth)

#reading in africa shapefile
afr <- readOGR('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/spatial_data/africa_0', 'africa_adm0')
afr <- fortify(afr, ADM0_CODE = 'GEO_ID')

#DRC
data_wfake <- fread( file = paste0(save_dir, '/all_data_4.5.17.csv'), stringsAsFactors = F)
drc <- data_wfake[point ==1 & country == 'COD',]
write.csv(drc, file = paste0(save_dir,'/drc.csv'))

plot(afr)
plot(data_wfake$latitude, data_wfake$longitude)

# plots
library(ggplot2)

plot(x=benin$year, y=benin$lf_prev)

ethdata <- data_wfake[country == 'ETH',]
plot(x=ethdata$year, y=ethdata$lf_prev)

gtbdata <- data_wfake[country=='GHA' | country == 'TGO'| country == 'BEN',]
plot(x=gtbdata$year, y=gtbdata$lf_prev)



data = all_poly_data
shp_path = '/home/j/WORK/11_geospatial/05_survey shapefile library/Shapefile directory'
ignore_warnings = TRUE
cores           = 1
indic           = 'lf_pos'
unique_vars     = NULL
density         = 10
perpixel        = TRUE
prob            = TRUE
use_1k_popraster= TRUE 
raster_dict     = rd
raster_col      = "sampling"
use_1k_popraster = TRUE 


# waitformodelstofinish() not stopping, still looking for sssa, but no data there, so no model

#read in each raster

ras.dir <- ('/share/geospatial/mbg/lf/lf_pos/output/2017_04_07_17_59_59') 

name_mean_raster <- raster(paste0(ras.dir, '/lf_pos_prediction_eb_bin0_name_0.grd'))
essa_mean_raster <- raster(paste0(ras.dir, '/lf_pos_prediction_eb_bin0_essa_0.grd'))
wssa_mean_raster <- raster(paste0(ras.dir, '/lf_pos_prediction_eb_bin0_wssa_0.grd'))
cssa_mean_raster <- raster(paste0(ras.dir, '/lf_pos_prediction_eb_bin0_cssa_0.grd'))

#merge

Regions=c('cssa','essa','name', 'wssa') #, 'sssa'

strata <- Regions

for(reg in strata){
  message(reg)
  
  load(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/2017_04_07_17_59_59/', indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))
  
  ## get aggregated estimates for all admin0. Aggregate to level you rake to
  simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = TRUE)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]

  mean_raster <- make_cell_pred_summary( draw_level_cell_pred = cell_pred,
                                         mask                 = simple_raster,
                                         return_as_raster     = TRUE,
                                         summary_stat         = 'median')

  assign(sprintf('%s_mean_raster',reg),mean_raster)
  rm(mean_raster); rm(cell_pred); #rm(pop_wts) #rm(raked_mean_raster); rm(raked_cell_pred); 
  
}




m = do.call(raster::merge,list(name_mean_raster,
                               essa_mean_raster,
                               wssa_mean_raster,
                               cssa_mean_raster))
                             
#write to a raster
writeRaster(
  m,
  file = paste0('/home/j/temp/stearns7/lf/mbg/lf_mbg_model_runs/2017_04_07_17_59_59/mostofssa_mean_raster'),
  format='GTiff',
  overwrite = TRUE)







