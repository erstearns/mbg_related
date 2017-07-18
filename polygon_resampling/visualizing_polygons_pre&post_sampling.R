###########################################################################################
# author: erin stearns
# objective: plot polygons before and after resampling
# date: 11 July 2017
###########################################################################################

## Set arguments for this run and user.
    repo <- '/share/code/geospatial/stearns7/mbg/'
    cores <- 30
    package_list <- c('mgcv', 'nlme', 'survey', 'foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr', 'foreach', 'doParallel')
    indicator <- 'lf_pos_all'
    indicator_group <- 'lf'

## load packages and custom functions.
## There is a difference depending on which nodes you are on.
root        <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                      paste0(root,'temp/geospatial/geos_packages'),
                      paste0(root,'temp/geospatial/packages'))
commondir   <- sprintf('/share/geospatial/mbg/common_inputs')
.libPaths(package_lib)
message(paste0('Loading packages from ',package_lib))
if(Sys.info()[1] == "Windows"){
  stop("STOP! you will overwrite these packages if you run from windows")
} else {
  for(package in c(package_list)) {
    message(paste0('loading ',package))
    library(package, lib.loc = package_lib, character.only=TRUE)
  }
}

# Load  data
collapsed <- fread(paste0(root, 'WORK/11_geospatial/10_mbg/input_data/lf_pos_all.csv'), stringsAsFactors=FALSE)

# format data
collapsed <- collapsed[, start_year := as.numeric(start_year)]
collapsed <- collapsed[grepl('KEN',country), country := 'KEN'] # Just replace IHME Kenya subnational iso3 with national iso3. We only use this to categorize data.
coverage_data <- copy(collapsed)
coverage_data <- coverage_data[(shapefile!="COL_DHS_2005" & shapefile!="GTM_DHS_1987") & shapefile!="PRY_central_wo_asuncion" | is.na(shapefile), ]
coverage_data <- coverage_data[, latitude := as.numeric(latitude)]
coverage_data <- coverage_data[, longitude := as.numeric(longitude)]
coverage_data <- coverage_data[point==1, lf_prev := lf_pos_all/N]

source(paste0(repo, 'mbg_central/graph_data_coverage.R'))

#Plot maps
coverage_maps <- graph_data_coverage_values(df = coverage_data,
                                            var = 'lf_prev',
                                            out_file = "/snfs1/temp/stearns7/lf/mbg/plots/test.png",
                                            title = 'LF',
                                            year_min = '1998',
                                            year_max = '2016',
                                            year_var = 'start_year',
                                            region = 'africa',
                                            sum_by = NULL,
                                            cores = cores,
                                            indicator = 'lf_pos_all',
                                            high_is_bad = TRUE,
                                            return_maps = TRUE,
                                            legend_title = "Prevalence",
                                            extra_file_tag = paste0('_', data_version),
                                            save_on_share = TRUE)


#1. Visualize polygons -- hack up polygon viz from collapse code

#2. Visualize polygons in annual and 5-year stills  - save as raster brick and do steps 3
#   & 4 in arc if short on time

#3. Layer resampled points on top of polygons

#4. Plot points weighted by
