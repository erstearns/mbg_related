#############################################
# Author: Erin Stearns
# Objective: Pull out BRT output for examination
# Date: 6/20/2017
############################################

rm(list = ls())

# -------------------------------------------------- set up environment -------------------------------------------
my_lib <- ("C:/Users/stearns7/Documents/R/win-library/3.3")
.libPaths(my_lib)
pacman::p_load(fields, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, 
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, seegMBG, sp, rgeos, maptools, rgeos, rgdal)

# ------------------------------------------------ load data --------------------------------------------------------------
load('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/mbg_output/2017_06_06_14_35_57/child_model_list_wssa_lf_0.RData')

#smo_cssa <- child_models
#smo_essa <- child_models
#smo_wssa <- child_models

# ------------------------------------------------ brt output -------------------------------------------------------------
gbm_cssa <- smo_cssa[['gbm']]
gbm_essa <- smo_essa[['gbm']]
gbm_wssa <- smo_wssa[['gbm']]

# ------------------------------------------------ examining output -------------------------------------------------------

this_config <- fread(paste0('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/mbg_output/2017_06_06_14_35_57/config.csv'))
stacker_list <- this_config[V1 == 'stacked_fixed_effects', V2]
stackers_used <- strsplit(stacker_list," ")
stackers_used <- stackers_used[[1]][stackers_used[[1]] != "+"]

load('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/mbg_output/2017_06_06_14_35_57_bin0_cssa_lf_0.RData')

fit.data <- as.data.frame(df)
## we want the part of the df after the stackers
fit.data <- fit.data[, -(1:max(match(stackers_used, colnames(fit.data))))]



