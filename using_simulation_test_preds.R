######################################
# Code author: Erin Stearns
# Code objective: Investigating prediction values for LF prevalence
# Date: 6.27.2017
#####################################

rm(list = ls())

# -------------------------------------------------- set up environment -------------------------------------------
my_lib <- ("C:/Users/stearns7/Documents/R/win-library/3.3")
.libPaths(my_lib)
pacman::p_load(fields, data.table, magrittr, stringr, reshape2, ggplot2, readbulk, 
               dplyr, Amelia, plyr, rgdal, raster,  seegSDM, seegMBG, sp, rgeos, maptools, rgeos, rgdal)

data_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/data/processed_data/latest/')
afr_dir <- 'C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/output_rasters/africa/2017_06_06_14_35_57/'
save_dir <- ('C:/Users/stearns7/OneDrive - UW Office 365/ntds/LF/mbg_stuff/investigations/prediction_v_raw/')

# ------------------------------------------------ load data --------------------------------------------------------------

#reading in csv of latest dataset
lf <- fread(paste0(data_dir, 'lf_data_w_gahi_nids.csv'), stringsAsFactors = F)

#reading in csv of africa only dataset
model2pts <- fread(file = paste0(afr_dir, 'raster_vals_to_original_data_pts.csv'), stringsAsFactors = F)

# ---------------------------------------------- clean up as necessary --------------------------------------------------
model2pts <- model2pts[RASTERVALU != -9999,]
model2pts[,lf_prev_o:=lf_pos/N]
afr <- model2pts[ country != 'FJI',]
afr[,pred_prev:=RASTERVALU]
afr[,obs_prev:=lf_prev_o]
# ------------------------------------------------ get sense of data to inform simulation ---------------------------------
summary(afr)

a <- copy(afr)

hist(a$N, breaks=20)

n_dens <- density(a$N)
plot(n_dens)

hist(a$N, breaks=40, probability=TRUE)
lines(n_dens, color='red')
rug(a$N)

ggplot(a, aes(N)) +
  geom_freqpoly(binwidth=500)



#get mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode(a$lf_pos)

#plot distributions
# plotvars <- names(a[,c(8:9, 16:17)])
# 
# for (p in plotvars) {
#   basic <- ggplot(a, aes(x=a[[p]])) +
#     geom_point(shape=16) + 
#     labs(x = (p))
#   loes <- basic + geom_smooth()
#   ggsave(basic, filename = paste0(save_dir, 'basic_', p, '.png'))
#   ggsave(loes, filename = paste0(save_dir, 'loes_', p, '.png'))
# }

# ------------------------------------------------ simulate similar dataset -----------------------------------------------
#500 draws from rbinom(n=30, p=(vector of ps))
#Median N == 100, mean prev == 0.06912

#p <- rep(c(0.02, 0.12,0.06912,0.05708), times = 25)
#p <- runif(100, min=0, max=0.1)


#get vector of probs
p <- (rbeta(100, 0.15, 5))

#get estimated x
x <- rbinom(500, size=50, prob = p)

#estimate phat (=x/30)
phat <- (x/100)

#plot phat vs p
sim <- data.table(x,p,phat)

ggplot(sim, aes(x=p, y=phat)) +
  geom_point(colour='red', alpha=0.25) +
  labs(x='True p', y='Estimated p') +
  geom_smooth()  
  # xlim(0,0.1) +
  # ylim(0,0.1)

ggsave(filename = paste0(save_dir,'n_50_p_range.png'))


# simulating varying denominators

#get a vector of denominators
denoms <- c(0,1,5,10,50,100,500)

for (d in denoms){
  p <- (rbeta(100, 0.15, 5))
  x1 <- rbinom(n=500, size = d, prob = p)
  phat1 <- (x1/d)
  sim1 <- data.table(x1,p,phat)
  plot <- ggplot(sim, aes(x=p, y=phat)) +
    geom_point(colour='red', alpha=0.25) +
    labs(x='True p', y='Estimated p') +
    geom_smooth()  +
    xlim(0,0.1) +
    ylim(0,0.1)
}



#plotting histo of N

barfill <- "#4271AE"
barlines <- "#1F3552"

p7 <- ggplot(a, aes(x = N)) +
  geom_histogram(aes(y = ..count..), binwidth = 5,
                 colour = barlines, fill = barfill) +
  scale_x_continuous(name = "Denominators",
                     breaks = seq(0, 2000, 250),
                     limits=c(0, 2000)) +
  scale_y_continuous(name = "Count") +
  ggtitle("Frequency histogram of median LF denominator") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 9),
        axis.text.y=element_text(colour="black", size = 9)) +
  geom_vline(xintercept = 100, size = 1, colour = "#FF3721",
             linetype = "dashed")
p7

#look at cameroon
cmr <- a[country=='CMR',]


