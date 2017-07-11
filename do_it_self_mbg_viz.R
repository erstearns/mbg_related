# plot
peru.shp <- shapefile("/snfs1/WORK/11_geospatial/05_survey shapefile library/Shapefile directory/PER_DHS_2010.shp")
cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)
brks <- seq(0, 1, length = 266)

years <- 2000:2016

pdf(file = sprintf("%s/MSS_Stunting_Peru_2000_2015.pdf", sharedir),
    width = 20, height = 20)
par(mfrow = c(4, 4))
zr <- c(0, max(values(range(mean_raster)), na.rm = T))

# TO-DO ----   make sure # of pixels in raster is less than maxpixels value
for(i in 1:16){
  raster::plot(mean_raster[[i]], main = sprintf("MSS Stunt: %i", years[i]), maxpixels = 1e9, col = mypalette, breaks = brks)
  plot(peru.shp, add = TRUE)
}
dev.off()

png(file = sprintf("%s/MSS_Stunting_Peru_2000_2015.png", sharedir),
    width = 20, height = 20, units = "in", res = 350)
par(mfrow = c(4, 4))
zr <- c(0, max(values(range(mean_raster)), na.rm = T))
for(i in 1:16){
  raster::plot(mean_raster[[i]], main = sprintf("MSS Stunt: %i", years[i]), maxpixels = 1e9, col = mypalette, breaks = brks,
               legend.width = 2,
               axis.args=list(at=seq(0, 1, length = 11),
                              labels=seq(0, 1, length = 11),
                              cex.axis=0.6))
  plot(peru.shp, add = TRUE)
}
dev.off()