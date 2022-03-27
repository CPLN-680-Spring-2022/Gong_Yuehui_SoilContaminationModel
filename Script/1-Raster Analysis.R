R_LIBS_SITE="C:\\Program Files\\R\\R-4.1.2\\library"

install.packages("exact extractr")
install.packages("raster")
library ("exact extractr")
library("raster")
install.packages("RStoolbox")
install.packages("devtools")
library("devtools")
install_github("bleutner/RStoolbox")
library("sp")
library("raster")
library("ggplot2")
library("RStoolbox")
library(sf)



DIR<-"D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Data-satellite\\LC08_L2SP_014032_20160202_20200907_02_T1\\LC08_L2SP_014032_20160202_20200907_02_T1_SR_B"
data_path<-function(file)paste0(DIR,file)

DIR2<-"D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Data-satellite\\LE07_L2SP_014032_20020219_20200917_02_T1\\"
data_path2<-function(file)paste0(DIR2,file)

DIR3<-"D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\Output\\Raster\\"
data_path3<-function(file)paste0(DIR3,file)




#import bands in 2016
b2 <- raster(data_path("2.tif"))
b3 <- raster(data_path("3.tif"))
b4 <- raster(data_path("4.tif"))
b5 <- raster(data_path("5.tif"))
b6 <- raster(data_path("6.tif"))
b7 <- raster(data_path("7.tif"))


r2016 <- paste0(DIR, 1:7, ".tif")
landsat2016 <- stack(r2016)
names(landsat2016) <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7')



#import study area and sample points in 2016
area <- shapefile("D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\Output\\studyarea\\studyarea.shp")
pt2017<-shapefile("D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\Output\\sample2017\\sample2017.shp")

#transform vectors' Coordinate reference system to bands' CRS
areaT<-spTransform(area,crs(b5))
pt2017T<-spTransform(pt2017,crs(b5))

plot(b5, maxpixels=500000, main = "NIR", col = gray(0:100 / 100))
points(pt2017T, col='red', pch=20, cex=1)
plot(areaT, border='white', lwd=2,add=TRUE)


#crop bands using areaT
b2c <- crop(b2, areaT)
b3c <- crop(b3, areaT)
b4c <- crop(b4, areaT)
b5c <- crop(b5, areaT)
b6c <- crop(b6, areaT)
b7c <- crop(b7, areaT)


#plot bands2-7
par(mfrow = c(2,3))
plot(b2c, main = "Band2", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(b3c, main = "Band3", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(b4c, main = "Band4", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(b5c, main = "Band5", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(b6c, main = "Band6", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(b7c, main = "Band7", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)




#spectral indices
#NVDI
vi <- function(x, y) {
  (x - y) / (x + y)
}
ndvi <- overlay(landsat2016[[5]], landsat2016[[4]], fun=vi)
ndvi <-crop(ndvi, areaT)
plot(ndvi, col=rev(terrain.colors(10)), main="2016-Landsat-NDVI")

#hist(ndvi2,main = "Distribution of NDVI values",xlab = "NDVI",ylab= "Frequency",col = "wheat", xlim = c(-0.5, 1),breaks = 30, xaxt = 'n')
#axis(side=1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

#try to filter "open space"
veg <- reclassify(ndvi, cbind(-Inf, 0.2, NA))
plot(veg, main='open space')
plot(areaT, border='dark grey', lwd=2,add=TRUE)


#CMR
FCMR <- function(x, y) {
  x/ y
}
CMR <- overlay(landsat2016[[6]], landsat2016[[7]], fun=FCMR)
CMR <-crop(CMR, areaT)
plot(CMR, col=rev(terrain.colors(10)), main="2016-Landsat-CMR")


#Brightness
Fbright <- function(a, b,c,d,e,f) {
  0.303*a+0.279*b+0.473*c+0.56*d+0.508*e+0.187*f
}
brightness <- overlay(landsat2016[[2]], landsat2016[[3]], landsat2016[[4]],landsat2016[[5]], landsat2016[[6]], landsat2016[[7]],  fun=Fbright)
brightness <-crop(brightness, areaT)
plot(brightness, col=rev(terrain.colors(10)), main="2016-Landsat-brightness")


#Greenness
Fgreen <- function(a, b,c,d,e,f) {
  -0.294*a-0.243*b-0.542*c+0.728*d+0.071*e-0.161*f
}
greenness <- overlay(landsat2016[[2]], landsat2016[[3]], landsat2016[[4]],landsat2016[[5]], landsat2016[[6]], landsat2016[[7]],  fun=Fgreen)
greenness <-crop(greenness, areaT)
plot(greenness, col=rev(terrain.colors(10)), main="2016-Landsat-greenness")
plot(areaT, border='white', lwd=1,add=TRUE)


#MNDWI
FMNDWI <- function(a, b) {
  (a-b)/(a+b)
}
MNDWI <- overlay(landsat2016[[3]], landsat2016[[6]],fun=FMNDWI)
MNDWI <-crop(MNDWI, areaT)
plot(MNDWI, col=rev(terrain.colors(10)), main="2016-Landsat-MNDWI")
plot(areaT, border='white', lwd=1,add=TRUE)


#EVI
FEVI <- function(a, b,c) {
  2.5*(a-b)/(a+6*b-7.5*c+1)
}
EVI <- overlay(landsat2016[[5]], landsat2016[[4]],landsat2016[[2]],fun=FEVI)
EVI <-crop(EVI, areaT)
plot(EVI, col=rev(terrain.colors(10)), zlim=c(-10,10),main="2016-Landsat-EVI")
plot(areaT, border='white', lwd=1,add=TRUE)



#other raster
#DEM.Slope
dem<- raster(data_path3("DEM.tif"))
dem<-crop(dem,areaT)
slope <- raster(data_path3("Slope.tif"))
slope<-crop(slope,areaT)
disroad<-raster(data_path3("dis_road.tif"))
disroad<-crop(disroad,areaT)
dismili<-raster(data_path3("dis_mili.tif"))
dismili<-crop(dismili,areaT)
diswaste<-raster(data_path3("dis_waste.tif"))
diswaste<-crop(diswaste,areaT)
brown<-raster(data_path3("brownfielddensity.tif"))
brown<-crop(brown,areaT)
brown2<-raster(data_path3("brownfielddensity_re.tif"))
brown2<-crop(brown2,areaT)


plot(slope, col=rev(terrain.colors(10)), main="Slope")
plot(dem, col=rev(terrain.colors(10)), main="Elevation")
plot(disroad, col=rev(terrain.colors(10)), main="Distance to Roads")
plot(dismili, col=rev(terrain.colors(10)), main="Distance to Military brownfields")
plot(diswaste, col=rev(terrain.colors(10)), main="Distance to Waste management brownfields")
plot(brown, col=rev(terrain.colors(20)), main="Brownfields Density")

plot(areaT, border='black', lwd=2,add=TRUE)


#Extract raster value into Sample points
pt2017T$b2<-extract(b2c,pt2017T)
pt2017T$b3<-extract(b3c,pt2017T)
pt2017T$b4<-extract(b4c,pt2017T)
pt2017T$b5<-extract(b5c,pt2017T)
pt2017T$b6<-extract(b6c,pt2017T)
pt2017T$b7<-extract(b7c,pt2017T)
pt2017T$ndvi<-extract(ndvi,pt2017T)
pt2017T$CMR<-extract(CMR,pt2017T)
pt2017T$brightness<-extract(brightness,pt2017T)
pt2017T$greenness<-extract(greenness,pt2017T)
pt2017T$MNDWI<-extract(MNDWI,pt2017T)
pt2017T$EVI<-extract(EVI,pt2017T)

pt2017T$dem<-extract(dem,pt2017T)
pt2017T$slope<-extract(slope,pt2017T)
pt2017T$disroad<-extract(disroad,pt2017T)
pt2017T$dismili<-extract(dismili,pt2017T)
pt2017T$diswaste<-extract(diswaste,pt2017T)
pt2017T$brown<-extract(brown,pt2017T)

pt2017test<- as.data.frame(pt2017T)
pt2017test<-pt2017test[-c(1:7,9,11,25,26),] #delete points not in the study area


write.table(pt2017test, file="D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\Output\\2017extracted.CSV", sep=",")






-----------------------------------------------------------------------------------------------
  
  
  
  
#2002
#Import bands1-7 in 2002
cb1 <- raster(data_path2("LE07_L2SP_014032_20020219_20200917_02_T1_SR_B1.TIF"))
cb2 <- raster(data_path2("LE07_L2SP_014032_20020219_20200917_02_T1_SR_B2.TIF"))
cb3 <- raster(data_path2("LE07_L2SP_014032_20020219_20200917_02_T1_SR_B3.TIF"))
cb4 <- raster(data_path2("LE07_L2SP_014032_20020219_20200917_02_T1_SR_B4.TIF"))
cb5 <- raster(data_path2("LE07_L2SP_014032_20020219_20200917_02_T1_SR_B5.TIF"))
cb6 <- raster(data_path2("LE07_L2SP_014032_20020219_20200917_02_T1_ST_B6.TIF"))
cb7 <- raster(data_path2("LE07_L2SP_014032_20020219_20200917_02_T1_SR_B7.TIF"))


#Import Sample Points in 2001 
pt2001<-shapefile("D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\Output\\sample2001\\sample2001.shp")

#CRS transformation
areaT<-spTransform(area,crs(cb5)) #same area
pt2001T<-spTransform(pt2001,crs(cb5))

points(pt2001T, col='red', pch=20, cex=1)
plot(cb5, maxpixels=500000, main = "NIR", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=2,add=TRUE)


#Crop bands
cb1c <- crop(cb1, areaT)
cb2c <- crop(cb2, areaT)
cb3c <- crop(cb3, areaT)
cb4c <- crop(cb4, areaT)
cb5c <- crop(cb5, areaT)
cb6c <- crop(cb6, areaT)
cb7c <- crop(cb7, areaT)

par(mfrow = c(2,3))

plot(cb1c, main = "Band1", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(cb2c, main = "Band2", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(cb3c, main = "Band3", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(cb4c, main = "Band4", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(cb5c, main = "Band5", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)
plot(cb7c, main = "Band7", col = gray(0:100 / 100))
plot(areaT, border='white', lwd=1,add=TRUE)


#Spectral Indices
#NDVI
cndvi <- overlay(cb4c, cb3c, fun=vi)
plot(cndvi, col=rev(terrain.colors(10)), main="2001-Landsat-NDVI")
plot(areaT, border='white', lwd=2,add=TRUE)


#CMR
cCMR <- overlay(cb5c, cb7c, fun=FCMR)
plot(cCMR, col=rev(terrain.colors(10)), main="2001-Landsat-CMR")
plot(areaT, border='white', lwd=2,add=TRUE)


#Brightness
cbrightness <- overlay(cb1c, cb2c, cb3c, cb4c, cb5c, cb7c,  fun=Fbright)
plot(cbrightness, col=rev(terrain.colors(10)), zlim=c(20000,50000),main="2001-Landsat-brightness")
plot(areaT, border='white', lwd=2,add=TRUE)

#Greenness
cgreenness <- overlay(cb1c, cb2c, cb3c, cb4c, cb5c, cb7c, fun=Fgreen)
plot(cgreenness, col=rev(terrain.colors(10)), main="2001-Landsat-greenness")
plot(areaT, border='white', lwd=2,add=TRUE)


#MNDWI
cMNDWI <- overlay(cb2c, cb5c,fun=FMNDWI)
plot(cMNDWI, col=rev(terrain.colors(10)), main="2001-Landsat-MNDWI")
plot(areaT, border='white', lwd=2,add=TRUE)

#EVI
cEVI <- overlay(cb4c, cb3c,cb1c,fun=FEVI)
plot(cEVI, col=rev(terrain.colors(10)), zlim=c(-10,10),main="2001-Landsat-EVI")
plot(areaT, border='white', lwd=2,add=TRUE)


#Raster from above


#Extract raster value for 2001 sample points
#band1 is labeled as "b2" since band1 in Landsat7 have the similar wave length with band2 in Landsat8 and it makes the two Landsat source comparable
pt2001T$b2<-extract(cb1c,pt2001T)
pt2001T$b3<-extract(cb2c,pt2001T)
pt2001T$b4<-extract(cb3c,pt2001T)
pt2001T$b5<-extract(cb4c,pt2001T)
pt2001T$b6<-extract(cb5c,pt2001T)
pt2001T$b7<-extract(cb7c,pt2001T)
pt2001T$ndvi<-extract(cndvi,pt2001T)
pt2001T$CMR<-extract(cCMR,pt2001T)
pt2001T$brightness<-extract(cbrightness,pt2001T)
pt2001T$greenness<-extract(cgreenness,pt2001T)
pt2001T$MNDWI<-extract(cMNDWI,pt2001T)
pt2001T$EVI<-extract(cEVI,pt2001T)

pt2001T$dem<-extract(dem,pt2001T)
pt2001T$slope<-extract(slope,pt2001T)
pt2001T$disroad<-extract(disroad,pt2001T)
pt2001T$dismili<-extract(dismili,pt2001T)
pt2001T$diswaste<-extract(diswaste,pt2001T)
pt2001T$brown<-extract(brown2,pt2001T)


pt2001test<- as.data.frame(pt2001T)
write.table(pt2001test, file="D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\Output\\2001extracted.CSV", sep=",")
