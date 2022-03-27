install.packages("GGally")
install.packages("outliers")
install.packages("caTools")
install.packages("terra")
install.packages("rgdal")

library(outliers)
library("GGally")
library("caTools")
library("ggplot2")
library("raster")
library("terra")
library(raster)
library(rgdal)

#Import Sample points with extracted values
pt2017test<-read.csv(file="D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\Output\\2017extracted.CSV", sep=",")
pt2001test<-read.csv(file="D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\Output\\2001extracted.CSV", sep=",")



#As

As<-rbind(pt2017test[,c(11,33:50)],pt2001test[,c(4,29:46)])

#test outliers
hist(As$Arsenic,
     main = "Arsenic Histogram",
     xlab = "Arsenic Amount",
     ylab= "Frequency",
     col = "grey",
)

test <- grubbs.test(As$Arsenic)
test <- grubbs.test(As$Arsenic,opposite = TRUE)
test
#35.6 (p-value = 1.558e-06) is an outlier

#outlier removed
As<-rbind(pt2017test[-c(1),c(11,33:50)],pt2001test[,c(4,29:46)])


plot(As, col="purple", main="Plotting Pairs Against Each Other")

#ggpairs see correlation
ggpairs(As)

#select 80% data for training, 20% for testing
sample = sample.split(As$Arsenic, SplitRatio = 0.8)
train = subset(As, sample == TRUE)
test = subset(As, sample == FALSE)
train_size = dim(train)
test_size = dim(test)

ggpairs(train)

Model <- lm(Arsenic ~ b4 + b3 + b2 + b5 + brown + brightness + MNDWI + greenness, data = train)
summary(Model)


#Predicting
data_size = dim(As)
pred <- predict(Model, test)
numx <- data_size[1]*(1 - 0.8)+0.4
x_axis <- seq(numx)
df <- data.frame(x_axis,pred,test$Arsenic)

#Plotting the predicted values against the actual values
g <- ggplot(df, aes(x=x_axis))
g <- g + geom_line(aes(y=pred, colour="Predicted"))
g <- g + geom_point(aes(x=x_axis, y=pred, colour="Predicted"))
g <- g + geom_line(aes(y=test$Arsenic, colour="Actual"))
g <- g + geom_point(aes(x=x_axis, y=test$Arsenic, colour="Actual"))
g <- g + scale_colour_manual("", values = c(Predicted="red", Actual="blue"))
g

original = test$Arsenic
predicted = pred
d = original-predicted
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((original-mean(original))^2))
cat(" MAE:", mae, "\n", "MSE:", mse, "\n", 
    "RMSE:", rmse, "\n", "R-squared:", R2)

#extract the parameters of the estimated regression equation
coeffs = coefficients(Model); coeffs 



#use model to visualize As distribution
#plot As in 2016
Asmodel <- function(a,b,c,d,e,f,g,h) {
  coeffs[1]+coeffs[2]*a+coeffs[3]*b+ coeffs[4]*c +coeffs[5]*d + coeffs[6]*e +coeffs[7]*f +coeffs[8]*g +coeffs[9]*h
}
Arsenic <- overlay(b4c,b3c,b2c,b5c,brown2,brightness, MNDWI,greenness,fun=Asmodel)

breakpoints <- c(-100,0,3,6,10,20,200)
colors <- c("black","dark green","green","yellow","orange","red")
par(mfrow = c(1,1))
plot(Arsenic,breaks=breakpoints,col=colors,main="Arsenic Distribution in 2016")
plot(areaT, border='white', lwd=2,add=TRUE)

#plot As in 2002
Arsenic <- overlay(cb3c,cb2c,cb1c,cb4c,brown2,cbrightness, cMNDWI,cgreenness,fun=Asmodel)
Arsenic <-crop(Arsenic, areaT)
breakpoints <- c(-100,0,3,6,10,20,50)
colors <- c("black","dark green","green","yellow","orange","red")
par(mfrow = c(1,1))
plot(Arsenic,breaks=breakpoints,col=colors,main="Arsenic Distribution in 2002")
plot(areaT, border='white', lwd=2,add=TRUE)

#plot As in 2022
DIR4<-"D:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Data-satellite\\LC09_L2SP_014032_20220210_20220212_02_T1\\LC09_L2SP_014032_20220210_20220212_02_T1_SR_B"
data_path4<-function(file)paste0(DIR4,file)
#import bands in 2022
tb1 <- raster(data_path4("1.tif"))
tb2 <- raster(data_path4("2.tif"))
tb3 <- raster(data_path4("3.tif"))
tb4 <- raster(data_path4("4.tif"))
tb5 <- raster(data_path4("5.tif"))
tb6 <- raster(data_path4("6.tif"))
tb7 <- raster(data_path4("7.tif"))
tb1c <- crop(tb1, areaT)
tb2c <- crop(tb2, areaT)
tb3c <- crop(tb3, areaT)
tb4c <- crop(tb4, areaT)
tb5c <- crop(tb5, areaT)
tb6c <- crop(tb6, areaT)
tb7c <- crop(tb7, areaT)

plot(tb2c)
plot(areaT, border='white', lwd=2,add=TRUE)

tMNDWI <- overlay(tb3c, tb6c,fun=FMNDWI)
tbrightness <- overlay(tb2c, tb3c, tb4c,tb5c, tb6c, tb7c,  fun=Fbright)
tgreenness <- overlay(tb2c, tb3c, tb4c,tb5c, tb6c, tb7c,  fun=Fgreen)

Arsenic <- overlay(tb4c,tb3c,tb2c,tb5c,brown2,tbrightness,tMNDWI,tgreenness,fun=Asmodel)
breakpoints <- c(-100,0,3,6,10,20,200)
colors <- c("black","dark green","green","yellow","orange","red")
par(mfrow = c(1,1))
plot(Arsenic,breaks=breakpoints,col=colors,main="Arsenic Distribution in 2022")
plot(areaT, border='white', lwd=2,add=TRUE)


#03/26/2022
--------------------------------------------------------------------------------------------------------


#Pb
#Pb<-pt2017test[-c(9),c(20,33:44)]
Pb<-rbind(pt2017test[-c(9),c(20,33:50)],pt2001test[-c(10),c(19,29:46)])

plot(Pb, col="purple", main="Plotting Pairs Against Each Other")

test <- grubbs.test(Pb$Lead)
test

boxplot(Pb$Lead,
        ylab = "Lead"
)
hist(Pb$Lead,
     main = "Lead Histogram",
     xlab = "Lead Amount",
     ylab= "Frequency",
     col = "grey",
)
axis(side=1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

ggpairs(Pb)



#Cu
Cu<-pt2017test[,c(18,33:44)]

plot(Cu, col="purple", main="Plotting Pairs Against Each Other")

test <- grubbs.test(Cu$Copper)
test <- grubbs.test(Cu$Copper,opposite = TRUE)
test

boxplot(Cu$Copper,
        ylab = "Cu"
)
hist(Cu$Copper,
     main = "Cu Histogram",
     xlab = "Cu Amount",
     ylab= "Frequency",
     col = "grey",
)
axis(side=1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

ggpairs(Cu[-c(18),])


#Cr
Cr<-pt2017test[,c(16,33:44)]

plot(Cr, col="purple", main="Plotting Pairs Against Each Other")

test <- grubbs.test(Cr$Chromium)
test <- grubbs.test(Cr$Chromium,opposite = TRUE)
test

boxplot(Cr$Chromium,
        ylab = "Cr"
)
hist(Cr$Chromium,
     main = "Cr Histogram",
     xlab = "Cr Amount",
     ylab= "Frequency",
     col = "grey",
)
axis(side=1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

ggpairs(Cr)

