library(readxl)
Sample2017verified<-read_excel("C:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\RawData\\PAH_sample_2017\\PAH Sampling Location Information.xlsx")
Sample2017metal<-read_excel("C:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\RawData\\PAH_sample_2017\\Metal results by county_all_data 02-15-2022 1338.xlsx", sheet = "CAMD")
Sample2017location<-subset(Sample2017verified, select = c("Sample ID","Easting","Northing"))
Sample2017metal<-Sample2017metal[,2:50]
samplemetal_transpose = t(Sample2017metal)
library(janitor)
samplemetal_transpose<- row_to_names(samplemetal_transpose,row_number = 1)
colnames(samplemetal_transpose)[1] <- "Sample ID"
sample2017CAMD<-merge(samplemetal_transpose,Sample2017location[48:61,],by="Sample ID")
#sample2017CAMD includes sample sites location and metal measurement for each site

sample2017CAMDshallow<-subset(sample2017CAMD,sample2017CAMD$"Sample Depth" =="shallow")
for (t in (4:5,9:30)) {
sample2017CAMDshallow[,t] <- as.numeric(sample2017CAMDshallow[,t] )
}

write.table(sample2017CAMDshallow, file="C:\\001\\PennMLA\\2022 spring\\Advanced GIS\\Git\\Gong_Yuehui_SoilContaminationModel\\ChangedData\\sample2017CAMDshallow.csv", sep=",")

###
library("tidyverse")
ggplot(sample2017CAMDshallow)+
  geom_point(data = sample2017CAMDshallow, 
             aes(x = "Cobalt", y = "Cadmium"))

library(tidycensus)
library(sf)
ggplot()+
  geom_sf(data = sample2017CAMDshallow.sf)
