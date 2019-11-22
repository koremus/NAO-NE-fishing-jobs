#******************************  METADATA *********************************
# generateFigure_1_2.do
# Script Description: This R script creates Fig. 1B in the main text. 
# Code developed for R version 3.5.3
# 
# Author:  Kimberly L. Oremus
# Date: June 1, 2019
# 
# Data from: 
# NOAA OI SST V2 High Resolution Dataset 
# daily SST values in degrees Celsius over a 0.25 degree latitude by 0.25 degree longitude grid 
# from 1/1/1981 to the present
# https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html
# data located in input folder "sst_Reynolds.nc"
# accessed on January 30, 2018
# 
# NCAR ClimateDataGuide
# DJFM annual NAO index
# https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based
# data located in input folder "NAO_slp_djfm_1864_2017.csv"
# accessed on January 10, 2019
#**************************************************************************/
library(png)
library(grid)
library(sf)
library(maptools)
library(raster)
library(rgdal)
library(dplyr)
library(ncdf4)
library(lattice)
library(RColorBrewer)
library(maps)
library(reshape2)
library(rworldmap)
library(rasterVis)
library(sp)
library(ggplot2)
library("PBSmapping")
library("data.table")
library(ggpubr)

#set working directory
setwd("")
rm(list=ls())
setwd("input")
#load data
ncfile<-nc_open("sst_Reynolds.nc")
print(ncfile)
sst<-ncvar_get(ncfile, "sst")
#lon, lat, time

time<-ncvar_get(ncfile, "T")
lat<-ncvar_get(ncfile,"Y")
lon<-ncvar_get(ncfile,"X")

nlon<-dim(lon)
nlat<-dim(lat)
nt<-dim(time)

for(i in 1:nlon) {
  if(lon[i]>=180){
    lon[i]<-(lon[i]-360);
  }
}

#Base date is Nov 1981
j<-11;
k<-1981;

#assign month and year to each time step
month<-rep(0,nt);
year<-rep(0,nt);
for(i in 1:nt) {
  month[i]<-j;
  year[i]<-k;
  if(j==12){
    j<-1;
    k<-k+1;
  }
  else{j<-j+1;}
}

#create sst average for DJFM by year
yr<-unique(year);
nyr<-length(yr);

sst_DJFM<-array(NA, dim=c(nlon, nlat, nyr-1));
for(i in 1:(nyr-1)){
  sst_DJFM[,,i]<-(sst[,,12*i-10]+sst[,,12*i+1-10]+sst[,,12*i+2-10]+sst[,,12*i+3-10])/4;
}

#load NAO data
nao<-read.csv("nao_pc_djfm_1899_2017.csv")
min_yr<-max(yr[2], nao[1,1])
max_yr<-min(yr[nyr], nao[1,length(nao[,1])])
NAO_min_yr<-which(nao[,1]==min_yr)
NAO_max_yr<-which(nao[,1]==max_yr)


#resize NAO to match sst time series
nao_DJFM<-nao[NAO_min_yr:NAO_max_yr,2]

#creating time trends
t1<-seq(1,length(nao_DJFM),1);
t2<-t1*t1;

#detrending nao
detrendNAO<-lm(nao_DJFM~t1);
nao_dt<-residuals(detrendNAO);

#East Coast lat lon bounds
lon_max<-which(lon==-59.5);
lon_min<-which(lon==-82.5);
lat_min<-which(lat==24.5);
lat_max<-which(lat==49.5);
nlon<-lon_max-lon_min;
nlat<-lat_max-lat_min;
dim<-nlon*nlat;

#correlating NAO to sst
ttdjfm<-data.frame(lon=rep(NA,dim),lat=rep(NA,dim),z=rep(NA,dim))
ct<-0;

for(i in 1:nlon){
  for(j in 1:nlat){
    ct<-ct+1;
    sst_temp<-sst_DJFM[i+lon_min,j+lat_min,];
    detrendSST<-lm(sst_temp~t1);
    sst_dt<-residuals(detrendSST);
    lm_djfm<-lm(sst_dt~nao_dt);
    ttdjfm$z[ct]<-summary(lm_djfm)$coefficients[2,1];
    ttdjfm$lon[ct]<-lon[i+lon_min];
    ttdjfm$lat[ct]<-lat[j+lat_min];
  }
}
    
#MAP
coordinates(ttdjfm)<- ~lon+lat;
gridded(ttdjfm)<-TRUE;
rdjfm<-raster(ttdjfm);
projection(rdjfm)<-crs("+proj=longlat +datum=WGS84")

#plot limits
xlim=c(-82,-59)
ylim=c(25, 49)

#grab administrative land boundaries
worldmap=map_data("world")
setnames(worldmap, c("X", "Y", "PID", "POS", "region", "subregion"))
worldmap<-clipPolys(worldmap, xlim=xlim, ylim=ylim, keepExtra=TRUE)

statemap = map_data("state")
setnames(statemap, c("X","Y","PID","POS","region","subregion"))
statemap = clipPolys(statemap, xlim=xlim,ylim=ylim, keepExtra=TRUE)


#treatment and control regions
states_ne  <-c("maine", "new hampshire", "massachusetts", "rhode island", "connecticut")
states_sa  <- c("south carolina", "georgia", "florida", "north carolina")
statemap2 <- statemap %>%
  mutate(State = case_when(region %in% states_ne ~ "NE",
                           region %in% states_sa ~ "SA"))

#PLOT
gplot(rdjfm)+geom_tile(aes(fill=value))+
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  coord_map(xlim=xlim, ylim=ylim) +
  geom_polygon(data=worldmap,aes(X,Y,group=PID), fill = "grey85",color="white") +
  geom_polygon(data=statemap,aes(X,Y,group=PID),fill = "grey85",color="white") +
  geom_polygon(data=statemap2,aes(X,Y,group=PID, color = State), fill = NA) +
  scale_colour_manual(values = c("NE" = "darkorange", "SA" = "forestgreen")) +
  labs(y="",x="") +
  theme_bw() +
  labs(fill="SST NAO correlation", title="", x="", y="") 

gplot(rdjfm)+geom_tile(aes(fill=value))+
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  coord_map(xlim=xlim, ylim=ylim) +
  geom_polygon(data=worldmap,aes(X,Y,group=PID), fill = "grey85",color="white") +
  geom_polygon(data=statemap,aes(X,Y,group=PID),fill = "grey85",color="white") +
  geom_polygon(data = filter(statemap, region %in% states_ne),aes(X,Y,group=PID),fill = "grey85",color="orange") +
  labs(y="",x="") +
  geom_polygon(data = filter(statemap, region %in% states_sa),aes(X,Y,group=PID),fill = "grey85",color="green3") +
  labs(y="",x="") +
  theme_bw() +
  labs(fill="correlation", title="", x="", y="")

setwd("..")
ggsave("figures/fig1b_winter_NAO_sst_corr.pdf")