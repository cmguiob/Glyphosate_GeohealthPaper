#-------------------------------------------------------------------------------
# Clip raster to a polygon of "reserva" study area and generate K-Means sampling 
# design based on covariates. Based on the tutorial and paperfrom Dick Brus 
#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 27-08-2020
# Input: Covariates (vegetation and mineral indices)
# Output: K-means centered stratified sampling design (not random)
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script creates a sampling design based on K-means centers. 
# The covariates were selected to represent the geomorphology and vegetation of the study area. 
# Covariates were calculated from Sentinel 2 images aqcuired on April the 6th 2019, during the 
# summer season.
#----------------------------------------------------------------------------------------
#RESOURCES
# Based on the tutorial and paper  from Dick Brus 
#-------------------------------------------------------------------------------
#PRELIMINARY WORK ON THE POLYGONS

# 1. Create .csv file from coordinate points
# 2. Open .csv file in QGIS. Export as .shp in desired CRS
# 3. Open .shp of points in google earth
# 3. Draw polygon in google earth and save as .kml file. 
# 4. Load .kml polygon and points in QGIS and export polygon as .shp in desired CRS
#    - sometimes points need to be opened and exported from R to be read in QGIS

#PRELIMINARY WORK ON THE IMAGES
# 1. Acquisition of sentinel L2A images  at sci-hub from copernicus
# 2. Exploration, resampling and subsetting with SNAP software
# 3. Calculation of vegetation indices, mineral indices and distance raster in R (see 01_Covariates_Sampling.R)

#-------------------------------------------------------------------------------
# PROCEDURES:
# 1.	Save a subset of  images as .tif file to a local directory
# 2.	Go to "Session", "Set working directory" and "Choose directory" in the 
#     menu bar of Rstudio. 
#3.   Choose the back the directory where the several data from the project are sttored.
# 4.	Run entire script

library(raster)
library(sf)
library(sp)
library(rgdal)
library(fields)
library(ggplot2)
library(rasterVis)

#-------------------------------------------------------------------------------
# Load indices rasters
#-------------------------------------------------------------------------------

# Find the .tif files in the directory
file_name <- list.files(path = "./INPUT_02_KMeans/",pattern = 'COVARIATES_.*tif$')
file_list <- paste("./INPUT_02_KMeans/",file_name, sep = '')

# Load sentinel 2A images
COV_0419 <- stack(x = file_list)
names(COV_0419) <- c("NDVI_0419","PSRI_0419", "SIER_0419", "KIER_0419")
plot(COV_0419)


#-------------------------------------------------------------------------------
# Clip raster to polygon: reserva and contaminada1
#-------------------------------------------------------------------------------

# Reserva, when saved with another CRS
Reserva_shp<-readOGR(dsn="./INPUT_02_KMeans/Poligono_Reserva_8419N.shp") # read shape file of polygon

plot(COV_0419[[1]])
plot(Reserva_shp, add = TRUE)

# Clip raster to reserva
extent_Res<- extent(Reserva_shp) + 2000
COV_Res_bbox <- crop(COV_0419, extent_Res)
COV_Res <- mask(COV_Res_bbox, Reserva_shp)

#-------------------------------------------------------------------------------
# PLOT AND EXPORT A MINERAL AND A VEGETATION INDEX
#-------------------------------------------------------------------------------

#Convert to df for ggplot
NDVI_bb_df <- as.data.frame(COV_Res_bbox[[1]], xy = TRUE)
names(NDVI_bb_df) <- c("long", "lat", "NDVI")

SIER_bb_df <- as.data.frame(COV_Res_bbox[[3]], xy = TRUE)
names(SIER_bb_df) <- c("long", "lat", "SIER")

# Extract single layer for ggplot
NDVI <- COV_Res_bbox[[1]]
SIER <- COV_Res_bbox[[3]]

#Convert polygon to sf
Res_sf <- st_as_sf(Reserva_shp, coords = c("long", "lat"), crs = 32619)#transform to sf object


#Plot NDVI
png(filename="Figure 2a.png", 
    type="cairo",
    units="in", 
    width=4, 
    height=6.5, 
    pointsize=1, 
    res=300)

ggplot()+
  geom_raster(data = NDVI_bb_df, aes(fill = NDVI, x=long, y=lat))+
  scale_fill_gradient2(high= '#00AFBB', mid='#FDF5E6', low='#FC4E07', midpoint = 0.35)+
  theme_minimal()+
  geom_sf(data = Res_sf, fill = "white", alpha = 0.2, color = "#8B8682")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.ticks.length = unit(0, "cm"),
        axis.text = element_text(size= 8),
        axis.title= element_text(size = 9),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

dev.off() # Finish saving plot

#Plot SIER

png(filename="Figure 2b.png", 
    type="cairo",
    units="in", 
    width=4, 
    height=6.5, 
    pointsize=1, 
    res=300)

ggplot()+
  geom_raster(data = SIER_bb_df, aes(fill = SIER, x=long, y=lat))+
  scale_fill_gradient2(low= '#00AFBB', mid='#FDF5E6', high='#FC4E07', midpoint = 1.27)+
  theme_minimal()+
  geom_sf(data = Res_sf, fill = NA, alpha = 0.2, color = "white")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.ticks.length = unit(0, "cm"),
        axis.text = element_text(size= 8),
        axis.title= element_text(size = 9),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

dev.off() # Finish saving plot


#-------------------------------------------------------------------------------
# Generate sampling points for "reserva"
#-------------------------------------------------------------------------------

#Set number of sampling locations to be selected for the contaminated an the blank area
n <- 15

#Compute clusters
set.seed(300)
myClusters <- kmeans(scale(Cov_Res_df[,3:6]), centers=n, iter.max=100,nstart=30)
Cov_Res_df$clusters <- myClusters$cluster


#Select locations closest to the centers of the clusters
rdist.out <- rdist(x1=myClusters$centers,x2=scale(Cov_Res_df[,3:6]))
ids.mindist <- apply(rdist.out,MARGIN=1,which.min)
mySample <- Cov_Res_df[ids.mindist,]

#Plot clusters and sampling points

ggplot(Cov_Res_df) +
  geom_tile(mapping = aes(x = x, y = y, fill = factor(clusters))) +
  scale_fill_discrete(name = "cluster") +
  geom_point(data=mySample,mapping=aes(x=x,y=y),size=2) +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "") +
  coord_fixed() +
  theme(legend.position="none")


ggplot(Cov_Res_df) +
  geom_point(mapping=aes(y=SFeox,x=PSRI,colour=factor(clusters))) +
  geom_point(data=mySample,mapping=aes(y=SFeox,x=PSRI),size=2) +
  scale_y_continuous(name = "SFeox") +
  scale_x_continuous(name = "PSRI") +
  theme(legend.position="none")


#-------------------------------------------------------------------------------
# Extract sampling points
#-------------------------------------------------------------------------------

KM_xy <- mySample[,1:2]
KM_xy$Name <- paste("PARS",seq(01,15,by =1),sep = '') # To export later as .gpx, ID must be named "Name" 


# Before writing sampling coordinates, review fields. For the .gpx, just x, y, and name
write.csv(KM_xy, file = "KM_xy.csv", row.names = FALSE)



