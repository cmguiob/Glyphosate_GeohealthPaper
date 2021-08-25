#-------------------------------------------------------------------------------
# Clip raster to a polygon of study area and generate stratified random sampling design based on covariates
#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 28-08-2020
# Input: Covariates (vegetation and mineral indices)
# Output: K-means based stratified random sampling design
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script creates a K-means stratified (K-means just used to make strata) random sampling design. 
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
# 2. Exploration, resampling and subsetting (exporting just a set of bands) with SNAP software
# 3. Calculation of vegetation indices, mineral indices and distance raster in R (see 01_Covariates_Sampling.R)

#-------------------------------------------------------------------------------
# PROCEDURES:
# 1.	Save a subset of  images as .tif file to a local directory
# 2.	Go to "Session", "Set working directory" and "Choose directory" in the 
#     menu bar of Rstudio. 
#3.   Choose the back the directory where the several data from the project are sttored.
# 4.	Run entire script

#-------------------------------------------------------------------------------
# Load indices rasters
#-------------------------------------------------------------------------------
library(raster)
library(rgdal)
library(sf)

# Find the .tif files in the directory
file_name <- list.files(path = "./Datos/SD/",pattern = 'COVARIATES_.*tif$')
file_list <- paste("./Datos/SD/",file_name[1], sep = '')

# Load indices which were already cropped to bounding box of reserve
COV_RES_0419 <- stack(x = file_list)
names(COV_RES_0419) <- c("NDVI_0419","PSRI_0419", "SIER_0419", "KIER_0419")
plot(COV_RES_0419)


# Reserva, when saved with another CRS
RES_shp<-readOGR(dsn="./Datos/SD/Poligono_Reserva_8419N.shp") # read shape file of polygon
RES_sf <- st_as_sf(RES_shp, coords = c("long", "lat"), crs = 32619)#transform to sf object

#-------------------------------------------------------------------------------
# Load points
#-------------------------------------------------------------------------------

library(sf)
library(rgdal)

Loc <- read.csv("Datos/Locations_all.csv")  

NWconv <- function(i,j,k){i+j/60+k/3600} # create new variable with coordinates
Loc$lat <- NWconv(Loc$NG,Loc$NM,Loc$NS)
Loc$long <- NWconv(Loc$WG,Loc$WM,Loc$WS)*-1
Loc <-Loc[Loc$SITUATION == "NR" & Loc$DAT == "Mai",] #Points sampled in Natural Reserve


# Create spatial points
pt_sf <- st_as_sf(Loc, coords = c("long", "lat"), crs = 4326, agr = "constant") #transform df to sf
pt_sf_trans <- st_transform(pt_sf, crs = 32619,"+proj=longlat +datum=WGS84+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs") # transform coordinate system to fit the raster
pt_sp <- as_Spatial(pt_sf)

# Export as shape
writeOGR(obj=pt_sp, dsn="Locations_all",layer = "Locations_all", driver="ESRI Shapefile", overwrite_layer = TRUE) # this is in geographical projection


#-------------------------------------------------------------------------------
# Unsupervised classification K-Means: create polygons, interpret them and classify land
#-------------------------------------------------------------------------------

# #execute the kmeans function on the rasterStack values and search for 3 to 4 clusters (e.g. centers = 4)
# corresponding to esteros, bancos, madreviejas and land use.

set.seed(200)
km_result <- kmeans(COV_RES_0419[], centers= 6, iter.max=200, nstart = 50) # the brackets are necessary!!


#create a dummy raster using the first layer of our image 
#and replace the values of the dummy raster with the clusters (classes) of the kMeans classification
classmap <- raster(COV_RES_0419[[1]])
classmap <- setValues(classmap, km_result$cluster)

#-------------------------------------------------------------------------------
# PLOT CLASSIFIED IMAGE
#-------------------------------------------------------------------------------
library(ggplot2)

#create df to plot in ggplot
classmap_df <- as.data.frame(classmap, xy = TRUE)
names(classmap_df) <- c("long", "lat", "STRATA")

## Set up color gradient
breaks <- seq(0,6, by=1)
mypal <- colorRampPalette(c("lightblue4","olivedrab3", "lightblue2","orange","olivedrab4","lightgoldenrod2"))(length(breaks)-1)


#Plot with ggplot
png(filename="Figure 2c.png", 
    type="cairo",
    units="in", 
    width=4, 
    height=6.5, 
    pointsize=1, 
    res=400)

ggplot()+
  geom_raster(data = classmap_df, aes(fill = STRATA, x=long, y=lat), alpha = 0.65)+
  scale_fill_gradientn(colors = mypal)+
  geom_sf(data=pt_sf_trans, aes(colour= TARGET), size = 1)+
  scale_colour_manual(values = c("#bc5090", "#ff6361","#003f5c"))+
  theme_minimal(base_size = 9)+
  geom_sf(data = RES_sf, fill = NA, alpha = 0.2, color = "white")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.ticks.length = unit(0, "cm"),
        axis.text = element_text(size= 9),
        axis.title= element_text(size = 10),
        legend.position="top",
        legend.justification="center",
        legend.key.size = unit(0.2, "in"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(1,-10,-10,-10),
        legend.text = element_text(size = 11),
        legend.title = element_blank())+
  guides(fill = FALSE)+
  scale_y_continuous(position = "right")+
  coord_sf(label_axes = list(bottom = "E",right = "N"))


dev.off()

# To plot with rasterVis
#levelplot(classmap, at=breaks, col.regions=mypal, margin = FALSE)+
#layer(sp.polygons(RES_shp, col = "white"))


#Barplot
png(filename="Figure 2d.png", 
    type="cairo",
    units="in", 
    width=4, 
    height=6.5, 
    pointsize=1, 
    res=400)

barplot(classmap,
        main = "Environmental classes",
        ylab= "Pixel frequency",
        col = c("lightblue4","olivedrab3", "lightblue2","orange","olivedrab4","lightgoldenrod2"),
        names.arg = c(  "D-Basin","Schrubs","S-Basin","Levee","Trees", "Splay"),
        las=2)

# Export classified image
writeRaster(classmap, filename = "Kmeans classmap.tif", format="GTiff", overwrite=TRUE )


#-------------------------------------------------------------------------------
# ASSES JUST FOR RESERVE
#-------------------------------------------------------------------------------

# Clip raster to reserva
classmap_RES <- crop(classmap, extent(RES_shp)) 
classmap_RES <- mask(classmap_RES, RES_shp)


#create df to plot in ggplot
classmap_RES_df <- as.data.frame(classmap_RES, xy = TRUE)
names(classmap_RES_df) <- c("long", "lat", "STRATA")

#Plot
#Plot with ggplot
ggplot()+
  geom_tile(data = classmap_RES_df, aes(fill = STRATA, x=long, y=lat))+
  scale_fill_gradientn(colors = mypal)+
  theme_minimal()+
  geom_sf(data = RES_sf, fill = NA, alpha = 0.2, color = "white")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.ticks.length = unit(0, "cm"),
        axis.text = element_text(size= 8),
        axis.title= element_text(size = 9),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))


#Barplot
barplot(classmap_RES,
        main = "Environmental classes",
        xlab = "Classes",
        ylab= "Pixel frequency",
        col = c("orange","olivedrab3", "lightblue2","olivedrab4", "lightgoldenrod2","lightblue4"),
        names.arg = c(  "Levee","Schrubs","S-Basin","Forest","Splay","D-Basin"))

#
#-------------------------------------------------------------------------------
# REMOVE CLASSES THAT WONT BE SAMPLED
#-------------------------------------------------------------------------------

grd_RES <- as(classmap_RES,"SpatialPixelsDataFrame")
gridded(grd_RES) <-F # turn gridded structure to false, means, transform to spatial points df
grd_RES <- as(grd_RES,"data.frame") 
names(grd_RES)[1] <- "env" # geo is a category corresponding to environmental class


for(i in c(2,4,6)){
  ids_RES <- which(grd_RES$env == i) #identify rows that fulfill the condition
  grd_RES <- grd_RES[-ids_RES,] # remove those rows: correspond to shrubs, forest, deep basin
  grd_RES
}


#-------------------------------------------------------------------------------
#select number of samples for each stratum
#-------------------------------------------------------------------------------

#first sort units on geo
unique(grd_RES$env) 
grd_RES <- grd_RES[order(grd_RES$env),] 


#compute stratum sample sizes for proportional allocation
Nc_RES <- tapply(grd_RES$x,INDEX=grd_RES$env,FUN=length) # for each geo class, calculate the number of x coordinates
Pc_RES <- Nc_RES/sum(Nc_RES) # calculates the proportion of each geo class


#set total sample size
n <- 15
Nsamp <- round(Pc_RES*n) # number of samples for each geo class
sum(Nsamp)

#-------------------------------------------------------------------------------
# Generate points
#-------------------------------------------------------------------------------
library(sampling)

# Points are generated at the center of each pixel. Using jitter to randomly assign withing each pixel was problematic.
set.seed(102)
units <-strata(grd_RES,stratanames="env",size=Nsamp,method="srswr") # method is simple random sampling with replacement

mysample <-getdata(grd_RES,units) # getdata extracts the observed data from a data frame

#-------------------------------------------------------------------------------
# PLot
#-------------------------------------------------------------------------------
library(ggrepel)
library(ggplot2)

# Add ID to the samples
mysample$Name <- paste("102S",seq(1,15,by =1),sep = '')

# Set colors
mycolors <- c("orange","lightblue2", "lightgoldenrod2")

#Plot clusters and sampling points

ggplot(grd_RES) +
  geom_tile(mapping = aes(x = x, y = y, fill = factor(env))) +
  scale_fill_discrete(name = "env") +
  scale_fill_manual(values=mycolors) +
  geom_point(data=mysample,mapping=aes(x=x,y=y),size=2) +
  geom_text_repel(aes(x= x, y = y,label = Name), data = mysample, size =3)+
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "") +
  coord_fixed() +
  theme(legend.position="none")

#Barplot
barplot(Nsamp,
        main = "Environmental classes",
        xlab = "Classes",
        ylab= "Sample number",
        col = c("orange","lightblue2", "lightgoldenrod2"),
        names.arg = c("Levee","Basin","Splay"))


#-------------------------------------------------------------------------------
# Write
#-------------------------------------------------------------------------------
library(plotKML)

SRS_xy <- mysample[,c(1,2,7)]

write.csv(SRS_xy,file="103SRS_xy.csv",row.names=F)
plotKML(classmap_RES)
writeRaster(classmap_RES, filename = "KM6_Clases_Reserva", format = "GTiff")
