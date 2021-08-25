#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 27-08-2020
# Input: Subset of Sentinel-2 resampled with SNAP for bands with res higher than 10m
# Output: MINERAL AND VEGATATION INDICES
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script demonstrates how to generate specific mineral and vegetation indices from 
# from Sentinel 2 images from 6.April.2019 (preprocesses with SNAP) to use them as environmental 
# covariates for the sampling design.
#-------------------------------------------------------------------------------
# OBJECTIVES:
# 1. Load and confirm data 
# 2. Calculate mineral and vegetation indices
# 3. Create and export correlation plot for Anex
# 4. Plot and export a mineral and a vegetation index for conceptual presentation
# 5. Export indices to use them in the next step (K-means clustering)
#-------------------------------------------------------------------------------
#PRELIMINARY WORK ON IMAGES
# 1. Acquisition of sentinel L2A images  at sci-hub from copernicus
# 2. Exploration, resampling and subsetting with SNAP: Sentinel images were resampled to 10m 
#    with Neares Neighbours method and subsetted to a smaller area and specific bands with SNAP 
#    software, before exporting them as GEOTIFF-BigTIF.
# NOTE: The complete output from SNAP, including images from September and December, resampled with
#        bilinear splines are stored in the GIS folder of this project.
#-------------------------------------------------------------------------------
# LOAD DATA 
#-------------------------------------------------------------------------------
library(raster)

# Find the .tif files in the directory
file_name <- list.files(path = "./Datos/SD/",pattern = 'subset_.*tif$')
file_list <- paste("./Datos/SD/",file_name[1], sep = '')

# Load sentinel 2A images
S2A_NN_0419 <- stack(x = file_list)
names(S2A_NN_0419)
names(S2A_NN_0419) <- c("NN_0419_B2","NN_0419_B3","NN_0419_B4","NN_0419_B5","NN_0419_B6",
                   "NN_0419_B7","NN_0419_B8", "NN_0419_B8A", "NN_0419_B9","NN_0419_B11",
                   "NN_0419_B12", "NN_0419_AOT","NN_0419_WV","NN_0419_CC")

#-------------------------------------------------------------------------------
#Corroborate that the bands were correctly assigned plotting an RGB natural colors
#-------------------------------------------------------------------------------# 
plotRGB(S2A_NN_0419, r=3, g=2, b=1, stretch = 'lin') 

#-------------------------------------------------------------------------------
# VEGETATION INDEX NDVI: from NIR B8 (842nm) or NIR B8a (865nm), and  red B4 (665nm)
#-------------------------------------------------------------------------------

# Calculate NDVI from reflectance. For k = NIR  and i = Red bands in the stack
VI <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi) / (bk + bi)
  return(vi)
}

NDVI_0419 <- VI(S2A_NN_0419, 7, 3)
names(NDVI_0419) <- "NDVI"
plot(NDVI_0419, main = "NDVI")

NDVI865_0419 <- VI(S2A_NN_0419, 8, 3)
names(NDVI865_0419) <- "NDVI865"
plot(NDVI865_0419, main = "NDVI 865nm")

#-------------------------------------------------------------------------------
# PLANT STRESS INDEX Chlorophyll red-edge: bands NIR red-edge: B7 (783nm), B5(705nm)
# Low CI indicates low chlorophyll content and severe stress. 
# Source Zhang et al. 2018
#-------------------------------------------------------------------------------

# Calculate chlorofile index from reflectance. For k = B7  and i = B5 red-edge NIR bands in the stack
VI2 <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi2 <- (bk/bi) -1
  return(vi2)
}

CI_0419 <- VI2(S2A_NN_0419, 6, 4)
names(CI_0419) <- "CI"
plot(CI_0419, main = "CI-red edge")

#-------------------------------------------------------------------------------
# PLANT STRESS INDEX PSRI red-edge B6(740nm), red B4(665nm)  and blue B2 (490nm)
# This index maximizes the sensitivity to the ration of bulk carothenoids to chlorophyll.
# An increased PSRI indicates more canopy stress.
# Source: Zhang et al. 2018.
#-------------------------------------------------------------------------------

# Calculate PSRI index from reflectance. For k = red B4, i = blue B2  and j = B6 NIR bands in the stack

VI3 <- function(img, k, i, j) {
  bk <- img[[k]]
  bi <- img[[i]]
  bj <- img[[j]]
  vi3 <- (bk/bi)/bj
  return(vi3)
}

PSRI_0419 <- VI3(S2A_NN_0419, 3, 1, 5)
names(PSRI_0419) <- "PSRI"
plot(PSRI_0419, main = "PSRI")


#-------------------------------------------------------------------------------
# SABINS IRON ENHANCEMENT RATIO -SIER (Sabins, 1999): original TM-5,  van der Werff  and van der Meer (2016) - adaptated to OLI-8 
# It uses the bands red B4 (665nm) and blue B2 (490nm), both are originally 10m resolution
# Sabins suggests also a ratio for claya, but it uses the bands B11 (1610) and B12 (2190), which are originally 20m res
# Likewise, the Soil Survey Manual (2017) recommends several normalized ratios (Table 5-2),
# which mostly use SWIR bands. For Sentinel 2,these bands are 20m resolution, and therefore excluded.
#-------------------------------------------------------------------------------

# Calculate Sabine's iron oxide index from reflectance. For k = red and i = blue bands  in the stack
FeI <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  feox <- bk/bi
  return(feox)
}

SIER_0419 <- FeI(S2A_NN_0419, 3, 1)
names(SIER_0419) <- "SIER"
plot(SIER_0419, main = "Sabins' Fe3 Enhancement Ratio")

#-------------------------------------------------------------------------------
# KAUFMANN IRON ENHANCEMENT RATIO -KIER (Kaufmann, 1988) uses the bands SWIR B12 (2190nm) and NIR B8A (865nm)
# This bands could introduce an error since bands 12 and 8a are originally 20m resolution and
# were resampled to 10m with Near Neighbours.
#-------------------------------------------------------------------------------

# Calculate Kaufmann's iron oxide index from reflectance. For k = SWIR and i = NIR bands  in the stack

KIER_0419 <- FeI(S2A_NN_0419, 9, 8)
names(KIER_0419) <- "KIER"
plot(KIER_0419, main = "Kaufmann's  Fe3 Enhancement Ratio")

#-------------------------------------------------------------------------------
# DISTANCES: A raster of distances between suspected punctual pollution sources
# could be a covariate for point source pollution (e.g. from a dump), but it was not 
# used here since glyphosate is non-point source pollution.
#-------------------------------------------------------------------------------

# Calculate the shortest distances in meters to  points
# DIST <- distanceFromPoints(object = BASE_RASTER, xy = POINTS)


#-------------------------------------------------------------------------------
## Clip raster to bounding box
#-------------------------------------------------------------------------------
library(rgdal)

# Reserva, when saved with another CRS
RES_shp<-readOGR(dsn="./Datos/SD/Poligono_Reserva_8419N.shp") # read shape file of polygon

# Stack covariates
COV_stack <- stack(NDVI_0419,NDVI865_0419,CI_0419,PSRI_0419,SIER_0419,KIER_0419)

# Clip raster to bounding box reserva
extent_RES <- extent(RES_shp) + 2000 # extend reserva by 2km
COV_RES_bbox <- crop(COV_stack, extent_RES)
COV_RES <- mask(COV_RES_bbox, RES_shp)

#-------------------------------------------------------------------------------
# PLOT AND EXPORT INDICES IMAGES (for presentation)
#-------------------------------------------------------------------------------
library(ggplot2)
library(sf)

## NOTE: Ideally, all covariates should be plotted with levelplot, but it has an error at the moment

#Convert to df for ggplot
NDVI_bb_df <- as.data.frame(COV_RES_bbox[[1]], xy = TRUE)
names(NDVI_bb_df) <- c("long", "lat", "NDVI")

SIER_bb_df <- as.data.frame(COV_RES_bbox[[5]], xy = TRUE)
names(SIER_bb_df) <- c("long", "lat", "SIER")

#Convert polygon to sf
RES_sf <- st_as_sf(RES_shp, coords = c("long", "lat"), crs = 32619)#transform to sf object


#Plot NDVI
png(filename="Figure 2a_NDVI.png", 
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
  geom_sf(data = RES_sf, fill = "white", alpha = 0.2, color = "#8B8682")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.ticks.length = unit(0, "cm"),
        axis.text = element_text(size= 9),
        axis.title= element_text(size = 9),
        legend.position="top",
        legend.justification="center",
        legend.key.size = unit(0.4, "in"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(1,-10,-10,-10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))

dev.off() # Finish saving plot

#Plot SIER

png(filename="Figure 2b_SIER.png", 
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
  geom_sf(data = RES_sf, fill = NA, alpha = 0.2, color = "white") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.ticks.length = unit(0, "cm"),
        axis.text = element_text(size= 9),
        axis.title= element_text(size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.position="top",
        legend.justification="center",
        legend.key.size = unit(0.4, "in"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(1,-10,-10,-10))

dev.off() # Finish saving plot

#-------------------------------------------------------------------------------
# Write the rasters as .tif files
#-------------------------------------------------------------------------------

# Write rasters of large area
COV_stack4_a <- stack(NDVI_0419, PSRI_0419,SIER_0419,KIER_0419) # covariate stack without distance to classify geoforms later
COV_stack4_b <- stack(NDVI865_0419,PSRI_0419,SIER_0419,KIER_0419) # covariate stack without distance to classify geoforms later


writeRaster(NDVI_0419, filename="NDVI_20190406.tif", format="GTiff", overwrite=TRUE)
writeRaster(NDVI865_0419, filename="NDVI865_20190406.tif", format="GTiff", overwrite=TRUE)
writeRaster(PSRI_0419, filename="PSRI_20190406.tif", format="GTiff", overwrite=TRUE)
writeRaster(SIER_0419, filename="SIER_20190406.tif", format="GTiff", overwrite=TRUE)
writeRaster(KIER_0419, filename="KIER_20190406.tif", format="GTiff", overwrite=TRUE)
writeRaster(COV_stack4_a, filename = "COVARIATES_stack4_a.tif", format = "GTiff", overwrite = TRUE)
writeRaster(COV_stack4_b, filename = "COVARIATES_stack4_b.tif", format = "GTiff", overwrite = TRUE)


# Write rasters of bbox
writeRaster(COV_RES_bbox[[1]], filename="NDVI_20190406_bbox.tif", format="GTiff", overwrite=TRUE)
writeRaster(COV_RES_bbox[[2]], filename="NDVI865_20190406_bbox.tif", format="GTiff", overwrite=TRUE)
writeRaster(COV_RES_bbox[[4]], filename="PSRI_20190406_bbox.tif", format="GTiff", overwrite=TRUE)
writeRaster(COV_RES_bbox[[5]], filename="SIER_20190406_bbox.tif", format="GTiff", overwrite=TRUE)
writeRaster(COV_RES_bbox[[6]], filename="KIER_20190406_bbox.tif", format="GTiff", overwrite=TRUE)
writeRaster(COV_RES_bbox[[c(1,4,5,6)]], filename = "COVARIATES_stack4a_bbox.tif", format = "GTiff", overwrite = TRUE)
writeRaster(COV_RES_bbox[[c(4,4,5,6)]], filename = "COVARIATES_stack4b_bbox.tif", format = "GTiff", overwrite = TRUE)
