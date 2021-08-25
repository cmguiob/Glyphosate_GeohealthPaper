#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 27-08-2020
# Input: Subset of Sentinel-2 resampled with SNAP for bands with res higher than 10m
# Output: CORRELATION PLOT
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script demonstrates how to generate create a correlation plot of covariates
#-------------------------------------------------------------------------------
#PRELIMINARY WORK ON THE IMAGES
# 1. Acquisition of sentinel L2A images  at sci-hub from copernicus
# 2. Exploration, resampling and subsetting (exporting just a set of bands) with SNAP software
# 3. Calculation of vegetation indices, mineral indices and distance raster in R (see 01_Covariates_Sampling.R)

#-------------------------------------------------------------------------------------
# Load indices rasters
#-------------------------------------------------------------------------------------
library(raster)
library(rgdal)
library(sf)

# Find the .tif files in the directory
file_name <- list.files(path = "./Datos/",pattern = 'COVARIATES_.*tif$')
file_list <- paste("./Datos/",file_name, sep = '')

# Load indices which were already cropped to bounding box of reserve
COV_RES_0419 <- stack(x = file_list)
names(COV_RES_0419) <- c("NDVI_0419","PSRI_0419", "SIER_0419", "KIER_0419")
plot(COV_RES_0419)

#-------------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------------

library(ggcorrplot)
library(wesanderson)

# Create stack as df to analyse
COV_RES_df <- as.data.frame(COV_RES_0419, xy= FALSE, na.rm = TRUE)
summary(COV_RES_df)
cor_mx <- cor(COV_RES_df) # CI was left out because of its high similarity with NDVI corrrelating other variables

# Scale range and color
scale_limits <- range(cor_mx, na.rm = TRUE) + c(-0.02, 0.02)
wespal <- rev(wes_palette(name = "Darjeeling1", type = "continuous", n = 3))

png(filename="Anex 2b_Correlation_Covariates.png", 
    type="cairo",
    units="in", 
    width=2.5, 
    height=2.5, 
    pointsize=1, 
    res=400)

ggcorrplot(cor_mx, type = "upper", hc.order = TRUE)+
  scale_fill_gradientn(limits = scale_limits, colors = wespal)+
  theme_bw(base_size = 7)+
  theme(axis.title= element_blank(),
        axis.text = element_text(size= 5.5),
        legend.key.size = unit(0.1, "in"),
        legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 7.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

dev.off()

