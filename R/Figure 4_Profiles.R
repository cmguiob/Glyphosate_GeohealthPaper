#-------------------------------------------------------------------------------
# SOIL PROFILE MODELING - PAZ DE ARIPORO
#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 10-10-2019
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script creates soil profile collection, groupes it by landform  and 
# plots it with glyphosate depth locaitons.
#-------------------------------------------------------------------------------
# OBJECTIVES:
# 1. Munsell property modelling
# 2. Plotting soil profiles
# 3. Calssification


#-------------------------------------------------------------------------------
#PRELIMINARY WORK ON DATA
# Data was collected on paper formats and digitalized to a .csv table following the format
# of a soil profile collection.

#-------------------------------------------------------------------------------
# LOAD DATA 
#-------------------------------------------------------------------------------

# Load soil horizon data
PARS <- read.csv("./Datos/TCIPA_Suelos_Hz.csv")  

# Load soil site data
PARS_Site <- read.csv("Datos/TCIPA_Suelos_Site.csv")


#-------------------------------------------------------------------------------
# FORMAT DATA AND CREATE GLOBAL VARIABLES 
#-------------------------------------------------------------------------------

library(aqp)

# Create color variables per horizon - hex colors need to be as characters
PARS$RGBmx <- munsell2rgb(PARS$MX_H, PARS$MX_V, PARS$MX_C)
PARS$RGBmo <- munsell2rgb(PARS$MOT1_H, PARS$MOT1_V, PARS$MOT1_C)

# Create long-lat variables
NWconv <- function(i,j,k){i+j/60+k/3600} # create new variable with coordinates
PARS_Site$lat <- NWconv(PARS_Site$NG,PARS_Site$NM,PARS_Site$NS)
PARS_Site$long <- NWconv(PARS_Site$WG,PARS_Site$WM,PARS_Site$WS)*-1

# Extract NF and glyphosate variables to plot as brackets later
PARS_NF1 <- PARS_Site[,c("ID","NF1_top","NF1_base")]
PARS_NF2 <- PARS_Site[,c("ID","NF2_top","NF2_base")]
PARS_GLY <- PARS_Site[,c("ID","Gly_top","Gly_base")]
names(PARS_NF1) <- c("ID", "top", "bottom")
names(PARS_NF2) <- c("ID", "top", "bottom")
names(PARS_GLY) <- c("ID", "top", "bottom")

#-------------------------------------------------------------------------------
# SET SOIL PROFILE COLLECTION
#-------------------------------------------------------------------------------
library(aqp)
depths(PARS) <- ID ~ TOPE + BASE

PARS$LANDFORM <- PARS_Site$LANDFORM # Alternative, create a site level variable
PARS$LANDFORM <- factor(PARS$LANDFORM, levels = c("LEVEE", "SPLAY","UPPER-BASIN", "BASIN"))
PARS$lat <- PARS_Site$lat # first add it to site data to be able to convert them to profile coordinates
PARS$long <- PARS_Site$long
print(PARS)



#Subset for profiles in the natural reserve (NR)

PARS$SITUATION <- PARS_Site$SITUATION #Add situation to the profile collection
idx <- which(PARS$SITUATION == 'NR') # explicit string matching
PARS<- PARS[idx, ] # perform subset based on index

#initialize coordinates
coordinates(PARS) <- ~ long + lat
coordinates(PARS)
print(PARS)

# set spatial reference system
proj4string(PARS) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
proj4string(PARS) 
print(PARS)

# create profile - landform label
df.label <- site(PARS)
df.label$label <- paste(df.label$ID,df.label$LANDFORM, sep = "-" )
site(PARS) <- df.label

# See soil profile locations
PARS.sp <- as(PARS, 'SpatialPointsDataFrame')

# plot the fake coordinates
plot(PARS.sp)
box()

#-------------------------------------------------------------------------------
# VISUALIZATION OF SOIL PROFILES: 
# Profiles should visualize water tables, glyphosate occurrence and horizon should be 
# labeled with texture/structure index. Below Infiltration rates ... could validate interflow transport
# where valeus are low and hanging water horizons occur.

# Profiles classified based on surface geo-chemical properties (dendrogram) are more or less consistent 
# with biogeographical model (NDVI shown below), which determined the identification of landforms for sampling representativeness. 
# Nevetheless: below the surface, complexity is such, that it has little relevance to the 
# occurrence of shallow-groundwater and glyphosate... so, it makes no sense to show the dendrogram,
# nor the landforms.
#-------------------------------------------------------------------------------
library(aqp)

# ***********Plot grouped version**********

# Save as tiff
tiff("Figure 4.tiff", type = "cairo", width=7.5, height=4.5, units = "in", res = 400)

par(mar=c(0.2,0.2,0,0.7))
groupedProfilePlot(PARS,  groups = "LANDFORM",group.name.offset = -12, group.line.col = "lightgrey", group.name.cex = 0.65, 
                   name = "TEXT_USDA",cex.names = 0.6, color = 'RGBmx',divide.hz=FALSE, label = 'ID',  cex.depth.axis = 0.65, )
addVolumeFraction(PARS, "PRC_MOT1",pch = 18, col = PARS@horizons$RGBmo, cex.min = 0.2, cex.max = 0.3)
addBracket(PARS_NF1, col = '#3aa6b4', lwd = 1.5, missing.bottom.depth = 25)
addBracket(PARS_NF2, col = '#3aa6b4', lwd = 1.5, missing.bottom.depth = 25)
addBracket(PARS_GLY,  lwd=2, col = "#FC4E07",  missing.bottom.depth = 25)


dev.off() # Finish saving plot

# Save as .svg
svg("Figure 4.svg", width=7.5, height=4.5)

par(mar=c(0.2,0.2,0,0.7))
groupedProfilePlot(PARS,  groups = "LANDFORM",group.name.offset = -12, group.line.col = "lightgrey", group.name.cex = 0.65, 
                   name = "TEXT_USDA",cex.names = 0.6, color = 'RGBmx',divide.hz=FALSE, label = 'ID',  cex.depth.axis = 0.65, )
addVolumeFraction(PARS, "PRC_MOT1",pch = 18, col = PARS@horizons$RGBmo, cex.min = 0.2, cex.max = 0.3)
addBracket(PARS_NF1, col = '#3aa6b4', lwd = 1.5, missing.bottom.depth = 25)
addBracket(PARS_NF2, col = '#3aa6b4', lwd = 1.5, missing.bottom.depth = 25)
addBracket(PARS_GLY,  lwd=2, col = "#FC4E07",  missing.bottom.depth = 25)


dev.off() # Finish saving plot

