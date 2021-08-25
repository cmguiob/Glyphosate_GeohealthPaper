#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 13-07-2020
# Output: MODELS OF 4 SOIL PROFILES WHERE GLYPHOSATE OCCURRED
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script plots 4 soil profiles with location of glyphosate and groundwater,
# showing horizon labels corresponding to field texture.
# This version is a simplification of the profiles studied, summarized for the
# case study presented in the report for DeJusticia.
#-------------------------------------------------------------------------------
# OBJECTIVES:
# 1. Munsell property modelling
# 2. Plotting soil profiles
# 3. Location of glyphosate at depth and size proportional to concentration
# 4. Location of groundwater levels
# 5. Label horizons with field texture
#-------------------------------------------------------------------------------
#PRELIMINARY WORK ON DATA
# Data was collected on paper formats and digitalized to a .csv table following the format
# of a soil profile collection.

#-------------------------------------------------------------------------------
# LOAD DATA 
#-------------------------------------------------------------------------------
# Load soil horizon data
PARS <- read.csv("Datos/TCIPA_Suelos_Hz.csv")  

# Subset just the profiles that have glyphosate
PARS <- PARS[PARS$GLY > 0,]

# Load soil site data
PARS_Site <- read.csv("Datos/TCIPA_Suelos_Site.csv")

# Subset just the profiles that have glyphosate
PARS_Site <- PARS_Site[PARS_Site$GLY_CONC > 0.1,]


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

# Create meaan depth value for plotting glyphosate occurence
calcmean <- function(x,y){(x+y)/2}
PARS_Site$Gly_med <- calcmean(PARS_Site$Gly_top, PARS_Site$Gly_base)

# Extract NF and glyphosate variables to plot as brackets later
PARS_NF1 <- PARS_Site[,c("ID","NF1_top","NF1_base")]
PARS_NF2 <- PARS_Site[,c("ID","NF2_top","NF2_base")]
PARS_GLY <- PARS_Site[,c("ID","Gly_med", "GLY_CONC")]
names(PARS_NF1) <- c("ID", "top", "bottom")
names(PARS_NF2) <- c("ID", "top", "bottom")



#-------------------------------------------------------------------------------
# SET SOIL PROFILE COLLECTION
#-------------------------------------------------------------------------------
library(aqp)


depths(PARS) <- ID ~ TOPE + BASE

PARS$GEOFORMA <- PARS_Site$GEOFORMA # Alternative, create a site level variable
PARS$lat <- PARS_Site$lat # first add it to site data to be able to convert them to profile coordinates
PARS$long <- PARS_Site$long
print(PARS)

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
df.label$label <- paste(df.label$ID,df.label$GEOFORMA, sep = "-" )
site(PARS) <- df.label


#-------------------------------------------------------------------------------
# VISUALIZATION OF SOIL PROFILES
#-------------------------------------------------------------------------------

#Create points x axis
x.pos <- 1:length(PARS)

#Create vector fot sizes
GLY_size <- c(0.38, 0.5, 0.68, 0.40)

pdf("Alternative_Figure 4_ Perfiles DeJusticia-PECIG.pdf", height = 6, width = 6)
par(mar=c(1,0.5,1,1.5))
plot(PARS, name = "TEXTURA", width = 0.15, color = 'RGBmx',divide.hz=FALSE, label = 'label', cex.depth.axis = 0.9, cex.names = 0.9)
addVolumeFraction(PARS, "PRC_MOT1",pch = 18, col = PARS$RGBmo, cex.min = 0.35, cex.max = 0.45)
addBracket(PARS_NF1, lwd = 2, col = '#3aa6b4', label.cex = 0.75, missing.bottom.depth = 25)
points(x.pos, PARS_GLY$Gly_med, col='white', bg="#E96149",pch=21, cex= GLY_size*6)

dev.off()
