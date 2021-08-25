#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 23-09-2019
# Output: LOCATION MAP
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script  plots a map with inset that shows the region and study
# area, as well as the distribution of sampling points in the area. This
# Figure wasn't used in the end for the paper, but it still can be useful
# for very short presentations.
#-------------------------------------------------------------------------------
# OBJECTIVES:
# 1. Polygon of colombia
# 2. Polygon of casanare
# 3. Locations as points
# 4. Sentinel data from study area -> compute NDVI
# 5. Polygon of study area
# 6. DEM of colombian mountains
# 7. Overlay raster on polygon ir polygon (of casanare) on raster: map of colombia with mountains and region of study
# 8. Overlay locaitons on the study area
# 9. Create inset map of 6 in 7

#-------------------------------------------------------------------------------
# LOAD DATA : example: readRDS("path/to/file/xxx.rds")
#-------------------------------------------------------------------------------

library(sp)

# 1. Polygon od Colombia
spl1_Col <- readRDS("Datos/GADM/gadm36_COL_1_sp.rds") # Level 1 (departments) data from GADM.org provided by the university davis
spl0_Col <- readRDS("Datos/GADM/gadm36_COL_0_sp.rds")

sfl1_Col <- readRDS("Datos/GADM/gadm36_COL_1_sf.rds") # Level 1 (departments) data from GADM.org provided by the university davis
sfl0_Col <- readRDS("Datos/GADM/gadm36_COL_0_sf.rds")

#Check plot
spl1_Col@data
regionalValues <- runif(26, min = 0.5, max = 0.9)  # Simulate random  values for each region between 0.5 and 0.9
plot(spl1_Col, col = gray(regionalValues), border = 0)

#-------------------------------------------------------------------------------
# 2. Polygon of Casanare

sp_csn <- spl1_Col[9,] # extraer el poligono de casanare
sf_csn <- sfl1_Col[sfl1_Col$NAME_1 == "Casanare",]

#Check plot
plot(sp_csn)

#-------------------------------------------------------------------------------
# 3. Location points

library(sf)

Loc <- read.csv("Datos/Locations_all.csv")  
NWconv <- function(i,j,k){i+j/60+k/3600} # create new variable with coordinates
Loc$lat <- NWconv(Loc$NG,Loc$NM,Loc$NS)
Loc$long <- NWconv(Loc$WG,Loc$WM,Loc$WS)*-1

gly_Loc <-Loc[Loc$SAMPLED == "Yes" & Loc$SITUATION == "NR",] #Points sampled for the glyphosate study in NR
gly_Loc$TYPE <- as.factor(gly_Loc$TYPE)
gly_Loc$PRESENCE <- as.factor(gly_Loc$GLYPHOSATE_P)


pt_sf <- st_as_sf(gly_Loc, coords = c("long", "lat"), crs = 4326, agr = "constant") #transform df to sf

pt_sp <- gly_Loc # rename object to keep the gly_Loc as df
coordinates(pt_sp) <- ~  long + lat # transform  df to sp
proj4string(pt_sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "

# Point for general map
Landmark_sp <- pt_sp[17,] # extraer un punto cerca de la casa de Ramiro
Landmark_sf <- pt_sf[17,]

#Create labels for samples with glyphoste
gly_Loc$label <- gly_Loc$ID_PAPER
gly_Loc$label[gly_Loc$PRESENCE == "No"] <- ""


#-------------------------------------------------------------------------------
# 4. Polygon of area

library(rgdal)
library(sf)

rn_sp <- readOGR(dsn ="./Datos/Poligono_Reserva.shp") # load polygon from the area
proj4string(rn_sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "

rn_sf <- st_as_sf(rn_sp, coords = c("long", "lat"), crs = 4326)#transform to sf object

#Check plot
plot(rn_sp, axes = TRUE )

#-------------------------------------------------------------------------------
# 5. DEM of Colombian mountain ranges

library(raster)

diva_DEM_Col <- raster("./Datos/DEM_SouthAmerica/COL_alt/COL_alt.grd" ) # Descargado de diva-gis.org
Col_DEM <- mask(diva_DEM_Col,spl0_Col)

diva_DEM_df <- as.data.frame(diva_DEM_Col, xy = TRUE, na.rm = TRUE)
Col_DEM_df <- as.data.frame(Col_DEM, xy = TRUE, na.rm = TRUE)
names(diva_DEM_df) <- c("long", "lat", "Elevation")
names(Col_DEM_df) <- c("long", "lat", "Elevation")

#Check plot
plot(diva_DEM_Col)


#-------------------------------------------------------------------------------

#7. Linea de corte para perfil DEM

library(TTR) # para moving window

#Create line as sp
x_line <- c(-77.3,-67.5)
y_line <- c(5.741941,5.741941)
line_corte <- SpatialLines(list(Lines(Line(cbind(x_line,y_line)), ID="a")))
proj4string(line_corte) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "

#Convert line to sf
line_corte_sf <- st_as_sf(line_corte,coords = c("long", "lat"), crs = 4326)

#Extract elevations from raster
elevations <- as.data.frame(extract(diva_DEM_Col, line_corte, cellnumbers = TRUE, na.rm = TRUE))
xy_line <- xyFromCell(diva_DEM_Col, elevations[,1])
xy_corte <- cbind(elevations,xy_line)

#Calcular moving window average for the altitudes
movingN <- 30 # define the n for the moving average calculations
xy_corte$SMA <- SMA(xy_corte$COL_alt, 
                    n = movingN)
corte <- xy_corte[movingN:length(xy_corte$SMA), ]

#Ubicacion del punto sobre el perfil
xy_corte$punto <- ifelse(xy_corte$x == -71.41250,"yes",NA)

#-------------------------------------------------------------------------------
# 6. NDVI of the area of study

NDVI <- raster("./Datos/NDVI_850_20190406_WGS84.tif") # raster was first opened in QGIS and then saved as CRS WGS84
ext_bb <- extent(rn_sp) + 0.1 #extednded bounding box
NDVI_bb <- crop(NDVI,ext_bb)
NDVI_rn <- mask(NDVI_bb,rn_sp)

#Convert to df for ggplot
NDVI_bb_df <- as.data.frame(NDVI_bb, xy = TRUE)
names(NDVI_bb_df) <- c("long", "lat", "NDVI")

#Check plot
plot(NDVI_bb)
plot(rn_sp, add = TRUE)

#-------------------------------------------------------------------------------
# Mapping with base: just a try
#-------------------------------------------------------------------------------

library(prettymapr) #to add scalebar and arrow to baseplots

# Main map
#Graphical parameters 
myColours3 <- rep("black", 26) # Los puntos PARPOZ04 y PAO22 quedan por fuera de la zona
myColours3[c(3:5,8,9,12:14)] <- "#737373" #soil+W
myPCH1 <- rep(19,26) #Soil
myPCH1[c(17:20)]<- 17 #Deep water
myPCH1[c(21:26)] <- 2 #Surface water
#Plot

plot(rn, axes = TRUE, col = "#EDEDED", yaxt="none")
axis(4)
addscalebar(plotepsg = 4326)
addnortharrow(pos = "topleft", scale = 0.5)
plot(gly_Loc, add = TRUE, pch = myPCH1, col = myColours3, cex = 0.9) 
legend("topright", legend=c("Soil", "Soil + W","S Water","D Water"),col=c("black","#737373","black", "black"), pch=c(19,19,2,17), cex=0.7,
       title="Locations", text.font=4, bg='#EDEDED')
plot_RN <- recordPlot()
# as_gtable in case of need


# Inset map
myColours <- rep("lightgrey", 25) # graphical parameters
myColours[9] <- "darkgrey"
plot(spl1_Col, col = myColours, border = "lightgrey")


#-------------------------------------------------------------------------------
# Mapping with rasterVis: just a try 
#-------------------------------------------------------------------------------

library(rasterVis)
library(colorspace)

myTheme <- rasterTheme(region=sequential_hcl(10, palette= "Light Grays"))
p_Col_DEM <- levelplot(Col_DEM, par.settings = myTheme, scales=list(draw=FALSE), xlab = "Elevation (m.a.s.l.)", ylab = NULL) + 
  layer(sp.polygons(sp_csn, col = "#4D4D4D", fill = "#363636")) + 
  layer(sp.points(pt_Loc,pch=20, cex=0.5, col="White")) 

p_Col_DEM$legend$right <- NULL
p_Col_DEM


#-------------------------------------------------------------------------------
# Mapping with individual maps with ggplot: final versions
#-------------------------------------------------------------------------------

library(ggplot2)
library(ggsn) # to add scalebar
library(ggrepel)


#Color version of main map
main_color <- ggplot()+
  geom_raster(data = NDVI_bb_df, aes(fill = NDVI,x=long, y=lat))+
  scale_fill_gradient2(low='#C67171', mid='#FFFACD', high='#9ACD32', midpoint = 0.3)+ #olivedrab, eggshell and salmon colors
  geom_sf(data = rn_sf, fill = "white", alpha = 0.2, color = "#8B8682")+
  geom_point(data=gly_Loc, aes(colour= TYPE, size = GLYPHOSATE_C, pch= PRESENCE, x=long, y=lat))+
  geom_text_repel(data = gly_Loc, aes(label = label,  x=long, y=lat, colour = TYPE), direction = "y", nudge_x = -71.42, hjust = 1, segment.size = 0.2, segment.color = '#8B8682',size = 3, fontface = "bold")+
  annotate("text", x = -71.41, y = 5.725, label= "F1", size = 2.5, color = "#757575")+
  annotate("text", x = -71.416, y = 5.732, label= "F2", size = 2.5, color = "#757575")+
  annotate("text", x = -71.42, y = 5.737, label= "F3", size = 2.5, color = "#757575")+
  annotate("text", x = -71.418, y = 5.748, label= "F4", size = 2.5, color = "#757575")+
  annotate("text", x = -71.415, y = 5.757, label= "F5", size = 2.5, color = "#757575")+
  annotate("text", x = -71.404, y = 5.76, label= "F6", size = 2.5, color = "#757575")+
  annotate("text", x = -71.401, y = 5.753, label= "F7", size = 2.5, color = "#757575")+
  annotate("text", x = -71.399, y = 5.749, label= "F8", size = 2.5, color = "#757575")+
  annotate("text", x = -71.396, y = 5.744, label= "F9", size = 2.5, color = "#757575")+
  annotate("text", x = -71.394, y = 5.733, label= "F10", size = 2.5, color = "#757575")+
  annotate("text", x = -71.395, y = 5.72, label= "F11", size = 2.5, color = "#757575")+
  scale_color_manual(name = "", labels = c("Groundwater", "Soil", "Soil + Water", "Surface water"),values = c("#2F4F4F","#8B4513","#68838B","#33A1C9"))+
  scale_size(name = "Glyphosate [mg/L]", range = c(1.5, 4))+
  scale_shape_manual(name = "", labels = c("Undetected","Detected"), values = c(1, 19))+
  coord_sf(xlim=c(-71.44,-71.385), ylim=c(5.715, 5.77))+ #to set the limits
  theme(panel.border=element_rect(fill = NA,colour="black",size=0.8),
        axis.text=element_text(colour='black'), 
        legend.key = element_rect(fill = "white", color = NA),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        legend.position = "left",
        legend.spacing.y = unit(0.05, 'cm'),
        legend.text=element_text(size=11),
        text=element_text(size=11,  family="Times"))+
  scalebar(NDVI_bb_df, transform = TRUE, dist = 1, border.size = 0.5, dist_unit = "km", model = "WGS84", anchor = c(x = -71.42, y = 5.715), height = 0.006, st.dist = 0.008, st.size = 3)+
  guides(shape = FALSE, size = guide_legend(order = 2) ,colour = guide_legend(override.aes = list(size = 2), order = 1),fill=guide_colourbar(barheight = 3.5, order = 0))


#Color version of inset
inset_color <- ggplot()  +
  geom_raster(data = diva_DEM_df, aes(fill = Elevation, x=long, y=lat))+
  scale_fill_gradient2(low='white', mid='#FFF5EE', high='#8B8682', midpoint = 1000)+
  geom_path(data = spl0_Col, col = "#CDC5BF", size = 0.4, aes(group = group, x=long, y=lat))+
  geom_rect(data = NDVI_bb_df, aes(xmin = -71.51723, xmax = -71.29574, 
                                   ymin = 5.621294,  ymax = 5.862588), 
            color = "#B22222", fill = NA)+
  geom_sf(data = line_corte_sf, linetype = "dashed", color = "#4D4D4D")+
  coord_sf(xlim=c(-80,-67.4), ylim=c(-3.53, 13), label_axes = list(top = "E", right = "N"))+
  annotate("text", x = -72.3, y = 3.44, label= "COLOMBIA", size = 4, color = "#8B8378")+
  annotate("text", x = -74.072, y = 4.711, label= "Bogota", size = 3, color = "#4D4D4D")+
  annotate("text", x = -75.574, y = 6.2486, label= "Medellin", size = 3, color = "#4D4D4D")+
  annotate("text", x = -76.5320, y = 3.4516, label= "Cali", size = 3, color = "#4D4D4D")+
  annotate("text", x = -71.8929, y = 5.8806, label= "Paz de \n Ariporo", size = 2.5, color = "#4D4D4D")+
  geom_vline(xintercept = -71.40649, colour = "#4D4D4D", linetype = "dashed")+
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#AEEEEE",colour = "black"),
        axis.title =element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_rect(fill = NA,colour="black",size=0.8),
        text=element_text(size=11,  family="Times"))


#Corte
inset_corte <- ggplot(xy_corte,aes(x = x, y = SMA)) +
  geom_ribbon(aes(ymin = 0,  ymax = SMA),fill = "#8B8682")+
  geom_point(aes(pch = punto), col = "#B22222", fill = "#B22222", cex = 1.5)+
  scale_shape_manual(values = 25)+
  geom_vline(xintercept = -71.40649, colour = "#4D4D4D", linetype = "dashed")+
  theme(legend.position = "none",
        panel.background = element_rect(fill = NA,colour = "black"),
        axis.title =element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(colour = "#4D4D4D", size = 0.2, linetype = "dashed"),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.border=element_rect(fill = NA,colour="black",size=0.8),
        text=element_text(size=11,  family="Times"))+
        scale_y_continuous(limits = c(0,3500), position = "right")+ 
        scale_x_continuous(limits = c(-80,-67.4))

                           
#-------------------------------------------------------------------------------
# Putting maps together
#-------------------------------------------------------------------------------

#Utilizando funciones de cowplot pero sin llamar la libreria porque modifica los bordes de las graficas
Figure1 <- cowplot::ggdraw(xlim = c(0, 8.5), ylim = c(0, 4)) + 
  cowplot::draw_plot(main_color, x = 0, y = 0, width = 6, height = 4)+
  cowplot::draw_plot(inset_corte, x = 5.65, y = 0.15, width = 2.5, height = 0.75)+
  cowplot::draw_plot(inset_color, x = 5.65, y = 0.77, width = 2.5, height = 3.3)

#-------------------------------------------------------------------------------
# Save
#-------------------------------------------------------------------------------

ggsave("Figure1.pdf", width = 9, height = 4)
