#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 25-08-2020
# Output: Hydrochemical analysis with PCA
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script performs a screening on a hydrochemical dataset of major ions plots
# and plots a biplot of PCA for publication.
#-------------------------------------------------------------------------------
# RESOURCES:
# Log transform adding 1: https://www.youtube.com/watch?v=aaN_atlM1kM
# Data quality assurance following Güler et al. 2002
# Format Conversions: http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/

#-------------------------------------------------------------------------------
# OBJECTIVES:
# 1. Load and confirm data 
# 2. Perform quality assurance: location check and charge balance error.
# 3. PLot with ggplot2

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------

WPA <- read.csv("./Datos/Agua_TCI_PA.csv")

#-------------------------------------------------------------------------------
# QUALITY ASSURANCE
#-------------------------------------------------------------------------------

#************** Spatial check ************************************
library(sp)

# Check that samples locations make sense and there are no missing data

# reate lat-long data
NWconv <- function(i,j,k){i+j/60+k/3600} # create new variable with coordinates
WPA$lat <- NWconv(WPA$NG,WPA$NM,WPA$NS)
WPA$long <- NWconv(WPA$WG,WPA$WM,WPA$WS)*-1

# Missing locations
summary(is.na(WPA$lat)) 
summary(is.na(WPA$long)) 

#initialize coordinates
WPAsp <- WPA
coordinates(WPAsp) <- ~ long + lat
coordinates(WPAsp)

# set spatial reference system
proj4string(WPAsp) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
proj4string(WPAsp) 
plot(WPAsp)

#************** Charge balance error  ************************************

library(smwrGraphs)
WPA_ions <- WPA[,c(20:29)]

meq_PA <- transform(WPA_ions, Ca.meq = conc2meq(Ca, "calcium"),
                    K.meq = conc2meq(K, "potassium"),
                    Mg.meq = conc2meq(Mg, "magnesium"),
                    Na.meq = conc2meq(Na, "sodium"),
                    NH4.meq = conc2meq(NH4, "ammonium"),
                    Cl.meq = conc2meq(Cl, "chloride"),
                    SO4.meq = conc2meq(SO4, "sulfate"),
                    HCO3.meq = conc2meq(HCO3, "bicarbonate"),
                    PO4.meq = conc2meq(PO4, "phosphorus as p"),
                    NO3.meq = conc2meq(NO3,"nitrate as n"))

meq_PA$CatSum <- meq_PA$Ca.meq + meq_PA$K.meq + meq_PA$Na.meq + meq_PA$Mg.meq + meq_PA$NH4.meq
meq_PA$AniSum <- meq_PA$Cl.meq + meq_PA$SO4.meq + meq_PA$NO3.meq + meq_PA$PO4.meq + meq_PA$HCO3.meq

meq_PA$CB <- (meq_PA$CatSum - meq_PA$AniSum)/(meq_PA$CatSum + meq_PA$AniSum)
plot(meq_PA$CB)
abline(h = 0.1, col = "red")
abline(h=-0.1, col = "red")

# CONCLUSION: Some important cations were not analysed: Fe2+?

#-------------------------------------------------------------------------------
# DATA SCREENING
#-------------------------------------------------------------------------------

#*********************************** Outliers ******************************************
library(dplyr)
library(ggrepel)

is_above75 <- function(x) {
  return(x > quantile(x, 0.75) + 1 * IQR(x))
}

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

WPA_outliers <- group_by(WPA_mstat_long, variable) %>% mutate(outlier = ifelse(is_outlier(value), as.character(ID_PAPER), as.character(NA)))

WPA_above75 <- group_by(WPA_mstat_long, variable) %>% mutate(above75 = ifelse(is_above75(value), as.character(ID_PAPER), as.character(NA)))


ggplot(WPA_above75, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill =  variable), lwd = 0.3, outlier.size = 1.5) +
  geom_text_repel(aes(label = above75), na.rm = TRUE, size = 2.5, nudge_y = 0.01, nudge_x = 0.15, force = 0.5, segment.color = "grey")+
  facet_wrap(~variable, nrow = 2, scales = "free_x")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = "Major ions", y = "ppm")+
  theme(strip.background = element_blank())

ggplot(WPA_outliers, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill =  variable), lwd = 0.3, outlier.size = 1.5) +
  geom_text_repel(aes(label = outlier), na.rm = TRUE, size = 2.5, nudge_y = 0.01, nudge_x = 0.15, force = 0.5, segment.color = "grey")+
  facet_wrap(~variable, nrow = 2, scales = "free_x")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = "Major ions", y = "ppm")+
  theme(strip.background = element_blank())

#******************************* Check bias  ************************************

library(reshape2)
library(ggplot2)

#Change to long form for plotting in facets
WPA_mstat <- WPA[,c(4,20:29, 33)]
WPA_mstat_long <- melt(WPA_mstat, id.vars = c("ID_PAPER", "RESUL"), measure.vars = c("HCO3","Cl","SO4","NO3", "PO4","Ca","K","Mg","Na","NH4"))
WPA_mstat_long$variable <- as.factor(WPA_mstat_long$variable)

# Histograms 
ggplot(WPA_mstat_long, aes(x = value)) +
  geom_histogram(aes(fill =  variable)) +
  facet_wrap(~variable, nrow = 2)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_flip()+
  labs(y = "Counts", x = "ppm")

#Log transform due to bias
WPA_mstat_long$logT <- log(WPA_mstat_long$value + 1) # To avoid the problem with zeros

# HIstogram after tranform
ggplot(WPA_mstat_long, aes(x = logT)) +
  geom_histogram(aes(fill =  variable)) +
  facet_wrap(~variable, nrow = 2)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_flip()+
  labs(y = "Counts", x = "ppm")


# Tansform back to wide format
WPA_mstat_long$value <- NULL
WPA_mstat_logT <- dcast(WPA_mstat_long, formula = ID_PAPER + RESUL ~ variable , value.var = "logT")

# CONCLUSION: Data were positively skewed: data contained small numbers of hugh values

#-------------------------------------------------------------------------------
#PCA ANALYSIS
# How was PCA used in papers?
# Should correlated variables be removed?
# What is the criteria to select the input variables?
# Should variables of different units be used as input? e.g. pH vs. [PO4]
# Should data be transformed beyond scaling?
# How to interpret the direction of the arrows? and the grouping of samples around the arrows?
# Which Principal COmponents should be used? Can any combiantion be used?
# What is the minimum explained variance to be accepted on a biplot?
# What is the meaning of the different outputs from the analysis?
#-------------------------------------------------------------------------------

# ************Fast exploratory analysis ***********

library(ggbiplot)

# pca with non-transformed data
WPA_pca <- prcomp(WPA_mstat[,2:11], center = TRUE,scale. = TRUE) # centered and scaled
ggbiplot(WPA_pca, labels=WPA_mstat$ID_PAPER, groups=WPA_mstat$RESUL,choices=c(1,2),obs.scale = 1, var.scale = 1, ellipse = TRUE)

WPA_logT_pca <- prcomp(WPA_mstat_logT[,3:12], center = TRUE,scale. = TRUE)
ggbiplot(WPA_logT_pca, labels=WPA_mstat_logT$ID_PAPER, groups=WPA_mstat_logT$RESUL,choices=c(1,2),obs.scale = 1, var.scale = 1, ellipse = TRUE)

# Conclusion: log transformed data show a better partitioning than not transformed. 
# P1-P2: Samples with glyphosate exhibit higher contents of anions
# P2-P3: S5 and S13 are more related to each other due to relatively higher NO3 and NH4 contents,
# while S9 and S12 exhibit the highest SO4 contents and low N. 

# ************ Plot for paper ***********

# Label in red the ones that have glyphosate detected and in gray the others.
# Add group circle with alpha
# Color lines according to contribution to PCAs?

library(factoextra)
library(FactoMineR)
library(wesanderson)

rownames(WPA_mstat_logT) <- WPA_mstat_logT$ID_PAPER #change this to get labeled data on biplot
WPA.pca <- PCA(WPA_mstat_logT[,3:12], scale.unit = TRUE, graph = FALSE)
mypal <- wes_palette("Zissou1", 2, type = "continuous")
# Save as png
tiff(filename="Figure 3a.tiff", 
    type="cairo",
    units="in", 
    width=2.8, 
    height=3, 
    pointsize='1', 
    res=400)

fviz_pca_biplot(WPA.pca,
                fill.ind = WPA_mstat_logT$RESUL, col.ind = "#848484",
                labelsize = 2, 
                arrowsize = 0.4,
                pointshape = 21, pointsize = 1.7,
                addEllipses = TRUE, ellipse.level = 0.6,ellipse.type = "norm", invisible="quali", #remove centroid
                col.var = "contrib",
                gradient.cols = c("#ffa600", "#bc5090"),
                legend.title = list(fill = "Glyphosate", color = "Contribution"))+  
  scale_fill_manual(values=c("#848484", "#FC4E07"))+
  theme_bw(base_size = 8)+
  guides(color = FALSE, fill = FALSE)+
  theme(legend.position="top",
        legend.key = element_rect(fill = NA, color = NA),
        legend.text = element_text(size = 7.5),
        legend.margin=margin(0,0,0,0),
        axis.text = element_text(size= 7),
        axis.title= element_text(size = 7.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("")+ xlab("PC1 (36.4%)") + ylab("PC2 (23.6%)")

dev.off() # Finish saving plot









