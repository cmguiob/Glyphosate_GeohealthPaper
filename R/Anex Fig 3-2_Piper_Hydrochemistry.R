#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 14-09-2019
# Output: PIPER PLOTS OF SAMPLES WITH WATER CHEMISTRY ANALYSES
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script  plots piper diagrams of hydrochemical data from Paz de Ariporo
# using two packages for exploratory analysis
#-------------------------------------------------------------------------------
# EXPLORATORY ANALYSIS WITH "Hydrogeo"
#-------------------------------------------------------------------------------

library(hydrogeo)

water_PA <- read.csv("./Datos/Agua_TCI_PA.csv")
select_PA <- water_PA[,c(18, 19, 20, 21, 24, 25, 26, 27)]
perc_PA <- toPercent(select_PA) 

piper_PA <- piper(perc_PA)
piper_PA@pt.pch = c(6,6,6,2,15,2,2,15,15,1,1,19,19,19,1,19) 
# Circles are pits, squares streams,  upward triangles wetlands, downward triangles wells. Filled are samples that showd unusual values
piper_PA@pt.col = c(1,1,1,1,1,1,1,1,1,1,1,8,1,1,1,8) # Grey are samples with unusual relaiton on the diagram.
plot(main = "Diagrama Piper Aves d'Jah", piper_PA, cex = 1, cex.axis = 1.5)

# *****************COMENTARIOS*************************
# Hay una muestra de calicata con un contenido relativo de cloruro muy alto 

#-------------------------------------------------------------------------------
# EXPLORATORY ANALYSIS WITH "smwrGraphs" (from the USGS repository)
#-------------------------------------------------------------------------------
# ***********TO INSTALL**************************
# 1. check if USGS repository is added to the Rprofile using setRepositories()
# 2. If not, add it:

#rprofile_path = file.path(Sys.getenv("HOME"), ".Rprofile")
#write('\noptions(repos=c(getOption(\'repos\'),
#      CRAN=\'https://cloud.r-project.org\',
#      USGS=\'https://owi.usgs.gov/R\'))\n',
#      rprofile_path, 
#      append =  TRUE)

# 3. Restart R session
# 4. install.packages("smwrGraphs")
# 5. install.packages("smwrData") # for data examples

# ***********LOAD**************************

library(smwrGraphs)

water_PA <- read.csv("./Datos/Agua_TCI_PA.csv")
select1_PA <- water_PA[,c(10,18, 19, 20, 21, 23,24, 25, 26, 27, 28, 31, 32)]

# ***********PREPARE VARIABLES **************************

# Transform to milliequivalents
meq_PA <- transform(select1_PA, Ca.meq = conc2meq(Ca, "calcium"),
                    K.meq = conc2meq(K, "potassium"),
                    Mg.meq = conc2meq(Mg, "magnesium"),
                    Na.meq = conc2meq(Na, "sodium"),
                    NH4.meq = conc2meq(NH4, "ammonium"),
                    Cl.meq = conc2meq(Cl, "chloride"),
                    SO4.meq = conc2meq(SO4, "sulfate"),
                    HCO3.meq = conc2meq(HCO3, "bicar"))


# Own calculation of milliequivalents to check
meq_PA$Na_omeq <- (meq_PA$Na/22.989769)*1 #not used
meq_PA$PO4.meq <- (meq_PA$PO4/94.9714)*3 #not used

# ***********PLOTS **************************

# Set graphic convetions
meq_PA$SYMBOL <- c("downtri","downtri","downtri","uptri","square","uptri","uptri","square","square","circle","circle","circle","circle","circle","circle","circle")
meq_PA$COL <- c("grey50","grey50","grey50","grey50","black","grey50","grey50","black","black","grey50","grey50","grey50","black","black","grey50","grey50")

# Plot with Sodium
setSweave("Anex 1_Piper_Na", 7, 7) # Start plot

AA.pl <- with(meq_PA, piperPlot(Ca.meq, Mg.meq, Na.meq, Cl.meq, HCO3.meq, SO4.meq,
                                Plot=list(name=RESUL, symbol = SYMBOL, color = COL),
                                xCat.title = "Calcium",
                                yCat.title = "Magnesium",
                                zCat.title = "Sodium",
                                xAn.title = "Chloride",
                                yAn.title = "Bicarbonate",
                                zAn.title = "Sulfate"))

addExplanation(AA.pl, where="ul", title="")
graphics.off() #End plot

# Plot with Potassium
setSweave("Anex 1_Piper_K", 7, 7) # Start plot

AA.pl <- with(meq_PA, piperPlot(Ca.meq, Mg.meq, Na.meq, Cl.meq, HCO3.meq, SO4.meq,
                                Plot=list(name=RESUL, symbol = SYMBOL, color = COL),
                                xCat.title = "Calcium",
                                yCat.title = "Magnesium",
                                zCat.title = "Potassium",
                                xAn.title = "Chloride",
                                yAn.title = "Bicarbonate",
                                zAn.title = "Sulfate"))

addExplanation(AA.pl, where="ul", title="")
graphics.off() #End plot


# Plot with Potassium and size proportional to glyphosate
PA.size <- 0.04+log(meq_PA$Glypho+1)/5

setSweave("Piper_K_GlyphoSized", 7, 7) # Start plot
AA.pl <- with(meq_PA, piperPlot(Ca.meq, Mg.meq, Na.meq, Cl.meq, HCO3.meq, SO4.meq,
                                Plot=list(what = "none"),
                                ticks = TRUE,
                                xCat.title = "Calcium",
                                yCat.title = "Magnesium",
                                zCat.title = "Potassium",
                                xAn.title = "Chloride",
                                yAn.title = "Bicarbonate",
                                zAn.title = "Sulfate"))

with(AA.pl, addPiper(xCat=cations$x, yCat=cations$y, xAn=anions$x, yAn=anions$y,
                     xPip=piper$x, yPip=piper$y,
                     Plot=list(size=PA.size, filled=FALSE), current=AA.pl))

addExplanation(AA.pl, where="ul", title="")
graphics.off() # End plot

# Plot with Ammonium

# Add the bases according to valence
meq_PA$Na_K.meq <- meq_PA$Na.meq  + meq_PA$K.meq 
meq_PA$Ca_Mg.meq <- meq_PA$Ca.meq  + meq_PA$Mg.meq 

meq_PA$SYMBOL <- c("downtri","downtri","downtri","uptri","square","uptri","uptri","square","square","circle","circle","circle","circle","circle","circle","circle")
meq_PA$COL <- c("grey50","grey50","grey50","grey50","black","grey50","grey50","black","black","grey50","grey50","grey50","black","black","grey50","grey50")

setSweave("Anex 1_Piper_Ammonium", 7, 7)

AA.pl <- with(meq_PA, piperPlot(Ca.meq, Mg.meq, NH4.meq, Cl.meq, HCO3.meq, SO4.meq,
                                Plot=list(name=RESUL, symbol = SYMBOL, color = COL),
                                xCat.title = "Calcium",
                                yCat.title = "Magnesium",
                                zCat.title = "Ammonium",
                                xAn.title = "Chloride",
                                yAn.title = "Bicarbonate",
                                zAn.title = "Sulphate"))

addExplanation(AA.pl, where="ul", title="")
graphics.off() # End plot

