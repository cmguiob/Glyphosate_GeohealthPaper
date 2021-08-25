#-------------------------------------------------------------------------------
# Author: Carlos Guio
# Contact: macguiob@gmail.com 
# Organization: Terrae
# Project: TCI-Paz de Ariporo-Glifosato
# Date last modified: 24-08-2020
# Output: CALIBRATION CURVE
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script  plots a calibration curve showing calibration points with
# error bars and the 95% confidence interval for the regression line.
#-------------------------------------------------------------------------------
# RESOURCES:
# https://aosmith.rbind.io/2018/11/16/plot-fitted-lines/
# https://plotly.com/ggplot2/geom_ribbon/
# https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
# https://stackoverflow.com/questions/29554796/meaning-of-band-width-in-ggplot-geom-smooth-lm

#-------------------------------------------------------------------------------
# OBJECTIVES:
# 1. Load and confirm data of means and standard deviations
# 2. Calculate regression line
# 3. PLot with ggplot2

#-------------------------------------------------------------------------------
# LOAD AND CONFIRM
#-------------------------------------------------------------------------------

cali <- read.csv("./Datos/Glifo_cali.csv")

library(dplyr)
grouped <- group_by(cali, ID)
summarise(grouped, mean=mean(Abs), sd=sd(Abs))

#-------------------------------------------------------------------------------
#CALCULATE REGRESSION 
#-------------------------------------------------------------------------------

fit <- lm(Abs ~ Glifo, data = cali)
summary(fit)

#-------------------------------------------------------------------------------
# PLOT AND SAVE
#-------------------------------------------------------------------------------

library(ggplot2)

# Labels to add in the plot
labform <- as.character(expression("y == 0.7514 + 0.2465*x"))
labform2 <- as.character(expression("R^2 == 0.9732"))

# Save as png
png(filename="Figure 3.png", 
    type="cairo",
    units="in", 
    width=2.5, 
    height=2.5, 
    pointsize=1, 
    res=300)

# Plot in device
ggplot(cali, aes(x = Glifo, y = Abs) ) +
  geom_smooth(method = "lm", lwd = 0.4, col = "#bc5090", fill = "#CDC9C9")+
  geom_errorbar(aes(ymin=Mean-Stdev, ymax=Mean+Stdev), width=.03, lwd = 0.25, col = "#8B8989")+
  geom_point(aes(x = Glifo, y = Mean ),col = "#8B8989", size = 0.35)+
  geom_vline(xintercept = 0.25, linetype="dashed", 
             color = "#CDC9C9", size=0.5)+
  theme_bw(base_size = 7.5)+
  labs(x= "Glyphosate (ppm)", y = "Absorbance 265 (nm)")+
  annotate(geom="text", x=2.5, y=0.25, label= labform, color="darkgray", hjust = 1,size = 2.5, face = "bold", parse= TRUE)+
  annotate(geom="text", x=2.5, y=0.1, label= labform2, color="darkgray", hjust = 1, size = 2.5, face = "bold", parse= TRUE)+
  theme(axis.text = element_text(size= 5.5),
        axis.title= element_text(size = 6.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  

dev.off() # Finish saving plot
