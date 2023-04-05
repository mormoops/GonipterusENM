### - Gonipterus distribution models - ###
### - script by: Angelo Soto-Centeno - ###
### - revised date: 5 pril 2023 - ###

# load packages
library(sp) # mapping
library(raster) # raster management
library(maptools) # mapping
library(rgdal) # raster management
library(dismo) # sp data mining, gis
library(ecospat) # niche analysis
library(tidyverse) # data wrangling
library(maps) # mapping
library(rJava) # java controls
library(ENMeval) # enm tuning
library(rasterVis) # raster management
library(RColorBrewer) # color palettes
library(rgeos)
library(mapdata) # maps

######################################################
######## 1. ### import & wrangle climate data ########
######################################################
raster_files <- list.files("/Volumes/ANGELO4/ASC_GIS/Layers/Climate/Climate/Present/wc2-2.5m/", full.names = T, pattern = ".tif")
# create a raster stack using the list you created
predictors <- stack(raster_files)
# inspect full clim data
plot(predictors$bio_1)

## create geographic extent of region of interest
# geo.extEC for Ecuador
# geo.extSA for South America
geo.extEC <- extent(x = c(-81, -75, -5, 1.5))
geo.extSA <- extent(x = c(-81, -34, -56, 12))

# trim clim data by region/country 
# 1. create predictors for the Calibration extent
envSA <- crop(x = predictors, y = geo.extSA)
# inspect trimmed clim data
plot(envSA$bio_1)


#################################################################
######## 2. ### mine Gonipterus platensis data from GBIF ########
#################################################################
# wrangled species observation data included in Dryad (doi:XXXXXXX)
goplSA.3a <- read.csv("~/Documents/Working_On/WorkingPubs/GonipterusBeetles_Pinto_etal/SpatialData/NEW_data/goplSA.3a.csv", header = T)

# plot to inspect
plot(envSA$bio_1, col = viridis::viridis(99), main = "G. platensis localities", xlab = "longitude", ylab = "latitude")
points(goplSA.3a[, 2:3], col = "black")


##############################################
######## 3. crate calibration extents ########
##############################################

# Make our occs into a sf object
occs.gopl.sf <- sf::st_as_sf(goplSA.3a[, 2:3], coords = c("lon", "lat"), crs = raster::crs(envSA)) 

# Re-project the point data
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.gopl.sf <- sf::st_transform(occs.gopl.sf, crs = eckertIV)

# Buffer all occurrences by 500 km
occs.buf <- sf::st_buffer(occs.gopl.sf, dist = 500000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envSA))
plot(envSA$bio_1, main = names(envSA)[1])
points(goplSA.3a[, 2:3])
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)

# Crop environmental rasters to match the study extent
env.bg.gopl <- raster::crop(envSA, occs.buf)
# Mask the rasters to the shape of the buffers
env.bg.gopl <- raster::mask(env.bg.gopl, occs.buf)
plot(env.bg.gopl[[1]], main = names(envSA)[1])
points(goplSA.3a[, 2:3])
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)

# Randomly sample 10,000 background points from one background extent raster 
bg <- dismo::randomPoints(env.bg.gopl[[9]], n = 10000) %>% as.data.frame()
colnames(bg) <- colnames(goplSA.3a[, 2:3])

# Inspect cell coverage
plot(env.bg.gopl[[1]])
points(bg, pch = 20, cex = 0.2)


##########################################
######## 4. ENMeval model testing ########
##########################################

# set the parameter tunings
tune.args <- list(fc = c("L", "LQ", "H", "LQH", "LQHP"), rm = 1:5)
# G platensis
eval.gopl <- ENMevaluate(occs = goplSA.3a[, 2:3], envs = envSA, bg = bg, 
                         algorithm = 'maxnet', partitions = 'block', 
                         tune.args = tune.args)
# save the results
write.csv(eval.gopl@results, file = "enmeval.res.csv", row.names = F)
# best tuning parameters found in the *.csv file
# LQH rm = 1

#################################################################
######## 5.  Default ENMs with Maxent                    ########
#################################################################

# create a Maxent model with default parameters
mx.dflt <- maxent(env.bg.gopl, goplSA.3a[ , 2:3] , path = "./gopl.dflt")

# create a raster of your Maxent model prediction
mx.dflt.dist <- predict(mx.dflt, envSA, args = c("outputformat=cloglog"), progress = 'text')
# check variable contribution
plot(mx.dflt, main = "G platensis: variable contribution, default model")
# make response curves
response(mx.dflt, main = "G platensis: response curves, default model")

# test the prediction accuracy 
ecospat.boyce(mx.dflt.dist, goplSA.3a[ , 2:3], window.w = "default", res = 100, PEplot = T)
# statistic is the Spearman.cor value (= 0.882)

# examine the predicted distribution
plot(mx.dflt.dist, main = "G platensis (Default - Train AUC 0.966, BI = 0.882)", 
                    xlab = "longitude", ylab = "latitude")
points(x = goplSA.3a$lon, y = goplSA.3a$lat, col = "black", pch = 1, cex = 0.75)

# write the modeled distribution raster
writeRaster(mx.dflt.dist, filename = "./mx.dflt.dist.tif")


##################################################################
######## 6. improved Maxent model using custom parameters ########
##################################################################

## create Maxent model using ALL species points and best parameters: fc.LQHP_rm.3
mx.cust <- maxent(env.bg.gopl, goplSA.3a[ , 2:3], args = c("linear=true", "quadratic=true", 
                                                          "product=true", "hinge=true", "threshold=false", 
                                                          "betamultiplier=3"), path = "./Gonipterus3")

# check variable contribution
plot(mx.cust, main = "G platensis: variable contribution, custom model")
# check response curves
response(mx.cust, main = "G plat: var cont, custom")

# make present model prediction
mx.cust.dist <- predict(mx.cust, envSA, args = c("outputformat=cloglog"), progress = "text")

# calculate Boyce index
ecospat.boyce(mx.cust.dist, goplSA.3a[ , 2:3], window.w = "default", res = 100, PEplot = T)

# examine the predicted distribution
mx.cust.mod <- plot(mx.cust.dist, main = "G platensis (Custom - Train AUC 0.959, BI = 0.862)", 
                    xlab = "longitude", ylab = "latitude")
points(x = goplSA.3a$lon, y = goplSA.3a$lat, col = "black", pch = 1, cex = 0.75)

# write model distribution raster
writeRaster(mx.cust.dist, filename = "Gonipterus3/mx.cust.dist.tif")


#######################################################################
######## 7. Characterization of climate envelopes              ########
#######################################################################

## create a boxplot of Australia (native) vs S America (invasive) ranges for important environmental variables
# make a colors vector
clrs = c("#8c510a", "#35978f")

# create the boxplot for bio1
bio1vp <- ggplot(goplSA.3a, aes(x = colabel, y = bio_1)) +
  geom_boxplot(outlier.alpha = 1, coef = 1, color = "#525252", fill = clrs, width = 0.5, alpha = 0.75) +
  geom_jitter(color = "#525252", size = 2, alpha = 0.4) +
  labs(x = "Locality", y = "Temperature seasonality") +
  theme_light() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12))

## NOTE: repeat boxplot code changing each of the climate variables, for example:
##       y = bio_4 ... etc.

## statistics ##
# inspect country labels
unique(goplSA.3a$colabel)

# Shapiro test of normality per variable. the two non-independent groups are Ecuador and S America
with(goplSA.3a, shapiro.test(bio_4[colabel == "South America"])) # W = 0.9065, p-value = 0.01406 **NOT normal**
with(goplSA.3a, shapiro.test(bio_4[colabel == "Ecuador"])) # W = 0.71628, p-value = 2.293e-05 **NOT normal**

with(goplSA.3a, shapiro.test(bio_1[colabel == "South America"])) # W = 0.95074, p-value = 0.1913 **YES normal**
with(goplSA.3a, shapiro.test(bio_1[colabel == "Ecuador"])) # W = 0.89871, p-value = 0.02378 **NOT normal**

with(goplSA.3a, shapiro.test(bio_19[colabel == "South America"])) # W = 0.91293, p-value = 0.02026 **NOT normal**
with(goplSA.3a, shapiro.test(bio_19[colabel == "Ecuador"])) # W = 0.8665, p-value = 0.005499 **NOT normal**

with(goplSA.3a, shapiro.test(bio_6[colabel == "South America"])) # W = 0.87658, p-value = 0.002818 **NOT normal**
with(goplSA.3a, shapiro.test(bio_6[colabel == "Ecuador"])) # W = 0.9588, p-value = 0.4395 **YES normal**

with(goplSA.3a, shapiro.test(elev[colabel == "South America"])) # W = 0.5912, p-value = 8.056e-08 **NOT normal**
with(goplSA.3a, shapiro.test(elev[colabel == "Ecuador"])) # W = 0.963, p-value = 0.5266 **YES normal**


# Wilcoxon ranked sum test (non-parametric)
wilcox.test(bio_4 ~ colabel, data = goplSA.3a, exact = F) # W = 9, p-value = 2.38e-09 ** REJECT 
wilcox.test(bio_1 ~ colabel, data = goplSA.3a, exact = F) # W = 366, p-value = 0.5555 ** NOT REJECT
wilcox.test(bio_19 ~ colabel, data = goplSA.3a, exact = F) # W = 228, p-value = 0.05303 ** NOT REJECT
wilcox.test(bio_6 ~ colabel, data = goplSA.3a, exact = F) # W = 543, p-value = 0.0001178 ** REJECT
wilcox.test(elev ~ colabel, data = goplSA.3a, exact = F) # W = 626, p-value = 7.449e-08 ** REJECT


## correct p-values for the 5 comparisons
# create a vector of your p-value results from statistical tests
p.val <- c(0.00000000238, 0.5555, 0.05303, 0.0001178, 0.00000007449)

#  use p.adjust to adjust the p-values based on a specific correction method (e.g. Bonferoni). see ?p.adjust
p.adjust(p.val, method = "bonferroni", n = 5)
## this produces the following list of corrected p-values based on Bonferroni correction
# [1] 1.1900e-08 1.0000e+00 2.6515e-01 5.8900e-04 3.7245e-07
## written in non sci form
# 0.0000000119 - BIO4 = Temperature Seasonality (standard deviation ×100)
# 1.0 - BIO1 = Annual Mean Temperature 
# 0.26515 - BIO19 = Precipitation of Coldest Quarter
# 0.000589 - BIO6 = Mean Temperature of Coldest Month
# 0.00000037245 - ELEV = Elevation



## ¡muerto el pollo! ##

# clear environment & plotting dev
rm(list = ls())
dev.off()