# Code for landscape biomass map

### Install required packages ----
library(raster)
library(sf)
library(rgdal)

### Load data ----
# Load rasters
CHM <- raster("data/site_rasters/CHM.tif")                                      # Read in canopy height model at 0.010 m grain.
NDVI_018 <- raster("data/site_rasters/NDVI_018.tif")                            # Read in NDVI raster at 0.018 m grain.
NDVI_047 <- raster("data/site_rasters/NDVI_047.tif")                            # Read in NDVI raster at 0.047 m grain.
NDVI_119 <- raster("data/site_rasters/NDVI_119.tif")                            # Read in NDVI raster at 0.119 m grain.
NDVI_121 <- raster("data/site_rasters/NDVI_121.tif")                            # Read in NDVI raster at 0.121 m grain.
RGB <- stack("data/site_rasters/RGB_040.tif")                                   # Read in multi-band RGB orthmosaic at 0.040 m grain.

# Load polygon
AOI <- st_read("data/site_rasters/monitoring_plot.geojson", crs = 32607)        # Import polygon of monitoirng plot as sf object, using st_read to allow the non-standard CRS to be specified.


# Clip rasters to the monitoring extent
AOI_CHM <- crop(CHM, AOI)
AOI_NDVI_018 <- crop(NDVI_018, AOI)
AOI_NDVI_047 <- crop(NDVI_047, AOI)
AOI_NDVI_119 <- crop(NDVI_119, AOI)
AOI_NDVI_121 <- crop(NDVI_121, AOI)
AOI_RGB <- crop(RGB, AOI)




# compute mean canopy height
meanCH <- raster::extract(AOI_CHM, AOI, fun=mean)


# plot results
plot(AOI_NDVI_121)

test <- plot(AOI_NDVI_121)

par("mar")


exact_extract(CHM, AOI, 'mean')



# # view histogram of data
hist(AOI_NDVI_121,
     main="Distribution of CHM Values",
     xlab="CHM Elevation Value (m)",
     ylab="Frequency",
     col="wheat")

# Set common limits to sclars for plotting?
plot(CHM, 
     main="canopy height model")


#   # Calculate mean caopy height
#   plots$mean_NDVI_018 <- exact_extract(CHM, AOI, 'mean')
#  
#
#   # Tidy dataframe
#   plots_df <- st_drop_geometry(plots)                                           # Create dataframe from simple features object (dropping geometry).
#   plots_df$EPSG <- NULL                                                         # Remove unecessary EPSG column from dataframe.
#   plots_df <- plots_df[order(plots_df$PlotID),]                                 # Order dataframe by PlotID.
#
#   # Export NDVI values
#   write.csv(plots_df,"data/Extracted_NDVI.csv", row.names = FALSE)            # extracted NDVI values were added to the main_database file. ndvi_data <- read.csv("data/Extracted_NDVI.csv", header = T)                  # Read in NDVI values from Exact Extract pipeline.




### Apply  models to estiamte biomass


### Calcualte biomass difference rasters
https://www.neonscience.org/dc-raster-calculations-r 

### Tabluate biomass estimates




### Plot rasters