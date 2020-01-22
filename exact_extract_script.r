# Raster analysis with exactextractr 

# Establish operating environment ----
home <- "C:/workspace/OrcaManuscript/"

### Load required packages ----
library(geojsonsf)                                                              # Converts GeoJSON to sf object.
library(sf)                                                                     # Package for handeling simple features.
library(raster)                                                                 # 
library(exactextractr)                                                          # Efficient and exact extraction of raster statistics. 
library(rgeos)                                                                  # for validating which hole belongs to which exterior ring


### Load data ----
feature_filename <- paste0(home, "data/20160725_AC_ORC - formated for exact extractr.geojson")
plots <- st_read(feature_filename, crs = 32607)                                 # Import geoJSON as sf object, using st_read to allow the non-standerd CRS to be specified.


# Import rasters (Available from NERC Polar Data Centre - see readme) ----
raster_018 <- raster(paste0(home, "inputs_NDVI/NDVI_019m_from_20160726.tif")) # import plot shapefiles.
raster_047 <- raster(paste0(home, "inputs_NDVI/NDVI_050m_from_20160730.tif")) # import plot shapefiles.
raster_119 <- raster(paste0(home, "inputs_NDVI/NDVI_120m_from_20160730.tif")) # import plot shapefiles.
raster_121 <- raster(paste0(home, "inputs_NDVI/NDVI_120m_from_20160726.tif")) # import plot shapefiles.


# Calculate mean NDVI for each polygon ----
plots$mean_NDVI_018 <- exact_extract(raster_018, plots, 'mean')                        
plots$mean_NDVI_047 <- exact_extract(raster_047, plots, 'mean')                        
plots$mean_NDVI_119 <- exact_extract(raster_119, plots, 'mean')                        
plots$mean_NDVI_121 <- exact_extract(raster_121, plots, 'mean')       

# 
plots_df <- st_drop_geometry(plots)                                             # Create dataframe from simple features object (dropping geometry).
plots_df$EPSG <- NULL                                                           # Remove unecessary EPSG column from dataframe.
plots_df <- plots_df[order(plots_df$PlotID),]                                   # Order dataframe by PlotID.

# Export NDVI values
write.csv(plots_df,"data/Extracted_NDVI.csv", row.names = FALSE)