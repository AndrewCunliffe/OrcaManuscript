# Raster analysis with exactextractr 
# https://github.com/isciences/exactextractr 

# Establish operating environment ----
home <- "C:/workspace/Geospatial_Pipeline/Geospatial_Pipeline_OrcaProject/"

#### Load required packages ####
# install.packages("rgeos")
library(geojsonsf)                                                              # Converts GeoJSON to sf object.
library(sf)                                                                     # Package for handeling simple features.
library(raster)                                                                 # 
library(exactextractr)                                                          # Efficient and exact extraction of raster statistics. 
library(rgeos) # for validating which hole belongs to which exterior ring



# # Example from https://github.com/isciences/exactextractr
# brazil <- st_as_sf(getData('GADM', country='BRA', level=2))
# brazil2 <- brazil
# prec <- getData('worldclim', var='prec', res=10)[[12]]
# brazil2$mean_prec <- exact_extract(prec, brazil, 'mean')
# head(brazil)


#### Load data ----
feature_filename <- paste0(home, "outputs_04_geoJSON/20160725_AC_ORC - formated for exact extractr.geojson")
plots <- st_read(feature_filename, crs = 32607) # Import geoJSON as sf object. NB. st_read rather than geojson_sf is necesasary to allow the non-standerd CRS to be specified.
# head(plots)

# Import rasters
raster_018 <- raster(paste0(home, "inputs_4_NDVI/20160725_AC_ORC_NDVI_019m_from_20160726.tif")) # import plot shapefiles.
raster_047 <- raster(paste0(home, "inputs_4_NDVI/20160725_AC_ORC_NDVI_050m_from_20160730.tif")) # import plot shapefiles.
raster_119 <- raster(paste0(home, "inputs_4_NDVI/20160725_AC_ORC_NDVI_120m_from_20160730.tif")) # import plot shapefiles.
raster_121 <- raster(paste0(home, "inputs_4_NDVI/20160725_AC_ORC_NDVI_120m_from_20160726.tif")) # import plot shapefiles.
# raster_018        # Summarise raster
# summary(raster)   # Summarise a sample of raster values.


# Calculate mean NDVI for each polygon.
# NB: 'exact_extract' does not like Z values in the features ("ParseException: Unknown WKB type 235"!).
# NB. polygons must be CLOSED linear rings.
plots$mean_NDVI_018 <- exact_extract(raster_018, plots, 'mean')                        
plots$mean_NDVI_047 <- exact_extract(raster_047, plots, 'mean')                        
plots$mean_NDVI_119 <- exact_extract(raster_119, plots, 'mean')                        
plots$mean_NDVI_121 <- exact_extract(raster_121, plots, 'mean')       

# Example syntax
# plots[, c('min_ndvi', 'max_ndvi')] <- exact_extract(raster, plots, c('min', 'max')) # Find min and max values in a single pass

plots_df <- st_drop_geometry(plots)                                             # Create dataframe from simple features object (dropping geometry).
plots_df$EPSG <- NULL                                                           # Remove unecessary EPSG column from dataframe.
plots_df <- plots_df[order(plots_df$PlotID),]                                   # Order dataframe by PlotID.

write.csv(plots_df,"outputs_01_main/Extracted NDVI.csv", row.names = FALSE)