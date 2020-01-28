# Clear working envrionmnt
rm(list = ls(all.names = TRUE))

# Establish operating environment ----
home <- "C:/workspace/OrcaManuscript/"

# Install required packages ----
library(tidyverse)
library(viridis)
library(grid)                                                                   # required for plot annotation 
library(gridExtra)                                                              # for arranging multi-panel plots
library(xlsx)
library(ggpubr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(DescTools)                                                              # for computing CCC
library(ggpmisc)                                                                # for adding model parameters to plots
library(geojsonsf)                                                              # converts GeoJSON to sf object.
library(sf)                                                                     # package for handeling simple features.
library(raster)                                                                 # 
library(exactextractr)                                                          # efficient and exact extraction of raster statistics. 
library(rgeos)                                                                  # for validating which hole belongs to which exterior ring
library(lme4)                                                                   # For linear mixed effects models.
# library(lmerTest)                                                             # For extracting p values from linear mixed effects models (but can sometimes conflict with other packages, so use carfeully and be aware of errors!).
library(ggeffects)                                                              # For plotting mixed effects models.
library(sjPlot)                                                                 # For plotting mixed effects models.
library(raster)                                                                 # For extracting raster values
library(rgdal)
library(sp)


# #### Extract NDVI values from rasters ----
#   # Load data for raster extraction
#   feature_filename <- paste0(home, "data/20160725_AC_ORC - formated for exact extractr.geojson")
#   plots <- st_read(feature_filename, crs = 32607)                               # Import geoJSON as sf object, using st_read to allow the non-standerd CRS to be specified.
#   
#   
#   # Import rasters (Available from NERC Polar Data Centre - see readme) 
#   raster_018 <- raster(paste0(home, "inputs_NDVI/NDVI_019m_from_20160726.tif")) # import plot shapefiles.
#   raster_047 <- raster(paste0(home, "inputs_NDVI/NDVI_050m_from_20160730.tif")) # import plot shapefiles.
#   raster_119 <- raster(paste0(home, "inputs_NDVI/NDVI_120m_from_20160730.tif")) # import plot shapefiles.
#   raster_121 <- raster(paste0(home, "inputs_NDVI/NDVI_120m_from_20160726.tif")) # import plot shapefiles.
#   
#   
#   # Calculate mean NDVI for each polygon 
#   plots$mean_NDVI_018 <- exact_extract(raster_018, plots, 'mean')                        
#   plots$mean_NDVI_047 <- exact_extract(raster_047, plots, 'mean')                        
#   plots$mean_NDVI_119 <- exact_extract(raster_119, plots, 'mean')                        
#   plots$mean_NDVI_121 <- exact_extract(raster_121, plots, 'mean')       
#   
#   # Tidy dataframe
#   plots_df <- st_drop_geometry(plots)                                           # Create dataframe from simple features object (dropping geometry).
#   plots_df$EPSG <- NULL                                                         # Remove unecessary EPSG column from dataframe.
#   plots_df <- plots_df[order(plots_df$PlotID),]                                 # Order dataframe by PlotID.
#   
#   # Export NDVI values
#   write.csv(plots_df,"data/Extracted_NDVI.csv", row.names = FALSE)            # extracted NDVI values were added to the main_database file. ndvi_data <- read.csv("data/Extracted_NDVI.csv", header = T)                  # Read in NDVI values from Exact Extract pipeline.

  

# Load datasets ----
  dataset <- read.csv("data/main_database.csv", header = T)                     # Read in summary  data
  PF_observations <- read.csv("data/point_framing_observations.csv")            # Read in canopy height from point framing

# Data preparation ----
  # Point framing observations
  # remove unwanted index position.
    PF_observations <- PF_observations[,-1]                                
  
  # Pointframe data - Extracting only Height and PlotN from the pointframe data, and omit NAs
    PF_observations2 <- dplyr::select(PF_observations, Height, PlotN) %>%
        na.omit(PF_observations)

    # Convert heights from cm into m. Could be integrated into the above pipe.
    PF_observations2$Height <- PF_observations2$Height/100
    
    # Compute canopy height summary metrics (min, max, median, mean and quantiles) 
    # for each plot.
    PF_HAG_summary <- PF_observations2 %>%
        dplyr::select(Height, PlotN) %>%
        group_by(PlotN) %>%
        mutate(mean = mean(Height, na.rm = TRUE),
               min = min(Height, na.rm = TRUE),
               lower_IQR = quantile(Height, 0.25, na.rm = TRUE),
               median = median(Height, na.rm = TRUE),
               upper_ICR = quantile(Height, 0.75, na.rm = TRUE),
               max = max(Height, na.rm = TRUE)
               )

    # Return only unique row summarising canopy heights for each plot. This could be integrated into the preceeding pipe...
    PF_HAG_summary <- PF_HAG_summary[!duplicated(PF_HAG_summary$PlotN), ]
    
    # Remove unwanted/erronious individual height. This could be integrated into the preceding pipe...
    PF_HAG_summary <- subset(PF_HAG_summary, select = -Height)
    
    # Add point framing HAGs into main dataframe
    dataset$PF_HAG_min <- PF_HAG_summary$min
    dataset$PF_HAG_max <- PF_HAG_summary$max
    dataset$PF_HAG_median <- PF_HAG_summary$median
    dataset$PF_HAG_mean <- PF_HAG_summary$mean
    dataset$PF_HAG_upper_ICR <- PF_HAG_summary$upper_ICR
    dataset$PF_HAG_lower_IQR <- PF_HAG_summary$lower_IQR
    
  # Biomass observations
    # Calculate spatially normalised values for photosynthetic tissue
    dataset$leaf_biomass <- dataset$Comp2_Mass / dataset$plot_area_from_field_length  # Leaves.
    dataset$herbacious_biomass <- dataset$Comp3_Mass / dataset$plot_area_from_field_length  # Herbacious
    dataset$phytomass <- dataset$leaf_biomass + dataset$herbacious_biomass 
    
        
  # Create long form dataset for plotting canopy heights in SI
    PlotID <- rep(dataset$PlotID, times = 2)
    method <- rep(c('SfM', 'Point Framing'), each = 36)
    min <- c(dataset$HAG_plotmin_of_cellmax_m, dataset$PF_HAG_min)
    lower <- c(dataset$HAG_plot25percentile_of_cellmax_m, dataset$PF_HAG_lower_IQR)
    median <- c(dataset$HAG_plotmedian_of_cellmax_m, dataset$PF_HAG_median)
    upper <- c(dataset$HAG_plot75percentile_of_cellmax_m, dataset$PF_HAG_upper_ICR)
    max <- c(dataset$HAG_plotmax_of_cellmax_m, dataset$PF_HAG_max)
    dataset_long <- data.frame(PlotID, method, min, lower, median, upper, max)  # Create vectors into new dataframe
    rm(PlotID, method, min, lower, median, upper, max)

    
  ### Analysis of moss cover effect on NDVI-biomass relationships ###
    # Count number of locations sampled within each plot
    PF_count <- PF_observations %>%
      dplyr::select(-Species, -Status, -Tissue, -Count, -Height) %>% 
      group_by(PlotN) %>%
      distinct() %>% 
      count() %>% 
      rename(PF_obs = n)
    
    # Count the number of locations where moss was observed within each plot
    PF_moss <- PF_observations %>% 
      group_by(PlotN) %>% 
      filter(Species == "XXXothermoss") %>% 
      count() %>% 
      rename(moss_obs = n)
    
    # Add zero for plots containing no moss
    PF_moss[nrow(PF_moss) + 1,] = list("SP04", 0)
    PF_moss[nrow(PF_moss) + 1,] = list("SP07", 0)
    PF_moss[nrow(PF_moss) + 1,] = list("SP11", 0)
    PF_moss[nrow(PF_moss) + 1,] = list("SP34", 0)
    PF_moss <- dplyr::arrange(PF_moss, PlotN)  # re-order plots
    
    # Determine the proportion of PF locations that were moss for each plot
    PF_moss_summary <- PF_count
    PF_moss_summary$moss_obs <- PF_moss$moss_obs
    PF_moss_summary$moss_prop <- PF_moss_summary$moss_obs / PF_moss_summary$PF_obs
    
    # Add moss proportion to main dataframe
    dataset$moss_prop <- PF_moss_summary$moss_prop
    
  ### Extract distributions of canopy heights from SfM rasters ###
    
    # List of elevation rasters 
    raster_list <- list.files(path = "data/Height_rasters/",
                              full.names=TRUE)
    
    # Create dataframe
    SFM_heights <- data.frame()
    
    # Extract height values
    for (file in raster_list){
      PlotN <- substr(file, 37, 39)  # Extract PlotN from file
      
      # Read raster
      height_rast <- raster(file, band=2)  # NB. Band 2 specifies the maximum height recorded in each cell.
      
      # # Plot the distribution of values in the raster
      # hist(height_rast, main=paste0("Distribution of heights in plot ", PlotN, sep=""), 
      #      col= "purple", 
      #      maxpixels=10000,
      #      xlab = "height (m)", 
      #      xlim = c(0, 1.2))
      
      heights <- getValues(height_rast)  # Summarise height values
      
      heights <- na.exclude(heights)  # Remove NAs.
      
      height_distributions <- data.frame(heights)  # Convert to dataframe
      
      height_distributions$PlotN <- PlotN  # Add PlotN
      
      SFM_heights <- rbind(SFM_heights, height_distributions)  # Integrate heights to dataframe. 
      
      rm(heights, height_distributions)  # Tidy up
    }
    
    
    
    
### Data Analysis ----

  ### Analysis of canopy heights ####
    # Concordence correlation coefficient
    ccc.test <- CCC(dataset$HAG_plotmean_of_cellmax_m,
                dataset$PF_HAG_mean,
                ci = "z-transform",
                conf.level = 0.95,
                na.rm = FALSE)
    
    ccc.test$rho.c  # The concordance correlation coefficient
    
    # Fit model
    model_heights_pow <- nls(HAG_plotmean_of_cellmax_m ~ a*PF_HAG_mean^b, data = dataset,
                             start = list(a =1, b =1), na.action=na.exclude)
    summary(model_heights_pow)
    
    
    # Summarise bias
    HAG_rediduals <- dataset$HAG_plotmean_of_cellmax_m - dataset$PF_HAG_mean
    HAG_bias_mean <- round(mean(HAG_rediduals), 3)
    HAG_bias_mean
    
    HAG_bias_mean_SD <- round(sd(HAG_rediduals), 3)
    HAG_bias_mean_SD
    


    
    
    
  #### Analysis of canopy height as predictor of biomass ####
    # unconstrained intercept
    model_PFu <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$PF_HAG_mean)
    summary(model_PFu)
    
    # constrained intercept  
    model_PF <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$PF_HAG_mean + 0)
    summary(model_PF)
    
    # unconstrained intercept
    model_SfMu <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$HAG_plotmean_of_cellmax_m)
    summary(model_SfMu)
    
    # constrained intercept  
    modelSfM <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$HAG_plotmean_of_cellmax_m + 0)
    summary(modelSfM)

    
    
  ### Analysis of NDVI - biomass relationships ####
  # Exponential models
    exp_model_total_NDVI_018 <- nls(AGB_spatially_normalised_g_m2 ~ a*exp(b*mean_NDVI_018), data=dataset, start = list(a=30, b=4), na.action=na.exclude)
    exp_model_total_NDVI_047 <- nls(AGB_spatially_normalised_g_m2 ~ a*exp(b*mean_NDVI_047), data=dataset, start = list(a=30, b=4), na.action=na.exclude)
    exp_model_total_NDVI_119 <- nls(AGB_spatially_normalised_g_m2 ~ a*exp(b*mean_NDVI_119), data=dataset, start = list(a=30, b=4), na.action=na.exclude)
    exp_model_total_NDVI_121 <- nls(AGB_spatially_normalised_g_m2 ~ a*exp(b*mean_NDVI_121), data=dataset, start = list(a=30, b=4), na.action=na.exclude)
    exp_model_photo_NDVI_018 <- nls(phytomass ~ a*exp(b*mean_NDVI_018), data=dataset, start = list(a=30, b=1), na.action=na.exclude)
    exp_model_photo_NDVI_047 <- nls(phytomass ~ a*exp(b*mean_NDVI_047), data=dataset, start = list(a=30, b=1), na.action=na.exclude)
    exp_model_photo_NDVI_119 <- nls(phytomass ~ a*exp(b*mean_NDVI_119), data=dataset, start = list(a=30, b=1), na.action=na.exclude)
    exp_model_photo_NDVI_121 <- nls(phytomass ~ a*exp(b*mean_NDVI_121), data=dataset, start = list(a=30, b=1), na.action=na.exclude)
    exp_model_leaf_NDVI_018 <- nls(leaf_biomass ~ a*exp(b*mean_NDVI_018), data=dataset, start = list(a=30, b=1), na.action=na.exclude)
    exp_model_leaf_NDVI_047 <- nls(leaf_biomass ~ a*exp(b*mean_NDVI_047), data=dataset, start = list(a=30, b=1), na.action=na.exclude)
    exp_model_leaf_NDVI_119 <- nls(leaf_biomass ~ a*exp(b*mean_NDVI_119), data=dataset, start = list(a=30, b=1), na.action=na.exclude)
    exp_model_leaf_NDVI_121 <- nls(leaf_biomass ~ a*exp(b*mean_NDVI_121), data=dataset, start = list(a=30, b=1), na.action=na.exclude)
    
    summary(exp_model_total_NDVI_018)
    summary(exp_model_total_NDVI_047)
    summary(exp_model_total_NDVI_119)
    summary(exp_model_total_NDVI_121)
    summary(exp_model_photo_NDVI_018)
    summary(exp_model_photo_NDVI_047)
    summary(exp_model_photo_NDVI_119)
    summary(exp_model_photo_NDVI_121)
    summary(exp_model_leaf_NDVI_018)
    summary(exp_model_leaf_NDVI_047)
    summary(exp_model_leaf_NDVI_119)
    summary(exp_model_leaf_NDVI_121)

    
    

    
### Analysis of moss cover effect on NDVI-biomass relationships ####
  # Quantitative evaluation of the Figure S2 plots (below) suggests that plots with low
  # proportions of moss cover (observed in the point framing) may have lower NDVI values
  # for a given value of biomass.
    
  # Mixed effects model to test the effect of moss_prop on the NDVI-biomass relationships #
    # biomass as a function of NDVI, with moss_prop as a fixed effect.
    #### !!!!! Mixed effects model to be refined / confirmed !!!!! ####
    model <- lmer(phytomass ~ mean_NDVI_121 + (1|moss_prop), data = dataset)
    summary(model)
    
    
  # Q: Should we ignore moss under litter?
    
  # Q: Should we consider moss differently dpeending on whether it is the *first hit*?
    
    
    
    
    
    
    
# Visualisation ----
# Plotting themes ----
theme_coding <- function(){
    theme_bw()+
        theme(axis.text = element_text(size = 8),
              axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
              axis.title = element_text(size = 10),
              panel.grid = element_blank(),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
              plot.title = element_text(size = 12, vjust = 1, hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 6, face = "italic"),
              legend.key.size = unit(0.9,"line"),
              legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"),
              legend.position = c(0.9, 0.9))
}

# Set scaling parameters ----
max_agb <- 1.1*max(dataset$AGB_spatially_normalised_g_m2, na.rm = TRUE)
max_hag <- 1.1*max(max(dataset$HAG_plotmax_of_cellmax_m, na.rm = TRUE), max(dataset$PF_HAG_max))
max_mean_hag <- 1.1*max(max(dataset$HAG_plotmean_of_cellmax_m, na.rm = TRUE), max(dataset$PF_HAG_mean))
max_ndvi <- 0.91
min_ndvi <- 0.6
photo_biomass_max <- 1.12*max(dataset$phytomass)
spacing <- 2


# Figure 1. Site map ----
# Produced in ArcGIS

# Figure 2. Canopy height comparison ----
# mean point framing canopy height versus mean structure-from-motion canopy height

  # Create plot
    (Canopy_heights_plot <- ggplot(data = dataset,
                                  aes(x = PF_HAG_mean, y = HAG_plotmean_of_cellmax_m)) +
                geom_point(shape = 1) +
                labs(x = "Point Frame - Canopy Height (m)",
                     y = "SfM - Canopy height (m)") +
                theme_coding() +
                coord_cartesian(ylim = c(0, 1), xlim = c(0, 1), expand=FALSE) +
                geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
             stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
                          formula = y ~ x, na.rm = TRUE, coef.digits = 4, rr.digits = 3, size = 3, parse = TRUE,
                          label.x.npc = 0.90, label.y.npc = 0.10) +
                geom_smooth(method="lm", formula= y ~ x, se=FALSE, size=0.5, na.rm = TRUE))
  
  
    (Canopy_heights_plot <- ggplot(data = dataset,
                                   aes(x = PF_HAG_mean, y = HAG_plotmean_of_cellmax_m)) +
        geom_point(shape = 1) +
        labs(x = "Point Frame - Canopy Height (m)",
             y = "SfM - Canopy height (m)") +
        theme_coding() +
        coord_cartesian(ylim = c(0, 1), xlim = c(0, 1), expand=FALSE) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
        stat_function(fun = function(x) (coef(summary(model_heights_pow))[, "Estimate"])[1]*(x)^(coef(summary(model_heights_pow))[, "Estimate"])[2],
                      aes(), size = 1, lty = "solid", colour="Black") +
        annotate("text", x = 0.8, y = 0.1, label = "italic(Y)==1.085~italic(X)^0.715", parse = TRUE,
                 color = "black", size = 4, family = "serif")
      )
    
  # Export plot
    png(filename = "plots/Figure 2 - Canopy Heights2.png", width = 10, height = 10, units = "cm", res = 400)
    plot(Canopy_heights_plot)
    dev.off()



# Figure 3. Biomass prediction ----
# illustrating relationships between (SfM and PF) canopy height and aboveground biomass.

# Create plot
  (biomass_CH_SfM <- ggplot(data = dataset,
                              aes(x = HAG_plotmean_of_cellmax_m,
                                  y = AGB_spatially_normalised_g_m2)) + 
     geom_point(shape = 1, na.rm = TRUE) +
     theme_coding() +
     coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 1), expand=FALSE) +
     labs(x = expression("Canopy height (m)"),
          y = expression("Dry biomass (g m"^"-2"*")"),
          title = "SfM") +
     stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
                  formula = y ~ x-1, na.rm = TRUE, coef.digits = 4, rr.digits = 3, size = 3, parse = TRUE,
                  label.x.npc = 0.05, label.y.npc = 0.95) +
     theme(legend.position = c(0.15, 0.9)) +
     geom_smooth(method="lm", formula= y ~ x-1, se=TRUE, size=0.5, na.rm = TRUE))
  
  
  (biomass_CH_PF <- ggplot(data = dataset,
                             aes(x = PF_HAG_mean,
                                 y = AGB_spatially_normalised_g_m2)) + 
      geom_point(shape = 1, na.rm = TRUE) +
      theme_coding() +
      coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 1), expand=FALSE) +
      labs(x = expression("Canopy height (m)"),
           y = expression("Dry biomass (g m"^"-2"*")"),
           title = "Point Framing") +
      stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
                   formula = y ~ x-1, na.rm = TRUE, coef.digits = 4, rr.digits = 3, size = 3, parse = TRUE,
                   label.x.npc = 0.05, label.y.npc = 0.95) +
      theme(legend.position = c(0.15, 0.9)) +
      geom_smooth(method="lm", formula= y ~ x-1, se=TRUE, size=0.5, na.rm = TRUE))
  
# Combine plots
  biomass_plots <- ggpubr::ggarrange(biomass_CH_PF, biomass_CH_SfM,
                                  heights = c(9,9),
                                  labels = c("(a)", "(b)"),
                                  ncol = 2, nrow = 1,
                                  align = "h")
  
  
# Export figure
  png(filename="plots/Figure 3 - Height vesus biomass2.png", width=14, height=7, units="cm", res=400)
  plot(biomass_plots)
  dev.off()
  


  # Figure 4. NDVI vs. biomass ----
  # Illustrating relationships between NDVI and various biomass components
  
  # Create plots
  # Total biomass
  NDVI_vs_total_biomass_121 <- ggplot(data = dataset, aes(x = mean_NDVI_121, y = AGB_spatially_normalised_g_m2)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(
      # x = expression("mean NDVI"),
      x = expression(""),
      y = expression(atop("Total", paste ("biomass (g m"^"-2"*")"))),
      title = "0.121 m grain") +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_total_NDVI_121))[, "Estimate"])[1]*exp((coef(summary(exp_model_total_NDVI_121))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_total_biomass_119 <- ggplot(data = dataset, aes(x = mean_NDVI_119, y = AGB_spatially_normalised_g_m2)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(
      x = expression(""),
      y = expression(""),
      title = "0.119 m grain") +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_total_NDVI_119))[, "Estimate"])[1]*exp((coef(summary(exp_model_total_NDVI_119))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_total_biomass_047 <- ggplot(data = dataset, aes(x = mean_NDVI_047, y = AGB_spatially_normalised_g_m2)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(
      x = expression(""),
      y = expression(""),
      title = "0.047 m grain") +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_total_NDVI_047))[, "Estimate"])[1]*exp((coef(summary(exp_model_total_NDVI_047))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_total_biomass_018 <- ggplot(data = dataset, aes(x = mean_NDVI_018,  y = AGB_spatially_normalised_g_m2)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(
      x = expression(""),
      y = expression(""),
      title = "0.018 m grain") +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_total_NDVI_018))[, "Estimate"])[1]*exp((coef(summary(exp_model_total_NDVI_018))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  # Photosynthetic biomass
  NDVI_vs_photo_biomass_121 <- ggplot(data = dataset, aes(x = mean_NDVI_121,  y = phytomass)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(
      x = expression(""),
      y = expression(atop("Photosynthetic", paste ("biomass (g m"^"-2"*")")))) +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_photo_NDVI_121))[, "Estimate"])[1]*exp((coef(summary(exp_model_photo_NDVI_121))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_photo_biomass_119 <- ggplot(data = dataset, aes(x = mean_NDVI_119,  y = phytomass)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(
      x = expression(""),
      y = expression("")
    ) +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_photo_NDVI_119))[, "Estimate"])[1]*exp((coef(summary(exp_model_photo_NDVI_119))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_photo_biomass_047 <- ggplot(data = dataset, aes(x = mean_NDVI_047, y = phytomass)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(
      x = expression(""),
      y = expression("")
    ) +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_photo_NDVI_047))[, "Estimate"])[1]*exp((coef(summary(exp_model_photo_NDVI_047))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_photo_biomass_018 <- ggplot(data = dataset, aes(x = mean_NDVI_018,  y = phytomass)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(
      x = expression(""),
      y = expression("")
    ) +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_photo_NDVI_018))[, "Estimate"])[1]*exp((coef(summary(exp_model_photo_NDVI_018))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  # Leaf biomass
  # Create plots
  NDVI_vs_leaf_biomass_121 <- ggplot(data = dataset, aes(x = mean_NDVI_121, y = leaf_biomass)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(x = expression("NDVI"),
         y = expression(atop("Leaf", paste ("biomass (g m"^"-2"*")")))
    ) +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_leaf_NDVI_121))[, "Estimate"])[1]*exp((coef(summary(exp_model_leaf_NDVI_121))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_leaf_biomass_119 <- ggplot(data = dataset, aes(x = mean_NDVI_119, y = leaf_biomass)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(x = expression("NDVI"),
         y = expression("")
    ) +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_leaf_NDVI_119))[, "Estimate"])[1]*exp((coef(summary(exp_model_leaf_NDVI_119))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_leaf_biomass_047 <- ggplot(data = dataset, aes(x = mean_NDVI_047, y = leaf_biomass)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(x = expression("NDVI"),
         y = expression("")
    ) +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_leaf_NDVI_047))[, "Estimate"])[1]*exp((coef(summary(exp_model_leaf_NDVI_047))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  NDVI_vs_leaf_biomass_018 <- ggplot(data = dataset, aes(x = mean_NDVI_018, y = leaf_biomass)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    labs(x = expression("NDVI"),
         y = expression("")
    ) +
    # geom_smooth(method='lm', formula= y~x, se=FALSE) +
    stat_function(fun = function(x) (coef(summary(exp_model_leaf_NDVI_018))[, "Estimate"])[1]*exp((coef(summary(exp_model_leaf_NDVI_018))[, "Estimate"])[2]*x),
                  aes(), size = 1, lty = "solid") +
    coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
    theme_coding() +
    theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
  
  
  # Combine plots
  NDVI_biomass <- ggpubr::ggarrange(NDVI_vs_total_biomass_121,
                                    NDVI_vs_total_biomass_119,
                                    NDVI_vs_total_biomass_047, 
                                    NDVI_vs_total_biomass_018,
                                    NDVI_vs_photo_biomass_121,
                                    NDVI_vs_photo_biomass_119,
                                    NDVI_vs_photo_biomass_047,
                                    NDVI_vs_photo_biomass_018,
                                    NDVI_vs_leaf_biomass_121,
                                    NDVI_vs_leaf_biomass_119,
                                    NDVI_vs_leaf_biomass_047,
                                    NDVI_vs_leaf_biomass_018,
                                    heights = c(9,9,9),
                                    ncol = 4, nrow = 3,
                                    align = "v",
                                    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"),
                                    font.label = list(size = 11, face = "bold"))
  
  # Export figure
  png(filename="plots/Figure 4 - NDVI vs biomass.png", width=22, height=13, units="cm", res=400)
  plot(NDVI_biomass)
  dev.off()
  

  
# Figure S1. Boxplot of canopy height observations ---- 
# Create boxplot
  (HAG_boxplot <- ggplot(data = dataset_long,
                         aes(x = reorder(PlotID, median, FUN = "median"),
                             y = median, fill = method)) +
          geom_boxplot(aes(fill = method,
                           ymin = min,
                           lower = lower,
                           middle = median,
                           upper = upper,
                           ymax = max),
                       stat = "identity", colour = "grey60",
                       outlier.colour = "white", outlier.size = 0,
                       position = position_dodge(0.7), width = 0.7) +
          geom_boxplot(aes(fill = method,
                           upper = upper,
                           lower = lower,
                           middle = median,
                           ymin = ..lower..,
                           ymax = ..upper..),
                       stat = "identity",
                       position = position_dodge(0.7), width = 0.7) +  # Added again to include black boarder.
          scale_y_continuous(lim = c(0, 1)) +
          scale_color_manual(values = c("white", "darkgrey")) +
          scale_fill_manual(values = c("white", "darkgrey")) +
          labs(x = "Plot ID", y = "Canopy height (m)") +
          theme_coding() +
          theme(legend.position = c(0.15, 0.9), 
                axis.line.x = element_line(color="black", size = 0.5),
                axis.line.y = element_line(color="black", size = 0.5),
                axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5)
          ))
        
# Export boxplot
  png(filename = "plots/Figure S1 - Height boxplot.png", width = 16, height = 11, units = "cm", res = 400)
  plot(HAG_boxplot)
  dev.off()


### Figure S2. Effect of moss on NDVI-Biomass relationship ####
  # What is the effect of moss_prop on the NDVI-biomass relationships?
  
  # Total biomass
    plot1 <- ggplot(data = dataset, aes(mean_NDVI_121,
                                       AGB_spatially_normalised_g_m2,
                                       color=moss_prop,
                                       size=moss_prop)) +
      geom_point(na.rm = TRUE) + 
      coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
      labs(
        x = expression("mean NDVI"),
        y = expression("Total biomass (g m"^"-2"*")")
        # title = "Total biomass (0.121 m grain)"
        ) +
      theme_coding() +
      theme(legend.position = c(0.9, 0.5))
  

  # Phytomas
    plot2 <- ggplot(data = dataset, aes(mean_NDVI_121,
                                       phytomass,
                                       color=moss_prop,
                                       size=moss_prop)) +
      geom_point(na.rm = TRUE) + 
      coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
      labs(
        x = expression("mean NDVI"),
        y = expression("Phytomass (g m"^"-2"*")")
        # title = "Phytomass (0.121 m grain)"
        ) +
      theme_coding() +
      theme(legend.position = c(0.9, 0.5))

  
  # Leaf biomass
    plot3 <- ggplot(data = dataset, aes(mean_NDVI_121,
                                       leaf_biomass,
                                       color=moss_prop,
                                       size=moss_prop)) +
      geom_point(na.rm = TRUE) + 
      coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
      labs(
        x = expression("mean NDVI"),
        y = expression("Leaf biomass (g m"^"-2"*")")
        # title = "Leaf biomass (0.121 m grain)"
        ) +
      theme_coding() +
      theme(legend.position = c(0.9, 0.5))

  
  # Combine plots
  moss_plots <- ggpubr::ggarrange(plot1, plot2, plot3,
                                     heights = c(6),
                                     labels = c("(a)", "(b)", "(c)"),
                                     ncol = 1, nrow = 3,
                                     align = "v")
  
  
  # Export figure
  png(filename="plots/Figure S2 - Moss interaction with NDVI vesus biomass.png", width=11, height=24.5, units="cm", res=300)
  plot(moss_plots)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  


  
  
      






