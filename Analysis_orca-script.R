# Clear working environment ----
rm(list = ls(all.names = TRUE))

# Establish operating environment ----
home <- "C:/workspace/OrcaManuscript/"

# Install required packages ----
library(suncalc)
library(tidyverse)
library(viridis)
library(grid)                                                                   # required for plot annotation 
library(gridExtra)                                                              # for arranging multi-panel plots
library(xlsx)
library(ggpubr)
library(ggplot2)
library(DescTools)                                                              # for computing CCC
library(ggpmisc)                                                                # for adding model parameters to plots
library(geojsonsf)                                                              # converts GeoJSON to sf object.
library(sf)                                                                     # package for handeling simple features.
library(raster)                                                                 # 
library(exactextractr)                                                          # efficient and exact extraction of raster statistics. 
library(rgeos)                                                                  # for validating which hole belongs to which exterior ring
library(lme4)                                                                   # For linear mixed effects models.
library(ggeffects)                                                              # For plotting mixed effects models.
library(sjPlot)                                                                 # For plotting mixed effects models.
library(raster)                                                                 # For extracting raster values
library(rgdal)
library(sp)
library(miscTools)


#### Extract NDVI values from rasters ----
  #   # Load data for raster extraction
  #   feature_filename <- paste0(home, "data/20160725_AC_ORC - formated for exact extractr.geojson")
  #   plots <- st_read(feature_filename, crs = 32607)                               # Import geoJSON as sf object, using st_read to allow the non-standerd CRS to be specified.
  # 
  # 
  #   # Import rasters (Available from the NERC Polar Data Centre - see readme for details)
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



#### Load Data ----
    loc_time_data <- read.csv("data/location_time_for_suncalc.csv", header = T)         # Read in locations and survey times for the suncalc
    dataset <- read.csv("data/main_database.csv", header = T)                     # Read in summary  data
    PF_observations <- read.csv("data/point_framing_observations.csv")            # Read in canopy height from point framing


# Extract solar angles for the surveys
    datetime_UTC <- as.character(loc_time_data$datetime_UTC)
    SunPosition <- getSunlightPosition(date = datetime_UTC, lat = 69.57133, lon = -138.8909)
    
    # Function to convert from radians to degrees
    rad2deg <- function(rad) {(rad * 180) / (pi)}
    
    SunPosition$altitude_degrees <- rad2deg(SunPosition$altitude)
    SunPosition



#### Data Preparation ----
    # Point framing observations
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
    
    # Return only unique row summarising canopy heights for each plot.
    PF_HAG_summary <- PF_HAG_summary[!duplicated(PF_HAG_summary$PlotN), ]
    
    # Remove unwanted/erronious individual height.
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
    
    # Create long form dataset for NDVI boxplot
    # NB. It is not straightforward to extract summary statistics like IQR for 
    # the NDVI rasters, because these metrics are not readily accessible via the
    # exactextractr package. Mean, mode, max and min are available.
    PlotID <- rep(dataset$PlotID, times = 4)
    grain <- rep(c('0.018', '0.047', '0.119', '0.121'), each = 36)
    mean_NDVI <- c(dataset$mean_NDVI_018, 
                   dataset$mean_NDVI_047, 
                   dataset$mean_NDVI_119, 
                   dataset$mean_NDVI_121)
    NDVI_data_long <- data.frame(PlotID, grain, mean_NDVI)  # Create vectors into new dataframe
    rm(PlotID, grain, mean_NDVI)  # Tidy up
    
  ### Summarise moss cover within each plot ###
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
    raster_list <- list.files(path = "data/Height_rasters/", full.names=TRUE)   # List of elevation rasters
    
    SFM_heights <- data.frame()                                                 # Create dataframe
    
    # Extract height values
    for (file in raster_list){
      PlotN <- substr(file, 37, 39)                                             # Extract PlotN from file
      
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


#### Plotting themes and scaling parameters ----
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
    
  ## Set scaling parameters
    max_agb <- 1.1*max(dataset$AGB_spatially_normalised_g_m2, na.rm = TRUE)
    max_hag <- 1.1*max(max(dataset$HAG_plotmax_of_cellmax_m, na.rm = TRUE), max(dataset$PF_HAG_max))
    max_mean_hag <- 1.1*max(max(dataset$HAG_plotmean_of_cellmax_m, na.rm = TRUE), max(dataset$PF_HAG_mean))
    max_ndvi <- 0.91
    min_ndvi <- 0.6
    phyto_biomass_max <- 1.12*max(dataset$phytomass)
    spacing <- 2

    
#### Data Analysis ----
#### Comparison of canopy height measurements ----    
  # Analysis
    # Concordence correlation coefficient
      ccc.test <- CCC(dataset$HAG_plotmean_of_cellmax_m,
                      dataset$PF_HAG_mean,
                      ci = "z-transform",
                      conf.level = 0.95,
                      na.rm = FALSE)
      
      ccc.test$rho.c                                                            # The concordance correlation coefficient
      
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
      
  # Visualistation (Figure 2)
    # mean point framing canopy height versus mean structure-from-motion canopy height
    # Create plot
      {
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
      }

      
#### Predictors of biomass (height and NDVI) ---- 
  # Testing canopy height (from point framing and photogrammetry) and NDVI as predictors of biomass.
  # Analysis
    # Point framing    
      model_PFu <- lm(AGB_spatially_normalised_g_m2 ~ PF_HAG_mean, data=dataset)  # Unconstrained intercept
      model_PF <- lm(AGB_spatially_normalised_g_m2 ~ PF_HAG_mean + 0, data=dataset)  # constrained intercept
      model_SfMu <- lm(AGB_spatially_normalised_g_m2 ~ HAG_plotmean_of_cellmax_m, data=dataset)  # unconstrained intercept
      model_SfM <- lm(AGB_spatially_normalised_g_m2 ~ HAG_plotmean_of_cellmax_m + 0, data=dataset)  # constrained intercept
      
      # Tabulate and export model parameters
      {
      model_list <- list(model_PFu)
      model_results1 <- bind_rows(lapply(model_list, function(model){
        model_glance <- broom::glance(model)
        model_tidy <- broom::tidy(model)
        return(data.frame(Canopy_height_method = "Point Framing",
                          model_form = "Y = a X + b",
                          a = paste0(round(model_tidy$estimate[2], 1), " ± ", 
                                     round(model_tidy$std.error[2], 1)),
                          b = paste0(round(model_tidy$estimate[1], 1), " ± ", 
                                     round(model_tidy$std.error[1], 1)),
                          r2 = round(model_glance$r.squared, 2),
                          resid = round(model_glance$sigma, 3), stringsAsFactors = F))
      }))
      
      
      model_list <- list(model_PF)
      model_results2 <- bind_rows(lapply(model_list, function(model){
        model_glance <- broom::glance(model)
        model_tidy <- broom::tidy(model)
        return(data.frame(Canopy_height_method = "Point Framing",
                          model_form = "Y = a X",
                          a = paste0(round(model_tidy$estimate[1], 1), " ± ", 
                                     round(model_tidy$std.error[1], 1)),
                          b = paste0("0"),
                          r2 = round(model_glance$r.squared, 2),
                          resid = round(model_glance$sigma, 3), stringsAsFactors = F))
      }))
      
      
      model_list <- list(model_SfMu)
      model_results3 <- bind_rows(lapply(model_list, function(model){
        model_glance <- broom::glance(model)
        model_tidy <- broom::tidy(model)
        return(data.frame(Canopy_height_method = "Photogrammetry",
                          model_form = "Y = a X + b",
                          a = paste0(round(model_tidy$estimate[2], 1), " ± ", 
                                     round(model_tidy$std.error[2], 1)),
                          b = paste0(round(model_tidy$estimate[1], 1), " ± ", 
                                     round(model_tidy$std.error[1], 1)),
                          r2 = round(model_glance$r.squared, 2),
                          resid = round(model_glance$sigma, 3), stringsAsFactors = F))
      }))
      
      
      model_list <- list(model_SfM)
      model_results4 <- bind_rows(lapply(model_list, function(model){
        model_glance <- broom::glance(model)
        model_tidy <- broom::tidy(model)
        return(data.frame(Canopy_height_method = "Photogrammetry",
                          model_form = "Y = a X",
                          a = paste0(round(model_tidy$estimate[1], 1), " ± ", 
                                     round(model_tidy$std.error[1], 1)),
                          b = paste0("0"),
                          r2 = round(model_glance$r.squared, 2),
                          resid = round(model_glance$sigma, 3), stringsAsFactors = F))
      }))
      
      # Combine tables
      blended <- rbind(model_results1,
                       model_results2,
                       model_results3,
                       model_results4)
      
      # Export model parameters to table
      write.csv(blended, file = "tables/Table S1 model fits.csv", row.names = F)
      }
      
  # Visualisation
      # Create plots
      {
      (biomass_CH_SfM <- ggplot(data = dataset,
                                aes(x = HAG_plotmean_of_cellmax_m,
                                    y = AGB_spatially_normalised_g_m2)) + 
         geom_point(shape = 1, na.rm = TRUE) +
         theme_coding() +
         coord_cartesian(ylim = c(0, 3000), xlim = c(0, 1), expand=FALSE) +
         labs(x = expression("Canopy height (m)"),
              y = expression("Dry biomass (g m"^"-2"*")"),
              title = "SfM") +
         stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
                      formula = y ~ x-1, na.rm = TRUE, coef.digits = 4, rr.digits = 2, size = 2.5, parse = TRUE,
                      label.x.npc = 0.03, label.y.npc = 0.99) +
         theme(legend.position = c(0.15, 0.9)) +
         geom_smooth(method="lm", formula= y ~ x-1, se=TRUE, size=0.5, na.rm = TRUE))
      
      
      (biomass_CH_PF <- ggplot(data = dataset,
                               aes(x = PF_HAG_mean,
                                   y = AGB_spatially_normalised_g_m2)) + 
          geom_point(shape = 1, na.rm = TRUE) +
          theme_coding() +
          coord_cartesian(ylim = c(0, 3000), xlim = c(0, 1), expand=FALSE) +
          labs(x = expression("Canopy height (m)"),
               y = expression("Dry biomass (g m"^"-2"*")"),
               title = "Point Framing") +
          stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
                       formula = y ~ x-1, na.rm = TRUE, coef.digits = 4, rr.digits = 2, size = 2.5, parse = TRUE,
                       label.x.npc = 0.03, label.y.npc = 0.99) +
          theme(legend.position = c(0.15, 0.9)) +
          geom_smooth(method="lm", formula= y ~ x-1, se=TRUE, size=0.5, na.rm = TRUE))
      
      (biomass_NDVI <- ggplot(data = dataset,
                              aes(x = mean_NDVI_121,
                                  y = AGB_spatially_normalised_g_m2)) + 
          geom_point(shape = 1, na.rm = TRUE) +
          theme_coding() +
          coord_cartesian(ylim = c(4.5, 3000), xlim = c(0.65, 0.85), expand=FALSE) +
          labs(x = expression("NDVI (0.121 m grain)"),
               y = expression("Dry biomass (g m"^"-2"*")"),
               title = "NDVI") +
          stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
                       formula = y ~ x, na.rm = TRUE, coef.digits = 4, rr.digits = 2, size = 2.5, parse = TRUE,
                       label.x.npc = 0.03, label.y.npc = 0.99) +
          theme(legend.position = c(0.15, 0.9)) +
          geom_smooth(method="lm", formula= y ~ x, se=TRUE, size=0.5, na.rm = TRUE))
      
      # ln version
      # (biomass_NDVI <- ggplot(data = dataset,
      #                         aes(x = mean_NDVI_121,
      #                             y = log(AGB_spatially_normalised_g_m2))) + 
      #     geom_point(shape = 1, na.rm = TRUE) +
      #     theme_coding() +
      #     coord_cartesian(ylim = c(4.5, log(5000)), xlim = c(0.65, 0.85), expand=FALSE) +
      #     labs(x = expression("NDVI (0.121 m grain)"),
      #          y = expression("ln ( Dry biomass (g m"^"-2"*"))"),
      #          title = "NDVI") +
      #     stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
      #                  formula = y ~ x, na.rm = TRUE, coef.digits = 4, rr.digits = 2, size = 2.5, parse = TRUE,
      #                  label.x.npc = 0.03, label.y.npc = 0.99) +
      #     theme(legend.position = c(0.15, 0.9)) +
      #     geom_smooth(method="lm", formula= y ~ x, se=TRUE, size=0.5, na.rm = TRUE))
      
      # Combine plots
      biomass_plots <- ggpubr::ggarrange(biomass_CH_PF, biomass_CH_SfM, biomass_NDVI,
                                         heights = c(10,10,10),
                                         labels = c("(a)", "(b)", "(c)"),
                                         ncol = 3, nrow = 1,
                                         align = "h")
 
      # Export figure
      png(filename="plots/Figure 3 - Biomass Predictions.png", width=20, height=7, units="cm", res=400)
      plot(biomass_plots)
      dev.off()
      }
      
      

### Analysis of NDVI - biomass relationships ----
# Testing relationships between NDVI and various biomass components

  # Analysis
    # Logarithmic models
      log_model_total_NDVI_121 <- lm(log(AGB_spatially_normalised_g_m2) ~ mean_NDVI_121, data=dataset, na.action=na.exclude)
      log_model_total_NDVI_119 <- lm(log(AGB_spatially_normalised_g_m2) ~ mean_NDVI_119, data=dataset, na.action=na.exclude)
      log_model_total_NDVI_047 <- lm(log(AGB_spatially_normalised_g_m2) ~ mean_NDVI_047, data=dataset, na.action=na.exclude)
      log_model_total_NDVI_018 <- lm(log(AGB_spatially_normalised_g_m2) ~ mean_NDVI_018, data=dataset, na.action=na.exclude)
      log_model_phyto_NDVI_121 <- lm(log(phytomass) ~ mean_NDVI_121, data=dataset, na.action=na.exclude)
      log_model_phyto_NDVI_119 <- lm(log(phytomass) ~ mean_NDVI_119, data=dataset, na.action=na.exclude)
      log_model_phyto_NDVI_047 <- lm(log(phytomass) ~ mean_NDVI_047, data=dataset, na.action=na.exclude)
      log_model_phyto_NDVI_018 <- lm(log(phytomass) ~ mean_NDVI_018, data=dataset, na.action=na.exclude)
      log_model_leaf_NDVI_121 <- lm(log(leaf_biomass) ~ mean_NDVI_121, data=dataset, na.action=na.exclude)
      log_model_leaf_NDVI_119 <- lm(log(leaf_biomass) ~ mean_NDVI_119, data=dataset, na.action=na.exclude)
      log_model_leaf_NDVI_047 <- lm(log(leaf_biomass) ~ mean_NDVI_047, data=dataset, na.action=na.exclude)
      log_model_leaf_NDVI_018 <- lm(log(leaf_biomass) ~ mean_NDVI_018, data=dataset, na.action=na.exclude)

  
    # Compile model objects
      log_models <- list(log_model_total_NDVI_121,
                         log_model_total_NDVI_119,
                         log_model_total_NDVI_047,
                         log_model_total_NDVI_018,
                         log_model_phyto_NDVI_121,
                         log_model_phyto_NDVI_119,
                         log_model_phyto_NDVI_047,
                         log_model_phyto_NDVI_018,
                         log_model_leaf_NDVI_121,
                         log_model_leaf_NDVI_119,
                         log_model_leaf_NDVI_047,
                         log_model_leaf_NDVI_018)
  
    # Tabulate model parameters
      log_model_results <- bind_rows(lapply(log_models, function(model){
        model_glance <- broom::glance(model)
        model_tidy <- broom::tidy(model)
        return(data.frame(ndvi_grain = gsub("mean_NDVI_", "0\\.", model_tidy$term[2]),
                          model_form = "ln(Y) = a X + b",
                          a = paste0(round(model_tidy$estimate[2], 3), " ± ", 
                                     round(model_tidy$std.error[2], 3)),
                          b = paste0(round(model_tidy$estimate[1], 3), " ± ", 
                                     round(model_tidy$std.error[1], 3)),
                          r2 = round(model_glance$r.squared, 2),
                          resid = round(model_glance$sigma, 3), stringsAsFactors = F))
      }))
  
    # Export model parameters to table
      write.csv(log_model_results, file = "tables/Table 2 log model fits.csv", row.names = F)

      
  # Visualisation
      {
        # Create plots
        # Total biomass
        # Set parameters for the log plots.
        phyto_biomass_max_log <- ceiling(max(log(dataset$phytomass)))*1.1
        max_agb_log <- ceiling(max(log(dataset$AGB_spatially_normalised_g_m2)))*1.1
        phyto_biomass_min_log <- floor(min(log(dataset$phytomass)))*0.9
        min_agb_log <- floor(min(log(dataset$AGB_spatially_normalised_g_m2)))*0.9
        
        NDVI_vs_total_biomass_121 <- ggplot(data = dataset, aes(x = mean_NDVI_121, y = log(AGB_spatially_normalised_g_m2))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(
            # x = expression("mean NDVI"),
            x = expression(""),
            y = expression(atop("ln (Total", paste ("biomass (g m"^"-2"*"))"))),
            title = "0.121 m grain") +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_total_NDVI_121))[, "Estimate"])[1]*exp((coef(summary(exp_model_total_NDVI_121))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(min_agb_log, max_agb_log), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_total_biomass_119 <- ggplot(data = dataset, aes(x = mean_NDVI_119, y = log(AGB_spatially_normalised_g_m2))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(
            x = expression(""),
            y = expression(""),
            title = "0.119 m grain") +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_total_NDVI_119))[, "Estimate"])[1]*exp((coef(summary(exp_model_total_NDVI_119))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(min_agb_log, max_agb_log), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_total_biomass_047 <- ggplot(data = dataset, aes(x = mean_NDVI_047, y = log(AGB_spatially_normalised_g_m2))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(
            x = expression(""),
            y = expression(""),
            title = "0.047 m grain") +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_total_NDVI_047))[, "Estimate"])[1]*exp((coef(summary(exp_model_total_NDVI_047))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(min_agb_log, max_agb_log), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_total_biomass_018 <- ggplot(data = dataset, aes(x = mean_NDVI_018,  y = log(AGB_spatially_normalised_g_m2))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(
            x = expression(""),
            y = expression(""),
            title = "0.018 m grain") +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_total_NDVI_018))[, "Estimate"])[1]*exp((coef(summary(exp_model_total_NDVI_018))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(min_agb_log, max_agb_log), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        # Photosynthetic biomass
        NDVI_vs_phyto_biomass_121 <- ggplot(data = dataset, aes(x = mean_NDVI_121,  y = log(phytomass))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(
            x = expression(""),
            y = expression(atop("ln (Photosynthetic", paste ("biomass (g m"^"-2"*"))")))) +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_phyto_NDVI_121))[, "Estimate"])[1]*exp((coef(summary(exp_model_phyto_NDVI_121))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(phyto_biomass_min_log, phyto_biomass_max_log), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_phyto_biomass_119 <- ggplot(data = dataset, aes(x = mean_NDVI_119,  y = log(phytomass))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(
            x = expression(""),
            y = expression("")
          ) +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_phyto_NDVI_119))[, "Estimate"])[1]*exp((coef(summary(exp_model_phyto_NDVI_119))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(phyto_biomass_min_log, phyto_biomass_max_log), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_phyto_biomass_047 <- ggplot(data = dataset, aes(x = mean_NDVI_047, y = log(phytomass))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(
            x = expression(""),
            y = expression("")
          ) +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_phyto_NDVI_047))[, "Estimate"])[1]*exp((coef(summary(exp_model_phyto_NDVI_047))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(phyto_biomass_min_log, phyto_biomass_max_log), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_phyto_biomass_018 <- ggplot(data = dataset, aes(x = mean_NDVI_018,  y = log(phytomass))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(
            x = expression(""),
            y = expression("")
          ) +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_phyto_NDVI_018))[, "Estimate"])[1]*exp((coef(summary(exp_model_phyto_NDVI_018))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(phyto_biomass_min_log, phyto_biomass_max_log), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        # Leaf biomass
        # Create plots
        NDVI_vs_leaf_biomass_121 <- ggplot(data = dataset, aes(x = mean_NDVI_121, y = log(leaf_biomass))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(x = expression("NDVI"),
               y = expression(atop("ln (Leaf", paste ("biomass (g m"^"-2"*"))")))
          ) +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_leaf_NDVI_121))[, "Estimate"])[1]*exp((coef(summary(exp_model_leaf_NDVI_121))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(3, 5.5), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_leaf_biomass_119 <- ggplot(data = dataset, aes(x = mean_NDVI_119, y = log(leaf_biomass))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(x = expression("NDVI"),
               y = expression("")
          ) +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_leaf_NDVI_119))[, "Estimate"])[1]*exp((coef(summary(exp_model_leaf_NDVI_119))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(3, 5.5), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_leaf_biomass_047 <- ggplot(data = dataset, aes(x = mean_NDVI_047, y = log(leaf_biomass))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(x = expression("NDVI"),
               y = expression("")
          ) +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_leaf_NDVI_047))[, "Estimate"])[1]*exp((coef(summary(exp_model_leaf_NDVI_047))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(3, 5.5), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        NDVI_vs_leaf_biomass_018 <- ggplot(data = dataset, aes(x = mean_NDVI_018, y = log(leaf_biomass))) + 
          geom_point(shape = 1, na.rm = TRUE) +
          labs(x = expression("NDVI"),
               y = expression("")
          ) +
          geom_smooth(method='lm', formula= y~x, se=TRUE, size = 1, lty = "solid", col="black") +
          # geom_smooth(method='lm', formula= y~x, se=FALSE) +
          # stat_function(fun = function(x) (coef(summary(exp_model_leaf_NDVI_018))[, "Estimate"])[1]*exp((coef(summary(exp_model_leaf_NDVI_018))[, "Estimate"])[2]*x),
          #               aes(), size = 1, lty = "solid") +
          coord_cartesian(ylim = c(3, 5.5), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
          theme_coding() +
          theme(plot.margin = margin(t = spacing, r = spacing, b = spacing, l = spacing, unit = "pt"))
        
        
        # Combine plots
        NDVI_biomass <- ggpubr::ggarrange(NDVI_vs_total_biomass_121,
                                          NDVI_vs_total_biomass_119,
                                          NDVI_vs_total_biomass_047, 
                                          NDVI_vs_total_biomass_018,
                                          NDVI_vs_phyto_biomass_121,
                                          NDVI_vs_phyto_biomass_119,
                                          NDVI_vs_phyto_biomass_047,
                                          NDVI_vs_phyto_biomass_018,
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
        png(filename="plots/Figure 4 - NDVI vs biomass - log_2.png", width=22, height=13, units="cm", res=400)
        plot(NDVI_biomass)
        dev.off()
      }
      
      
    
      
      
    
### Analysis of moss cover effect on NDVI-biomass relationships ####
    # Looking for the inteaction between the proportion of moss cover ('moss_prop') with NDVI and Biomass.
    # looking for an mean_NDVI_121 interaction with moss_prop. The coefficient of that interaction effect 
    # indicates how moss_prop influences the phytomass/NDVI relationship

    # Check the distribution of phytomass. Suggests phytomass should be transformed to normalise the distribution.
    hist(dataset$phytomass)
    hist(log(dataset$phytomass))
    
    model3 <- lm(phytomass ~ mean_NDVI_121*moss_prop, data = dataset)
    summary(model3)
    
    # It is easier to interpret an interaction with a graph
    interaction_NDVI_moss <- ggpredict(model3, terms = c("mean_NDVI_121", "moss_prop")) %>% plot() + theme_coding()
    # The interaction effect is for two continuous variables (NDVI and moss prop), but for the sake
    # of visualisation, ggpredict() takes the second continuous variable and splits it into three levels
    # so it's like low, medium and high moss cover.
    
    # Export plot
    png(filename = "plots/interaction_NDVI_moss.png", width = 10, height = 10, units = "cm", res = 400)
    plot(interaction_NDVI_moss)
    dev.off()



# Visualisation ----




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


# Figure XXX. Testing Leaf mass & Phytomass as predictors of biomass ----
# Create plot
(phytomass_biomass <- ggplot(data = dataset,
                             aes(x = phytomass,
                                 y = AGB_spatially_normalised_g_m2)) + 
   geom_point(shape = 1, na.rm = TRUE) +
   theme_coding() +
   coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 600), expand=FALSE) +
   labs(x = expression("Phytomass (g m"^"-2"*")"),
        y = expression("Biomass (g m"^"-2"*")"),
        title = "Phytomass - Biomass") +
   stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
                formula = y ~ x, na.rm = TRUE, coef.digits = 4, rr.digits = 2, size = 3, parse = TRUE,
                label.x.npc = 0.6, label.y.npc = 0.95) +
   geom_smooth(method="lm", formula= y ~ x, se=TRUE, size=0.5, na.rm = TRUE))

# Create plot
(leafmass_biomass <- ggplot(data = dataset,
                            aes(x = leaf_biomass,
                                y = AGB_spatially_normalised_g_m2)) + 
    geom_point(shape = 1, na.rm = TRUE) +
    theme_coding() +
    coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 150), expand=FALSE) +
    labs(x = expression("Leaf biomass (g m"^"-2"*")"),
         y = expression("Biomass (g m"^"-2"*")"),
         title = "Leaf biomass - Biomass") +
    stat_poly_eq(aes(label = paste("atop(", ..eq.label.., ",", ..rr.label.., ")", sep="")),
                 formula = y ~ x, na.rm = TRUE, coef.digits = 4, rr.digits = 2, size = 3, parse = TRUE,
                 label.x.npc = 0.05, label.y.npc = 0.95) +
    geom_smooth(method="lm", formula= y ~ x, se=TRUE, size=0.5, na.rm = TRUE))

# Combine plots
combined_biomass_parts <- ggpubr::ggarrange(phytomass_biomass, leafmass_biomass,
                                            heights = c(10, 10),
                                            labels = c("(a)", "(b)"),
                                            ncol = 2, nrow = 1,
                                            align = "h")


# Export figure
png(filename="plots/Figure XXX - leaf and phytomass vesus biomass.png", width=14, height=7, units="cm", res=400)
plot(combined_biomass_parts)
dev.off()




# Figure XXX. Comparison of NDVI observations  ----
# Consider using geom_errorline to add error bars indicating the maximum and minimum NDVI values observed...

# Create barplot
(NDVI_barplot <- ggplot(data = NDVI_data_long,
                        aes(x = reorder(PlotID, mean_NDVI, FUN = "median"),
                            y = mean_NDVI,
                            fill = grain)) +
   geom_col(aes(y=mean_NDVI, fill = grain), width=0.75, position = position_dodge2(preserve = "single")) +
   scale_fill_grey(start = 0.2, end = 0.8, aesthetics = "fill", guide = guide_legend(direction = "horizontal")) +
   # scale_y_continuous(lim = c(0, 0.9)) +
   coord_cartesian(ylim = c(0, 0.9), expand=FALSE) +
   labs(x = "Plot ID", 
        y = expression("mean NDVI"),
        title = "Comparison of mean NDVI by grain") +
   theme_coding() +
   theme(legend.position = c(0.19, 0.96), 
         axis.line.x = element_line(color="black", size = 0.5),
         axis.line.y = element_line(color="black", size = 0.5),
         axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5))
)

# Export barplot
png(filename = "plots/Figure XXX - NDVI barplot.png", width = 16, height = 20, units = "cm", res = 400)
plot(NDVI_barplot)
dev.off()


