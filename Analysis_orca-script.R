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
library(patchwork)                                                              # For arranging multi-panel plots
library(modelr)
library(rasterVis)                                                              # For visualising rasters
library(gridExtra)                                                              # for arranging pulti-panel raster plot
library(colorspace)                                                             # Generate colour ramps with the colorspace package


# Plotting theme
theme_fancy <- function() {
  theme_bw() +
    theme(
      text = element_text(family = "Helvetica"),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      axis.line.x = element_line(size = 0.3, color = "black"),
      axis.line.y = element_line(size = 0.3, color = "black"),
      axis.ticks = element_line(size = 0.3, color = "black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
      plot.title = element_text(
        size = 10,
        vjust = 1,
        hjust = 0.5,
        color = "black"
      ),
      legend.text = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 10, color = "black"),
      legend.position = c(0.9, 0.9),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 2,
        linetype = "blank"
      )
    )
}

#### Extract NDVI values from rasters ----
#   # Load data for raster extraction
#   feature_filename <- paste0(home, "data/20160725_AC_ORC - formated for exact extractr.geojson")
#   plots <- st_read(feature_filename, crs = 32607)                               # Import geoJSON as sf object, using st_read to allow the non-standard CRS to be specified.
#
#
#   # Import rasters (Available from the NERC Polar Data Centre - see readme for details)
#   raster_018 <- raster(paste0(home, "inputs_NDVI/NDVI_019m_from_20160726.tif"))
#   raster_047 <- raster(paste0(home, "inputs_NDVI/NDVI_050m_from_20160730.tif"))
#   raster_119 <- raster(paste0(home, "inputs_NDVI/NDVI_120m_from_20160730.tif"))
#   raster_121 <- raster(paste0(home, "inputs_NDVI/NDVI_120m_from_20160726.tif"))
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

PF_observations <- read.csv("data/point_framing_observations.csv")            # Read in canopy height from point intercept

DTM_accuracy <- read.csv("data/DTM_accuracy_assessment.csv")                    # Read in DTM accuracy assessment data.


# Compute averge NDVI across the four rasters
dataset$NDVImeans <-
  rowMeans(subset(
    dataset,
    select = c(mean_NDVI_018, mean_NDVI_047, mean_NDVI_119, mean_NDVI_121),
    na.rm = TRUE
  ))
dataset$NDVImedians <-
  rowMedians(subset(
    dataset,
    select = c(mean_NDVI_018, mean_NDVI_047, mean_NDVI_119, mean_NDVI_121),
    na.rm = TRUE
  ))

# Extract solar angles for the surveys
datetime_UTC <- as.character(loc_time_data$datetime_UTC)
SunPosition <-
  getSunlightPosition(date = datetime_UTC,
                      lat = 69.57133,
                      lon = -138.8909)

# Function to convert from radians to degrees
rad2deg <- function(rad) {
  (rad * 180) / (pi)
}

SunPosition$altitude_degrees <- rad2deg(SunPosition$altitude)
SunPosition


#### Data Preparation ----
# Point intercept observations
# Pointframe data - Extracting only Height and PlotN from the pointframe data, and omit NAs
PF_observations2 <-
  dplyr::select(PF_observations, Height, PlotN) %>%
  na.omit(PF_observations) %>%
  mutate(Height = Height / 100) # Convert heights from cm into m. Could be integrated into the above pipe.

# Taxanomic ID
unique_species <- unique(PF_observations$Species)
unique_species

# Compute canopy height summary metrics (min, max, median, mean and quantiles)
# for each plot.
PF_HAG_summary <- PF_observations2 %>%
  dplyr::select(Height, PlotN) %>%
  group_by(PlotN) %>%
  mutate(
    mean = mean(Height, na.rm = TRUE),
    min = min(Height, na.rm = TRUE),
    lower_IQR = quantile(Height, 0.25, na.rm = TRUE),
    median = median(Height, na.rm = TRUE),
    upper_ICR = quantile(Height, 0.75, na.rm = TRUE),
    max = max(Height, na.rm = TRUE)
  )

# Return only unique row summarising canopy heights for each plot.
PF_HAG_summary <-
  PF_HAG_summary[!duplicated(PF_HAG_summary$PlotN),]

# Remove unwanted/erronious individual height.
PF_HAG_summary <- subset(PF_HAG_summary, select = -Height)

# Add point intercept heights into main dataframe
dataset$PF_HAG_min <- PF_HAG_summary$min
dataset$PF_HAG_max <- PF_HAG_summary$max
dataset$PF_HAG_median <- PF_HAG_summary$median
dataset$PF_HAG_mean <- PF_HAG_summary$mean
dataset$PF_HAG_upper_ICR <- PF_HAG_summary$upper_ICR
dataset$PF_HAG_lower_IQR <- PF_HAG_summary$lower_IQR

# Biomass observations
# Calculate spatially normalised values for photosynthetic tissue
dataset$leaf_biomass <-
  dataset$Comp2_Mass / dataset$plot_area_from_field_length  # Leaves.
dataset$herbacious_biomass <-
  dataset$Comp3_Mass / dataset$plot_area_from_field_length  # Herbacious
dataset$phytomass <-
  dataset$leaf_biomass + dataset$herbacious_biomass

# Create long form dataset for plotting canopy heights in SI
PlotID <- rep(dataset$PlotID, times = 2)
method <- rep(c('SfM', 'Point Framing'), each = 36)
min <- c(dataset$HAG_plotmin_of_cellmax_m, dataset$PF_HAG_min)
lower <-
  c(dataset$HAG_plot25percentile_of_cellmax_m,
    dataset$PF_HAG_lower_IQR)
median <-
  c(dataset$HAG_plotmedian_of_cellmax_m, dataset$PF_HAG_median)
upper <-
  c(dataset$HAG_plot75percentile_of_cellmax_m,
    dataset$PF_HAG_upper_ICR)
max <- c(dataset$HAG_plotmax_of_cellmax_m, dataset$PF_HAG_max)
dataset_long <-
  data.frame(PlotID, method, min, lower, median, upper, max)  # Create vectors into new dataframe
rm(PlotID, method, min, lower, median, upper, max)

# Create long form dataset for NDVI boxplot
# NB. It is not straightforward to extract summary statistics like IQR for
# the NDVI rasters, because these metrics are not readily accessible via the
# exactextractr package. Mean, mode, max and min are available.
PlotID <- rep(dataset$PlotID, times = 4)
grain <- rep(c('0.018', '0.047', '0.119', '0.121'), each = 36)
mean_NDVI <- c(
  dataset$mean_NDVI_018,
  dataset$mean_NDVI_047,
  dataset$mean_NDVI_119,
  dataset$mean_NDVI_121
)
NDVI_data_long <-
  data.frame(PlotID, grain, mean_NDVI)  # Create vectors into new dataframe
rm(PlotID, grain, mean_NDVI)  # Tidy up

### Summarise moss cover within each plot ###
# Count number of locations sampled within each plot
PF_count <- PF_observations %>%
  dplyr::select(-Species,-Status,-Tissue,-Count,-Height) %>%
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
PF_moss[nrow(PF_moss) + 1, ] = list("P04", 0)
PF_moss[nrow(PF_moss) + 1, ] = list("P07", 0)
PF_moss[nrow(PF_moss) + 1, ] = list("P11", 0)
PF_moss[nrow(PF_moss) + 1, ] = list("P34", 0)
PF_moss <- dplyr::arrange(PF_moss, PlotN)  # re-order plots

# Determine the proportion of PF locations that were moss for each plot
PF_moss_summary <- PF_count
PF_moss_summary$moss_obs <- PF_moss$moss_obs
PF_moss_summary$moss_prop <-
  PF_moss_summary$moss_obs / PF_moss_summary$PF_obs

# Add moss proportion to main dataframe
dataset$moss_prop <- PF_moss_summary$moss_prop

### Extract distributions of canopy heights from SfM rasters ###
raster_list <-
  list.files(path = "data/Height_rasters/", full.names = TRUE)   # List of elevation rasters

SFM_heights <-
  data.frame()                                                 # Create dataframe

# Extract height values
for (file in raster_list) {
  PlotN <-
    substr(file, 37, 39)                                             # Extract PlotN from file
  
  # Read raster
  height_rast <-
    raster(file, band = 2)  # NB. Band 2 specifies the maximum height recorded in each cell.
  
  # # Plot the distribution of values in the raster
  # hist(height_rast, main=paste0("Distribution of heights in plot ", PlotN, sep=""),
  #      col= "purple",
  #      maxpixels=10000,
  #      xlab = "height (m)",
  #      xlim = c(0, 1.2))
  
  heights <- getValues(height_rast)  # Summarise height values
  
  heights <- na.exclude(heights)  # Remove NAs.
  
  height_distributions <-
    data.frame(heights)  # Convert to dataframe
  
  height_distributions$PlotN <- PlotN  # Add PlotN
  
  SFM_heights <-
    rbind(SFM_heights, height_distributions)  # Integrate heights to dataframe.
  
  rm(heights, height_distributions)  # Tidy up
}


#### Scaling parameters ----

## Set scaling parameters
max_agb <-
  1.1 * max(dataset$AGB_spatially_normalised_g_m2, na.rm = TRUE)
max_hag <- 0.85
max_mean_hag <-
  1.1 * max(max(dataset$HAG_plotmean_of_cellmax_m, na.rm = TRUE),
            max(dataset$PF_HAG_mean))
max_ndvi <- 0.91
min_ndvi <- 0.6
phyto_biomass_max <- 1.12 * max(dataset$phytomass)
spacing <- 2
# Set parameters for the log plots.
phyto_biomass_max_log <- ceiling(max(log(dataset$phytomass))) * 1.1
max_agb_log <-
  ceiling(max(log(dataset$AGB_spatially_normalised_g_m2))) * 1.1
phyto_biomass_min_log <- floor(min(log(dataset$phytomass))) * 0.9
min_agb_log <-
  floor(min(log(dataset$AGB_spatially_normalised_g_m2))) * 0.9


#### Comparison of canopy height measurements ----
# Analysis
# Concordence correlation coefficient
ccc.test <- CCC(
  dataset$HAG_plotmean_of_cellmax_m,
  dataset$PF_HAG_mean,
  ci = "z-transform",
  conf.level = 0.95,
  na.rm = FALSE
)

ccc.test$rho.c                                                            # The concordance correlation coefficient

# Fit model
model_heights_pow <-
  nls(
    HAG_plotmean_of_cellmax_m ~ a * PF_HAG_mean ^ b,
    data = dataset,
    start = list(a = 1, b = 1),
    na.action = na.exclude
  )
summary(model_heights_pow)


# Summarise bias
HAG_rediduals <-
  dataset$HAG_plotmean_of_cellmax_m - dataset$PF_HAG_mean
HAG_bias_mean <- round(mean(HAG_rediduals), 3)
HAG_bias_mean

HAG_bias_mean_SD <- round(sd(HAG_rediduals), 3)
HAG_bias_mean_SD


# Canopy height mixed effects model

# Create dataset
height_data <- dataset %>%
  dplyr::select(contributor_assigned_plot_id, HAG_plotmean_of_cellmax_m) %>%
  rename(PlotN = contributor_assigned_plot_id) %>%
  inner_join(PF_observations2) %>%
  na.omit()

# Plot raw data
plot(HAG_plotmean_of_cellmax_m ~ Height, data = height_data)

# Fit model
model_heights_mixed <-
  glmer(
    HAG_plotmean_of_cellmax_m ~ Height + (1|PlotN),
    family = Gamma(link = "log"),
    data = height_data
  )

# Warnings - model  didn't converge properly
summary(model_heights_mixed)

# Predictions plot flat line - probably because model didn't converge
ggpredict(model_heights_mixed, terms = "Height") %>% plot()


# Plot New Figure 2 - using exponential line from below
(
  Canopy_heights_mixed_plot <- ggplot(data = dataset,
                                      aes(x = PF_HAG_mean, y = HAG_plotmean_of_cellmax_m)) +
    geom_point(shape = 1, colour = "white") +
    geom_point(
      data = height_data,
      aes(x = Height, y = HAG_plotmean_of_cellmax_m),
      shape = 1,
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    labs(x = "Point Frame - Canopy Height (m)",
         y = "SfM - Canopy height (m)") +
    theme_fancy() +
    coord_cartesian(
      ylim = c(0, 1),
      xlim = c(0, 1),
      expand = FALSE
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed"
    ) +
    stat_function(
      fun = function(x)
        (coef(summary(model_heights_pow))[, "Estimate"])[1] * (x) ^ (coef(summary(model_heights_pow))[, "Estimate"])[2],
      aes(),
      size = 1,
      lty = "solid",
      colour = "Black"
    ) +
    annotate(
      "text",
      x = 0.8,
      y = 0.1,
      label = "italic(Y)==1.085~italic(X)^0.715",
      parse = TRUE,
      color = "black",
      size = 3,
      family = "sans"
    )
)

ggsave(
  Canopy_heights_mixed_plot,
  filename = "plots/Figure 2 - Canopy Heights - alternative.pdf",
  width = 10,
  height = 10,
  units = "cm"
)

ggsave(
  Canopy_heights_mixed_plot,
  filename = "plots/Figure 2 - Canopy Heights - alternative.png",
  width = 10,
  height = 10,
  units = "cm"
)


# Visualistation (Figure 2)
# mean point intercept canopy height versus mean structure-from-motion canopy height
# Create plot
(
  Canopy_heights_plot <- ggplot(data = dataset,
                                aes(x = PF_HAG_mean, y = HAG_plotmean_of_cellmax_m)) +
    geom_point(shape = 1) +
    labs(x = "Point Frame - Canopy Height (m)",
         y = "SfM - Canopy height (m)") +
    theme_fancy() +
    coord_cartesian(
      ylim = c(0, 1),
      xlim = c(0, 1),
      expand = FALSE
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed"
    ) +
    stat_poly_eq(
      aes(label = paste(
        "atop(", ..eq.label.., ",", ..rr.label.., ")", sep = ""
      )),
      formula = y ~ x,
      na.rm = TRUE,
      coef.digits = 4,
      rr.digits = 3,
      size = 3,
      parse = TRUE,
      label.x.npc = 0.90,
      label.y.npc = 0.10
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      size = 0.5,
      na.rm = TRUE
    )
)


(
  Canopy_heights_plot <- ggplot(data = dataset,
                                aes(x = PF_HAG_mean, y = HAG_plotmean_of_cellmax_m)) +
    geom_point(shape = 1) +
    labs(x = "Point Frame - Canopy Height (m)",
         y = "SfM - Canopy height (m)") +
    theme_fancy() +
    coord_cartesian(
      ylim = c(0, 1),
      xlim = c(0, 1),
      expand = FALSE
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed"
    ) +
    stat_function(
      fun = function(x)
        (coef(summary(model_heights_pow))[, "Estimate"])[1] * (x) ^ (coef(summary(model_heights_pow))[, "Estimate"])[2],
      aes(),
      size = 1,
      lty = "solid",
      colour = "Black"
    ) +
    annotate(
      "text",
      x = 0.8,
      y = 0.1,
      label = "italic(Y)==1.085~italic(X)^0.715",
      parse = TRUE,
      color = "black",
      size = 3,
      family = "sans"
    )
)

# Export plot
ggsave(
  Canopy_heights_plot,
  filename = "plots/Figure 2 - Canopy Heights.pdf",
  width = 10,
  height = 10,
  units = "cm"
)

ggsave(
  Canopy_heights_plot,
  filename = "plots/Figure 2 - Canopy Heights.png",
  width = 10,
  height = 10,
  units = "cm"
)


#### Predictors of biomass (height and NDVI) ----
# Testing canopy height (from point intercept and photogrammetry) and NDVI as predictors of biomass.
# Analysis
# Point intercept
model_PFu <-
  lm(AGB_spatially_normalised_g_m2 ~ PF_HAG_mean, data = dataset)  # Unconstrained intercept
model_PF <-
  lm(AGB_spatially_normalised_g_m2 ~ PF_HAG_mean + 0, data = dataset)  # constrained intercept
model_SfMu <-
  lm(AGB_spatially_normalised_g_m2 ~ HAG_plotmean_of_cellmax_m,
     data = dataset)  # unconstrained intercept
model_SfM <-
  lm(AGB_spatially_normalised_g_m2 ~ HAG_plotmean_of_cellmax_m + 0,
     data = dataset)  # constrained intercept

summary(model_PFu)
summary(model_SfMu)


# Tabulate and export model parameters
model_list <- list(model_PFu)
model_results1 <- bind_rows(lapply(model_list, function(model) {
  model_tidy <- broom::tidy(model)
  model_argument <- broom::augment(model)
  model_glance <- broom::glance(model)
  return(
    data.frame(
      Canopy_height_method = "Point Framing",
      model_form = "Y = a X + b",
      a = paste0(
        round(model_tidy$estimate[2], 1),
        " ± ",
        round(model_tidy$std.error[2], 1)
      ),
      b = paste0(
        round(model_tidy$estimate[1], 1),
        " ± ",
        round(model_tidy$std.error[1], 1)
      ),
      r2 = round(model_glance$r.squared, 2),
      RMSE = round(sqrt(c(crossprod(model_argument$.resid)) / length(model_argument$.resid)), 1),
      # resid = round(model_glance$sigma, 3),
      stringsAsFactors = F
    )
  )
}))


model_list <- list(model_PF)
model_results2 <- bind_rows(lapply(model_list, function(model) {
  model_tidy <- broom::tidy(model)
  model_argument <- broom::augment(model)
  model_glance <- broom::glance(model)
  return(
    data.frame(
      Canopy_height_method = "Point Framing",
      model_form = "Y = a X",
      a = paste0(
        round(model_tidy$estimate[1], 1),
        " ± ",
        round(model_tidy$std.error[1], 1)
      ),
      b = paste0("0"),
      r2 = round(model_glance$r.squared, 2),
      RMSE = round(sqrt(c(crossprod(model_argument$.resid)) / length(model_argument$.resid)), 1),
      # resid = round(model_glance$sigma, 3),
      stringsAsFactors = F
    )
  )
}))


model_list <- list(model_SfMu)
model_results3 <- bind_rows(lapply(model_list, function(model) {
  model_tidy <- broom::tidy(model)
  model_argument <- broom::augment(model)
  model_glance <- broom::glance(model)
  return(
    data.frame(
      Canopy_height_method = "Photogrammetry",
      model_form = "Y = a X + b",
      a = paste0(
        round(model_tidy$estimate[2], 1),
        " ± ",
        round(model_tidy$std.error[2], 1)
      ),
      b = paste0(
        round(model_tidy$estimate[1], 1),
        " ± ",
        round(model_tidy$std.error[1], 1)
      ),
      r2 = round(model_glance$r.squared, 2),
      RMSE = round(sqrt(c(crossprod(model_argument$.resid)) / length(model_argument$.resid)), 1),
      # resid = round(model_glance$sigma, 3),
      stringsAsFactors = F
    )
  )
}))


model_list <- list(model_SfM)
model_results4 <- bind_rows(lapply(model_list, function(model) {
  model_tidy <- broom::tidy(model)
  model_argument <- broom::augment(model)
  model_glance <- broom::glance(model)
  return(
    data.frame(
      Canopy_height_method = "Photogrammetry",
      model_form = "Y = a X",
      a = paste0(
        round(model_tidy$estimate[1], 1),
        " ± ",
        round(model_tidy$std.error[1], 1)
      ),
      b = paste0("0"),
      r2 = round(model_glance$r.squared, 2),
      RMSE = round(sqrt(c(crossprod(model_argument$.resid)) / length(model_argument$.resid)), 1),
      # resid = round(model_glance$sigma, 3),
      stringsAsFactors = F
    )
  )
}))

# Combine tables
blended <- rbind(model_results1,
                 model_results2,
                 model_results3,
                 model_results4)

# Export model parameters to table
write.csv(blended, file = "tables/Table S1 model fits.csv", row.names = F)


# Visualisation

# Create plots
(
  biomass_CH_SfM <- ggplot(
    data = dataset,
    aes(x = HAG_plotmean_of_cellmax_m,
        y = AGB_spatially_normalised_g_m2)
  ) +
    geom_point(shape = 1, na.rm = TRUE) +
    theme_fancy() +
    coord_cartesian(
      ylim = c(0, 3000),
      xlim = c(0, max_hag),
      expand = FALSE
    ) +
    labs(
      x = expression("Canopy height (m)"),
      y = expression("Dry biomass (g m" ^ "-2" * ")"),
      title = "SfM"
    ) +
    annotate(geom='text', x = (0.04*max_hag), y = (0.97*3000), size =2.6, label = 'y = 2522 x', hjust = 0) +
    annotate(geom='text', x = (0.04*max_hag), y = (0.87*3000), size =2.6, label = paste('R^2 == 0.90'), hjust = 0, parse = TRUE) +
    theme(legend.position = c(0.15, 0.9)) +
    geom_smooth(
      method = "lm",
      formula = y ~ x - 1,
      se = TRUE,
      size = 0.5,
      color = "black",
      na.rm = TRUE
    )
)


(
  biomass_CH_PF <- ggplot(data = dataset,
                          aes(x = PF_HAG_mean,
                              y = AGB_spatially_normalised_g_m2)) +
    geom_point(shape = 1, na.rm = TRUE) +
    theme_fancy() +
    coord_cartesian(
      ylim = c(0, 3000),
      xlim = c(0, max_hag),
      expand = FALSE
    ) +
    labs(
      x = expression("Canopy height (m)"),
      y = expression("Dry biomass (g m" ^ "-2" * ")"),
      title = "Point Intercept"
    ) +
    annotate(geom='text', x = (0.04*max_hag), y = (0.97*3000), size =2.6, label = 'y = 3623 x', hjust = 0) +
    annotate(geom='text', x = (0.04*max_hag), y = (0.87*3000), size =2.6, label = paste('R^2 == 0.92'), hjust = 0, parse = TRUE) +
    theme(legend.position = c(0.15, 0.9)) +
    geom_smooth(
      method = "lm",
      formula = y ~ x - 1,
      se = TRUE,
      size = 0.5,
      color = "black",
      na.rm = TRUE
    )
)

(
  biomass_NDVI <- ggplot(data = dataset,
                         aes(
                           x = mean_NDVI_119,
                           y = log(AGB_spatially_normalised_g_m2)
                         )) +
    geom_point(shape = 1, na.rm = TRUE) +
    theme_fancy() +
    coord_cartesian(
      ylim = c(4.5, 9),
      xlim = c(0.6, 0.8),
      expand = FALSE
    ) +
    labs(
      x = expression("NDVI (0.119 m grain)"),
      y = expression("ln(Dry biomass (g m" ^ "-2" * "))"),
      title = "NDVI"
    ) +
    annotate(geom='text', x = 0.66, y = 8.85, size =2.6, label = 'y = -0.372 +9.902 x', hjust = 0) +
    annotate(geom='text', x = 0.66, y = 8.35, size =2.6, label = paste('R^2 == 0.23'), hjust = 0, parse = TRUE) +
    theme(legend.position = c(0.15, 0.9)) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      size = 0.5,
      color = "black",
      na.rm = TRUE
    )
)


# Combine plots
biomass_plots <-
  ggpubr::ggarrange(
    biomass_CH_PF,
    biomass_CH_SfM,
    biomass_NDVI,
    heights = c(10, 10, 10),
    labels = c("(a)", "(b)", "(c)"),
    font.label = list(size = 10, face = "bold"),
    ncol = 3,
    nrow = 1,
    align = "h"
  )

# Export figure
ggsave(
  biomass_plots,
  filename = "plots/Figure 3 - Biomass Predictions.pdf",
  width = 20,
  height = 7,
  units = "cm"
)

ggsave(
  biomass_plots,
  filename = "plots/Figure 3 - Biomass Predictions.png",
  width = 20,
  height = 7,
  units = "cm"
)


### Analysis of NDVI - biomass relationships ----
# Testing relationships between NDVI and various biomass components

### Exponential models fitted to  biomass values ###
# Analysis
# Exponential models
exp_model_total_NDVI_121 <-
  nls(
    AGB_spatially_normalised_g_m2 ~ a * exp(b * mean_NDVI_121),
    data = dataset,
    start = list(a = 30, b = 4),
    na.action = na.exclude
  )
exp_model_total_NDVI_119 <-
  nls(
    AGB_spatially_normalised_g_m2 ~ a * exp(b * mean_NDVI_119),
    data = dataset,
    start = list(a = 30, b = 4),
    na.action = na.exclude
  )
exp_model_total_NDVI_047 <-
  nls(
    AGB_spatially_normalised_g_m2 ~ a * exp(b * mean_NDVI_047),
    data = dataset,
    start = list(a = 30, b = 4),
    na.action = na.exclude
  )
exp_model_total_NDVI_018 <-
  nls(
    AGB_spatially_normalised_g_m2 ~ a * exp(b * mean_NDVI_018),
    data = dataset,
    start = list(a = 30, b = 4),
    na.action = na.exclude
  )
exp_model_phyto_NDVI_121 <-
  nls(
    phytomass ~ a * exp(b * mean_NDVI_121),
    data = dataset,
    start = list(a = 30, b = 1),
    na.action = na.exclude
  )
exp_model_phyto_NDVI_119 <-
  nls(
    phytomass ~ a * exp(b * mean_NDVI_119),
    data = dataset,
    start = list(a = 30, b = 1),
    na.action = na.exclude
  )
exp_model_phyto_NDVI_047 <-
  nls(
    phytomass ~ a * exp(b * mean_NDVI_047),
    data = dataset,
    start = list(a = 30, b = 1),
    na.action = na.exclude
  )
exp_model_phyto_NDVI_018 <-
  nls(
    phytomass ~ a * exp(b * mean_NDVI_018),
    data = dataset,
    start = list(a = 30, b = 1),
    na.action = na.exclude
  )
exp_model_leaf_NDVI_121 <-
  nls(
    leaf_biomass ~ a * exp(b * mean_NDVI_121),
    data = dataset,
    start = list(a = 30, b = 1),
    na.action = na.exclude
  )
exp_model_leaf_NDVI_119 <-
  nls(
    leaf_biomass ~ a * exp(b * mean_NDVI_119),
    data = dataset,
    start = list(a = 30, b = 1),
    na.action = na.exclude
  )
exp_model_leaf_NDVI_047 <-
  nls(
    leaf_biomass ~ a * exp(b * mean_NDVI_047),
    data = dataset,
    start = list(a = 30, b = 1),
    na.action = na.exclude
  )
exp_model_leaf_NDVI_018 <-
  nls(
    leaf_biomass ~ a * exp(b * mean_NDVI_018),
    data = dataset,
    start = list(a = 30, b = 1),
    na.action = na.exclude
  )

# Compile model objects
exp_models <- list(
  exp_model_total_NDVI_121,
  exp_model_total_NDVI_119,
  exp_model_total_NDVI_047,
  exp_model_total_NDVI_018,
  exp_model_phyto_NDVI_121,
  exp_model_phyto_NDVI_119,
  exp_model_phyto_NDVI_047,
  exp_model_phyto_NDVI_018,
  exp_model_leaf_NDVI_121,
  exp_model_leaf_NDVI_119,
  exp_model_leaf_NDVI_047,
  exp_model_leaf_NDVI_018
)

summary(exp_model_total_NDVI_121)

# Tabulate model parameters
exp_model_results <- bind_rows(lapply(exp_models, function(model) {
  model_glance <- broom::glance(model)
  model_tidy <- broom::tidy(model)
  return(
    data.frame(
      ndvi_grain = gsub("mean_NDVI_", "0\\.", model_tidy$term[2]),
      model_form = "Y = a e(b X)",
      a = paste0(
        round(model_tidy$estimate[2], 3),
        " ± ",
        round(model_tidy$std.error[2], 3)
      ),
      b = paste0(
        round(model_tidy$estimate[1], 3),
        " ± ",
        round(model_tidy$std.error[1], 3)
      ),
      resid = round(model_glance$sigma, 3),
      stringsAsFactors = F
    )
  )
}))



# Tabulate model parameters
exp_model_results <- bind_rows(lapply(exp_models, function(model) {
  model_tidy <- broom::tidy(model)
  model_argument <- broom::augment(model)
  model_glance <- broom::glance(model)
  return(
    data.frame(
      ndvi_grain = gsub("mean_NDVI_", "0\\.", model_tidy$term[2]),
      model_form = "Y = a e(b X)",
      a = paste0(
        round(model_tidy$estimate[2], 3),
        " ± ",
        round(model_tidy$std.error[2], 3)
      ),
      b = paste0(
        round(model_tidy$estimate[1], 3),
        " ± ",
        round(model_tidy$std.error[1], 3)
      ),
      RMSE = round(sqrt(c(crossprod(model_argument$.resid)) / length(model_argument$.resid)), 1),
      # resid = round(model_glance$sigma, 3),
      stringsAsFactors = F
    )
  )
}))

# Export model parameters to table
write.csv(exp_model_results, file = "tables/Table S2 exp model fits.csv", row.names = F)


# Visualisation
# Create plots
# Total biomass
NDVI_vs_total_biomass_121 <-
  ggplot(data = dataset, aes(x = mean_NDVI_121, y = AGB_spatially_normalised_g_m2)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(# x = expression("mean NDVI"),
    x = expression(""),
    y = expression(atop("Total", paste (
      "biomass (g m" ^ "-2" * ")"
    ))),
    title = "0.121 m grain") +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_total_NDVI_121
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_total_NDVI_121)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, max_agb),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_total_biomass_119 <-
  ggplot(data = dataset, aes(x = mean_NDVI_119, y = AGB_spatially_normalised_g_m2)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression(""),
       title = "0.119 m grain") +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_total_NDVI_119
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_total_NDVI_119)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, max_agb),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_total_biomass_047 <-
  ggplot(data = dataset, aes(x = mean_NDVI_047, y = AGB_spatially_normalised_g_m2)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression(""),
       title = "0.047 m grain") +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_total_NDVI_047
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_total_NDVI_047)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, max_agb),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_total_biomass_018 <-
  ggplot(data = dataset, aes(x = mean_NDVI_018,  y = AGB_spatially_normalised_g_m2)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression(""),
       title = "0.018 m grain") +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_total_NDVI_018
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_total_NDVI_018)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, max_agb),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))


# Photosynthetic biomass
NDVI_vs_phyto_biomass_121 <-
  ggplot(data = dataset, aes(x = mean_NDVI_121,  y = phytomass)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression(atop("Phytomass", paste("(g m" ^ "-2" * ")")))) +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_phyto_NDVI_121
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_phyto_NDVI_121)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, phyto_biomass_max),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_phyto_biomass_119 <-
  ggplot(data = dataset, aes(x = mean_NDVI_119,  y = phytomass)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression("")) +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_phyto_NDVI_119
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_phyto_NDVI_119)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, phyto_biomass_max),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_phyto_biomass_047 <-
  ggplot(data = dataset, aes(x = mean_NDVI_047, y = phytomass)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression("")) +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_phyto_NDVI_047
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_phyto_NDVI_047)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, phyto_biomass_max),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_phyto_biomass_018 <-
  ggplot(data = dataset, aes(x = mean_NDVI_018,  y = phytomass)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression("")) +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_phyto_NDVI_018
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_phyto_NDVI_018)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, phyto_biomass_max),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

# Leaf biomass
# Create plots
NDVI_vs_leaf_biomass_121 <-
  ggplot(data = dataset, aes(x = mean_NDVI_121, y = leaf_biomass)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression("NDVI"),
       y = expression(atop("Leaf", paste (
         "biomass (g m" ^ "-2" * ")"
       )))) +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_leaf_NDVI_121
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_leaf_NDVI_121)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, 200),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_leaf_biomass_119 <-
  ggplot(data = dataset, aes(x = mean_NDVI_119, y = leaf_biomass)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression("NDVI"),
       y = expression("")) +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_leaf_NDVI_119
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_leaf_NDVI_119)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, 200),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_leaf_biomass_047 <-
  ggplot(data = dataset, aes(x = mean_NDVI_047, y = leaf_biomass)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression("NDVI"),
       y = expression("")) +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_leaf_NDVI_047
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_leaf_NDVI_047)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, 200),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_leaf_biomass_018 <-
  ggplot(data = dataset, aes(x = mean_NDVI_018, y = leaf_biomass)) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression("NDVI"),
       y = expression("")) +
  stat_function(
    fun = function(x)
      (coef(summary(
        exp_model_leaf_NDVI_018
      ))[, "Estimate"])[1] * exp((coef(
        summary(exp_model_leaf_NDVI_018)
      )[, "Estimate"])[2] * x),
    aes(),
    size = 1,
    lty = "solid"
  ) +
  coord_cartesian(
    ylim = c(0, 200),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))


# Combine plots
NDVI_biomass <- ggpubr::ggarrange(
  NDVI_vs_total_biomass_121,
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
  heights = c(9, 9, 9),
  ncol = 4,
  nrow = 3,
  align = "v",
  labels = c(
    "(a)",
    "(b)",
    "(c)",
    "(d)",
    "(e)",
    "(f)",
    "(g)",
    "(h)",
    "(i)",
    "(j)",
    "(k)",
    "(l)"
  ),
  font.label = list(size = 10, face = "bold")
)

# Export figure
ggsave(
  NDVI_biomass,
  filename = "plots/Figure S5 - NDVI vs biomass - exp.pdf",
  width = 25,
  height = 15,
  units = "cm"
)

ggsave(
  NDVI_biomass,
  filename = "plots/Figure S5 - NDVI vs biomass - exp.png",
  width = 25,
  height = 15,
  units = "cm"
)


### Exploration of average NDVI
(
  NDVImedians_vs_total_biomass <-
    ggplot(data = dataset, aes(
      x = NDVImedians, y = log(AGB_spatially_normalised_g_m2)
    )) +
    geom_point(shape = 1, na.rm = TRUE) +
    geom_text(aes(label = PlotID), hjust = 0, vjust = 0) +
    labs(
      # x = expression("median mean NDVI"),
      x = expression(""),
      y = expression(atop(
        "ln (Total", paste ("biomass (g m" ^ "-2" * "))")
      )),
      title = "median NDVI value across grains"
    ) +
    geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = TRUE,
      size = 1,
      lty = "solid",
      col = "black"
    ) +
    coord_cartesian(
      ylim = c(min_agb_log, max_agb_log),
      xlim = c(min_ndvi, max_ndvi),
      expand = FALSE
    ) +
    theme_fancy() +
    theme(
      plot.margin = margin(
        t = spacing,
        r = spacing,
        b = spacing,
        l = spacing,
        unit = "pt"
      )
    )
)

(
  NDVImedians_vs_leaf_biomass <-
    ggplot(data = dataset, aes(
      x = NDVImedians, y = log(leaf_biomass)
    )) +
    geom_point(shape = 1, na.rm = TRUE) +
    geom_text(aes(label = PlotID), hjust = 0, vjust = 0) +
    labs(x = expression("NDVI"),
         y = expression("")) +
    geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = TRUE,
      size = 1,
      lty = "solid",
      col = "black"
    ) +
    coord_cartesian(
      ylim = c(3, 5.5),
      xlim = c(min_ndvi, max_ndvi),
      expand = FALSE
    ) +
    theme_fancy() +
    theme(
      plot.margin = margin(
        t = spacing,
        r = spacing,
        b = spacing,
        l = spacing,
        unit = "pt"
      )
    )
)




### Linear models fitted to logged biomass values ###
# Analysis
# Logarithmic models
log_model_total_NDVI_121 <- lm(log(AGB_spatially_normalised_g_m2) ~ mean_NDVI_121, data = dataset, na.action = na.exclude)
log_model_total_NDVI_119 <- lm(log(AGB_spatially_normalised_g_m2) ~ mean_NDVI_119, data = dataset, na.action = na.exclude)
log_model_total_NDVI_047 <- lm(log(AGB_spatially_normalised_g_m2) ~ mean_NDVI_047, data = dataset, na.action = na.exclude)
log_model_total_NDVI_018 <- lm(log(AGB_spatially_normalised_g_m2) ~ mean_NDVI_018, data = dataset, na.action = na.exclude)
log_model_phyto_NDVI_121 <- lm(log(phytomass) ~ mean_NDVI_121, data = dataset, na.action = na.exclude)
log_model_phyto_NDVI_119 <- lm(log(phytomass) ~ mean_NDVI_119, data = dataset, na.action = na.exclude) 
log_model_phyto_NDVI_047 <- lm(log(phytomass) ~ mean_NDVI_047, data = dataset, na.action = na.exclude)
log_model_phyto_NDVI_018 <- lm(log(phytomass) ~ mean_NDVI_018, data = dataset, na.action = na.exclude)
log_model_leaf_NDVI_121 <- lm(log(leaf_biomass) ~ mean_NDVI_121, data = dataset, na.action = na.exclude)
log_model_leaf_NDVI_119 <- lm(log(leaf_biomass) ~ mean_NDVI_119, data = dataset, na.action = na.exclude)
log_model_leaf_NDVI_047 <- lm(log(leaf_biomass) ~ mean_NDVI_047, data = dataset, na.action = na.exclude)
log_model_leaf_NDVI_018 <- lm(log(leaf_biomass) ~ mean_NDVI_018, data = dataset, na.action = na.exclude)

# Compile model objects
log_models <- list(
  log_model_total_NDVI_121,
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
  log_model_leaf_NDVI_018
)


# Tabulate model parameters
log_model_results <- bind_rows(lapply(log_models, function(model) {
  model_tidy <- broom::tidy(model)
  model_argument <- broom::augment(model)
  model_glance <- broom::glance(model)
  return(
    data.frame(
      ndvi_grain = gsub("mean_NDVI_", "0\\.", model_tidy$term[2]),
      model_form = "ln(Y) = a X + b",
      a = paste0(
        round(model_tidy$estimate[2], 3),
        " ± ",
        round(model_tidy$std.error[2], 3)
      ),
      b = paste0(
        round(model_tidy$estimate[1], 3),
        " ± ",
        round(model_tidy$std.error[1], 3)
      ),
      r2 = round(model_glance$r.squared, 2),
      RMSE = round(sqrt(c(crossprod(model_argument$.resid)) / length(model_argument$.resid)), 4),
      # resid = round(model_glance$sigma, 3),
      stringsAsFactors = F
    )
  )
}))

# Export model parameters to table
write.csv(log_model_results, file = "tables/Table 2 log model fits.csv", row.names = F)


# Visualisation
# Create plots
# Total biomass

NDVI_vs_total_biomass_121 <-
  ggplot(data = dataset, aes(x = mean_NDVI_121, y = log(AGB_spatially_normalised_g_m2))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(# x = expression("mean NDVI"),
    x = expression(""),
    y = expression(atop(
      "ln (Total", paste ("biomass (g m" ^ "-2" * "))")
    )),
    title = "0.121 m grain") +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(min_agb_log, max_agb_log),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_total_biomass_119 <-
  ggplot(data = dataset, aes(x = mean_NDVI_119, y = log(AGB_spatially_normalised_g_m2))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression(""),
       title = "0.119 m grain") +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(min_agb_log, max_agb_log),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_total_biomass_047 <-
  ggplot(data = dataset, aes(x = mean_NDVI_047, y = log(AGB_spatially_normalised_g_m2))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression(""),
       title = "0.047 m grain") +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(min_agb_log, max_agb_log),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_total_biomass_018 <-
  ggplot(data = dataset, aes(x = mean_NDVI_018,  y = log(AGB_spatially_normalised_g_m2))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression(""),
       title = "0.018 m grain") +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(min_agb_log, max_agb_log),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

# Phytomass biomass
NDVI_vs_phyto_biomass_121 <-
  ggplot(data = dataset, aes(x = mean_NDVI_121,  y = log(phytomass))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression(atop("ln (Phytomass", paste("(g m" ^ "-2" * "))")))) +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(phyto_biomass_min_log, phyto_biomass_max_log),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_phyto_biomass_119 <-
  ggplot(data = dataset, aes(x = mean_NDVI_119,  y = log(phytomass))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression("")) +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(phyto_biomass_min_log, phyto_biomass_max_log),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_phyto_biomass_047 <-
  ggplot(data = dataset, aes(x = mean_NDVI_047, y = log(phytomass))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression("")) +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(phyto_biomass_min_log, phyto_biomass_max_log),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_phyto_biomass_018 <-
  ggplot(data = dataset, aes(x = mean_NDVI_018,  y = log(phytomass))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression(""),
       y = expression("")) +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(phyto_biomass_min_log, phyto_biomass_max_log),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

# Leaf biomass
# Create plots
NDVI_vs_leaf_biomass_121 <-
  ggplot(data = dataset, aes(x = mean_NDVI_121, y = log(leaf_biomass))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression("NDVI"),
       y = expression(atop(
         "ln (Leaf", paste ("biomass (g m" ^ "-2" * "))")
       ))) +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(3, 5.5),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_leaf_biomass_119 <-
  ggplot(data = dataset, aes(x = mean_NDVI_119, y = log(leaf_biomass))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression("NDVI"),
       y = expression("")) +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(3, 5.5),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_leaf_biomass_047 <-
  ggplot(data = dataset, aes(x = mean_NDVI_047, y = log(leaf_biomass))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression("NDVI"),
       y = expression("")) +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(3, 5.5),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))

NDVI_vs_leaf_biomass_018 <-
  ggplot(data = dataset, aes(x = mean_NDVI_018, y = log(leaf_biomass))) +
  geom_point(shape = 1, na.rm = TRUE) +
  labs(x = expression("NDVI"),
       y = expression("")) +
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = TRUE,
    size = 1,
    lty = "solid",
    col = "black"
  ) +
  coord_cartesian(
    ylim = c(3, 5.5),
    xlim = c(min_ndvi, max_ndvi),
    expand = FALSE
  ) +
  theme_fancy() +
  theme(plot.margin = margin(
    t = spacing,
    r = spacing,
    b = spacing,
    l = spacing,
    unit = "pt"
  ))


# Combine plots
NDVI_biomass <- ggpubr::ggarrange(
  NDVI_vs_total_biomass_121,
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
  heights = c(9, 9, 9),
  ncol = 4,
  nrow = 3,
  align = "v",
  labels = c(
    "(a)",
    "(b)",
    "(c)",
    "(d)",
    "(e)",
    "(f)",
    "(g)",
    "(h)",
    "(i)",
    "(j)",
    "(k)",
    "(l)"
  ),
  font.label = list(size = 10, face = "bold")
)

# Export figure
ggsave(
  NDVI_biomass,
  filename = "plots/Figure 4 - NDVI vs biomass - log.pdf",
  width = 25,
  height = 15,
  units = "cm"
)

ggsave(
  NDVI_biomass,
  filename = "plots/Figure 4 - NDVI vs biomass - log.png",
  width = 25,
  height = 15,
  units = "cm"
)


# Figure 5. Leaf mass & Phytomass as predictors of biomass ----
# Create plot
(
  phytomass_biomass <- ggplot(data = dataset,
                              aes(x = phytomass,
                                  y = AGB_spatially_normalised_g_m2)) +
    geom_point(shape = 1, na.rm = TRUE) +
    theme_fancy() +
    coord_cartesian(
      ylim = c(0, max_agb),
      xlim = c(0, 600),
      expand = FALSE
    ) +
    labs(
      x = expression("Phytomass (g m" ^ "-2" * ")"),
      y = expression("Total Biomass (g m" ^ "-2" * ")")
    ) +
    stat_poly_eq(
      aes(label = paste(
        "atop(", ..eq.label.., ",", ..rr.label.., ")", sep = ""
      )),
      formula = y ~ x,
      na.rm = TRUE,
      coef.digits = 4,
      rr.digits = 2,
      size = 3,
      parse = TRUE,
      label.x.npc = 0.95,
      label.y.npc = 0.95
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      size = 0.5,
      na.rm = TRUE,
      colour = "black"
    )
)

# Create plot
(
  leafmass_biomass <- ggplot(data = dataset,
                             aes(x = leaf_biomass,
                                 y = AGB_spatially_normalised_g_m2)) +
    geom_point(shape = 1, na.rm = TRUE) +
    theme_fancy() +
    coord_cartesian(
      ylim = c(0, max_agb),
      xlim = c(0, 150),
      expand = FALSE
    ) +
    labs(
      x = expression("Shrub leaf biomass (g m" ^ "-2" * ")"),
      y = expression("Total Biomass (g m" ^ "-2" * ")")
    ) +
    stat_poly_eq(
      aes(label = paste(
        "atop(", ..eq.label.., ",", ..rr.label.., ")", sep = ""
      )),
      formula = y ~ x,
      na.rm = TRUE,
      coef.digits = 4,
      rr.digits = 2,
      size = 3,
      parse = TRUE,
      label.x.npc = 0.05,
      label.y.npc = 0.95
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      size = 0.5,
      na.rm = TRUE,
      colour = "black"
    )
)

# Combine plots
  combined_biomass_parts <- ggpubr::ggarrange(leafmass_biomass, phytomass_biomass,
                                              heights = c(10, 10),
                                              labels = c("(a)", "(b)"),
                                              font.label = list(size = 10, face = "bold"),
                                              ncol = 2, nrow = 1,
                                              align = "h")

# Export figure
  ggsave(
    combined_biomass_parts,
    filename = "plots/Figure 5 - leaf and phytomass vesus biomass.pdf",
    width = 18,
    height = 9,
    units = "cm"
  )
  
  ggsave(
    combined_biomass_parts,
    filename = "plots/Figure 5 - leaf and phytomass vesus biomass.png",
    width = 18,
    height = 9,
    units = "cm"
  )
  

dataset$herb_prop <- dataset$herbacious_biomas / dataset$phytomass
hist(dataset$herb_prop)
summary(dataset$herb_prop)


### Figure 6. Analysis of moss cover effect on NDVI-biomass relationships ####
# Looking at the inteaction between the proportion of moss cover ('moss_prop') with NDVI and Biomass.
# The coefficient of that interaction effect indicates how moss_prop influences the phytomass/NDVI relationship.

# Check the distribution of phytomass. Suggests phytomass should be transformed to normalise the distribution. But that leaf and total biomass are better without normalisation.
hist(dataset$moss_prop)
hist(dataset$phytomass)
hist(log(dataset$phytomass))
hist(dataset$leaf_biomass)
hist(log(dataset$leaf_biomass))
hist(dataset$AGB_spatially_normalised_g_m2)
hist(log(dataset$AGB_spatially_normalised_g_m2))

# Create models
mod_biomass_121 <- lm(AGB_spatially_normalised_g_m2 ~ mean_NDVI_121 * moss_prop, data = dataset)
mod_biomass_119 <- lm(AGB_spatially_normalised_g_m2 ~ mean_NDVI_119 * moss_prop, data = dataset)
mod_biomass_047 <- lm(AGB_spatially_normalised_g_m2 ~ mean_NDVI_047 * moss_prop, data = dataset)
mod_biomass_018 <- lm(AGB_spatially_normalised_g_m2 ~ mean_NDVI_018 * moss_prop, data = dataset)

mod_phytomass_121 <- lm(phytomass ~ mean_NDVI_121 * moss_prop, data = dataset)
mod_phytomass_119 <- lm(phytomass ~ mean_NDVI_119 * moss_prop, data = dataset)
mod_phytomass_047 <- lm(phytomass ~ mean_NDVI_047 * moss_prop, data = dataset)
mod_phytomass_018 <- lm(phytomass ~ mean_NDVI_018 * moss_prop, data = dataset)

mod_leafmass_121 <- lm(leaf_biomass ~ mean_NDVI_121 * moss_prop, data = dataset)
mod_leafmass_119 <- lm(leaf_biomass ~ mean_NDVI_119 * moss_prop, data = dataset)
mod_leafmass_047 <- lm(leaf_biomass ~ mean_NDVI_047 * moss_prop, data = dataset)
mod_leafmass_018 <- lm(leaf_biomass ~ mean_NDVI_018 * moss_prop, data = dataset)


summary(mod_biomass_121)
summary(mod_biomass_119)
summary(mod_biomass_047)
summary(mod_biomass_018)

summary(mod_phytomass_121)
summary(mod_phytomass_119)
summary(mod_phytomass_047)
summary(mod_phytomass_018)

summary(mod_leafmass_121)
summary(mod_leafmass_119)
summary(mod_leafmass_047)
summary(mod_leafmass_018)


# Compile model objects
moss_models <- list(
  mod_biomass_121,
  mod_biomass_119,
  mod_biomass_047,
  mod_biomass_018,
  mod_phytomass_121,
  mod_phytomass_119,
  mod_phytomass_047,
  mod_phytomass_018,
  mod_leafmass_121,
  mod_leafmass_119,
  mod_leafmass_047,
  mod_leafmass_018
)


# Tabulate interaction model parameters
moss_model_results <- bind_rows(lapply(moss_models, function(model) {
  model_glance <- broom::glance(model)
  model_tidy <- broom::tidy(model)
  return(
    data.frame(
      dependent_variable = as.character(formula(model)[2]),
      NDVI_Grain_m = gsub("mean_NDVI_", "0\\.", model_tidy$term[2]),
      Term = model_tidy$term,
      estimate = round(model_tidy$estimate, 2),
      Std_error = round(model_tidy$std.error, 2),
      Statistic = round(model_tidy$statistic, 3),
      P_value = round(model_tidy$p.value, 3),
      stringsAsFactors = F
    )
  )
}))

# Export model parameters to table
write.csv(moss_model_results, file = "tables/Table S3 interaction models.csv", row.names = F)



# Visualising the moss interaction
# The interaction effect is for two continuous variables (NDVI and moss prop), but for the sake of visualisation, ggpredict() takes the second continuous variable and  splits it into three levels of moss cover
moss_levels <- "moss_prop[0.25, 0.50, 0.90]"  # set levels
legend_loc <- c(0.2, 0.8)

preds_biomass_121 <- ggpredict(mod_biomass_121, terms = c("mean_NDVI_121", moss_levels))
preds_biomass_119 <- ggpredict(mod_biomass_119, terms = c("mean_NDVI_119", moss_levels))
preds_biomass_047 <- ggpredict(mod_biomass_047, terms = c("mean_NDVI_047", moss_levels))
preds_biomass_018 <- ggpredict(mod_biomass_018, terms = c("mean_NDVI_018", moss_levels))

preds_phytomass_121 <- ggpredict(mod_phytomass_121, terms = c("mean_NDVI_121", moss_levels))
preds_phytomass_119 <- ggpredict(mod_phytomass_119, terms = c("mean_NDVI_119", moss_levels))
preds_phytomass_047 <- ggpredict(mod_phytomass_047, terms = c("mean_NDVI_047", moss_levels))
preds_phytomass_018 <- ggpredict(mod_phytomass_018, terms = c("mean_NDVI_018", moss_levels))

preds_leafmass_121 <- ggpredict(mod_leafmass_121, terms = c("mean_NDVI_121", moss_levels))
preds_leafmass_119 <- ggpredict(mod_leafmass_119, terms = c("mean_NDVI_119", moss_levels))
preds_leafmass_047 <- ggpredict(mod_leafmass_047, terms = c("mean_NDVI_047", moss_levels))
preds_leafmass_018 <- ggpredict(mod_leafmass_018, terms = c("mean_NDVI_018", moss_levels))

# Visualisations
# Biomass
(
  plot_biomass_121 <- ggplot() +
    geom_line(
      data = preds_biomass_121,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_biomass_121,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.121 m)",
      y = expression("Biomass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_biomass_119 <- ggplot() +
    geom_line(
      data = preds_biomass_119,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_biomass_119,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.119 m)",
      y = expression("Biomass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_biomass_047 <- ggplot() +
    geom_line(
      data = preds_biomass_047,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_biomass_047,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.047 m)",
      y = expression("Biomass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_biomass_018 <- ggplot() +
    geom_line(
      data = preds_biomass_018,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_biomass_018,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.018 m)",
      y = expression("Biomass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)


# Phytomass
(
  plot_phytomass_121 <- ggplot() +
    geom_line(
      data = preds_phytomass_121,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_phytomass_121,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.121 m)",
      y = expression("Phytomass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_phytomass_119 <- ggplot() +
    geom_line(
      data = preds_phytomass_119,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_phytomass_119,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.119 m)",
      y = expression("Phytomass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_phytomass_047 <- ggplot() +
    geom_line(
      data = preds_phytomass_047,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_phytomass_047,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.047 m)",
      y = expression("Phytomass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_phytomass_018 <- ggplot() +
    geom_line(
      data = preds_phytomass_018,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_phytomass_018,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.018 m)",
      y = expression("Phytomass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)


# Leaf biomass
(
  plot_leafmass_121 <- ggplot() +
    geom_line(
      data = preds_leafmass_121,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_leafmass_121,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.121 m)",
      y = expression("Leaf mass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_leafmass_119 <- ggplot() +
    geom_line(
      data = preds_leafmass_119,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_leafmass_119,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.119 m)",
      y = expression("Leaf mass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_leafmass_047 <- ggplot() +
    geom_line(
      data = preds_leafmass_047,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_leafmass_047,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.047 m)",
      y = expression("Leaf mass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

(
  plot_leafmass_018 <- ggplot() +
    geom_line(
      data = preds_leafmass_018,
      aes(x = x, y = predicted, colour = group),
      size = 1
    ) +
    geom_ribbon(
      data = preds_leafmass_018,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.2
    ) +
    theme_fancy() +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6, face = "italic"),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 4,
        linetype = "blank"
      ),
      legend.position = legend_loc
    ) +
    labs(
      x = "Mean NDVI (0.018 m)",
      y = expression("Leaf mass (g m" ^ "-2" * ")"),
      fill = "Moss cover",
      colour = "Moss cover"
    ) +
    scale_colour_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    ) +
    scale_fill_viridis_d(
      option = "magma",
      direction = -1,
      end = 0.8
    )
)

# Combine plots
all_interactions <-
  ggpubr::ggarrange(
    plot_biomass_121,
    plot_biomass_119,
    plot_biomass_047,
    plot_biomass_018,
    plot_phytomass_121,
    plot_phytomass_119,
    plot_phytomass_047,
    plot_phytomass_018,
    plot_leafmass_121,
    plot_leafmass_119,
    plot_leafmass_047,
    plot_leafmass_018,
    heights = c(12),
    labels = c(
      "(a)",
      "(b)",
      "(c)",
      "(d)",
      "(e)",
      "(f)",
      "(g)",
      "(h)",
      "(i)",
      "(j)",
      "(k)",
      "(l)"
    ),
    font.label = list(size = 10, face = "bold"),
    ncol = 4,
    nrow = 3,
    align = "h"
  )


# Export figures
ggsave(
  all_interactions,
  filename = "plots/Figure S6 - moss interactions.pdf",
  width = 30,
  height = 20,
  units = "cm"
)

ggsave(
  all_interactions,
  filename = "plots/Figure S6 - moss interactions.png",
  width = 30,
  height = 20,
  units = "cm"
)


# Combine plots with patchwork
(
  Figure_5 <-
    leafmass_biomass + phytomass_biomass + plot_phytomass_121 +
    plot_annotation(
      tag_levels = 'a',
      tag_prefix = '(',
      tag_suffix = ')'
    )
)

# Export figure
ggsave(
  Figure_5,
  filename = "plots/Figure 5 - leaf and phytomass vesus biomass and interaction.pdf",
  width = 25,
  height = 8,
  units = "cm"
)

ggsave(
  Figure_5,
  filename = "plots/Figure 5 - leaf and phytomass vesus biomass and interaction.png",
  width = 25,
  height = 8,
  units = "cm"
)




# Figure S1. Boxplot of canopy height observations ----
# Create boxplot
(
  HAG_boxplot <- ggplot(data = dataset_long,
                        aes(
                          x = reorder(PlotID, median, FUN = "median"),
                          y = median,
                          fill = method
                        )) +
    geom_boxplot(
      aes(
        fill = method,
        ymin = min,
        lower = lower,
        middle = median,
        upper = upper,
        ymax = max
      ),
      stat = "identity",
      colour = "grey60",
      outlier.colour = "white",
      outlier.size = 0,
      position = position_dodge(0.7),
      width = 0.7
    ) +
    geom_boxplot(
      aes(
        fill = method,
        upper = upper,
        lower = lower,
        middle = median,
        ymin = ..lower..,
        ymax = ..upper..
      ),
      stat = "identity",
      position = position_dodge(0.7),
      width = 0.7
    ) +  # Added again to include black boarder.
    scale_y_continuous(lim = c(0, 1)) +
    scale_color_manual(values = c("white", "darkgrey")) +
    scale_fill_manual(values = c("white", "darkgrey")) +
    labs(x = "Plot ID", y = "Canopy height (m)") +
    theme_fancy() +
    theme(
      legend.position = c(0.15, 0.9),
      axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 0.5
      )
    )
)

# Export boxplot
ggsave(
  HAG_boxplot,
  filename = "plots/Figure S1 - Height boxplot.pdf",
  width = 16,
  height = 11,
  units = "cm"
)

ggsave(
  HAG_boxplot,
  filename = "plots/Figure S1 - Height boxplot.png",
  width = 16,
  height = 11,
  units = "cm"
)


# Figure S4. Comparison of mean NDVIs  ----
# Boxplot
(
  NDVI_boxplot <- ggplot(data = NDVI_data_long,
                         aes(
                           x = as.factor(grain),
                           y = mean_NDVI,
                           fill = grain
                         )) +
    geom_boxplot() +
    scale_fill_grey(
      start = 0.3,
      end = 0.9,
      aesthetics = "fill"
    ) +
    coord_cartesian(ylim = c(0.6, 0.9)) +
    labs(
      x = "Spatial grain",
      y = expression("Distribution of plot-level mean NDVI"),
      title = "Comparison of NDVI by spatial grain"
    ) +
    theme_fancy() +
    theme(legend.position = "none")
)

# Create barplot
(
  NDVI_barplot <- ggplot(data = NDVI_data_long,
                         aes(
                           x = reorder(PlotID, mean_NDVI, FUN = "median"),
                           y = mean_NDVI,
                           fill = grain
                         )) +
    geom_col(
      aes(y = mean_NDVI, fill = grain),
      width = 0.75,
      position = position_dodge2(preserve = "single")
    ) +
    scale_fill_grey(
      start = 0.3,
      end = 0.8,
      aesthetics = "fill",
      guide = guide_legend(direction = "horizontal")
    ) +
    coord_cartesian(ylim = c(0.5, 0.9), expand = FALSE) +
    labs(
      x = "Plot ID",
      y = expression("mean NDVI"),
      title = "Comparison of mean NDVI by plot"
    ) +
    theme_fancy() +
    theme(
      legend.position = c(0.35, 0.96),
      axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 0.5
      )
    )
)

# Combine plots
NDVI_plots <- ggpubr::ggarrange(
  NDVI_boxplot,
  NDVI_barplot,
  heights = c(10, 15),
  labels = c("(a)", "(b)"),
  font.label = list(size = 10, face = "bold"),
  ncol = 1,
  nrow = 2,
  align = "h"
)

# Export barplot
ggsave(
  NDVI_plots,
  filename = "plots/Figure S4 - NDVI comparison.pdf",
  width = 16,
  height = 18,
  units = "cm"
)

ggsave(
  NDVI_plots,
  filename = "plots/Figure S4 - NDVI comparison.png",
  width = 16,
  height = 18,
  units = "cm"
)



# SI Figure of plant functional groups across the 36 plots ----
# Calculate percentage cover of different veg covers per plot
PF_observations3 <- PF_observations %>%
  drop_na(Count) %>%
  group_by(PlotN) %>%
  mutate(total_hits = sum(Count)) %>%
  ungroup() %>%
  group_by(PlotN, Species) %>%
  mutate(species_hits = sum(Count)) %>%
  ungroup() %>%
  mutate(percentage_cover = species_hits/total_hits)

PF_observations3$Species <- factor(PF_observations3$Species,
                                   levels = c("Dryas integrifolia", "Equisetum", "Salix arctica", "Salix richardsonii",
                                              "XXXbareground",  "XXXfungus", "XXXlitter", "XXXotherforb",
                                              "XXXothergram", "XXXothermoss", "XXXvegetatedground"),
                                   labels = c("Dryas integrifolia", "Equisetum spp.",
                                   "Salix arctica", "Salix richardsonii", "Bare ground",
                                   "Fungus", "Leaf litter", "Forbs", "Graminoids",
                                   "Moss", "Vegetated ground"))

(plot_cover_fig <- ggplot(PF_observations3, aes(x = PlotN, y = percentage_cover, 
                             colour = Species, fill = Species)) +
    geom_bar(stat = "identity", position = "fill",
             width = 0.5) +
    coord_flip() +
    labs(x = "Plot ID\n", y = "\nProportion of contacts") +
    theme_fancy() +
    scale_fill_manual(values = c("#fcba03", "#9de650",
                                 "#649c28", "#025c03", "#545954",
                                 "#aa84ab", "#946d41", "#cda3d9", "#73b4ba",
                                 "#b3ba73", "#346078")) +
    scale_colour_manual(values = c("#fcba03", "#9de650",
                                   "#649c28", "#025c03", "#545954",
                                   "#aa84ab", "#946d41", "#cda3d9", "#73b4ba",
                                   "#b3ba73", "#346078")) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 8, face = "italic"),
          axis.line.x = element_line(),
          axis.line.y = element_line()))



ggsave(
  plot_cover_fig,
  filename = "plots/Figure S2 - taxa proportion.png",
  width = 16,
  height = 18,
  units = "cm"
)

ggsave(
  plot_cover_fig,
  filename = "plots/Figure S2 - taxa proportion.pdf",
  width = 16,
  height = 18,
  units = "cm"
)




### SI figure of biomass components ----
# prepare data
df_biomass <- data.frame("plotID" = rep(dataset$PlotID,3),
                         "biomass" = c((dataset$Comp1_Mass*4), (dataset$Comp2_Mass*4), (dataset$Comp3_Mass*4)),
                         "partition" = c(rep("woody",36),rep("leaves",36),rep("herbaceous",36))
                         )

# create plot
(plot_biomass_components <- ggplot(df_biomass, aes(x = plotID, y = biomass, 
                                                colour = partition, fill = partition)) +
   geom_bar(stat = "identity",
            width = 0.5) +
   coord_flip() +
   labs(x = "Plot ID\n", 
        y = expression("Dry biomass (g m" ^ "-2" * ")")) +
   theme_fancy() +
   scale_fill_manual(values = c("#aaeb07", "#2c8232","#755504")) +
   scale_colour_manual(values = c("#aaeb07", "#2c8232","#755504")) +
   theme(legend.position = "bottom",
         legend.title = element_blank(),
         legend.text = element_text(size = 8),
         axis.line.x = element_line(),
         axis.line.y = element_line())
  )

# save plot
ggsave(
  plot_biomass_components,
  filename = "plots/Figure S3 - biomass components.png",
  width = 16,
  height = 18,
  units = "cm"
)

ggsave(
  plot_biomass_components,
  filename = "plots/Figure S3 - biomass components.pdf",
  width = 16,
  height = 18,
  units = "cm"
)




### SI Figure of terrain model error ----
# prepare data
DTM_accuracy$z_error <- DTM_accuracy$DTM_elevation - DTM_accuracy$GNSS_elevation

DTM_error_mean <- round(mean(DTM_accuracy$z_error), 3)
DTM_error_SD <- round(sd(DTM_accuracy$z_error), 3)

# create plot
(plot_terrain_error <- ggplot(DTM_accuracy, aes(x = z_error)) +
    geom_density(fill="lightgrey") +
    labs(x = "\nZ error (m)\n (DTM elevation - GNSS elevation)",
         y = "Density\n") +
    theme_fancy() +
    # scale_x_continuous(breaks=c(-0.1, 0.1, 0.01)) +
    scale_x_continuous(lim = c(-0.15, 0.15)) +
    
    
    geom_vline(aes(xintercept=mean(z_error)),
               color="black", linetype="dashed", size=1)
)

# save plot
ggsave(
  plot_terrain_error,
  filename = "plots/Figure S7 - DTM error.png",
  width = 16,
  height = 18,
  units = "cm"
)

ggsave(
  plot_terrain_error,
  filename = "plots/Figure S7 - DTM error.pdf",
  width = 16,
  height = 18,
  units = "cm"
)







#### Upscaling for landscape biomass estimate #### 

# Load rasters
rast_CHM <- raster("data/site_rasters/CHM.tif")                                      # Read in canopy height model at 0.010 m grain.
rast_NDVI_018 <- raster("data/site_rasters/NDVI_018.tif")                            # Read in NDVI raster at 0.018 m grain.
rast_NDVI_047 <- raster("data/site_rasters/NDVI_047.tif")                            # Read in NDVI raster at 0.047 m grain.
rast_NDVI_119 <- raster("data/site_rasters/NDVI_119.tif")                            # Read in NDVI raster at 0.119 m grain.
rast_NDVI_121 <- raster("data/site_rasters/NDVI_121.tif")                            # Read in NDVI raster at 0.121 m grain.
rast_RGB <- stack("data/site_rasters/RGB_040.tif")                                   # Read in multi-band RGB orthmosaic at 0.040 m grain.


# Load polygon of monitoring area (corner marker locations)
AOI <- st_read("data/site_rasters/monitoring_plot.geojson", crs = 32607)        # Import polygon of monitoirng plot as sf object, using st_read to allow the non-standard CRS to be specified.
st_area(AOI)                                                                    # Return area of the AOI (with units)
AOI_area <- as.numeric(st_area(AOI))                                                # Return area of AOI (without unit attributes)


# Clip rasters to the monitoring extent
rast_AOI_CHM <- crop(rast_CHM, AOI)
rast_AOI_NDVI_018 <- crop(rast_NDVI_018, AOI)
rast_AOI_NDVI_047 <- crop(rast_NDVI_047, AOI)
rast_AOI_NDVI_119 <- crop(rast_NDVI_119, AOI)
rast_AOI_NDVI_121 <- crop(rast_NDVI_121, AOI)
rast_AOI_RGB <- crop(rast_RGB, AOI)

# Load rasters
# rast_AOI_CHM <- raster("data/site_rasters/rast_AOI_CHM.tif")
# rast_AOI_RGB <- brick("data/site_rasters/rast_AOI_RGB.tif")
# rast_AOI_NDVI_018 <- raster("data/site_rasters/rast_AOI_NDVI_018.tif")
# rast_AOI_NDVI_047 <- raster("data/site_rasters/rast_AOI_NDVI_047.tif")
# rast_AOI_NDVI_119 <- raster("data/site_rasters/rast_AOI_NDVI_119.tif")
# rast_AOI_NDVI_121<- raster("data/site_rasters/rast_AOI_NDVI_121.tif")

# Review value distributions
hist(rast_AOI_CHM,
     main="Distribution of canopy height Values",
     xlab="Canopy height (m)",
     ylab="Frequency",
     col="grey")

hist(rast_AOI_NDVI_018,
     main="Distribution of NDVI Values (0.018m)",
     xlab="NDVI",
     ylab="Frequency",
     col="grey",
     xlim = c(0,1))

hist(rast_AOI_NDVI_047,
     main="Distribution of NDVI Values (0.047m)",
     xlab="NDVI",
     ylab="Frequency",
     col="grey",
     xlim = c(0,1))

hist(rast_AOI_NDVI_119,
     main="Distribution of NDVI Values (0.119m)",
     xlab="NDVI",
     ylab="Frequency",
     col="grey",
     xlim = c(0,1))

hist(rast_AOI_NDVI_121,
     main="Distribution of NDVI Values (0.121m)",
     xlab="NDVI",
     ylab="Frequency",
     col="grey",
     xlim = c(0,1))


# compute mean canopy height
meanCH <- exact_extract(rast_AOI_CHM, AOI, 'mean')                              # USed exact_extract because raster::extract is much slower.


# calculate biomass maps (units in g m^2 per pixel, to facilitate comparison between different spatial grain protucts)
rast_Biomass_CHM <- rast_AOI_CHM * model_SfM$coefficients[1] 
rast_Biomass_NDVI_018 <- coef(exp_model_total_NDVI_018)[1] * exp(coef(exp_model_total_NDVI_018)[2] * rast_AOI_NDVI_018)
rast_Biomass_NDVI_047 <- coef(exp_model_total_NDVI_047)[1] * exp(coef(exp_model_total_NDVI_047)[2] * rast_AOI_NDVI_047)
rast_Biomass_NDVI_119 <- coef(exp_model_total_NDVI_119)[1] * exp(coef(exp_model_total_NDVI_119)[2] * rast_AOI_NDVI_119)
rast_Biomass_NDVI_121 <- coef(exp_model_total_NDVI_121)[1] * exp(coef(exp_model_total_NDVI_121)[2] * rast_AOI_NDVI_121)


# Calculate total biomass in each raster, and convert to standard units (Mg ha-1)
Biomass_CHM_gm2 <- 
  round(cellStats(rast_Biomass_CHM, 'sum') / ncell(rast_Biomass_CHM), 1)                  # mean biomass in g m-2
Biomass_CHM_Mgha <- 
  round(cellStats(rast_Biomass_CHM, 'sum') / ncell(rast_Biomass_CHM) / 100, 2)                                                       # mean biomass in Mg ha-1

Biomass_NDVI_018_gm2 <- 
  round(cellStats(rast_Biomass_NDVI_018, 'sum') / ncell(rast_Biomass_NDVI_018), 1)        # mean biomass in g m-2
Biomass_NDVI_018_Mgha <- 
  round(cellStats(rast_Biomass_NDVI_018, 'sum') / ncell(rast_Biomass_NDVI_018) / 100, 2)                                                   # mean biomass in Mg ha-1

Biomass_NDVI_047_gm2 <- 
  round(cellStats(rast_Biomass_NDVI_047, 'sum') / ncell(rast_Biomass_NDVI_047), 1)        # mean biomass in g m-2
Biomass_NDVI_047_Mgha <- 
  round(cellStats(rast_Biomass_NDVI_047, 'sum') / ncell(rast_Biomass_NDVI_047) / 100, 2)                                                  # mean biomass in Mg ha-1

Biomass_NDVI_119_gm2 <- 
  round(cellStats(rast_Biomass_NDVI_119, 'sum') / ncell(rast_Biomass_NDVI_119), 1)        # mean biomass in g m-2
Biomass_NDVI_119_Mgha <- 
  round(cellStats(rast_Biomass_NDVI_119, 'sum') / ncell(rast_Biomass_NDVI_119) / 100, 2)                                                  # mean biomass in Mg ha-1

Biomass_NDVI_121_gm2 <- 
  round(cellStats(rast_Biomass_NDVI_121, 'sum') / ncell(rast_Biomass_NDVI_121), 1)        # mean biomass in g m-2
Biomass_NDVI_121_Mgha <- 
  round(cellStats(rast_Biomass_NDVI_121, 'sum') / ncell(rast_Biomass_NDVI_121) / 100, 2)                                                  # mean biomass in Mg ha-1


# create dataframe of total biomass estimates
df_biomass_est <- data.frame("Raster" = c("CHM",
                                          "NDVI 0.018 m",
                                          "NDVI 0.047 m",
                                          "NDVI 0.119 m",
                                          "NDVI 0.121 m"),
                             "Biomass_Mg_ha1" = c(Biomass_CHM_Mgha,
                                                  Biomass_NDVI_018_Mgha,
                                                  Biomass_NDVI_047_Mgha,
                                                  Biomass_NDVI_119_Mgha,
                                                  Biomass_NDVI_121_Mgha),
                             "Biomass_g_m2" = c(Biomass_CHM_gm2,
                                                Biomass_NDVI_018_gm2,
                                                Biomass_NDVI_047_gm2,
                                                Biomass_NDVI_119_gm2, 
                                                Biomass_NDVI_121_gm2)
                             )

# Export total biomass estimates
write.csv(df_biomass_est,"tables/Table 3 Biomass estimates.csv", row.names = FALSE)            # extracted NDVI values were added to the main_database file. ndvi_data <- read.csv("data/Extracted_NDVI.csv", header = T)                  # Read in NDVI values from Exact Extract pipeline.

## write clipped rasters to faciliate sharing
# writeRaster(rast_AOI_CHM, "data/site_rasters/rast_AOI_CHM.tif")
# writeRaster(rast_AOI_RGB, "data/site_rasters/rast_AOI_RGB.tif")
# writeRaster(rast_AOI_NDVI_018, "data/site_rasters/rast_AOI_NDVI_018.tif")
# writeRaster(rast_AOI_NDVI_047, "data/site_rasters/rast_AOI_NDVI_047.tif")
# writeRaster(rast_AOI_NDVI_119, "data/site_rasters/rast_AOI_NDVI_119.tif")
# writeRaster(rast_AOI_NDVI_121, "data/site_rasters/rast_AOI_NDVI_121.tif")


### Calculate biomass difference rasters
# Calcualting difference maps requires rasters with identical resolution and 
# alignment. Our aim here is illustrate the spatial distribution of differences
# in estimated biomass (which are on the order of 25%-70%, so are pretty large 
# relative to errors expected from resampling. Bi-linear resampling introduces 
# less uncertainty if we also aggregate a bit and we wanna make sure that the 
# new resolution is almost a multiple of the smaller resolutions.
# Checking out the relative rest of the division (using modulo)
(200%%c(18,47,119,121)) / c(18,47,119,121)  # 0.2 m
(250%%c(18,47,119,121)) / c(18,47,119,121)  # 0.25 m
# 0.25 m would end up with less bilinear resampling for the coarser resolution,
# the relative remainder for the smaller resolution is larger, but we also 
# aggregate across more pixels of those, so the error will be smaller, so we're
# going for 0.25 m. 


# Create vector with raster names for convenience in handling later
raster_names <- c("rast_Biomass_CHM", 
                  "rast_Biomass_NDVI_018",
                  "rast_Biomass_NDVI_047",
                  "rast_Biomass_NDVI_119",
                  "rast_Biomass_NDVI_121")

# check the extent of all rasters matches (just as I don't know them)
lapply(raster_names, function(x) extent(get(x)))
# They match optically, but R rounds up the display with low precision. 
# So the actual floats don't  match:
lapply(raster_names, function(x) extent(get(x))@xmin == 581926.7) 
# and: 
lapply(raster_names, function(x) sprintf("%.10f",extent(get(x))@xmin)) 
lapply(raster_names, function(x) sprintf("%.10f",extent(get(x))@xmax))
lapply(raster_names, function(x) sprintf("%.10f",extent(get(x))@ymin))
lapply(raster_names, function(x) sprintf("%.10f",extent(get(x))@ymax))

# We solved this by trimming the edges, this introduces more error in the
# re-sampling, but can't easily be avoided. (looking at those differences, we
# could also just re-sample to 0.2 m, but let's stick to the 0.25 m planned).
# Ceate an empty target raster with the desired properties and then resample to that. 
target_raster <- raster(xmn = 581926.75,
                        xmx = 582040.50,
                        ymn = 7719550.50,
                        ymx = 7719597.00,
                        crs = crs(rast_AOI_CHM),
                        res = 0.25)

list2env(
  lapply(setNames(raster_names,
                paste0(raster_names, "_coarse")),
       function(x) raster::resample(get(x), target_raster)),
       envir = .GlobalEnv)

# Compute difference maps
list2env(
  lapply(setNames(paste0(raster_names, "_coarse")[-1],
                  paste0(raster_names, "_coarse_diff")[-1]),
         function(x) rast_Biomass_CHM_coarse - get(x)),
  envir = .GlobalEnv)

# Check them quickly:
# library(rasterVis)
# lapply(paste0(raster_names, "_coarse_diff")[-1],
#          function(x) levelplot(get(x)))
# levelplot(rast_Biomass_CHM_coarse)
# levelplot(rast_Biomass_NDVI_018_coarse)
# levelplot(rast_Biomass_NDVI_018_coarse_diff)

png("plots/Figure S8.png")
par(mfrow = c(2, 2)) 

hist(rast_Biomass_NDVI_018_coarse_diff,
     main="Biomass residuals (NDVI 0.018m)",
     xlab="Diff. Biomass CHM - NDVI (g m2)",
     ylab="Frequency",
     col="grey",
     xlim = c(-3000,2000))

hist(rast_Biomass_NDVI_047_coarse_diff,
     main="Biomass residuals (NDVI 0.047m)",
     xlab="Diff. Biomass CHM - NDVI (g m2)",
     ylab="Frequency",
     col="grey",
     xlim = c(-3000,2000))

hist(rast_Biomass_NDVI_119_coarse_diff,
     main="Biomass residuals (NDVI 0.119m)",
     xlab="Diff. Biomass CHM - NDVI (g m2)",
     ylab="Frequency",
     col="grey",
     xlim = c(-3000,2000))

hist(rast_Biomass_NDVI_121_coarse_diff,
     main="Biomass residuals (NDVI 0.121m)",
     xlab="Diff. Biomass CHM - NDVI (g m2)",
     ylab="Frequency",
     col="grey",
     xlim = c(-3000,2000))
dev.off()
par(mfrow = c(1, 1)) 


### Final figure for publication with  15 rasters together (3 cols x 5 rows)
# Include the RGB image in the empty space left because the CHM is not differenced against itself.
# ggplot can handle rasters well, but it's ability to do so well is still developing
# Most people use the rasterVis package which uses lattice plots that can be arranged into multiple grobs with gridExtra
# This can all be nicely done with levelplot from the rasterVis package

# Create dataframe of rasters to plot
rasters_to_plot <- data.frame(
  raster_name = c("rast_AOI_CHM",
                  "rast_Biomass_CHM",
                  "rast_AOI_RGB",
                  "rast_AOI_NDVI_018",
                  "rast_Biomass_NDVI_018",
                  "rast_Biomass_NDVI_018_coarse_diff",
                  "rast_AOI_NDVI_047",
                  "rast_Biomass_NDVI_047",
                  "rast_Biomass_NDVI_047_coarse_diff",
                  "rast_AOI_NDVI_119",
                  "rast_Biomass_NDVI_119",
                  "rast_Biomass_NDVI_119_coarse_diff",
                  "rast_AOI_NDVI_121",
                  "rast_Biomass_NDVI_121",
                  "rast_Biomass_NDVI_121_coarse_diff"),
  plot_type = c("chm",
                 "biomass",
                 "rgb",
                 rep(c("ndvi",
                       "biomass",
                       "diff"), 4)),
  col_ramp_title = c("'Canopy Height (m)'",
                    "'Biomass (g m' ^ -2 * ')'",
                    "'none'",
                    rep(c("'NDVI'",
                          "'Biomass (g m' ^ -2 * ')'",
                          "'Δ Biomass (g m' ^ -2 * ')'"), 4)),
  col_ramp_min = c(0,0,NA, 
                   rep(c(0,0,-3000), 4)), # Minimum scale values
  col_ramp_max = c(1.4,3500,NA,
                   rep(c(1,3500,2000), 4)), # Maximum scale values
  col_ramp_step = c(0.2,500,NA,
                    rep(c(0.2,500,1000), 4)), # Color ramp klegend steps
  panel_label = paste0("(",
                       letters[1:15],
                       ") ",
                       c("SfM Canopy Height",
                         "Biomass SfM",
                         "RGB",
                         "NDVI (0.018 m)",
                         "Biomass NDVI (0.018 m)",
                         "Diff. Biomass CHM - NDVI (0.018 m)",
                         "NDVI (0.047 m)",
                         "Biomass NDVI (0.047 m)",
                         "Diff. Biomass CHM - NDVI (0.047 m)",
                         "NDVI (0.119 m)",
                         "Biomass NDVI (0.119 m)",
                         "Diff. Biomass CHM - NDVI (0.119 m)",
                         "NDVI (0.121 m)",
                         "Biomass NDVI (0.121 m)",
                         "Diff. Biomass CHM - NDVI (0.121 m)")),
  panel_label_xpos = rep(c(0.07,0.07, 0.06),5),
  stringsAsFactors = F)

# Generate colour ramps with the colorspace package
chm_col <- sequential_hcl(100, palette = "Blues")
ndvi_col <- sequential_hcl(100, palette = "Oranges")
biomass_col <- sequential_hcl(100, palette = "Greens")
diff_col <- diverging_hcl(600, palette = "Purple-Brown")[1:500] # Calculate 600 values, only use 500 due to imbalance aorund 0 


plot_pretty_raster <- function(raster_name) {
  # Get raster object
  raster_to_plot <- get(raster_name)
  # get plot type
  plot_type <- rasters_to_plot$plot_type[rasters_to_plot$raster_name == raster_name]
  
   # If RGB raster do the following
  if(plot_type == "rgb") {
    raster_to_plot <- raster_to_plot[[1:3]]  # throw out alpha band
    # It is huge, so when trying layout changes aggregate to a sensible res first
    # raster_to_plot <- raster::aggregate(raster_to_plot, 5)
    # To use lattice (for compatibility with the other raster plots) we
    # need to create the RGB colourspace ourselves.
    
    # Create RGB color for cell values
    cols <- factor(rgb(raster_to_plot[], maxColorValue=255))
    
    # Creat single band temp raster
    temp_raster <- raster(raster_to_plot)
    # re-assign cell values
    temp_raster[] <- cols
    
    # Plot
    pretty_plot <- levelplot(temp_raster, 
                             margin=FALSE,  # don't plot margins
                             scales=list(draw=FALSE), # suppress axis labels
                             col.regions=as.character(levels(cols)),
                             colorkey=FALSE,
                             xlab.top = list("ffadsf", col = "white", cex = 0.8), # Quick work around to make sure spacing is the same
                             #xlab = list("dfsafasdf", col = "white", cex = 1), # Quick work around to make sure spacing is the same
                             ylab.right = list("sdfsfa", col = "white", cex = 6.75) # This is a quick work around to replace the color key bar...
                             ) + 
      # Add scale bar
      latticeExtra::layer({
        ## Scale bar
        # Determine position of scale bar (bottom left corner)
        scale_bar_length <- 20 # Scale bar length 20 m
        scale_bar_height <- 2.5 # scale bar height 2.5 m (seeme like a good start for 50 m height raster)
        scale_bar_nsegments <- 2 # 4 segments
        scale_bar_xpos <- 10 # x position from bottom left corner in percent
        scale_bar_ypos <- 10 # y position from bottom left corner in percent
        scale_bar_col <- rasters_to_plot$scale_bar_col[rasters_to_plot$raster_name == raster_name]
        
        # calculate derived parameters:
        # (there is a bug in lattice so we need to shove it into the global enviornment)
        scale_bar_xmin <- extent(raster_to_plot[[1]])@xmin + 
          ((extent(raster_to_plot[[1]])@xmax -
              extent(raster_to_plot[[1]])@xmin) / 100 * scale_bar_xpos)
        scale_bar_xmax <- scale_bar_xmin + scale_bar_length 
        scale_bar_ymin <- extent(raster_to_plot[[1]])@ymin + 
          ((extent(raster_to_plot[[1]])@ymax -
              extent(raster_to_plot[[1]])@ymin) / 100 * scale_bar_ypos)
        scale_bar_step <- scale_bar_length / scale_bar_nsegments
        
        xs <- seq(scale_bar_xmin, scale_bar_xmax, scale_bar_step)
        grid.rect(x = xs[1:(length(xs)-1)],
                  y = scale_bar_ymin,
                  width = scale_bar_step, height=scale_bar_height,
                  gp= gpar(fill = rep(c('transparent', "white"),
                                      2),
                           col = "white"),
                  default.units='native')
        grid.text(x = xs - (scale_bar_step / 2), 
                  y = scale_bar_ymin + scale_bar_height * 1.5,
                  paste(seq(0, scale_bar_length, scale_bar_step), "m"),
                  gp=gpar(cex=0.8, col = "white"),
                  default.units='native')
      }, data = list(raster_name = raster_name,
                     raster_to_plot = raster_to_plot)) # + 
    #### Add north arrow (can be replaced with any polygon, this one is ugly)
    # latticeExtra::layer({
    #   SpatialPolygonsRescale(layout.north.arrow(),
    #                          offset = c(scale_bar_xmin + scale_bar_length * 1.1,
    #                                     scale_bar_ymin),
    #                          scale = scale_bar_height * 2)
    #   theme = list(col.regions = "white")
    # }) 
  } else { # else....
    # Set color ramp
    col_ramp <- get(paste0(plot_type, "_col"))
    # Set color key label
    col_key_label <- rasters_to_plot$col_ramp_title[rasters_to_plot$raster_name == raster_name]
    
    # Set col ramp limits 
    col_scale_min <- rasters_to_plot$col_ramp_min[rasters_to_plot$raster_name == raster_name]
    col_scale_max <- rasters_to_plot$col_ramp_max[rasters_to_plot$raster_name == raster_name]
    col_scale_step <- rasters_to_plot$col_ramp_step[rasters_to_plot$raster_name == raster_name]
    
    pretty_plot <- levelplot(raster_to_plot, 
                             margin=FALSE,                       # don't plot margins
                             colorkey=list(
                               space='right',                   # right
                               labels=list(at=seq(col_scale_min,
                                                  col_scale_max,
                                                  col_scale_step), 
                                           font=1, cex = 1)  # Colorkey axis text
                             ),    
                             legend=list(
                               top=list(
                                 fun=grid::textGrob(
                                   parse(text = col_key_label),
                                   y=1.1, 
                                   x=1.07),
                                 cex = 1)),
                             scales=list(draw=FALSE),            # suppress axis labels
                             col.regions=col_ramp,                   # colour ramp
                             at=seq(col_scale_min, col_scale_max, length.out = 100))
  }
  # Add label to top left
  panel_label <- rasters_to_plot$panel_label[rasters_to_plot$raster_name == raster_name]
  panel_label_xpos <- rasters_to_plot$panel_label_xpos[rasters_to_plot$raster_name == raster_name]
  pretty_plot_with_label <- arrangeGrob(pretty_plot,
                         top = textGrob(panel_label,
                                        x = unit(panel_label_xpos, "npc"), 
                                        y   = unit(0, "npc"), 
                                        just=c("left","top"),
                                        gp=gpar(col="black", 
                                                fontsize=18, 
                                                fontfamily="Arial")))
  return(pretty_plot_with_label)
}

# Execute for all plots
pretty_plot_list <- lapply(rasters_to_plot$raster_name,
                       plot_pretty_raster)

# Export PNG for plot grid
png("plots/Figure 6 Biomass Maps.png", 
    width = 7 * 3,
    height = 3 * 5,
    units = "in",
    res = 300)
print(grid.arrange(grobs = pretty_plot_list,
                   ncol = 3))
dev.off()




