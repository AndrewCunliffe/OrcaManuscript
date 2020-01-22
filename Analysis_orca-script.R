###  Script for analysing observatiosn from the Orca project                    #
###  Written by Andrew Cunliffe (andrewmcunliffe@gmail.com) & Gergana Daskalova #
###                                                                             #

# Establish operating environment ----
path <- getwd()
# options(repos = c(CRAN = "http://cran.rstudio.com"))  # set default CRAN mirror

# Load required packages ----
# # Download
# install.packages("exactextractr")

 
# Install libr
library(tidyverse)
library(viridis)
library(grid) # required for plot annotation 
library(gridExtra) # for arranging multi-panel plots
library(xlsx)
library(ggpubr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(DescTools)  # for computing CCC
library(ggpmisc) # for adding model parameters to plots


# Load data ----
    # Primary dataet
dataset <- read.csv("outputs_01_main/processed_master_database_combined.csv", header = T)  # Read in summary  data
ndvi_data <- read.csv("outputs_01_main/Extracted NDVI.csv", header = T)  # Read in NDVI values from Exact Extract pipeline.
    





    # Canopy height from point framing
    PF_observations <- read.csv("outputs_01_main/PointFraming_2016_tidy.csv")
    PF_observations <- PF_observations[,-1]                                     # remove unwanted index position.

# Data preparation/manipulation ----
    # Pointframe data - Extracting only Height and PlotN from the pointframe data, and using a pipe to omit NAs
    PF_observations2 <- dplyr::select(PF_observations, Height, PlotN) %>%
        na.omit(PF_observations)

    # Convert height ovservations from cm into m. Could be integrated into the above pipe.
    PF_observations2$Height <- PF_observations2$Height/100
    
    # Pipe to compute canopy height summary metrics (min, max, median, mean and quantiles) 
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
    # Remoce unwanted/erronious individual height. This could be integrated into the preceeding pipe...
    PF_HAG_summary <- subset(PF_HAG_summary, select = -Height)
    
    # Add point framing HAGs into main dataframe
    dataset$PF_HAG_min <- PF_HAG_summary$min
    dataset$PF_HAG_max <- PF_HAG_summary$max
    dataset$PF_HAG_median <- PF_HAG_summary$median
    dataset$PF_HAG_mean <- PF_HAG_summary$mean
    dataset$PF_HAG_upper_ICR <- PF_HAG_summary$upper_ICR
    dataset$PF_HAG_lower_IQR <- PF_HAG_summary$lower_IQR
    
# Visualisation ----
# Plotting themes ----
theme_coding <- function(){
    theme_bw()+
        theme(axis.text = element_text(size = 6),
              axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
              axis.title = element_text(size = 8),
              panel.grid = element_blank(),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
              plot.title = element_text(size = 10, vjust = 1, hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 5, face = "italic"),
              legend.key.size = unit(0.9,"line"),
              legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"),
              legend.position = c(0.9, 0.9))
}

# Compute scaling parameters from minimum and maximum values ----
max_agb <- 1.1*max(dataset$AGB_spatially_normalised_g_m2, na.rm = TRUE)
max_hag <- 1.1*max(max(dataset$HAG_plotmax_of_cellmax_m, na.rm = TRUE), max(dataset$PF_HAG_max))
max_mean_hag <- 1.1*max(max(dataset$HAG_plotmean_of_cellmax_m, na.rm = TRUE), max(dataset$PF_HAG_mean))
# max_ndvi <- max(max(dataset$mean_NDVI_018), max(dataset$mean_NDVI_047), max(dataset$mean_NDVI_119), max(dataset$mean_NDVI_121))
max_ndvi <- 0.91
# min_ndvi <- min(min(dataset$mean_NDVI_018), min(dataset$mean_NDVI_047), min(dataset$mean_NDVI_119), min(dataset$mean_NDVI_121))
min_ndvi <- 0.6


# SfM versus PF HAG (comparison of height above ground estimates) ----
# (Structure from motion versus point framing)
    # mean HAG scatter plot
        # Create plot
        PF_Vs_SfM_HAG_mean <- ggplot(data = dataset,
                                     aes(x = PF_HAG_mean, y = HAG_plotmean_of_cellmax_m)) +
                                geom_point() +
                                labs(x = "Point Frame - Canopy Height (m)",
                                     y = "SfM - Canopy height (m)",
                                     title = "Mean") +
                                theme_coding() +
                                coord_cartesian(ylim = c(0, 1), xlim = c(0, 1), expand=FALSE) +
                                geom_abline(intercept = 0, slope = 1, linetype = "dashed")
            
        # Export plot
        png(filename = "outputs_01_main/PF_HAG_Vs_SfM_HAG_mean.png", width = 10, height = 10, units = "cm", res = 400)
        plot(PF_Vs_SfM_HAG_mean)
        dev.off()

    # median HAG scatter plot
        # Create plot
        PF_Vs_SfM_HAG_median <- ggplot(data = dataset,
                                       aes(x = PF_HAG_median, y = HAG_plotmedian_of_cellmax_m)) +
            geom_point() +
            labs(x = "Point Frame - Canopy Height (m)",
                 y = "SfM - Canopy height (m)",
                 title = "Median") +
            theme_coding() +
            coord_cartesian(ylim = c(0, 1), xlim = c(0, 1), expand=FALSE) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed")

        # Export plot
        png(filename = "outputs_01_main/PF_Vs_FsM_HAG_median.png", width = 10, height = 10, units = "cm", res = 400)
        plot(PF_Vs_SfM_HAG_median)
        dev.off()

    # Boxplot of canopy height distributions.
        # The intention here is to illustrate the similarity (and bias) between these two methods of observing canopy height.
        # Create combined long form dataset for plotting
            PlotID <- rep(dataset$PlotID, times = 2)
            method <- rep(c('SfM', 'Point Framing'), each = 36)
            min <- c(dataset$HAG_plotmin_of_cellmax_m, dataset$PF_HAG_min)
            lower <- c(dataset$HAG_plot25percentile_of_cellmax_m, dataset$PF_HAG_lower_IQR)
            median <- c(dataset$HAG_plotmedian_of_cellmax_m, dataset$PF_HAG_median)
            upper <- c(dataset$HAG_plot75percentile_of_cellmax_m, dataset$PF_HAG_upper_ICR)
            max <- c(dataset$HAG_plotmax_of_cellmax_m, dataset$PF_HAG_max)
        dataset_long <- data.frame(PlotID, method, min, lower, median, upper, max)  # Create vectors into new dataframe
        
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
                      axis.line.y = element_line(color="black", size = 0.5)))
        
        # Export boxplot
            png(filename = "outputs_01_main/HAG_boxplot.png", width = 16, height = 11, units = "cm", res = 400)
            plot(HAG_boxplot)
            dev.off()
        
        
    # Statistical analysis of HAG agreement.
    # Bias summary
        HAG_rediduals <- dataset$HAG_plotmean_of_cellmax_m - dataset$PF_HAG_mean
        HAG_bias_mean <- mean(HAG_rediduals)
        HAG_bias_mean
        HAG_bias_mean_SD <- sd(HAG_rediduals)
        HAG_bias_mean_SD
        
        HAG_rediduals <- dataset$HAG_plotmedian_of_cellmax_m - dataset$PF_HAG_median
        HAG_bias_median <- mean(HAG_rediduals)
        HAG_bias_median
        HAG_bias_median_SD <- sd(HAG_rediduals)
        HAG_bias_median_SD

        # Inference
        # There is fractionally better agreement between the median rather than the mean observations of canopy height.
        
    # Concordence correlation coefficient
        # I'm still trying to identify which is the most popular package to use...
        # in any case, if we report CCC we should report the package used in the methods.
     # mean CCC
        test <- CCC(dataset$HAG_plotmean_of_cellmax_m,
            dataset$PF_HAG_mean,
            ci = "z-transform",
            conf.level = 0.95,
            na.rm = FALSE)

      test$rho.c
      
      # median CCC
      test <- CCC(dataset$HAG_plotmedian_of_cellmax_m,
                  dataset$PF_HAG_median,
                  ci = "z-transform",
                  conf.level = 0.95,
                  na.rm = FALSE)
      
      test$rho.c
      
      
      names(dataset)
# HAG as predictor of AGB ----
      # median # NOT USED!
      (biomass_SfM <- ggplot(data = dataset,
                            aes(x = HAG_plotmedian_of_cellmax_m,
                                y = AGB_spatially_normalised_g_m2)) +
                            geom_point(na.rm = TRUE) +
                            theme_coding() +
                            coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 1), expand=FALSE) +
                            labs(x = expression("Canopy height (m)"),
                                 y = expression("Dry biomass (g m"^"-2"*")"),
                                 title = "SfM") +
                            theme(legend.position = c(0.15, 0.9)) +
                            geom_smooth(method="lm", formula= y ~ x-1, se=TRUE, size=0.5, na.rm = TRUE))

      (biomass_PF <- ggplot(data = dataset,
                            aes(x = PF_HAG_median,
                                y = AGB_spatially_normalised_g_m2)) +
              geom_point(na.rm = TRUE) +
              theme_coding() +
              coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 1), expand=FALSE) +
              labs(x = expression("Canopy height (m)"),
                   y = expression("Dry biomass (g m"^"-2"*")"),
                   title = "Point Framing") +
              theme(legend.position = c(0.15, 0.9)) +
              geom_smooth(method="lm", formula= y ~ x-1, se=TRUE, size=0.5, na.rm = TRUE))

      # Combine plots with ggarrange()[in ggpubr]
      biomass_plots_median <- ggarrange(biomass_PF, biomass_SfM,
                                    heights = c(10,10),
                                    labels = c("(a)", "(b)"),
                                    ncol = 2, nrow = 1,
                                    align = "h")

      # Export figure
      png(filename="outputs_01_main/HAG Vs biomass_median.png", width=16, height=9, units="cm", res=400)
      plot(biomass_plots_median)
      dev.off()
      
      # mean
      # (biomass_SfM_mean <- ggplot(data = dataset,
      #                        aes(x = HAG_plotmean_of_cellmax_m,
      #                            y = AGB_spatially_normalised_g_m2)) + 
      #         geom_point(na.rm = TRUE) +
      #         theme_coding() +
      #         coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 1), expand=FALSE) +
      #         labs(x = expression("Canopy height (m)"),
      #              y = expression("Dry biomass (g m"^"-2"*")"),
      #              title = "SfM") +
      #     annotate("text", x = 0.16, y = 2550, label = expression("Y==2522.3~X"),
      #              parse = TRUE, color = "black", size = 3.5, family = "serif") +
      #     annotate(geom = "text", x = 0.05, y = 2400, label = expression("r"^"2"*" = 0.899"),
      #              family = "serif", hjust = 0, parse = TRUE, size = 3.5, colour = "black") +
      #         theme(legend.position = c(0.15, 0.9)) +
      #         geom_smooth(method="lm", formula= y ~ x-1, se=TRUE, size=0.5, na.rm = TRUE))
      
      # (biomass_PF_mean <- ggplot(data = dataset,
      #                            aes(x = PF_HAG_mean,
      #                                y = AGB_spatially_normalised_g_m2)) + 
      #     geom_point(na.rm = TRUE) +
      #     theme_coding() +
      #     coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 1), expand=FALSE) +
      #     labs(x = expression("Canopy height (m)"),
      #          y = expression("Dry biomass (g m"^"-2"*")"),
      #          title = "Point Framing") +
      #     annotate("text", x = 0.16, y = 2550, label = expression("Y==3622.8~X"),
      #              parse = TRUE, color = "black", size = 3.5, family = "serif") +
      #     annotate(geom = "text", x = 0.05, y = 2400, label = expression("r"^"2"*" = 0.922"),
      #              family = "serif", hjust = 0, parse = TRUE, size = 3.5, colour = "black") +
      #     theme(legend.position = c(0.15, 0.9)) +
      #     geom_smooth(method="lm", formula= y ~ x-1, se=TRUE, size=0.5, na.rm = TRUE))
      
      (biomass_SfM_mean <- ggplot(data = dataset,
                                  aes(x = HAG_plotmean_of_cellmax_m,
                                      y = AGB_spatially_normalised_g_m2)) + 
          geom_point(na.rm = TRUE) +
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
      
      
      (biomass_PF_mean <- ggplot(data = dataset,
                            aes(x = PF_HAG_mean,
                                y = AGB_spatially_normalised_g_m2)) + 
              geom_point(na.rm = TRUE) +
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
      
      # Combine plots with ggarrange()[in ggpubr]
      biomass_plots_mean <- ggarrange(biomass_PF_mean, biomass_SfM_mean,
                                 heights = c(9,9),
                                 labels = c("(a)", "(b)"),
                                 ncol = 2, nrow = 1,
                                 align = "h")
      
      
      # Export figure
      png(filename="outputs_01_main/HAG Vs biomass mean.png", width=16, height=8, units="cm", res=400)
      plot(biomass_plots_mean)
      dev.off()
      
      # height biomass regression models
        # means
        model_PF <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$PF_HAG_mean + 0)
        summary(model_PF)
        model_PFu <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$PF_HAG_mean)
        summary(model_PFu)
        
        modelSfM <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$HAG_plotmean_of_cellmax_m + 0)
        summary(modelSfM)
        model_SfMu <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$HAG_plotmean_of_cellmax_m)
        summary(model_SfMu)
        
        #median
        # model_PF <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$PF_HAG_median + 0)
        # summary(model_PF)
        # modelSfM <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$HAG_plotmedian_of_cellmax_m + 0)
        # summary(modelSfM)
      

      
      
      mean_model_info <- list()
      mean_model_info <- c(percentile = 'mean',
                           slope = round(summary(model)$coefficients[1], digits = 2),
                           slope_error = round(summary(model)$coefficients[2], digits = 2),
                           r2 = round(summary(model)$r.squared, digits = 3),
                           p = round(summary(model)$coefficients[1,4], digits = 5),
                           df = summary(model)$df[2]
      )
      
    # Probably not used ----
    # # Create combined long dataset for plotting
    #   PlotID <- rep(dataset$PlotID, times = 2)
    #   method <- rep(c('SfM', 'PF'), each = 36)
    #   AGB <- rep(dataset$AGB_spatially_normalised_g_m2, times = 2)
    #   median_HAG <- c(dataset$HAG_plotmedian_of_cellmax_m, dataset$PF_HAG_median)
    #   mean_HAG <- c(dataset$HAG_plotmean_of_cellmax_m, dataset$PF_HAG_mean)
    #   dataset_long <- data.frame(PlotID, method, AGB, median_HAG, mean_HAG)
    #   
    #   # plot height as predictor of biomass
    #   # Median
    #   (biomass_plot_median <- ggplot(data = dataset_long,
    #                          aes(x = median_HAG,
    #                              y = AGB,
    #                              colour = method)) + 
    #       geom_point(na.rm = TRUE) +
    #       theme_coding() +
    #       scale_colour_viridis_d() +
    #       coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 1), expand=FALSE) +
    #       labs(x = expression("Median canopy height (m)"),
    #            y = expression("Dry biomass (g m"^"-2"*")")) +
    #       theme(legend.position = c(0.15, 0.9)) +
    #       geom_smooth(method="lm", formula= y ~ x-1, aes(group=method, colour=method), se=TRUE, size=0.5, na.rm = TRUE))
    # 
    #     # Export figure
    #         png(filename="outputs_01_main/HAG Vs biomass_median.png", width=12, height=12, units="cm", res=400)
    #         plot(biomass_plot_median)
    #         dev.off()
    # 
    #     # mean
    #     (biomass_plot_mean <- ggplot(data = dataset_long,
    #                                        aes(x = mean_HAG,
    #                                            y = AGB,
    #                                            colour = method)) + 
    #                 geom_point(na.rm = TRUE) +
    #                 theme_coding() +
    #                 scale_colour_viridis_d() +
    #                 coord_cartesian(ylim = c(0, max_agb), xlim = c(0, 1), expand=FALSE) +
    #                 labs(x = expression("Mean canopy height (m)"),
    #                      y = expression("Dry biomass (g m"^"-2"*")")) +
    #                 theme(legend.position = c(0.15, 0.9)) +
    #                 geom_smooth(method="lm", formula= y ~ x-1, aes(group=method, colour=method), se=TRUE, size=0.5, na.rm = TRUE))
    #         
    #         # Export figure
    #         png(filename="outputs_01_main/HAG Vs biomass_mean.png", width=12, height=12, units="cm", res=400)
    #         plot(biomass_plot_mean)
    #         dev.off()
                       
            
            
    


### HAG sensitvitiy analysis of HAG for different percentile canopy heights ### ----
    # Create long form dataframe for plotting
        temp1 <- rep(c(0, 05, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100), each = 36)
        
        temp2 = rep(dataset$AGB_spatially_normalised_g_m2, 21)
        
        temp3 <- c(dataset$HAG_plot00percentile_of_cellmax_m,
                   dataset$HAG_plot05percentiel_of_cellmax_m,
                   dataset$HAG_plot10percentile_of_cellmax_m,
                   dataset$HAG_plot15percentile_of_cellmax_m,
                   dataset$HAG_plot20percentile_of_cellmax_m,
                   dataset$HAG_plot25percentile_of_cellmax_m,
                   dataset$HAG_plot30percentile_of_cellmax_m,
                   dataset$HAG_plot35percentile_of_cellmax_m,
                   dataset$HAG_plot40percentile_of_cellmax_m,
                   dataset$HAG_plot45percentile_of_cellmax_m,
                   dataset$HAG_plot50percentile_of_cellmax_m,
                   dataset$HAG_plot55percentile_of_cellmax_m,
                   dataset$HAG_plot60percentile_of_cellmax_m,
                   dataset$HAG_plot65percentile_of_cellmax_m,
                   dataset$HAG_plot70percentile_of_cellmax_m,
                   dataset$HAG_plot75percentile_of_cellmax_m,
                   dataset$HAG_plot80percentile_of_cellmax_m,
                   dataset$HAG_plot85percentile_of_cellmax_m,
                   dataset$HAG_plot90percentiel_of_cellmax_m,
                   dataset$HAG_plot95percentile_of_cellmax_m,
                   dataset$HAG_plot100percentile_of_cellmax_m)
        
        longdata <- data.frame(temp1, temp2, temp3)
        colnames(longdata) <- c("percentile", "agb", "hag")
        remove(temp1, temp2, temp3)
        longdata$percentile <- as.factor(longdata$percentile)

    # Scatter plot
        height_mass_all_percentiles <- ggplot(data = longdata,
                    aes(x = hag,
                        y = agb,
                        colour = percentile)) +
                geom_point(na.rm = TRUE) +
                scale_colour_viridis_d() +
                labs(x = expression("canopy height (m)"),
                     y = expression("Dry biomass (g m"^"-2"*")")) +
                theme_coding() +
                geom_smooth(method="lm", formula= y ~ x-1, aes(group=percentile, colour=percentile), se=TRUE, size=0.5, na.rm = TRUE)
        height_mass_all_percentiles    

        # Export figure
            png(filename="outputs_01_main/HAG vs mass for all percentiles.png", width=16, height=16, units="cm", res=400)
            plot(height_mass_all_percentiles)
            dev.off()
        

    # Create empty dataframe for model summary statistics.
        model_summary_table <- data.frame(percentile=character(),
                                          slope=double(),
                                          slope_error=double(),
                                          r2=double(),
                                          p=double(),
                                          df=integer(),
                                          stringsAsFactors=FALSE)

    # Summarise linear regression parameters for HAG versus biomass for all percentiels.
        # Extract summary metrics.
        for (i in sort(unique(longdata$percentile))){
            subset_df<- subset(longdata, longdata==i)
        
            # define key variables
            x <- subset_df$hag
            y <- subset_df$agb
         
            # Compute and summarise model statistics
            model <- lm(y ~ x + 0) # Define model: linear with constrained intercept.
            model_info <- list()
            model_info <- c(percentile = i,
                                # as.character(unique(subset_df$percentile)),
                            slope = round(summary(model)$coefficients[1], digits = 2),
                            slope_error = round(summary(model)$coefficients[2], digits = 2),
                            r2 = round(summary(model)$r.squared, digits = 3),
                            p = round(summary(model)$coefficients[1,4], digits = 5),
                            df = summary(model)$df[2]
                )
                
            # Convert to dataframe
            model_info <- data.frame(lapply(model_info, type.convert), stringsAsFactors=FALSE)
                
            # Add new record to table.
            model_summary_table <- rbind(model_summary_table, model_info)
            
            # Remove temporary subset dataframe
            rm(subset_df)
            }

        # add parameters from meanHAG model
        mean_model_info <- data.frame(lapply(mean_model_info, type.convert), stringsAsFactors=FALSE)
        model_summary_table <- rbind(model_summary_table, mean_model_info)
        
        # prescribe order to levels
        model_summary_table$percentile <- factor(model_summary_table$percentile, as.character(model_summary_table$percentile))
        
        # plot summary model parameters.
            # r2 plot
            HAG_percentiles_r2 <- ggplot(data = model_summary_table, 
                        aes(x = percentile,
                            y = r2)) +
                geom_point(na.rm = TRUE) +
                labs(x = expression("Percentile canopy height"),
                     y = expression("r"^"2")) +
                coord_cartesian(xlim = NULL, ylim = c(0, 1)) +
                theme_coding() + 
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
                

            # slope plot
            HAG_percentiles_slope <- ggplot(data = model_summary_table, 
                        aes(x = percentile,
                            y = slope)) +
                geom_point(na.rm = TRUE) +
                coord_cartesian(xlim = NULL, ylim = c(1500, 3500)) +
                labs(x = expression("Percentile canopy height"),
                     y = expression("Slope")) +
                theme_coding() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


        # Combine plots with ggarrange()[in ggpubr]
            HAG_percentiles_summary <- ggarrange(HAG_percentiles_r2, HAG_percentiles_slope,
                                          heights = c(10,10),
                                          labels = c("(a)", "(b)"),
                                          ncol = 2,
                                          align = "h")            
            
        # Export figure
            png(filename="outputs_01_main/HAG percentiles summary.png", width=14, height=8, units="cm", res=400)
            plot(HAG_percentiles_summary)
            dev.off()
            
        
        

### Analysis of NDVI relationships ####
  # Calculate spatially normalised values for photosynthetic tissue
      dataset$leaf_biomass <- dataset$Comp2_Mass / dataset$plot_area_from_field_length  # Leaves.
      dataset$herbacious_biomass <- dataset$Comp3_Mass / dataset$plot_area_from_field_length  # Herbacious
      dataset$photosynthetic_biomass <- dataset$leaf_biomass + dataset$herbacious_biomass 
      
  # Fitting models
      model_NDVI_Biomass_018 <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$mean_NDVI_018)
      model_NDVI_Biomass_047 <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$mean_NDVI_047)
      model_NDVI_Biomass_119 <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$mean_NDVI_119)
      model_NDVI_Biomass_121 <- lm(dataset$AGB_spatially_normalised_g_m2 ~ dataset$mean_NDVI_121)
      summary(model_NDVI_Biomass_018)
      summary(model_NDVI_Biomass_047)
      summary(model_NDVI_Biomass_119)
      summary(model_NDVI_Biomass_121)
      
      model_NDVI_018_photo <- lm(dataset$photosynthetic_biomass ~ dataset$mean_NDVI_018)
      model_NDVI_047_photo <- lm(dataset$photosynthetic_biomass ~ dataset$mean_NDVI_047)
      model_NDVI_119_photo <- lm(dataset$photosynthetic_biomass ~ dataset$mean_NDVI_119)
      model_NDVI_121_photo <- lm(dataset$photosynthetic_biomass ~ dataset$mean_NDVI_121)
      summary(model_NDVI_018_photo)
      summary(model_NDVI_047_photo)
      summary(model_NDVI_119_photo)
      summary(model_NDVI_121_photo)
      
      model_NDVI_018_leaf <- lm(dataset$leaf_biomass ~ dataset$mean_NDVI_018)
      model_NDVI_047_leaf <- lm(dataset$leaf_biomass ~ dataset$mean_NDVI_047)
      model_NDVI_119_leaf <- lm(dataset$leaf_biomass ~ dataset$mean_NDVI_119)
      model_NDVI_121_leaf <- lm(dataset$leaf_biomass ~ dataset$mean_NDVI_121)
      summary(model_NDVI_018_leaf)
      summary(model_NDVI_047_leaf)
      summary(model_NDVI_119_leaf)
      summary(model_NDVI_121_leaf)
      
      
  # Plotting NDVI vs Biomass
    # Total biomass
      NDVI_vs_total_biomass_121 <- ggplot(data = dataset, 
                                    aes(x = mean_NDVI_121,
                                        y = AGB_spatially_normalised_g_m2)) + 
                            geom_point(shape = 1, na.rm = TRUE) +
                            labs(
                                # x = expression("mean NDVI"),
                                x = expression(""),
                                 y = expression("Total biomass (g m"^"-2"*")"),
                                 title = "0.121 m grain") +
                            # annotate(geom = "text", x = 0.61, y = 2500, label = expression("r"^"2"*" = 0.193"), 
                            #          family = "serif", hjust = 0, parse = TRUE, size = 3, colour = "black") +
                            geom_smooth(method='lm', formula= y~x) +
                            coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                            theme_coding()
                          
      NDVI_vs_total_biomass_119 <- ggplot(data = dataset, 
                                      aes(x = mean_NDVI_119,
                                          y = AGB_spatially_normalised_g_m2)) + 
                              geom_point(shape = 1, na.rm = TRUE) +
                              labs(
                                # x = expression("mean NDVI"),
                                x = expression(""),
                                # y = expression("Dry biomass (g m"^"-2"*")"),
                                y = expression(""),
                                title = "0.119 m grain") +
                              # annotate(geom = "text", x = 0.61, y = 2500, label = expression("r"^"2"*" = 0.160"), 
                              #          family = "serif", hjust = 0, parse = TRUE, size = 3, colour = "black") +
                              geom_smooth(method='lm', formula= y~x) +
                              coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                              theme_coding()
      
      NDVI_vs_total_biomass_047 <- ggplot(data = dataset, 
                                      aes(x = mean_NDVI_047,
                                          y = AGB_spatially_normalised_g_m2)) + 
                              geom_point(shape = 1, na.rm = TRUE) +
                              labs(
                                # x = expression("mean NDVI"),
                                x = expression(""),
                                # y = expression("Dry biomass (g m"^"-2"*")"),
                                y = expression(""),
                                title = "0.047 m grain") +
                              # annotate(geom = "text", x = 0.61, y = 2500, label = expression("r"^"2"*" = 0.137"), 
                              #          family = "serif", hjust = 0, parse = TRUE, size = 3, colour = "black") +
                              geom_smooth(method='lm', formula= y~x) +
                              coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                              theme_coding()
      
      NDVI_vs_total_biomass_018 <- ggplot(data = dataset, 
                                     aes(x = mean_NDVI_018,
                                         y = AGB_spatially_normalised_g_m2)) + 
                              geom_point(shape = 1, na.rm = TRUE) +
                              labs(
                                # x = expression("mean NDVI"),
                                x = expression(""),
                                # y = expression("Dry biomass (g m"^"-2"*")"),
                                y = expression(""),
                                title = "0.018 m grain") +
                              # annotate(geom = "text", x = 0.61, y = 2500, label = expression("r"^"2"*" = 0.138"), 
                              #          family = "serif", hjust = 0, parse = TRUE, size = 3, colour = "black") +
                              geom_smooth(method='lm', formula= y~x) +
                              coord_cartesian(ylim = c(0, max_agb), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                              theme_coding()
      
    # Photosynthetic biomass
      # Create plots 
        photo_biomass_max <- 1.12*max(dataset$photosynthetic_biomass)
        NDVI_vs_photo_biomass_121 <- ggplot(data = dataset, 
                                           aes(x = mean_NDVI_121,
                                               y = photosynthetic_biomass)) + 
                                    geom_point(shape = 1, na.rm = TRUE) +
                                    labs(
                                      # x = expression("mean NDVI"),
                                      x = expression(""),
                                      y = expression("Photosynthetic biomass (g m"^"-2"*")")
                                      # title = "0.121 m grain"
                                         ) +
                                    geom_smooth(method='lm', formula= y~x) +
                                    coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                                    theme_coding()
        
        NDVI_vs_photo_biomass_119 <- ggplot(data = dataset, 
                                             aes(x = mean_NDVI_119,
                                                 y = photosynthetic_biomass)) + 
                                    geom_point(shape = 1, na.rm = TRUE) +
                                    labs(
                                      # x = expression("mean NDVI"),
                                      x = expression(""),
                                      # y = expression("Photosynthetic dry biomass (g m"^"-2"*")"),
                                      y = expression("")
                                      # title = "0.119 m grain"
                                         ) +
                                    geom_smooth(method='lm', formula= y~x) +
                                    coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                                    theme_coding()
        
        NDVI_vs_photo_biomass_047 <- ggplot(data = dataset,
                                             aes(x = mean_NDVI_047,
                                                 y = photosynthetic_biomass)) + 
                                      geom_point(shape = 1, na.rm = TRUE) +
                                      labs(
                                        # x = expression("mean NDVI"),
                                        x = expression(""),
                                        # y = expression("Photosynthetic dry biomass (g m"^"-2"*")"),
                                        y = expression("")
                                        # title = "0.047 m grain"
                                           ) +
                                      geom_smooth(method='lm', formula= y~x) +
                                      coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                                      theme_coding()
        
        NDVI_vs_photo_biomass_018 <- ggplot(data = dataset, 
                                            aes(x = mean_NDVI_018,
                                                y = photosynthetic_biomass)) + 
                                    geom_point(shape = 1, na.rm = TRUE) +
                                    labs(
                                      # x = expression("mean NDVI"),
                                      x = expression(""),
                                      # y = expression("Photosynthetic dry biomass (g m"^"-2"*")"),
                                      y = expression("")
                                      # title = "0.018 m grain"
                                         ) +
                                    geom_smooth(method='lm', formula= y~x) +
                                    coord_cartesian(ylim = c(0, photo_biomass_max), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                                    theme_coding()  
        
    # Leaf biomass
        # Create plots
        NDVI_vs_leaf_biomass_121 <- ggplot(data = dataset, 
                                          aes(x = mean_NDVI_121,
                                              y = leaf_biomass)) + 
                                    geom_point(shape = 1, na.rm = TRUE) +
                                    labs(x = expression("mean NDVI"),
                                         y = expression("Leaf biomass (g m"^"-2"*")")
                                         # title = "0.121 m grain"
                                         ) +
                                    geom_smooth(method='lm', formula= y~x) +
                                    coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                                    theme_coding()
        
        NDVI_vs_leaf_biomass_119 <- ggplot(data = dataset, 
                                            aes(x = mean_NDVI_119,
                                                y = leaf_biomass)) + 
                                    geom_point(shape = 1, na.rm = TRUE) +
                                    labs(x = expression("mean NDVI"),
                                         # y = expression("Leaf dry biomass (g m"^"-2"*")"),
                                         y = expression("")
                                         # title = "0.119 m grain"
                                         ) +
                                    geom_smooth(method='lm', formula= y~x) +
                                    coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                                    theme_coding()
        
        NDVI_vs_leaf_biomass_047 <- ggplot(data = dataset, 
                                            aes(x = mean_NDVI_047,
                                                y = leaf_biomass)) + 
                                    geom_point(shape = 1, na.rm = TRUE) +
                                    labs(x = expression("mean NDVI"),
                                         # y = expression("Leaf dry biomass (g m"^"-2"*")"),
                                         y = expression("")
                                         # title = "0.047 m grain"
                                         ) +
                                    geom_smooth(method='lm', formula= y~x) +
                                    coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                                    theme_coding()
        
        NDVI_vs_leaf_biomass_018 <- ggplot(data = dataset, 
                                           aes(x = mean_NDVI_018,
                                               y = leaf_biomass)) + 
                                    geom_point(shape = 1, na.rm = TRUE) +
                                    labs(x = expression("mean NDVI"),
                                         # y = expression("Leaf dry biomass (g m"^"-2"*")"),
                                         y = expression("")
                                         # title = "0.018 m grain"
                                         ) +
                                    geom_smooth(method='lm', formula= y~x) +
                                    coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
                                    theme_coding()
             
    # Combine into multipanel plots 
        # Total biomass
        NDVI_total_biomass <- ggpubr::ggarrange(NDVI_vs_total_biomass_121, NDVI_vs_total_biomass_119, NDVI_vs_total_biomass_047, NDVI_vs_total_biomass_018,
                                  heights = c(10,10),
                                  labels = c("(a)", "(b)", "(c)", "(d)"),
                                  ncol = 2, nrow = 2,
                                  align = "v")
          # Export figure
            png(filename="outputs_01_main/NDVI vs total biomass.png", width=16, height=16, units="cm", res=400)
            plot(NDVI_total_biomass)
            dev.off()
        
        # Photosynthetic biomass
        NDVI_photo_biomass <- ggpubr::ggarrange(NDVI_vs_photo_biomass_121, NDVI_vs_photo_biomass_119, NDVI_vs_photo_biomass_047, NDVI_vs_photo_biomass_018,
                                       heights = c(10,10),
                                       labels = c("(a)", "(b)", "(c)", "(d)"),
                                       ncol = 2, nrow = 2,
                                       align = "v")
          # Export figure
            png(filename="outputs_01_main/NDVI vs total biomass.png", width=16, height=16, units="cm", res=400)
            plot(NDVI_photo_biomass)
            dev.off()
            
        # Leaf biomass
        NDVI_leaf_biomass <- ggpubr::ggarrange(NDVI_vs_leaf_biomass_121, NDVI_vs_leaf_biomass_119, NDVI_vs_leaf_biomass_047, NDVI_vs_leaf_biomass_018,
                                          heights = c(10,10),
                                          labels = c("(a)", "(b)", "(c)", "(d)"),
                                          ncol = 2, nrow = 2,
                                          align = "v")
          # Export figure
            png(filename="outputs_01_main/NDVI vs leaf biomass.png", width=16, height=16, units="cm", res=400)
            plot(NDVI_leaf_biomass)
            dev.off()
            
        # All together
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
                                          heights = c(10,10,10),
                                          ncol = 4, nrow = 3,
                                          align = "v",
                                          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"),
                                          font.label = list(size = 11, face = "bold")
                                          )
                                            
        
            
        # Export figure
            png(filename="outputs_01_main/NDVI vs biomass.png", width=24, height=14, units="cm", res=400)
            plot(NDVI_biomass)
            dev.off()

                      
            
# Developing non-linar models
            PowModel_Biomass_NDVI_047 <- nls(leaf_biomass ~ a*mean_NDVI_047^b, data = dataset,
                                         start = list(a =1, b =1), na.action=na.exclude)
            summary(PowModel_Biomass_NDVI_047)

            ExpModel_Biomass_NDVI_047 <- nls(leaf_biomass ~ a*exp(b*mean_NDVI_047), data = dataset,
                                             start = list(a =1, b =1), na.action=na.exclude)
            summary(ExpModel_Biomass_NDVI_047)
            
            
      # plot
            test_pow <- ggplot(data = dataset,
                                               aes(x = mean_NDVI_047,
                                                   y = leaf_biomass)) +
              geom_point(shape = 1, na.rm = TRUE) +
              labs(x = expression("mean NDVI"),
                   # y = expression("Leaf dry biomass (g m"^"-2"*")"),
                   y = expression("")
                   # title = "0.018 m grain"
              ) +
              stat_function(fun = function(x) (coef(summary(PowModel_Biomass_NDVI_047))[, "Estimate"])[1]*(x)^(coef(summary(PowModel_Biomass_NDVI_047))[, "Estimate"])[2],
                            aes(), size = 1, lty = "solid") +
              geom_smooth(method='lm', formula= y~x) +
              coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
              theme_coding()

            test_pow
            
            
            # plot
            test_exp <- ggplot(data = dataset,
                               aes(x = mean_NDVI_047,
                                   y = leaf_biomass)) +
              geom_point(shape = 1, na.rm = TRUE) +
              labs(x = expression("mean NDVI"),
                   # y = expression("Leaf dry biomass (g m"^"-2"*")"),
                   y = expression("")
                   # title = "0.018 m grain"
              ) +
              stat_function(fun = function(x) (coef(summary(ExpModel_Biomass_NDVI_047))[, "Estimate"])[1]*exp((coef(summary(ExpModel_Biomass_NDVI_047))[, "Estimate"])[2]*x),
                            aes(), size = 1, lty = "solid") +
              geom_smooth(method='lm', formula= y~x) +
              coord_cartesian(ylim = c(0, 200), xlim = c(min_ndvi, max_ndvi), expand=FALSE) +
              theme_coding()
            
            test_exp
