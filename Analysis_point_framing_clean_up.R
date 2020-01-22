#################################################
#                                               #
#  Point Framing Data Clean Up & Analysis Prep  #
#                                               #
#    Written by Gergana Daskalova 2017-02-02    #
#             gndaskalova@gmail.com             #
#                                               #
#################################################

# Packages
library(data.table)
library(plyr) # Load plyr first to avoid clashes with dplyr
library(dplyr)
library(tidyr)
library(ggplot2)

# Setting the working directory to where the csv files are and loading them in
# Careful not to have any other csvs in that folder

setwd("D:/TeamShrub/Point Framing/CSVs")  # All the csvs are in the DataStore, but I couldn't figure out how to link to it
list <- list.files(pattern = ".csv")

# A function to add a filename (i.e. plot number) column
read_csv_filename <- function(filename){
  point_framing <- read.csv(filename)
  point_framing$Plot <- filename
  point_framing
}

# Loading all the csv files into one dataframe
pointfr <- ldply(list, read_csv_filename)

# Turning SP01.csv, SP02.csv, etc., into just SP01

pointfr <- pointfr %>% separate(Plot, c("PlotN", "filename"), sep="\\.") 
pointfr <- pointfr %>% select (-filename)
head(pointfr)

# Checking the data

# Does each of the 36 sub-plots have 36 height observations? ------NO------

# There should be 1296 (36*36) HEIGHT observations

summary(pointfr$Height) # 5542 observations in total, of which 4265 NAs
5542-4265    # 1277
1296-1277 # 19 height observations missing

# Do all observations have spatial XY coordinates (multiple rows can share the same XY coordinates)
unique(pointfr$X) 
unique(pointfr$Y) # No NAs, all fine

# summary(pointfr$X) # 11 NAs in the X column
# which(is.na(pointfr$X)) # Rows 3836-3846 - went back and deleted them, they were just extra blank rows

# Do all unique locations have a height reading? -----NO------
    # Fixed the ones which have hit the surface (added height=0)
    # A few have hit Salix richardsonii, but don't have anything written down for Height
          # I imagine the salix did have some height, wasn't just zero? Especially the ones that were stems...

# A reality check on the height values (no >100 cm)
# Corrected a 196 cm tall gramminoid to 19.6cm (SP10 1x5)
# Highest height 103 cm. This is plausible.
summary(pointfr$Height)

# How many subplots are missing observations

# SP19 1x5 no height
# SP22 3x6, 4x1, 4x2, 4x3, 4x4 MISSING; 4x5 missing height for Salix r.
# SP23 1x5 missing height salix r.
# SP27 1x6 missing height, and 1x2 has TWO height measurements
# SP32 4x6 missing height
# SP33 6x1, 6x2, 6x3, 6x4, 6x5, 6x6 MISSING
# SP36 6x1, 6x2, 6x3, 6x4, 6x5, 6x6 MISSING

# Ensuring that all row/observations with a value of 'XXXothermoss' in column Tissue have 'N/A' in the Tissue variable 
         # and  'Live' in the Status variable.

# Need to make the variables characters to make the change
pointfr$Status <- as.character(pointfr$Status)
pointfr$Tissue <- as.character(pointfr$Tissue)
pointfr[pointfr$Species == 'XXXothermoss',]$Tissue <- "NA"
pointfr[pointfr$Species == 'XXXothermoss',]$Status <- "Live"

# Ensuring that all row/observations with a value of 'XXXlitter' in column C have 'Dead' in the Status variable.

pointfr[pointfr$Species == 'XXXlitter',]$Status <-  "Dead"
  
# Standardising the case (live vs Live) used in the Status column
# Making the variables factors again
pointfr$Status <- as.factor(pointfr$Status)
pointfr$Tissue <- as.factor(pointfr$Tissue)


pointfr$Status <- recode(pointfr$Status, live = "Live")
levels(pointfr$Status)
levels(pointfr$Tissue) # Need to make blank cells, N/A and NA all NAs

# Empty rows to NAs

# Define a helper function
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

## Transform all blank cells to NAs
pointfr <- pointfr %>% mutate_each(funs(empty_as_na)) 
str(pointfr)

# As a side effect of transforming spaces and blank cells to NAs, R now thinks factors are characters, changing them back
pointfr$Status <- as.factor(pointfr$Status)
pointfr$Species <- as.factor(pointfr$Species)
pointfr$Tissue <- as.factor(pointfr$Tissue)
pointfr$PlotN <- as.factor(pointfr$PlotN)
str(pointfr)

# Changing them to NA
pointfr$Status <- recode(pointfr$Status, "N/A" = "NA")
levels(pointfr$Status)

pointfr$Tissue <- recode(pointfr$Tissue, "N/A" = "NA")
levels(pointfr$Tissue)

# Standardising Equisetum and Equisetum spp. to be the same
levels(pointfr$Species)
pointfr$Species <- recode(pointfr$Species, "Equisetum spp." = "Equisetum")

write.csv(pointfr, file = "pointfr.csv") #Note file name will need to change for the analysis script. 
