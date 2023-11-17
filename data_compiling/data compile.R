#data compiling
library(tidyverse)
rawdata <- read_csv(paste(datpath, "/All_Areas_combined_plotnames_adjusted.csv", sep = ""))
plot_info <- read_csv(paste(datpath, "/Soda_Fire_Sites.csv", sep = ""))
plantcomp <- read_csv("C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Soda_plantcomp/Plant_composition.csv")


