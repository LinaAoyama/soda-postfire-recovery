#data compiling
library(tidyverse)
data <- read_csv(paste(datpath, "/All_Areas_combined_plotnames_adjusted.csv", sep = ""))
plot_info <- read_csv(paste(datpath, "/Soda_Fire_Sites.csv", sep = ""))
plantcomp <- read_csv("C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Soda_plantcomp/Plant_composition.csv")


#count number of unique values
#data_9018 <- read_csv(paste(datpath, "/Area E/9018.csv", sep = ""))
library(dbplyr)
#print(data_9018 %>%
#  pivot_longer(Allele1:Allele2)%>%
#  group_by(value) %>%
#  summarize(n_unique = length((value))) %>%
#  ungroup(), n = 21)


