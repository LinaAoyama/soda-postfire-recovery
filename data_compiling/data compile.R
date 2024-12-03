#data compiling
library(tidyverse)
library(gstudio)

plot_info <- read_csv(paste(datpath, "/Genotyping/Cleaned/Soda_Fire_Sites.csv", sep = ""))
plot_info <- plot_info %>%
  mutate(Area = case_when(Area == 'Anatone' ~ 'Anatone', 
                          Area == 'C' ~ 'Salmon',
                          Area == 'D' ~ 'West',
                          Area == 'E' ~ 'Rockville'))

plantcomp <- read_csv(paste(datpath, "/Soda_plantcomp/Plant_composition.csv", sep = ""))

genomic_data <- read_population(paste(datpath, "/Genotyping/Cleaned/soda_fire_genomic_data_cleaned_v2.csv",sep =""),
                                type = "column", locus.columns = 10:29)
genomic_data_raw <- read_csv(paste(datpath, "/Genotyping/Cleaned/soda_fire_genomic_data_cleaned_v2.csv",sep ="")) %>%
  mutate(Area = case_when(Area == 'Anatone' ~ 'Anatone',
                          Area == 'C' ~ 'Salmon',
                          Area == 'D' ~ 'West',
                          Area == 'E' ~ 'Rockville'))

genomic_data <- genomic_data %>%
  mutate(Area = case_when(Area == 'Anatone' ~ 'Anatone',
                          Area == 'C' ~ 'Salmon',
                          Area == 'D' ~ 'West',
                          Area == 'E' ~ 'Rockville'))

psspdensity <- read_csv(paste(datpath, "/Soda_plantcomp/Bluebunch density.csv", sep = ""))
