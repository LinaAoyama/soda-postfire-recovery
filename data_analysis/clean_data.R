#Reorganize genotype data to use gstudio

#load packages
library(tidyverse)
library(dplyr)

#load data
rawdata <- read_csv(paste(datpath, "/All_Areas_combined.csv", sep = ""))
data <- rawdata

#remove NA entries (165 rows)
data <- data[!(is.na(data$Allele1)),]
data <- data[!(is.na(data$Allele2)),]
data <- data[!(is.na(data$SampleID)),]

#change "none" in Allele1 and Allele2 to 0
data$Allele1[data$Allele1 == 'none'] <- '0'
data$Allele2[data$Allele2 == 'none'] <- '0'

#combine Allele1 and Allele2 as Genotype and separate with :
data$Genotype <- paste(data$Allele1, ":", data$Allele2)

#change loci (primer) names
data$Primer[data$Primer == 'FF347548.1 HEX' | data$Primer ==  'FF347548.1'] <- 'Loc548'
data$Primer[data$Primer == 'FF343209.1 HEX' | data$Primer == 'FF343209.1'] <- 'Loc209'
data$Primer[data$Primer == 'FF340831.1 FAM' | data$Primer == 'FF340831.1'] <- 'Loc831'
data$Primer[data$Primer == 'FF340262.1 HEX' | data$Primer == 'FF340262.1'] <- 'Loc262'
data$Primer[data$Primer == 'FF343025.1 FAM' | data$Primer == 'FF343025.1'] <- 'Loc025'
data$Primer[data$Primer == 'FF344307.1 HEX' | data$Primer == 'FF344307.1'] <- 'Loc307'
data$Primer[data$Primer == 'FF347040.1 HEX' | data$Primer == 'FF347040.1'] <- 'Loc040'
data$Primer[data$Primer == 'FF344338.1 FAM' | data$Primer == 'FF344338.1'] <- 'Loc338'
data$Primer[data$Primer == 'FF344396.1 FAM' | data$Primer == 'FF344396.1'] <- 'Loc396'
data$Primer[data$Primer == 'FF342618.1 HEX' | data$Primer == 'FF342618.1'] <- 'Loc618'

#pivot_wider to change loci as columns
data_wide <- data %>% pivot_wider(names_from = Primer, values_from = Genotype)
