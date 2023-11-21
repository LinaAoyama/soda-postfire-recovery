#Reorganize genomic data to use gstudio

#load packages
library(tidyverse)
library(dplyr)

#load data
source("data_compiling/data compile.R")
data <- rawdata

#remove NA entries (165 rows)
data <- data[!(is.na(data$Allele1)),]
data <- data[!(is.na(data$Allele2)),]
data <- data[!(is.na(data$SampleID)),]

#change "none" in Allele1 and Allele2 to 0
data$Allele1[data$Allele1 == 'none'] <- '0'
data$Allele2[data$Allele2 == 'none'] <- '0'

#combine Allele1 and Allele2 as Genotype and separate with :
data$Genotype <- as.character(paste(data$Allele1, ":", data$Allele2))

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
data_wide <- data %>% 
  select(SampleID, Treatment, Plot, Rep, Area, Primer, Genotype) %>%
  pivot_wider(names_from = Primer, values_from = Genotype) %>%
  separate(Loc548, c("Loc548", "Loc548A2")) %>% #split genotype to two columns
  separate(Loc209, c("Loc209", "Loc209A2")) %>%
  separate(Loc831, c("Loc831", "Loc831A2")) %>%
  separate(Loc262, c("Loc262", "Loc262A2")) %>%
  separate(Loc025, c("Loc025", "Loc025A2")) %>%
  separate(Loc307, c("Loc307", "Loc307A2")) %>%
  separate(Loc040, c("Loc040", "Loc040A2")) %>%
  separate(Loc338, c("Loc338", "Loc338A2")) %>%
  separate(Loc396, c("Loc396", "Loc396A2")) %>%
  separate(Loc618, c("Loc618", "Loc618A2")) 
  
#combine genomic data and plot info 
full_data <- full_join(plot_info, data_wide) 

#Bin alleles by loci: +2 bp (round down), at least 3 bp apart
binned_genomic_data <- full_data

#84, 87, 141, 144, 147, 150, 153, 156, 160, 164, 167
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '85'| binned_genomic_data$Loc025 == '86'] <- '84'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '88'] <- '87'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '142'| binned_genomic_data$Loc025 == '143'] <- '141'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '145'| binned_genomic_data$Loc025 == '146'] <- '144'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '148'| binned_genomic_data$Loc025 == '149'] <- '147'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '151'| binned_genomic_data$Loc025 == '152'] <- '150'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '154'| binned_genomic_data$Loc025 == '155'] <- '153'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '157'| binned_genomic_data$Loc025 == '158'] <- '156'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '161'| binned_genomic_data$Loc025 == '162'] <- '160'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '165'] <- '164'
binned_genomic_data$Loc025[binned_genomic_data$Loc025 == '169'] <- '167'

binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '85'| binned_genomic_data$Loc025A2 == '86'] <- '84'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '88'] <- '87'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '142'| binned_genomic_data$Loc025A2 == '143'] <- '141'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '145'| binned_genomic_data$Loc025A2 == '146'] <- '144'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '148'| binned_genomic_data$Loc025A2 == '149'] <- '147'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '151'| binned_genomic_data$Loc025A2 == '152'] <- '150'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '154'| binned_genomic_data$Loc025A2 == '155'] <- '153'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '157'| binned_genomic_data$Loc025A2 == '158'] <- '156'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '161'| binned_genomic_data$Loc025A2 == '162'] <- '160'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '165'] <- '164'
binned_genomic_data$Loc025A2[binned_genomic_data$Loc025A2 == '169'] <- '167'

#169, 172, 175, 179, 182, 185, 188, 191, 194, 197
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '170'| binned_genomic_data$Loc040 == '171'] <- '169'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '173'| binned_genomic_data$Loc040 == '174'] <- '172'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '177'] <- '175'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '180'| binned_genomic_data$Loc040 == '181'] <- '179'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '183'| binned_genomic_data$Loc040 == '184'] <- '182'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '186'| binned_genomic_data$Loc040 == '187'] <- '185'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '189'| binned_genomic_data$Loc040 == '190'] <- '188'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '192'| binned_genomic_data$Loc040 == '193'] <- '191'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '195'| binned_genomic_data$Loc040 == '196'] <- '194'
binned_genomic_data$Loc040[binned_genomic_data$Loc040 == '198'| binned_genomic_data$Loc040 == '199'] <- '197'

binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '170'| binned_genomic_data$Loc040A2 == '171'] <- '169'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '173'| binned_genomic_data$Loc040A2 == '174'] <- '172'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '177'] <- '175'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '180'| binned_genomic_data$Loc040A2 == '181'] <- '179'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '183'| binned_genomic_data$Loc040A2 == '184'] <- '182'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '186'| binned_genomic_data$Loc040A2 == '187'] <- '185'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '189'| binned_genomic_data$Loc040A2 == '190'] <- '188'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '192'| binned_genomic_data$Loc040A2 == '193'] <- '191'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '195'| binned_genomic_data$Loc040A2 == '196'] <- '194'
binned_genomic_data$Loc040A2[binned_genomic_data$Loc040A2 == '198'| binned_genomic_data$Loc040A2 == '199'] <- '197'

#167, 170, 173, 179, 182, 185, 188, 191, 196, 199
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == ' 0'] <- '0'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == 'noe'] <- '0'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '168'] <- '167'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '171' | binned_genomic_data$Loc209 == '172'] <- '170'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '175' ] <- '173'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '180' | binned_genomic_data$Loc209 == '181'] <- '179'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '183' | binned_genomic_data$Loc209 == '184'] <- '182'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '187'] <- '185'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '189'] <- '188'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '192' | binned_genomic_data$Loc209 == '193'] <- '191'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '197'] <- '196'
binned_genomic_data$Loc209[binned_genomic_data$Loc209 == '201'] <- '199'

binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == ' 0'] <- '0'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == 'noe'] <- '0'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '168'] <- '167'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '171' | binned_genomic_data$Loc209A2 == '172'] <- '170'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '175' ] <- '173'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '180' | binned_genomic_data$Loc209A2 == '181'] <- '179'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '183' | binned_genomic_data$Loc209A2 == '184'] <- '182'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '187'] <- '185'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '189'] <- '188'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '192' | binned_genomic_data$Loc209A2 == '193'] <- '191'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '197'] <- '196'
binned_genomic_data$Loc209A2[binned_genomic_data$Loc209A2 == '201'] <- '199'

#473, 476, 479, 482, 485, 488
binned_genomic_data$Loc262[binned_genomic_data$Loc262 == '474' | binned_genomic_data$Loc262 == '475'] <- '473'
binned_genomic_data$Loc262[binned_genomic_data$Loc262 == '478'] <- '476'
binned_genomic_data$Loc262[binned_genomic_data$Loc262 == '480' | binned_genomic_data$Loc262 == '481'] <- '479'
binned_genomic_data$Loc262[binned_genomic_data$Loc262 == '483' | binned_genomic_data$Loc262 == '484'] <- '482'
binned_genomic_data$Loc262[binned_genomic_data$Loc262 == '486' | binned_genomic_data$Loc262 == '487'] <- '485'
binned_genomic_data$Loc262[binned_genomic_data$Loc262 == '490'] <- '488'

binned_genomic_data$Loc262A2[binned_genomic_data$Loc262A2 == '474' | binned_genomic_data$Loc262A2 == '475'] <- '473'
binned_genomic_data$Loc262A2[binned_genomic_data$Loc262A2 == '478'] <- '476'
binned_genomic_data$Loc262A2[binned_genomic_data$Loc262A2 == '480' | binned_genomic_data$Loc262A2 == '481'] <- '479'
binned_genomic_data$Loc262A2[binned_genomic_data$Loc262A2 == '483' | binned_genomic_data$Loc262A2 == '484'] <- '482'
binned_genomic_data$Loc262A2[binned_genomic_data$Loc262A2 == '486' | binned_genomic_data$Loc262A2 == '487'] <- '485'
binned_genomic_data$Loc262A2[binned_genomic_data$Loc262A2 == '490'] <- '488'

#162, 165, 168, 171, 174, 177, 197
binned_genomic_data$Loc307[binned_genomic_data$Loc307 == '163' | binned_genomic_data$Loc307 == '164'] <- '162'
binned_genomic_data$Loc307[binned_genomic_data$Loc307 == '166' | binned_genomic_data$Loc307 == '167'] <- '165'
binned_genomic_data$Loc307[binned_genomic_data$Loc307 == '169' | binned_genomic_data$Loc307 == '170'] <- '168'
binned_genomic_data$Loc307[binned_genomic_data$Loc307 == '172' | binned_genomic_data$Loc307 == '173'] <- '171'
binned_genomic_data$Loc307[binned_genomic_data$Loc307 == '176'] <- '174'
binned_genomic_data$Loc307[binned_genomic_data$Loc307 == '178' | binned_genomic_data$Loc307 == '179'] <- '177'
binned_genomic_data$Loc307[binned_genomic_data$Loc307 == '199'] <- '197'

binned_genomic_data$Loc307A2[binned_genomic_data$Loc307A2 == '163' | binned_genomic_data$Loc307A2 == '164'] <- '162'
binned_genomic_data$Loc307A2[binned_genomic_data$Loc307A2 == '166' | binned_genomic_data$Loc307A2 == '167'] <- '165'
binned_genomic_data$Loc307A2[binned_genomic_data$Loc307A2 == '169' | binned_genomic_data$Loc307A2 == '170'] <- '168'
binned_genomic_data$Loc307A2[binned_genomic_data$Loc307A2 == '172' | binned_genomic_data$Loc307A2 == '173'] <- '171'
binned_genomic_data$Loc307A2[binned_genomic_data$Loc307A2 == '176'] <- '174'
binned_genomic_data$Loc307A2[binned_genomic_data$Loc307A2 == '178' | binned_genomic_data$Loc307A2 == '179'] <- '177'
binned_genomic_data$Loc307A2[binned_genomic_data$Loc307A2 == '199'] <- '197'

#169, 172, 175, 179, 182, 185, 188, 192, 206
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '170' | binned_genomic_data$Loc338 == '171'] <- '169'
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '173' | binned_genomic_data$Loc338 == '174'] <- '172'
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '176' | binned_genomic_data$Loc338 == '177'] <- '175'
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '180' ] <- '179'
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '183' | binned_genomic_data$Loc338 == '184'] <- '182'
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '186' | binned_genomic_data$Loc338 == '187'] <- '185'
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '189' | binned_genomic_data$Loc338 == '190'] <- '188'
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '194' ] <- '192'
binned_genomic_data$Loc338[binned_genomic_data$Loc338 == '207' ] <- '206'

binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '170' | binned_genomic_data$Loc338A2 == '171'] <- '169'
binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '173' | binned_genomic_data$Loc338A2 == '174'] <- '172'
binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '176' | binned_genomic_data$Loc338A2 == '177'] <- '175'
binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '180' ] <- '179'
binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '183' | binned_genomic_data$Loc338A2 == '184'] <- '182'
binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '186' | binned_genomic_data$Loc338A2 == '187'] <- '185'
binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '189' | binned_genomic_data$Loc338A2 == '190'] <- '188'
binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '194' ] <- '192'
binned_genomic_data$Loc338A2[binned_genomic_data$Loc338A2 == '207' ] <- '206'

#164, 177, 183, 228, 231, 234, 237, 240
binned_genomic_data$Loc396[binned_genomic_data$Loc396 == '179'] <- '177'
binned_genomic_data$Loc396[binned_genomic_data$Loc396 == '185'] <- '183'
binned_genomic_data$Loc396[binned_genomic_data$Loc396 == '229' | binned_genomic_data$Loc396 == '230'] <- '228'
binned_genomic_data$Loc396[binned_genomic_data$Loc396 == '232' | binned_genomic_data$Loc396 == '233'] <- '231'
binned_genomic_data$Loc396[binned_genomic_data$Loc396 == '236'] <- '234'
binned_genomic_data$Loc396[binned_genomic_data$Loc396 == '238' | binned_genomic_data$Loc396 == '239'] <- '237'
binned_genomic_data$Loc396[binned_genomic_data$Loc396 == '241' | binned_genomic_data$Loc396 == '242'] <- '240'

binned_genomic_data$Loc396A2[binned_genomic_data$Loc396A2 == '179'] <- '177'
binned_genomic_data$Loc396A2[binned_genomic_data$Loc396A2 == '185'] <- '183'
binned_genomic_data$Loc396A2[binned_genomic_data$Loc396A2 == '229' | binned_genomic_data$Loc396A2 == '230'] <- '228'
binned_genomic_data$Loc396A2[binned_genomic_data$Loc396A2 == '232' | binned_genomic_data$Loc396A2 == '233'] <- '231'
binned_genomic_data$Loc396A2[binned_genomic_data$Loc396A2 == '236'] <- '234'
binned_genomic_data$Loc396A2[binned_genomic_data$Loc396A2 == '238' | binned_genomic_data$Loc396A2 == '239'] <- '237'
binned_genomic_data$Loc396A2[binned_genomic_data$Loc396A2 == '241' | binned_genomic_data$Loc396A2 == '242'] <- '240'

#39, 71, 84, 92, 95, 98, 101, 104, 107, 110, 113, 116, 119, 122, 125, 128, 132
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '40'] <- '39'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '73'] <- '71'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '85'] <- '84'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '94'] <- '92'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '96' | binned_genomic_data$Loc548 == '97'] <- '95'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '99' | binned_genomic_data$Loc548 == '100'] <- '98'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '102' | binned_genomic_data$Loc548 == '103'] <- '101'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '105' | binned_genomic_data$Loc548 == '106'] <- '104'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '108' | binned_genomic_data$Loc548 == '109'] <- '107'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '111'] <- '110'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '115'] <- '113'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '117' | binned_genomic_data$Loc548 == '118'] <- '116'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '120' | binned_genomic_data$Loc548 == '121'] <- '119'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '123' | binned_genomic_data$Loc548 == '124'] <- '122'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '126' | binned_genomic_data$Loc548 == '127'] <- '125'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '129' | binned_genomic_data$Loc548 == '130'] <- '128'
binned_genomic_data$Loc548[binned_genomic_data$Loc548 == '134' ] <- '132'

binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '40'] <- '39'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '73'] <- '71'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '85'] <- '84'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '94'] <- '92'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '96' | binned_genomic_data$Loc548A2 == '97'] <- '95'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '99' | binned_genomic_data$Loc548A2 == '100'] <- '98'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '102' | binned_genomic_data$Loc548A2 == '103'] <- '101'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '105' | binned_genomic_data$Loc548A2 == '106'] <- '104'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '108' | binned_genomic_data$Loc548A2 == '109'] <- '107'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '111'] <- '110'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '115'] <- '113'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '117' | binned_genomic_data$Loc548A2 == '118'] <- '116'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '120' | binned_genomic_data$Loc548A2 == '121'] <- '119'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '123' | binned_genomic_data$Loc548A2 == '124'] <- '122'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '126' | binned_genomic_data$Loc548A2 == '127'] <- '125'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '129' | binned_genomic_data$Loc548A2 == '130'] <- '128'
binned_genomic_data$Loc548A2[binned_genomic_data$Loc548A2 == '134' ] <- '132'

#146, 159, 165, 168, 172, 178, 181, 184, 187, 190, 195, 236
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '147'] <- '146'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '161'] <- '159'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '166' | binned_genomic_data$Loc618 == '167'] <- '165'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '169'] <- '168'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '173' | binned_genomic_data$Loc618 == '174'] <- '172'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '179' | binned_genomic_data$Loc618 == '180'] <- '178'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '182' | binned_genomic_data$Loc618 == '183'] <- '181'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '185' | binned_genomic_data$Loc618 == '186'] <- '184'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '188' | binned_genomic_data$Loc618 == '189'] <- '187'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '191'] <- '190'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '196'] <- '195'
binned_genomic_data$Loc618[binned_genomic_data$Loc618 == '237'] <- '236'

binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '147'] <- '146'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '161'] <- '159'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '166' | binned_genomic_data$Loc618A2 == '167'] <- '165'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '169'] <- '168'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '173' | binned_genomic_data$Loc618A2 == '174'] <- '172'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '179' | binned_genomic_data$Loc618A2 == '180'] <- '178'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '182' | binned_genomic_data$Loc618A2 == '183'] <- '181'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '185' | binned_genomic_data$Loc618A2 == '186'] <- '184'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '188' | binned_genomic_data$Loc618A2 == '189'] <- '187'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '191'] <- '190'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '196'] <- '195'
binned_genomic_data$Loc618A2[binned_genomic_data$Loc618A2 == '237'] <- '236'

#235, 238, 241, 252, 255, 258, 261, 264, 267, 272
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '236' | binned_genomic_data$Loc831 == '237'] <- '235'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '239' | binned_genomic_data$Loc831 == '240'] <- '238'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '242' | binned_genomic_data$Loc831 == '243'] <- '241'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '253' | binned_genomic_data$Loc831 == '254'] <- '252'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '256' | binned_genomic_data$Loc831 == '257'] <- '255'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '259' | binned_genomic_data$Loc831 == '260'] <- '258'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '262' | binned_genomic_data$Loc831 == '263'] <- '261'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '265' | binned_genomic_data$Loc831 == '266'] <- '264'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '268' | binned_genomic_data$Loc831 == '269'] <- '267'
binned_genomic_data$Loc831[binned_genomic_data$Loc831 == '273' ] <- '272'

binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '236' | binned_genomic_data$Loc831A2 == '237'] <- '235'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '239' | binned_genomic_data$Loc831A2 == '240'] <- '238'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '242' | binned_genomic_data$Loc831A2 == '243'] <- '241'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '253' | binned_genomic_data$Loc831A2 == '254'] <- '252'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '256' | binned_genomic_data$Loc831A2 == '257'] <- '255'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '259' | binned_genomic_data$Loc831A2 == '260'] <- '258'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '262' | binned_genomic_data$Loc831A2 == '263'] <- '261'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '265' | binned_genomic_data$Loc831A2 == '266'] <- '264'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '268' | binned_genomic_data$Loc831A2 == '269'] <- '267'
binned_genomic_data$Loc831A2[binned_genomic_data$Loc831A2 == '273' ] <- '272'

write.csv(binned_genomic_data, "C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_cleaned.csv", row.names = FALSE)

