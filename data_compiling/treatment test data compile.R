#data compiling
library(tidyverse)
data <- read_csv(paste(datpath, "/Area D and E/Analysis Area E half D adjusted combined.csv", sep = ""))


#count number of unique values
#data_9018 <- read_csv(paste(datpath, "/Area E/9018.csv", sep = ""))
library(dbplyr)
#print(data_9018 %>%
#  pivot_longer(Allele1:Allele2)%>%
#  group_by(value) %>%
#  summarize(n_unique = length((value))) %>%
#  ungroup(), n = 21)


