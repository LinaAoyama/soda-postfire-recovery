#Load data
source("data_compiling/data compile.R")

#Load library
library(tidyverse)
library(ggplot2)

# this is a function for calculating standard error
se<-function(x){
  sd(x, na.rm = TRUE)/sqrt(length(x))
} 

#PSSP density by treatment
summary_density <- psspdensity %>%
  group_by(Treatment) %>%
  summarise(mean = mean(PSSP6_density), 
            se = se(PSSP6_density))

summary(aov(PSSP6_density ~ Treatment, psspdensity))
TukeyHSD(aov(PSSP6_density ~ Treatment, psspdensity))
