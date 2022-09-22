#Treatment test 
#Set data pathway
#Data in Area E
source("data_compiling/treatment test data compile.R")

library(ggplot2)
library(multcomp)

#Clean data
clean_data <- data %>%
  filter(Allele1 != "NA") 

# Function for calculating standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
}

#Metrics of genetic diversity:
#Allele richness - number of unique alleles 
#Heterozygosity - proportion of heterozygotes (two different peaks) per population

#Total diversity across treatments:
richness_trt <- clean_data %>%
  group_by(Treatment) %>%
  summarize(richness = length(unique(c(Allele1, Allele2))))

hetero_trt <- clean_data %>%
  group_by(Treatment) %>%
  summarize(heterozygosity = length(which(Allele2!="none"))/length(Allele1))
  
#Diversity within populations:
richness_pop <- clean_data %>%
  group_by(Treatment, Plot, Primer) %>%
  summarize(richness = length(unique(c(Allele1, Allele2))))

hetero_pop <- clean_data %>%
  group_by(Treatment, Plot, Primer) %>%
  summarize(heterozygosity = length(which(Allele2!="none"))/length(Allele1))

#Summary
mean_richness <- richness_pop %>%
  group_by(Treatment, Primer) %>%
  summarize(mean = mean(richness), se = se(richness))

mean_hetero <- hetero_pop %>%
  group_by(Treatment, Primer) %>%
  summarize(mean = mean(heterozygosity), se = se(heterozygosity))

#Plot them
ggplot(mean_richness, aes(x = Treatment, y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
  ylab(bquote(Mean~allelic~richness))+
  theme_bw()+
  facet_wrap(~Primer, ncol = 5)+
  geom_jitter(data = richness_pop, aes(x = Treatment, y = richness ))

ggplot(mean_hetero, aes(x = Treatment, y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
  ylab(bquote(Mean~heterozygosity))+
  theme_bw()+
  facet_wrap(~Primer, ncol = 5)+
  geom_jitter(data = hetero_pop, aes(x = Treatment, y = heterozygosity))

#anova test
anova(lm(heterozygosity~Treatment, hetero_pop))

TukeyHSD(aov(lm(heterozygosity~Treatment, hetero_pop)))
