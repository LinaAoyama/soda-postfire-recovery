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

#Number of plots per treatment
num_plots <- clean_data %>%
  group_by(Area, Treatment) %>%
  summarise(num_plots = length(unique(Plot)))

#Metrics of genetic diversity:
#Allele richness - number of unique alleles 
#Heterozygosity - proportion of heterozygotes (two different peaks) per population

#Total diversity across treatments:
richness_trt <- clean_data %>%
  group_by(Treatment, Area) %>%
  summarize(richness = length(unique(c(Allele1, Allele2))))

hetero_trt <- clean_data %>%
  group_by(Treatment, Area) %>%
  summarize(heterozygosity = length(which(Allele2!="none"))/length(Allele1))
  
#Diversity within populations:
richness_pop <- clean_data %>%
  group_by(Treatment, Plot, Primer, Area) %>%
  summarize(richness = length(unique(c(Allele1, Allele2))))

hetero_pop <- clean_data %>%
  group_by(Treatment, Plot, Primer, Area) %>%
  summarize(heterozygosity = length(which(Allele2!="none"))/length(Allele1))

#Summary by treatment and area
mean_richness <- richness_pop %>%
  group_by(Treatment, Area) %>%
  summarize(mean = mean(richness), se = se(richness)) 
mean_hetero <- hetero_pop %>%
  group_by(Treatment, Area) %>%
  summarize(mean = mean(heterozygosity), se = se(heterozygosity))

#anova test: Effect of Treatment and Area on genetic diversity
anova(lm(richness~Treatment*Area, richness_pop)) # sig area diff but not treatment diff 
TukeyHSD(aov(lm(richness~Treatment*Area, richness_pop)))
anova(lm(heterozygosity~Treatment*Area, hetero_pop)) # sig area diff but not treatment diff 
TukeyHSD(aov(lm(heterozygosity~Treatment*Area, hetero_pop)))

#Plot diversity by treatment and area
ggplot(mean_richness, aes(x = Treatment, y = mean, col = Area))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1)+
  ylab(bquote(Mean~allelic~richness))+
  theme_bw()+
  #facet_wrap(~Area, ncol = 5)+
  geom_jitter(data = richness_pop, aes(x = Treatment, y = richness ))

ggplot(mean_hetero, aes(x = Treatment, y = mean, col = Area))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1)+
  ylab(bquote(Mean~heterozygosity))+
  theme_bw()+
  #facet_wrap(~Area, ncol = 5)+
  geom_jitter(data = hetero_pop, aes(x = Treatment, y = heterozygosity ))

#Summary by each primer
mean_richness_primer <- richness_pop %>%
  group_by(Treatment, Primer, Area) %>%
  summarize(mean = mean(richness), se = se(richness))

mean_hetero_primer <- hetero_pop %>%
  group_by(Treatment, Primer, Area) %>%
  summarize(mean = mean(heterozygosity), se = se(heterozygosity))

#Plot diversity by each primer
ggplot(mean_richness_primer, aes(x = Treatment, y = mean, col = Area))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1)+
  ylab(bquote(Mean~allelic~richness))+
  theme_bw()+
  facet_wrap(~Primer, ncol = 5)+
  geom_jitter(data = richness_pop, aes(x = Treatment, y = richness ))

ggplot(mean_hetero_primer, aes(x = Treatment, y = mean, col = Area))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1)+
  ylab(bquote(Mean~heterozygosity))+
  theme_bw()+
  facet_wrap(~Primer, ncol = 5)+
  geom_jitter(data = hetero_pop, aes(x = Treatment, y = heterozygosity))

