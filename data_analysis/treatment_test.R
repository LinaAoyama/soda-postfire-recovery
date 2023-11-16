#Treatment test 
#Set data pathway
#Load data
source("data_compiling/treatment test data compile.R")

library(ggplot2)
library(multcomp)
library(tidyverse)

#Clean data
clean_data <- data %>%
  filter(Allele1 != "NA") %>%
  filter(Treatment != "NA")

# Function for calculating standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
}

####Metric of species diversity: species richness
sp_richness_distance <- left_join(plantcomp, plot_info) #combine plant comp data and plot info

#Total sp diversity across treatments:
sprichness_trt <- sp_richness_distance %>%
  group_by(Treatment, Area) %>%
  summarize(richness = length(unique(Species)))

#Sp Diversity within populations:
sprichness_pop <- sp_richness_distance %>%
  group_by(Treatment, Plot, Area, Distance_m) %>%
  summarize(richness = length(unique(Species)))

#Summary by treatment and area
mean_sprichness <- sprichness_pop %>%
  group_by(Treatment, Area) %>%
  summarize(mean = mean(richness), se = se(richness)) 

#Plot sp diversity by treatment and area
ggplot(mean_sprichness, aes(x = Treatment, y = mean, col = Area))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1)+
  ylab(bquote(Mean~species~richness))+
  theme_bw()+
  facet_wrap(~Area, ncol = 5)+
  geom_jitter(data = sprichness_pop, aes(x = Treatment, y = richness ))

anova(lm(richness~Treatment*Area, sprichness_pop))
TukeyHSD(aov(lm(richness~Treatment*Area, sprichness_pop)))

#Sp Richness by distance from fire edge
ggplot(sprichness_pop %>% filter(Treatment%in%c("BS", "BU")), aes(x = Distance_m, y = richness, col = Treatment))+
  geom_jitter()+
  theme_bw()+
  geom_smooth(method = lm)+
  ylab(bquote(Mean~Species~Richness))+
  xlab("Distance from Fire Edge (m)")+
  annotate("text", x = 3500, y = 18, label = "Burn-Seeded (BS): y = - 0.0003x + 12.49, R2 = 0.02, p = 0.53")+
  annotate("text", x = 3500, y = 17, label = "Burn-Unseeded (BU): y = - 0.0006x + 13.10, R2 = 0.11, p = 0.22")

anova(lm(richness~Distance_m*Treatment, sprichness_pop%>%filter(Treatment%in%c("BS", "BU"))))
summary(lm(richness~Distance_m, sprichness_pop%>%filter(Treatment%in%c("BS"))))
summary(lm(richness~Distance_m, sprichness_pop%>%filter(Treatment%in%c("BU"))))

####Metrics of genetic diversity:
#Allele richness - number of unique alleles 
#Heterozygosity - proportion of heterozygotes (two different peaks) per population

#Number of plots per treatment
num_plots <- clean_data %>%
  group_by(Area, Treatment) %>%
  summarise(num_plots = length(unique(Plot)))

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

#Combine plot info with richness_pop and hetero_pop
richness_distance <- left_join(richness_pop, plot_info)
hetero_distance <- left_join(hetero_pop, plot_info)

#Richness by distance from fire edge
ggplot(richness_distance %>% filter(Treatment%in%c("BS", "BU")), aes(x = Distance_m, y = richness, col = Treatment))+
  geom_jitter()+
  theme_bw()+
  geom_smooth(method = lm)+
  ylab(bquote(Mean~Allelic~Richness))+
  xlab("Distance from Fire Edge (m)")+
  annotate("text", x = 3900, y = 20, label = "Burn-Seeded (BS): y = - 0.0003x + 6.98, R2 = 0.01, p = 0.19")+
  annotate("text", x = 3900, y = 17.5, label = "Burn-Unseeded (BU): y = - 0.0004x + 6.50, R2 = 0.03, p = 0.02")

anova(lm(richness~Distance_m*Treatment, richness_distance%>%filter(Treatment%in%c("BS", "BU"))))
summary(lm(richness~Distance_m, richness_distance%>%filter(Treatment%in%c("BS"))))
summary(lm(richness~Distance_m, richness_distance%>%filter(Treatment%in%c("BU"))))


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

#