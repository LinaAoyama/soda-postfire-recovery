#Treatment test
#Set data pathway
#Data
source("data_compiling/treatment test data compile.R")

#Clean data
clean_data <- data %>%
  filter(Allele1 != "none") 

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
  group_by(Treatment, Plot) %>%
  summarize(richness = length(unique(c(Allele1, Allele2))))

hetero_pop <- clean_data %>%
  group_by(Treatment, Plot) %>%
  summarize(heterozygosity = length(which(Allele2!="none"))/length(Allele1))

#Summary
mean_richness <- richness_pop %>%
  group_by(Treatment) %>%
  summarize(mean = mean(richness), se = se(richness))

mean_hetero <- hetero_pop %>%
  group_by(Treatment) %>%
  summarize(mean = mean(heterozygosity), se = se(heterozygosity))

#Plot them
ggplot(mean_richness, aes(x = Treatment, y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
  ylab(bquote(Mean~allelic~richness))+
  theme_bw()

ggplot(mean_hetero, aes(x = Treatment, y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, alpha = 0.9, size = 1,position = position_dodge(width = 0.5))+
  ylab(bquote(Mean~heterozygosity))+
  theme_bw()
