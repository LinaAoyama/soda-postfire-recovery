#Population structure

#load packages
library(gstudio)
library(tidyverse)
library(ggplot2)

#load data
genomic_data <- read_population("C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_cleaned.csv",
                                type = "column", locus.columns = 10:29)

#PCA of raw genotypes
matrix_raw_data <- to_mv(genomic_data[,10:19], drop.allele = TRUE)
fit.pca <- princomp(matrix_raw_data, cor = TRUE)
summary(fit.pca)
pred <- predict(fit.pca)
PCA_df <- data.frame(PC1 = pred[,1], PC2 = pred[,2], Treatment = genomic_data$Treatment,
                     Area = genomic_data$Area, Pop = genomic_data$Plot)
ggplot(PCA_df)+
  geom_point(aes(x = PC1, y = PC2, shape = Treatment, color = Area), size = 3)

#take out anatone cultivar
wild_only_data <- genomic_data[9:760,]
matrix_wild_data <- to_mv(wild_only_data[,10:19], drop.allele = TRUE)
fit.pca.wild <- princomp(matrix_wild_data, cor = TRUE)
summary(fit.pca.wild)
pred_wild <- predict(fit.pca.wild)
PCA_df_wild <- data.frame(PC1 = pred_wild[,1], PC2 = pred_wild[,2], Treatment = wild_only_data$Treatment,
                     Area = wild_only_data$Area, Pop = wild_only_data$Plot)
ggplot(PCA_df_wild)+
  geom_point(aes(x = PC1, y = PC2, shape = Treatment, color = Area), size = 3)+
  theme_bw()

#create data file in structure format 
write_population(genomic_data, "C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_cleaned.str", row.names = TRUE, mode = "structure", stratum = "Plot")
write_population(wild_only_data, "C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_no_anatone.str", row.names = TRUE, mode = "structure", stratum = "Plot")
