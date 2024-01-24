#Population structure

#load packages
library(gstudio)
library(tidyverse)
library(ggplot2)
library(adegenet)
library(ggpubr)

#load data
source("data_compiling/data compile.R")

#PCA of raw genotypes
matrix_raw_data <- to_mv(genomic_data[,10:19], drop.allele = TRUE)
fit.pca <- princomp(matrix_raw_data, cor = TRUE)
summary(fit.pca)
pred <- predict(fit.pca)
PCA_df <- data.frame(PC1 = pred[,1], PC2 = pred[,2], Treatment = genomic_data$Treatment,
                     Area = genomic_data$Area, Pop = genomic_data$Plot)
fig_pca_all <- ggplot(PCA_df)+
                  geom_point(aes(x = PC1, y = PC2, shape = Treatment, color = Area), size = 3)+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                        axis.title = element_text(size = 15),
                        legend.position = "right")+
                  xlab("PC1 (3.9%)")+
                  ylab("PC2 (6.9%)")

#PCA subset by area to see Anatone and seeded trts better
rockville_raw_data <- genomic_data %>%
  filter(Area == "Rockville"| Area == "Anatone")
matrix_rockville <- to_mv(rockville_raw_data[,10:19], drop.allele = TRUE)
fit.pca.rockville <- princomp(matrix_rockville, cor = TRUE)
summary(fit.pca.rockville)
pred.rockville <- predict(fit.pca.rockville)
PCA_df_rockville <- data.frame(PC1 = pred.rockville[,1], PC2 = pred.rockville[,2], Treatment = rockville_raw_data$Treatment,
                     Area = rockville_raw_data$Area, Pop = rockville_raw_data$Plot)
fig_pca_rockville <- ggplot(PCA_df_rockville)+
                        geom_point(aes(x = PC1, y = PC2, color = Treatment, shape = Treatment), size = 3)+
                        theme(text = element_text(size=15),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              axis.line = element_line(colour = "black"),
                              panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                              axis.title = element_text(size = 15),
                              legend.position = "right")+
                        xlab("PC1 (3.8%)")+
                        ylab("PC2 (6.9%)")+
                        scale_color_manual(values = c("#F8766D",
                                                      "#7CAE00",
                                                      "#999999",
                                                      "#E69F00"))

salmon_raw_data <- genomic_data %>%
  filter(Area == "Salmon"| Area == "Anatone")
matrix_salmon <- to_mv(salmon_raw_data[,10:19], drop.allele = TRUE)
fit.pca.salmon <- princomp(matrix_salmon, cor = TRUE)
summary(fit.pca.salmon)
pred.salmon <- predict(fit.pca.salmon)
PCA_df_salmon <- data.frame(PC1 = pred.salmon[,1], PC2 = pred.salmon[,2], Treatment = salmon_raw_data$Treatment,
                               Area = salmon_raw_data$Area, Pop = salmon_raw_data$Plot)
fig_pca_salmon <- ggplot(PCA_df_salmon)+
                    geom_point(aes(x = PC1, y = PC2, color = Treatment, shape = Treatment), size = 3)+
                    theme(text = element_text(size=15),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                          axis.title = element_text(size = 15),
                          legend.position = "right")+
                    xlab("PC1 (7.4%)")+
                    ylab("PC2 (13.1%)")+
                    scale_color_manual(values = c("#F8766D",
                                                  "#00BFC4",
                                                  "#999999",
                                                  "#E69F00"))

west_raw_data <- genomic_data %>%
  filter(Area == "West"| Area == "Anatone")
matrix_west <- to_mv(west_raw_data[,10:19], drop.allele = TRUE)
fit.pca.west <- princomp(matrix_west, cor = TRUE)
summary(fit.pca.west)
pred.west <- predict(fit.pca.west)
PCA_df_west <- data.frame(PC1 = pred.west[,1], PC2 = pred.west[,2], Treatment = west_raw_data$Treatment,
                            Area = west_raw_data$Area, Pop = west_raw_data$Plot)
fig_pca_west <- ggplot(PCA_df_west)+
                  geom_point(aes(x = PC1, y = PC2, color = Treatment, shape = Treatment), size = 3)+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                        axis.title = element_text(size = 15),
                        legend.position = "right")+
                  xlab("PC1 (9.1%)")+
                  ylab("PC2 (16.8%)")+
                  scale_color_manual(values = c("#F8766D",
                                                "#C77CFF",
                                                "#999999",
                                                "#E69F00"))

ggarrange(fig_pca_all, fig_pca_rockville, fig_pca_salmon, fig_pca_west, labels = c("a)", "b)", "c)", "d)"), ncol =2, nrow =2)
# #take out anatone cultivar
# wild_only_data <- genomic_data[9:760,]
# matrix_wild_data <- to_mv(wild_only_data[,10:19], drop.allele = TRUE)
# fit.pca.wild <- princomp(matrix_wild_data, cor = TRUE)
# summary(fit.pca.wild)
# pred_wild <- predict(fit.pca.wild)
# PCA_df_wild <- data.frame(PC1 = pred_wild[,1], PC2 = pred_wild[,2], Treatment = wild_only_data$Treatment,
#                      Area = wild_only_data$Area, Pop = wild_only_data$Plot)
# ggplot(PCA_df_wild)+
#   geom_point(aes(x = PC1, y = PC2, shape = Treatment, color = Area), size = 3)+
#   theme_bw()

#create data file in STRUCTURE format 
write_population(genomic_data, "C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_cleaned.str", row.names = TRUE, mode = "structure", stratum = "Plot")
#write_population(wild_only_data, "C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_no_anatone.str", row.names = TRUE, mode = "structure", stratum = "Plot")

#Run STRUCTURE externally on soda_fire_genomic_data_cleaned.str
#100k burn-in period, 100k iterations, K=1 to 6, run 10 times

#Run STRUCTURE Harvester externally on STRUCTURE output
#Visualize clumpp files with pophelper: https://www.royfrancis.com/pophelperShiny/index.html
# install dependencies and remotes
install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes",
                   "colourpicker","DT","highcharter","htmlwidgets","magrittr",
                   "markdown","RColorBrewer","shiny","shinyAce","shinyBS",
                   "shinythemes","shinyWidgets","viridisLite","writexl"),
                 repos = "http://cran.us.r-project.org")

# install pophelper package from GitHub
remotes::install_github('royfrancis/pophelper')

# install the package from GitHub
remotes::install_github('royfrancis/pophelperShiny')

# load library for use
library(pophelperShiny)

# launch app
runPophelper()

###Are closer populations more similar to each other than more distant populations?
#Calculate pair-wise Gst (Genetic differentiation between populations)
Gst <- genetic_structure(genomic_data, stratum = "Plot", mode = "Gst", pairwise = TRUE) #0-1; large values mean more differentiated (no sharing of genetic mat)
#Calculate pair-wise distance matrix
plot_info[order(plot_info$Plot),]
GeoDistance <- dist(plot_info[-38,4:5], method = "euclidean",
                    diag = TRUE, upper = TRUE) #omit UU13
#Plot pairwise Gst and distance


map <- population_map(plot_info)
ggplot2::ggmap(map)+
  geom_point(aes(x = Longitude, y = Latitude, col = Area), data = genomic_data)
