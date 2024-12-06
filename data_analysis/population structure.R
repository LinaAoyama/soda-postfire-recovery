#Population structure

#load packages
library(gstudio)
library(tidyverse)
library(ggplot2)
library(adegenet)
library(ggpubr)
library(vegan)
library(reshape2)

#load data
source("data_compiling/data compile.R")

#PCA of raw genotypes
#matrix_raw_data <- to_mv(genomic_data[,10:19], drop.allele = FALSE)
#matrix_raw_data_test <- genomic_data[,10:19]
matrix_raw_data <- genomic_data_raw[,10:29]
fit.pca <- princomp(matrix_raw_data, cor = TRUE)
summary(fit.pca)
pred <- predict(fit.pca)
PCA_df <- data.frame(PC1 = pred[,1], PC2 = pred[,2], Treatment = genomic_data_raw$Treatment,
                     Area = genomic_data_raw$Area, Pop = genomic_data_raw$Plot)
fig_pca_all <- ggplot(PCA_df)+
                  geom_point(aes(x = PC1, y = PC2, color = Area, shape = Treatment), size = 3)+
                  theme(text = element_text(size=15),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                        axis.title = element_text(size = 15),
                        legend.position = "right")+
                  xlab("PC1 (27.3%)")+
                  ylab("PC2 (13.4%)")+
                  scale_shape_manual(labels = c("Anatone", "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), values = c(16, 3, 17, 15))

#PERMANOVA
# genomic_data_long <- genomic_data %>%
#   pivot_longer(cols = Loc025:Loc831, names_to = "Locus", values_to = "Alleles")
# genomic_data_long$Alleles <- as.character(genomic_data_long$Alleles)
# genomic_data_long <- as.data.frame(genomic_data_long)
# genomic_data_long <- separate(genomic_data_long, col = Alleles, into = c('Allele1', 'Allele2'), sep = ':')
# #genomic_data_long <- genomic_data_long %>%
# #  pivot_longer(cols = Allele1:Allele2, names_to = "Allele", values_to = "values")
# genomic_data_long$Allele1 <- as.numeric(genomic_data_long$Allele1)
# permanova <- adonis(Allele1~Area, data = genomic_data_long, perm = 99, method = "euclidean") #PERMANOVA results: Significant treatment effect p = 0.001
genomic_data_env <- genomic_data_raw[,1:9]
adonis2(matrix_raw_data ~ Area*Treatment, data = genomic_data_env)  #PERMANOVA results: Significant area and treatment effects p = 0.001

#Try removing loci 307 and 396 (>50% missing values)
matrix_raw_data_truncated <- genomic_data_raw[,10:29]
matrix_raw_data_truncated <- matrix_raw_data_truncated[,-c(9,10, 13, 14)]
fit.pca_truncated <- princomp(matrix_raw_data_truncated, cor = TRUE)
summary(fit.pca_truncated)
pred_truncated <- predict(fit.pca_truncated)
PCA_df_truncated <- data.frame(PC1 = pred_truncated[,1], PC2 = pred_truncated[,2], Treatment = genomic_data_raw$Treatment,
                     Area = genomic_data_raw$Area, Pop = genomic_data_raw$Plot)
fig_pca_all_truncated <- ggplot(PCA_df_truncated)+
  geom_point(aes(x = PC1, y = PC2, color = Area, shape = Treatment), size = 3)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.position = "right")+
  xlab("PC1 (29.9%)")+
  ylab("PC2 (14.6%)")+
  scale_shape_manual(labels = c("Anatone", "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), values = c(16, 3, 17, 15))


#PCA subset by area to see Anatone and seeded trts better
rockville_raw_data <- genomic_data_raw %>%
  filter(Area == "Rockville"| Area == "Anatone")
matrix_rockville <- rockville_raw_data[,10:29]
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
                              legend.position = "bottom")+
                        xlab("PC1 (18.0%)")+
                        ylab("PC2 (13.3%)")+
                        scale_color_manual(values = c("#F8766D", 
                                                      "#f7b2d8",
                                                      "#999999",
                                                      "#b2ffc3"),                                                      
                                           labels = c("Anatone", "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded" ))+
                        scale_shape_manual(labels = c( "Anatone", "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), 
                                           values = c( 16, 3, 17, 15))+
                        geom_point(data = PCA_df_rockville %>%
                                     filter( Area == "Anatone"), 
                                   aes(x = PC1, y = PC2, color = Treatment, shape = Treatment), size = 3)

                        #annotate("text", label = "Rockville", size = 5, x = -6, y = )

salmon_raw_data <- genomic_data_raw %>%
  filter(Area == "Salmon"| Area == "Anatone")
matrix_salmon <- salmon_raw_data[,10:29]
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
                          legend.position = "bottom")+
                    xlab("PC1 (26.6%)")+
                    ylab("PC2 (14.6%)")+
                    scale_color_manual(values = c("#F8766D",
                                                  "#f7b2d8",
                                                  "#999999",
                                                  "#b2ffc3"), 
                                       labels = c("Anatone", "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"))+
                    scale_shape_manual(labels = c("Anatone", "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), 
                                       values = c(16, 3, 17, 15))

west_raw_data <- genomic_data_raw %>%
  filter(Area == "West"| Area == "Anatone")
matrix_west <- west_raw_data[,10:29]
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
                        legend.position = "bottom")+
                  xlab("PC1 (26.2%)")+
                  ylab("PC2 (15.1%)")+
                  scale_color_manual(values = c("#F8766D",
                                                "#f7b2d8",
                                                "#999999",
                                                "#b2ffc3"), 
                                     labels = c("Anatone", "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"))+
                  scale_shape_manual(labels = c("Anatone", "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), 
                                     values = c(16, 3, 17, 15))

ggarrange(fig_pca_all, ggarrange(fig_pca_rockville, fig_pca_salmon, fig_pca_west, labels = c("b)", "c)", "d)"), 
                                 ncol =3, nrow =1, common.legend = TRUE), labels = c("a)"), ncol = 1, nrow =2, heights = c(1.5, 1))
######################################################################
#PCA without Anatone
#All areas
data_no_anatone <- genomic_data_raw%>%
                        filter( Area != "Anatone")
matrix_raw_data_no_anatone <- data_no_anatone[,10:29]
fit.pca_no_anatone <- princomp(matrix_raw_data_no_anatone, cor = TRUE)
summary(fit.pca_no_anatone)
pred_no_anatone <- predict(fit.pca_no_anatone)
PCA_df_no_anatone <- data.frame(PC1 = pred_no_anatone[,1], PC2 = pred_no_anatone[,2], Treatment = data_no_anatone$Treatment,
                     Area = data_no_anatone$Area, Pop = data_no_anatone$Plot)
fig_pca_no_anatone <- ggplot(PCA_df_no_anatone)+
  geom_point(aes(x = PC1, y = PC2, color = Area, shape = Treatment), size = 3)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.position = "right")+
  xlab("PC1 (26.7%)")+
  ylab("PC2 (13.5%)")+
  scale_color_manual(values = c("#7CAE00",
                                "#00BFC4",
                                "#C77CFF"), 
                     labels = c( "Rockville", "Salmon", "West"))+
  scale_shape_manual(labels = c("Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), values = c(3, 17, 15))

genomic_data_env_no_anatone <- data_no_anatone[,1:9]
adonis2(matrix_raw_data_no_anatone ~ Area*Treatment, data = genomic_data_env_no_anatone)  #PERMANOVA results: Significant area and treatment effects p = 0.001

#Rockville no anatone
rockville_raw_data_no_anatone <- data_no_anatone %>%
  filter(Area == "Rockville")
matrix_rockville_no_anatone <- rockville_raw_data_no_anatone[,10:29]
fit.pca.rockville_no_anatone <- princomp(matrix_rockville_no_anatone, cor = TRUE)
summary(fit.pca.rockville_no_anatone)
pred.rockville_no_anatone <- predict(fit.pca.rockville_no_anatone)
PCA_df_rockville_no_anatone <- data.frame(PC1 = pred.rockville_no_anatone[,1], PC2 = pred.rockville_no_anatone[,2], Treatment = rockville_raw_data_no_anatone$Treatment,
                               Area = rockville_raw_data_no_anatone$Area, Pop = rockville_raw_data_no_anatone$Plot)
fig_pca_rockville_no_anatone <- ggplot(PCA_df_rockville_no_anatone)+
  geom_point(aes(x = PC1, y = PC2, color = Treatment, shape = Treatment), size = 3)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.position = "bottom")+
  xlab("PC1 (18.3%)")+
  ylab("PC2 (13.5%)")+
  scale_color_manual(values = c(
                                "#f7b2d8",
                                "#999999",
                                "#b2ffc3"),                                                      
                     labels = c( "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded" ))+
  scale_shape_manual(labels = c( "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), 
                     values = c( 3, 17, 15))

adonis2(matrix_rockville_no_anatone ~ Treatment, data = genomic_data_env_no_anatone%>%filter(Area == "Rockville"))  #PERMANOVA results: Significant treatment effects p = 0.001

#Salmon no anatone
salmon_raw_data_no_anatone <- data_no_anatone %>%
  filter(Area == "Salmon")
matrix_salmon_no_anatone <- salmon_raw_data_no_anatone[,10:29]
fit.pca.salmon_no_anatone <- princomp(matrix_salmon_no_anatone, cor = TRUE)
summary(fit.pca.salmon_no_anatone)
pred.salmon_no_anatone <- predict(fit.pca.salmon_no_anatone)
PCA_df_salmon_no_anatone <- data.frame(PC1 = pred.salmon_no_anatone[,1], PC2 = pred.salmon_no_anatone[,2], Treatment = salmon_raw_data_no_anatone$Treatment,
                            Area = salmon_raw_data_no_anatone$Area, Pop = salmon_raw_data_no_anatone$Plot)
fig_pca_salmon_no_anatone <- ggplot(PCA_df_salmon_no_anatone)+
  geom_point(aes(x = PC1, y = PC2, color = Treatment, shape = Treatment), size = 3)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.position = "bottom")+
  xlab("PC1 (25.3%)")+
  ylab("PC2 (14.7%)")+
  scale_color_manual(values = c(
                                "#f7b2d8",
                                "#999999",
                                "#b2ffc3"), 
                     labels = c( "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"))+
  scale_shape_manual(labels = c( "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), 
                     values = c( 3, 17, 15))

adonis2(matrix_salmon_no_anatone ~ Treatment, data = genomic_data_env_no_anatone%>%filter(Area == "Salmon"))  #PERMANOVA results: Significant treatment effects p = 0.001

#West no anatone
west_raw_data_no_anatone <- data_no_anatone %>%
  filter(Area == "West")
matrix_west_no_anatone <- west_raw_data_no_anatone[,10:29]
fit.pca.west_no_anatone <- princomp(matrix_west_no_anatone, cor = TRUE)
summary(fit.pca.west_no_anatone)
pred.west_no_anatone <- predict(fit.pca.west_no_anatone)
PCA_df_west_no_anatone <- data.frame(PC1 = pred.west_no_anatone[,1], PC2 = pred.west_no_anatone[,2], Treatment = west_raw_data_no_anatone$Treatment,
                          Area = west_raw_data_no_anatone$Area, Pop = west_raw_data_no_anatone$Plot)
fig_pca_west_no_anatone <- ggplot(PCA_df_west_no_anatone)+
  geom_point(aes(x = PC1, y = PC2, color = Treatment, shape = Treatment), size = 3)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.position = "bottom")+
  xlab("PC1 (24.8%)")+
  ylab("PC2 (15.2%)")+
  scale_color_manual(values = c(
                                "#f7b2d8",
                                "#999999",
                                "#b2ffc3"), 
                     labels = c( "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"))+
  scale_shape_manual(labels = c( "Burned-Seeded", "Burned-Unseeded", "Unburned-Unseeded"), 
                     values = c( 3, 17, 15))

adonis2(matrix_west_no_anatone ~ Treatment, data = genomic_data_env_no_anatone%>%filter(Area == "West"))  #PERMANOVA results: Significant treatment effects p = 0.001

ggarrange(fig_pca_no_anatone, ggarrange(fig_pca_rockville_no_anatone, fig_pca_salmon_no_anatone, fig_pca_west_no_anatone, labels = c("b)", "c)", "d)"), 
                                 ncol =3, nrow =1, common.legend = TRUE), labels = c("a)"), ncol = 1, nrow =2, heights = c(1.5, 1))
###################################################################################
#PCA subset by treatments
BS_raw_data <- genomic_data_raw %>%
  filter(Treatment == "BS"| Treatment == "Anatone")
matrix_BS <- BS_raw_data[,10:29]
fit.pca.BS <- princomp(matrix_BS, cor = TRUE)
summary(fit.pca.BS)
pred.BS <- predict(fit.pca.BS)
PCA_df_BS <- data.frame(PC1 = pred.BS[,1], PC2 = pred.BS[,2], Treatment = BS_raw_data$Treatment,
                            Area = BS_raw_data$Area, Pop = BS_raw_data$Plot)
fig_pca_BS <- ggplot(PCA_df_BS)+
  geom_point(aes(x = PC1, y = PC2, color = Area, shape = Treatment), size = 3)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.position = "none")+
  xlab("PC1 (32.9%)")+
  ylab("PC2 (14.7%)")+
  scale_shape_manual(labels = c("Anatone", "Burned-Seeded"), values = c(16, 3))

BU_raw_data <- genomic_data_raw %>%
  filter(Treatment == "BU"| Treatment == "Anatone")
matrix_BU <- BU_raw_data[,10:29]
fit.pca.BU <- princomp(matrix_BU, cor = TRUE)
summary(fit.pca.BU)
pred.BU <- predict(fit.pca.BU)
PCA_df_BU <- data.frame(PC1 = pred.BU[,1], PC2 = pred.BU[,2], Treatment = BU_raw_data$Treatment,
                        Area = BU_raw_data$Area, Pop = BU_raw_data$Plot)
fig_pca_BU <- ggplot(PCA_df_BU)+
  geom_point(aes(x = PC1, y = PC2, color = Area, shape = Treatment), size = 3)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.position = "none")+
  xlab("PC1 (32.9%)")+
  ylab("PC2 (14.7%)")+
  scale_shape_manual(labels = c("Anatone", "Burned-Unseeded"), values = c(16, 17))
  
UU_raw_data <- genomic_data_raw %>%
  filter(Treatment == "UU"| Treatment == "Anatone")
matrix_UU <- UU_raw_data[,10:29]
fit.pca.UU <- princomp(matrix_UU, cor = TRUE)
summary(fit.pca.UU)
pred.UU <- predict(fit.pca.UU)
PCA_df_UU <- data.frame(PC1 = pred.UU[,1], PC2 = pred.UU[,2], Treatment = UU_raw_data$Treatment,
                        Area = UU_raw_data$Area, Pop = UU_raw_data$Plot)
fig_pca_UU <- ggplot(PCA_df_UU)+
  geom_point(aes(x = PC1, y = PC2, color = Area, shape = Treatment), size = 3)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 15),
        legend.position = "none")+
  xlab("PC1 (25.0%)")+
  ylab("PC2 (15.8%)")+
  scale_shape_manual(labels = c("Anatone", "Unburned-Unseeded"), values = c(16, 15))

ggarrange(fig_pca_all, ggarrange(fig_pca_BS, fig_pca_BU, fig_pca_UU, labels = c("b)", "c)", "d)"), 
                                 ncol =3, nrow =1), labels = c("a)"), ncol = 1, nrow =2, heights = c(1.5, 1))

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

# #REMOVE loci with >50% missing data (L262, L307, L396)
# reduced_data <-  subset(genomic_data, select = -c(Loc262,Loc307, Loc396) )
#   
# matrix_red_data <- to_mv(reduced_data[,10:16], drop.allele = TRUE)
# fit.pca.red <- princomp(matrix_red_data, cor = TRUE)
# summary(fit.pca.red)
# pred.red <- predict(fit.pca.red)
# PCA_df_red <- data.frame(PC1 = pred.red[,1], PC2 = pred.red[,2], Treatment = reduced_data$Treatment,
#                      Area = reduced_data$Area, Pop = reduced_data$Plot)
# fig_pca_red <- ggplot(PCA_df_red)+
#   geom_point(aes(x = PC1, y = PC2, shape = Treatment, color = Area), size = 3)+
#   theme(text = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         axis.title = element_text(size = 15),
#         legend.position = "right")+
#   xlab("PC1 (3.9%)")+
#   ylab("PC2 (2.8%)")
# 
# salmon_red_data <- reduced_data %>%
#   filter(Area == "Salmon"| Area == "Anatone")
# matrix_salmon_red <- to_mv(salmon_red_data[,10:16], drop.allele = TRUE)
# fit.pca.salmon.red <- princomp(matrix_salmon_red, cor = TRUE)
# summary(fit.pca.salmon.red)
# pred.salmon.red <- predict(fit.pca.salmon.red)
# PCA_df_salmon_red <- data.frame(PC1 = pred.salmon.red[,1], PC2 = pred.salmon.red[,2], Treatment = salmon_red_data$Treatment,
#                             Area = salmon_red_data$Area, Pop = salmon_red_data$Plot)
# fig_pca_salmon_red <- ggplot(PCA_df_salmon_red)+
#   geom_point(aes(x = PC1, y = PC2, color = Treatment, shape = Treatment), size = 3)+
#   theme(text = element_text(size=15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         axis.title = element_text(size = 15),
#         legend.position = "bottom")+
#   xlab("PC1 (7.2%)")+
#   ylab("PC2 (5.7%)")+
#   scale_color_manual(values = c("#F8766D",
#                                 "#f7b2d8",
#                                 "#999999",
#                                 "#b2ffc3"))

#Remove anatone from genomic_data
genomic_data_no_anatone <- genomic_data%>%
  filter( Area != "Anatone")

#create data file in STRUCTURE format 
write_population(genomic_data, "C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_cleaned.str", row.names = TRUE, mode = "structure", stratum = "Plot")
write_population(genomic_data_no_anatone, "C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_no_anatone.str", row.names = TRUE, mode = "structure", stratum = "Plot")

#Run STRUCTURE externally on soda_fire_genomic_data_cleaned.str
#100k burn-in period, 100k iterations, K=1 to 6, run 10 times

#Run STRUCTURE Harvester externally on STRUCTURE output
#Visualize clumpp files with pophelper: https://www.royfrancis.com/pophelperShiny/index.html
# install dependencies and remotes
#install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes",
#                   "colourpicker","DT","highcharter","htmlwidgets","magrittr",
#                   "markdown","RColorBrewer","shiny","shinyAce","shinyBS",
#                   "shinythemes","shinyWidgets","viridisLite","writexl"),
#                repos = "http://cran.us.r-project.org")

# install pophelper package from GitHub
#remotes::install_github('royfrancis/pophelper')

# install the package from GitHub
#remotes::install_github('royfrancis/pophelperShiny')


# load library for use R version 4.3
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
