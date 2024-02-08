#Estimate allelic richness and heterozygosity with gstudio!

#load packages
library(gstudio)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

#load data
source("data_compiling/data compile.R")
column_class(genomic_data, class = "locus")

# this is a function for calculating standard error
se<-function(x){
  sd(x, na.rm = TRUE)/sqrt(length(x))
} 

#visualize allele frequencies
plot(genomic_data$Loc025) #84, 87, 141, 144, 147, 150, 153, 156, 160, 164, 167
L025 <- ggplot() + geom_locus(aes(x = Loc025, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L025", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 144, 150 OK
plot(genomic_data$Loc040) #169, 172, 175, 179, 182, 185, 188, 191, 194, 197
L040 <- ggplot() + geom_locus(aes(x = Loc040, fill = Treatment), data = genomic_data)+ 
  annotate("text", label = "L040", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 182, 188, 191 OK
plot(genomic_data$Loc209) #167, 170, 173, 179, 182, 185, 188, 191, 196, 199
L209 <- ggplot() + geom_locus(aes(x = Loc209, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L209", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 182 OK
plot(genomic_data$Loc262) #473, 476, 479, 482, 485, 488 >50% missing data
L262 <- ggplot() + geom_locus(aes(x = Loc262, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L262", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 225, 292 Weird
plot(genomic_data$Loc307) #162, 165, 168, 171, 174, 177, 197 >50% missing data
L307 <- ggplot() + geom_locus(aes(x = Loc307, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L307", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 177 OK
plot(genomic_data$Loc338) #169, 172, 175, 179, 182, 185, 188, 192, 206
L338 <- ggplot() + geom_locus(aes(x = Loc338, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L338", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 175, 182 Weird
plot(genomic_data$Loc396) #164, 177, 183, 228, 231, 234, 237, 240 >50% missing data
L396 <- ggplot() + geom_locus(aes(x = Loc396, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L396", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 164, 228, 231 Weird
plot(genomic_data$Loc548) #39, 71, 84, 92, 95, 98, 101, 104, 107, 110, 113, 116, 119, 122, 125, 128, 132
L548 <- ggplot() + geom_locus(aes(x = Loc548, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L548", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 84 Weird
plot(genomic_data$Loc618) #146, 159, 165, 168, 172, 178, 181, 184, 187, 190, 195, 236
L618 <- ggplot() + geom_locus(aes(x = Loc618, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L618", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 178, 187 OK
plot(genomic_data$Loc831) #235, 238, 241, 252, 255, 258, 261, 264, 267, 272
L831 <- ggplot() + geom_locus(aes(x = Loc831, fill = Treatment), data = genomic_data)+
  annotate("text", label = "L831", x = 2, y = 1, size = 5)+
  ylim(c(0,1))#Anatone: 252, 261 OK

ggarrange(L025, L040, L209, L262, L307, L338, L396, L548, L618, L831,
          ncol = 2, nrow = 5, common.legend = TRUE)

# freqs.loci.strata <- frequencies(genomic_data, stratum = "Treatment")
# ggplot(freqs.loci.strata) +
#   geom_frequencies(freqs.loci.strata) +
#   facet_grid(Stratum ~.)+ theme(legend.position = "none")

#number of alleles per treatment per loci 
#Allelic richness (total number of alleles)
Arichness.diversity <- genetic_diversity(genomic_data, stratum = "Treatment", mode = "A")
colnames(Arichness.diversity) <- c('Treatment', 'Locus', 'A')
A_plot <- ggplot(Arichness.diversity, aes(y = A, x = Treatment, fill = Treatment))+
  geom_col()+
  facet_wrap(~Locus, ncol = 2)+
  scale_fill_manual(values = c("#F8766D",
                                "#f7b2d8",
                                "#999999",
                                "#b2ffc3"))+
  theme(legend.position="top")+
  ylab("Total number of alleles (A)")

#Allelic richness (Effective Number of Alleles)
A.diversity <- genetic_diversity(genomic_data, stratum = "Plot", mode = "Ae")
colnames(A.diversity) <- c('Plot', 'Locus', 'Ae')

#Observed heterozygosity 
Ho.diversity <- genetic_diversity(genomic_data, stratum = "Plot", mode = "Ho")
colnames(Ho.diversity) <- c('Plot', 'Locus', 'Ho')

#Expected heterozygosity
He.diversity <- genetic_diversity(genomic_data, stratum = "Plot", mode = "He")
colnames(He.diversity) <- c('Plot', 'Locus', 'He')

#Combine diversity metrics in one table
genetic.diversity <- left_join(A.diversity, Ho.diversity) %>%
  left_join(., He.diversity) %>%
  left_join(., plot_info) 

genetic.diversity.pop <- genetic.diversity %>%
  group_by(Plot, Treatment, Area) %>%
  summarize(mean_Ae = mean(Ae), 
            se_Ae = se(Ae),
            mean_Ho = mean(Ho), 
            se_Ho = se(Ho),
            mean_He = mean(He),
            se_He = se(He)) %>%
  left_join(., plot_info) 

anova(lm(Ae ~ Treatment * Area, genetic.diversity %>% filter(Treatment != "Anatone")))
TukeyHSD(aov(Ae ~ Treatment * Area, genetic.diversity))
anova(lm(Ho ~ Treatment * Area, genetic.diversity %>% filter(Treatment != "Anatone")))
TukeyHSD(aov(Ho ~  Area, genetic.diversity))
anova(lm(He ~ Treatment * Area, genetic.diversity %>% filter(Treatment != "Anatone")))
TukeyHSD(aov(He ~ Treatment * Area, genetic.diversity))
#Summary by area
mean.genetic.diveristy <- genetic.diversity.pop %>%
  group_by(Area) %>%
  summarize(Ae = mean(mean_Ae), se_Ae = se(mean_Ae),
            Ho = mean(mean_Ho), se_Ho = se(mean_Ho),
            He = mean(mean_He), se_He = se(mean_He))

#Summary by treatment and area
mean.genetic.diveristy <- genetic.diversity.pop %>%
  group_by(Treatment, Area) %>%
  summarize(Ae = mean(mean_Ae), se_Ae = se(mean_Ae),
            Ho = mean(mean_Ho), se_Ho = se(mean_Ho),
            He = mean(mean_He), se_He = se(mean_He))

#Plot diversity by treatment and area
Ae_plot <- ggplot(mean.genetic.diveristy, aes(x = Treatment, y = Ae, col = Area))+
  geom_point()+
  geom_errorbar(aes(ymin = Ae-se_Ae, ymax = Ae+se_Ae), width = 0.2, alpha = 0.9, size = 1)+
  geom_jitter(data = genetic.diversity.pop %>% filter(Treatment != "Anatone"), aes(x = Treatment, y = mean_Ae))+
  ylab(bquote(Ae))+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "top")
  #facet_wrap(~Area, ncol = 5)+
  
Ho_plot <- ggplot(mean.genetic.diveristy, aes(x = Treatment, y = Ho, col = Area))+
  geom_point()+
  geom_errorbar(aes(ymin = Ho-se_Ho, ymax = Ho+se_Ho), width = 0.2, alpha = 0.9, size = 1)+
  geom_jitter(data = genetic.diversity.pop %>% filter(Treatment != "Anatone"), aes(x = Treatment, y = mean_Ho))+
  ylab(bquote(Ho))+
  ylim(0, 1)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 12))
  #facet_wrap(~Area, ncol = 5)+
  
He_plot <- ggplot(mean.genetic.diveristy, aes(x = Treatment, y = He, col = Area))+
  geom_point()+
  geom_errorbar(aes(ymin = He-se_He, ymax = He+se_He), width = 0.2, alpha = 0.9, size = 1)+
  geom_jitter(data = genetic.diversity.pop %>% filter(Treatment != "Anatone"), aes(x = Treatment, y = mean_He))+
  ylab(bquote(He))+
  ylim(0, 1)+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank())
  #ggtitle("Mean expected heterozygosity")
  #facet_wrap(~Area, ncol = 5)+
  
ggarrange(A_plot, ggarrange(Ae_plot,He_plot,Ho_plot, ncol = 1, common.legend = TRUE, align = "hv", 
          legend = "top", labels = c("b)", "c)", "d)")), ncol = 2, labels = c("a)"))

Fis <- genetic_diversity(genomic_data, stratum = "Plot", mode = "Fis")
colnames(Fis) <- c('Plot', 'Locus', 'Fis')
Fis.pop <- left_join(Fis, plot_info)  %>%
  filter(Plot != "Anatone") %>%
  group_by(Area) %>%
  summarize(mean_Fis = mean(Fis, na.rm = TRUE),
            se_Fis = se(Fis)) 
Fis.overall <- Fis %>%
  summarize(mean_Fis = mean(Fis, na.rm = TRUE),
            se_Fis = se(Fis))
# #Summary by each primer
# mean.genetic.diveristy.primer <- genetic.diversity %>%
#   group_by(Treatment, Area, Locus) %>%
#   summarize(mean_Ae = mean(Ae), se_Ae = se(Ae),
#             mean_Ho = mean(Ho), se_Ho = se(Ho),
#             mean_He = mean(He), se_He = se(He))
# 
# #Plot diversity by each primer
# ggplot(mean.genetic.diveristy.primer, aes(x = Treatment, y = mean_Ae, col = Area))+
#   geom_point()+
#   geom_errorbar(aes(ymin = mean_Ae-se_Ae, ymax = mean_Ae+se_Ae), width = 0.2, alpha = 0.9, size = 1)+
#   ylab(bquote(Mean~allelic~richness))+
#   theme_bw()+
#   facet_wrap(~Locus, ncol = 5)+
#   geom_jitter(data = genetic.diversity, aes(x = Treatment, y = Ae))
# 
# ggplot(mean.genetic.diveristy.primer, aes(x = Treatment, y = mean_Ho, col = Area))+
#   geom_point()+
#   geom_errorbar(aes(ymin = mean_Ho-se_Ho, ymax = mean_Ho+se_Ho), width = 0.2, alpha = 0.9, size = 1)+
#   ylab(bquote(Mean~Observed~Heterozygosity))+
#   theme_bw()+
#   facet_wrap(~Locus, ncol = 5)+
#   geom_jitter(data = genetic.diversity, aes(x = Treatment, y = Ho))
# 
# ggplot(mean.genetic.diveristy.primer, aes(x = Treatment, y = mean_He, col = Area))+
#   geom_point()+
#   geom_errorbar(aes(ymin = mean_He-se_He, ymax = mean_He+se_He), width = 0.2, alpha = 0.9, size = 1)+
#   ylab(bquote(Mean~Expected~Heterozygosity))+
#   theme_bw()+
#   facet_wrap(~Locus, ncol = 5)+
#   geom_jitter(data = genetic.diversity, aes(x = Treatment, y = He))

#Richness by distance from fire edge
Ae_distance <- ggplot(genetic.diversity.pop %>% filter(Treatment%in%c("BS", "BU")), aes(x = Distance_m, y = mean_Ae, col = Area))+
  geom_jitter()+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank())+
  geom_smooth(method = lm)+
  ylab(bquote(Mean~Allelic~Richness))+
  xlab("Distance from Fire Edge (m)")+
  facet_grid(~Treatment)+
  scale_color_manual(values = c("#7CAE00",
                                "#00BFC4",
                                "#C77CFF"))
  #annotate("text", x = 3900, y = 20, label = "Burn-Seeded (BS): y = - 0.0003x + 6.98, R2 = 0.01, p = 0.19")+
  #annotate("text", x = 3900, y = 17.5, label = "Burn-Unseeded (BU): y = - 0.0004x + 6.50, R2 = 0.03, p = 0.02")

anova(lm(Ae~Distance_m*Treatment, genetic.diversity%>%filter(Treatment%in%c("BS", "BU"))))
summary(lm(Ae~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "Rockville")))
summary(lm(Ae~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "Salmon")))
summary(lm(Ae~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "West")))
summary(lm(Ae~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "Rockville")))
summary(lm(Ae~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "Salmon")))
summary(lm(Ae~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "West")))

Ho_distance <- ggplot(genetic.diversity.pop %>% filter(Treatment%in%c("BS", "BU")), aes(x = Distance_m, y = mean_Ho, col = Area))+
  geom_jitter()+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 12))+
  geom_smooth(method = lm)+
  ylab(bquote(Mean~Obs~Heterozygosity))+
  xlab("Distance from Fire Edge (m)")+
  facet_grid(~Treatment)+
  scale_color_manual(values = c("#7CAE00",
                                "#00BFC4",
                                "#C77CFF"))

summary(lm(Ho~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "Rockville")))
summary(lm(Ho~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "Salmon")))
summary(lm(Ho~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "West")))
summary(lm(Ho~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "Rockville")))
summary(lm(Ho~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "Salmon")))
summary(lm(Ho~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "West")))

He_distance <- ggplot(genetic.diversity.pop %>% filter(Treatment%in%c("BS", "BU")), aes(x = Distance_m, y = mean_He, col = Area))+
  geom_jitter()+
  theme(text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank())+
  geom_smooth(method = lm)+
  ylab(bquote(Mean~Exp~Heterozygosity))+
  xlab("Distance from Fire Edge (m)")+
  facet_grid(~Treatment)+
  scale_color_manual(values = c("#7CAE00",
                                "#00BFC4",
                                "#C77CFF"))

summary(lm(He~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "Rockville")))
summary(lm(He~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "Salmon")))
summary(lm(He~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BS"))%>%filter(Area == "West")))
summary(lm(He~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "Rockville")))
summary(lm(He~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "Salmon")))
summary(lm(He~Distance_m, genetic.diversity%>%filter(Treatment%in%c("BU"))%>%filter(Area == "West")))

ggarrange(Ae_distance, He_distance, Ho_distance, ncol = 1, align = "hv",
          common.legend = TRUE, legend = "right", labels = c("a)", "b)", "c)"))

#Fst
Fst <- Fst(genomic_data, stratum = "Plot")
