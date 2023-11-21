#Estimate H0 and A with gstudio!

library(gstudio)
library(ggplot2)

#load data
genomic_data <- read_population("C:/Users/Lina/Dropbox/Academics/Projects/Soda_Fire/Data/Genotyping/Cleaned/soda_fire_genomic_data_cleaned.csv",
                                type = "column", locus.columns = 10:29)
column_class(genomic_data, class = "locus")

#visualize allele frequencies
plot(genomic_data$Loc025) #84, 87, 141, 144, 147, 150, 153, 156, 160, 164, 167
ggplot() + geom_locus(aes(x = Loc025, fill = Treatment), data = genomic_data) #Anatone: 144, 150
plot(genomic_data$Loc040) #169, 172, 175, 179, 182, 185, 188, 191, 194, 197
ggplot() + geom_locus(aes(x = Loc040, fill = Treatment), data = genomic_data) #Anatone: 182, 188, 191
plot(genomic_data$Loc209) #167, 170, 173, 179, 182, 185, 188, 191, 196, 199
ggplot() + geom_locus(aes(x = Loc209, fill = Treatment), data = genomic_data) #Anatone: 182
plot(genomic_data$Loc262) #473, 476, 479, 482, 485, 488
ggplot() + geom_locus(aes(x = Loc262, fill = Treatment), data = genomic_data) #Anatone: 225, 292
plot(genomic_data$Loc307) #162, 165, 168, 171, 174, 177, 197
ggplot() + geom_locus(aes(x = Loc307, fill = Treatment), data = genomic_data) #Anatone: 177
plot(genomic_data$Loc338) #169, 172, 175, 179, 182, 185, 188, 192, 206
ggplot() + geom_locus(aes(x = Loc338, fill = Treatment), data = genomic_data) #Anatone: 175, 182
plot(genomic_data$Loc396) #164, 177, 183, 228, 231, 234, 237, 240
ggplot() + geom_locus(aes(x = Loc396, fill = Treatment), data = genomic_data) #Anatone: 164, 228, 231
plot(genomic_data$Loc548) #39, 71, 84, 92, 95, 98, 101, 104, 107, 110, 113, 116, 119, 122, 125, 128, 132
ggplot() + geom_locus(aes(x = Loc548, fill = Treatment), data = genomic_data) #Anatone: 84 
plot(genomic_data$Loc618) #146, 159, 165, 168, 172, 178, 181, 184, 187, 190, 195, 236
ggplot() + geom_locus(aes(x = Loc618, fill = Treatment), data = genomic_data) #Anatone: 178, 187
plot(genomic_data$Loc831) #235, 238, 241, 252, 255, 258, 261, 264, 267, 272
ggplot() + geom_locus(aes(x = Loc831, fill = Treatment), data = genomic_data) #Anatone: 252, 261

freqs.loci.strata <- frequencies(genomic_data, stratum = "Treatment")
ggplot(freqs.loci.strata) +
  geom_frequencies(freqs.loci.strata) +
  facet_grid(Stratum ~.)+ theme(legend.position = "none")

#Allelic diversity
genetic_diversity(genomic_data$Loc025, mode = "A95")
