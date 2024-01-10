#UPDATED ANALYSES




library(ggplot2)
library(ggfortify)
require(FactoMineR)
require(factoextra)
require(tidyr)
require(dplyr)
require(MASS)
require(reshape2)
require(cowplot)
require(tidyverse)
require(forcats)
library(RColorBrewer)
library(poppr)


# Resource: https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
#           http://rstudio-pubs-static.s3.amazonaws.com/53162_cd16ee63c24747459ccd180f69f07810.html

setwd("C:/Users/nikki/Michigan State University/PSM.Hausbecklab - Nikki Lukasko - Nikki Lukasko/Botrytis/Molecular/Plates123 Population Genomics/Population Structure Analysis/Plink PCA")

#----------PCAs -------------


### IMPORTANT: LOOK AT .EIGENVAL FILES ON HPCC FOR THE ACTUAL PCA PERCENTAGES (below)

# For no maf filter (1,122,956 var): 43.07, 19.01, 15.43, 14.87, 13.71
# For 0.05 maf filter (81,220 var): 37.41, 24.20, 22.81, 20.44, 15.56
# for 0.1 maf filter (48,188 var): 36.89, 27.16, 25.69, 20.12, 17.11



#Basic PCA
pca_table <- read.table("BcinMAF05_PCA_results.eigenvec", header = TRUE, comment.char = "")
plot(pca_table[, c("PC1", "PC2", "PC3", "PC4", "PC5")])


#Detailed PCA
library(readxl)
Metadata_eigenvec <- read_excel("Metadata_eigenvec_current005and01.xlsx", 
                                col_types = c("text", "text", "text", 
                                              "text", "text", "text", "text", "numeric", 
                                              "text", "text", "text", "numeric", 
                                              "text", "text", "text", "text",
                                              "text", "text", "text", 
                                              "text", "text", "numeric", 
                                              "numeric", "numeric", "numeric",
                                              "numeric", "numeric", "numeric",
                                              "numeric", "numeric", "numeric",
                                              "numeric", "numeric", "numeric",
                                              "numeric", "numeric", "numeric",
                                              "numeric", "numeric", "numeric", "numeric",
                                              "numeric", "numeric"))
View(Metadata_eigenvec)

plot(Metadata_eigenvec[, c("PC1_maf005", "PC2_maf005", "PC3_maf005", "PC4_maf005", "PC5_maf005"), color="Group S" ])


dfnomaf <- Metadata_eigenvec[25:29]
pca_nomaf <- prcomp(dfnomaf, scale. = TRUE)
autoplot(pca_nomaf)
summary(pca_nomaf)

df005 <- Metadata_eigenvec[30:34]
pca_005 <- prcomp(df005, scale. = TRUE)
autoplot(pca_005)
summary(pca_005)

df01 <- Metadata_eigenvec[35:39]
pca_01 <- prcomp(df01, scale. = TRUE)
autoplot(pca_01)
summary(pca_01)


autoplot(pca_005, data = Metadata_eigenvec, colour = 'CCR')
autoplot(pca_005, data = Metadata_eigenvec, colour = 'Group_S')
autoplot(pca_005, data = Metadata_eigenvec, colour = 'Growing_Cycle')
autoplot(pca_005, data = Metadata_eigenvec, colour = 'Crop')
autoplot(pca_005, data = Metadata_eigenvec, colour = 'Greenhouse_simplified')
autoplot(pca_005, data = Metadata_eigenvec, colour = 'Region')
autoplot(pca_005, data = Metadata_eigenvec, colour = 'Fenhexamid')
autoplot(pca_005, data=Metadata_eigenvec, colour='Group_S', frame = TRUE, frame.type = 'norm')
autoplot(pca_005, data=Metadata_eigenvec, colour='CCR', frame = TRUE)
autoplot(pca_01, data = Metadata_eigenvec, colour = 'CCR')
autoplot(pca_01, data = Metadata_eigenvec, colour = 'Group_S')
autoplot(pca_01, data = Metadata_eigenvec, colour = 'Growing_Cycle')
autoplot(pca_01, data = Metadata_eigenvec, colour = 'Crop')
autoplot(pca_01, data = Metadata_eigenvec, colour = 'Greenhouse_simplified')
autoplot(pca_01, data = Metadata_eigenvec, colour = 'Region')
autoplot(pca_01, data = Metadata_eigenvec, colour = 'Fenhexamid')
autoplot(pca_nomaf, data=Metadata_eigenvec, colour='Group_S', frame = TRUE, frame.type = 'norm')
autoplot(pca_nomaf, data=Metadata_eigenvec, colour="Group_S", frame = TRUE)






#Using ggplot may be easier to modify

#PC1 = 20.20%; PC2= 20.14%

library(pals)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)

library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()

# Standard with sample names

#consider removing BU10 in nomaf pca's (outlier)
nomaf_BU10rm <- Metadata_eigenvec[Metadata_eigenvec$Sample != "BU10", ]

p_na_BU10rm <- ggplot(nomaf_BU10rm,aes(x=PC3_nomaf,y=PC2_nomaf,col=Mating_type))+
  geom_text(aes(label = Sample)) +
  scale_color_brewer(palette="Dark2")+ 
  theme_classic()
p_na_BU10rm

p_01 <- ggplot(Metadata_eigenvec,aes(x=PC3_maf01,y=PC4_maf01,col=Mating_type))+
  geom_text(aes(label = Sample)) +
  scale_color_brewer(palette="Dark2")+ 
  theme_classic()
p_01


p_005 <- ggplot(Metadata_eigenvec,aes(x=PC1_maf005,y=PC2_maf005,col=Mating_type))+
  geom_text(aes(label = Sample)) +
  scale_color_brewer(palette="Dark2")+ 
  theme_classic()
p_005


cd 


# Group S with points or sample names
ggplot(Metadata_eigenvec,aes(x=PC1_maf005,y=PC2_maf005,col=Group_S))+
  geom_point(size=3,alpha=0.5)+ 
  scale_color_brewer(palette="Dark2")+ 
  theme_classic()
ggplot(Metadata_eigenvec,aes(x=PC1_maf005,y=PC2_maf005,col=Group_S))+
  geom_text(aes(label = Sample)) +
  scale_color_brewer(palette="Dark2")+ 
  theme_classic()

# Change color scheme
ggplot(Metadata_eigenvec,aes(x=PC1_maf005,y=PC2_maf005,col=Greenhouse_simplified))+
  geom_point(size=3,alpha=0.8)+ 
  scale_fill_manual(values=as.vector(watlington(12)))+
  #scale_color_brewer(palette="RdYlGn")+ 
  theme_classic()






p_framed <- p_005 + stat_ellipse(geom="polygon", aes(fill = Mating_type),
                             alpha = 0.2, show.legend = FALSE, level = 0.975) +
  xlab("PC 1") + 
  ylab("PC 2") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent"))
p_framed




------------------
  
# ------------ Population structure ---------------

Metadata_K4 <- data.frame(Metadata_eigenvec)
Metadata_K4[Metadata_K4==""] <- NA

Group_S_labels <- c("sensu stricto", "Group S")
names(Group_S_labels) <- c("0", "1")

ggsave("Structure_GroupS.tiff" , Structure_GroupS, device='tiff', dpi=400)


#facet wrap vs facet grid

Structure_Matingtype <- Metadata_K4 %>%
  #arrange(Cluster_3, desc(Value)) %>%
  mutate(Sample = fct_reorder(Sample, K4_C3, .desc=T)) %>%
  gather(key=Cluster, value=Value, c(K4_C1, K4_C2, K4_C3, K4_C4)) %>%
  ggplot(aes(x=Sample, y=Value, fill=Cluster)) + 
  facet_wrap(~Mating_type, scales="free", nrow=3)+
  geom_col() +
  theme_light() +
  theme(strip.background = element_rect(fill="gray86", linewidth =1, linetype = "solid"),
        text = element_text(size = 10, color = "black"), 
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 10, color = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_fill_manual(values = c("#669966", "#666666", "#006699", "#CC9900"))


Structure_Matingtype
ggsave("Structure_Matingtype.tiff" , Structure_Matingtype, device='tiff', dpi=400)



Structure_GroupS <- Metadata_K4 %>%
  #arrange(Cluster_3, desc(Value)) %>%
  mutate(Sample = fct_reorder(Sample, K4_C1, .desc=T)) %>%
  gather(key=Cluster, value=Value, c(K4_C1, K4_C2, K4_C3, K4_C4)) %>%
  ggplot(aes(x=Sample, y=Value, fill=Cluster)) + 
  facet_wrap(~Group_S, scales="free", nrow=2, labeller=labeller(Group_S= Group_S_labels))+
  geom_col() +
  theme_light() +
  theme(strip.background = element_rect(fill="gray86", linewidth =1, linetype = "solid"),
        text = element_text(size = 11, color = "black"), 
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 11, color = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_fill_manual(values = c("#669966", "#666666", "#006699", "#CC9900"), labels= c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))
Structure_GroupS
ggsave("StructureK4_GroupS_C2.tiff" , Structure_GroupS, device='tiff', dpi=400)


Metadata_K4_phenotyped <- subset(Metadata_K4, subset=Phenotyped %in% c("1"))
Structure_CCR <- Metadata_K4_phenotyped %>%
  #arrange(Cluster_3, desc(Value)) %>%
  mutate(Sample = fct_reorder(Sample, K4_C3, .desc=T)) %>%
  gather(key=Cluster, value=Value, c(K4_C1, K4_C2, K4_C3, K4_C4)) %>%
  ggplot(aes(x=Sample, y=Value, fill=Cluster)) + 
  facet_wrap(~CCR, scales="free", nrow=7)+
  geom_col() +
  theme_light() +
  theme(strip.background = element_rect(fill="gray86", linewidth =1, linetype = "solid"),
        text = element_text(size = 11, color = "black"), 
        axis.title.y = element_blank(),
        #axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 11, color = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_fill_manual(values = c("#669966", "#666666", "#006699", "#CC9900"), labels= c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))
Structure_CCR
ggsave("StructureK4_CCR_C1.tiff" , Structure_CCR, device='tiff', dpi=400)




#look at log4.out to see Fst values between clusters
#look at pong results in browser for propotions of each cluster





#---------- Genetic Distance Tree ---------
  
#  Genetic distance tree, statistics, and more

# https://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html

library(vcfR)
library(poppr)
library(ape)


setwd("C:/Users/nikki/Michigan State University/PSM.Hausbecklab - Nikki Lukasko - Nikki Lukasko/Botrytis/Molecular/Plates123 Population Genomics/Population Structure Analysis")


# read in the vcf file ------------UPDATE TO VCF CONTAINING 15k SNPs INSTEAD OF 32k---------------DONE---

vcf_Bcin <- read.vcfR("SNV_plinkMAF05_FINAL_v2.vcf.gz")
head(vcf_Bcin)
vcf_Bcin@meta
write.vcf(vcf_Bcin, "SNV_plinkMAF05_FINAL_v2.vcf.gz")
colnames(vcf_Bcin@gt)

#convert file type to a genlight object for downstream compatibility
vcf_Bcin_genlight <- vcfR2genlight(vcf_Bcin)

#check a few genotypes from the genlight object
gt <- extract.gt(vcf_Bcin, element = "GT")
gt[c(2,6,18), 1:3]

ploidy(vcf_Bcin_genlight) <- 1

#Add population data (can change according to interest)

pop(vcf_Bcin_genlight) <- Metadata_K4$Crop
vcf_Bcin_genlight


myDiff <- genetic_diff(vcf_Bcin, pops = pop(vcf_Bcin_genlight), method = 'nei')
knitr::kable(head(myDiff[,1:13]))
knitr::kable(round(colMeans(myDiff[,c(3:11)], na.rm = TRUE), digits = 3))
dpf <- melt(myDiff[,c(3:6,10,11)], varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
violin <- ggplot(dpf, aes(x=variable, y=Depth)) + geom_violin(fill="#2ca25f", adjust = 1.2)
violin <- violin + xlab("")
violin <- violin + ylab("")
violin <- violin + theme_bw()
violin

# Plot Crop (or other 'pop' variable) by genomic position (Manhattan plot)
Crop_Gst_manhat <- plot(getPOS(vcf_Bcin), myDiff$Gprimest,  pch = 20, col = "#1E90FF44", xlab = "", ylab = "", ylim = c(0, 1), xaxt = "n") +
axis(side = 1, at = seq(0, 43000000, by = 430), labels = seq(0, 43000000, by = 430)) +
title(xlab='Genomic position (Kbp)') +
title(ylab = expression(italic("G'"["ST"])))



# Subset example if wanted: https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
vcf_Bcin_CCR <- vcf_Bcin[,c(TRUE, Metadata_K4$CCR == c("0", "1", "2", "3", "4", "5", "6", "7"))] #this would remove CCR=0 or 7)
vcf_Bcin_CCR #doesnt seem right, check later.

#Create distance tree
?poppr::bitwise.dist
?dist
#may want to consider other distances matrices?
#vcf_Bcin_genlight.dist <- dist(vcf_Bcin_genlight)
vcf_Bcin_genlight.dist <- poppr::bitwise.dist(vcf_Bcin_genlight, scale_missing = TRUE)

# change bootstrap (sample) to 10 instead of 100 for quicker viewing
tree <- aboot(vcf_Bcin_genlight, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = TRUE, cutoff = 70, quiet = F)
?aboot


#-------- run tree with 100 bootstraps, need to rerun and save stuff below ----------


pop(vcf_Bcin_genlight) <- Metadata_K4$Mating_type
pop(vcf_Bcin_genlight) <- Metadata_K4$CCR
pop(vcf_Bcin_genlight) <- Metadata_K4$Group_S


cols <- brewer.pal(n = nPop(vcf_Bcin_genlight), name = "Set3")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(vcf_Bcin_genlight)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("Group S", "sensu stricto"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("sensu stricto", "Group S"), fill = cols, border = FALSE, bty = "n", cex = 0.8)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

# Saving tree as a tiff

pop(vcf_Bcin_genlight) <- Metadata_K4$Group_S
tiff('GeneticDistanceTree_GroupS.tiff', units="in", width=4, height=35, res=300, compression = 'lzw')
plot.phylo(tree, cex = 0.8, font = 1, adj = 0, tip.color =  c("blue", "red")[pop(vcf_Bcin_genlight)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 1, xpd = TRUE)
#legend(35,10,c("Group S", "sensu stricto"),cols, border = FALSE, bty = "n")
#legend('topleft', legend = c("sensu stricto", "Group S"), fill = cols, border = FALSE, bty = "n", cex = 0.8)
axis(side = 1)
title(xlab = "Genetic distance \n(proportion of loci that are different)")
dev.off()

pop(vcf_Bcin_genlight) <- Metadata_K4$CCR
tiff('GeneticDistanceTree_CCR.tiff', units="in", width=4, height=35, res=300, compression = 'lzw')
plot.phylo(tree, cex = 0.8, font = 1, adj = 0, tip.color =  c("darkblue", "blue", "blue", "blue", "green", "green", "red", "red")[pop(vcf_Bcin_genlight)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 1, xpd = TRUE)
#legend(35,10,c("Group S", "sensu stricto"),cols, border = FALSE, bty = "n")
#legend('topleft', legend = c("sensu stricto", "Group S"), fill = cols, border = FALSE, bty = "n", cex = 0.8)
axis(side = 1)
title(xlab = "Genetic distance \n(proportion of loci that are different)")
dev.off()


#-----trying---
  
pop(vcf_Bcin_genlight) <- Metadata_K4$CCR
tiff('GeneticDistanceTree_test.tiff', units="in", width=4, height=35, res=300, compression = 'lzw')
plot.phylo(tree, cex = 0.8, font = 1, adj = 0, tip.color =  c("blue", "red")[pop(vcf_Bcin_genlight)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 1, xpd = TRUE)
pop(vcf_Bcin_genlight) <- Metadata_K4$Crop
tiplabels(pch=19, frame="circle", adj=-5, bg=c("blue", "red", "green")[pop(vcf_Bcin_genlight)])
dev.off()


plot.phylo(tree, cex = 0.8, font = 1, adj = 0, tip.color =  c("blue", "red")[pop(vcf_Bcin_genlight)])
plot(tree)
phydataplot(Fungicides, tree, style = "bars", offset=3)



Metadata_K4 %>% gather(key=Fungicide, value=Value, c(Fenhexamid, Iprodione, Fludioxonil, 
                                     Thiophanate.methyl, Boscalid, Fluopyram,
                                     Cyprodinil, Pyraclostrobin))


#legend(35,10,c("Group S", "sensu stricto"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("sensu stricto", "Group S"), fill = cols, border = FALSE, bty = "n", cex = 0.8)
axis(side = 1)
title(xlab = "Genetic distance \n(proportion of loci that are different)")
plot(bird.orders, "c", FALSE, font = 1, label.offset = 3,
     x.lim = 31, no.margin = TRUE)
tiplabels(pch = 1, bg = gray(1:23/23), cex = 2, adj = 1.4)
tiplabels(pch = 19, col = c("yellow", "red", "blue"), adj = 2.5, cex = 2)

-----



require(phytools)
trait_values <- Metadata_K4$Crop[match(tree$tip.label, Metadata_K4$CCR)]

# Create a color palette for the 'Crop' trait (adjust the number of unique values if needed)
colors <- rainbow(length(unique(trait_values)))

# Assign colors to trait values
tip_colors <- colors[as.factor(trait_values)]

# Plot the tree with tip labels colored by the 'Crop' trait
plot(tree, tip.color = tip_colors, show.tip.label = TRUE, cex = 0.6)
legend("topright", legend = levels(factor(trait_values)), fill = colors, title = "Crop")



# Extract the CCR column from genlight_popdata
trait_values <- Metadata_K4$CCR[match(tree$tip.label, Metadata_K4$Sample)]

# Define colors based on CCR values
tip_colors <- rep(" ", length(trait_values))  # Set label to empty for NA values

# Assign colors based on CCR values
tip_colors[trait_values == 0] <- "blue"  # CCR = 0, set color to blue
tip_colors[trait_values %in% c(1, 2)] <- "blue"  # CCR = 1 or 2, set color to blue
tip_colors[trait_values %in% c(3, 4, 5, 6)] <- "green"  # CCR = 3, 4, 5, or 6, set color to green
tip_colors[trait_values == 7] <- "red"  # CCR = 7, set color to red

# Plot the tree with tip labels colored by the 'CCR' trait
plot(tree, tip.color = tip_colors, show.tip.label = TRUE, cex = 0.6)



# Distance tree colored by CCR (leaving NA as blank) ---------------

tree2 <- aboot(vcf_Bcin_genlight_CCR, tree = "upgma", distance = bitwise.dist, sample = 10, showtree = TRUE, cutoff = 50, quiet = F)
tree2

# Create a vector of tip labels for the filtered data
tip_labels <- genlight_popdata_CCR$Sample

# Create a vector of trait values for the filtered data
trait_values <- as.numeric(genlight_popdata_CCR$CCR)

# Define colors for different CCR values
colors <- c("blue", "purple", "green", "red")

# Assign colors based on CCR values
tip_colors <- colors[cut(trait_values, breaks = c(-1, 0.5, 2.5, 6.5, 8), labels = FALSE)]

# Ensure row names of 'tree' match the values in the 'Sample' column of filtered_data
tree2$tip.label <- tip_labels[match(tree2$tip.label, tip_labels)]

# Plot the modified tree with colored tips
plot(tree2, tip.color = tip_colors, show.tip.label = TRUE, cex = 0.6)






library(phytools)
dotTree(tree,genlight_popdata,colors=cols[pop(vcf_Bcin_genlight)])
pop(vcf_Bcin_genlight) <- genlight_popdata$CCR
color.plot.phylo(tree, genlight_popdata, trait = genlight_popdata$Crop)



# -------------Discriminant analysis of principal components (DAPC)----------------
?dapc
Bcin.dapc <- dapc(vcf_Bcin_genlight, n.pca = NULL, n.da = NULL) # chose 100 and 5

Bcin.dapc1 <- dapc(vcf_Bcin_genlight, n.pca = 20, n.da = 3, scale = TRUE)
Bcin.dapc2 <- dapc(vcf_Bcin_genlight, n.pca = 20, n.da = 3, scale = FALSE)

Bcin.dapc_100 <- dapc(vcf_Bcin_genlight, n.pca = 100, n.da = 5)

scatter(Bcin.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topright", cleg = 0.75)




## ------------RUN AMOVA ------------------

package?poppr


Metadata_K4
strata(vcf_Bcin_genlight) <- Metadata_K4
table(strata(vcf_Bcin_genlight, ~CCR))

# without clone correction

amova.result_GH <- poppr.amova(vcf_Bcin_genlight, ~Greenhouse)
amova.result_GH
amova.test_GH <- randtest(amova.result_GH, nrepet = 999) #significance testing
amova.test_GH
plot(amova.test_GH)

#individual greenhouse visits
amova.result_code <- poppr.amova(vcf_Bcin_genlight, ~Code)
amova.result_code
amova.test_code <- randtest(amova.result_code, nrepet = 999) #significance testing
amova.test_code
plot(amova.test_code)

amova.result_Crop <- poppr.amova(vcf_Bcin_genlight, ~Crop)
amova.result_Crop
amova.test_Crop <- randtest(amova.result_Crop, nrepet = 999) #significance testing
amova.test_Crop
plot(amova.test_Crop)

amova.result_cycle <- poppr.amova(vcf_Bcin_genlight, ~Growing_Cycle)
amova.result_cycle
amova.test_cycle <- randtest(amova.result_cycle, nrepet = 999) #significance testing
amova.test_cycle
plot(amova.test_cycle)

amova.result_groupS <- poppr.amova(vcf_Bcin_genlight, ~Group_S)
amova.result_groupS
amova.test_groupS <- randtest(amova.result_groupS, nrepet = 999) #significance testing
amova.test_groupS
plot(amova.test_groupS)



# Run AMOVA for each variable
amova_result_greenhouse <- poppr.amova(vcf_Bcin_genlight, ~Greenhouse)
amova_result_cycle <- poppr.amova(vcf_Bcin_genlight, ~Growing_Cycle)
amova_result_group <- poppr.amova(vcf_Bcin_genlight, ~Group_S)
amova_result_code <- poppr.amova(vcf_Bcin_genlight, ~Code)

# Perform significance testing for each variable
amova_test_greenhouse <- randtest(amova_result_greenhouse, nrepet = 999)
amova_test_cycle <- randtest(amova_result_cycle, nrepet = 999)
amova_test_group <- randtest(amova_result_group, nrepet = 999)
amova_test_code <- randtest(amova_result_code, nrepet = 999)

# Compile results into a list (excluded p-values from above)
amova_results_list <- list(
  "Greenhouse", amova_result_greenhouse,
  "Growing Cycle", amova_result_cycle,
  "Group S", amova_result_group,
  "Individual Greenhouse Visit", amova_result_code
)
#above is missing crop so added manually to manuscript data

amova_results_list

#NEED TO SUBSET TO REMOVE "BOTH" AS A MATING TYPE POPULATION
pop(vcf_Bcin_genlight) <- Metadata_K4$Mating_type
vcf_Bcin_genlight_mat <- popsub(vcf_Bcin_genlight, exclude = "Both") #Removes all samples without fungicide phenotype data
Metadata_K4_mat <- Metadata_K4[!grepl("Both", Metadata_K4$Mating_type), ]
popNames(vcf_Bcin_genlight_mat)
strata(vcf_Bcin_genlight_mat) <- Metadata_K4_mat


amova.result_mat <- poppr.amova(vcf_Bcin_genlight_mat, ~Mating_type)
amova.result_mat
amova.test_mat <- randtest(amova.result_mat, nrepet = 999) #significance testing
amova.test_mat
plot(amova.test_mat)


# AMOVA on subset of phenotyped data (Can run on CCR or each fungicide individually)
#make sure pop is CCR before continuing
vcf_Bcin_genlight_CCR <- popsub(vcf_Bcin_genlight, exclude = "NA") #Removes all samples without fungicide phenotype data
popNames(vcf_Bcin_genlight_CCR)
pop(vcf_Bcin_genlight_CCR) <- Metadata_K4_phenotyped$CCR
Metadata_K4_phenotyped <- subset(Metadata_K4, subset=Phenotyped %in% c("1"))

strata(vcf_Bcin_genlight_CCR) <- Metadata_K4_phenotyped
amova.result_CCR <- poppr.amova(vcf_Bcin_genlight_CCR, ~CCR)
amova.result_CCR
amova.test_CCR <- randtest(amova.result_CCR, nrepet = 999)
amova.test_CCR
plot(amova.test_CCR)

# Rerun for each fungicide
amova.result_fungicide <- poppr.amova(vcf_Bcin_genlight_CCR, ~Thiophanate.methyl)
amova.result_fungicide
amova.test_fungicide <- randtest(amova.result_fungicide, nrepet = 999)
amova.test_fungicide
plot(amova.test_fungicide)

# List of fungicides
fungicides <- c("Thiophanate.methyl", "Fenhexamid", "Iprodione", "Boscalid", "Fluopyram", "Pyraclostrobin", "Cyprodinil", "Fludioxonil")

# Creating an empty list to store results
results_list <- list()

# Loop through each fungicide
for (fungicide in fungicides) {
  # Run the commands for each fungicide
  amova_result <- poppr.amova(vcf_Bcin_genlight_CCR, as.formula(paste("~", fungicide)))
  print(amova_result)
  
  amova_test <- randtest(amova_result, nrepet = 999)
  
  # Store the results in the list
  results_list[[fungicide]] <- list(amova_result = amova_result, amova_test = amova_test)
}

# Access the results for a specific fungicide, for example, "Iprodione"
results_list[["Iprodione"]]











# Could group CCR 0-1, 2-5, 6-7 or otherwise for AMOVA (just need to change names in excel sheet and reload strata)
#vcf_Bcin_genlight_CCR0167 <- popsub(vcf_Bcin_genlight, exclude = c("2", "3", "4", "5")) #Removes all samples without fungicide phenotype data
#the following is grouped into "0-2" and "3-7"
#genlight_popdata_CCR_GROUPED <- read_excel("genlight_popdata_CCR_GROUPED.xlsx", 
#                                           col_types = c("text", "text", "text", 
#                                                         "numeric", "text", "text", "text", 
#                                                         "text", "text", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "text"))
#strata(vcf_Bcin_genlight_CCR) <- genlight_popdata_CCR_GROUPED
#pop(vcf_Bcin_genlight_CCR) <- genlight_popdata_CCR_GROUPED$CCR
#popNames(vcf_Bcin_genlight_CCR)
#amova.result_CCR_group <- poppr.amova(vcf_Bcin_genlight_CCR, ~CCR)
#amova.result_CCR_group
#amova.test_CCR_group <- randtest(amova.result_CCR_group, nrepet = 999)
#amova.test_CCR_group
#plot(amova.test_CCR_group)



#--- Should I clone correct before running AMOVA? B. cinerea reproduces sexually, so probably not. ---

#clone correcting does not change-I believe there are no clones.
amova.result_GroupScc <- poppr.amova(vcf_Bcin_genlight, ~Group_S, clonecorrect = TRUE)
amova.result_GroupScc #clone-corrected
amova.result_groupS #not clone-corrected

# Unused
amova.result_CrGH <- poppr.amova(vcf_Bcin_genlight, ~Growing_Cycle/Greenhouse)
amova.result_CrGH
amova.test_CrGH <- randtest(amova.result_CrGH, nrepet = 999)
amova.test_CrGH



# https://grunwaldlab.github.io/Population_Genetics_in_R/clustering_plot.html


#This will take a while to run- load the environment (open in R studio) from the same folder to skip re-running
#check correct SNP subset first
library(adegenet)
maxK <- 5
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(vcf_Bcin_genlight, n.pca = 50, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1



my_k <- 2:5
# ------------ TRY CHANGE MY_K TO 2:3 INSTEAD BECAUSE IT APPEARS THAT IS MORE ACCURATE? ------------

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(133)
  grp_l[[i]] <- find.clusters(vcf_Bcin_genlight, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(vcf_Bcin_genlight, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(vcf_Bcin_genlight, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

# Write the dataframe to an Excel file so I can change colors according to CCR in DAPC graph

write.csv(my_df, "DAPC_coord.csv")




Comparison_for_DAPC_Groups <- read_excel("Comparison for DAPC Groups.xlsx", 
                                         col_types = c("text", "numeric", "numeric", 
                                                       "numeric", "numeric", "text", "text", 
                                                       "text", "text", "numeric", "numeric", 
                                                       "text", "text", "numeric", "text", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "text"))


# CCR only
scatter_CCR <- ggplot(Comparison_for_DAPC_Groups, aes(x = LD1, y = LD2, color = DAPC_Group, size = CCR)) +
  geom_point(alpha=0.4) +
  labs(x = "LD1", y = "LD2", color = "DAPC Group", size = "CCR") +
  theme_minimal()

pairs(Comparison_for_DAPC_Groups[c("LD1", "LD2", "LD3", "LD4")], pch = 16, col = Comparison_for_DAPC_Groups$DAPC_Group)



# CCR and Group S
ggplot(Comparison_for_DAPC_Groups, aes(x = LD2, y = LD4, color = DAPC_Group, size = CCR, shape = as.factor(GroupS))) +
  geom_point(alpha = 0.5) +
  labs(x = "LD2", y = "LD4", color = "DAPC Group", size = "CCR") +
  scale_shape_manual(values = c(16, 1), name = "GroupS", labels = c("0", "1")) +
  theme_minimal()

# Group S only
scatter_GroupS <- ggplot(Comparison_for_DAPC_Groups, aes(x = LD1, y = LD2, color = DAPC_Group, shape = as.factor(GroupS))) +
  geom_point(alpha=0.5, size=5) +
  labs(x = "LD1", y = "LD2", color = "DAPC Group", shape = "Group S") +
  scale_shape_manual(values = c(16, 12), name = "GroupS", labels = c("0", "1")) +
  theme_minimal()
ggplot(Comparison_for_DAPC_Groups, aes(x = LD4, y = LD1, color = DAPC_Group, size = GroupS)) +
  geom_point(alpha=0.4) +
  labs(x = "LD4", y = "LD1", color = "DAPC Group", size = "GroupS") +
  theme_minimal()


DAPC_frame <- data.frame(Comparison_for_DAPC_Groups)
DAPC_frame_phenotyped <- subset(DAPC_frame, subset = Phenotyped == 1)

# Table with Counts of each CCR within each DAPC Group
dapc_ccr_table <- table(DAPC_frame_phenotyped$CCR, DAPC_frame_phenotyped$DAPC_Group)
print(dapc_ccr_table)



# Create a stacked bar chart
barchart_CCR <- ggplot(data = DAPC_frame_phenotyped, aes(x = DAPC_Group, fill = as.factor(CCR))) +
  geom_bar() +
  labs(x = "DAPC Group", y = "Count") +
  scale_fill_viridis_d() +
  theme_minimal()

barchart_GroupS <- ggplot(data = DAPC_frame_phenotyped, aes(x = DAPC_Group, fill = as.factor(GroupS))) +
  geom_bar() +
  labs(x = "DAPC Group", y = "Count") +
  scale_fill_viridis_d() +
  theme_minimal()



# Combining for output

# Load required libraries
library(patchwork)
library(viridisLite)


# Combine plots using patchwork
combined_plots_GroupS <- scatter_GroupS / barchart_GroupS
combined_plots_GroupS



# write.table(my_df,file='DAPC_groups.csv')

#--------- View my_df_ and pull out group 2, which separated horizontally in p2 graph -----------

my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

p2 <- ggplot(my_df, aes(x = LD3, y = LD4, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2


# --------TRY AGAIN WITH LD2 AND LD3 --------------------------







# Calculating nucleotide diversity ----- https://www.york.ac.uk/res/dasmahapatra/teaching/MBiol_sequence_analysis/workshop4_2019.html#plotting_nucleotide_diversity

pi.all <- read.table("NucDiv_10kb.windowed.pi",header=T)
head(pi.all)
hist(pi.all$PI,br=20)

boxplot(pi.all$PI,ylab="diversity")

pi.chr5 <- subset(pi.all, CHROM == "5")
boxplot(pi.chr5$PI,ylab="diversity")

#find out how many rows your objects have
nrow(pi.chr5)
nrow(pi.all)

#see thr first and last lines
head(pi.chr5)
head(pi.all)

tail(pi.chr5)
tail(pi.all)

#see a summary of the data frame
summary(pi.chr5)
summary(pi.all)

plot(pi.chr5$BIN_START,pi.chr5$PI,xlab="position",ylab="diversity")
plot(pi.all$BIN_START,pi.all$PI,xlab="position",ylab="diversity") #not sure this shows all chromosomes or maybe overlaps them?


# Graphing Tajima's D values
# using VCF that originated from the hwe filtered file in newindeppairwise folder, before taking LD subset

taj.hwe <- read.table("Tajima_hwe_10kb.Tajima.D",header=T) 
taj.ld <- read.table("Tajima_all_samples_10kb.Tajima.D",header=T) 

# note- this was run on LD subset so its bias - try on entire vcf
hist(taj.hwe$TajimaD,br=20)
hist(taj.ld$TajimaD,br=20)


# Create a list of datasets
data <- list("Before LD subset" = taj.hwe$TajimaD, "After LD subset" = taj.ld$TajimaD)

# Create side-by-side boxplots
boxplot(data, col = c("skyblue", "salmon"),
        main = "Tajima's D Comparison",
        names = c("Before LD subset (570,854)", "After LD subset (15,207)"),
        xlab = "Dataset (SNVs)", ylab = "Tajima's D")


# Comparison of Tajima's D before and after LD filtering on chromosome 5

taj.chr5_hwe <- subset(taj.hwe, CHROM == "5")
taj.chr5_ld <- subset(taj.ld, CHROM == "5")

data2 <- list("Before LD subset" = taj.chr5_hwe$TajimaD, "After LD subset" = taj.chr5_ld$TajimaD)
boxplot(data2, col = c("skyblue", "salmon"),
        main = "Tajima's D Comparison (Chr 5 only)",
        names = c("Before LD subset (570,854)", "After LD subset (15,207)"),
        xlab = "Dataset (SNVs)", ylab = "Tajima's D")


plot(taj.chr5_hwe$BIN_START,taj.chr5_hwe$TajimaD,xlab="position",ylab="Tajima's D")
plot(taj.chr5_ld$BIN_START,taj.chr5_ld$TajimaD,xlab="position",ylab="Tajima's D")

#----------- per site Tajima's D -------------


# read in the file
head(taj.ld)
# capdown the headers
names(taj.ld) <- tolower(taj.ld)
names(taj.ld)[grep("taj.ld")] <- "taj.ld"

TajimaD_chr1 <- subset(taj.ld, subset=CHROM %in% c("1"))
TajimaD_chr2 <- subset(taj.ld, subset=CHROM %in% c("2"))
TajimaD_chr3 <- subset(taj.ld, subset=CHROM %in% c("3"))
TajimaD_chr4 <- subset(taj.ld, subset=CHROM %in% c("4"))
TajimaD_chr5 <- subset(taj.ld, subset=CHROM %in% c("5"))
TajimaD_chr6 <- subset(taj.ld, subset=CHROM %in% c("6"))
TajimaD_chr7 <- subset(taj.ld, subset=CHROM %in% c("7"))
TajimaD_chr8 <- subset(taj.ld, subset=CHROM %in% c("8"))
TajimaD_chr9 <- subset(taj.ld, subset=CHROM %in% c("9"))
TajimaD_chr10 <- subset(taj.ld, subset=CHROM %in% c("10"))
TajimaD_chr11 <- subset(taj.ld, subset=CHROM %in% c("11"))
TajimaD_chr12 <- subset(taj.ld, subset=CHROM %in% c("12"))
TajimaD_chr13 <- subset(taj.ld, subset=CHROM %in% c("13"))
TajimaD_chr14 <- subset(taj.ld, subset=CHROM %in% c("14"))
TajimaD_chr15 <- subset(taj.ld, subset=CHROM %in% c("15"))
TajimaD_chr16 <- subset(taj.ld, subset=CHROM %in% c("16"))

# CHROM 5
ggplot(TajimaD_chr5, aes(BIN_START, TajimaD)) + geom_point()
quantile(TajimaD_chr5$TajimaD, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(TajimaD_chr5$TajimaD, 0.975, na.rm = T)
# make an outlier column in the data.frame
TajimaD_chr5_outliers <- TajimaD_chr5 %>% mutate(outlier = ifelse(TajimaD > my_threshold, "outlier", "background"))
TajimaD_chr5_outliers %>% tally()
TajimaD_chr5_outliers %>% group_by(outlier) %>% tally()
TajimaD_plot_chr5 <- ggplot(TajimaD_chr5_outliers, aes(BIN_START, TajimaD, colour = outlier)) + geom_point()
TajimaD_plot_chr5

# All

ggplot(taj.ld, aes(BIN_START, TajimaD)) + geom_point()
quantile(taj.ld$TajimaD, c(0.975, 0.995), na.rm = T)
quantile(taj.ld$TajimaD, c(0.025, 0.005), na.rm = T)

# identify the 95% percentile
my_threshold <- quantile(taj.ld$TajimaD, 0.975, na.rm = T)
my_threshold_low <- quantile(taj.ld$TajimaD, 0.005, na.rm = T)

# make an outlier column in the data.frame
TajimaD_all_outliers <- taj.ld %>% mutate(outlier = ifelse(TajimaD > my_threshold, "outlier", "background"))
TajimaD_all_outliers %>% tally()
TajimaD_all_outliers %>% group_by(outlier) %>% tally()
taj.ld_highoutlier <- ggplot(TajimaD_all_outliers, aes(BIN_START, TajimaD, colour = outlier)) + geom_point()
taj.ld_highoutlier

outliers_all <- subset(TajimaD_all_outliers, subset=outlier %in% c("outlier"))
outliers_all

# make an outlier column in the data.frame
TajimaD_all_outliers_low <- taj.ld %>% mutate(outlier = ifelse(TajimaD < my_threshold_low, "outlier", "background"))
TajimaD_all_outliers_low %>% tally()
TajimaD_all_outliers_low %>% group_by(outlier) %>% tally()
taj.ld_lowoutliers <- ggplot(TajimaD_all_outliers_low, aes(BIN_START, TajimaD, colour = outlier)) + geom_point()
taj.ld_lowoutliers

outliers_low <- subset(TajimaD_all_outliers_low, subset=outlier %in% c("outlier"))
outliers_low







#-------Relationship between cluster membership and fungicide resistance ------------


setwd("C:/Users/nikki/Michigan State University/PSM.Hausbecklab - Nikki Lukasko - Nikki Lukasko/Botrytis/Molecular/Plates123 Population Genomics/Population Structure Analysis"

library(readxl)
library(dplyr)
library(tidyr)


membership_fungicide <- read_excel("membership_fungicide.xlsx", 
                                          col_types = c("text", "numeric", "numeric", 
                                                                  "numeric", "numeric", "numeric", 
                                                                  "numeric", "numeric", "numeric", 
                                                                  "numeric", "numeric", "numeric", 
                                                                  "numeric"))
View(membership_fungicide)
membership_fungicide <- as.data.frame(membership_fungicide)

