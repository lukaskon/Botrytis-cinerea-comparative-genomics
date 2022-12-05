#Botrytis cinerea Population Structure Analysis 
#Admixture v1.3
#10/31/22
#Author: Nikki Lukasko

library(ggplot2)
library(readxl)
library(tidyverse)
library(forcats)
library(dplyr)


###NOTE: Samples need to be reordered before analysis



setwd("C:/Users/nikki/Michigan State University/PSM.Hausbecklab - Nikki Lukasko - Nikki Lukasko/Botrytis/Molecular/Population Structure Analysis")



#####Including Indels


tbl=read.table("BcinSNPIND_ann.3.Q")
view(tbl)
barplot(t(as.matrix(tbl)), col=rainbow(3), xlab="Individual #", ylab="Ancestry", border=NA)

K3_pop <- read_excel("Snp_Indel_K3_pop_admixture.xlsx", 
                                                col_types = c("text", "text", "text", 
                                                                       "text", "text", "text", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric"))

K3_framed <- data.frame(K3_pop)
#option 1
K3_framed %>%
  #group_by() %>%
  #summarise_at(vars(Cluster_1, Cluster_2, Cluster_3)) %>%
  gather(key=Cluster, value=Value, c(Cluster_1, Cluster_2, Cluster_3)) %>%
  ggplot(aes(x=Sample, y=Value, fill=Cluster)) + facet_wrap(~Region, scales="free")+
  geom_col() +
  theme(text = element_text(size = 10), 
        axis.text.x.bottom = element_text(size=10, angle=45, vjust = 0.85, hjust = 0.8),
        axis.title.y = element_text(vjust = 2.5))

#option 2
K3_framed %>%
  #group_by() %>%
  #summarise_at(vars(Cluster_1, Cluster_2, Cluster_3)) %>%
  gather(key=Cluster, value=Value, c(Cluster_1, Cluster_2, Cluster_3)) %>%
  ggplot(aes(x=Sample, y=Value, fill=Cluster)) + facet_grid(~Greenhouse, scales="free_x", space ="free_x", switch="x")+
  geom_col() +
  theme()




#####SNPs only



K4_pop <- read_excel("Snp_K4_pop_admixture.xlsx", 
                     col_types = c("text", "text", "text", 
                                   "text", "text", "text", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", "numeric"))

K4_framed <- data.frame(K4_pop)





#option 1

#lvls <- names(sort(table(K4_framed == "Cluster_3", "Sample")))





K4_framed %>%
  #arrange(Cluster_3, desc(Value)) %>%
  mutate(Sample = fct_reorder(Sample, Cluster_4, .desc=T)) %>%
  gather(key=Cluster, value=Value, c(Cluster_1, Cluster_2, Cluster_3, Cluster_4)) %>%
  ggplot(aes(x=Sample, y=Value, fill=Cluster)) + facet_wrap(~CCR, scales="free")+
  geom_col() +
  theme(text = element_text(size = 10), 
        axis.text.x.bottom = element_text(size=10, angle=45, vjust = 0.85, hjust = 0.8),
        axis.title.y = element_text(vjust = 2.5))


K4_framed %>%
  mutate(Sample = fct_reorder(Sample, Cluster_4)) %>%
  gather(key=Cluster, value=Value, c(Cluster_1, Cluster_2, Cluster_3, Cluster_4)) %>%
  ggplot(aes(x=Sample, y=Value, fill=Cluster)) + facet_wrap(~CCR, scales="free")+
  geom_col() +
  theme(text = element_text(size = 10), 
        axis.text.x.bottom = element_text(size=10, angle=45, vjust = 0.85, hjust = 0.8),
        axis.title.y = element_text(vjust = 2.5))


K4_framed %>%
  mutate(Sample = fct_reorder(Sample, Cluster_4)) %>%
  gather(key=Cluster, value=Value, c(Cluster_1, Cluster_2, Cluster_3, Cluster_4)) %>%
  ggplot(aes(x=Sample, y=Value, fill=Cluster)) + facet_grid(~Region, scales="free_x", space ="free_x", switch="x")+
  geom_col() +
  theme(text = element_text(size = 10), 
        axis.text.x.bottom = element_text(size=10, angle=45, vjust = 0.85, hjust = 0.8),
        axis.title.y = element_text(vjust = 2.5))




###https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html#population-genetic-analyses-for-gbs-data


library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)


#Read in VCF
Bcin.VCF <- read.vcfR("Bcin1_snps_pass_names_ann.vcf.gz")
Bcin.VCF
queryMETA(Bcin.VCF)


#Read in population data
pop.data= read.table("populationdata_plate1.txt", sep = "\t", header=TRUE)
view(pop.data)


#Check that all samples are in VCF and population data
all(colnames(Bcin.VCF@gt)[-1] == pop.data$Sample)

#Convert vcf to genlight object
Bcin_gl <- vcfR2genlight(Bcin.VCF)
Bcin_gl
ploidy(Bcin_gl) <- 1
pop(Bcin_gl) <- pop.data$Crop
Bcin_gl

#Distance matrices (two options)
?dist
x.dist <- dist(Bcin_gl) #option 1
x.dist
?poppr::bitwise.dist(x) #option 2
x.dist <- poppr::bitwise.dist(Bcin_gl)



#Distance tree (takes too long)
tree <- aboot(Bcin_gl, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
?aboot
cols <- brewer.pal(n = nPop(Bcin_gl), name = "Dark2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(Bcin_gl)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("CA","OR","WA"), fill = cols, border = FALSE, bty = "n", cex = 2)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")




#PCA
?glPca
Bcin.pca <- glPca(Bcin_gl, nf = 4)
#Start=1:02PM
barplot(100*Bcin.pca$eig/sum(Bcin.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

Bcin.pca.scores <- as.data.frame(Bcin.pca$scores)
Bcin.pca.scores$pop <- pop(Bcin.pca)

library(ggplot2)
set.seed(9)
p <- ggplot(Bcin.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()

p


#Subsetting a vcfR object to X variants

alphabet <- c("A","B","C","D")
alphabet[c(1,2)]

Bcin.VCF[c(1:10),]

subset.1 <- sample(size = 200, x= c(1:nrow(Bcin.VCF)))
subset.2 <- sample(size = 200, x= c(1:nrow(Bcin.VCF)))
identical(subset.1, subset.2) #should be false to show that subsets are different

Bcin.VCF.sub1 <- Bcin.VCF[subset.1,]
Bcin.VCF.sub2 <- Bcin.VCF[subset.2,]
Bcin.VCF.sub1
Bcin.VCF.sub2


# Creating a list object to save our subsets in.
Bcin.variant.subset <- vector(mode = "list", length = 50)

# Using a for loop to generate 50 subsets of 200 random variants from the rubi.VCF vcfR object.
for (i in 1:50){
  Bcin.variant.subset[[i]] <- Bcin.VCF[sample(size = 200, x= c(1:nrow(Bcin.VCF)))]
}

# Checking we have 50 vcfR objects:
length(Bcin.variant.subset)

head(Bcin.variant.subset, n=2)


# Creating the GenLight object
Bcin.gl.subset <- lapply(Bcin.variant.subset, function (x) suppressWarnings(vcfR2genlight(x)))
for (i in 1:length(Bcin.gl.subset)){
  ploidy(Bcin.gl.subset[[i]]) <- 1
}

# Creating a simple UPGMA tree per object
install.packages("phangorn")
library(phangorn)
Bcin.trees <- lapply(Bcin.gl.subset, function (x) upgma(bitwise.dist(x)))
class(Bcin.trees) <- "multiPhylo"

# Overlapping the trees
densiTree(Bcin.trees, scaleX = T, show.tip.label = F, alpha = 0.1)
title(xlab = "Proportion of variants different")




