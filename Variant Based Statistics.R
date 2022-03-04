#Botrytis cinerea variant based statistics
#Data output from vcftools
#Author: Nikki Lukasko
#Date: 3-1-22
#Resource: https://speciationgenomics.github.io/filtering_vcfs/


library (tidyverse)
library (dplyr)





############## MAF filter = 0.1 ###############



# Variant Quality (Phred score) 


Bcin1_var_quality <- read_delim("./Bcin_filteredMAF.vcf.gz.lqual", delim = "\t",
                          col_names = c("chr", "pos", "qual"), skip = 1)
Bcin1_var_quality_hist <- ggplot(Bcin1_var_quality, aes(qual)) + geom_density(fill="dodgerblue1", color = "black",
                                                                          alpha=0.3) + theme_light()
Bcin1_var_quality_hist 




# Mean Depth

Bcin1_var_depth <- read_delim("./Bcin_filteredMAF.vcf.gz.ldepth.mean", delim = "\t",
                              col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
Bcin1_var_depth_hist <- ggplot(Bcin1_var_depth, aes(mean_depth)) + geom_density(fill="dodgerblue1", color = "black",
                                                                          alpha=0.3) + theme_light() + xlim(0,50)
Bcin1_var_depth_hist
summary(Bcin1_var_depth$mean_depth)




# Missingness (measure of how many individuals lack a genotype at a call site)

Bcin1_var_miss <- read_delim("./Bcin_filteredMAF.vcf.gz.lmiss", delim = "\t",
                            col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
Bcin1_var_miss_hist <- ggplot(Bcin1_var_miss, aes(fmiss)) + geom_density(fill="dodgerblue1", color = "black",
                                                                            alpha=0.3) + theme_light()
Bcin1_var_miss_hist
summary(Bcin1_var_miss$fmiss)



# Minor allele frequency (maf)

Bcin1_var_freq <- read_delim("./Bcin_filteredMAF.vcf.gz.frq", delim = "\t",
                          col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
Bcin1_var_freq$maf <- Bcin1_var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
Bcin1_var_freq_hist <- ggplot(Bcin1_var_freq, aes(maf)) + geom_density(fill="dodgerblue1", color = "black",
                                                                   alpha=0.3) + theme_light()
Bcin1_var_freq_hist
summary(Bcin1_var_freq$maf)



# Mean depth per individual


Bcin1_ind_depth <- read_delim("./Bcin_filteredMAF.vcf.gz.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
Bcin1_ind_depth_hist <- ggplot(Bcin1_ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", 
                                                colour = "black", alpha = 0.3, bins=50) + theme_light()
Bcin1_ind_depth_hist



# Proportion of missingness per individual

Bcin1_ind_miss  <- read_delim("./Bcin_filteredMAF.vcf.gz.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered",
                                      "nmiss", "fmiss"), skip = 1)
Bcin1_ind_miss_hist <- ggplot(Bcin1_ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1",
                          colour = "black", alpha = 0.3, bins=50) + theme_light()
Bcin1_ind_miss_hist
#individuals with a high amount of missingness may represent another species? or contamination?


# Heterozygosity and inbreeding per individual

Bcin1_ind_het <- read_delim("./Bcin_filteredMAF.vcf.gz.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
Bcin1_ind_het_hist <- ggplot(Bcin1_ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", 
                            colour = "black", alpha = 0.3, bins=50) + theme_light()
Bcin1_ind_het_hist
#expect this to be close to 1 in asexual populations






########### No MAF Filter ##################


setwd("C:/Users/nikki/Michigan State University/PSM.Hausbecklab - Nikki Lukasko - Nikki Lukasko/Botrytis/Molecular/Sequencing Stats for R from HPCC/No MAF Filter")


# Variant Quality (Phred score)

Bcin2_var_quality <- read_delim("./Bcin_filterednoMAF.vcf.gz.lqual", delim = "\t",
                                col_names = c("chr", "pos", "qual"), skip = 1)
Bcin2_var_quality_hist <- ggplot(Bcin2_var_quality, aes(qual)) + geom_density(fill="dodgerblue1", color = "black",
                                                                              alpha=0.3) + theme_light()
Bcin2_var_quality_hist 




# Mean Depth

Bcin2_var_depth <- read_delim("./Bcin_filterednoMAF.vcf.gz.ldepth.mean", delim = "\t",
                              col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
Bcin2_var_depth_hist <- ggplot(Bcin2_var_depth, aes(mean_depth)) + geom_density(fill="dodgerblue1", color = "black",
                                                                                alpha=0.3) + theme_light() + xlim(0,50)
Bcin2_var_depth_hist
summary(Bcin2_var_depth$mean_depth)




# Missingness (measure of how many individuals lack a genotype at a call site)

Bcin2_var_miss <- read_delim("./Bcin_filterednoMAF.vcf.gz.lmiss", delim = "\t",
                             col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
Bcin2_var_miss_hist <- ggplot(Bcin2_var_miss, aes(fmiss)) + geom_density(fill="dodgerblue1", color = "black",
                                                                         alpha=0.3) + theme_light()
Bcin2_var_miss_hist
summary(Bcin2_var_miss$fmiss)



# Minor allele frequency (maf)

Bcin2_var_freq <- read_delim("./Bcin_filterednoMAF.vcf.gz.frq", delim = "\t",
                             col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
Bcin2_var_freq$maf <- Bcin2_var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
Bcin2_var_freq_hist <- ggplot(Bcin2_var_freq, aes(maf)) + geom_density(fill="dodgerblue1", color = "black",
                                                                       alpha=0.3) + theme_light()
Bcin2_var_freq_hist
summary(Bcin2_var_freq$maf)



# Mean depth per individual


Bcin2_ind_depth <- read_delim("./Bcin_filterednoMAF.vcf.gz.idepth", delim = "\t",
                              col_names = c("ind", "nsites", "depth"), skip = 1)
Bcin2_ind_depth_hist <- ggplot(Bcin2_ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", 
                                                                             colour = "black", alpha = 0.3, bins=50) + theme_light()
Bcin2_ind_depth_hist



# Proportion of missingness per individual

Bcin2_ind_miss  <- read_delim("./Bcin_filterednoMAF.vcf.gz.imiss", delim = "\t",
                              col_names = c("ind", "ndata", "nfiltered",
                                            "nmiss", "fmiss"), skip = 1)
Bcin2_ind_miss_hist <- ggplot(Bcin2_ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1",
                                                                           colour = "black", alpha = 0.3, bins=50) + theme_light()
Bcin2_ind_miss_hist
#individuals with a high amount of missingness may represent another species? or contamination?


# Heterozygosity and inbreeding per individual

Bcin2_ind_het <- read_delim("./Bcin_filterednoMAF.vcf.gz.het", delim = "\t",
                            col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
Bcin2_ind_het_hist <- ggplot(Bcin2_ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", 
                                                                     colour = "black", alpha = 0.3, bins=50) + theme_light()
Bcin2_ind_het_hist
#expect this to be close to 1 in asexual populations





############ Grunwald statistics ##########


#Resource: https://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html



#download required packages
library('vcfR')
#version 2.12.0
library(ggplot2)
install.packages('Rcpp')
library(Rcpp)
install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library(adegenet)

?read.vcfR

#might need to increase R ram before reading in the vcf?
memory.limit
memory.limit(24000)
gc() #might help with large files?

#read in the vcf file
vcf_Bcin <- read.vcfR("Bcin_filtered_noMAFfilter.vcf.gz")
head(vcf_Bcin)
write.vcf(vcf_Bcin, "Bcin_filtered_noMAFfilter.vcf.gz")

#histogram of quality scores
qplot(getQUAL(vcf), geom = "histogram")
colnames(vcf@gt)

#convert file type to a genlight object for downstream compatibility
vcf_Bcin_genlight <- vcfR2genlight(vcf)

#check a few genotypes from the genlight object
gt <- extract.gt(vcf_Bcin, element = "GT")
gt[c(2,6,18), 1:3]

t(as.matrix(vcf_Bcin_genlight))[c(1,5,17), 1:3]


#VCF does not include population info, so create a vector that is the same length the number of samples.
pop()







