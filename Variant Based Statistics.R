#Botrytis cinerea variant based statistis
#data from HPCC
#Author: Nikki Lukasko
#Date: 3-1-22
#Resource: https://speciationgenomics.github.io/filtering_vcfs/


library (tidyverse)
library (dplyr)


# Variant Quality (Phred score)

Bcin1_var_quality <- read_delim("./BU8_S96_aln_bwamem2_sort_rmvdups_calls_subset.lqual.lqual", delim = "\t",
                          col_names = c("chr", "pos", "qual"), skip = 1)
Bcin1_var_quality_hist <- ggplot(Bcin1_var_quality, aes(qual)) + geom_density(fill="dodgerblue1", color = "black",
                                                                          alpha=0.3) + theme_light()
Bcin1_var_quality_hist 




# Mean Depth

Bcin1_var_depth <- read_delim("./BU8.ldepth.mean", delim = "\t",
                              col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
Bcin1_var_depth_hist <- ggplot(Bcin1_var_depth, aes(mean_depth)) + geom_density(fill="dodgerblue1", color = "black",
                                                                          alpha=0.3) + theme_light() + xlim(0,100)
Bcin1_var_depth_hist
summary(Bcin1_var_depth$mean_depth)




# Missingness (measure of how many individuals lack a genotype at a call site)

Bcin1_var_miss <- read_delim("./BU8.lmiss", delim = "\t",
                            col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
Bcin1_var_miss_hist <- ggplot(Bcin1_var_miss, aes(fmiss)) + geom_density(fill="dodgerblue1", color = "black",
                                                                            alpha=0.3) + theme_light()
Bcin1_var_miss_hist
summary(Bcin1_var_miss$fmiss)



# Minor allele frequency (maf)

Bcin1_var_freq <- read_delim("./BU8.frq", delim = "\t",
                          col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
Bcin1_var_freq$maf <- Bcin1_var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
Bcin1_var_freq_hist <- ggplot(Bcin1_var_freq, aes(maf)) + geom_density(fill="dodgerblue1", color = "black",
                                                                   alpha=0.3) + theme_light()
Bcin1_var_freq_hist
summary(Bcin1_var_freq$maf)



# Mean depth per individual


Bcin1_ind_depth <- read_delim("./cichlid_subset.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
Bcin1_ind_depth_hist <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", 
                                                colour = "black", alpha = 0.3) 
                                                + theme_light()
Bcin1_ind_depth_hist



# Proportion of missingness per individual

Bcin1_ind_miss  <- read_delim("./cichlid_subset.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered",
                                      "nmiss", "fmiss"), skip = 1)
Bcin1_ind_miss_hist <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1",
                          colour = "black", alpha = 0.3) + theme_light()
Bcin1_ind_miss_hist



# Heterozygosity and inbreeding per individual

Bcin1_ind_het <- read_delim("./cichlid_subset.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
Bcin1_ind_het_hist <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", 
                            colour = "black", alpha = 0.3) + theme_light()
Bcin1_ind_het_hist


















