#SNPRelate
#Author: Nikki Lukasko
#Date: 3-8-2022
#Resource: https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#format-conversion-from-vcf-files
#Resuorce: http://www.bioconductor.org/packages/devel/bioc/vignettes/SeqArray/inst/doc/SeqArray.html#pca-r-implementation


library(ggplot2)
#if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("gdsfmt")
require(gdsfmt)
#BiocManager::install("SNPRelate")
require(SNPRelate)
  BiocManager::install("SeqArray")
library(SeqArray)



# the VCF file, using the example in the SeqArray package
vcf.fn <- "C:/Users/nikki/Michigan State University/PSM.Hausbecklab - Nikki Lukasko - Nikki Lukasko/Botrytis/Molecular/Sequencing Stats for R from HPCC/VCFs/Bcin_filtered.vcf"


# convert, save in "tmp.gds" with the default lzma compression algorithm
seqVCF2GDS(vcf.fn, "Bcin.gds") #take about 15m to run


genofile <- seqOpen("Bcin.gds")
# take out sample id
head(samp.id <- seqGetData(genofile, "sample.id"))

# take out variant id
head(variant.id <- seqGetData(genofile, "variant.id"))

# get "chromosome"
table(seqGetData(genofile, "chromosome"))

# get "allele"
head(seqGetData(genofile, "allele"))

# get "annotation/info/GP"
#head(seqGetData(genofile, "annotation/info/GP"))

# get "sample.annotation/family"
#head(seqGetData(genofile, "sample.annotation/family"))


set.seed(1000)
# may try different LD thresholds for sensitivity analysis
snpset1 <- snpgdsLDpruning(genofile, ld.threshold=0.2)
names(snpset1)
head(snpset1$chr3)
snpset.id <- unlist(snpset1)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
plot(pca, eig=1:4, pch=20, cex=0.5)


#add population data in the 4th column here:
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")


#adding population data
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
metadata <- read.delim("samplemetadata_GRCY.txt")
#can add fungicide resistance data in later



# PCA with crop coloring

tab2 <- data.frame(sample.id = pca$sample.id,
                  Crop = factor(metadata$Crop)[match(pca$sample.id, metadata$sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab2)
tail(tab2)

plot(tab2$EV2, tab2$EV1, col=as.integer(tab2$Crop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab2$Crop), pch="o", col=1:nlevels(tab2$Crop))

lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab2$Crop, labels=lbls)




# PCA with greenhouse coloring

tab3 <- data.frame(sample.id = pca$sample.id,
                   Greenhouse = factor(metadata$Greenhouse)[match(pca$sample.id, metadata$sample.id)],
                   EV1 = pca$eigenvect[,1],    # the first eigenvector
                   EV2 = pca$eigenvect[,2],    # the second eigenvector
                   stringsAsFactors = FALSE)
head(tab3)
tail(tab3)

plot(tab3$EV2, tab3$EV1, col=as.integer(tab3$Greenhouse), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab3$Greenhouse), pch="o", col=1:nlevels(tab3$Greenhouse))

lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab3$Greenhouse, labels=lbls)





# PCA with region coloring

tab4 <- data.frame(sample.id = pca$sample.id,
                   Region = factor(metadata$Region)[match(pca$sample.id, metadata$sample.id)],
                   EV1 = pca$eigenvect[,1],    # the first eigenvector
                   EV2 = pca$eigenvect[,2],    # the second eigenvector
                   stringsAsFactors = FALSE)
head(tab4)
tail(tab4)

plot(tab4$EV2, tab4$EV1, col=as.integer(tab4$Region), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab4$Region), pch="o", col=1:nlevels(tab4$Region))

lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab4$Region, labels=lbls)





----
  
  
# Relatedness analysis
  
CFG.id <- sample.id[metadata$Greenhouse == "CFG"]
# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, sample.id=CFG.id, snp.id=snpset.id, num.thread=2)
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1", main="CFG samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

TEB.id <- sample.id[metadata$Greenhouse == "TEB"]
# Estimate IBD coefficients
ibd2 <- snpgdsIBDMoM(genofile, sample.id=TEB.id, snp.id=snpset.id, num.thread=2)
ibd2.coeff <- snpgdsIBDSelection(ibd2)
head(ibd2.coeff)
plot(ibd2.coeff$k0, ibd2.coeff$k1, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1", main="TEB samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

#Identity by state analaysis


pop.code <- factor(seqGetData(genofile, "sample.annotation/Population"))
#need to integrate population data correctly...















