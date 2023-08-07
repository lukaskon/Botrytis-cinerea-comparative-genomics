#PCA of 276 Botrytis cinerea isolates


library(ggplot2)
library(ggfortify)

# Resource: https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
#           http://rstudio-pubs-static.s3.amazonaws.com/53162_cd16ee63c24747459ccd180f69f07810.html

setwd("C:/Users/nikki/Michigan State University/PSM.Hausbecklab - Nikki Lukasko - Nikki Lukasko/Botrytis/Molecular/Population Structure Analysis/Plink PCA")


#Basic PCA
pca_table <- read.table("BcinMAF05_PCA_results.eigenvec", header = TRUE, comment.char = "")
plot(pca_table[, c("PC1", "PC2", "PC3", "PC4", "PC5")])


#Detailed PCA
library(readxl)
Metadata_eigenvec <- read_excel("Metadata_eigenvec.xlsx", 
                                col_types = c("text", "text", "text", 
                                              "text", "text", "text", "text", "numeric", 
                                              "text", "text", "text", "numeric", 
                                              "text", "text", "text", 
                                              "text", "text", "text", 
                                              "text", "text", "numeric", 
                                              "numeric", "numeric", "numeric", 
                                              "numeric"))
View(Metadata_eigenvec)

plot(Metadata_eigenvec[, c("PC1", "PC2", "PC3", "PC4", "PC5"), color="Group S" ])

df <- Metadata_eigenvec[21:25]
pca_res <- prcomp(df, scale. = TRUE)
autoplot(pca_res)

summary(pca_res)

#wrong pca percentages- need to manually input

autoplot(pca_res, data = Metadata_eigenvec, colour = 'CCR')
autoplot(pca_res, data = Metadata_eigenvec, colour = 'Group_S')
autoplot(pca_res, data = Metadata_eigenvec, colour = 'Growing_Cycle')
autoplot(pca_res, data = Metadata_eigenvec, colour = 'Crop')
autoplot(pca_res, data = Metadata_eigenvec, colour = 'Greenhouse_simplified')
autoplot(pca_res, data = Metadata_eigenvec, colour = 'Region')
autoplot(pca_res, data = Metadata_eigenvec, colour = 'Fenhexamid')


autoplot(pca_res, data=Metadata_eigenvec, colour='Group_S', frame = TRUE, frame.type = 'norm')

# GroupS appears more 'diverse' than sensu stricto.


#Using ggplot may be easier to modify

#PC1 = 20.20%; PC2= 20.14%

install.packages("pals")
library(pals)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)

library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()

ggplot(Metadata_eigenvec,aes(x=PC1,y=PC2,col=Crop))+
  geom_point(size=3,alpha=0.5)+ 
  scale_color_brewer(palette="Dark2")+ 
  theme_classic()
ggplot(Metadata_eigenvec,aes(x=PC1,y=PC2,col=Group_S))+
  geom_point(size=3,alpha=0.5)+ 
  scale_color_brewer(palette="Dark2")+ 
  theme_classic()
ggplot(Metadata_eigenvec,aes(x=PC1,y=PC2,col=Greenhouse_simplified))+
  geom_point(size=3,alpha=0.8)+ 
  scale_fill_manual(values=as.vector(watlington(12)))+
  #scale_color_brewer(palette="RdYlGn")+ 
  theme_classic()
ggplot(Metadata_eigenvec,aes(x=PC1,y=PC2,col=Region))+
  geom_point(size=3,alpha=0.5)+ 
  scale_color_brewer(palette="Dark2")+ 
  theme_classic()