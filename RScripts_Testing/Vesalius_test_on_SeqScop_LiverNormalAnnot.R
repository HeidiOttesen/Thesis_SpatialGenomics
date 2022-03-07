#------------------------------------------------------------------------------#
# Testing out seq scope
#------------------------------------------------------------------------------#

library(imagerExtra)
library(imager)
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tvR)
library(sp)
library(grid)
library(parallel)
library(devtools)

## If using the published version of Vesalius
#install_github("patrickCNMartin/Vesalius")
library(vesalius)

##For using the develop-version of Vesalius
#source("~/VesaliusDev/quickload.R")


### Annotated seurat object
livn.seq <-readRDS("~/group/seqScope/Liver_normal_10um_annotated.rds")
cells <- as.character(Idents(livn.seq)) ##Hep_MT or _PC etc
#seq <- rgbPCA(SO =seq, slices=1,conserveSparse=F)

## This takes time! convert expression in seurat object form to R+G+B values (complex matrix -> simple table)
livn.seq <- rgbUMAP(SO = livn.seq, conserveSparse=F)

#Reduce pixel size? X and Y = location ()
livn.seq$x <- livn.seq$x/10
livn.seq$y <- livn.seq$y/10

## Add column with cell type (OBS repeat cell-types 3 times (R+G+B values))
cells <- data.frame(livn.seq,cells=rep(cells,times =3))

## What does your slide look like now? How many slides are inclluded?
x <- cells$x
y <- cells$y
plot(x, y, col='purple', xlab='X', ylab='Y', las=1)

## Keep only UMIs within a certain slide of the full picture 
cells1 <- cells %>% filter(cc==1 & x > 1 & x < 2900)
cells2 <- cells %>% filter(cc==1 & x > 2900 & x < 5800)
cells3 <- cells %>% filter(cc==1 & x > 5800 & x < 8700)
cells4 <- cells %>% filter(cc==1 & x > 8700 & x < 11600)
cells5 <- cells %>% filter(cc==1 & x > 11600 & x < 14500)
cells6 <- cells %>% filter(cc==1 & x > 14500 & x < 17400)

## Show chosen section/slide:
## x <- cells1$x
## y <- cells1$y
## plot(x, y, col='purple', xlab='X', ylab='Y', las=1)


## Why do you pick 8 colours from the Accent palettte?
ter_pal <- colorRampPalette(brewer.pal(8, "Accent")) #Function to colour the different types

## ter_col <- length(unique(cells1$cells)) # no of different Cell-types
## ter_col <- ter_pal(ter_col) # Use function to colour the different cell types
## Merge the two above lines:
# Use function to pick put colours for the different cell types
ter_col <- ter_pal(length(unique(cells1$cells))) 

## What is this??? Size?
## "override.aes	A list specifying aesthetic parameters of legend key"
cex <- 10

pdf("cellTypeLiverNormalSeqScope_2.pdf", width = 6, height= 5)
an <- ggplot(cells6, aes(x,y,col = cells))+ geom_point(size=1, alpha = 0.65)+
  theme_void() +
  scale_color_manual(values = ter_col)+
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size=12),
        plot.title = element_text(size =15),
        legend.position = "left")+
  guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
  labs(colour = "Cell Type", title = "Cell Types, slide 6")
print(an)
dev.off()


## I need to study these following steps much more!!!

## SliceID confusion - in livn.seq they are all set to slice 1. How does it know which of the 6 "slides"
## to work on? Or is it handling all 6 of them?


##Creates image array based on coloured coordinates
## Error in .checkVesalius(coordinates, sliceID, verbose) with dev version of vesalius: 
## this takes time! about 15 mins with ~900MB original dataset
image <- buildImageArray(livn.seq, sliceID=1, resolution=30, filterThreshold=0.99, cores= 10)
## resolution 30% of original image, 
## filterThreshold "numeric (range 0 -1) describing the quantile threshold at which barcodes and tiles 
## should be retained (seed details)"????
## 10 cores is a lot?

##image details are enhanced/amplified via colour histogram equalization (especially borders?)
tmpImage <- equalizeHistogram(image, sleft = 2.5 , sright = 2.5)
tmpImage <- regulariseImage(tmpImage, lambda = 10,invert =F)
tmpImage <- iterativeSegmentation.array(tmpImage,colDepth = 9,
                                        smoothIter = 10,
                                        method = c("box","iso"),
                                        acrossLevels = "mean",
                                        sigma = seq(1,3),
                                        box = 10,
                                        useCenter = T)
tmpImage <- isolateTerritories.array(tmpImage,minBar=5,captureRadius = 0.035)
plotTmp <- filter(tmpImage, x > 800 & x < 1700)
## what is this PAtrick??? Honestly you can do better.
## 12 hours straight of work is not an excuse for this kind of shitty code
for(i in seq_len(length(unique(plotTmp$territory)))){
  print(i)
  plotTmp$territory[plotTmp$territory == unique(plotTmp$territory)[i]] <- i
}

pdf("seqScopeTerrLiverNormal_ann.pdf", width = 7, height=6)
g <-territoryPlot(plotTmp,cex = 15,cex.pt=1)+
  theme_void()+
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size=12),
        plot.title = element_text(size =15),
        legend.position = "right")+
  labs(title= "Vesalius - Territories", colour= "Territory nr.")
print(g)
dev.off()

pdf("seqScopeLiverNormal_combined.pdf",width = 10,height=4.5)
print(an +g )
dev.off()

