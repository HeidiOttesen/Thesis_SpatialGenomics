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

#install.packages("devtools")
library(devtools)

## If using the published version of Vesalius
#install_github("patrickCNMartin/Vesalius")
#library(vesalius)

##For using the develop-version of Vesalius
source("~/R/quickload.R")


### Annotated seurat object
seq <-readRDS("~/group/seqScope/Liver_normal_10um_annotated.rds")
cells <- as.character(Idents(seq)) ##Hep_MT or _PC etc
#seq <- rgbPCA(SO =seq, slices=1,conserveSparse=F)

## This takes time! convert expression to RGB colors (complex matrix -> simple table)
seq <- rgbUMAP(SO = seq, conserveSparse=F)
## Holds UMI (collapsed?), coordinates, color channel, value (of colour), slice?

#Reduce number of pixels
seq$x <- seq$x/10
seq$y <- seq$y/10

## Add column with cell type (OBS 49392 observations in seq, 16464 items in cells.. 
## 3 times for the three different color channels)
cells <- data.frame(seq,cells=rep(cells,times =3))
## Keep only UMIs with counts 2800<5700?
cells <- cells %>% filter(cc==1 & x > 2800 & x < 5700)


ter_col <- length(unique(cells$cells)) #13 Cell-types
ter_pal <- colorRampPalette(brewer.pal(8, "Accent")) #Function
ter_col <- ter_pal(ter_col) #Gives colour to the different cell types
# Define text and and symbol sizes. Scaling factor: (see more in ?par)
cex <- 10


an <- ggplot(cells, aes(x,y,col = cells))+ geom_point(size=1, alpha = 0.65)+
theme_void() +
scale_color_manual(values = ter_col)+
theme(legend.text = element_text(size = 12),
      legend.title = element_text(size=12),
      plot.title = element_text(size =15),
      legend.position = "left")+
guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
labs(colour = "Cell Type", title = "Cell Types")
pdf("cellTypeLiverSeqScope.pdf", width = 6, height= 5)
print(an)
dev.off()

## Error in .checkVesalius(coordinates, sliceID, verbose) : 
image <- buildImageArray(seq, sliceID=1, resolution=30, filterThreshold=0.99, cores= 10)
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

pdf("seqScopeTerLiver.pdf", width = 7, height=6)
g <-territoryPlot(plotTmp,cex = 15,cex.pt=1)+
   theme_void()+
   theme(legend.text = element_text(size = 12),
         legend.title = element_text(size=12),
         plot.title = element_text(size =15),
         legend.position = "right")+
  labs(title= "Vesalius - Territories", colour= "Territory nr.")
print(g)
dev.off()

pdf("seqScopeLiver_combined.pdf",width = 10,height=4.5)
print(an +g )
dev.off()

### Annotated seurat object
seq <-readRDS("~/group/seqScope/Colon_10um_annotated.rds")
cells <- as.character(Idents(seq))
#seq <- rgbPCA(SO =seq, slices=1,conserveSparse=F)
seq <- rgbUMAP(SO = seq, conserveSparse=F)


seq$x <- seq$x/10
seq$y <- seq$y/10

cells <- data.frame(seq,cells=rep(cells,times =3))
cells <- cells %>% filter(cc==1 & x > 6000 & x < 8500 & y<3000)


ter_col <- length(unique(cells$cells))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- ter_pal(ter_col)
cex <- 10

an1 <- ggplot(cells, aes(x,y,col = cells))+ geom_point(size=1, alpha = 0.65)+
theme_void() +
scale_color_manual(values = ter_col)+
theme(legend.text = element_text(size = 12),
      legend.title = element_text(size=12),
      plot.title = element_text(size =15),
      legend.position = "left")+
guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
labs(colour = "Cell Type", title = "Cell Types")
pdf("cellTypeColonSeqScope.pdf", width = 6, height= 5)
print(an1)
dev.off()


image <- buildImageArray(seq,sliceID=1,resolution=30,filterThreshold=0.99, cores= 10)
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


plotTmp <- filter(tmpImage, x > 1600 & x < 2500 & y <1000)
for(i in seq_len(length(unique(plotTmp$territory)))){
  print(i)
    plotTmp$territory[plotTmp$territory == unique(plotTmp$territory)[i]] <- i
}

pdf("seqScopeTerColon.pdf", width = 18, height=4.5)
g1 <-territoryPlot(plotTmp,cex = 15,cex.pt=1)+
   theme_void()+
   theme(legend.text = element_text(size = 12),
         legend.title = element_text(size=12),
         plot.title = element_text(size =15),
         legend.position = "right")
print(g1)
dev.off()


pdf("seqScopeColon_combined.pdf",width = 10,height=4.5)
print(an1 +g1 )
dev.off()


### Not annotated seurat object
seq <-readRDS("~/group/seqScope/Liver_normal_Segmentation_Hepatocytes.rds")

cells <- as.character(Idents(seq))
#seq <- rgbPCA(SO =seq, slices=1,conserveSparse=F)
seq <- rgbUMAP(SO = seq, conserveSparse=F)

seq_rgbUMAP <- seq 

seq$x <- seq$x/10
seq$y <- seq$y/10

#plot(seq$x, seq$y)

cells <- data.frame(seq,cells=rep(cells,times =3))


ter_col <- length(unique(cells$cells))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- ter_pal(ter_col)
cex <- 10

#an2 <- 
ggplot(cells, aes(x,y,col = cells))+ geom_point(size=1, alpha = 0.65)+
theme_void() +
scale_color_manual(values = ter_col)+
theme(legend.text = element_text(size = 12),
      legend.title = element_text(size=12),
      plot.title = element_text(size =15),
      legend.position = "left")+
guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
labs(colour = "Cell Type", title = "Cell Types")
pdf("cellTypeLiverCellsSeqScope.pdf", width = 6, height= 5)
print(an2)
dev.off()


image <- buildImageArray(seq,sliceID=1,resolution=30,filterThreshold=0.99, cores= 10)
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


plotTmp <- filter(tmpImage, x > 1600 & x < 2500 & y <1000)
for(i in seq_len(length(unique(plotTmp$territory)))){
  print(i)
    plotTmp$territory[plotTmp$territory == unique(plotTmp$territory)[i]] <- i
}

pdf("seqScopeTerColon.pdf", width = 18, height=4.5)
g1 <-territoryPlot(plotTmp,cex = 15,cex.pt=1)+
   theme_void()+
   theme(legend.text = element_text(size = 12),
         legend.title = element_text(size=12),
         plot.title = element_text(size =15),
         legend.position = "right")
print(g1)
dev.off()


pdf("seqScopeColon_combined.pdf",width = 10,height=4.5)
print(an1 +g1 )
dev.off()
