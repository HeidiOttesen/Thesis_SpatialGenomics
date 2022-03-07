#------------------------------------------------------------------------------#
# Testing out seq scope segmented + DCT (+FFT)
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
#library(vesalius)

source("~/R/quickload.R")


### Annotated seurat object
seq <-readRDS("~/group/seqScope/Liver_normal_Segmentation_Hepatocytes.rds")
cells.id <- as.character(Idents(seq))
#seq <- rgbPCA(SO =seq, slices=1,conserveSparse=F)
seq <- rgbUMAP(SO = seq, conserveSparse=F)


seq$x <- jitter(seq$x/10,factor=0.1)
seq$y <- jitter(seq$y/10,factor=0.1)

cells <- data.frame(seq,cells=rep(cells.id,times =3))
cells2 <- cells %>% filter(cc==1 & x > 2800 & x < 5700)
cells1 <- cells %>% filter(cc==1 & x < 2800)
cells3 <- cells %>% filter(cc==1 & x > 5700 & x < 8600)
cells4 <- cells %>% filter(cc==1 & x > 8600)


##Heidi extra plotting:
layout(matrix(c(1,3,2,4), 2, 2))
an <- ggplot(cells2, aes(x,y,col = cells))+ geom_point(size=1, alpha = 0.65)+
  theme_void() +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size=12),
        plot.title = element_text(size =15),
        legend.position = "left")+
  guides(colour = guide_legend(override.aes = list(size=10 * 0.3)))+
  labs(colour = "Cell Type", title = "Cell Types")
print(an1 + an2 + an3 + an4)
print(an)


#ter_col <- length(unique(cells$cells))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- ter_pal(length(unique(cells$cells)))
cex <- 10


gg.an2 <- ggplot(cells2, aes(x,y,col = cells))+ geom_point(size=1, alpha = 0.65)+
  theme_void() +
  scale_color_manual(values = ter_col)+
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size=12),
        plot.title = element_text(size =15),
        legend.position = "left")+
  guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
  labs(colour = "Cell Type", title = "Cell Types")
#pdf("cellTypeLiverSeqScope.pdf", width = 6, height= 5)
layout(matrix(c(1,3,2,4), 2, 2))
#print(an.sect1 + an.sect2 + an1 + an2)
print(an.sect2 + gg.an2)
#dev.off()

##Removes any pixel overlaps (can create problem with the x,y divided by 10) Add jitter function before dividing to fix (get unique x,y).
image <- buildImageArray(seq,sliceID=1,resolution=30,filterThreshold=0.99, cores= 10)
im.eq <- equalizeHistogram(image, sleft = 2.5 , sright = 2.5)
im.reg <- regulariseImage(im.eq, lambda = 10,invert =F)
im.itSeg <- iterativeSegmentation.array(im.reg,colDepth = 9,
                                        smoothIter = 10,
                                        method = c("box","iso"),
                                        acrossLevels = "mean",
                                        sigma = seq(1,3),
                                        box = 10,
                                        useCenter = T)
im.terr <- isolateTerritories.array(im.itSeg,minBar=5,captureRadius = 0.035)



## Heidi Adding stuff
ss1 <- subset(image, select=c(x, y, cc, value))
ss.eq <- subset(im.eq, select=c(x, y, cc, value))
ss.reg <- subset(im.reg, select=c(x, y, cc, value))
ss.itSeg <- subset(im.itSeg, select=c(x, y, cc, value))




#### DCT ImagerExtra
pdf("DCT_SegmentedLiverSeqScope_comb.pdf", width = 6, height= 5)
layout(matrix(c(1,3,2,4), 2, 2))
g.or <- as.cimg(ss1) %>% plot(main="Original Array as.cimg")
g <- grayscale(g.or) %>% plot(main = "Original, gray")
gg1 <- DCT2D(g) %>% plot(main = "DCT Transformed")
gg2 <- IDCT2D(gg1) %>% plot(main = "DCT Inversed")
mean((g - gg2)^2)

g.eq <- as.cimg(ss.eq) %>% plot(main="Equalized as.cimg")
g <- grayscale(g.eq) %>% plot(main = "Equalized, gray")
gg1 <- DCT2D(g) %>% plot(main = "DCT Transformed")
gg2 <- IDCT2D(gg1) %>% plot(main = "DCT Inversed")
mean((g - gg2)^2)

g.reg <- as.cimg(ss.reg) %>% plot(main="Regularized as.cimg")
g <- grayscale(g.reg) %>% plot(main = "Regularized, gray")
gg1 <- DCT2D(g) %>% plot(main = "DCT Transformed")
gg2 <- IDCT2D(gg1) %>% plot(main = "DCT Inversed")
mean((g - gg2)^2)

g.seg <- as.cimg(ss.itSeg) %>% plot(main="Iteratively segmented as.cimg")
g <- grayscale(g.seg) %>% plot(main = "Segmented, gray")
gg1 <- DCT2D(g) %>% plot(main = "DCT Transformed")
gg2 <- IDCT2D(gg1) %>% plot(main = "DCT Inversed")
mean((g - gg2)^2)


dev.off()

## Min Max normalization:
normalize <- function(x) {return ((x - min(x)) / (max(x) - min(x)))}
itSeg.dct <- as.data.frame(gg1)
itSeg.dct$value <- normalize(itSeg.dct$value)
gg1 <- as.cimg(itSeg.dct) %>% plot(main="DCT Transformed, Normalized")

itSeg.dct.inv <- as.data.frame(gg2)
itSeg.dct.inv$value <- normalize(itSeg.dct.inv$value)
gg2 <- as.cimg(itSeg.dct.inv) %>% plot(main="DCT Inversed, Normalized")

## Log transformation
#logTrans <- function(mat){
#  cons <- 255/log(1+abs(max(mat)-min(mat)))
#  mat <- cons * log(mat)
#  return(mat)
#}
logTrans <- function(mat){
  mat <- mat + abs(min(mat)) + 1
  cons <- 255/log(1+abs(max(mat)-min(mat)))
  mat <- cons * log(mat)
  return(mat)
}
itSeg.dct <- as.data.frame(gg1)
itSeg.dct$value <- logTrans(itSeg.dct$value)
gg1 <- as.cimg(itSeg.dct) %>% plot(main="DCT Transformed, log-Normalized")




##FFT Imager
ss1.cim <- as.cimg(ss1)
tmp.fft <- FFT(ss1.cim)
layout(matrix(c(1,3,2,4), 2, 2))
plot(tmp.fft$real,main="Real part of the transform")
plot(tmp.fft$imag,main="Imaginary part of the transform")
sqrt(tmp.fft$real^2+tmp.fft$imag^2) %>% plot(main="Power spectrum")
# Do we get our image back?:
check <- FFT(tmp.fft$real,tmp.fft$imag,inverse=TRUE)$real #Should be the same as original
mean((check-im)^2)





## Vesalius resume:

plotTmp <- im.terr


for(i in seq_len(length(unique(plotTmp$territory)))){
  print(i)
  plotTmp$territory[plotTmp$territory == unique(plotTmp$territory)[i]] <- i
}

#pdf("seqScopeTerLiver.pdf", width = 7, height=6)
g <- territoryPlot(plotTmp,cex = 15,cex.pt=1)+
  theme_void()+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size=10),
        plot.title = element_text(size =12),
        legend.position = "right")+
  labs(title= "Vesalius - Territories", colour= "Territory nr.")
print(g)
#dev.off()

#pdf("seqScopeLiver_combined.pdf",width = 10,height=4.5)
print(an.sect2 + g)
dev.off()




