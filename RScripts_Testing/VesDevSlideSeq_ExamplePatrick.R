#### Test Runs for vesalius
library(Matrix)
library(parallel)
library(imager)
library(imagerExtra)
library(Seurat)
library(tidyverse)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(tvR)
library(sp)
library(grid)

library(SuperpixelImageSegmentation)



# Standard libraries - just in case you are a command line fiend
library(utils)
library(stats)
library(graphics)
library(grDevices)

source("quickLoad.R")


#### Loading raw data
slideBeads <-c("~/group/slide_seqV2/Puck_200115_08_bead_locations.csv")
#slideBeads <-c("~/group/slide_seqV2/Puck_190926_03_bead_locations.csv")

slideBeads <- read.csv(slideBeads, header=TRUE)

slideCounts <- c("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz")
#slideCounts <- c("~/group/slide_seqV2/Puck_190926_03.digital_expression.txt.gz")
slideCounts <- read.table(slideCounts, header = TRUE )
rownames(slideCounts) <- slideCounts[,1]
slideCounts <- slideCounts[,-1]

### Loading vesalius object for quick test
load("~/Vesalius/vesalius.rda")
slideBeads <- GetTissueCoordinates(vesalius) %>% select(c("cells","x","y"))
colnames(slideBeads) <- c("barcodes","xcoord","ycoord")
slideCounts <- GetAssayData(vesalius, slot = "counts")


## Testing conversion and multi version integration
## this is older version approach
count <- CreateSeuratObject(slideCounts, assay ="Spatial")
bead <- ReadSlideSeq("~/group/slide_seqV2/Puck_200115_08_bead_locations.csv")
bead <- bead[Cells(x = count)]
DefaultAssay(object = bead) <- "Spatial"
count[["slice1"]] <- bead

count <- NormalizeData(count)
count <- FindVariableFeatures(count, nfeatures = 2000)
count <- ScaleData(count)
count <- rgbUMAP(count, pcs =30, conserveSparse = FALSE)



#### Building Vesalius Object
vesalius <- buildVesaliusObject(slideBeads,slideCounts)


ves <- buildVesaliusEmbeddings(vesalius,
                               method ="PCA",pcs = 30,
                               tensorResolution = 0.1,cores =10)

ves <- buildVesaliusEmbeddings(ves,
                               method ="PCA_L",pcs = 3,
                               tensorResolution = 0.1,cores =1)
ves <- buildVesaliusEmbeddings(ves,
                               norm = "SCT",
                               method ="UMAP",pcs = 30,
                               tensorResolution = 0.1,cores =1)
ves <- buildVesaliusEmbeddings(ves,
                               norm = "log",
                               method ="UMAP",pcs = 30,
                               tensorResolution = 0.1,cores =1)
#ves <- equalizeHistogram(ves, embedding = "PCA", dims = 1:30)
ves <- imageSegmentation(ves,dims=1:30, colDepth = 5,smoothIter=5,smoothType="iso",sigma =2,embedding = "PCA")

ves <- isolateTerritories(ves,trial = "last")

ves <- extractMarkers(ves,seed=7)

ves <- territoryMorphing(ves,c(5,6))

ves <- layerTerritory(ves,c(28,34,32),morphologyFactor=0,layerDepth=5)

viewGeneExpression(ves,gene= "Cryba4",as.layer=T,norm=F)


################# OLD Depreciated for now ######################################
mos <- buildMosaic(ves,dims = seq(1,30),
                   box=seq(3,11),
                   sigma=c(0.5,1,1.5),
                   iter=3,
                   maskThreshold=0.9,
                   maskType ="double",
                   method = "Kmeans",
                   k=5,
                   lambda =5,
                   smoothingMethod=c("iso","median"))

ter <- mos@territories
coord <- mos@tiles

ter <- right_join(coord, ter,by = "barcodes") %>% distinct(barcodes,.keep_all=TRUE)
sorted_labels <- order(levels(as.factor(ter$territory)))
sorted_labels[length(sorted_labels)] <- "isolated"
ter$territory <- factor(ter$territory, levels = sorted_labels)
ter_col <- length(levels(ter$territory))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- sample(ter_pal(ter_col),ter_col)
pdf("compositeEM.pdf", width = 16, height=12)
image <- ggplot(ter,aes(x,y,col = as.factor(territory)))+
  geom_point()+
  scale_color_manual(values = ter_col)
print(image)
dev.off()

pdf("embedsEM.pdf",height=32,width=18)
plotCols = 5
plotRows = 7

grid.newpage()
pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
vplayout <- function(x, y)
  viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
for (i in 1:30) {
  curRow = ceiling(i/plotCols)
  curCol = (i-1) %% plotCols + 1
  
  p <- imagePlotNew(ves, dims = i, cex = 12)
  print(p, vp = vplayout(curRow, curCol ))
  
}
dev.off()

pdf("Loadings.pdf",height=14,width=16)
plotCols = 4
plotRows = 4

grid.newpage()
pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
vplayout <- function(x, y)
  viewport(layout.pos.row = x, layout.pos.col = y)

# Make each plot, in the correct location
for (i in 1:16) {
  curRow = ceiling(i/plotCols)
  curCol = (i-1) %% plotCols + 1
  
  p <- imagePlotNew(vesl, dims = i, cex = 12)
  print(p, vp = vplayout(curRow, curCol ))
  
}
dev.off()




#### poly test