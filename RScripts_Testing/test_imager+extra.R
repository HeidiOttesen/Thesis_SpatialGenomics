## Testing Imager + ImagerExtra

library(imager)
plot(boats)

b <- boats
class(boats) ##cimg - class type of normal image - all numbers

boats
grayscale(boats)
dim(boats) ## show the 4 dimensions of a cimg-class (W*H/pixels, depth/frames, col.channels)

#Do arithmetic to show it's an array with numbers:
log(boats)+3*sqrt(boats)
mean(boats)
sd(boats)

b2 <- b/2

layout(t(1:2)) # change plot layout to 1 row of 2 objects
plot(boats)
plot(boats/2) # plot function rescales RGB values as default

layout(t(1:2))
plot(boats,rescale=FALSE)
plot(boats/2,rescale=FALSE)
rgb(0,1,0)

# swap values between red and green
cscale <- function(r,g,b) rgb(g,r,b)
plot(boats,colourscale=cscale,rescale=FALSE)

#Map grayscale values to blue
cscale <- function(v) rgb(0,0,v)
grayscale(boats) %>% plot(colourscale=cscale,rescale=FALSE)

# function for returning colour codes from colour scale:
cscale <- scales::gradient_n_pal(c("red","purple","lightblue"),c(0,.5,1))
#cscale is now a function returning colour codes
cscale(0)
grayscale(boats) %>% plot(colourscale=cscale,rescale=FALSE)


## 2nd image:
fpath <- system.file('extdata/parrots.png',package='imager')
parrots <- load.image(fpath)
plot(parrots)

grayscale(boats) %>% hist(main="Luminance values in boats picture")
R(boats) %>% hist(main="Red channel values in boats picture")
G(boats) %>% hist(main="Green channel values in boats picture")
B(boats) %>% hist(main="Blue channel values in boats picture")
channel(boats,1) %>% hist(main="Red channel values in boats picture")

library(ggplot2)
library(dplyr)
bdf <- as.data.frame(boats)
head(bdf,3)
## Make new factor column based on the channels
bdf <- mutate(bdf,channel=factor(cc,labels=c('R','G','B')))
head(bdf,3)
ggplot(bdf,aes(value,col=channel))+geom_histogram(bins=30)+facet_wrap(~ channel)

## Normalization (Histogram Equalization)
## Random example:
## sample 100 random values (normal distribution) with mean=0 and SD=+/-1
x <- rnorm(100)
#layout(t(1:2))
hist(x,main="Histogram of x")
## empirical cumulative distribution function
f <- ecdf(x)
## Each pillar represent 10 of the samples in their empirical order
hist(f(x),main="Histogram of ecdf(x)")

## Normalization with the boat picture:
boats.g <- grayscale(boats)
boats.g ##still array
## 98304 pixels total
f <- ecdf(boats.g) ##order based on grayscale values
f
## Values from low 0.0845.. to high 0.997. gray-values or empirical?
plot(f,main="Empirical CDF of luminance values")

f(boats.g) %>% hist(main="Transformed luminance values")
f(boats.g) %>% str #not an image array, just a vector

f(boats.g) %>% as.cimg(dim=dim(boats.g)) %>% plot(main="With HistEq")
boats.g %>% plot(main="Witout HistEq")


## Histogram Equalizer for all color channels:
#Hist. equalisation for grayscale
hist.eq <- function(im) as.cimg(ecdf(im)(im),dim=dim(im))

#Split across colour channels,
cn <- imsplit(boats,"c")
cn #we now have a list of images
cn.eq <- map_il(cn,hist.eq) #run hist.eq on each
imappend(cn.eq,"c") %>% plot(main="All channels equalised") #recombine and plot

library(purrr)
#Convert to HSV, reduce saturation, convert back
RGBtoHSV(boats) %>% imsplit("c") %>%
  modify_at(2,~ . / 2) %>% imappend("c") %>%
  HSVtoRGB %>% plot(rescale=FALSE, main="RGB to HSV All channels")


## Make it a function:
desat <- function(im) RGBtoHSV(im) %>% imsplit("c") %>%
  modify_at(2,~ . / 2) %>% imappend("c") %>%
  HSVtoRGB

#Split image into 3 blocks, reduce saturation in 2nd block, recombine
im <- load.example("parrots")
imsplit(im,"x",3) %>% modify_at(2,desat) %>%
  imappend("x") %>% plot(rescale=FALSE)



## Edge detection
## search for gradients along the x and y axis
gr <- imgradient(boats.g,"xy")
gr
plot(gr,layout="row")
dx <- imgradient(boats.g,"x")
dy <- imgradient(boats.g,"y")
grad.mag <- sqrt(dx^2+dy^2)
plot(grad.mag,main="Gradient magnitude")

## Same as above in only one line: (using enorm - euclidian norm pixel-wise)
imgradient(boats.g,"xy") %>% enorm %>% plot(main="Gradient magnitude (again)")






## GG-plotting imager dataframes:
#convert image array to data.frame/table
df <- grayscale(boats) %>% as.data.frame

p <- ggplot(df,aes(x,y))+geom_raster(aes(fill=value))
p
## Gives a flipped boats image as plot
## ggplot plots in blue as default

## Reverse the Y-axis:
p + scale_y_continuous(trans=scales::reverse_trans())
## remove gray margin:
p <- p+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())
p

# Use grayscale instead of blur with ggplot
p+scale_fill_gradient(low="black",high="white")

df <- as.data.frame(boats)
## Plot all 3 channels separately
p <- ggplot(df,aes(x,y))+geom_raster(aes(fill=value))+facet_wrap(~ cc)
## Remove grayscale
p+scale_y_reverse()

## Adds columns: c.1 for channel 1 etc.. 
## so we have the values of each pixel next to each other in the same row
as.data.frame(boats,wide="c") %>% head

#Add new column with the translated colour HexCodes
df <- as.data.frame(boats,wide="c") %>% mutate(rgb.val=rgb(c.1,c.2,c.3))
head(df,3)

p <- ggplot(df,aes(x,y))+geom_raster(aes(fill=rgb.val))+scale_fill_identity()
p+scale_y_reverse()















## ImagerExtra
library(imagerExtra)
g <- grayscale(boats)
gd <- grayscale(dogs)
layout(matrix(1:2,1,2))
plot(g, main = "boats")
plot(gd, main = "dogs")


## Contrast Enhancement:
layout(matrix(1:4, 2, 2))
plot(g, main = "Original")
#EqualizePiecewise(im, N, smax, smin, range)
# N = how the gray level will be mapped
# smax = the upper limit of contrast stretching
EqualizePiecewise(g, 2) %>% plot(main = "N = 2")
EqualizePiecewise(g, 10) %>% plot(main = "N = 10")
EqualizePiecewise(g, 1000) %>% plot(main = "N = 1000")


## Color Balance (saturation)
layout(matrix(1:4, 2, 2))
plot(g, main = "Original")
## Saturates a percentage of the pixels on the left vs the right side of the "histogram"??
BalanceSimplest(g, 1, 1) %>% plot(main = "sleft = 1, sright = 1")
BalanceSimplest(g, 1, 20) %>% plot(main = "sleft = 1, sright = 20")



## SPE (Screened Poisson Equation) - Background correction
layout(matrix(1:2, 1, 2))
plot(papers, main = "Original")
SPE(papers, 0.1) %>% plot(main = "SPE (lamda = 0.1)")



## Denoising - remove noise (sharper image)
noisy <- g + imnoise(dim = dim(g), sd = 0.1)
layout(matrix(c(1,3,2,4), 2, 2))
plot(g, main = "Original")
plot(noisy, main = "Noisy Boats")
DenoiseDCT(noisy, 0.1) %>% plot(., main = "Denoised (8x8 window)")
DenoiseDCT(noisy, 0.1, flag_dct16x16 = TRUE) %>% plot(., main = "Denoised (16x16 window)")



## ThresholdAdaptive (Local adaptive thresholding)
# Extract objects from image

layout(matrix(1:2,1,2))
plot(papers, main = "Original")
hello <- ThresholdAdaptive(papers, 0.1, windowsize = 17, range = c(0,1))
plot(hello, main = "Binarizesd")





