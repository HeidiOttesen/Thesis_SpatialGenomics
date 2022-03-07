## Imager/FFT

im <- as.cimg(function(x,y) sin(x/5)+cos(x/4)*sin(y/2),128,128)
ff <- FFT(im)
plot(ff$real,main="Real part of the transform")
plot(ff$imag,main="Imaginary part of the transform")
sqrt(ff$real^2+ff$imag^2) %>% plot(main="Power spectrum")
#Check that we do get our image back
check <- FFT(ff$real,ff$imag,inverse=TRUE)$real #Should be the same as original
mean((check-im)^2)



## ImagerExtra/DCT
DCT2D(imormat, returnmat = FALSE)

IDCT2D(imormat, returnmat = FALSE)

## imormat - image or matrix - a grayscale image of class cimg or a numeric matrix

## returnmat	- if returnmat is TRUE, returns numeric matrix. if FALSE, returns a grayscale image of class cimg.


## Example:
g <- grayscale(boats)
layout(matrix(1:2, 1, 2))
plot(g, main = "Original")
gg <- DCT2D(g) %>% plot(main = "Transformed")
gg %>% IDCT2D() %>% plot(main = "Inversed")
mean((g - gg)^2)
