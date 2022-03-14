library(list)
files <- list.files(path = '~/Thesis_SpatialGenomics/VesaliusDev/RSource/', pattern = 'R')
files <- paste0('~/Thesis_SpatialGenomics/VesaliusDev/RSource/',files)
for(i in seq_along(files)){
  source(files[i])
}
