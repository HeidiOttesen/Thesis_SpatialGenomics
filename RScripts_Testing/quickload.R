library(list)
files <- list.files(path = '~/R/VesaliusDev/RSource/', pattern = 'R')
files <- paste0('~/R/VesaliusDev/RSource/',files)
for(i in seq_along(files)){
  source(files[i])
}
