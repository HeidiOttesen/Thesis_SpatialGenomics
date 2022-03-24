Knn(
  m,
  k,
  nThreads = 10L,
  verbose = TRUE,
  indexType = "angular",
  M = 12L,
  MaxM0 = 24L,
  ef_search_multiplier = 50,
  quiet = FALSE
)


library(sdpep)
#Let's make some sample data so people can run it. First a matrix of point coordinates:
#coords <- cbind(runif(100), runif(100))

coords <- bead@coordinates
coords <- cbind(coords$xcoord, coords$ycoord)

#coords <- cbind(ctmp[,2],ctmp[,1])
#coords <- coords[1:sub, ]

#Then suppose we have A measured at each of the 100 points:
# counts / weights to each coordinate:
#A <- runif(100)
A <- ctmp[,3]
a <- A[1:sub]

#And some ID values which we might need:

#IDs = paste0("ID-",1:100)
IDs <- bins[1:sub]
#Then we compute the nearest neighbours (I'm doing 5 nearest for simplicity):
ctmp.knearneigh <- knearneigh(coords, k = 5)

#knn50 <- knn2nb(knn, row.names = IDs)
#Each element of knn50 is a vector of the index in coords of the neighbours:
knn50[[1]]


#[1] 12 21 55 87 92
#So for example we can plot all the points, plot the neighbours of point 5 in green, then point 5 itself in red:
    
plot(coords)
points(coords[knn50[[5]],],col="green",pch=19)
points(coords[5,,drop=FALSE],col="red",pch=19)


#Then if we want to average over some measure at each of those points, here A, we can use sapply and a function that subsets from A according to those nearest elements:
  
sapply(1:length(knn50), function(N){mean(A[N])})

#giving the mean of A for the 5 nearest neighbours of each of the 100 points. Note the A value of the point itself isn't included.










# https://github.com/patrickCNMartin/VesaliusDev/blob/main/R/embeddings.R#:~:text=%23%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D%2D-,%23,-%23%20creating%20new%20merged


library("dbscan")
data(iris)
x <- iris[, -5]
nn <- kNN(x, k = 5)
ctmp.knn <- kNN(coords, k = 5)

