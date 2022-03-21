
min_cov <- 100
pass_cov <- sum(sparse(:,1:end-1),2) > min_cov

# NOT RUN {
data(iris)
iris_df = data.matrix(iris[-5]) ## convert to a numeric matrix 
Knn(m=iris_df, 4)

# }

knn_smoothing <- function(M, SNN, k = 10, method = "ave", strict = TRUE) {
  mat <- as.matrix( M ) # Give up SparseMatrix here.
  cname <- colnames( mat )
  gname <- rownames( mat )
  S <- mat
  if ( method == "ave" ) S <- .smootherAverageNN( mat, SNN, k, strict )
  else S <- .smootherSumNN( mat, SNN, k, strict )
  if (! is.null(cname)) colnames(S) <- cname
  if (! is.null(gname)) rownames(S) <- gname
  return ( S )
}