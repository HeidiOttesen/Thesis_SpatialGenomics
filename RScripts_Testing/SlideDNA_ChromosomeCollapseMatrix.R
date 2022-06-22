
chromosomeCounts <- function(bead, df){
  barcodes <- rownames(bead@coordinates)
  bead.coords <- bead@coordinates
  coords <- cbind(bead.coords$xcoord, bead.coords$ycoord) # I also need to include the chromosome info somewhere?!
  
  bc <- vector(mode="list", length=length(barcodes)) 
  cat(paste(Sys.time()," Change to barcode names (instead of barcode index numbers)", "\n"))
  for(i in seq_along(barcodes)){
    bc[[i]] <- barcodes[coords[i,]]
  }
  
  
  
  
 
  ## Summing counts over knn-groups:
  #- Iterate over knn groups - take their barcodes 
  #- take the counts for those barcodes for all bins from the sparsematrix df
  #- df - Sparse matrix vector from the DNA_vesalius markdown script - Barcodes as columns, bins as rows.
  #- Sum the counts for each group and each bin 
  grMtrx <- matrix(0, nrow = nrow(df), ncol = ncol(df))
  cat( paste(Sys.time()," Summing the counts for each KNN group and genomic bin", "\n"))
  for(i in seq_along(knn.bc)){
    grMtrx[,i] <- rowSums(df[,knn.bc[[i]]])
  }
  
  # If creating bead info and count matrix as object in list instead of exporting to files:
  colnames(grMtrx) <- rownames(bead.coords)
  rownames(grMtrx) <- bins
  #knnSpMtx <- Matrix(grMtrx, sparse = TRUE)
  return(list(grMtrx, bead.coords)) 
}






#' KNN.DNA - Perform K nearest neighbor on count matrix with coordinates
#' @param bead SlideSeq object containing barcode coordinate information
#' @param df Matrix containing the sparse counts of bins
#' @param k Number of neighbors to sum counts of

KNN.DNA <- function(bead, df, k){
  barcodes <- rownames(bead@coordinates)
  bead.coords <- bead@coordinates
  coords <- cbind(bead.coords$xcoord, bead.coords$ycoord)
  
  #Perform KNN:
  cat(paste(Sys.time()," Running K-Nearest Neighbor with k =",k, "number of neighbors \n on", length(barcodes), "barcoded bead locations", "\n"))
  coords.knearneigh <- knearneigh(coords, k = k)
  knnIx <- coords.knearneigh$nn
  
  #Replace barcode index number with barcode string
  knn.bc <- vector(mode="list", length=length(barcodes)) #holds unique barcodes in groups of up to k
  cat(paste(Sys.time()," Change to barcode names (instead of barcode index numbers)", "\n"))
  for(i in seq_along(barcodes)){
    knn.bc[[i]] <- barcodes[knnIx[i,]]
  }
  
  
  ## Summing counts over knn-groups:
  #- Iterate over knn groups - take their barcodes 
  #- take the counts for those barcodes for all bins from the sparsematrix df
  #- df - Sparse matrix vector from the DNA_vesalius markdown script - Barcodes as columns, bins as rows.
  #- Sum the counts for each group and each bin 
  grMtrx <- matrix(0, nrow = nrow(df), ncol = ncol(df))
  cat( paste(Sys.time()," Summing the counts for each KNN group and genomic bin", "\n"))
  for(i in seq_along(knn.bc)){
    grMtrx[,i] <- rowSums(df[,knn.bc[[i]]])
  }
  
  # If creating bead info and count matrix as object in list instead of exporting to files:
  colnames(grMtrx) <- rownames(bead.coords)
  rownames(grMtrx) <- bins
  #knnSpMtx <- Matrix(grMtrx, sparse = TRUE)
  return(list(grMtrx, bead.coords))
}

