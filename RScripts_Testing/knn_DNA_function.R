# Copied from KNN_DNA markdown file - see for code comments and more info

KNN.DNA <- function(bead, df, k){
  barcodes <- rownames(bead@coordinates)
  bead.coords <- bead@coordinates
  coords <- cbind(bead.coords$xcoord, bead.coords$ycoord)
  l <- length(barcodes)
  coords.knearneigh <- knearneigh(coords, k = k)
  knn50 <- knn2nb(coords.knearneigh, row.names = barcodes)
  knn5 <- knn50[1:l] #Just the grouped indexes (to the barcodes) - not the other attributes
  
  bc <- barcodes
  cat( paste(Sys.time()," Running K-Nearest Neighbor with k=",k, " number of neighbors", "\n"))
  knn5.bc <- vector(mode="list", length=l) #holds unique barcodes in groups of up to k
  doubles <- vector(mode="list", length=l) #Probably don't need this..
  knn5.bc1 <- vector(mode="list", length=l) #unique again but without NA elements (Still null groups)
  
  i <- 1
  
  cat(paste(Sys.time()," Using only unique barcodes in KNN groups", "\n"))
  while(i <= l){
    j <- 1
    while(j <= k){
      ix <- knn5[[i]][j]
      b <- barcodes[ix]
      if(b %in% bc){
        knn5.bc[[i]][j] <- b
      }else{
        doubles[[i]][j] <- b
      }
      j <- j + 1
      bc <- bc[!bc %in% b]
    }
    x <- knn5.bc[[i]]
    knn5.bc1[[i]] <- x[!is.na(x)]
    i <- i + 1
  }
  cat(paste(Sys.time()," Number of barcodes not included: ", length(bc), "\n"))
  i <- 1
  l <- length(knn5.bc1)
  comb <- matrix(0, nrow = l, ncol = ncol(bead.coords))
  bc.g <- vector(length=l)
  
  
  cat(paste(Sys.time()," Finding median coordinate for each knn group", "\n"))
  while(i <= l){
    tmp <- knn5.bc1[[i]]
    t <- bead.coords[tmp,]
    ltmp <- length(tmp)
    comb[i,] <- apply(t, 2, median)
    if(ltmp > 0){
      bc.g[[i]] <- tmp[[1]]
    }
    i <- i + 1 
  }
  
  rownames(comb) <- bc.g
  colnames(comb) <- c("xcoord", "ycoord")
  
  row.has.na <- apply(comb, 1, function(x){any(is.na(x))})
  comb <- comb[!row.has.na,]
  knn.bc <- knn5.bc1[!sapply(knn5.bc1,is.null)]
  spMtrx <- df
  i <- 1
  
  l <- length(knn.bc)
  cat(paste(Sys.time()," Number of KNN groups with only unique barcodes: ", l, "\n"))
  grMtrx <- matrix(0, nrow = nrow(spMtrx), ncol = l)
  cat( paste(Sys.time()," Summing the counts for each group and bin", "\n"))
  while(i <= l){
    tmp <- knn.bc[[i]]
    t <- spMtrx[,tmp]
    ltmp <- length(tmp)
    if(ltmp > 1){
      grMtrx[,i] <- rowSums(t)
    }else if(ltmp == 1){
      grMtrx[,i] <- t
    }
    i <- i + 1 
  }
  bc.comb <- rownames(comb)
  
  colnames(grMtrx) <- bc.comb
  rownames(grMtrx) <- bins
  knnSpMtx <- Matrix(grMtrx, sparse = TRUE)
  fname <- paste0(in.path, "knn/", alias, "_knn", k)
  write.csv(comb, file = paste0(fname, ".bead_locations.csv"), row.names = TRUE, quote = FALSE)
  cat(paste(Sys.time()," New bead location file saved as ", paste0(fname, ".bead_locations.csv"), "\n"))
  return(list(knnSpMtx, fname))
}
