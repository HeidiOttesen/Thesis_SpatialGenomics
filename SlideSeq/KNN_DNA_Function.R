# Copied from KNN_DNA markdown file

KNN.DNA <- function(bead, df, k){
  ## Start
  #- Make sure you have the count df and bead info from the DNA_vesalius markdown..
  #- Barcodes and coordinates from bead_locations
  barcodes <- rownames(bead@coordinates)
  bead.coords <- bead@coordinates
  coords <- cbind(bead.coords$xcoord, bead.coords$ycoord)
  l <- length(barcodes)
  
  #Perform KNN:
  cat( paste(Sys.time()," Running K-Nearest Neighbor with k =",k, "number of neighbors", "\n"))
  coords.knearneigh <- knearneigh(coords, k = k)
  knn50 <- knn2nb(coords.knearneigh, row.names = barcodes)
  knn5 <- knn50[1:l] #Just the grouped indexes (to the barcodes) - not the other attributes
  
  ## Unique barcode knn list
  #- take the first barcode from the list and all of its neighbours and put them in the list.
  #- Remove those barcodes from your barcodes vector so that each barcode is only used once
  #- There are also sometimes leftovers not placed in any group
  
    
  
  bc <- barcodes
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
  
  ## Find median coordinate within knn group
  #- iterate over groups of knn 
  #- stores their corresponding x and y coordinates in vector t
  #- Stores the median x and y for each column in new list - comb
  #- makes a new list of new barcode names, bc.g - the first barcode in each group
  
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
  
  #Remove NA values and empty groups
  
  row.has.na <- apply(comb, 1, function(x){any(is.na(x))})
  comb <- comb[!row.has.na,]
  knn.bc <- knn5.bc1[!sapply(knn5.bc1,is.null)]
  
  ## Summing counts over knn-groups:
  #- Iterate over knn groups 
  #- take the counts for those barcodes for all bins from the sparsematrix
  #- df - Sparse matrix vector from the DNA vesalius script - Barcodes as columns, bins as rows.
  #- Sum the counts for each group and each bin 
  
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
  
  #Naming rows using new list of barcode names (first barcode of each knn group)
  bc.comb <- rownames(comb)
  

  # If reading slideseq bead info from file:
  #fname <- paste0(in.path, "knn/", alias, "_knn", k)
  #write.csv(comb, file = paste0(fname, ".bead_locations.csv"), row.names = TRUE, quote = FALSE)
  #cat(paste(Sys.time()," New bead location file saved as ", paste0(fname, ".bead_locations.csv"), "\n"))
  #ss <- ReadSlideSeq(paste0(fname, ".bead_locations.csv"))
  
  # If creating slideseq bead info as object in list instead:
  comb <- as.data.frame(comb)
  ss <- new(Class = 'SlideSeq',assay = "Spatial",
            coordinates = comb[,c("xcoord","ycoord")])
  rownames(ss@coordinates) <- rownames(comb)
  
  colnames(grMtrx) <- bc.comb
  rownames(grMtrx) <- bins
  knnSpMtx <- Matrix(grMtrx, sparse = TRUE)
  return(list(knnSpMtx, ss))
}
