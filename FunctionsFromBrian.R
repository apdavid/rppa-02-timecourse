makeOneWayAnovaContrastMatrix <- function(design){
  cont.matrix <- matrix(nrow = length(colnames(design)), ncol = (length(colnames(design)) - 1)*(length(colnames(design)))/2) # levels X Contrasts
  prevNumContrasts <- 1
  for ( i in 1:nrow(cont.matrix)){
    numContrasts <- nrow(cont.matrix) - i
    endContrasts <- prevNumContrasts + numContrasts - 1
    if (numContrasts > 0){
      cont.matrix[i,c(prevNumContrasts:endContrasts)] <- 1
    }
    k <- i
    l <- i + numContrasts - 1
    if ( k <= l){
      for ( j in k:l){
        cont.matrix[j+1,c(prevNumContrasts:endContrasts)[j - i + 1]] <- -1      
      }
    }
    prevNumContrasts <- prevNumContrasts + numContrasts
    
  }  
  cont.matrix[is.na(cont.matrix)] <- 0
  row.names(cont.matrix) <- colnames(design)
  colnames(cont.matrix) <- paste("Contrast", c(1:ncol(cont.matrix)), sep = "")
  return(cont.matrix)
}

scaleHeatMapColor <- function(fcMatrix, lower, upper){
  if(missing(lower)){
    LOWER <- "10%"
  }else{
    LOWER <- lower
  }
  if(missing(upper)){
    UPPER <- "90%"
  }else{
    UPPER <- upper
  }
  quantile.range <- quantile(fcMatrix, probs = seq(0, 1, 0.01))
  palette.breaks.neg <- seq(quantile.range[LOWER], 0, (0 - quantile.range[LOWER])/25)
  palette.breaks.pos <- seq(0, quantile.range[UPPER], (quantile.range[UPPER] - 0)/25)
  
  # use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
  color.palette  <- colorRampPalette(c("blue","white", "red"))(50)
  
  col <- color.palette
  breaks <- c(palette.breaks.neg, palette.breaks.pos[2:length(palette.breaks.pos)])
  return(list("col" = col, "breaks" = breaks))
}
makeRowd <- 
  function(intensities){
    rowcor <- as.dist(1-cor(t(intensities)))
    rowhc <- hclust(rowcor, method="average")
    rowd <- as.dendrogram(rowhc)
    return(rowd)
  }