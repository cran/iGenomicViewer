#
# convert chromosome location to genomic location
#

# chromosome should be numeric only
# X and Y also numberic, ie. X=23, Y=24

# x should be numeric only,
# either numeric vector or numeric matrix

convertCloc <- function(x, chr, row=TRUE, bandobj=NA){

  if(is.na(bandobj[[1]][1])){
    data(Band.Info)
  }else{
    band.info = bandobj
  }

  # shouldnt need this b/c if X and Y NOT numeric 
  if(length(which(chr == "X")) != 0) chr[which(chr=="X")] = 23 
  if(length(which(chr == "Y")) != 0) chr[which(chr=="Y")] = 24
  
  if(is.null(dim(x))){
    if(length(chr)==1) chr = rep(chr, length(x))
    x = x + band.info$offset[chr]
  }else{
    
    if(row){
      if(length(chr)==1) chr = rep(chr, dim(x)[2])
      x[,1:dim(x)[2]] = x[,1:dim(x)[2]] + band.info$offset[chr]
    }
    if(!row){
      if(length(chr)==1) chr = rep(chr, dim(x)[1])
      x = t(x)
      x[,1:dim(x)[2]] = x[,1:dim(x)[2]] + band.info$offset[chr]
      x=t(x)
    }    
  }

  return(x)

}



convertGloc <- function(x, chr, row=TRUE, bandobj=NA){

  if(is.na(bandobj[[1]][1])){
    data(Band.Info)
  }else{
    band.info = bandobj
  }
  
  # shouldnt need this b/c if X and Y NOT numeric 
  if(length(which(chr == "X")) != 0) chr[which(chr=="X")] = 23 
  if(length(which(chr == "Y")) != 0) chr[which(chr=="Y")] = 24
  
  if(is.null(dim(x))){
    if(length(chr)==1) chr = rep(chr, length(x))
    x = x - band.info$offset[chr]
  }else{
    
    if(row){
      if(length(chr)==1) chr = rep(chr, dim(x)[2])
      x[,1:dim(x)[2]] = x[,1:dim(x)[2]] - band.info$offset[chr]
    }
    if(!row){
      if(length(chr)==1) chr = rep(chr, dim(x)[1])
      x = t(x)
      x[,1:dim(x)[2]] = x[,1:dim(x)[2]] - band.info$offset[chr]
      x=t(x)
    }    
  }

  return(x)

}
