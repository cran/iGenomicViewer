
makeTrack <- function(Broad.Band=NA,
                      Fine.Band=NA,
                      genomicLoc=NA,
                      geneName=NA){

  trackRegion = list()
  if(!is.na(Broad.Band[1])){
    if(class(Broad.Band) != "character") warning("Broad.Band should be character vector \n corresponding to broad.band label in associated mapObj Continuing with NA\n \n", immediate.=TRUE)
    if(class(Broad.Band) == "character") trackRegion$Broad.Band = Broad.Band
  }
  if(!is.na(Fine.Band[1])){
    if(class(Fine.Band) != "character") warning("Fine.Band should be character vector \n corresponding to fine.band label in associated mapObj \n Continuing with NA\n", immediate.=TRUE)
    if(class(Fine.Band) == "character") trackRegion$Fine.Band = Fine.Band
  }
  if(!is.na(geneName[1])){
    if(class(geneName) != "character") warning("geneName should be character vector \n corresponding to geneName label in associated annObj \n Continuing with NA\n", immediate.=TRUE)
    if(class(geneName) == "character") trackRegion$geneName = geneName
  }
  if(!is.na(genomicLoc[1])){
    if((class(genomicLoc) == "numeric") | (class(genomicLoc) == "integer")){
      warning("genomicLoc given as numeric \n Assuming in order of start, stop, start, stop, etc.\n  Placing in matrix \n", immediate.=TRUE)
      if((length(genomicLoc)/2 - round(length(genomicLoc)/2)) != 0){
        warning("Not even length, missing a location\n Continuing with NA\n")
      }else{
        genomicLoc = matrix(genomicLoc, ncol=2, byrow=TRUE)
      }
    }    
    if(class(genomicLoc) == "matrix"){
      trackRegion$genomicLoc = genomicLoc
    }
  } 
  
  return(trackRegion)
}

