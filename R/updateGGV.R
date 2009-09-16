
updateGGV <- function(GGV,
                      trackRegions,
                      appendTo=TRUE,
                      returnVl=TRUE,
                      saveFlag=FALSE,
                      saveName="GGVobj.RData"){


  if(appendTo){
    otr = GGV$values$trackRegions 
    if(!is.null(trackRegions$Broad.Band)){
      if(is.null(otr$Broad.Band)) otr$Broad.Band = trackRegions$Broad.Band
      if(!is.null(otr$Broad.Band)) otr$Broad.Band = c(otr$Broad.Band, trackRegions$Broad.Band)
    }
    if(!is.null(trackRegions$Fine.Band)){
      if(is.null(otr$Fine.Band)) otr$Fine.Band = trackRegions$Fine.Band
      if(!is.null(otr$Fine.Band)) otr$Fine.Band = c(otr$Fine.Band, trackRegions$Fine.Band)
    }
    if(!is.null(trackRegions$geneName)){
      if(is.null(otr$geneName)) otr$geneName = trackRegions$geneName
      if(!is.null(otr$geneName)) otr$geneName = c(otr$geneName, trackRegions$geneName)
    }

    if(!is.null(trackRegions$genomicLoc)){
      if(is.null(otr$genomicLoc)) otr$geneName = trackRegions$genomicLoc
      if(!is.null(otr$genomicLoc)) otr$geneName = rbind(otr$genomicLoc, trackRegions$genomicLoc)
    }    

    GGV$values$trackRegions = otr
    
  }else{
    GGV$values$trackRegions = trackRegions
  }
  return(GGV)
}

