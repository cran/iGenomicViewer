

##########################################
# function to make MapObj from file
##########################################


mappingObj <- function(file,
                   spot.ID,
                   chrom,
                   chrom.levels,
                   loc=NA,
                   loc.start=NA,
                   loc.stop=NA,
                   file.sep="\t",
                   additional=NA,
                   names.additional=NA,
                   links=NA,
                   names.links=NA,
                   images=NA,
                   names.images=NA,
                   band.info = NA,
                   returnVl = TRUE,
                   saveFile = FALSE,
                   saveName="MapObj.RData",
                   ...){

  # reads file
  raw = read.table(file, sep=file.sep, ...)

  # a genomic location for each mapping piece must be given
  if(is.na(loc) & (is.na(loc.start) | is.na(loc.stop)) ) {
    stop("Must specify spot location \n")
  }

  # Retrieve labels and chromosome location for each piece of
  # mapping  [i.e BACs, Spots]
  # this may be numeric indication of column or column header name
  cnames = names(raw)
  if(class(spot.ID)=="character") spot.ID = which(cnames == spot.ID)
  if(class(chrom)=="character") chrom = which(cnames == chrom)
  subC = c(spot.ID, chrom)
  mapping.info = raw[,subC]
  names(mapping.info) = c("Spot.ID", "Chrom")

  #
  # genomic location must be given
  #   either
  #     a single location in loc argument - which will be taken as a central location
  #       assumed location is given over entire genome not within chromosome
  #   or
  #     a start and stop location in loc.start and loc.stop arguments 
  # this may be numeric indication of column or column header name
  if(!is.na(loc)){
    if(class(loc)=="character") loc = which(cnames == loc)
    subC = c(subC, loc)
    mapping.info$loc.center = raw[,loc]    
  }else{
    if(class(loc.start)=="character") loc.start = which(cnames == loc.start)
    if(class(loc.stop)=="character") loc.stop = which(cnames == loc.stop)
    subC = c(subC, loc.start, loc.stop)
    mapping.info$loc.start = raw[,loc.start]
    mapping.info$loc.stop = raw[,loc.stop]
    mapping.info$loc.center = (as.numeric(raw[,loc.start])+as.numeric(raw[,loc.stop]))/2
  }

  # additiional represents any additional columns in the mapping file besides
  #  label, chromosome, and genomic location that the user wishes to include in the
  #  mapping object
  # NA will include all additional columns
  # 0 will include no additional columns 
  # this may be numeric indication of column or column header name
  if(!is.na(additional[1])){
    # if additional is 0 then no additional columns are added
    if(additional[1] == 0){
      cat("\n")
    # columns have been specified by name or numeric indications
    }else{
      nm = names(mapping.info)
      if(class(additional) == "character"){
        mm = match(additional, names(raw))
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0) mapping.info = cbind(mapping.info, raw[,mm])
      }
      if(class(additional) == "numeric" | class(additional) == "integer"){
        mapping.info = cbind(mapping.info, raw[,additional])  
      }
      # optional way to change column header information
      if(!is.na(names.additional[1])){
        names(mapping.info) = c(nm, names.additional)
      }else{
        names(mapping.info) = c(nm, names(raw)[mm])
      }
    }
  # if additional is left as NA
  }else{
    # if additional is NA automatically adds all extra columns
    cdx = setdiff(1:dim(raw)[2], subC)
    adnam = names(mapping.info)

    mapping.info = cbind(mapping.info, raw[,cdx])  

  
    if(!is.na(names.additional[1])){
      names(mapping.info) = c(adnam, names.additional)
    }else{
      names(mapping.info) = c(adnam, cnames[cdx])
    }



  }



  # way to include hyperlink information for annotation
  # if links is a table, assumes in correct order of annotation label
  if(!is.null(dim(links))){
    if(dim(links)[1] == dim(mapping.info)[1]){
      x.link=links
    }else{
      warning("links dimension does not match file length \n continuing with no links\n", immediate.=TRUE)
    }
  }else{
    # if links is not NA, links is [are] column[s] in file 
    #  this may be numeric indication of column or column header name
    if(!is.na(links[1])){
      nm = names(mapping.info)
      if(class(links) == "character"){
        mm = match(links, names(raw))
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0) x.link = as.data.frame(raw[,mm])
      }
      if(class(links) == "numeric" | class(additional) == "integer"){
        mm = links
        x.link = as.data.frame(raw[,links])  
      }
      #  optional way to change column header information
      if(!is.na(names.links[1])){
        names(x.link) =  names.links
      }else{
        names(x.link) = names(raw)[mm]
      }
    }else{
      # no links included
      x.link = NA
    }
  }

  names.links = names(x.link)
    




  # way to include image information for annotation
  # if imagess is a table, assumes in correct order of annotation label
  if(!is.null(dim(images))){
    if(dim(images)[1] == dim(mapping.info)[1]){
      x.images=images
    }else{
      warning("images dimension does not match file length \n continuing with no images\n", immediate.=TRUE)
    }
  }else{
    # if images is not NA, images is [are] column[s] in file 
    #  this may be numeric indication of column or column header name
    if(!is.na(images[1])){
      nm = names(mapping.info)
      if(class(images) == "character"){
        mm = match(images, names(raw))
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0) x.images = as.data.frame(raw[,mm])
      }
      if(class(images) == "numeric" | class(additional) == "integer"){
        mm = images
        x.images = as.data.frame(raw[,images])  
      }
      #  optional way to change column header information
      if(!is.na(names.images[1])){
        names(x.images) =  names.images
      }else{
        names(x.images) = names(raw)[mm]
      }
    }else{
      # no links included
      x.images = NA
    }
  }

  names.images = names(x.images)









  MapObj = list()
   
  ############################
  # band.info information 
  ############################

  # if NA use program default 
  if(class(band.info) == "logical"){
    # uses default object
    data(Band.Info,envir =environment())
    MapObj$band.info = band.info 
    
  }else{
    # check to make sure object is bandinfo object
    if(class(band.info) == "bandinfo"){
      MapObj$band.info = band.info
    }else{
      warning("Band.Info is not of correct format\n Must be of the class bandinfo.\n See makeBandInfo for details\n Or using band.info=NA to use default object\n", immediate.=TRUE)
      band.info=NA
    }
  }


  #
  # if locations are based on start and stoping genomic locations
  #
  if(is.na(loc)){
  #################################################
  # make genomic locations and reorder mapping 
  #################################################
    g.loc.start = mapping.info$loc.start
    g.loc.center =  mapping.info$loc.center
    g.loc.stop =  mapping.info$loc.stop
    
    # add buffers for previous chromosome so location is
    # over entire genome not just chromosome
    for(i in 1:length(chrom.levels)){
      idx = which(as.character(mapping.info$Chrom) == chrom.levels[i])
      g.loc.start[idx] = g.loc.start[idx] + band.info$offset[i]
      g.loc.stop[idx] = g.loc.stop[idx] + band.info$offset[i]
      g.loc.center[idx] = g.loc.center[idx] + band.info$offset[i]
    }
    mapping.info$g.loc.start = g.loc.start
    mapping.info$g.loc.center = g.loc.center
    mapping.info$g.loc.stop = g.loc.stop

    # reorder file listing
    # place in genomic order 
    reOrder = order(g.loc.center)
    MapObj$mapping.info = mapping.info[reOrder,]
    # hyperlinks
    if(!is.null(dim(x.link))){
      MapObj$links = as.data.frame(x.link[reOrder,])
      names(MapObj$links) = names.links
    }else{
      MapObj$links = NA
    }
    # images
    if(!is.null(dim(x.images))){
      MapObj$images = as.data.frame(x.images[reOrder,])
      names(MapObj$images) = names.images
    }else{
      MapObj$images = NA
    }
  # if locations are based on a single central location 
  }else{

    g.loc.center =  mapping.info$loc.center
    # add buffers for previous chromosome so location is
    # over entire genome not just chromosome
    for(i in 1:length(chrom.levels)){
      g.loc.center[idx] = g.loc.center[idx] + band.info$offset[i]
    }
    mapping.info$g.loc.center = g.loc.center
    mapping.info$g.loc.start = g.loc.center
    mapping.info$g.loc.stop = g.loc.center
    
    # reorder file listing
    # place in genomic order 
    reOrder = order(mapping.info$g.loc.center)
    MapObj$mapping.info = mapping.info[reOrder,]
    # hyperlinks
    if(!is.null(dim(x.link))){
      MapObj$links = as.data.frame(x.link[reOrder,])
      names(MapObj$links) = names.links
    }else{
      MapObj$links = NA
    }      
    # hyperlinks
    if(!is.null(dim(x.images))){
      MapObj$images = as.data.frame(x.images[reOrder,])
      names(MapObj$images) = names.images
    }else{
      MapObj$images = NA
    }      
  }

  class(MapObj) <- "mapobj"
  
  # save or return  
  if(saveFile) save(MapObj, file=saveName, compress=TRUE)
  if(returnVl) return(MapObj)
  
}











