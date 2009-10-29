 
mappingObjADF <- function(adf,
                   spot.ID=NA,
                   chrom,
                   locBy,
                   base.chrm=NA,
                   reg.exp=NA,           
                   loc=NA,
                   loc.start=NA,
                   loc.stop=NA,
                   additional=NA,
                   names.additional=NA,
                   links=NA,
                   names.links=NA,
                   images=NA,
                   names.images=NA,
                   band.info = NA,
                   returnVl = TRUE,
                   saveFile = FALSE,
                   saveName="MapObj.RData"){

  # first check to make sure a genomic location is identified
  loc.found = FALSE
  loc.s.found=TRUE
  if(!is.na(loc[1])){
    loc.found=TRUE
    loc.s.found=FALSE
  }
  if(is.na(loc.start[1]) | is.na(loc.stop[1])){
    loc.s.found=FALSE
  }
  
  # a genomic location for each mapping piece must be given
  #if(is.na(loc) & (is.na(loc.start) | is.na(loc.stop)) ) {
  if(!loc.found & !loc.s.found){
    stop("Must specify spot location \n")
  }
  
  # names of columns
  cnames = names(adf@data)



  # get IDs 
  # if NA or invalid use row names
  # this may be numeric indication of column or column header name
  if(length(spot.ID) > 1){
    spot.ID = spot.ID[1]
    warning("Multiple column selected for ID \n  Using only first value \n", immediate.=TRUE)
  }  
  if(!is.na(spot.ID)){

    if(class(spot.ID)=="character") spot.ID = which(cnames == spot.ID)

    probes = try(as.character(adf@data[,spot.ID]))

    if(class(probes) == "try-error"){
      warning("Invalid column selected for ID \n  Using labels \n", immediate.=TRUE)
      spot.ID = NA
    }
  } 
  if(is.na(spot.ID)){
    probes = rownames(adf@data) 
  }


  # chromosome
  if(length(chrom) > 1){
    chrom = chrom[1]
    warning("Multiple column selected for chromosome \n  Using only first value \n", immediate.=TRUE)
  }
  if(!is.na(chrom)){

    if(class(chrom)=="character") chrom = which(cnames == chrom)

    chr = try(as.character(adf@data[,chrom]))
    if(class(chr) == "try-error"){
      warning("Invalid chromosome column selected \n", immediate.=TRUE)
      chrom = NA
    }
  }
  if(is.na(chrom)){
    stop("Chromosome Information not given \n")
  }

  
  subC = c(spot.ID, chrom)
  mapping.info = cbind(probes, chr)
  mapping.info = as.data.frame(mapping.info)
  names(mapping.info) = c("Spot.ID", "Chrom")




  
  if((locBy != "within") & (locBy != "across")){
    stop("Specifiy 'within' chromosome or 'across' entire genome for location values \n")
  }


    
  # get chr in correct format = could be 'chr1', 'chrom1'
  # in case of having to convert


  if(!is.na(base.chrm)){
    
    if(is.na(reg.exp)) reg.exp = rep(FALSE, length(base.chrm))
    if(length(reg.exp)==1) reg.exp = rep(reg.exp[1], length(base.chrm))
    
    for(i in 1:length(base.chrm)){
      chr = gsub(chr, pattern=base.chrm[i], replacement="",perl=reg.exp[i])
    }
    newchr=chr
    
    if(length(grep(newchr, pattern="X")>0)){
      newchr[grep(newchr, pattern="X")]= (max(suppressWarnings(as.numeric(newchr)), na.rm=TRUE)+1)
    }
    if(length(grep(newchr, pattern="x")>0)){
      newchr[grep(newchr, pattern="x")]= (max(suppressWarnings(as.numeric(newchr)), na.rm=TRUE)+1)
    }
    if(length(grep(newchr, pattern="Y")>0)){
      newchr[grep(newchr, pattern="Y")]= (max(suppressWarnings(as.numeric(newchr)), na.rm=TRUE)+1)
    }
    if(length(grep(newchr, pattern="y")>0)){
      newchr[grep(newchr, pattern="y")]= (max(suppressWarnings(as.numeric(newchr)), na.rm=TRUE)+1)
    }
  }


 
  #
  # genomic location must be given
  #   either
  #     a single location in loc argument - which will be taken as a central location
  #       ?assumed location is given over entire genome not within chromosome
  #   or
  #     a start and stop location in loc.start and loc.stop arguments 
  # this may be numeric indication of column or column header name
  if(!is.na(loc)){

    if(length(loc) > 1){
      loc = loc[1]
      warning("Multiple column selected for starting location \n  Using only first value \n", immediate.=TRUE)
    }   
    if(class(loc)=="character") loc = which(cnames == loc)
    locdata = try(adf@data[,loc])
    if(class(locdata) == "try-error"){
      warning("Invalid starting location column selected \n", immediate.=TRUE)
      loc = NA
    }
    if(is.na(loc)){
      stop("Location Information not valide \n")
    }
    subC = c(subC, loc)

    mapping.info$loc.center = as.numeric(locdata)
    # convert location if needed 
    if(locBy == "within") mapping.info$g.loc.center = convertCloc(as.numeric(locdata), as.numeric(newchr))
    if(locBy == "across") mapping.info$g.loc.center = as.numeric(locdata)
    
    mapping.info$g.loc.stop= mapping.info$g.loc.center
    mapping.info$g.loc.start= mapping.info$g.loc.center
    

  }else{

    # starting location 
    if(length(loc.start) > 1){
      loc.start = loc.start[1]
      warning("Multiple column selected for starting location \n  Using only first value \n", immediate.=TRUE)
    }
    # ending location
    if(length(loc.stop) > 1){
      loc.stop = loc.stop[1]
      warning("Multiple column selected for ending location \n  Using only first value \n", immediate.=TRUE)
    }
    if(class(loc.start)=="character") loc.start = which(cnames == loc.start)
    if(class(loc.stop)=="character") loc.stop = which(cnames == loc.stop)

    startData = try(adf@data[,loc.start])
    endData = try(adf@data[,loc.stop])

    if(class(startData) == "try-error"){
      warning("Invalid starting location column selected \n", immediate.=TRUE)
      loc.start = NA
    }
    if(is.na(loc.start)){
      stop("Starting Location Information not given \n", immediate.=TRUE)
    }
    if(class(loc.stop) == "try-error"){
      warning("Invalid ending location column selected \n", immediate.=TRUE)
      loc.stop = NA
    }
    if(is.na(loc.stop)){
      stop("Ending Location Information not given \n")
    }
    subC = c(subC, loc.start, loc.stop)


    mapping.info$loc.start=as.numeric(startData)
    mapping.info$loc.stop=as.numeric(endData)
    mapping.info$loc.center= (as.numeric(startData)+as.numeric(endData))/2
    
    # convert location if needed 
    if(locBy == "across"){
      mapping.info$g.loc.start = convertGloc(as.numeric(startData), as.numeric(newchr))
      mapping.info$g.loc.stop = convertGloc(as.numeric(endData), newchr)
      mapping.info$g.loc.center = (as.numeric(convertGloc(startData, as.numeric(newchr)))+as.numeric(convertGloc(endData, as.numeric(newchr))))/2
    }
    if(locBy == "within"){
      mapping.info$g.loc.start = as.numeric(startData)
      mapping.info$g.loc.stop = as.numeric(endData)
      mapping.info$g.loc.center = (as.numeric(startData)+as.numeric(endData))/2
    }

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
        mm = match(additional, cnames)
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0){
          extraInc = intersect(1:dim(adf@data)[2], mm)
          if(length(extraInc) > 0)mapping.info = cbind(mapping.info, adf@data[,extraInc])
          mm = extraInc
        }
      }
      if(class(additional) == "numeric" | class(additional) == "integer"){
        extraInc = intersect(1:dim(adf@data)[2], additional)
        mapping.info = cbind(mapping.info, adf@data[,additional])
        mm = additional
      }
      # optional way to change column header information
      if(!is.na(names.additional[1])){
        names(mapping.info) = c(nm, names.additional)
      }else{
        names(mapping.info) = c(nm, cnames[mm])
      }
    }
  # if additional is left as NA
  }else{

    naloc = which(is.na(subC))
    if(length(naloc) > 0) subC = subC[-naloc]
    
    # if additional is NA automatically adds all extra columns
    cdx = setdiff(1:dim(adf@data)[2], subC)
     adnam = names(mapping.info)
    
    mapping.info = cbind(mapping.info, adf@data[,cdx])  

  
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
        mm = match(links, cnames)
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0){
          linksInc = intersect(1:dim(adf@data)[2], mm)
          mm = linksInc
          if(length(linksInc)>0) x.link = as.data.frame(adf@data[,linksInc])
          else warning("no valid column[s] given for links \n", immediate.=TRUE)
        }
      }
      if(class(links) == "numeric" | class(additional) == "integer"){
        mm = links
        linksInc = intersect(1:dim(adf@data)[2], links)
        if(length(linksInc) > 0) x.link = as.data.frame(adf@data[,linksInc])
        else warning("no valid column[s] given for links \n", immediate.=TRUE)
        mm = linksInc
      }  
      #  optional way to change column header information
      if(!is.na(names.links[1])){
        names(x.link) =  names.links
      }else{
        names(x.link) = cnames[mm]
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
        mm = match(images, cnames)
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0){
          imagesInc = intersect(1:dim(adf@data)[2], mm)
          mm = imagesInc
          if(length(imagesInc)>0)x.images = as.data.frame(adf@data[,imagesInc])
          else warning("no valid column[s] given for images \n", immediate.=TRUE)
          mm = imagesInc
        }
      }
      if(class(images) == "numeric" | class(additional) == "integer"){
        mm = images
        imagesInc = intersect(1:dim(adf@data)[2], images)
        if(length(imagesInc) > 0)x.images = as.data.frame(adf@data[,imagesInc])
        else warning("no valid column[s] given for images \n", immediate.=TRUE)
        mm = imagesInc       
      }
      #  optional way to change column header information
      if(!is.na(names.images[1])){
        names(x.images) =  names.images
      }else{
        names(x.images) = cnames[mm]
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
    data(Band.Info, envir =environment())
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
    
    g.loc.start = mapping.info$g.loc.start
    g.loc.center =  mapping.info$g.loc.center
    g.loc.stop =  mapping.info$g.loc.stop
    
    # add buffers for previous chromosome so location is
    # over entire genome not just chromosome
    #for(i in 1:length(chrom.levels)){
    #  idx = which(as.character(mapping.info$Chrom) == chrom.levels[i])
    #  g.loc.start[idx] = g.loc.start[idx] + band.info$offset[i]
    #  g.loc.stop[idx] = g.loc.stop[idx] + band.info$offset[i]
    #  g.loc.center[idx] = g.loc.center[idx] + band.info$offset[i]
    #}
    #mapping.info$g.loc.start = g.loc.start
    #mapping.info$g.loc.center = g.loc.center
    #mapping.info$g.loc.stop = g.loc.stop

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
    # reorder file listing
    # place in genomic order 
    
    g.loc.center =  mapping.info$g.loc.center
    reOrder = order(g.loc.center)
    #reOrder = order(mapping.info$loc.center)

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




