

#####################################################################################
# function to make annotation object from individual annotation information objects
#####################################################################################

# currently add only one at a time
annotationObj <- function(annotation,  # see makeAnnotation
                          annotationObj = NA, 
                          obj.name = NA,
                          returnVl = TRUE,
                          saveVl = FALSE,
                          saveName="AnnotationObj.RData"
                          ){

  # if no annotationObj is given start new object
  #   otherwise will add annotation to given object
  if(class(annotationObj)=="logical"){
    cat("Creating new annotation object \n")
    annotationObj = list()
  }

  # check if individual annotation information object already
  #    exists in annotation object
  # NOTE: If it does exist the function will replace
  if(eval.js(paste("is.null(annotationObj$", obj.name, ")", sep=""))){
   eval.js(paste("annotationObj$", obj.name, "= annotation", sep=""))
  }else{
    warning(paste("annotation object: ", obj.name, " already exisits \n replacing object now\n"), immediate.=TRUE) 
    eval.js(paste("annotationObj$", obj.name, "= annotation", sep=""))
  } 

  class(annotationObj) <- "annobj"
  
  # save or return 
  if(saveVl) save(annotationObj, file=saveName, compress=TRUE)
  if(returnVl) return(annotationObj)

}



############################################################
# function to make annotation information object from file
############################################################
makeAnnotation <- function(file,
                           label,
                           chrom,
                           chrom.levels,
                           band.info=NA,
                           loc=NA,
                           loc.start=NA,
                           loc.stop=NA,
                           file.sep="\t",
                           additional=NA,
                           names.additional = NA,
                           links=NA,
                           names.links=NA,
                           images=NA,
                           names.images=NA,
                           returnVl = TRUE,
                           saveVl = FALSE,
                           saveName="Annotation.RData",
                           ...
                           ){

  # read annotation file
  cg = read.table(file, sep=file.sep, ...)

  # a genomic location for each annotation piece must be given
  #  [otherwise cannot match up to genomic mapping] 
  if(is.na(loc) & (is.na(loc.start) | is.na(loc.stop)) ) {
    stop("Must specify spot location \n")
  }
  
  # Retrieve labels and chromosome location for each piece of
  # annotation  [i.e geneName, chromosome]
  # this may be numeric indication of column or column header name
  cnames = names(cg)
  if(class(label)=="character") label = which(cnames == label)
  if(class(chrom)=="character") chrom = which(cnames == chrom)
  subC = c(label, chrom)
  annotation = cg[,subC]
  names(annotation) = c("Label", "Chrom")

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
    annotation$loc.center = cg[,loc]    
  }else{
    if(class(loc.start)=="character") loc.start = which(cnames == loc.start)
    if(class(loc.stop)=="character") loc.stop = which(cnames == loc.stop)
    subC = c(subC, loc.start, loc.stop)
    annotation$loc.start = cg[,loc.start]
    annotation$loc.stop = cg[,loc.stop]
    annotation$loc.center = (as.numeric(cg[,loc.start])+as.numeric(cg[,loc.stop]))/2
  }

  # additiional represents any additional columns in the annotation file besides
  #  label, chromosome, and genomic location that the user wishes to include in the
  #  individual annotation information object
  # NA will include all additional columns
  # 0 will include no additional columns 
  # this may be numeric indication of column or column header name

  if(!is.na(additional[1])){
    # if additional is 0 than no extra information is included
    if(additional[1] == 0){
      additional=0
    # columns have been specified by name or numeric indications
    }else{
      nm = names(annotation)
      if(class(additional) == "character"){
        mm = match(additional, names(cg))
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0) annotation = cbind(annotation, cg[,mm])
      }
      if(class(additional) == "numeric"){
        annotation = cbind(annotation, cg[,additional])  
      }
      # optional way to change column header information
      if(!is.na(names.additional[1])) names(annotation) = c(nm, names.additional)
    }
  # if additional is left as NA
  }else{
      # this would automatically add additional columns
    cdx = setdiff(1:dim(cg)[2], subC)
    annotation = cbind(annotation, cg[,cdx])
  }
  
  # way to include hyperlink information for annotation
  # if links is a table, assumes in correct order of annotation label
  if(!is.null(dim(links))){
    if(dim(links)[1] == dim(annotation)[1]){
      x.link=links
    }else{
      warning("links dimension does not match file length \n continuing with no links\n", immediate.=TRUE)
    }
    
  }else{
    # if links is not NA, links is [are] column[s] in file 
    #  this may be numeric indication of column or column header name
    if(!is.na(links[1])){
      nm = names(cg)
      if(class(links) == "character"){
        mm = match(links, names(cg))
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0) x.link = as.data.frame(cg[,mm])
      }
      if(class(links) == "numeric"){
        mm = links
        x.link = as.data.frame(cg[,links])  
      }
      #  optional way to change column header information
      if(!is.na(names.links[1])){
        names(x.link) =  names.links
      }else{
        names(x.link) = names(cg)[mm]
      }
    }else{
      # no links included
      x.link = NA
    }
  }

  names.links = names(x.link)




  # way to include images information for annotation
  # if images is a table, assumes in correct order of annotation label
  if(!is.null(dim(images))){
    if(dim(images)[1] == dim(annotation)[1]){
      x.images=images
    }else{
      warning("images dimension does not match file length \n continuing with no images\n", immediate.=TRUE)
    }
    
  }else{
    # if links is not NA, images is [are] column[s] in file 
    #  this may be numeric indication of column or column header name
    if(!is.na(images[1])){
      nm = names(cg)
      if(class(images) == "character"){
        mm = match(images, names(cg))
        if(length(which(is.na(mm))) > 0) mm = mm[-which(is.na(mm))]
        if(length(mm) != 0) x.images = as.data.frame(cg[,mm])
      }
      if(class(images) == "numeric"){
        mm = images
        x.images = as.data.frame(cg[,images])  
      }
      #  optional way to change column header information
      if(!is.na(names.images[1])){
        names(x.images) =  names.images
      }else{
        names(x.images) = names(cg)[mm]
      }
    }else{
      # no links included
      x.images = NA
    }
  }

  names.images = names(x.images)

  

  
  ############################
  # Band.Info information 
  ############################
  # if NA use program default 
  if(class(band.info) == "logical"){
    # uses default object
    data(Band.Info)
  }

  annObj = list()

  #
  # if locations are based on start and stoping genomic locations
  #
  if(is.na(loc)){
  ####################################################
  # add genomic locations and reorder
  ###################################################
    g.loc.start = annotation$loc.start
    g.loc.center =  annotation$loc.center
    g.loc.stop =  annotation$loc.stop

    # add buffers for previous chromosome so location is
    # over entire genome not just chromosome
    for(i in 1:length(chrom.levels)){
      idx = which(as.character(annotation$Chrom) == chrom.levels[i])
      g.loc.start[idx] = g.loc.start[idx] + band.info$offset[i]
      g.loc.stop[idx] = g.loc.stop[idx] + band.info$offset[i]
      g.loc.center[idx] = g.loc.center[idx] + band.info$offset[i]
    }
    annotation$g.loc.start = g.loc.start
    annotation$g.loc.center = g.loc.center
    annotation$g.loc.stop = g.loc.stop

    # reorder file listing
    # place in genomic order 
    reOrder = order(g.loc.center)
    annObj$annotation = annotation[reOrder,]
    #hyperlinks
    if(!is.null(dim(x.link))){
      annObj$links = as.data.frame(x.link[reOrder,])
      names(annObj$links) = names.links
    }else{
      annObj$links = NA
    }
    #images
    if(!is.null(dim(x.images))){
      annObj$images = as.data.frame(x.images[reOrder,])
      names(annObj$images) = names.images
    }else{
      annObj$images = NA
    }    
    
  # if locations are based on a single central location 
  }else{
    # reorder file listing
    # place in genomic order 
    annotation$g.loc.center = loc
    reOrder = order(annotation$g.loc.center)
    annObj$annotation = annotation[reOrder,]
    #hyperlinks
    if(!is.null(dim(x.link))){
      annObj$links = as.data.frame(x.link[reOrder,])
      names(annObj$links) = names.links
    }else{
      annObj$links = NA
    }      
    #images
    if(!is.null(dim(x.images))){
      annObj$images = as.data.frame(x.images[reOrder,])
      names(annObj$images) = names.images
    }else{
      annObj$images = NA
    }      
    
  }

  class(annObj) <- "anninfo"
  
  # save or return  
  if(saveVl) save(annObj, file=saveName, compress=TRUE)
  if(returnVl) return(annObj)
  
}









# retrieveAnnotation ?  see retrieve band.info
