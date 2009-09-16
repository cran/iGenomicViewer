

##########################################
#
# bioinforamtics wrapper to sendplot
#
###########################################


#
# assumes TIplot locations are with respect to genome
#  not chromosome
#





#
# this add side plot
#

iGGVtiled <- function(TIplot,
                      
                  annObj,

                  x.labels=NA,
                  y.labels=NA,
                  xy.labels=NA,
                  x.links=NA,
                  y.links=NA,
                  xy.links=NA,
                  asLinks=NA,
                      
                  x.images=NA,
                  y.images=NA,
                  xy.images=NA,

                  mat=NA,    
                  mai.mat = NA,
                  mai.prc=FALSE,

                     
                  plot.extras=NA,    
                  smpLines=TRUE,
                  divCol="lightgrey",
                 
                  plot.call=NA,
                  plot.vec=NA,

                  lims = c(-0.5,0.5),
                  annotation = NA,
                  clrs=c("blue", "hotpink", "purple", "orange"),

                  mapObj.columns = NA,

                  fname.root="iGGV",
                  dir="./",
                  overwriteSourcePlot = NA,
                  makeInteractive=TRUE,
                  overrideInteractive=NA,
                  header="v3",
                  window.size = "800x1100",
                  image.size= "800x1100",

                  cleanDir=TRUE,

                  vrb = TRUE, 
                  ... # extra arguments to makeImap not in above 
                  ){

   

  mapObj = TIplot$map$mapObj
  smplDX = TIplot$vls$smplDX

  
  if(is.null(dim(mat))){
    #if(is.na(plot.call)) mat = matrix(c(rep(c(rep(1,7), rep(3,1)), 8), rep(c(rep(2,7), rep(0,1)), 2)), byrow=TRUE, ncol=8)
    if(is.na(plot.call)) mat = matrix(c(rep(c(rep(1,10), rep(3,1)), 10), rep(c(rep(2,10), rep(0,1)), 2)), byrow=TRUE, ncol=11)
    #if(!is.na(plot.call)) mat = matrix(c(rep(c(rep(1,7), rep(3,1), rep(4,2)), 9), rep(c(rep(2,7), rep(0,3)), 2)), byrow=TRUE, ncol=10)
    if(!is.na(plot.call)) mat = matrix(c(rep(c(rep(1,10), rep(3,1), rep(4,2)), 13), rep(c(rep(2,10), rep(0,3)), 2)), byrow=TRUE, ncol=13)
  }else{
    mat = mat
  }
  


  
  # heatmap
  plot1 = "makeTiled(TIplot, smpDiv=smpLines, divCol=divCol)"


  zlgnd=array(seq(from=TIplot$lims$zlim[1],to=TIplot$lims$zlim[2],length=1000),dim=c(1,1000))

  # legend
  plot2="image(x=as.vector(zlgnd),y=1,z=t(zlgnd),zlim=c(min(zlgnd,na.rm=TRUE),max(zlgnd, na.rm=TRUE)),col=c(hsv(h=2/6,v=seq(1,0,length=length(zlgnd))^1.15),hsv(h=0/6,v=seq(0,1,length=length(zlgnd))^1.15)), axes=FALSE,xlab='',ylab='')"
  plot2 = paste(plot2, paste("axis(1,seq(from=",min(zlgnd),",to=",max(zlgnd),",length=5),line=0)", sep=""), sep=";")



  # bug in following --
  #  do to overlapping BACs
  
  #g.range = TIplot$lims$ylim
  #scanDX = intersect(which(mapObj$mapping.info$g.loc.stop >= g.range[1]), which(mapObj$mapping.info$g.loc.start <= g.range[2]))
  #scanLoc = mapObj$mapping.info$g.loc.center[scanDX]
  #while(length(which(diff(scanLoc)<=0)) > 0){
  #  scanLoc[which(diff(scanLoc)<=0)]=scanLoc[which(diff(scanLoc)<=0)]-1
  #}
  
  scanDX = TIplot$map$bacDX
  scanLoc = mapObj$mapping.info$g.loc.center[scanDX]
  while(length(which(diff(scanLoc)<=0)) > 0){
    scanLoc[which(diff(scanLoc)<=0)]=scanLoc[which(diff(scanLoc)<=0)]-1
  }
  
 # annotation
 if(is.na(annotation[1])) annotation = 1:length(annObj)

  plot3 = "image(x=0:1,y=0:1,z=matrix(rep(NA,4),ncol=2),xlim=c(0,length(annotation)),ylim=range(scanLoc,na.rm=TRUE),zlim=c(0,1),axes=FALSE,xlab='',ylab='')"  

  annotationDX = list()
  
  for(a in 1:length(annotation)){
    obj = annObj[annotation[a]][[1]][[1]]
    if(!is.null(obj$g.loc.start)){
      #ann.idx = intersect(which(obj$g.loc.start > range(scanLoc,na.rm=TRUE)[1]), which(obj$g.loc.stop < range(scanLoc,na.rm=TRUE)[2]))
      ann.idx = intersect(which(obj$g.loc.stop >= range(scanLoc,na.rm=TRUE)[1]), which(obj$g.loc.start <= range(scanLoc,na.rm=TRUE)[2]))
      
      if(length(ann.idx) !=0){
        y.bottom = rep(NA, length(ann.idx))
        y.pos = rep(NA, length(ann.idx))
        for(i in 1:length(ann.idx)){
          plot3 = paste(plot3, paste("lines(x=c(",a-.5, ",", a-.5, "), y=c(",obj$g.loc.start[ann.idx[i]],",", obj$g.loc.stop[ann.idx[i]],"), col='",clrs[a],"',lwd=5)",sep=""), sep=";")
          y.bottom[i] = obj$g.loc.start[ann.idx[i]]
          y.pos[i] =  obj$g.loc.stop[ann.idx[i]]
        }
        x.pos = rep((a-.7), length(y.pos))
        x.right = x.pos+.4
        annDX = list()
        annDX$type = "rect"
        annDX$idx = ann.idx 
        annDX$y.bottom = y.bottom
        annDX$y.pos = y.pos
        annDX$x.pos = x.pos
        annDX$x.right = x.right
        #hyperlink
        if(!is.null(dim(annObj[annotation[a]][[1]][[2]]))){
          annDX$links = as.data.frame(annObj[annotation[a]][[1]][[2]][ann.idx,])
          names(annDX$links) = names(annObj[annotation[a]][[1]][[2]])
        }else{
          annDX$links = NA
        }
        # images
        if(!is.null(dim(annObj[annotation[a]][[1]][[3]]))){
          annDX$images = as.data.frame(annObj[annotation[a]][[1]][[3]][ann.idx,])
          names(annDX$images) = names(annObj[annotation[a]][[1]][[3]])
        }else{
          annDX$images = NA
        }

        
        eval.js(paste("annotationDX$annObjDX",a,"=annDX", sep=""))
      }else{
        eval.js(paste("annotationDX$annObjDX",a,"$type=NA", sep=""))
      }
    }else{
      ann.idx = intersect(which(obj$g.loc.center > range(scanLoc,na.rm=TRUE)[1]), which(obj$g.loc.center < range(scanLoc,na.rm=TRUE)[2]))
      if(length(ann.idx) !=0){
        x.pos = rep(NA, length(ann.idx))
        y.pos = rep(NA, length(ann.idx))
        for(i in ann.idx){
          plot3 = paste(plot3, paste("points(x=",a-.5,", y=",obj$g.loc.center[i],",col='",clrs[a],"',lwd=5)",sep=""), sep=";")
          x.pos[i] = a-.5
          y.pos[i] =  obj$g.loc.center[i]
        }
        annDX = list()
        annDX$type = "point"
        annDX$idx = ann.idx 
        annDX$x.pos = x.pos
        annDX$y.pos = y.pos    
        annDX$x.right = NA
        annDX$y.bottom = NA
        #hyperlink
        if(!is.null(dim(annObj[annotation[a]][[1]][[2]]))){
          annDX$links = annObj[annotation[a]][[1]][[2]][ann.idx,]
        }else{
          annDX$links = NA
        }
        #images
        if(!is.null(dim(annObj[annotation[a]][[1]][[3]]))){
          annDX$images = annObj[annotation[a]][[1]][[3]][ann.idx,]
        }else{
          annDX$images = NA
        }
        
        eval.js(paste("annotationDX$annObjDX",a,"=annDX", sep=""))
      }else{
        eval.js(paste("annotationDX$annObjDX",a,"$type=NA", sep=""))
      }
    }
  }
  plot3 =  paste(plot3, "axis(1, at=(1:length(annotation))-.5, labels=names(annObj)[annotation], cex.axis=1, las=2)", sep=";")

  # extra plot
  plot4 = plot.call
  plot.vec =  eval.js(plot.vec)



  
 if(is.na(mai.mat[1]) & is.na(plot4)) mai.mat = matrix(c(
                                           0.2,  1.0, 1.0, 0.05,
                                           0.5,  1.0, 0.3, 0.05,
                                           0.2, 0.05, 1.0, 0.2)
                                           , byrow=TRUE, ncol=4)
  
 if(is.na(mai.mat[1]) & !is.na(plot4)) mai.mat = matrix(c(
                                           0.2,  1.0, 1.0, 0.05,
                                           0.6,  1.0, 0.5, 0.05,
                                           0.2, 0.05, 1.0, 0.2,
                                           0.2,  0.2, 1.0, 0.2)
                                           , byrow=TRUE, ncol=4)


 
  
if(is.na(plot4))  plot.calls = list(plot1=plot1, plot2=plot2, plot3=plot3)
if(!is.na(plot4))  plot.calls = list(plot1=plot1, plot2=plot2, plot3=plot3, plot4=plot4)

  
    # set environments so can access local workspace
    environment(addBounding) <- environment()
    environment(addDefault) <- environment()
    environment(automapPts) <- environment()
    environment(getBounds) <- environment()
    environment(initSplot) <- environment()
    environment(makeCharacter) <- environment()
    environment(makeImageDF) <- environment()
    environment(makeImap) <- environment()
    environment(makePolyDF) <- environment()
    environment(makeRectDF) <- environment()
    environment(makeScatterDF) <- environment()
    environment(makeSplot) <- environment()
    environment(mapMethod) <- environment()
    environment(removeImap) <- environment()
    environment(writeDefault1) <- environment()
    environment(writeDefault2) <- environment()
    environment(writeCircle.1) <- environment()
    environment(writeRect.1) <- environment()
    environment(writePoly.1) <- environment()
    environment(writeCircle.2) <- environment()
    environment(writeRect.2) <- environment()
    environment(writePoly.2) <- environment()
    environment(writeToHTML1) <- environment()
    environment(writeToHTML2) <- environment()
    
  
if(is.na(plot4))  Splot = initSplot(mat=mat, plot.calls=plot.calls, Iflag=c(TRUE,FALSE,TRUE), figTypes=c("image", "image", "image" ), mai.mat=mai.mat,mai.prc=mai.prc, image.size=image.size, plot.extras=plot.extras) 
if(!is.na(plot4))  Splot = initSplot(mat=mat, plot.calls=plot.calls, Iflag=c(TRUE,FALSE,TRUE, TRUE), figTypes=c("image", "image", "image","image"), mai.mat=mai.mat,mai.prc=mai.prc, image.size=image.size, plot.extras=plot.extras) 

  
  #######################
  #######################
  # now add interactive
  #######################
  #######################

  
  # start with what is given through function calls (i.e. x.labels, y.labels, xy.labels...)
  #   goes with figure 1 
  # ASSUMES: already subset based on TIplot parameters
  #          if larger dimension will subset based on TIplot
  
  if(!is.null(dim(x.labels))){
    if(dim(x.labels)[1] >= length(smplDX)){
       xnm = names(x.labels)
      x.labels = as.data.frame(x.labels[smplDX,])
      names(x.labels) = xnm
      
    }
    if(dim(x.labels)[1] != length(smplDX)){
      if(vrb)  warning("xlabels of not correct dimension\n continuing with x.labels as NA\n", immediate.=TRUE)
      x.labels=NA
    }
  }
  if(!is.null(dim(y.labels))){
    if(dim(y.labels)[1] > length(scanDX)){
      if(vrb) warning("ylabels of not correct dimension\n subsetting based on scanDX\n", immediate.=TRUE)
      y.labels = y.labels[scanDX,]
    }
    if(dim(y.labels)[1] != length(scanDX)){
      if(vrb) warning("ylabels of not correct dimension\ncontinuing with y.labels as NA\n", immediate.=TRUE)
      y.labels=NA
    }
  }
  if(!is.na(xy.labels[1])){
    for(i in 1:length(xy.labels)){

      if(dim(xy.labels[[i]])[1] > length(scanDX)){
        if(vrb) warning(paste("Correcting xy.labels object ",i, " y dimension\nsusetting based on scanDX\n",sep=""), immediate.=TRUE)
        xy.labels[[i]] = xy.labels[[i]][scanDX,]
      }          
      if(dim(xy.labels[[i]])[2] >= length(smplDX)){
        xy.labels[[i]] = xy.labels[[i]][,smplDX]
      }
      if((dim(xy.labels[[i]])[1] != length(scanDX)) & (dim(xy.labels[[i]])[2] != length(smplDX))){
        if(vrb) warning(paste("xy.labels object ",i, " not of correct dimension\n continuing as NA \n",sep=""), immediate.=TRUE)
        xy.labels[[i]] = NA
      }
    }
  }
  ##############
  # hyperlinks
  ##############
  if(!is.null(dim(x.links))){
    if(dim(x.links)[1] >= length(smplDX)){
      xlnm = names(x.links)
      x.links = as.data.frame(x.links[smplDX,])
      names(x.links) = xlnm
      
    }
    if(dim(x.links)[1] != length(smplDX)){
      if(vrb) warning("x.links of not correct dimension\n  continuing with x.links as NA\n", immediate.=TRUE)
      x.links=NA
    }
  }
  if(!is.null(dim(y.links))){
    if(dim(y.links)[1] > length(scanDX)){
      if(vrb) warning("y.links of not correct dimension\n subsetting based on scanDX\n", immediate.=TRUE)
      y.links = y.links[scanDX,]
    }
    if(dim(y.links)[1] != length(scanDX)){
      if(vrb) warning("y.links of not correct dimension\n continuing with y.links as NA\n", immediate.=TRUE)
      y.links=NA
    }
  }
  if(!is.na(xy.links[1])){
    for(i in 1:length(xy.links)){

      if(dim(xy.links[[i]])[1] > length(scanDX)){
        if(vrb) warning(paste("Correcting xy.links object ",i," y dimension\n susetting based on scanDX\n",sep=""), immediate.=TRUE)
        xy.links[[i]] = xy.links[[i]][scanDX,]
      }          
      if(dim(xy.links[[i]])[2] >= length(smplDX)){
        xy.links[[i]] = xy.links[[i]][,1:TIplot$vls$nsmp]
      }
      if((dim(xy.links[[i]])[1] != length(scanDX)) & (dim(xy.links[[i]])[2] != length(smplDX))){
        if(vrb) warning(paste("xy.links object ",i, " not of correct dimension\n continuing as NA \n",sep=""), immediate.=TRUE)
        xy.links[[i]] = NA
      }
    }
  }
  if(length(asLinks) == 1){
    if(!is.na(asLinks)) asLinks = rep(asLinks, (TIplot$vls$nsmp*length(scanDX)))
  }else{
    if(!is.null(dim(asLinks))){
      if(dim(asLinks)[1] > length(scanDX)) asLinks = asLinks[scanDX,]
      if(dim(asLinks)[2] > length(smplDX)) asLinks = asLinks[,smplDX]
      if((dim(asLinks)[1] != length(scanDX)) & (dim(asLinks)[2] != length(smplDX))){
       if(vrb) warning("asLinks of not correct dimension\n continuing with asLinks as NA\n", immediate.=TRUE)
        asLinks=NA
      }
    }
  }
  NewXlinks = x.links
  NewYlinks = y.links
  
  ##############
  # images
  ##############
  if(!is.null(dim(x.images))){
    if(dim(x.images)[1] >= length(smplDX)){
      xinm = names(x.images)
      x.images = as.data.frame(x.images[smplDX,])
      names(x.images) = xinm
    }
    if(dim(x.images)[1] != length(smplDX)){
      if(vrb) warning("x.images of not correct dimension\n continuing with x.images as NA\n", immediate.=TRUE)
      x.images=NA
    }
  }
  if(!is.null(dim(y.images))){
    if(dim(y.images)[1] > length(scanDX)){
      if(vrb) warning("y.images of not correct dimension\n subsetting based on scanDX\n", immediate.=TRUE)
      y.images = y.images[scanDX,]
    }
    if(dim(y.images)[1] != length(scanDX)){
      if(vrb) warning("y.images of not correct dimension\n continuing with y.images as NA\n", immediate.=TRUE)
      y.images=NA
    }
  }
  if(!is.na(xy.images[1])){
    for(i in 1:length(xy.images)){

      if(dim(xy.images[[i]])[1] > length(scanDX)){
        if(vrb) warning(paste("Correcting xy.images object ",i, " y dimension\n susetting based on scanDX\n",sep=""), immediate.=TRUE)
        xy.images[[i]] = xy.images[[i]][scanDX,]
      }          
      if(dim(xy.images[[i]])[2] >= length(smplDX)){
        xy.image[[i]] = xy.images[[i]][,1:TIplot$vls$nsmp]
      }
      if((dim(xy.images[[i]])[1] != length(scanDX)) & (dim(xy.images[[i]])[2] != length(smplDX))){
        if(vrb) warning(paste("xy.image object ",i, " not of correct dimension\n continuing as NA \n",sep=""), immediate.=TRUE)
        xy.images[[i]] = NA
      }
    }
  }

  NewXimages=x.images
  NewYimages=y.images


  

  #
  # now also grab columns from set mapping obj
  #
  if(is.na(mapObj.columns[1])) mapObj.columns = 1:dim(mapObj$mapping.info)[2]
  if(mapObj.columns[1] == 0){
    y.labels = y.labels
  }else{
    if(class(mapObj.columns)=="character"){
      mapObj.columns = match(mapObj.columns, names(mapObj$mapping.info))
      idx = which(is.na(mapObj.columns))
      if(length(idx) != 0) mapObj.columns = mapObj.columns[-idx]
    }
    if((class(mapObj.columns) == "numeric" | class(mapObj.columns) == "integer") & (length(mapObj.columns) !=0)){
      map.lbls = mapObj$mapping.info[scanDX, mapObj.columns]
    }else{
      map.lbls = NA
    }

    if(is.null(dim(map.lbls)) & is.null(dim(y.labels))){
      y.labels = NA
    }else{
      if(is.null(dim(y.labels))){
        y.labels = map.lbls
      }
      if(!is.null(dim(map.lbls)) & !is.null(dim(y.labels))){
        y.labels = cbind(y.labels, map.lbls)
      }
    }
  }

  # add links is in mapping.obj
  # NOTE: as of right now...if links was set in mapping it will show up in interactive
  #       to control automatic display must be changed in mapping object

  if(!is.null(dim(mapObj$links))){
    if(is.null(dim(y.links))){
      y.links = as.data.frame(mapObj$links[scanDX,])
      names(y.links) = names(mapObj$links)
    }else{
      orig.names = names(y.links)
      y.links = cbind(y.links, mapObj$links[scanDX,])
      names(y.links) = c(orig.names, names(mapObj$links))
    }    
  }

  # add images is in mapping.obj
  # NOTE: as of right now...if images was set in mapping it will show up in interactive
  #       to control automatic display must be changed in mapping object

  if(!is.null(dim(mapObj$images))){
    if(is.null(dim(y.images))){
      y.images = as.data.frame(mapObj$images[scanDX,])
      names(y.images) = names(mapObj$images)
    }else{
      orig.names = names(y.images)
      y.images = cbind(y.images, mapObj$images[scanDX,])
      names(y.images) = c(orig.names, names(mapObj$images))
    }    
  }


  

  

  ############################
  #
  # rework interactive boxes
  #
  ############################
  
  #x.pos = 1:TIplot$vls$nsmp
  #y.pos = scanLoc
  
  #Splot = makeImap(Splot, figure=1, xy.type="image.midpoints", x.pos=x.pos, y.pos=y.pos, x.labels=x.labels, y.labels=y.labels, xy.labels=xy.labels, x.links=x.links, y.links=y.links, xy.links=xy.links, asLinks=asLinks, fname.root=fname.root, dir=dir, ...)
  


  # adjust x.labels, y.labels, etc. because using boxes
  if(!is.null(dim(x.labels))){
    x.labels$tempOrd = 1:dim(x.labels)[1]
    NewXlabels = rbind(x.labels, x.labels)
    times = length(scanDX)-2
    while(times > 0){
      NewXlabels = rbind(NewXlabels, x.labels)
      times=times-1
    }
    NewXlabels = NewXlabels[order(NewXlabels$tempOrd),1:(dim(NewXlabels)[2]-1)]
    x.labels = x.labels[,1:(dim(x.labels)[2]-1)]
  }
  if(!is.null(dim(x.links))){
    x.links$tempOrd = 1:dim(x.links)[1]
    NewXlinks = rbind(x.links, x.links)
    times = length(scanDX)-2
    while(times > 0){
      NewXlinks = rbind(NewXlinks, x.links)
      times=times-1
    }
    NewXlinks = NewXlinks[order(NewXlinks$tempOrd),1:(dim(NewXlinks)[2]-1)]
  }

  if(!is.null(dim(y.labels))){
    NewYlabels = rbind(y.labels, y.labels)
    times = TIplot$vls$nsmp-2 
    while(times > 0){
      NewYlabels = rbind(NewYlabels, y.labels)
      times = times-1
    }    
  }

  if(!is.null(dim(y.links))){
    NewYlinks = rbind(y.links, y.links)
    times = TIplot$vls$nsmp-2 
    while(times > 0){
      NewYlinks = rbind(NewYlinks, y.links)
      times = times-1
    }
  }

  nm = names(NewYlabels)
  if(!is.na(xy.labels[1])){
    for(i in 1:length(xy.labels)){
      NewYlabels = cbind(NewYlabels, as.vector(xy.labels[[i]]))
    }
    names(NewYlabels) = c(nm, names(xy.labels))
    xy.labels = NA
  }
 
  nm = names(NewYlinks)
  if(!is.na(xy.links[1])){
    for(i in 1:length(xy.links)){
      NewYlinks = cbind(NewYlinks, as.vector(xy.links[[i]]))
    }
    names(NewYlinks) = c(nm, names(xy.links))
    xy.links = NA
  }
  
  if(is.null(dim(asLinks))){
    if(length(asLinks) == TIplot$vls$nsmp) asLinks = rep(asLinks, each=length(scanDX))
    if(length(asLinks) == length(scanDX)) asLinks = rep(asLinks, times=TIplot$vls$nsmp)
  }else{
    asLinks = as.vector(asLinks)
  }
  

  
  
  if(!is.null(dim(x.images))){
    x.images$tempOrd = 1:dim(x.images)[1]
    NewXimages = rbind(x.images, x.images)
    times = length(scanDX)-2
    while(times > 0){
      NewXimages = rbind(NewXimages, x.images)
      times=times-1
    }
    NewXimages = NewXimages[order(NewXimages$tempOrd),1:(dim(NewXimages)[2]-1)]
  }

  if(!is.null(dim(y.images))){
    NewYimages = rbind(y.images, y.images)
    times = TIplot$vls$nsmp-2 
    while(times > 0){
      NewYimages = rbind(NewYimages, y.images)
      times = times-1
    }
  }
  nm = names(NewYimages)
  if(!is.na(xy.images[1])){
    for(i in 1:length(xy.images)){
      NewYimages = cbind(NewYimages, as.vector(xy.images[[i]]))
    }
    names(NewYimages) = c(nm, names(xy.images))
    xy.images = NA
  }









  

  # find each box boundary
  ttracks = (TIplot$vls$nsmp*TIplot$vls$H)
  tord = rep(1:TIplot$vls$H, ceiling(length(scanDX)/TIplot$vls$H))[1:length(scanDX)]
  xl = numeric(length(0))
  xr = numeric(length(0))
  yb = numeric(length(0))
  yt = numeric(length(0))
  
  for(s in 1:TIplot$vls$nsmp){
    for(i in 1:length(tord)){
      xl = c(xl,TIplot$vls$Xcoords[TIplot$tractBound[[tord[i]]]$xdx[s]] -.25)
      xr = c(xr,TIplot$vls$Xcoords[TIplot$tractBound[[tord[i]]]$xdx[s]] +.25)
      yb = c(yb,TIplot$vls$Y[i,1])
      yt = c(yt,TIplot$vls$Y[i,2])
    }
  }
  xl = xl[-1]
  xr = xr[-1]
  yb = yb[-1]
  yt = yt[-1]
  
  Splot = makeImap(Splot, figure=1, xy.type="rect",
                   x.pos=xl, y.pos=yt,
                   x.right.pos=xr, y.bottom.pos=yb,
                   x.labels=NewXlabels, y.labels=NewYlabels, xy.labels=xy.labels, x.links=NewXlinks, y.links=NewYlinks, xy.links=xy.links, asLinks=asLinks, fname.root=fname.root, xy.images=xy.images, x.images=NewXimages, y.images=NewYimages, dir=dir, cleanDir=cleanDir, ...)
  




  
  #
  # now add interactive data for annotation sets
  #

  for(i in (1:length(annotation))){
    if(!is.na(annotationDX[[i]]$type)){
      obj = annObj[[annotation[i]]]
      display.names = c("Label","Chrom", "loc.start", "loc.stop","loc.center")
      display = match( display.names, names(obj[[1]]))
      idx = which(is.na(display))
      if(length(idx) != 0){
        diplay = display[-idx]
        display.names = display.names[-idx]
      }
      exclude = match(c("g.loc.start", "g.loc.center", "g.loc.stop"), names(obj[[1]]))
      extra.dx = setdiff(c(exclude, display), 1:length(names(obj[[1]])))
      y.labels2 = obj[[1]][annotationDX[[i]]$idx,c(display,extra.dx)]
      # hyperlinks
      if(!is.null(dim(annotationDX[[i]]$links))){
        y.links2 =  annotationDX[[i]]$links
      }else{
        y.links2 = NA
      }
      # images
      if(!is.null(dim(annotationDX[[i]]$images))){
        y.images2 =  annotationDX[[i]]$images
      }else{
        y.images2 = NA
      }




      
      Splot = makeImap(Splot, figure=3, xy.type=annotationDX[[i]]$type, x.pos = annotationDX[[i]]$x.pos, y.pos = annotationDX[[i]]$y.pos,  x.right.pos = annotationDX[[i]]$x.right, y.bottom.pos = annotationDX[[i]]$y.bottom, y.labels = y.labels2, y.links=y.links2, y.images=y.images2, fname.root=fname.root, dir=dir,cleanDir=cleanDir, ...)
    }
  }



  #
  # now add interactive for extra plot if applicable
  #

  if(!is.na(plot4)){
    if(length(plot.vec) > length(scanDX)){
      y.vec.pos = rep(scanLoc, (length(plot.vec)/length(scanDX)))      
    }else{
      y.vec.pos = scanLoc
    }
    new.y.labels = y.labels
    new.links = y.links
    new.images=y.images
    if(length(plot.vec) > length(scanDX)){
      for(i in 1:((length(plot.vec)/length(scanDX))-1)){
        new.y.labels = rbind(new.y.labels, y.labels)
        new.links = rbind(new.links, y.links)
        new.images = rbind(new.images, y.images)
      }
    }
    
    nm = names(new.y.labels)
    new.y.labels = cbind(plot.vec, new.y.labels)
    names(new.y.labels) = c("plotValue", nm)

    Splot = makeImap(Splot, figure=4, xy.type="points", x.pos = plot.vec, y.pos = y.vec.pos, y.labels = new.y.labels, y.links=new.links, y.images=new.images, fname.root=fname.root, dir=dir, cleanDir=cleanDir,...)
 
  }


  Splot = makeSplot(Splot, fname.root=fname.root, dir=dir,overwriteSourcePlot =overwriteSourcePlot,makeInteractive=makeInteractive,overrideInteractive=overrideInteractive,header=header,window.size=window.size,returnObj=TRUE)

  
}


