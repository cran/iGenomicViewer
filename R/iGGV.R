

##########################################
#
# bioinforamtics wrapper to sendplot
#
###########################################


#
# this add side plot
#

iGGV <- function(vls,                  
                  mapObj,
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
                  maxLabels=25,
                  mai.mat = NA,
                  mai.prc=FALSE,
                  plot.x.index=NA,
                  smp.color = NA,

                  plot.y.index=NA,

                  goodDX=NA,
                 
                  genomic.start=NA,
                  genomic.stop=NA,

                  genomic.region=NA,   
                  region.type="chrom",

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
                  ... # extra arguments to makeImap not in above 
                  ){


  if(is.null(dim(mat))){
    if(is.na(plot.call)) mat = matrix(c(rep(c(rep(1,7), rep(3,1)), 8), rep(c(rep(2,7), rep(0,1)), 2)), byrow=TRUE, ncol=8)
    if(!is.na(plot.call)) mat = matrix(c(rep(c(rep(1,7), rep(3,1), rep(4,2)), 9), rep(c(rep(2,7), rep(0,3)), 2)), byrow=TRUE, ncol=10)
  }else{
    mat = mat
  }


  # set goodDX default
  if(is.na(goodDX[1])) goodDX = 1:dim(mapObj$mapping.info)[1]
  
  #
  # x axis -- samples 
  #  
  if(!is.na(plot.x.index[1])){
    smpls = plot.x.index
    if(is.null(colnames(vls))){
      x.names = paste("Smp", smpls, sep="")
    }else{
      x.names = colnames(vls)[plot.x.index]
    }    
  }else{
    smpls = 1:dim(vls)[2]
    if(is.null(colnames(vls))){
      x.names = paste("Smp", smpls, sep="")
    }else{
      x.names = colnames(vls)
    }        
  }
  nsmpls = length(smpls)
  if(is.na(smp.color[1])){
    smp.color = rep("black", nsmpls)
  }
  if(length(smp.color) == (dim(vls)[2])) smp.color = smp.color[smpls]
  if(length(smp.color) != length(smpls)){
    warning("Sample colors not of correct dimension\n Continuing with no color\n", immediate.=TRUE)
    smp.color = rep("black", nsmpls)
  }

  
  ###################
  # y axis 
  ##################
  
  #
  # by what method are they giving y-axis
  #

  if(dim(vls)[1] != dim(mapObj$mapping.info)[1]){
    stop("mapObj does not appear to correspond to value matrix \n")
  }
  
  if(!is.na(genomic.region) & ( (region.type != "chrom") & (region.type != "arm") & (region.type !="broad.band") & (region.type != "fine.band"))){
    stop("Error: region type not specified\n")
  }
  
  # if region is given 
  if(!is.na(genomic.region) & region.type=="chrom"){
    genomic.start = mapObj$band.info$Chrom$lower[which(mapObj$band.info$Chrom$Label == genomic.region)]
    genomic.stop = mapObj$band.info$Chrom$upper[which(mapObj$band.info$Chrom$Label == genomic.region)]
  }
  if(!is.na(genomic.region) & region.type=="arm"){
    genomic.start = mapObj$band.info$Arm$lower[which(mapObj$band.info$Arm$Arm == genomic.region)]
    genomic.stop = mapObj$band.info$Arm$upper[which(mapObj$band.info$Arm$Arm == genomic.region)]
  }
  if(!is.na(genomic.region) & region.type=="broad.band"){
    genomic.start = mapObj$band.info$Broad.Band$lower[which(mapObj$band.info$Broad.Band$Broad.Band == genomic.region)]
    genomic.stop = mapObj$band.info$Broad.Band$upper[which(mapObj$band.info$Broad.Band$Broad.Band == genomic.region)]
  }
  if(!is.na(genomic.region) & region.type=="fine.band"){
    genomic.start = mapObj$band.info$Fine.Band$lower[which(mapObj$band.info$Fine.Band$Fine.Band == genomic.region)]
    genomic.stop = mapObj$band.info$Fine.Band$upper[which(mapObj$band.info$Fine.Band$Fine.Band == genomic.region)]
  }
  # if genomic start and stop are given 
  if(!is.na(genomic.start) & !is.na(genomic.stop)){
    idx.start = which(mapObj$mapping.info$g.loc.start >= genomic.start)[1]
    idx.stop = which(mapObj$mapping.info$g.loc.stop <= genomic.stop)[length(which(mapObj$mapping.info$g.loc.stop <= genomic.stop))]
    y.index = intersect(idx.start:idx.stop, goodDX)
    y.names = as.character(mapObj$mapping.info$Spot.ID)[y.index]
    scanDX = y.index
  }
  # if index is given
  if(!is.na(plot.y.index[1])){
    plot.y.index = intersect(plot.y.index, goodDX)
    y.names = as.character(mapObj$mapping.info$Spot.ID)[plot.y.index]
    genomic.start = min(mapObj$mapping.info$g.loc.start[plot.y.index],na.rm=TRUE)
    genomic.stop = max(mapObj$mapping.info$g.loc.stop[plot.y.index],na.rm=TRUE)
    scanDX = plot.y.index
  }
  
  # subset vls
  z = vls[scanDX, smpls]
  zlgnd=array(seq(from=lims[1],to=lims[2],length=1000),dim=c(1,1000))
  z[z>lims[2]] = lims[2]
  z[z<lims[1]] = lims[1]

  #
  # For now go off spot.ID genomic locations for plotting
  #
  # eventually rework or extend this section for genomic location
  #    ?? does this make sense? add in breaks for regions not covered ?
  
  scanLoc = mapObj$mapping.info$g.loc.center[scanDX]
  while(length(which(diff(scanLoc)<=0)) > 0){
    scanLoc[which(diff(scanLoc)<=0)]=scanLoc[which(diff(scanLoc)<=0)]-1
  }
  
  # heatmap
  plot1 = "image(x=1:nsmpls,y=scanLoc,z=t(z),zlim=c(min(zlgnd,na.rm=TRUE),max(zlgnd, na.rm=TRUE)),ylim=range(scanLoc,na.rm=TRUE),col=c(hsv(h=2/6,v=seq(1,0,length=1000)^1.15), hsv(h=0/6,v=seq(0,1,length=1000)^1.15)),axes=FALSE,xlab='',ylab='')"

#  count.arm=sum((mapObj$band.info$Arm$upper>=min(scanLoc,na.rm=TRUE))
#    &(mapObj$band.info$Arm$lower<=max(scanLoc,na.rm=TRUE)))  
#  count.broadband=sum((mapObj$band.info$Broad.Band$upper>=min(scanLoc,na.rm=TRUE))
#    &(mapObj$band.info$Broad.Band$lower<=max(scanLoc,na.rm=TRUE)))
#  count.finband=sum((mapObj$band.info$Fine.Band$upper>=min(scanLoc,na.rm=TRUE))
#    &(mapObj$band.info$Fine.Band$lower<=max(scanLoc,na.rm=TRUE)))

#  if((count.arm == count.broadband) & (count.broadband == count.finband)){
#    ilbl = 4
#  }else{
#    ilbl=order(abs(c(count.arm,count.broadband,count.finband,length(scanDX))-maxLabels))[1]
#  }
#  if(ilbl<=3){
#    if(ilbl==1) bandDX=1:40
#    if(ilbl==2) bandDX=((sum(mapObj$band.info$Broad.Band$upper<=min(scanLoc,na.rm=TRUE),na.rm=TRUE)+1)
#         :(sum(mapObj$band.info$Broad.Band$lower<=max(scanLoc,na.rm=TRUE),na.rm=TRUE)))
#    if(ilbl==3) bandDX=((sum(mapObj$band.info$Fine.Band$upper<=min(scanLoc,na.rm=TRUE),na.rm=TRUE)+1)
#         :(sum(mapObj$band.info$Fine.Band$lower<=max(scanLoc,na.rm=TRUE),na.rm=TRUE)))
#    lbls=as.character(mapObj$band.info[ilbl+2][[1]][bandDX,1])
#    lbls.loc =mapObj$band.info[ilbl+2][[1]]$center[bandDX]
#    if(lbls.loc[1] < range(scanLoc, na.rm=TRUE)[1]) lbls.loc[1] = mean(c(range(scanLoc, na.rm=TRUE)[1],mapObj$band.info[ilbl+2][[1]]$lower[bandDX][2]))
#    if(lbls.loc[length(lbls.loc)] > range(scanLoc, na.rm=TRUE)[2]) lbls.loc[length(lbls.loc)] = mean(c(range(scanLoc, na.rm=TRUE)[2],mapObj$band.info[ilbl+2][[1]]$lower[bandDX][length(bandDX)]))
#    plot1 = paste(plot1,"axis(2,lbls.loc,tick=FALSE,labels=lbls,las=2,cex.axis=1)","axis(2,mapObj$band.info[ilbl+2][[1]]$lower[bandDX], labels=FALSE)", sep=";")
#  }
#  if(ilbl==4){   
#    lbls=as.character(mapObj$mapping.info$Spot.ID[scanDX])
#    plot1 = paste(plot1,"axis(2,mapObj$mapping.info$g.loc.center[scanDX],tick=FALSE,labels=lbls, las=2,cex.axis=1)", sep=";")
#  }

        # common reference from other programs 
        spotDX = scanDX

        gloc.start = mapObj$mapping.info$g.loc.start[spotDX[1]]
        gloc.stop = mapObj$mapping.info$g.loc.stop[spotDX[length(spotDX)]]
           
        # figure out labels
        idx.Chrom = intersect(which(mapObj$band.info$Chrom$lower <= gloc.stop), which(mapObj$band.info$Chrom$upper >= gloc.start))
        if(gloc.start < mapObj$band.info$Chrom$lower[idx.Chrom[1]]) idx.Chrom = c(idx.Chrom[1]-1, idx.Chrom)
        if(gloc.stop > mapObj$band.info$Chrom$upper[idx.Chrom[length(idx.Chrom)]]) idx.Chrom = c(idx.Chrom, idx.Chrom[length(idx.Chrom)]+1)
        idx.Arm = intersect(which(mapObj$band.info$Arm$lower <= gloc.stop), which(mapObj$band.info$Arm$upper >= gloc.start))
        if(gloc.start < mapObj$band.info$Arm$lower[idx.Arm[1]]) idx.Arm = c(idx.Arm[1]-1, idx.Arm)
        if(gloc.stop > mapObj$band.info$Arm$upper[idx.Arm[length(idx.Arm)]]) idx.Arm = c(idx.Arm, idx.Arm[length(idx.Arm)]+1)
        idx.Broad.Band = intersect(which(mapObj$band.info$Broad.Band$lower <= gloc.stop), which(mapObj$band.info$Broad.Band$upper >= gloc.start))
        if(gloc.start < mapObj$band.info$Broad.Band$lower[idx.Broad.Band[1]]) idx.Broad.Band = c(idx.Broad.Band[1]-1, idx.Broad.Band)
        if(gloc.stop > mapObj$band.info$Broad.Band$upper[idx.Broad.Band[length(idx.Broad.Band)]]) idx.Broad.Band = c(idx.Broad.Band, idx.Broad.Band[length(idx.Chrom)]+1)
        idx.Fine.Band = intersect(which(mapObj$band.info$Fine.Band$lower <= gloc.stop), which(mapObj$band.info$Fine.Band$upper >= gloc.start))
        if(gloc.start < mapObj$band.info$Fine.Band$lower[idx.Fine.Band[1]]) idx.Fine.Band = c(idx.Fine.Band[1]-1, idx.Fine.Band)
        if(gloc.stop > mapObj$band.info$Fine.Band$upper[idx.Fine.Band[length(idx.Fine.Band)]]) idx.Fine.Band = c(idx.Fine.Band, idx.Fine.Band[length(idx.Chrom)]+1)
        
        idx.spots = spotDX
        
        nChrom = length(idx.Chrom)
        nArm =  length(idx.Arm)
        nBroad.Band =  length(idx.Broad.Band)
        nFine.Band =  length(idx.Fine.Band)
        nSpots = length(spotDX)


 
            # determine which set of labels to use
            Spots.Used = FALSE
            Fine.Band.Used = FALSE
            Broad.Band.Used=FALSE
            Arm.Used=FALSE
            Chrom.Used=TRUE
            
            if(nSpots <= maxLabels){
              Spots.Used = TRUE
              Chrom.Used = FALSE
           
              # get section divides
              lwrDX = idx.spots[1]
              uprDX = idx.spots[length(idx.spots)]
              sectiondiv = mapObj$mapping.info$g.loc.start[lwrDX]
              for(i in lwrDX:(uprDX-1)){
                sectiondiv = c(sectiondiv, mean(c(mapObj$mapping.info$g.loc.stop[i], mapObj$mapping.info$g.loc.start[i+1])))
              }
              sectiondiv = c(sectiondiv, mapObj$mapping.info$g.loc.stop[uprDX])
              # get section centers
              sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              for(i in 2:(length(sectiondiv)-1)){
                sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
              }

              # back check if only on location
              if(lwrDX==uprDX){
                sectiondiv = c(mapObj$mapping.info$g.loc.start[lwrDX],mapObj$mapping.info$g.loc.stop[uprDX]) 
                sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              }                    
              # get sectio labels
              sectionlbls = as.character(mapObj$mapping.info$Spot.ID)[idx.spots]
              
            }else{
              if(nFine.Band <=maxLabels){
                Fine.Band.Used=TRUE
                Chrom.Used=FALSE
        
                # get section divides
                lwrDX = idx.Fine.Band[1]
                uprDX = idx.Fine.Band[length(idx.Fine.Band)]
                sectiondiv = mapObj$band.info$Fine.Band$lower[lwrDX]
                for(i in lwrDX:(uprDX-1)){
                  sectiondiv = c(sectiondiv, mean(c(mapObj$band.info$Fine.Band$upper[i], mapObj$band.info$Fine.Band$lower[i+1])))
                }
                sectiondiv = c(sectiondiv, mapObj$band.info$Fine.Band$upper[uprDX])
                # get section centers
                sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                for(i in 2:(length(sectiondiv)-1)){
                  sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
                }

                # back check if only one location
                if(lwrDX == uprDX){
                  sectiondiv = c(mapObj$band.info$Fine.Band$lower[lwrDX], mapObj$band.info$Fine.Band$upper[uprDX])
                  sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                }
                
                # get sectio labels
                sectionlbls = as.character(mapObj$band.info$Fine.Band$label)[idx.Fine.Band]
                
              }else{
                if(nBroad.Band<=maxLabels){
                  Broad.Band.Used=TRUE
                  Chrom.Used=FALSE
        
                  # get section divides
                  lwrDX = idx.Broad.Band[1]
                  uprDX = idx.Broad.Band[length(idx.Broad.Band)]
                  sectiondiv = mapObj$band.info$Broad.Band$lower[lwrDX]
                  for(i in lwrDX:(uprDX-1)){
                    sectiondiv = c(sectiondiv, mean(c(mapObj$band.info$Broad.Band$upper[i], mapObj$band.info$Broad.Band$lower[i+1])))
                  }
                  sectiondiv = c(sectiondiv, mapObj$band.info$Broad.Band$upper[uprDX])
                  # get section centers
                  sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                  for(i in 2:(length(sectiondiv)-1)){
                    sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
                  }
    
                  # back check if only on location
                  if(lwrDX==uprDX){
                    sectiondiv = c(mapObj$band.info$Broad.Band$lower[lwrDX], mapObj$band.info$Broad.Band$upper[uprDX])
                    sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))            
                  }
          
                  # get sectio labels
                  sectionlbls = as.character(mapObj$band.info$Broad.Band$label)[idx.Broad.Band]
                  
                }else{
                  if(nArm<=maxLabels){
                    Arm.Used=TRUE
                    Chrom.Used=FALSE
                
                    # get section divides
                    lwrDX = idx.Arm[1]
                    uprDX = idx.Arm[length(idx.Arm)]
                    sectiondiv = mapObj$band.info$Arm$lower[lwrDX]
                    for(i in lwrDX:(uprDX-1)){
                      sectiondiv = c(sectiondiv, mean(c(mapObj$band.info$Arm$upper[i], mapObj$band.info$Arm$lower[i+1])))
                    }
                    sectiondiv = c(sectiondiv, mapObj$band.info$Arm$upper[uprDX])
                    # get section centers
                    sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                    for(i in 2:(length(sectiondiv)-1)){
                      sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
                    }
                            
                    # back check if only one location
                    if(uprDX==lwrDX){
                      sectiondiv = c(mapObj$band.info$Arm$lower[lwrDX],mapObj$band.info$Arm$upper[uprDX])
                      sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                    }            
                    # get section labels
                    sectionlbls = as.character(mapObj$band.info$Arm$label)[idx.Arm]
                    
                  }
                }
              }
            }
            if(Chrom.Used){
          
              # get section divides
              lwrDX = idx.Chrom[1]
              uprDX = idx.Chrom[length(idx.Chrom)]
              sectiondiv = mapObj$band.info$Chrom$lower[lwrDX]
              for(i in lwrDX:(uprDX-1)){
                sectiondiv = c(sectiondiv, mean(c(mapObj$band.info$Chrom$upper[i], mapObj$band.info$Chrom$lower[i+1])))
              }
              sectiondiv = c(sectiondiv, mapObj$band.info$Chrom$upper[uprDX])
              # get section centers
              sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              for(i in 2:(length(sectiondiv)-1)){
                sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
              }
  
              # back check if only one location
              if(uprDX==lwrDX){
                sectiondiv = c(mapObj$band.info$Chrom$lower[lwrDX], mapObj$band.info$Chrom$upper[uprDX])
                sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              }                
              # get sectio labels
              sectionlbls = as.character(mapObj$band.info$Chrom$label)[idx.Chrom]
            }

            ylim.temp = range(mapObj$mapping.info$g.loc.center[scanDX],na.rm=TRUE)
            
            # fix sectiondiv and sectioncntr
            if(sectiondiv[1] < ylim.temp[1]){
              sectiondiv[1] = ylim.temp[1]
              sectioncntr[1] = mean(c(sectiondiv[1], sectiondiv[2]))
            }
            if(sectiondiv[length(sectiondiv)] > ylim.temp[2]){
              sectiondiv[length(sectiondiv)] = ylim.temp[2]
              sectioncntr[length(sectioncntr)] = mean(c(sectiondiv[length(sectiondiv)], sectiondiv[(length(sectiondiv)-1)]))
            }


  plot1 = paste(plot1,"axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=1,las=2)", sep=";")







  



  


  
  if(smpLines){
    plot1 = paste(plot1, paste("abline(v=(0:length(smpls))+1/2,col='",divCol,"',lty=1,lwd=1/3)",sep=""), sep=";")
  }
    
  # add sample lables
  for(c in 1:length(unique(smp.color))){
    iclr = unique(smp.color)[c]
    nm = paste("cdx", c, sep="")
    eval.js(paste(nm, "=which(smp.color == iclr)",sep=""))
    plot1 = paste(plot1, paste("axis(3, at=",nm, ",labels=x.names[",nm,"], cex.axis=1, las=2, col.axis='", iclr,"')",sep=""), sep=";")

  }
 
  # legend
  plot2="image(x=as.vector(zlgnd),y=1,z=t(zlgnd),zlim=c(min(zlgnd,na.rm=TRUE),max(zlgnd, na.rm=TRUE)),col=c(hsv(h=2/6,v=seq(1,0,length=length(zlgnd))^1.15),hsv(h=0/6,v=seq(0,1,length=length(zlgnd))^1.15)), axes=FALSE,xlab='',ylab='')"
  plot2 = paste(plot2, paste("axis(1,seq(from=",min(zlgnd),",to=",max(zlgnd),",length=5),line=0)", sep=""), sep=";")
  
  # annotation
  if(is.na(annotation[1])) annotation = 1:length(annObj)

  plot3 = "image(x=0:1,y=0:1,z=matrix(rep(NA,4),ncol=2),xlim=c(0,length(annotation)),ylim=range(scanLoc,na.rm=TRUE),zlim=c(0,1),axes=FALSE,xlab='',ylab='')"  

  annotationDX = list()
  
  for(a in 1:length(annotation)){
    obj = annObj[annotation[a]][[1]][[1]]
    if(!is.null(obj$g.loc.start)){
      #ann.idx = intersect(which(obj$g.loc.start > range(scanLoc,na.rm=TRUE)[1]), which(obj$g.loc.stop < range(scanLoc,na.rm=TRUE)[2]))
      ann.idx = intersect(which(obj$g.loc.start <= range(scanLoc,na.rm=TRUE)[2]), which(obj$g.loc.stop >= range(scanLoc,na.rm=TRUE)[1]))

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
        # hyperlinks
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

        # hyperlinks
        if(!is.null(dim(annObj[annotation[a]][[1]][[2]]))){
          annDX$links = annObj[annotation[a]][[1]][[2]][ann.idx,]
        }else{
          annDX$links = NA
        }
        # images
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
                                           1.3,  1.0, 0.3, 0.05,
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
  # we assume complete tables that will be subset based on x and y
  
  if(!is.null(dim(x.labels))){
    xnm = names(x.labels)
    x.labels = as.data.frame(x.labels[smpls,])
    names(x.labels) = xnm
  }
  if(!is.null(dim(y.labels))){
    y.labels = y.labels[scanDX,]
  }
  if(!is.na(xy.labels[1])){
    for(i in 1:length(xy.labels)){
      xy.labels[[i]] = xy.labels[[i]][scanDX,smpls]    
    }
  }
  # links
  if(!is.null(dim(x.links))){
    xlnm = names(x.links)
    x.links = as.data.frame(x.links[smpls,])
    names(x.links) = xlnm
    
  }
  if(!is.null(dim(y.links))){
    y.links = y.links[scanDX,]
  }
  if(!is.na(xy.links[1])){
    for(i in 1:length(xy.links)){
      xy.links[[i]] = xy.links[[i]][scanDX,smpls]    
    }
  }
  if(length(asLinks) == 1){
    if(!is.na(asLinks)) asLinks = rep(asLinks, (nsmpls*length(scanDX)))
  }else{
    if(!is.null(dim(asLinks))) asLinks = asLinks[scanDX,smpls]
  }
  # images
  if(!is.null(dim(x.images))){
    xinm = names(x.images)
    x.images = as.data.frame(x.images[smpls,])
    names(x.images) = xinm
  }
  if(!is.null(dim(y.images))){
    y.images = y.images[scanDX,]
  }
  if(!is.na(xy.images[1])){
    for(i in 1:length(xy.images)){
      xy.images[[i]] = xy.images[[i]][scanDX,smpls]    
    }
  }
 





  
 
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
      nms1 = names(mapObj$mapping.info)[mapObj.columns]
      map.lbls = as.data.frame(mapObj$mapping.info[scanDX, mapObj.columns])
      names(map.lbls) = nms1  
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
  
  x.pos = 1:nsmpls
  #   x.pos = c(.5, ((1:nsmpls) + .5))
  y.pos = scanLoc
  
  Splot = makeImap(Splot, figure=1, xy.type="image.midpoints", x.pos=x.pos, y.pos=y.pos, x.labels=x.labels, y.labels=y.labels, xy.labels=xy.labels, x.links=x.links, y.links=y.links, xy.links=xy.links, asLinks=asLinks, x.images=x.images, y.images=y.images, xy.images=xy.images, fname.root=fname.root, dir=dir, cleanDir=cleanDir, ...)
  
  
  #
  # now add interactive data for annotation sets
  #

  for(i in (1:length(annotation))){
    if(!is.na(annotationDX[[i]]$type)){
      obj = annObj[[annotation[i]]]
      display.names = c("Label","Chrom", "loc.start", "loc.stop")
      display = match( display.names, names(obj[[1]]))
      idx = which(is.na(display))
      if(length(idx) != 0){
        display = display[-idx]
        display.names = display.names[-idx]
      }
      exclude = match(c("g.loc.start", "g.loc.stop"), names(obj[[1]]))
      extra.dx = setdiff(c(exclude, display), 1:length(names(obj[[1]])))

      y.labels2 = obj[[1]][annotationDX[[i]]$idx,c(display,extra.dx)]

      # hyperlinks
      if(!is.null(dim(annotationDX[[i]]$links))){
        y.links2 =  annotationDX[[i]]$links
      }else{
        y.links2 = NA
      }
      #images
      if(!is.null(dim(annotationDX[[i]]$images))){
        y.images2 =  annotationDX[[i]]$images
      }else{
        y.images2 = NA
      }

      
      Splot = makeImap(Splot, figure=3, xy.type=annotationDX[[i]]$type, x.pos = annotationDX[[i]]$x.pos, y.pos = annotationDX[[i]]$y.pos,  x.right.pos = annotationDX[[i]]$x.right, y.bottom.pos = annotationDX[[i]]$y.bottom, y.labels = y.labels2, y.links=y.links2,y.images=y.images2, fname.root=fname.root, dir=dir,cleanDir=cleanDir, ...)
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
        if(!is.null(dim(new.y.labels))){
          new.y.labels = rbind(new.y.labels, y.labels)
        }else{
          new.y.labels=NA
        }
        if(!is.null(dim(new.links))){
          new.links = rbind(new.links, y.links)
        }else{
          new.links=NA
        }
        if(!is.null(dim(new.images))){
          new.images=rbind(new.images, y.images)
        }else{
          new.images=NA
        }
      }
    }
    nm = names(new.y.labels)
    new.y.labels = cbind(plot.vec, new.y.labels)
    names(new.y.labels) = c("plotValue", nm)

    Splot = makeImap(Splot, figure=4, xy.type="points", x.pos = plot.vec, y.pos = y.vec.pos, y.labels = new.y.labels, y.links=new.links, y.images=new.images,  fname.root=fname.root, dir=dir,cleanDir=cleanDir, ...)
 
  }


  Splot = makeSplot(Splot, fname.root=fname.root, dir=dir,overwriteSourcePlot =overwriteSourcePlot,makeInteractive=makeInteractive,overrideInteractive=overrideInteractive,header=header,window.size=window.size,returnObj=TRUE)

  
}
