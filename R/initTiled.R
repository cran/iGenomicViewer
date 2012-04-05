#
# newTiled.R
#   ititTile
#   makeTiled
#


# the mapObj must be made with start and stop locations
#  not a central location
initTile <- function(Z,
                     bacDX,
                     goodDX=NA,
                     mapObj=NA,
                     H=2,                     
                     zlims=c(-0.5,0.5),
                     smplDX =NA,
                     ylabels=NA,
                     xlabels=NA,
                     x.axis.cex =0.5,
                     y.axis.cex =0.5,
                     xlab="samples",
                     ylab="BAC location",
                     ttl=NA,
                     returnVl=TRUE,
                     saveFlag=FALSE,
                     saveName="TIplot.RData"){

  if(is.na(smplDX[1])) smplDX = 1:dim(Z)[2]


  # set goodDX default
  if(is.na(goodDX[1])) goodDX = 1:dim(mapObj$mapping.info)[1]
  origDX=bacDX
  bacDX = intersect(bacDX, goodDX)
  
  # check Z dimensions
  if(dim(Z)[1] > length(bacDX)) Z = as.matrix(Z[bacDX,], nrows=length(bacDX))
  if(dim(Z)[2] >= length(smplDX)) Z = as.matrix(Z[,smplDX],nrows=length(bacDX))
  if((dim(Z)[1] != length(bacDX)) & (dim(Z)[2] != length(smplDX))) {
    stop("Error: dimensions of Z not correct \n       dimension of Z must be equal to length of bacDX by length of smplDX \n       or a larger matrix to be subset based on bacDX and smplDX\n") 
  }

 
  # start and stop matrix
  Y = matrix(NA, length(bacDX),2)
  Y[,1] = mapObj$mapping.info$g.loc.start[bacDX]
  Y[,2] = mapObj$mapping.info$g.loc.stop[bacDX]

  # check mapping object - or use default program mapObj
  if(is.na(mapObj[[1]][[1]][1])){
    data(mapping.info)
    mapping.info=mapping.info
    mapObj = mapping.info
  }

  # check additiional xlabels and ylabels 
  if(!is.na(ylabels[1])){
    if((length(ylabels) != length(bacDX))){
      if(length(ylabels) >= length(origDX)){
        #ylabels = ylabels[bacDX]
        if(length(which(is.na(match(origDX, bacDX)))) != 0)   ylabels = ylabels[-which(is.na(match(origDX, bacDX)))]
      }else{
        ylabels = as.character(1:dim(Y)[1])
      }
    }
  }
  if(!is.na(xlabels[1])){
    if((length(xlabels) != length(smplDX))){
      if(length(xlabels) >= length(smplDX)){
        xlabels = xlabels[smplDX]
      }else{
        xlabels = as.character(1:dim(Z)[2])
      }
    }
  }
  
  if(sum(is.na(zlims))==0){
    Z[Z>zlims[2]]=zlims[2]
    Z[Z<zlims[1]]=zlims[1]
  }
    
  Zcol=c(hsv(h=2/6,v=seq(1,0,length=100)),
              hsv(h=0/6,v=seq(0,1,length=100)))


# X is I x J
# where I = # BACs
#       J = # samples
# Y is I x 2

  I=dim(Z)[1]
  J=dim(Z)[2]

  Xtract=rep(1:H,ceiling(I/H))[1:I]
# Xtract is I long and indicates which tract each BAC falls in

  Ylim=range(as.vector(Y))
  
  # now get boundaries for each tract..

  tractBound=pairlist()
  for(h in 1:H){
    hdx=which(Xtract==h)

    # check for boundary overlap..
    Ysegs=Y[hdx,]
    ovdx=which(Ysegs[1:length(hdx)-1,2]>Ysegs[2:length(hdx),1])
    if(length(ovdx)>0){
      cnd=TRUE
      while(cnd){
        ovdx=which(Ysegs[1:length(hdx)-1,2]>Ysegs[2:length(hdx),1])
        if(length(ovdx)==0) cnd=FALSE
        if(length(ovdx)>0){
          for(odx in ovdx){
            midpt=(Ysegs[odx+1,1]+Ysegs[odx,2])/2
            Ysegs[odx+1,1]=midpt
            Ysegs[odx,2]=midpt
          }
        }
      }# end while(cnd)
    }# end if(length(ovdx)>0)

    Ubreaks=sort(unique(c(Ylim,as.vector(Ysegs))))
    Ymap=match(Ysegs[,1],Ubreaks)
    tractBound=c(tractBound,list(list(tbound=Ubreaks,Ymap=Ymap,
                   xdx=matrix(1:(H*J),ncol=J,byrow=FALSE)[h,],ydx=hdx)))
  }# end loop over H


  # displays values for tracks
  #   tractBound
  
  Xcoords=rep((1:H)/(H+1),J)+as.vector(mapply(rep,1:J,MoreArgs=list(times=H)))-1
  Xlim=c(0,J)

  TIplot = list()
  TIplot$tractBound = tractBound
  TIplot$vls = list(Y=Y, Z=Z, H=H, Xcoords=Xcoords, nsmp=J, nBAC=I, smplDX=smplDX) 
  TIplot$lims = list(xlim=Xlim, ylim=Ylim, zlim=zlims)
  TIplot$labels = list(xlab=xlab, ylab=ylab, ttl=ttl, xlabels=xlabels, ylabels=ylabels)
  TIplot$Zcol = Zcol
  TIplot$cex = list(xcex = x.axis.cex, ycex = y.axis.cex)
  TIplot$map = list(bacDX=bacDX, mapObj = mapObj)

  # specify class 
  class(TIplot) <- "TIplot"

  if(saveFlag) save(TIplot, file=saveName, compress=TRUE)
    if(returnVl) return(TIplot)
  
}




