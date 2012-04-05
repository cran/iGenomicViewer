
  
makeGGV <- function(GGV,
                    goodDX=NA,
                    smplDX=NA,
                    smp.color=NA,

                    break.num = 125,
                    tileNum = 2, 
                    buffer = 5,
                    
                    makeWinArms=TRUE, 

                    tiledMat=NA,
                    tiledMai.mat = NA,
                    tiledMai.prc=FALSE,
                    
                    fname.root="iGGV",
                    dir="GGV/",
                    overwriteSourcePlot = NA,
                    header="v3",
                    window.size = "800x1100",
                    image.size = "800x1100",
                    tiled.window.size = "800x1100",
                    tiled.image.size = "1200x1100",
                    
                    cleanDir = TRUE
  
                    ){

  # set and create directory for files
  direct = dir
  if(eval.js(paste("length(which(dir() =='", substr(direct,start=1, stop=nchar(direct)-1), "'))", sep="")) == 0) system(paste("mkdir", direct)) 
  subdirect = "GGVArms/"
  if(eval.js(paste("length(which(dir(direct) =='", substr(subdirect,start=1, stop=nchar(subdirect)-1), "'))", sep="")) == 0) system(paste("mkdir ", direct,"/", subdirect,sep="")) 
  subdirect = "GGVwinArms/"
  if(eval.js(paste("length(which(dir(direct) =='", substr(subdirect,start=1, stop=nchar(subdirect)-1), "'))", sep="")) == 0) system(paste("mkdir ", direct,"/", subdirect,sep="")) 
  subdirect = "GGVtiled/"
  if(eval.js(paste("length(which(dir(direct) =='", substr(subdirect,start=1, stop=nchar(subdirect)-1), "'))", sep="")) == 0) system(paste("mkdir ", direct,"/", subdirect,sep="")) 


  mapObj = GGV$values$mapObj
  annObj = GGV$values$annObj
  
  # set goodDX default
  if(is.na(goodDX[1])) goodDX = 1:dim(mapObj$mapping.info)[1]
  #origDX=bacDX
  #bacDX = intersect(bacDX, goodDX)

  if(is.na(smplDX[1])) smplDX = 1:dim(GGV$values$vls)[2]
  if(is.na(smp.color[1])) smp.color = rep("black", length(smplDX))
  

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

    environment(iGGV) <- environment()
    environment(iGGVtiled) <- environment()


  
  ##########################
  #
  # make cover pages 
  #
  ##########################

#  cat("Making Index file...\n")
  
  # index file
  chrArms=GGV$values$chrArms




###################################################################################
###################################################################################
###################################################################################
  








  # genomic plot if applicable

  if(!is.na(GGV$values$plot.vec[1])){
    
    cat("Genomic plot index file...\n")
 
    # if plot.dx is NA assume across all spot.DX
    if(is.na(GGV$values$plot.dx[1])){
      wasNA = TRUE
      GGV$values$plot.dx= intersect(1:dim(mapObj$mapping.info)[1], goodDX)
    }else{
      GGV$values$plot.dx= intersect(GGV$values$plot.dx, goodDX)
      wasNA = FALSE
    }

    # get genomic start and stop locations for additional plot
    gloc.start = min(mapObj$mapping.info$g.loc.start[GGV$values$plot.dx],na.rm=TRUE)
    gloc.stop = max(mapObj$mapping.info$g.loc.stop[GGV$values$plot.dx],na.rm=TRUE)

    
    # if chrArm does not include arms in the genomic plot, add them
    if(!wasNA) chrArms = unique(c(chrArms, as.character(mapObj$band.info$Arm$Arm[max(which(mapObj$band.info$Arm$lower <= gloc.start)):min(which(mapObj$band.info$Arm$upper >= gloc.stop))])))
 
    
    # go off genomic centers
    g.loc.spots = mapObj$mapping.info$g.loc.center[GGV$values$plot.dx]
    arm = rep(NA, length(g.loc.spots))

    for(i in 1:dim(mapObj$band.info$Arm)[1]){
      arm[intersect(which(g.loc.spots >= mapObj$band.info$Arm$lower[i]), which(g.loc.spots <= mapObj$band.info$Arm$upper[i]))] = as.character(mapObj$band.info$Arm$Arm[i])      
    }
    
    # make interactive to GGVArms
    asLinks=paste("GGVArms/", arm, ".html",sep="")
    asLinks[which(is.na(arm))] = NA



#--------------------------------------------------------------------------------


    #
    # set up blank plot
    #

    
    # common reference from other programs 
    spotDX = GGV$values$plot.dx
    
    # figure out labels
    idx.Chrom = intersect(which(GGV$values$mapObj$band.info$Chrom$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Chrom$upper >= gloc.start))
    if(gloc.start < GGV$values$mapObj$band.info$Chrom$lower[idx.Chrom[1]]) idx.Chrom = c(idx.Chrom[1]-1, idx.Chrom)
    if(gloc.stop > GGV$values$mapObj$band.info$Chrom$upper[idx.Chrom[length(idx.Chrom)]]) idx.Chrom = c(idx.Chrom, idx.Chrom[length(idx.Chrom)]+1)
    idx.Arm = intersect(which(GGV$values$mapObj$band.info$Arm$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Arm$upper >= gloc.start))
    if(gloc.start < GGV$values$mapObj$band.info$Arm$lower[idx.Arm[1]]) idx.Arm = c(idx.Arm[1]-1, idx.Arm)
    if(gloc.stop > GGV$values$mapObj$band.info$Arm$upper[idx.Arm[length(idx.Arm)]]) idx.Arm = c(idx.Arm, idx.Arm[length(idx.Arm)]+1)
    idx.Broad.Band = intersect(which(GGV$values$mapObj$band.info$Broad.Band$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Broad.Band$upper >= gloc.start))
    if(gloc.start < GGV$values$mapObj$band.info$Broad.Band$lower[idx.Broad.Band[1]]) idx.Broad.Band = c(idx.Broad.Band[1]-1, idx.Broad.Band)
    if(gloc.stop > GGV$values$mapObj$band.info$Broad.Band$upper[idx.Broad.Band[length(idx.Broad.Band)]]) idx.Broad.Band = c(idx.Broad.Band, idx.Broad.Band[length(idx.Chrom)]+1)
    idx.Fine.Band = intersect(which(GGV$values$mapObj$band.info$Fine.Band$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Fine.Band$upper >= gloc.start))
    if(gloc.start < GGV$values$mapObj$band.info$Fine.Band$lower[idx.Fine.Band[1]]) idx.Fine.Band = c(idx.Fine.Band[1]-1, idx.Fine.Band)
    if(gloc.stop > GGV$values$mapObj$band.info$Fine.Band$upper[idx.Fine.Band[length(idx.Fine.Band)]]) idx.Fine.Band = c(idx.Fine.Band, idx.Fine.Band[length(idx.Chrom)]+1)

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

    if(nSpots <= GGV$info$maxLabels){
      Spots.Used = TRUE
      Chrom.Used = FALSE


      # get section divides
      lwrDX = idx.spots[1]
      uprDX = idx.spots[length(idx.spots)]
      sectiondiv = GGV$values$mapObj$mapping.info$g.loc.start[lwrDX]
      bndDX = lwrDX:uprDX
      bndDX = intersect(bndDX,goodDX)
           
      #for(i in lwrDX:(uprDX-1)){
      for(i in bndDX[-length(bndDX)]){
        sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$mapping.info$g.loc.stop[i], GGV$values$mapObj$mapping.info$g.loc.start[i+1])))
      }
      sectiondiv = c(sectiondiv, GGV$values$mapObj$mapping.info$g.loc.stop[uprDX])
      # get section centers
      sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
      for(i in 2:(length(sectiondiv)-1)){
        sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
      }

      # back check if only on location
      if(lwrDX==uprDX){
        sectiondiv = c(GGV$values$mapObj$mapping.info$g.loc.start[lwrDX],GGV$values$mapObj$mapping.info$g.loc.stop[uprDX]) 
        sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
      }      

      
      # get sectio labels
      sectionlbls = as.character(GGV$values$mapObj$mapping.info$Spot.ID)[idx.spots]
                
        
      
    }else{
      if(nFine.Band <=GGV$info$maxLabels){
        Fine.Band.Used=TRUE
        Chrom.Used=FALSE
        
        # get section divides
        lwrDX = idx.Fine.Band[1]
        uprDX = idx.Fine.Band[length(idx.Fine.Band)]
        sectiondiv = GGV$values$mapObj$band.info$Fine.Band$lower[lwrDX]
        for(i in lwrDX:(uprDX-1)){
          sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Fine.Band$upper[i], GGV$values$mapObj$band.info$Fine.Band$lower[i+1])))
        }
        sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Fine.Band$upper[uprDX])
        # get section centers
        sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
        for(i in 2:(length(sectiondiv)-1)){
          sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
        }

        # back check if only one location
        if(lwrDX == uprDX){
          sectiondiv = c(GGV$values$mapObj$band.info$Fine.Band$lower[lwrDX], GGV$values$mapObj$band.info$Fine.Band$upper[uprDX])
          sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
        }
       
        # get section labels
        sectionlbls = as.character(GGV$values$mapObj$band.info$Fine.Band$label)[idx.Fine.Band]
 
      }else{
        if(nBroad.Band<=GGV$info$maxLabels){
          Broad.Band.Used=TRUE
          Chrom.Used=FALSE
        
          # get section divides
          lwrDX = idx.Broad.Band[1]
          uprDX = idx.Broad.Band[length(idx.Broad.Band)]
          sectiondiv = GGV$values$mapObj$band.info$Broad.Band$lower[lwrDX]
          for(i in lwrDX:(uprDX-1)){
            sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Broad.Band$upper[i], GGV$values$mapObj$band.info$Broad.Band$lower[i+1])))
          }
          sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Broad.Band$upper[uprDX])
          # get section centers
          sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
          for(i in 2:(length(sectiondiv)-1)){
            sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
          }

          # back check if only on location
          if(lwrDX==uprDX){
            sectiondiv = c(GGV$values$mapObj$band.info$Broad.Band$lower[lwrDX], GGV$values$mapObj$band.info$Broad.Band$upper[uprDX])
            sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))            
          }
          
          # get section labels
          sectionlbls = as.character(GGV$values$mapObj$band.info$Broad.Band$label)[idx.Broad.Band]
 
        }else{
          if(nArm<=GGV$info$maxLabels){
            Arm.Used=TRUE
            Chrom.Used=FALSE
        
          # get section divides
            lwrDX = idx.Arm[1]
            uprDX = idx.Arm[length(idx.Arm)]
            sectiondiv = GGV$values$mapObj$band.info$Arm$lower[lwrDX]
            for(i in lwrDX:(uprDX-1)){
              sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Arm$upper[i], GGV$values$mapObj$band.info$Arm$lower[i+1])))
            }
            sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Arm$upper[uprDX])
            # get section centers
            sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
            for(i in 2:(length(sectiondiv)-1)){
              sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
            }

            # back check if only one location
            if(uprDX==lwrDX){
              sectiondiv = c(GGV$values$mapObj$band.info$Arm$lower[lwrDX],GGV$values$mapObj$band.info$Arm$upper[uprDX])
              sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
            }            
            # get section labels
            sectionlbls = as.character(GGV$values$mapObj$band.info$Arm$label)[idx.Arm]
            
          }
        }
      }
    }
 
   if(Chrom.Used){
        
      # get section divides
      lwrDX = idx.Chrom[1]
      uprDX = idx.Chrom[length(idx.Chrom)]
      sectiondiv = GGV$values$mapObj$band.info$Chrom$lower[lwrDX]
      for(i in lwrDX:(uprDX-1)){
        sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Chrom$upper[i], GGV$values$mapObj$band.info$Chrom$lower[i+1])))
      }
      sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Chrom$upper[uprDX])
      # get section centers
      sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
      for(i in 2:(length(sectiondiv)-1)){
        sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
      }

      # back check if only one location
      if(uprDX==lwrDX){
        sectiondiv = c(GGV$values$mapObj$band.info$Chrom$lower[lwrDX], GGV$values$mapObj$band.info$Chrom$upper[uprDX])
        sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
   
      }
      
      # get section labels
      sectionlbls = as.character(GGV$values$mapObj$band.info$Chrom$label)[idx.Chrom]
     
    }

  
 

    #-------------------------------------------------

    #
    # add points to plot
    #

    plot.vec = eval.js(GGV$values$plot.vec)  # should be across entire genome
    bacDX = GGV$values$plot.dx

    # if plot.vec is longer than number spots assume multiple values
    # for spots 
    if(length(eval.js(plot.vec)) > dim(mapObj$mapping.info)[1]){
      nt = length(eval.js(plot.vec))/(dim(mapObj$mapping.info)[1])
      new.bacDX = bacDX
      for(i in 2:nt){
        new.bacDX = c(new.bacDX, (bacDX + ((dim(mapObj$mapping.info)[1])*(nt -1))))
      }
      new.plot.vec = plot.vec[new.bacDX]
    }else{
      nt = 1
      new.bacDX = bacDX
      new.plot.vec = plot.vec[new.bacDX]
    }

    

    #abline(h=sectiondiv,col='gray57',lty=2);
    side.plot.call = "plot(0,0,pch=' ',xlab=' ',ylab=' ',main=' ',ylim=range(sectiondiv,na.rm=TRUE),xlim=range(eval.js(new.plot.vec),na.rm=TRUE),axes=FALSE);    axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=.7,las=2);axis(1); points(x=new.plot.vec, y=rep(mapObj$mapping.info$g.loc.center[bacDX],nt),pch=21, bg='black') "
    

       
    # plot.call is like plot.extras now -- add color points, lines, etc. through this call
    # default only plots interactive points 

    
    # initialize and make interactive
    #Splot = initSplot(mat=1, plot.calls= side.plot.call, Iflag=TRUE, plot.extras=GGV$values$side.plot.extras)   
    Splot = initSplot(mat=1, plot.calls= side.plot.call, Iflag=TRUE, plot.extras=GGV$values$side.plot.extras)   

    Splot = makeImap(Splot, xy.type='points', x.pos=new.plot.vec,
                     y.pos=rep(mapObj$mapping.info$g.loc.center[bacDX],nt),
                     spot.radius=10, dir=dir, fname.root=fname.root,
                     x.labels=rep(arm, nt),
                     asLinks=rep(asLinks, nt), cleanDir=cleanDir)
 
    
    # make plot
    makeSplot(Splot, fname.root=fname.root, dir=dir,
              overwriteSourcePlot=overwriteSourcePlot,header=header,
              window.size=window.size, returnObj=FALSE)              



  }

  
###################################################################################
###################################################################################
###################################################################################
  


  
 ##################################
 #
 # Add Known Tracks to annObj
 #
 ##################################

  TR = GGV$values$trackRegions
  TRdo = TRUE
  annObj = GGV$values$annObj
  len = length(annObj)
  clrs = GGV$info$clrs

  TRdo = TRUE
  if(length(TR) == 1){
    if(is.na(TR)) TRdo=FALSE
  }

  if(TRdo){
  
  #
  # four ways to get location: broad, fine, g.loc, geneName
  #  cycle over to get locations
  #

  #
  # broad band
  #
  if(!is.null(TR$Broad.Band)){
    broad.band = TR$Broad.Band
    if(!is.na(broad.band[1])){
      # match broad band to band.info object
      mm = match(broad.band,GGV$values$mapObj$band.info$Broad.Band$Broad.Band)
      # subset based on location
      bb.DF = GGV$values$mapObj$band.info$Broad.Band[mm,c(1,3:5, 3:5)]
      Chrm = rep(NA, dim(bb.DF)[1])
      # determine chromosome of region 
      for(bb.ch  in 1:length(Chrm)){
        Chrm[bb.ch] = as.character(GGV$values$mapObj$band.info$Chrom$Chrom[intersect(which(GGV$values$mapObj$band.info$Chrom$lower <= bb.DF[bb.ch,3]),  which(GGV$values$mapObj$band.info$Chrom$upper >= bb.DF[bb.ch,3]))])
      }
      bb.DF$Chrom = Chrm
      # adjust labels 
      names(bb.DF) = c("Label","loc.start", "loc.center", "loc.stop","g.loc.start", "g.loc.center", "g.loc.stop", "Chrom")
      bb.DF = bb.DF[,c(1,8,2:7)]
      # make region interactive 
      bb.DF.links = paste("../GGVtiled/", as.character(bb.DF$Label),".html", sep="")
    }else{
      bb.DF = rep(NA,8)
      bb.DF.links=NA
    }    
  }else{
    bb.DF = rep(NA,8)
    bb.DF.links=NA
   }
  #
  # fine band
  #
  if(!is.null(TR$Fine.Band)){
    fine.band = TR$Fine.Band
    if(!is.na(fine.band[1])){
      # match fine band to band.info object
      mm = match(fine.band,GGV$values$mapObj$band.info$Fine.Band$Fine.Band)
      # subset based on location
      fb.DF = GGV$values$mapObj$band.info$Fine.Band[mm,c(1,3:5, 3:5)]
      Chrm = rep(NA, dim(fb.DF)[1])
      # determine chromosome of region 
      for(fb.ch  in 1:length(Chrm)){
        Chrm[fb.ch] = as.character(GGV$values$mapObj$band.info$Chrom$Chrom[intersect(which(GGV$values$mapObj$band.info$Chrom$lower <= fb.DF[fb.ch,3]),  which(GGV$values$mapObj$band.info$Chrom$upper >= fb.DF[fb.ch,3]))])
      }
      fb.DF$Chrom = Chrm
     # adjust labels 
      names(fb.DF) = c("Label","loc.start", "loc.center", "loc.stop","g.loc.start", "g.loc.center", "g.loc.stop", "Chrom")
      fb.DF = fb.DF[,c(1,8,2:7)]
      # make region interactive 
      fb.DF.links = paste("../GGVtiled/", as.character(fb.DF$Label),".html", sep="")
    }else{
      fb.DF = rep(NA,8)
      fb.DF.links = NA
    }    
  }else{
    fb.DF =rep(NA,8)
    fb.DF.links = NA
  }
  #
  # genomic locations
  #
  if(!is.null(TR$genomicLoc)){
    if(!is.null(dim(TR$genomicLoc))){
      gl.DF = as.data.frame(TR$genomicLoc)
      names(gl.DF) = c("loc.start", "loc.stop")
      new.DF = as.data.frame(matrix(NA, ncol=2))
      names(new.DF) = c("loc.start", "loc.stop")
      k.region = numeric(length=0)
      for(tch in 1:dim(gl.DF)[1]){
        # if span more than one chromosome, break into seperate regions
        # should be small regions so will never cross more than one other chromosome
        if(length(intersect(which(GGV$values$mapObj$band.info$Chrom$lower <= gl.DF$loc.start[tch]),which(GGV$values$mapObj$band.info$Chrom$upper >= gl.DF$loc.stop[tch]))) == 0){
          new.DF = rbind(new.DF, c(gl.DF[tch,1], GGV$values$mapObj$band.info$Chrom$upper[which(GGV$values$mapObj$band.info$Chrom$lower <= gl.DF$loc.start[tch])[length(which(GGV$values$mapObj$band.info$Chrom$lower <= gl.DF$loc.start[tch]))]]))
          new.DF = rbind(new.DF, c(GGV$values$mapObj$band.info$Chrom$lower[which(GGV$values$mapObj$band.info$Chrom$upper >= gl.DF$loc.stop[tch])[1]], gl.DF[tch,2]))
          k.region = c(k.region, paste(tch, "A",sep=""), paste(tch,"B",sep=""))
        }else{
          new.DF = rbind(new.DF, gl.DF[tch,])
          k.region = c(k.region,tch)
        }
      }
      new.DF = new.DF[-1,]
      orig.DF = gl.DF
      gl.DF = new.DF
      Chrm = character(length=0)
      # determine chromosome for region
      for(tch in 1:dim(gl.DF)[1]){
        Chrm = c(Chrm, as.character(GGV$values$mapObj$band.info$Chrom$Chrom[intersect(which(GGV$values$mapObj$band.info$Chrom$lower <= gl.DF$loc.start[tch]),which(GGV$values$mapObj$band.info$Chrom$upper >= gl.DF$loc.stop[tch]))]))
      }
      gl.DF$loc.center = (gl.DF$loc.start + gl.DF$loc.stop)/2
      gl.DF = gl.DF[,c(1,3,2)]
      k.region = paste("KnownRegion", k.region, sep="")
      gl.DF = cbind(k.region,Chrm, gl.DF, gl.DF)
      # adjust labels 
      names(gl.DF) = c("Label", "Chrom", "loc.start", "loc.center", "loc.stop", "g.loc.start", "g.loc.center", "g.loc.stop")
      # make regions interactive 
      gl.DF.links = paste("../GGVtiled/", as.character(gl.DF$Label),".html", sep="")
    }else{
      gl.DF = rep(NA,8)
      gl.DF.links=NA
    }
  }else{
    gl.DF = rep(NA,8)
    gl.DF.links=NA
  }
  #
  # gene names 
  #
  if(!is.null(TR$geneName)){
    if(!is.na(TR$geneName[1])){
      gn.DF = as.data.frame(matrix(NA, ncol=8, nrow=length(TR$geneName))) 
      gn.DF[,1] = TR$geneName
      # use all possible annotation object to find match
      for(adx in 1:length(annObj)){
        dx = match(TR$geneName,annObj[[adx]]$annotation$Label)
        gdx = which(!is.na(dx))
        if(length(which(is.na(dx)))!= 0) dx = dx[-(which(is.na(dx)))]
        if(length(dx) != 0 ){
          gn.DF[gdx,3:8] =  annObj[[adx]]$annotation[dx,3:8]
          gn.DF[gdx,2] = as.character(annObj[[adx]]$annotation$Chrom[dx])
        }
      }
      # adjust labels
      names(gn.DF) =  c("Label", "Chrom", "loc.start", "loc.center", "loc.stop", "g.loc.start", "g.loc.center", "g.loc.stop")
      # make regions interactive
      gn.DF.links = paste("../GGVtiled/", as.character(gn.DF$Label),".html", sep="")

      # if match not found remove from list
      notFound = which(is.na(gn.DF$Chrom))
      if(length(notFound) !=0){
        gn.DF = gn.DF[-(notFound),]
        gn.DF.links = gn.DF.links[-(notFound)]
      }
      if(dim(gn.DF)[1] == 0){
        gn.DF = rep(NA,8)
        gn.DF.links=NA
      }
    }else{
      gn.DF = rep(NA,8)
      gn.DF.links=NA
    }
  }else{
    gn.DF = rep(NA,8)
    gn.DF.links=NA
  }
  
  
  #
  # now combined into one object of known regions
  #
  master.DF = rbind(bb.DF, fb.DF, gl.DF, gn.DF)
  known.region= as.data.frame(list(knownlinks = c(bb.DF.links, fb.DF.links, gl.DF.links, gn.DF.links)))

  # remove regions not found
  notFound = which(is.na(master.DF$Label))
  if(length(notFound) !=0){
    master.DF = master.DF[-(notFound),]
    known.region=known.region[-(notFound),]
  }

  known.region = as.data.frame(known.region)

  
  #
  # set new annotation object with known regions
  #
  if(dim(master.DF)[1] != 0){

    annObj$KnownRegions = list(annotation = master.DF, links = known.region, images=NA) 
    annotation = GGV$info$annotation
    if(is.na(annotation[1])){
      annotation = 1:length(annObj)
    }else{
      annotation = c(annotation, length(annObj))
    }
    GGV$info$annotation = annotation
    GGV$info$clrs = c(clrs[1:len], "gray57")
    GGV$values$annObj = annObj
       
  }

  if(length(grep(GGV$mapping.info$Chrom, pattern="chr")) == 0){
    master.DF$Chrom = gsub(master.DF$Chrom, pattern="chr", replacement="")
  }
  


  
  ########################################
  #
  # make Tracks
  #   make known track regions  
  #   as tiled images
  #
  ########################################
  
  cat("..Working on Known Regions \n") 


  # check regions already made
  u = as.character(unique(master.DF$Label))
  listMade = dir(paste(dir, "GGVtiled/", sep=""))
  mMade = match(paste(u, ".html", sep=""), listMade)
  u = u[which(is.na(mMade))]

  # if updated object known regions may be different
  # if the count is different regenerate
  #
  # NOTE:  if when updating results in same number of knowRegions  but regions are different genomic locations
  #        must manually remove tracks in directory GGVtiled
  
  if(length(grep(as.character(master.DF$Label), pattern="KnownRegion")) != length(grep(pattern="KnownRegion", listMade))/2){
    if(length(grep(as.character(master.DF$Label), pattern="KnownRegion")) != 0){
      if(length(u) != 0)  u = unique(c(u,as.character(master.DF$Label)[grep(as.character(master.DF$Label), pattern="KnownRegion")]))
      if(length(u) == 0)  u = as.character(master.DF$Label)[grep(as.character(master.DF$Label), pattern="KnownRegion")]
    }
  }


}else{

  u = numeric(length=0)
}

  
  # if there are tiled plots to generate
  if(length(u) != 0){

    # cycle through all 
  for(ut in 1:length(u)){
    
    cat("  ..Working on Known Region ", ut, " of ", length(u), "\n") 
  
    # original tile number, adjusted if fails
    H = tileNum
    # get label and genomic start and stop locations
    udx = which(master.DF$Label == u[ut])
    g.start = master.DF$g.loc.start[min(udx,na.rm=TRUE)]
    g.stop = master.DF$g.loc.stop[max(udx,na.rm=TRUE)]

    # find spot index
    bacDX.tr = intersect(which(GGV$values$mapObj$mapping.info$g.loc.start <= g.stop), which(GGV$values$mapObj$mapping.info$g.loc.stop >= g.start))
    # if no bacs are in region...find closest region
    if(length(bacDX.tr) == 0) {

      bacDX.tr = c(max(which(GGV$values$mapObj$mapping.info$g.loc.start <= g.stop)[length(which(GGV$values$mapObj$mapping.info$g.loc.start <= g.stop))],1),min(which(GGV$values$mapObj$mapping.info$g.loc.stop >= g.start)[1], dim(GGV$values$mapObj$mapping.info)[1]))
      
    }

   
    bacDX.tr = intersect(bacDX.tr,goodDX)


    armDX = intersect(which(GGV$values$mapObj$band.info$Arm$lower <= GGV$values$mapObj$mapping.info$g.loc.stop[bacDX.tr[1]]),which(GGV$values$mapObj$band.info$Arm$upper >= GGV$values$mapObj$mapping.info$g.loc.start[bacDX.tr[length(bacDX.tr)]]))
    arm.start = GGV$values$mapObj$band.info$Arm$lower[armDX]
    arm.stop = GGV$values$mapObj$band.info$Arm$upper[armDX]
    
    
    # adjust so not to span 2 different chromosomes
    bacDX.tr = max( (bacDX.tr[1]-buffer),
      
                #    which(GGV$values$mapObj$mapping.info$g.loc.start == min(GGV$values$mapObj$mapping.info$g.loc.start[which(GGV$values$mapObj$mapping.info$Chrom == master.DF$Chrom[udx])]))
                        which(GGV$values$mapObj$mapping.info$g.loc.stop >= arm.start)[1]
      

                  ):min((bacDX.tr[length(bacDX.tr)]+buffer),

                 #       which(GGV$values$mapObj$mapping.info$g.loc.stop == max(GGV$values$mapObj$mapping.info$g.loc.stop[which(GGV$values$mapObj$mapping.info$Chrom == master.DF$Chrom[udx])]))
                        which(GGV$values$mapObj$mapping.info$g.loc.start <= arm.stop)[length(which(GGV$values$mapObj$mapping.info$g.loc.start <= arm.stop))]
                             
                  )
    

    bacDX.tr = intersect(bacDX.tr,goodDX)


    
    
    TIplot = try(initTile(Z = GGV$values$vls, bacDX=bacDX.tr, mapObj=GGV$values$mapObj, H=H,
                      ylabels=as.character(GGV$values$mapObj$mapping.info$Spot.ID[bacDX.tr]),
                      xlab="samples", smplDX=smplDX, 
                      ylab="",ttl=u[ut]), silent=TRUE)

    # if cannot tile because of dimension, take dimension down
    while((class(TIplot)=="try-error") & (H > 1)){
      H = H - 1
      TIplot = try(initTile(Z = GGV$values$vls, bacDX=bacDX.tr, mapObj=GGV$values$mapObj, H=H,
                      ylabels=as.character(GGV$values$mapObj$mapping.info$Spot.ID[bacDX.tr]),
                      xlab="samples",smplDX=smplDX, 
                      ylab="",
                      ttl=u[ut]), silent=TRUE)
    }











    #
    # set up blank plot
    #

    
    # common reference from other programs 
    spotDX = bacDX.tr
    
    # figure out labels
    idx.Chrom = intersect(which(GGV$values$mapObj$band.info$Chrom$lower <= g.stop), which(GGV$values$mapObj$band.info$Chrom$upper >= g.start))
    if(g.start < GGV$values$mapObj$band.info$Chrom$lower[idx.Chrom[1]]) idx.Chrom = c(idx.Chrom[1]-1, idx.Chrom)
    if(g.stop > GGV$values$mapObj$band.info$Chrom$upper[idx.Chrom[length(idx.Chrom)]]) idx.Chrom = c(idx.Chrom, idx.Chrom[length(idx.Chrom)]+1)
    idx.Arm = intersect(which(GGV$values$mapObj$band.info$Arm$lower <= g.stop), which(GGV$values$mapObj$band.info$Arm$upper >= g.start))
    if(g.start < GGV$values$mapObj$band.info$Arm$lower[idx.Arm[1]]) idx.Arm = c(idx.Arm[1]-1, idx.Arm)
    if(g.stop > GGV$values$mapObj$band.info$Arm$upper[idx.Arm[length(idx.Arm)]]) idx.Arm = c(idx.Arm, idx.Arm[length(idx.Arm)]+1)
    idx.Broad.Band = intersect(which(GGV$values$mapObj$band.info$Broad.Band$lower <= g.stop), which(GGV$values$mapObj$band.info$Broad.Band$upper >= g.start))
    if(g.start < GGV$values$mapObj$band.info$Broad.Band$lower[idx.Broad.Band[1]]) idx.Broad.Band = c(idx.Broad.Band[1]-1, idx.Broad.Band)
    if(g.stop > GGV$values$mapObj$band.info$Broad.Band$upper[idx.Broad.Band[length(idx.Broad.Band)]]) idx.Broad.Band = c(idx.Broad.Band, idx.Broad.Band[length(idx.Chrom)]+1)
    idx.Fine.Band = intersect(which(GGV$values$mapObj$band.info$Fine.Band$lower <= g.stop), which(GGV$values$mapObj$band.info$Fine.Band$upper >= g.start))
    if(g.start < GGV$values$mapObj$band.info$Fine.Band$lower[idx.Fine.Band[1]]) idx.Fine.Band = c(idx.Fine.Band[1]-1, idx.Fine.Band)
    if(g.stop > GGV$values$mapObj$band.info$Fine.Band$upper[idx.Fine.Band[length(idx.Fine.Band)]]) idx.Fine.Band = c(idx.Fine.Band, idx.Fine.Band[length(idx.Chrom)]+1)

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

    if(nSpots <= GGV$info$maxLabels){
      Spots.Used = TRUE
      Chrom.Used = FALSE

      # get section divides
      lwrDX = idx.spots[1]
      uprDX = idx.spots[length(idx.spots)]
      sectiondiv = GGV$values$mapObj$mapping.info$g.loc.start[lwrDX]
      for(i in lwrDX:(uprDX-1)){
        sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$mapping.info$g.loc.stop[i], GGV$values$mapObj$mapping.info$g.loc.start[i+1])))
      }
      sectiondiv = c(sectiondiv, GGV$values$mapObj$mapping.info$g.loc.stop[uprDX])
      # get section centers
      sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
      for(i in 2:(length(sectiondiv)-1)){
        sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
      }

      # back check if only on location
      if(lwrDX==uprDX){
        sectiondiv = c(GGV$values$mapObj$mapping.info$g.loc.start[lwrDX],GGV$values$mapObj$mapping.info$g.loc.stop[uprDX]) 
        sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
      }      
      # get section labels
      sectionlbls = as.character(GGV$values$mapObj$mapping.info$Spot.ID)[idx.spots]
                

     
    }else{
      if(nFine.Band <=GGV$info$maxLabels){
        Fine.Band.Used=TRUE
        Chrom.Used=FALSE
        
        # get section divides
        lwrDX = idx.Fine.Band[1]
        uprDX = idx.Fine.Band[length(idx.Fine.Band)]
        sectiondiv = GGV$values$mapObj$band.info$Fine.Band$lower[lwrDX]
        for(i in lwrDX:(uprDX-1)){
          sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Fine.Band$upper[i], GGV$values$mapObj$band.info$Fine.Band$lower[i+1])))
        }
        sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Fine.Band$upper[uprDX])
        
        # get section centers
        sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
        for(i in 2:(length(sectiondiv)-1)){
          sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
        }

        # back check if only one location
        if(lwrDX == uprDX){
          sectiondiv = c(GGV$values$mapObj$band.info$Fine.Band$lower[lwrDX], GGV$values$mapObj$band.info$Fine.Band$upper[uprDX])
          sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
        }
        
        
        # get sectionlabels
        sectionlbls = as.character(GGV$values$mapObj$band.info$Fine.Band$label)[idx.Fine.Band]

        
      }else{
        if(nBroad.Band<=GGV$info$maxLabels){
          Broad.Band.Used=TRUE
          Chrom.Used=FALSE
        
          # get section divides
          lwrDX = idx.Broad.Band[1]
          uprDX = idx.Broad.Band[length(idx.Broad.Band)]
          sectiondiv = GGV$values$mapObj$band.info$Broad.Band$lower[lwrDX]
          for(i in lwrDX:(uprDX-1)){
            sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Broad.Band$upper[i], GGV$values$mapObj$band.info$Broad.Band$lower[i+1])))
          }
          sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Broad.Band$upper[uprDX])
          # get section centers
          sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
          for(i in 2:(length(sectiondiv)-1)){
            sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
          }
 
          # back check if only on location
          if(lwrDX==uprDX){
            sectiondiv = c(GGV$values$mapObj$band.info$Broad.Band$lower[lwrDX], GGV$values$mapObj$band.info$Broad.Band$upper[uprDX])
            sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))            
          }
          
       
          # get section labels
          sectionlbls = as.character(GGV$values$mapObj$band.info$Broad.Band$label)[idx.Broad.Band]
 
        }else{
          if(nArm<=GGV$info$maxLabels){
            Arm.Used=TRUE
            Chrom.Used=FALSE
        
          # get section divides
            lwrDX = idx.Arm[1]
            uprDX = idx.Arm[length(idx.Arm)]
            sectiondiv = GGV$values$mapObj$band.info$Arm$lower[lwrDX]
            for(i in lwrDX:(uprDX-1)){
              sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Arm$upper[i], GGV$values$mapObj$band.info$Arm$lower[i+1])))
            }
            sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Arm$upper[uprDX])
            # get section centers
            sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
            for(i in 2:(length(sectiondiv)-1)){
              sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
            }
            # back check if only one location
            if(uprDX==lwrDX){
              sectiondiv = c(GGV$values$mapObj$band.info$Arm$lower[lwrDX],GGV$values$mapObj$band.info$Arm$upper[uprDX])
              sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
            }            

            # get sectio labels
            sectionlbls = as.character(GGV$values$mapObj$band.info$Arm$label)[idx.Arm]
            
          }
        }
      }
    }
    if(Chrom.Used){
        
      # get section divides
      lwrDX = idx.Chrom[1]
      uprDX = idx.Chrom[length(idx.Chrom)]
      sectiondiv = GGV$values$mapObj$band.info$Chrom$lower[lwrDX]
      for(i in lwrDX:(uprDX-1)){
        sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Chrom$upper[i], GGV$values$mapObj$band.info$Chrom$lower[i+1])))
      }
      sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Chrom$upper[uprDX])
      # get section centers
      sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
      for(i in 2:(length(sectiondiv)-1)){
        sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
      }

      # back check if only one location
      if(uprDX==lwrDX){
        sectiondiv = c(GGV$values$mapObj$band.info$Chrom$lower[lwrDX], GGV$values$mapObj$band.info$Chrom$upper[uprDX])
        sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
   
      }
         
      # get section labels
      sectionlbls = as.character(GGV$values$mapObj$band.info$Chrom$label)[idx.Chrom]
     
    }


    # fix sectiondiv and sectioncntr
    if(sectiondiv[1] < TIplot$lims$ylim[1]){
      sectiondiv[1] = TIplot$lims$ylim[1]
      sectioncntr[1] = mean(c(sectiondiv[1], sectiondiv[2]))
    }
    if(sectiondiv[length(sectiondiv)] > TIplot$lims$ylim[2]){
      sectiondiv[length(sectiondiv)] = TIplot$lims$ylim[2]
      sectioncntr[length(sectioncntr)] = mean(c(sectiondiv[length(sectiondiv)], sectiondiv[(length(sectiondiv)-1)]))
    }
    
    

    #
    # add points to plot
    #

    plot.vec.tr = eval.js(GGV$values$plot.vec)  # should be across entire genome
    

    doPlot=TRUE
    if(length(plot.vec.tr) == 1){
      if(is.na(plot.vec.tr))  doPlot = FALSE
    }


    if(doPlot){
      
    # if plot.vec.tr is longer than number spots assume multiple values
    # for spots 
    if(length(eval.js(plot.vec.tr)) > dim(mapObj$mapping.info)[1]){
      nt = length(eval.js(plot.vec.tr))/(dim(mapObj$mapping.info)[1])
      new.bacDX.tr = bacDX.tr
      for(i in 2:nt){
        new.bacDX.tr = c(new.bacDX.tr, (bacDX.tr + ((dim(mapObj$mapping.info)[1])*(nt -1))))
      }
      new.plot.vec.tr = plot.vec.tr[new.bacDX.tr]
    }else{
      nt = 1
      new.bacDX.tr = bacDX.tr
      new.plot.vec.tr = plot.vec.tr[new.bacDX.tr]
    }


        
#plot.call="image(x=0:1,y=0:1,z=matrix(rep(NA,4),ncol=2),xlim=c(range(plot.vec,na.rm=T)),ylim=range(mapObj$mapping.info$g.loc.center,na.rm=T),zlim=c(0,1),axes=F,xlab='',ylab='');points(x=plot.vec,y=mapObj$mapping.info$g.loc.center, pch=3, cex=0.5, col='red');axis(2);axis(1)"

#  image(x=TIplot$vls$Xcoords,y=TIplot$tractBound[[h]]$tbound,z=tmat, xlim=TIplot$lims$xlim,ylim=TIplot$lims$ylim,zlim=TIplot$lims$zlim,col=TIplot$Zcol,axes=FALSE,xlab=TIplot$labels$xlab,ylab=TIplot$labels$ylab, main=TIplot$labels$ttl)
    #tmat=array(NA,dim=c(length(TIplot$vls$Xcoords),
    #                   length(TIplot$tractBound[[1]]$tbound)-1))
   
#image(x=TIplot$vls$Xcoords,y=TIplot$tractBound[[1]]$tbound,z=tmat, xlim=TIplot$lims$xlim,ylim=TIplot$lims$ylim,zlim=TIplot$lims$zlim,axes=FALSE, add=TRUE);

 #plot(0,0,pch=' ',xlab=' ',ylab=' ',main=' ',ylim=range(sectiondiv,na.rm=TRUE),xlim=range(eval.js(new.plot.vec.tr),na.rm=TRUE),axes=FALSE)




    side.plot.call.tr = "image(x=0:1,y=0:1,z=matrix(rep(NA,4),ncol=2),xlim=c(range(eval.js(new.plot.vec.tr),na.rm=T)),ylim=TIplot$lims$ylim,zlim=c(0,1),axes=F,xlab='',ylab=''); axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=.6,las=2);axis(1); points(x=new.plot.vec.tr, y=rep(mapObj$mapping.info$g.loc.center[bacDX.tr],nt),pch=21, bg='black')"

  }else{
    side.plot.call.tr=NA
    new.plot.vec.tr=NA

  }
    
    #if(!is.na(GGV$values$side.plot.extras[1])) side.plot.call.tr = paste(side.plot.call.tr, GGV$values$side.plot.extras, sep=";")
    if(!is.na(GGV$values$side.plot.extras[1])) side.plot.call.tr = paste(side.plot.call.tr, GGV$values$side.plot.extras, sep=";")


   if(is.null(dim(tiledMat))){
      if(is.na(side.plot.call.tr[1]))  tiledMat =  matrix(c(rep(c(rep(1,12),3), 11), rep(c(rep(2,12),0),2)), ncol=13, byrow=TRUE)  
      if(!is.na(side.plot.call.tr[1])) tiledMat =  matrix(c(rep(c(rep(1,12),3,rep(4,2)), 11), rep(c(rep(2,12), rep(0,3)),2)), ncol=15, byrow=TRUE)
   
    }

    if(is.null(dim(tiledMai.mat))){
      if(is.na(side.plot.call.tr[1]))  tiledMat.mat = matrix(c(0.2,  0.8, 0.5, 0.05,1.2,  0.8, 0.3, 0.05,0.2, 0.05, 0.5, 0.2), byrow=TRUE, ncol=4)
      if(!is.na(side.plot.call.tr[1])) tiledMai.mat = matrix(c(0.2,  0.8, 0.5, 0.05, 1.2,  0.8, 0.3, 0.05,0.2, 0.05, 0.5, 0.2,0.2,  0.2, 0.5, 0.2), byrow=TRUE, ncol=4)
  
    }


  
  # if tiled object successfully made, attempt to  make plot
    if(class(TIplot) != "try-error"){

      suppressWarnings(
                       
      iGGVtiled(TIplot=TIplot,
                annObj=GGV$values$annObj,
                x.labels=GGV$interactive$x.labels,
                y.labels=GGV$interactive$y.labels,
                xy.labels=GGV$interactive$xy.labels,
                x.links=GGV$interactive$x.links,
                y.links=GGV$interactive$y.links,
                xy.links=GGV$interactive$xy.links,
                asLinks=GGV$interactive$asLinks,
                x.images=GGV$interactive$x.images,
                y.images=GGV$interactive$y.images,
                xy.images=GGV$interactive$xy.images,
                #mai.mat=GGV$info$mai.mat,
                #mai.prc=GGV$info$mai.prc,
                mat=tiledMat,
                mai.mat=tiledMai.mat,
                mai.prc=tiledMai.prc,
                plot.extras=GGV$info$plot.extras,
                smpLines=GGV$info$smpLines,
                divCol=GGV$info$divCol,        
                plot.call =side.plot.call.tr,
                plot.vec = new.plot.vec.tr,
                lims = GGV$info$lims,
                annotation = GGV$info$annotation,
                clrs = GGV$info$clrs,
                mapObj.columns=GGV$info$mapObj.columns,
                fname.root=u[ut],
                dir=paste(dir,"GGVtiled/", sep=""),
                overwriteSourcePlot=overwriteSourcePlot, 
                header=header,
                window.size=tiled.window.size,
                image.size=tiled.image.size,
                vrb=FALSE,
                cleanDir=cleanDir
                )

                       )

     }# if TIplot is valid
   
  } # loop over tracks

} # if there were tracks to do 







  ########################################
  #
  # make Level 1  
  #   Chromosome Arms with limit
  #     interactive ability 
  #
  ########################################
  
  # determine number of cuts for each chromosome arm 
  nBACs = rep(NA, dim(mapObj$band.info$Arm)[1])
  g.loc = mapObj$mapping.info$g.loc.center
  BACstart = rep(NA, dim(mapObj$band.info$Arm)[1])
  BACend = rep(NA, dim(mapObj$band.info$Arm)[1])
  for(i in 1:dim(mapObj$band.info$Arm)[1]){
    dx = intersect(which(g.loc >= mapObj$band.info$Arm$lower[i]),which(g.loc <= mapObj$band.info$Arm$upper[i]))
    nBACs[i] = length(dx)
    if(length(dx) != 0){
      BACstart[i] = dx[1]
      BACend[i] = dx[length(dx)]
    }
  }
  #ncut = ceiling(nBACs/100)
  locInfo = as.data.frame(list(Arms=mapObj$band.info$Arm$Arm, starting=BACstart, ending=BACend, nSpot=nBACs))
  

  indx.arms = chrArms


  # check which arms are already made
  listMade = dir(paste(dir, "GGVArms/", sep=""))
  mMade = match(paste(chrArms, ".html", sep=""), listMade)
  chrArms = chrArms[which(is.na(mMade))]

  # check which arms should be updated if new known tracks
  new.arms = rep(NA, length(u))
  if(length(u) !=0){
    udx = match(u, master.DF$Label)
    for(i in 1:length(udx)){
      chr = master.DF$Chrom[udx[i]]
      center = master.DF$g.loc.center[udx[i]]
      new.arms[i] = as.character(GGV$values$mapObj$band.info$Arm$Arm[intersect(which(GGV$values$mapObj$band.info$Arm$lower <= center), which(GGV$values$mapObj$band.info$Arm$upper >= center))])
    }    
  }
  if(length(which(is.na(new.arms))) > 0) new.arms = new.arms[-(which(is.na(new.arms)))]

  # update chrArms 
  if( (length(chrArms) != 0) & (length(new.arms) != 0) ){
    chrArms = unique(c(chrArms, new.arms))
    indx.arms = unique(c(indx.arms, new.arms))
  }
  if( (length(chrArms) == 0) & (length(new.arms) != 0) ){
    chrArms =  new.arms
    indx.arms = unique(c(indx.arms, new.arms))
 }
  if( (length(chrArms) != 0) & (length(new.arms) == 0) ){
    chrArms =  chrArms
    indx.arms = indx.arms
  }
  if( (length(chrArms) == 0) & (length(new.arms) == 0) ){
    chrArms =  NA
    indx.arms = indx.arms
  }



  ##########################
  #
  # make cover pages 
  #
  ##########################

  cat("Making Index file...\n")
  
  # index file
  #chrArms=GGV$values$chrArms
  sink(paste(dir,fname.root,".Index.html", sep=""))
  cat("<html> \n <head> \n <title>index</title> \n <style type=\"text/css\"> \n /* CSS GOES HERE */ \n .plot {border:1px solid \n  max-width:800px; max-height:1100px;} \n </style> \n </head> \n <body> \n")
  cat("<H3> Chromosome Arms <H3> ")

  for(i in 1:length(indx.arms)){
    cat(paste(" <A href= \" GGVArms/",indx.arms[i], ".html \"> ", indx.arms[i], "  </A> <BR> \n", sep=""))    
  }
  cat("</body>\n</html>\n")
  sink()


  
  # if chrArms is not empty
  if( (length(chrArms) != 1) | !is.na(chrArms[1])  ){

  # now only do chromosome arms listed in chrArms
  for(ch in chrArms){

    cat("..Working on ", ch, "\n") 
            
    
    # find chromsome arm in mapObj/info object
    idx = which(locInfo$Arms==ch)

    # if info for that arm exists
    if(length(idx) != 0){

      # make sure enough spots to make image
      if(locInfo$nSpot[idx] > 2){

        #  may need to alter 
        # add subsets chromosomes to annotation object 
        # add tracks to annotation
        
        annObj = GGV$values$annObj
        annotation = GGV$info$annotation
        clrs = GGV$info$clrs
        plot.extras=GGV$info$plot.extras
        
        spot.dx = locInfo$starting[idx]:locInfo$ending[idx]

        spot.dx = intersect(spot.dx, goodDX)
        
        
        # if making within arm breakdown
        # determine breaks for interactivity
        if(makeWinArms){
     

         
          
          # Breaks are 125 long
          if(locInfo$nSpot[idx] > break.num){
            cuts = c(1,(1:ceiling(locInfo$nSpot[idx]/break.num)*break.num))
            cuts[length(cuts)] = locInfo$nSpot[idx]
          }
          if(locInfo$nSpot[idx] <= break.num) cuts = c(1, locInfo$nSpot[idx])





          if(locInfo$nSpot[idx] > break.num){

            
            # find actual spot dx in mapObj
            t.cuts = cuts + (locInfo$starting[idx] - 1) 
            g.cuts = as.data.frame(list(startloc=rep(NA, length(cuts) - 1), centerloc=rep(NA, length(cuts)- 1), endloc=rep(NA, length(cuts) - 1)))
            g.cuts[1,1]= mapObj$mapping.info$g.loc.start[t.cuts[1]]
            if(dim(g.cuts)[1] > 1){
              g.cuts[2:dim(g.cuts)[1], 1] = mapObj$mapping.info$g.loc.start[t.cuts[2:(length(t.cuts)-1)]+1]
              g.cuts$endloc[1:(dim(g.cuts)[1]-1)] = g.cuts$startloc[2:dim(g.cuts)[1]]-1
            }
            g.cuts$endloc[dim(g.cuts)[1]] =  mapObj$mapping.info$g.loc.stop[locInfo$ending[idx]]
            g.cuts$centerloc = (g.cuts$startloc + g.cuts$endloc)/2
            
            cuts2 = cuts[-(length(cuts))] + ceiling(break.num/2)
            cuts2[length(cuts2)] = locInfo$nSpot[idx]

            t.cuts2 = cuts2 + (locInfo$starting[idx] - 1)
            g.cuts2 = as.data.frame(list(startloc=rep(NA, length(cuts2) - 1), centerloc=rep(NA, length(cuts2)- 1), endloc=rep(NA, length(cuts2) - 1)))
            g.cuts2[1,1]= mapObj$mapping.info$g.loc.start[t.cuts2[1]]
            if(dim(g.cuts2)[1] > 1){
              g.cuts2[2:dim(g.cuts2)[1], 1] = mapObj$mapping.info$g.loc.start[t.cuts2[2:(length(t.cuts2)-1)]+1]
              g.cuts2$endloc[1:(dim(g.cuts2)[1]-1)] = g.cuts2$startloc[2:dim(g.cuts2)[1]]-1
            }
            g.cuts2$endloc[dim(g.cuts2)[1]] =  mapObj$mapping.info$g.loc.stop[t.cuts2[length(t.cuts2)]]
            g.cuts2$centerloc = (g.cuts2$startloc + g.cuts2$endloc)/2

            if(g.cuts$endloc[dim(g.cuts)[1]] == g.cuts2$endloc[dim(g.cuts2)[1]]) g.cuts = g.cuts[1:(dim(g.cuts)[1] -1),]
            

            full.g = rbind(g.cuts, g.cuts2)
            #full.g = full.g[order(full.g$centerloc),]

            g.cuts$Region = 1:dim(g.cuts)[1]
            g.cuts2$Region = (dim(g.cuts)[1]+1):((dim(g.cuts)[1]+1)+(dim(g.cuts2)[1]-1))

          }else{
            t.cuts = cuts + (locInfo$starting[idx] - 1) 
            g.cuts = as.data.frame(list(startloc=rep(NA, length(cuts) - 1), centerloc=rep(NA, length(cuts)- 1), endloc=rep(NA, length(cuts) - 1)))
            g.cuts[1,1]= mapObj$mapping.info$g.loc.start[t.cuts[1]]
            g.cuts$endloc[dim(g.cuts)[1]] =  mapObj$mapping.info$g.loc.stop[locInfo$ending[idx]]
            g.cuts$centerloc = (g.cuts$startloc + g.cuts$endloc)/2

            full.g = g.cuts
            g.cuts$Region = 1:dim(g.cuts)[1]
            
          }
          
          
          spots.g.loc.c = mapObj$mapping.info$g.loc.center[spot.dx]

          regionDX = rep(NA, length(spots.g.loc.c))
          for(i in 1:length(spots.g.loc.c)){
            tempdx = intersect(which(spots.g.loc.c[i] >= full.g$startloc), which(spots.g.loc.c[i] <= full.g$endloc))
            if(length(tempdx) == 1){
              regionDX[i] = tempdx
            }else{
               regionDX[i] = tempdx[which(abs(full.g$centerloc[tempdx] - spots.g.loc.c[i]) == min(abs(full.g$centerloc[tempdx] - spots.g.loc.c[i]), na.rm=TRUE))]
            }
          }

          
           
          if(is.na(annotation[1])) annotation = 1:length(annObj)
                          
          # add tracks to annotation obj
          ann1 = as.data.frame(list(Label =paste("SubRegion",g.cuts$Region,sep=""), g.loc.start=g.cuts$startloc, g.loc.center=g.cuts$centerloc, g.loc.stop=g.cuts$endloc))          
          annlink1 = as.data.frame(list(SubRegion =  paste("../GGVwinArms/", ch, ".sub", g.cuts$Region,".html", sep="")))
          annObj$SubRegion1 = list(annotation = ann1, links=annlink1, images=NA)
          
           if(locInfo$nSpot[idx] > break.num){
             ann2 = as.data.frame(list(Label =paste("SubRegion",g.cuts2$Region,sep=""),g.loc.start=g.cuts2$startloc, g.loc.center=g.cuts2$centerloc, g.loc.stop=g.cuts2$endloc))
             annlink2 = as.data.frame(list(SubRegion = paste("../GGVwinArms/", ch, ".sub", g.cuts2$Region,".html", sep="")))

             
             annObj$SubRegion2 = list(annotation = ann2, links=annlink2, images=NA)

             annotation = c((length(annObj)-1), length(annObj), annotation)
             clrs = c("gray57", "gray57", clrs)

           }else{
             annotation = c(length(annObj), annotation)
             clrs = c("gray57", clrs)
           }
          
        

          # add breaks between otherwise continuous lines 
          plt = "lines(x=c(.5,.5), y=c((g.cuts$endloc[1]-50), (g.cuts$endloc[1]+50)),lwd=5, col='white')"
          for(ic in 2:(dim(g.cuts)[1]-1)){
            plt = paste(plt,paste("lines(x=c(.5,.5), y=c((",g.cuts$endloc[ic],"-50), (",g.cuts$endloc[ic],"+50)),lwd=5, col='white')",sep=""), sep=";")
          }

          if(locInfo$nSpot[idx] > break.num){
         
            for(ic in 1:(dim(g.cuts2)[1]-1)){
              plt = paste(plt,paste("lines(x=c(1.5,1.5), y=c((",g.cuts2$endloc[ic],"-50), (",g.cuts2$endloc[ic],"+50)),lwd=5, col='white')",sep=""), sep=";")
            }
            
          }
          
          if(is.na(plot.extras[3])){
            plot.extras[3] = plt
          }else{
            plt.dx = length(plot.extras[3])
            plot.extras[3][plt.dx+1] = plt
          }

          if(!is.na(GGV$values$plot.vec[1])){
            if(length(plot.extras) != 4) plot.extras[4]=NA
          }


              
          
        }# end if makeWinArms

  

        if(is.na(GGV$values$plot.vec[1])) overrideInteractive = c(FALSE, FALSE, TRUE)
        if(!is.na(GGV$values$plot.vec[1])) overrideInteractive = c(FALSE, FALSE, TRUE, TRUE)
        

        
        # common reference from other programs 
        spotDX = spot.dx

        gloc.start = GGV$values$mapObj$mapping.info$g.loc.start[spotDX[1]]
        gloc.stop = GGV$values$mapObj$mapping.info$g.loc.stop[spotDX[length(spotDX)]]
           
        # figure out labels
        idx.Chrom = intersect(which(GGV$values$mapObj$band.info$Chrom$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Chrom$upper >= gloc.start))
        if(gloc.start < GGV$values$mapObj$band.info$Chrom$lower[idx.Chrom[1]]) idx.Chrom = c(idx.Chrom[1]-1, idx.Chrom)
        if(gloc.stop > GGV$values$mapObj$band.info$Chrom$upper[idx.Chrom[length(idx.Chrom)]]) idx.Chrom = c(idx.Chrom, idx.Chrom[length(idx.Chrom)]+1)
        idx.Arm = intersect(which(GGV$values$mapObj$band.info$Arm$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Arm$upper >= gloc.start))
        if(gloc.start < GGV$values$mapObj$band.info$Arm$lower[idx.Arm[1]]) idx.Arm = c(idx.Arm[1]-1, idx.Arm)
        if(gloc.stop > GGV$values$mapObj$band.info$Arm$upper[idx.Arm[length(idx.Arm)]]) idx.Arm = c(idx.Arm, idx.Arm[length(idx.Arm)]+1)
        idx.Broad.Band = intersect(which(GGV$values$mapObj$band.info$Broad.Band$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Broad.Band$upper >= gloc.start))
        if(gloc.start < GGV$values$mapObj$band.info$Broad.Band$lower[idx.Broad.Band[1]]) idx.Broad.Band = c(idx.Broad.Band[1]-1, idx.Broad.Band)
        if(gloc.stop > GGV$values$mapObj$band.info$Broad.Band$upper[idx.Broad.Band[length(idx.Broad.Band)]]) idx.Broad.Band = c(idx.Broad.Band, idx.Broad.Band[length(idx.Chrom)]+1)
        idx.Fine.Band = intersect(which(GGV$values$mapObj$band.info$Fine.Band$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Fine.Band$upper >= gloc.start))
        if(gloc.start < GGV$values$mapObj$band.info$Fine.Band$lower[idx.Fine.Band[1]]) idx.Fine.Band = c(idx.Fine.Band[1]-1, idx.Fine.Band)
        if(gloc.stop > GGV$values$mapObj$band.info$Fine.Band$upper[idx.Fine.Band[length(idx.Fine.Band)]]) idx.Fine.Band = c(idx.Fine.Band, idx.Fine.Band[length(idx.Chrom)]+1)
        
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
        
        if(nSpots <= GGV$info$maxLabels){
          Spots.Used = TRUE
          Chrom.Used = FALSE
          
          # get section divides
          lwrDX = idx.spots[1]
          uprDX = idx.spots[length(idx.spots)]
          sectiondiv = GGV$values$mapObj$mapping.info$g.loc.start[lwrDX]
          for(i in lwrDX:(uprDX-1)){
            sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$mapping.info$g.loc.stop[i], GGV$values$mapObj$mapping.info$g.loc.start[i+1])))
          }
          sectiondiv = c(sectiondiv, GGV$values$mapObj$mapping.info$g.loc.stop[uprDX])
          # get section centers
          sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
          for(i in 2:(length(sectiondiv)-1)){
            sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
          }

          # back check if only on location
          if(lwrDX==uprDX){
            sectiondiv = c(GGV$values$mapObj$mapping.info$g.loc.start[lwrDX],GGV$values$mapObj$mapping.info$g.loc.stop[uprDX]) 
            sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
          }      
          
          # get sectio labels
          sectionlbls = as.character(GGV$values$mapObj$mapping.info$Spot.ID)[idx.spots]
                  
        }else{
          if(nFine.Band <=GGV$info$maxLabels){
            Fine.Band.Used=TRUE
            Chrom.Used=FALSE
        
            # get section divides
            lwrDX = idx.Fine.Band[1]
            uprDX = idx.Fine.Band[length(idx.Fine.Band)]
            sectiondiv = GGV$values$mapObj$band.info$Fine.Band$lower[lwrDX]
            for(i in lwrDX:(uprDX-1)){
              sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Fine.Band$upper[i], GGV$values$mapObj$band.info$Fine.Band$lower[i+1])))
            }
            sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Fine.Band$upper[uprDX])
            # get section centers
            sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
            for(i in 2:(length(sectiondiv)-1)){
              sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
            }

            # back check if only one location
            if(lwrDX == uprDX){
              sectiondiv = c(GGV$values$mapObj$band.info$Fine.Band$lower[lwrDX], GGV$values$mapObj$band.info$Fine.Band$upper[uprDX])
              sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
            }            
            # get section labels
            sectionlbls = as.character(GGV$values$mapObj$band.info$Fine.Band$label)[idx.Fine.Band]
            
          }else{
            if(nBroad.Band<=GGV$info$maxLabels){
              Broad.Band.Used=TRUE
              Chrom.Used=FALSE
        
              # get section divides
              lwrDX = idx.Broad.Band[1]
              uprDX = idx.Broad.Band[length(idx.Broad.Band)]
              sectiondiv = GGV$values$mapObj$band.info$Broad.Band$lower[lwrDX]
              for(i in lwrDX:(uprDX-1)){
                sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Broad.Band$upper[i], GGV$values$mapObj$band.info$Broad.Band$lower[i+1])))
              }
              sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Broad.Band$upper[uprDX])
              # get section centers
              sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              for(i in 2:(length(sectiondiv)-1)){
                sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
              }
 
              # back check if only on location
              if(lwrDX==uprDX){
                sectiondiv = c(GGV$values$mapObj$band.info$Broad.Band$lower[lwrDX], GGV$values$mapObj$band.info$Broad.Band$upper[uprDX])
                sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))            
              }          
                   
              # get section labels
              sectionlbls = as.character(GGV$values$mapObj$band.info$Broad.Band$label)[idx.Broad.Band]
              
            }else{
              if(nArm<=GGV$info$maxLabels){
                Arm.Used=TRUE
                Chrom.Used=FALSE
                
                # get section divides
                lwrDX = idx.Arm[1]
                uprDX = idx.Arm[length(idx.Arm)]
                sectiondiv = GGV$values$mapObj$band.info$Arm$lower[lwrDX]
                for(i in lwrDX:(uprDX-1)){
                  sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Arm$upper[i], GGV$values$mapObj$band.info$Arm$lower[i+1])))
                }
                sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Arm$upper[uprDX])
                # get section centers
                sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                for(i in 2:(length(sectiondiv)-1)){
                  sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
                }

                # back check if only one location
                if(uprDX==lwrDX){
                  sectiondiv = c(GGV$values$mapObj$band.info$Arm$lower[lwrDX],GGV$values$mapObj$band.info$Arm$upper[uprDX])
                  sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                }            
    
                # get section labels
                sectionlbls = as.character(GGV$values$mapObj$band.info$Arm$label)[idx.Arm]
                
              }
            }
          }
        }
        if(Chrom.Used){
          
          # get section divides
          lwrDX = idx.Chrom[1]
          uprDX = idx.Chrom[length(idx.Chrom)]
          sectiondiv = GGV$values$mapObj$band.info$Chrom$lower[lwrDX]
          for(i in lwrDX:(uprDX-1)){
            sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Chrom$upper[i], GGV$values$mapObj$band.info$Chrom$lower[i+1])))
          }
          sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Chrom$upper[uprDX])
          # get section centers
          sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
          for(i in 2:(length(sectiondiv)-1)){
            sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
          }
           
          # back check if only one location
          if(uprDX==lwrDX){
            sectiondiv = c(GGV$values$mapObj$band.info$Chrom$lower[lwrDX], GGV$values$mapObj$band.info$Chrom$upper[uprDX])
            sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
            
          }
      
          # get sectio labels
          sectionlbls = as.character(GGV$values$mapObj$band.info$Chrom$label)[idx.Chrom]
        }

        ylim.temp = range(mapObj$mapping.info$g.loc.center[spot.dx],na.rm=TRUE)
        
        # fix sectiondiv and sectioncntr
        if(sectiondiv[1] < ylim.temp[1]){
          sectiondiv[1] = ylim.temp[1]
          sectioncntr[1] = mean(c(sectiondiv[1], sectiondiv[2]))
        }
        if(sectiondiv[length(sectiondiv)] > ylim.temp[2]){
          sectiondiv[length(sectiondiv)] = ylim.temp[2]
          sectioncntr[length(sectioncntr)] = mean(c(sectiondiv[length(sectiondiv)], sectiondiv[(length(sectiondiv)-1)]))
        }
            
        #
        # add points to plot
        #

        plot.vec = eval.js(GGV$values$plot.vec)  # should be across entire genome

        doPlot=TRUE
        if(length(plot.vec) == 1){
          if(is.na(plot.vec))  doPlot = FALSE
        }


        if(doPlot){
        
       # if plot.vec is longer than number spots assume multiple values
       # for spots 
        if(length(eval.js(plot.vec)) > dim(mapObj$mapping.info)[1]){
          nt = length(eval.js(plot.vec))/(dim(mapObj$mapping.info)[1])
         
          new.bacDX = spot.dx
          for(i in 2:nt){
            new.bacDX = c(new.bacDX, (spot.dx + ((dim(mapObj$mapping.info)[1])*(nt -1))))
          }
          new.plot.vec = plot.vec[new.bacDX]
        }else{
          nt = 1
          new.bacDX = spot.dx
          new.plot.vec = plot.vec[new.bacDX]
        }





 #   side.plot.call.tr = "image(x=0:1,y=0:1,z=matrix(rep(NA,4),ncol=2),xlim=c(range(eval.js(new.plot.vec.tr),na.rm=T)),ylim=TIplot$lims$ylim,zlim=c(0,1),axes=F,xlab='',ylab=''); axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=.5,las=2);axis(1); points(x=new.plot.vec.tr, y=rep(mapObj$mapping.info$g.loc.center[bacDX.tr],nt),pch=21, bg='black') "

 #  side.plot.call = "plot(0,0,pch=' ',xlab=' ',ylab=' ',main=' ',ylim=range(sectiondiv,na.rm=TRUE),xlim=range(eval.js(new.plot.vec),na.rm=TRUE),axes=FALSE);    axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=.5,las=2);axis(1); points(x=new.plot.vec, y=rep(mapObj$mapping.info$g.loc.center[spot.dx],nt),pch=21, bg='black') "


        
        side.plot.call = "image(x=0:1,y=0:1,z=matrix(rep(NA,4),ncol=2),xlim=range(eval.js(new.plot.vec),na.rm=TRUE), ylim=range(mapObj$mapping.info$g.loc.center[spot.dx],na.rm=TRUE), zlim=c(0,1), axes=FALSE, xlab='', ylab=''); axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=.6,las=2);axis(1); points(x=new.plot.vec, y=rep(mapObj$mapping.info$g.loc.center[spot.dx],nt),pch=21, bg='black')"


       #  if(!is.na(GGV$values$side.plot.extras[1])) side.plot.call = paste(side.plot.call, GGV$values$side.plot.extras, sep=";")
        if(!is.na(GGV$values$side.plot.extras[1])) side.plot.call = paste(side.plot.call, GGV$values$side.plot.extras, sep=";")

      }else{
        side.plot.call=NA
        new.plot.vec=NA
      }
       
        
        bacDX = spot.dx 

        suppressWarnings(
                         
        iGGV(vls=GGV$values$vls,
             mapObj=GGV$values$mapObj,
             annObj=annObj,
             x.labels=GGV$interactive$x.labels,
             y.labels=GGV$interactive$y.labels,
             xy.labels=GGV$interactive$xy.labels,             
             x.links=GGV$interactive$x.links,
             y.links=GGV$interactive$y.links,
             xy.links=GGV$interactive$xy.links,
             asLinks=GGV$interactive$asLinks,
             x.images=GGV$interactive$x.images,
             y.images=GGV$interactive$y.images,
             xy.images=GGV$interactive$xy.images,
             plot.y.index=spot.dx,
             maxLabels=GGV$info$maxLabels,
             mat=GGV$values$mat,
             mai.mat=GGV$info$mai.mat,
             mai.prc=GGV$info$mai.prc,
             plot.x.index=smplDX,
             smp.color=smp.color,
             plot.extras=plot.extras,
             smpLines=GGV$info$smpLines,
             divCol=GGV$info$divCol,             
             plot.call= side.plot.call,
             plot.vec= new.plot.vec,
             lims=GGV$info$lims,
             annotation=annotation,
             clrs=clrs,
             mapObj.columns=GGV$info$mapObj.columns,
             dir=paste(dir,"GGVArms/",sep=""),
             fname.root=ch,
             overwriteSourcePlot=overwriteSourcePlot,
             overrideInteractive= overrideInteractive,
             header=header,
             window.size=window.size,
             image.size=image.size,
             cleanDir=cleanDir) 

                         )


        
        ########################################
        #
        # make Level 2  
        #   Sub Arms with full 
        #     interactive ability 
        #
        ########################################
 

        
        if(makeWinArms){

     
          
          if(is.na(GGV$values$plot.vec[1])) overrideInteractive = c(TRUE, FALSE, TRUE)
          if(!is.na(GGV$values$plot.vec[1])) overrideInteractive = c(TRUE, FALSE, TRUE, TRUE)
  
          
          if(locInfo$nSpot[idx] > break.num)  g.sub = rbind(g.cuts, g.cuts2)
          if(locInfo$nSpot[idx] <= break.num)  g.sub = g.cuts
          
          
          for(gs in 1:dim(g.sub)[1]){

            cat("  ..Making sub region ", gs, " of ", dim(g.sub)[1], "\n") 
            
            # get bac index off genomic locations
            fname = paste(ch, ".sub", g.sub$Region[gs], sep="")
            g.start = g.sub$startloc[gs]
            g.stop = g.sub$endloc[gs]            
            bacDX.sub = intersect(which(GGV$values$mapObj$mapping.info$g.loc.stop >= g.start), which(GGV$values$mapObj$mapping.info$g.loc.start <= g.stop))
            ttl = paste(fname, "\n ", g.start, " : ", g.stop, sep="")
            # if there is no bacs in genomic range
            # map close
            if(length(bacDX.sub) == 0){
              warning("No spots in range.  Making sub plot with next closest spots \n", immediate.=TRUE)             
              bacDX.sub = (which(GGV$values$mapObj$mapping.info$g.loc.start <= g.stop)[length(which(GGV$values$mapObj$mapping.info$g.loc.start <= g.stop))]-buffer):(which(GGV$values$mapObj$mapping.info$g.loc.stop >= g.start)[1]+buffer)
              g.start = GGV$values$mapObj$mapping.info$g.loc.start[bacDX.sub[1]]
              g.stop = GGV$values$mapObj$mapping.info$g.loc.stop[bacDX.sub[length(bacDX.sub)]]
              ttl = paste(fname, "\n No spots in orignal range \n New range: ",g.start, " : ", g.stop,sep="" )
            }
            if(length(bacDX.sub) < 2){
              warning("spot range too low \n Adjusting range with buffer \n", immediate.=TRUE)
              bacDX.sub = (bacDX.sub[1]-buffer):(bacDX.sub[length(bacDX.sub)]+buffer)
            }

            #adjust to min and max of chromosome
            if(bacDX.sub[1] < locInfo$starting[idx]) bacDX.sub[1] = locInfo$starting[idx]
            if(bacDX.sub[length(bacDX.sub)] > locInfo$ending[idx]) bacDX.sub[length(bacDX.sub)] = locInfo$ending[idx]
               
             if(locInfo$nSpot[idx] > break.num) annotation = 1:(length(annObj)-2)
            if(locInfo$nSpot[idx] <= break.num) annotation = 1:(length(annObj)-1)
            
            clrs = GGV$info$clrs



            bacDX.sub = intersect(bacDX.sub, goodDX)

       
            # common reference from other programs 
            spotDX = bacDX.sub
            
            gloc.start = GGV$values$mapObj$mapping.info$g.loc.start[spotDX[1]]
            gloc.stop = GGV$values$mapObj$mapping.info$g.loc.stop[spotDX[length(spotDX)]]
           

           
            # figure out labels
            idx.Chrom = intersect(which(GGV$values$mapObj$band.info$Chrom$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Chrom$upper >= gloc.start))
            if(gloc.start < GGV$values$mapObj$band.info$Chrom$lower[idx.Chrom[1]]) idx.Chrom = c(idx.Chrom[1]-1, idx.Chrom)
            if(gloc.stop > GGV$values$mapObj$band.info$Chrom$upper[idx.Chrom[length(idx.Chrom)]]) idx.Chrom = c(idx.Chrom, idx.Chrom[length(idx.Chrom)]+1)
            idx.Arm = intersect(which(GGV$values$mapObj$band.info$Arm$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Arm$upper >= gloc.start))
            if(gloc.start < GGV$values$mapObj$band.info$Arm$lower[idx.Arm[1]]) idx.Arm = c(idx.Arm[1]-1, idx.Arm)
            if(gloc.stop > GGV$values$mapObj$band.info$Arm$upper[idx.Arm[length(idx.Arm)]]) idx.Arm = c(idx.Arm, idx.Arm[length(idx.Arm)]+1)
            idx.Broad.Band = intersect(which(GGV$values$mapObj$band.info$Broad.Band$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Broad.Band$upper >= gloc.start))
            if(gloc.start < GGV$values$mapObj$band.info$Broad.Band$lower[idx.Broad.Band[1]]) idx.Broad.Band = c(idx.Broad.Band[1]-1, idx.Broad.Band)
            if(gloc.stop > GGV$values$mapObj$band.info$Broad.Band$upper[idx.Broad.Band[length(idx.Broad.Band)]]) idx.Broad.Band = c(idx.Broad.Band, idx.Broad.Band[length(idx.Chrom)]+1)
            idx.Fine.Band = intersect(which(GGV$values$mapObj$band.info$Fine.Band$lower <= gloc.stop), which(GGV$values$mapObj$band.info$Fine.Band$upper >= gloc.start))
            if(gloc.start < GGV$values$mapObj$band.info$Fine.Band$lower[idx.Fine.Band[1]]) idx.Fine.Band = c(idx.Fine.Band[1]-1, idx.Fine.Band)
            if(gloc.stop > GGV$values$mapObj$band.info$Fine.Band$upper[idx.Fine.Band[length(idx.Fine.Band)]]) idx.Fine.Band = c(idx.Fine.Band, idx.Fine.Band[length(idx.Chrom)]+1)
            
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
            
            if(nSpots <= GGV$info$maxLabels){
              Spots.Used = TRUE
              Chrom.Used = FALSE
           
              # get section divides
              lwrDX = idx.spots[1]
              uprDX = idx.spots[length(idx.spots)]
              sectiondiv = GGV$values$mapObj$mapping.info$g.loc.start[lwrDX]
              for(i in lwrDX:(uprDX-1)){
                sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$mapping.info$g.loc.stop[i], GGV$values$mapObj$mapping.info$g.loc.start[i+1])))
              }
              sectiondiv = c(sectiondiv, GGV$values$mapObj$mapping.info$g.loc.stop[uprDX])
              # get section centers
              sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              for(i in 2:(length(sectiondiv)-1)){
                sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
              }

              # back check if only on location
              if(lwrDX==uprDX){
                sectiondiv = c(GGV$values$mapObj$mapping.info$g.loc.start[lwrDX],GGV$values$mapObj$mapping.info$g.loc.stop[uprDX]) 
                sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              }                    
              # get sectio labels
              sectionlbls = as.character(GGV$values$mapObj$mapping.info$Spot.ID)[idx.spots]
              
            }else{
              if(nFine.Band <=GGV$info$maxLabels){
                Fine.Band.Used=TRUE
                Chrom.Used=FALSE
        
                # get section divides
                lwrDX = idx.Fine.Band[1]
                uprDX = idx.Fine.Band[length(idx.Fine.Band)]
                sectiondiv = GGV$values$mapObj$band.info$Fine.Band$lower[lwrDX]
                for(i in lwrDX:(uprDX-1)){
                  sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Fine.Band$upper[i], GGV$values$mapObj$band.info$Fine.Band$lower[i+1])))
                }
                sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Fine.Band$upper[uprDX])
                # get section centers
                sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                for(i in 2:(length(sectiondiv)-1)){
                  sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
                }

                # back check if only one location
                if(lwrDX == uprDX){
                  sectiondiv = c(GGV$values$mapObj$band.info$Fine.Band$lower[lwrDX], GGV$values$mapObj$band.info$Fine.Band$upper[uprDX])
                  sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                }
                
                # get sectio labels
                sectionlbls = as.character(GGV$values$mapObj$band.info$Fine.Band$label)[idx.Fine.Band]
                
              }else{
                if(nBroad.Band<=GGV$info$maxLabels){
                  Broad.Band.Used=TRUE
                  Chrom.Used=FALSE
        
                  # get section divides
                  lwrDX = idx.Broad.Band[1]
                  uprDX = idx.Broad.Band[length(idx.Broad.Band)]
                  sectiondiv = GGV$values$mapObj$band.info$Broad.Band$lower[lwrDX]
                  for(i in lwrDX:(uprDX-1)){
                    sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Broad.Band$upper[i], GGV$values$mapObj$band.info$Broad.Band$lower[i+1])))
                  }
                  sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Broad.Band$upper[uprDX])
                  # get section centers
                  sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                  for(i in 2:(length(sectiondiv)-1)){
                    sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
                  }
    
                  # back check if only on location
                  if(lwrDX==uprDX){
                    sectiondiv = c(GGV$values$mapObj$band.info$Broad.Band$lower[lwrDX], GGV$values$mapObj$band.info$Broad.Band$upper[uprDX])
                    sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))            
                  }
          
                  # get sectio labels
                  sectionlbls = as.character(GGV$values$mapObj$band.info$Broad.Band$label)[idx.Broad.Band]
                  
                }else{
                  if(nArm<=GGV$info$maxLabels){
                    Arm.Used=TRUE
                    Chrom.Used=FALSE
                
                    # get section divides
                    lwrDX = idx.Arm[1]
                    uprDX = idx.Arm[length(idx.Arm)]
                    sectiondiv = GGV$values$mapObj$band.info$Arm$lower[lwrDX]
                    for(i in lwrDX:(uprDX-1)){
                      sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Arm$upper[i], GGV$values$mapObj$band.info$Arm$lower[i+1])))
                    }
                    sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Arm$upper[uprDX])
                    # get section centers
                    sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                    for(i in 2:(length(sectiondiv)-1)){
                      sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
                    }
                            
                    # back check if only one location
                    if(uprDX==lwrDX){
                      sectiondiv = c(GGV$values$mapObj$band.info$Arm$lower[lwrDX],GGV$values$mapObj$band.info$Arm$upper[uprDX])
                      sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
                    }            
                    # get section labels
                    sectionlbls = as.character(GGV$values$mapObj$band.info$Arm$label)[idx.Arm]
                    
                  }
                }
              }
            }
            if(Chrom.Used){
          
              # get section divides
              lwrDX = idx.Chrom[1]
              uprDX = idx.Chrom[length(idx.Chrom)]
              sectiondiv = GGV$values$mapObj$band.info$Chrom$lower[lwrDX]
              for(i in lwrDX:(uprDX-1)){
                sectiondiv = c(sectiondiv, mean(c(GGV$values$mapObj$band.info$Chrom$upper[i], GGV$values$mapObj$band.info$Chrom$lower[i+1])))
              }
              sectiondiv = c(sectiondiv, GGV$values$mapObj$band.info$Chrom$upper[uprDX])
              # get section centers
              sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              for(i in 2:(length(sectiondiv)-1)){
                sectioncntr = c(sectioncntr, mean(c(sectiondiv[i], sectiondiv[i+1])))
              }
  
              # back check if only one location
              if(uprDX==lwrDX){
                sectiondiv = c(GGV$values$mapObj$band.info$Chrom$lower[lwrDX], GGV$values$mapObj$band.info$Chrom$upper[uprDX])
                sectioncntr = mean(c(sectiondiv[1], sectiondiv[2]))
              }                
              # get sectio labels
              sectionlbls = as.character(GGV$values$mapObj$band.info$Chrom$label)[idx.Chrom]
            }

            ylim.temp = range(mapObj$mapping.info$g.loc.center[bacDX.sub],na.rm=TRUE)
            
            # fix sectiondiv and sectioncntr
            if(sectiondiv[1] < ylim.temp[1]){
              sectiondiv[1] = ylim.temp[1]
              sectioncntr[1] = mean(c(sectiondiv[1], sectiondiv[2]))
            }
            if(sectiondiv[length(sectiondiv)] > ylim.temp[2]){
              sectiondiv[length(sectiondiv)] = ylim.temp[2]
              sectioncntr[length(sectioncntr)] = mean(c(sectiondiv[length(sectiondiv)], sectiondiv[(length(sectiondiv)-1)]))
            }
                
            #
            # add points to plot
            #

            plot.vec = eval.js(GGV$values$plot.vec)  # should be across entire genome
  
            doPlot=TRUE
            if(length(plot.vec) == 1){
              if(is.na(plot.vec))  doPlot = FALSE
            }


            if(doPlot){
              
            # if plot.vec is longer than number spots assume multiple values
            # for spots 
            if(length(eval.js(plot.vec)) > dim(mapObj$mapping.info)[1]){
              nt = length(eval.js(plot.vec))/(dim(mapObj$mapping.info)[1])
              
              new.bacDX = bacDX.sub
              for(i in 2:nt){
                new.bacDX = c(new.bacDX, (bacDX.sub + ((dim(mapObj$mapping.info)[1])*(nt -1))))
              }
              new.plot.vec = plot.vec[new.bacDX]
            }else{
              nt = 1
              new.bacDX = bacDX.sub
              new.plot.vec = plot.vec[new.bacDX]
            }
            
            
      #  side.plot.call = "image(x=0:1,y=0:1,z=matrix(rep(NA,4),ncol=2),xlim=range(eval.js(new.plot.vec),na.rm=TRUE), ylim=range(mapObj$mapping.info$g.loc.center[spot.dx],na.rm=TRUE), xlim=c(0,1), axes=FALSE, xlab='', ylab=''); axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=.5,las=2);axis(1); points(x=new.plot.vec, y=rep(mapObj$mapping.info$g.loc.center[spot.dx],nt),pch=21, bg='black') "


       #     side.plot.call = "plot(0,0,pch=' ',xlab=' ',ylab=' ',main=' ',ylim=range(sectiondiv,na.rm=TRUE),xlim=range(eval.js(new.plot.vec),na.rm=TRUE),axes=FALSE);    axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=.5,las=2);axis(1);points(x=new.plot.vec, y=rep(mapObj$mapping.info$g.loc.center[bacDX.sub],nt),pch=21, bg='black') "

            side.plot.call = "image(x=0:1,y=0:1,z=matrix(rep(NA,4),ncol=2),xlim=range(eval.js(new.plot.vec),na.rm=TRUE), ylim=range(mapObj$mapping.info$g.loc.center[bacDX.sub],na.rm=TRUE), zlim=c(0,1), axes=FALSE, xlab='', ylab=''); axis(2,at=sectiondiv,labels=NA);axis(2,tick=FALSE,at=sectioncntr,labels=sectionlbls,las=0,cex.axis=.6,las=2);axis(1); points(x=new.plot.vec, y=rep(mapObj$mapping.info$g.loc.center[bacDX.sub],nt),pch=21, bg='black')"


            
            # if(!is.na(GGV$values$plot.call[1])) side.plot.call = paste(side.plot.call, GGV$values$plot.call, sep=";")
            if(!is.na(GGV$values$side.plot.extras[1])) side.plot.call = paste(side.plot.call, GGV$values$side.plot.extras, sep=";")

          }else{
            side.plot.call=NA
            new.plot.vec=NA
            
          }
            
            
            bacDX = bacDX.sub 

           
            suppressWarnings(
            
            iGGV(vls=GGV$values$vls,
                 mapObj=GGV$values$mapObj,
                 annObj=annObj,
                 x.labels=GGV$interactive$x.labels,
                 y.labels=GGV$interactive$y.labels,
                 xy.labels=GGV$interactive$xy.labels,             
                 x.links=GGV$interactive$x.links,
                 y.links=GGV$interactive$y.links,
                 xy.links=GGV$interactive$xy.links,
                 asLinks=GGV$interactive$asLinks,
                 x.images=GGV$interactive$x.images,
                 y.images=GGV$interactive$y.images,
                 xy.images=GGV$interactive$xy.images, 
                 plot.y.index=bacDX.sub,
                 mat=GGV$values$mat,
                 maxLabels=GGV$info$maxLabels,
                 mai.mat=GGV$info$mai.mat,
                 mai.prc=GGV$info$mai.prc,
                 plot.x.index=smplDX,
                 smp.color=smp.color,
                 plot.extras=GGV$info$plot.extras,
                 smpLines=GGV$info$smpLines,
                 divCol=GGV$info$divCol,
                 plot.call= side.plot.call,
                 plot.vec= new.plot.vec,             
                 lims=GGV$info$lims,
                 annotation=annotation,
                 clrs=clrs,
                 mapObj.columns=GGV$info$mapObj.columns,
                 dir=paste(dir,"GGVwinArms/",sep=""),
                 fname.root=fname,
                 overwriteSourcePlot=overwriteSourcePlot,
                 overrideInteractive= overrideInteractive,
                 header=header,
                 window.size=window.size,
                 image.size=image.size,
                 cleanDir=cleanDir) 
               
                             )

            
          } # end for each sub region 

        
          
        }# end if withinarms
      
      }# if nSpots > 2
     
    }# if length of index is not empty
    
  } # for each chromosome arm listed
  
} # if there are new chromosome arms to make 

  
} # end function 
