
initGGV <- function(vls,
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
                    
                       chrArms=NA,
                       trackRegions=NA,
                       
                       side.plot.extras=NA, 
                       plot.vec=NA,   # should be across genome
                       plot.dx=NA,
                  
                    
                       maxLabels=25,
                       mat = NA, 
                       mai.mat = NA,
                       mai.prc=FALSE,

                       plot.extras=NA,    
                       smpLines=TRUE,
                       divCol="lightgrey",

                       lims = c(-0.5,0.5),
                       annotation = NA,
                       clrs=c("blue", "hotpink", "purple", "orange"),

                       mapObj.columns = NA,
        
                       returnVl=TRUE,
                       saveFlag=FALSE,
                       saveName="GGVobj.RData"                

                       ){


  if(is.na(chrArms[1]) & is.na(plot.vec[1])) chrArms = paste(rep(1:24, each=2),c("p","q"),sep="")


  if(length(plot.vec) == 1 ){
    if(is.null(dim(mai.mat))) mai.mat = matrix(c(0.2,  1.0, 1.0, 0.05,1.4,  1.0, 0.3, 0.05,0.2, 0.05, 1.0, 0.2), byrow=TRUE, ncol=4)
  }
    
  

  GGV = list()
  GGV$values = list(vls=vls, mapObj=mapObj, annObj=annObj, chrArms=chrArms,
                    trackRegions=trackRegions, mat=mat,
                    side.plot.extras=side.plot.extras, plot.vec=plot.vec, plot.dx = plot.dx)
  GGV$interactive = list(x.labels=x.labels, y.labels=y.labels, xy.labels=xy.labels,
                         x.links=x.links, y.links=y.links, xy.links=xy.links,
                         asLinks=asLinks, x.images=x.images, y.images=y.images, xy.images=xy.images)
  GGV$info = list(maxLabels=maxLabels, mai.mat=mai.mat, mai.prc=mai.prc,
                  plot.extras=plot.extras, smpLines=smpLines,divCol=divCol, lims=lims,
                  annotation=annotation, clrs=clrs, mapObj.columns=mapObj.columns)
  class(GGV) <- 'GGVobj'

  if(returnVl) return(GGV)
  if(saveFlag) save(GGV, file=saveName, compress=TRUE)
  

}


