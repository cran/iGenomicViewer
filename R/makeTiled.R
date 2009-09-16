

makeTiled <- function(TIplot, smpDiv=TRUE, divCol="lightgrey"){


 for(h in 1:TIplot$vls$H){

   tmat=array(NA,dim=c(length(TIplot$vls$Xcoords),
                       length(TIplot$tractBound[[h]]$tbound)-1))
   
   if(h==1){
     image(x=TIplot$vls$Xcoords,y=TIplot$tractBound[[h]]$tbound,z=tmat,
           xlim=TIplot$lims$xlim,ylim=TIplot$lims$ylim,zlim=TIplot$lims$zlim,
           col=TIplot$Zcol,axes=FALSE,xlab=TIplot$labels$xlab,
           ylab=TIplot$labels$ylab, main=TIplot$labels$ttl)
     #axis(1,at=1:J-0.5,las=2,1:J)
   }
     
   tmat[TIplot$tractBound[[h]]$xdx,TIplot$tractBound[[h]]$Ymap]=t(TIplot$vls$Z[TIplot$tractBound[[h]]$ydx,])

   image(x=TIplot$vls$Xcoords,y=TIplot$tractBound[[h]]$tbound,z=tmat,add=TRUE,
         zlim=TIplot$lims$zlim,col=TIplot$Zcol)


   if(!is.na(TIplot$labels$xlabels[1])) axis(1,at=1:TIplot$vls$nsmp-0.5,las=2,labels=TIplot$labels$xlabels, cex.axis=TIplot$cex$xcex)
 
   
 }

 if(!is.na(TIplot$labels$ylabels[1])) axis(2, at=((TIplot$vls$Y[,1]+TIplot$vls$Y[,2])/2), labels=TIplot$labels$ylabels,
       las=2, cex.axis=TIplot$cex$ycex )
  
 if(smpDiv)abline(v=0:TIplot$vls$nsmp,col="white",lwd=2)
 if(smpDiv)abline(v=0:TIplot$vls$nsmp,col=divCol,lwd=1)

  

}



