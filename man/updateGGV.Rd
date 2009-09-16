\name{updateGGV}
\alias{updateGGV}

\title{Updates a 'GGV' object with new known track regions}

\description{This fucntion will update a GGV object with new known
  regions of interest.}

\usage{
updateGGV(GGV,
          trackRegions,
          appendTo=TRUE,
          returnVl=TRUE,
          saveFlag=FALSE,
          saveName="GGVobj.RData")
}

\arguments{
  \item{GGV}{A GGV object. See initGGV}
  \item{trackRegions}{A object made by or with the same structure as the
    makeTrack function. See makeTrack}
  \item{appendTo}{Logical indicating if the regions should be added to
    the existing trackRegion object or if the new regions should totally
    replace the GGV object trackRegions}
  \item{returnVl}{Should GGV object be returned}
  \item{saveFlag}{Should GGV object be saved}
  \item{saveName}{If saveFlag, path file name to save object}

}

\details{
  This function acts on an exisiting GGV object. The user can either
  replace the known regions in the GGV with new regions, or append new
  regions to the existing GGV trackRegion object.

  For more information on trackRegion see makeTrack help file or package
  vignette. 
  
}

\value{If returnVl, a new updated GGV object is returned.}

\note{see vignette for more details}

\author{Lori A. Shepherd}

\seealso{\code{\link{initGGV}}, \code{\link{makeTrack}}, \code{\link{GGVobj}}}

\examples{

# load library
library("iGenomicViewer")

# load example data
data(iGGVex)
data(mapping.info)
data(annObj)


mapObj = mapping.info
chrArm = c("8p", "18p")
y.lbls$Pdisc = round(y.lbls$Pdisc,3)
pvls = rep(rep(rep(1:4, each=5), 960))[1:length(mapObj$mapping.info$g.loc.center)]
plot.vec = pvls[1:length(mapObj$mapping.info$g.loc.center)]
side.plot.extras="title(main='pvals'"

# make track Region
trackRegion = makeTrack(Fine.Band = c("8p11.22","18p11.21"),genomicLoc = NA, geneName = "FANCE")


# make GGV object
GGV = initGGV(vls = round(mat, 3),
              mapObj = mapObj,
              annObj = annObj,
              x.labels=x.lbls,
              y.labels=y.lbls,
              xy.labels=NA,
              chrArms=chrArm,
              trackRegions=trackRegion,
              side.plot.extras=side.plot.extras,
              plot.vec=plot.vec,
              plot.dx=which(mapObj$mapping.info$Chrom=="chr8"),
              mapObj.columns =c(2,8,10),
              smpLines=TRUE,
              divCol="lightgrey")


# make new track region
trackRegionUpdate = makeTrack(geneName="MYC")


# this will append new regions to old object regions
GGVnew = updateGGV(GGV=GGV,
                   trackRegions=trackRegionUpdate,
                   appendTo=TRUE)



# this will replace GGV object's trackregions with
# new object
GGVnew = updateGGV(GGV=GGV,
                   trackRegions=trackRegionUpdate,
                   appendTo=FALSE)



}

\keyword{methods}
