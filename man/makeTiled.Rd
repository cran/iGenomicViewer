\name{makeTiled}
\alias{makeTile}
\alias{makeTiled}


\title{Makes Tiled Image Display Plot}

\description{
  This function makes a tiled image display plot from a tiled image
  display 'TIplot' object.
}

\usage{
makeTiled(TIplot,
          smpDiv=TRUE,
          divCol="lightgrey")
}

\arguments{
\item{TIplot}{tiled display image object, of the class 'TIplot'}
\item{smpDiv}{Should vertical dividing lines be placed in between
  samples of plot}
\item{divCol}{If smpDiv, the color of the dividing line}
}

\details{
  This function take the TIplot object and creates a static tiled
  display image. All information for plotting is stored in the TIplot
  object. The user may choose to display lines separating samples. This
  and the color of the separation lines are controlled by smpDiv and
  divCol. 
  
}

\value{Generates a static tiled display plot}


\author{Lori A. Shepherd, Daniel P. Gaile}

\seealso{\code{\link{initTile}}, \code{\link{iGGVtiled}}}

\examples{

# load data objects
library("iGenomicViewer")
data(iGGVex)
data(mapping.info)

# create a subset range
bacDX = 103:112
smplDX = 1:10
Z = mat[bacDX,smplDX]

# make object
TIplot = initTile(Z=Z,
                  bacDX=bacDX, mapObj=mapping.info,
                  H=3,zlims=c(-0.5,0.5),
                  ylabels=paste("Spot",bacDX, sep=""),
                  xlabels=paste("smp",smplDX, sep=""),
                  xlab="Samples",
                  ylab="SpotID",
                  ttl="tiledImage",
                  returnVl=TRUE,
                  saveFlag=FALSE,
                  saveName="TIplot.RData")


# plot object
#      makeTiled(TIplot)



}

\keyword{methods}
