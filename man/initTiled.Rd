\name{initTile}
\alias{initTile}
\alias{initTiled}


\title{Creates A Tiled Display Image 'TIplot' Object}

\description{The initTile function creates a tiled display image
  'TIplot' object. A TIplot object holds all necessary elements to make
  a tiled display plot.}

\usage{
initTile(Z,
         bacDX,
         goodDX=NA,
         mapObj=NA,
         H=2,
         zlims=c(-0.5,0.5),
         smplDX=NA,
         ylabels=NA,
         xlabels=NA,
         x.axis.cex =0.5,
         y.axis.cex =0.5,
         xlab="samples",
         ylab="BAC location",
         ttl=NA,
         returnVl=TRUE,
         saveFlag=FALSE,
         saveName="TIplot.RData")

}

\arguments{
\item{Z}{matrix of values for image. The number of rows and columns should be equal
  to the lenghts of bacDX and smplDX. If the matrix is larger the matrix
  will be subset based on bacDX and smplDX }
\item{bacDX}{index of spot.ID's to graph. should correspond to index in
  mapObj. This will be used to determine the genomic start and stopping
  location of the plot}
\item{goodDX}{range of acceptable y values. This is used as a
  way to removed spots of known bad quality. should correspond to
  locations in the mapObj\$mapping.info object}
\item{mapObj}{an object of the class  mapobj. This object should be made
  with genomic starting and stopping locations not a central location}
\item{H}{number of tracks to display per sample}
\item{zlims}{minimum and maximum range for values in Z}
\item{smplDX}{subset for x axis. If the second dimension of Z is larger
  than or equal to the length of smplDX, Z will be subset based on this
  index}
\item{ylabels}{vector indicating labels for Y axis. Should be equal in
  length to the number of rows in Z [or Y]}
\item{xlabels}{vector indicating labels for X axis. Should be equal in
  length to the number of columns in Z}
\item{x.axis.cex}{display size of xlabels}
\item{y.axis.cex}{display size of ylabels}
\item{xlab}{main x axis label for plot}
\item{ylab}{main y axis label for plot}
\item{ttl}{main title for plot}
\item{returnVl}{Should TIplot object be returned}
\item{saveFlag}{Should TIplot object be saved}
\item{saveName}{If saveFlag, path file name to save object}
}

\details{

  The initTile function creates a tiled display image 'TIplot' object. A
  TIplot object holds all necessary elements to make a tiled display
  plot. The function takes in a matrix of values to display as a tiled
  image, Z,and a bacDX indicating the spot.IDs that should be used and
  correspond to the rows of Z. If Z has more rows than the length of
  bacDX, it is assumed that Z should be subset based on bacDX. The bacDX
  should refer to spot.IDs given in a mapObj. If no mapObj is used, the
  default package mapObj is used. It is also possible to subset Z based
  on a smplDX. If no smplDX is given it is assumed that all columns of Z
  should be used. The reasoning behind the tiled image is that there are
  spot.ID, BAC clones, etc. that span the genome. The different regions
  may overlap or not include certain regions across the entire
  genome. This viewer breaks the spot.IDs to show how they are covering
  a particular region.
  
}

\value{If returnVl, an object of the class 'TIplot'}

\note{Assumes genomic location values on Y axis and samples on X axis}

\author{Lori A. Shepherd, Daniel P. Gaile}

\seealso{\code{\link{convertLoc}}, \code{\link{TIplot}},
  \code{\link{makeTile}}, \code{\link{iGGVtiled}}}

\examples{

library("iGenomicViewer")

# load data objects
data(iGGVex)
data(mapping.info)
mapObj = mapping.info

# create a subset range
bacDX = 103:112
smplDX = 1:10
Z = mat[bacDX,smplDX]

# make object
TIplot = initTile(Z=Z,
                  bacDX=bacDX, mapObj=mapObj,
                  H=3,zlims=c(-0.5,0.5),
                  ylabels=paste("Spot",bacDX, sep=""),
                  xlabels=paste("smp",smplDX, sep=""),
                  xlab="Samples",
                  ylab="SpotID",
                  ttl="tiledImage",
                  returnVl=TRUE,
                  saveFlag=FALSE,
                  saveName="TIplot.RData")


}

\keyword{methods}
