\name{iGGVtiled}
\alias{iGGVtiled}

\title{Interactive Tiled Display}

\description{
  This function takes a tiled image display 'TIplot' object and creates
  an interactive layout of plots. The tiled display will be the main
  plot of the layout.
}

\usage{
iGGVtiled(TIplot,
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
          vrb=TRUE,
          ...)
}

\arguments{
  \item{TIplot}{an object of the class TIplot. This object's genomic
    locations must be with respect to the entire genome not location
    within chromosome.}
  \item{annObj}{an annotation object}
  \item{x.labels}{data frame of n x m which contains values relating to
    the x-axis values. This information is displayed in the interactive plot
    window.}
  \item{y.labels}{data frame of n x m which contains values relating to
    the y-axis values. This information is displayed in the interactive plot
    window.}
  \item{xy.labels}{list of matricies. All matricies should be of n x m.
    This information is displayed in the interactive plot window.}
  \item{x.links}{data frame of n x m which contains web addresses
    for links relating to the x-axis values. This information is displayed as
    hyperlinks in the interactive plot window.}
  \item{y.links}{data frame of n x m which contains web addresses
    for links relating to the y-axis values. This information is displayed as
    hyperlinks in the interactive plot window.}
  \item{xy.links}{list of matricies. All matricies should be of n x m.
    This information is displayed in the interactive plot window as
    hyperlinks. The values in these matricies should be complete web
    address.}
  \item{asLinks}{contains complete web address for points that should be
    treated as hyperlinks.}
  \item{x.images}{data frame of n x m which contains paths to images
    relating to the x-axis values. This information is displayed as
    images in the interactive plot window.}
  \item{y.images}{data frame of n x m which contains paths to images
    relating to the y-axis values. This information is displayed as
    images in the interactive plot window.}
  \item{xy.images}{list of matricies. All matricies should be of n x m.
    This information is displayed in the interactive plot window as
    images. The values in these matricies should be complete paths to images.}
  \item{mat}{matrix indicating layout. This argument will be passed into
    the graphics package layout call as mat.Each value in the matrix
    must be '0' or a positive integer.  If N is the largest positive
    integer in the matrix, then the integers \{1,...,N-1\} must also
    appear at least once in the matrix. '0' indicates region of no
    plotting. This may be left as NA, and a default will be used}
  \item{mai.mat}{n x 4 matrix of values to be passed in for each plots
    par mai. n will be 3 if plot.call is NA, and 4 if plot.calls is
    specified}
  \item{mai.prc}{logical indicating if mai mat values are percentages
    or hard coded values. If mai.proc is T, indicates percentage.}
  \item{plot.extras}{List of length equal to the number of
    plots: 3 if plot.call is NA, 4 if plot.call is specified.  This
    object is a list of lists. The sublists contain any additional
    plotting calls that should be executed for the plot. Each ntry must
    be a character vector. If no additional plotting is equired, an NA
    should be used}
  \item{smpLines}{logical indicating if vertical, yellow lines should be
    added between each sample of the heatmap.}  
  \item{divCol}{If smpLines, the color of the dividing lines}
  \item{plot.call}{character vector containing plot call that will be
    evaluated. This plot is added to the right of the annoation
    tracks. It is designed to add additional statistical analysis such
    as p.values, linear order statistics, etc. If NA, no plot will be
    added to the display}
  \item{plot.vec}{vector of values that will be plotted in
    plot.call. While this may be specified in the plot.call, it is also
    necessary in order to add interactive tool-tip to values}
  \item{lims}{Lower and Upper limit for vls}
  \item{annotation}{Numeric indication of which annotation information
    objects to include from the annObj. If NA all are used.}
  \item{clrs}{Character vector of colors to use for annotation tracks}
  \item{mapObj.columns}{Which columns from mapObj's mapping.info
    data.frame to include in tool-tip. May be numeric or header names}
  \item{fname.root}{Base name to use for all file created}
  \item{dir}{directory path to where files should be created}
  \item{overwriteSourcePlot}{character, should static image generated be
    a postscript, png,jpeg, or tiff. see \code{makeSplot} for more details}
  \item{makeInteractive}{logical, should interactive html file be
    created.see \code{makeSplot} for more details}
  \item{overrideInteractive}{Indicates which figures of the layout should be
    interactive. see \code{makeSplot} for more details}  
  \item{header}{May either be v1,v2, or v3. Determines which
    tooltip header will be in the html file. see \code{makeSplot} for more details}
  \item{window.size}{size of the html window. see \code{makeSplot} for more details}
  \item{image.size}{character indicating resize value of
    image,'width'x'height' see \code{initSplot} for more details}
  \item{cleanDir}{logical indicating if intermediate files created for
    mapping purposes should be deleted}
  \item{vrb}{logical indicating if status messages should be printed}
  \item{...}{extra arguments to makeImap function}


}

\details{

  This function is a wrapper to sendplot to create an interactive
  bioinformatics heatmap display. 

  The main plot of the layout is a tiled display image. The
  reasoning behind the tiled image is that there are spot.ID, BAC
  clones, etc. that span the genome. The different regions may overlap
  or not include certain regions across the entire genome. This viewer
  breaks the spot.IDs to show how they are covering a particular region.  

  This function requires an object of the class TIplot,  a mapping
  object, and an annotation object. For details on creating the mapping
  and annoation objects see \code{mappingObj} and \code{annotation}. For
  details on creating the TIplot object see \code{initTile}. The TIplot
  start and stop genomic locations must be with respect to the entire
  genome not with chromosome. If the TIplot object was created with the
  later, plese see \code{convertLoc} for help on converting these
  locations and recreate the object.  

  It is possible to control what is displayed from the annotation
  object so that any, all, or none of the annotation information may be
  displayed as tracks along side the heatmap. The annotation argument is
  a numeric corresponding to the order of the annotation information
  objects in the annObj. NA will display all. 0 will display none.


  For more details on interactive plotting capabilities, please refer to
  vignette.  
}

\value{
  .png [.jpeg,.ps]  and .html graphs. .html is an interactive
  display. The display will have a tiled heatmap, legend, annotation track, and,
  if utilized, an extra descriptive plot. 
  }

\note{

  Works properly on linux/unix operating systems only.

  The TIplot object must have Y locations with respect to entire genome
  not within chromosome. See \code{convertLoc} for conversion option.
}

\author{Lori A. Shepherd, Daniel P. Gaile}

\seealso{\code{\link{initTile}},\code{\link{makeTiled}}, \code{\link{iGGV}}}

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


# annotaiton object
# see vignette and help files for more details
data(annObj)

# make interactive plot

iGGVtiled(TIplot=TIplot,
          annObj=annObj,
          x.labels=as.data.frame(list(
                  sample.ID=paste("smp",1:TIplot$vls$nsmp,sep=""),
                  xla1=c("a","b","c","d","e","f","g","h","i","j"),
                  xla2=10:1)),
          y.labels=as.data.frame(list(
                  Spot.ID=paste("Spot",bacDX,sep=""))),
          xy.labels=list(lgr=round(Z,3)),
          mapObj.columns = c(2,3,7),
          fname.root="iGGVtiled")

   
  



}

\keyword{methods}
\keyword{iplot}
