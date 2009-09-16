\name{initGGV}

\alias{initGGV}

\title{Initialize GGV Object}

\description{
  This function initializes a GGV object to use with makeGGV
}

\usage{

initGGV(vls,
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
        plot.vec=NA,
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
        saveName="GGVobj.RData")
}

\arguments{
 \item{vls}{a matrix of values for the heatmap. assumes full set do not
  subset based on a desired sample or spot index, opportunities to
  subset will be given in executing the plots and the function will
  automatically subset. The number of rows should correspond to the number
  of spot IDs in the mapObj}
  \item{mapObj}{an object of the class mapobj}
  \item{annObj}{an annotation object}
  \item{x.labels}{data frame of n x m which contains values relating to
    the x-axis values. This information is displayed in the interactive plot
    window}
  \item{y.labels}{data frame of n x m which contains values relating to
    the y-axis values. This information is displayed in the interactive plot
    window}
  \item{xy.labels}{list of matricies. All matricies should be of n x m.
    This information is displayed in the interactive plot window}
  \item{x.links}{data frame of n x m which contains web addresses
    for links relating to the x-axis values. This information is displayed as
    hyperlinks in the interactive plot window}
  \item{y.links}{data frame of n x m which contains web addresses
    for links relating to the y-axis values. This information is displayed as
    hyperlinks in the interactive plot window}
  \item{xy.links}{list of matricies. All matricies should be of n x m.
    This information is displayed in the interactive plot window as
    hyperlinks. The values in these matricies should be complete web
    address}
  \item{asLinks}{contains complete web address for points that should be
    treated as hyperlinks}
  \item{x.images}{data frame of n x m which contains paths to images
    relating to the x-axis values. This information is displayed as
    images in the interactive plot window.}
  \item{y.images}{data frame of n x m which contains paths to images
    relating to the y-axis values. This information is displayed as
    images in the interactive plot window.}
  \item{xy.images}{list of matricies. All matricies should be of n x m.
    This information is displayed in the interactive plot window as
    images. The values in these matricies should be complete paths to images.}
  \item{chrArms}{Character vector of chromosome arms to be plotted. These
    names should match arm names in the mapObj's band.info. These
    Chromosomes will be diplayed in an index file}  
  \item{trackRegions}{A list containing character vectors for
    Broad.Band, Fine.Band, and geneName. Also may contain a n by 2 matrix
    of genomic start and stop values. These are known regions of
    interested that will be displayed as tiledImages. See makeTrack}
  \item{side.plot.extras}{character vector containing extra plotting calls that will be
    evaluated on the optional side plot. This plot is added to the right of the annoation
    tracks. It is designed to add additional statistical analysis such
    as p.values, linear order statistics, etc. }
  \item{plot.vec}{vector of values that will be plotted in the side
    plot. If NA, no plot will be added to the display. These values are the
    additional plots x-values. See details}
  \item{plot.dx}{Index to subset plot.vec on when creating
    an initial genomic plot. If NA, all values are used}
  \item{maxLabels}{maximum number of labels to appear on the y axis}
  \item{mat}{matrix indicating layout. This argument will be passed into
    the graphics package layout call as mat.Each value in the matrix
    must be '0' or a positive integer.  If N is the largest positive
    integer in the matrix, then the integers \{1,...,N-1\} must also
    appear at least once in the matrix. '0' indicates region of no
    plotting. This may be left as NA, and a default will be used. This
    matrix will be used for Chromosome Arm and Sub.Arm Plots}
  \item{mai.mat}{n x 4 matrix of values to be passed in for each plots
    par mai. n will be 3 if plot.vec is NA, and 4 if plot.vec is
    specified. This will be used for Chromosome Arm and Sub.Arm plots}
  \item{mai.prc}{logical indicating if mai mat values are percentages
    or hard coded values. If mai.prc is T, indicates percentage. This
    will be used for Chromosome Arm and Sub.Arm plots}
  \item{plot.extras}{List of length equal to the number of
    plots: 3 if plot.vec is NA, 4 if plot.vec is specified.  This
    object is a list of lists. The sublists contain any additional
    plotting calls that should be executed for the plot. Each entry must
    be a character vector. If no additional plotting is equired, an NA
    should be used}
  \item{smpLines}{logical indicating if vertical lines should be
    added between each sample of the heatmap}
  \item{divCol}{If smpLines, the color of the dividing lines}
  \item{lims}{Lower and Upper limit for vls}
  \item{annotation}{Numeric indication of which annotation information
    objects to include from the annObj. If NA all are used}
  \item{clrs}{Character vector of colors to use for annotation tracks}
  \item{mapObj.columns}{Which columns from mapObj's mapping.info
    data.frame to include in tool-tip. May be numeric or header names}
  \item{returnVl}{Should GGV object be returned}
  \item{saveFlag}{Should GGV object be saved}
  \item{saveName}{If saveFlag, path file name to save object}

}

\details{
  trackRegions list known regions of interest that should show more
  descriptive plots. Regions may be specified by Broad.Band, Fine.Band,
  genomic location, and geneName. See makeTrack for more details


 
  plot.vec will be subset based on plot.dx to make the genomic
  plot in the makeGGV function. It is important the the plot.vec are
  across the genome so the additional plot may be displayed for all
  chromosome arms. The vector contains x-values to plot. The y-values
  will be determined by plot.dx, or internal subsetting. If plot.vec is
  NA no side plot will be included

  side.plot.extras is a character string. It is similar to the plot.extras
  argument only it is a single variable. These calls will be evaluated
  for the additional side plot, if utilized. The calls should be
  separated by a semicolon. 
  

  See vignette for more details

}

\value{If returnVl, an object of the class 'GGVobj', a Genomic Viewer
  Object, is returned}

\note{see vignette for more details}

\author{Lori A. Shepherd}

\seealso{\code{\link{GGVobj}}, \code{\link{iGGV}}, \code{\link{makeTrack}}}

\examples{

# load library
library("iGenomicViewer")

# load example data
data(iGGVex)
data(mapping.info)
data(annObj)


mapObj = mapping.info
y.lbls$Pdisc = round(y.lbls$Pdisc,3)
chrArm = c("8p", "18p")
pvls = rep(rep(rep(1:4, each=5), 960))[1:length(mapObj$mapping.info$g.loc.center)]
plot.vec = pvls[1:length(mapObj$mapping.info$g.loc.center)]
side.plot.extras="title(main='pvals')"


# make trackregion object
trackRegion = makeTrack(Fine.Band = c("8p11.22","18p11.21"),genomicLoc = NA, geneName = "FANCE")

# initialize GGV object
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



}

\keyword{methods}
