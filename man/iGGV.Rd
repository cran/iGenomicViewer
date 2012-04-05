\name{iGGV}
\alias{iGGV}

\title{Interactive Genomic Heatmap Display}

\description{
  This function creates an interactive genomic heatmap with any number
  of additional annotation tracks.  
}

\usage{

iGGV(vls,                  
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
      ps.paper="letter",
      ps.width=8,
      ps.height=11,
      cleanDir=TRUE,
      ...)
}

\arguments{
  \item{vls}{a matrix of values for the heatmap. assumes full set do not
  subset based on plot.x.index or plot.y.index, the function will
  automatically subset. The number of rows should correspond to the number
  of spot IDs in the mapObj.}
  \item{mapObj}{an object of the class  mapobj}
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
  \item{maxLabels}{maximum number of labels to appear on the heatmap y axis}
  \item{mai.mat}{n x 4 matrix of values to be passed in for each plots
    par mai. n will be 3 if plot.call is NA, and 4 if plot.calls is
    specified}
  \item{mai.prc}{logical indicating if mai mat values are percentages
    or hard coded values. If mai.proc is T, indicates percentage.}
  \item{plot.x.index}{subset of x values to use}
  \item{smp.color}{color for x-axis header values}
  \item{plot.y.index}{subset of y values to use}
  \item{goodDX}{range of acceptable y values. This is used as a
    way to removed spots of known bad quality. should correspond to
    locations in the mapObj\$mapping.info object}
  \item{genomic.start}{numeric value indicating starting genomic
    location for region to plot. Note: will round to nearest fine.band
    region}
  \item{genomic.stop}{numeric value indicating ending genomic
    location for region to plot. Note: will round to nearest fine.band
    region}
  \item{genomic.region}{character value of desired region to plot. Given
    in chromosome, arm, broad.band, or fine.band. i.e 1, 2p, 1p13.1,
    etc.}
  \item{region.type}{If genomic region is used, the type of region
    given. May either be, chrom, arm, broad.band, or fine.band.}
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
    a postscript, png, jpeg, or tiff. see \code{makeSplot} for more details}
  \item{makeInteractive}{logical, should interactive html file be
    created.see \code{makeSplot} for more details}
  \item{overrideInteractive}{Indicates which figures of the layout should be
    interactive. see \code{makeSplot} for more details}  
  \item{header}{May either be v1,v2, or v3. Determines which
    tooltip header will be in the html file. see \code{makeSplot} for more details}
  \item{window.size}{size of the html window. see \code{makeSplot} for more details}
  \item{image.size}{character indicating resize value of
    image,'width'x'height' see \code{initSplot} for more details}
  \item{ps.paper}{character indicating paper for postscript function}
  \item{ps.width}{width argument for postscript function}
  \item{ps.height}{hieght argument for postscript function}
   \item{cleanDir}{logical indicating if intermediate files created for
    mapping purposes should be deleted}
  \item{...}{extra arguments to makeImap function}
}

\details{

  This function is a wrapper to sendplot to create an interactive
  bioinformatics heatmap display. 

  The function requires a matrix of values for the heatmap, a mapping
  object, and an annotation object. For details on creating the mapping
  and annoation objects see \code{mappingObj} and \code{annotation}.  

  The matrix of values should be a complete listing, do not subset based
  on plot.x.index or plot.y.index. The function will subset
  automatically.

  The matricies x.labels, y.labels, xy.labels, x.links, y.links,
  xy.links,asLinks, x.images, y.images, and xy.images should also be
  complete listings. Do not subset based on plot.x.index or
  plot.y.index. The function will subset automatically.

  There are three ways to indicate y values that should be plotted. They
  may be specified directly through the plot.y.index. This is a numeric
  vector which corresponds to the ordering in the mapping object. They
  may be determined through giving a genomic starting and ending
  location, genomic.start and genomic stop respectively. Both starting
  and ending locations must be given if this option is utilized. The
  genomic locations should be the genomic location with respect to the
  entire genome not within a chromosome. Finally, they may be specified
  by listing a single specific region to be plotted with
  genomic.region. If this option is used, the user must also indicate
  what type of region is listed in the region.type argument. The four
  options for this argument are chrom, arm, broad.band, fine.band.  The
  region given should match up to a region in the mapping object.

  It is possible to control what is displayed from the annotation
  object so that any, all, or none of the annotation information may be
  displayed as tracks along side the heatmap. The annotation argument is
  a numeric corresponding to the order of the annotation information
  objects in the annObj. NA will display all.

}

\value{
  .png [.jpeg,.ps]  and .html graphs. .html is an interactive
  display. The display will have a heatmap, legend, annotation track, and,
  if utilized, an extra descriptive plot. 
}

\note{Works properly on linux/unix operating systems only. The below
  example show one way to use iGGV as far as plotting and
  subsetting. Please see vignette for other examples.}

\author{Lori A. Shepherd}

\seealso{\code{sendplot}, \code{\link{annotation}},
  \code{\link{mappingObj}}, \code{\link{iGGVex}}}

\examples{

library("iGenomicViewer")

data(mapping.info)
data(annObj)
data(iGGVex)


iGGV(vls = mat,
     mapObj=mapping.info,
     annObj=annObj,
     x.labels=x.lbls,
     y.labels=y.lbls,
     genomic.region="11q",
     region.type="arm",
     mapObj.columns =c(1,2,3,4))

}

\keyword{methods}
\keyword{iplot}
