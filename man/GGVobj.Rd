\name{GGVobj-class}
\docType{class}

\alias{GGVobj}

\title{GGV 'generic genomic viewer' class}

\description{
  The generic genomic viewer object contains all necessary information
  to plot a series of interactive genomic plots. 
}

\section{GGVobj structure}{
  \describe{
  \item{values}{Contains values, objects, and matricies with main data
    source}
  \item{interactive}{contains matricies and vectors of interactive data
    for the plots}
  \item{info}{contains values of other arguments used in functions call}
  }
}

\section{GGVobj's values structure}{
  \describe{

  \item{vls}{a matrix of values for the heatmap. assumes full set do not
  subset based on plot.x.index or plot.y.index, the function will
  automatically subset. The number of rows should correspond to the number
  of spot IDs in the mapObj}
  \item{mapObj}{an object of the class mapobj}
  \item{annObj}{an annotation object}
  \item{chrArms}{Character vector of chromosome arms to be plotted. These
    names should match arm names in the mapObj's band.info. These
    Chromosomes will be diplayed in an index file}
  \item{trackRegions}{A list containing character vectors for
    Broad.Band, Fine.Band, and geneName. Also may contain a n by 2 matrix
    of genomic start and stop values. These are known regions of
    interested that will be displayed as tiledImages. See makeTrack}
  \item{mat}{matrix indicating layout. This argument will be passed into
    the graphics package layout call as mat.Each value in the matrix
    must be '0' or a positive integer.  If N is the largest positive
    integer in the matrix, then the integers {1,...,N-1} must also
    appear at least once in the matrix. '0' indicates region of no
    plotting. This may be left as NA, and a default will be used. This
    matrix will be used for Chromosome Arm and Sub.Arm Plots}
   \item{plot.call}{character vector containing extra plotting calls that will be
    evaluated on the optional side plot. This plot is added to the right of the annoation
    tracks. It is designed to add additional statistical analysis such
    as p.values, linear order statistics, etc. }
  \item{plot.vec}{vector of values that will be plotted in the side
    plot. If NA, no plot will be added to the display. These values are the
    additional plots x-values. See details}
  \item{plot.dx}{Index to subset plot.vec on when creating an initial
  genomic plot. If NA, all values are used}
  }
 
}

\section{GGVobj's interactive structure}{
  \describe{
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
    images in the interactive plot window}
  \item{y.images}{data frame of n x m which contains paths to images
    relating to the y-axis values. This information is displayed as
    images in the interactive plot window}
  \item{xy.images}{list of matricies. All matricies should be of n x m.
    This information is displayed in the interactive plot window as
    images. The values in these matricies should be complete paths to
    images.}
  }
}

\section{GGVobj's info structure}{
  \describe{
  \item{maxLabels}{maximum number of labels to appear on the heatmap y axis}
  \item{mai.mat}{n x 4 matrix of values to be passed in for each plots
    par mai. n will be 3 if plot.call is NA, and 4 if plot.calls is
    specified. This will be used for Chromosome Arm and Sub.Arm plots}
  \item{mai.prc}{logical indicating if mai mat values are percentages
    or hard coded values. If mai.prc is T, indicates percentage. This
    will be used for Chromosome Arm and Sub.Arm plots}
  \item{plot.extras}{List of length equal to the number of
    plots: 3 if plot.call is NA, 4 if plot.call is specified.  This
    object is a list of lists. The sublists contain any additional
    plotting calls that should be executed for the plot. Each entry must
    be a character vector. If no additional plotting is equired, an NA
    should be used}
  \item{smpLines}{logical indicating if vertical  lines should be
    added between each sample of the heatmap}
  \item{divCol}{If smpLines, the color of the dividing lines}
  \item{lims}{Lower and Upper limit for vls}
  \item{annotation}{Numeric indication of which annotation information
    objects to include from the annObj. If NA all are used}
  \item{clrs}{Character vector of colors to use for annotation track}
  \item{mapObj.columns}{Which columns from mapObj's mapping.info
    data.frame to include in tool-tip. May be numeric or header names}
  }
}



\note{Users are encouraged to view the vignette for more detail}

\author{Lori A. Shepherd}

\seealso{\code{\link{makeGGV}}, \code{\link{initGGV}}, \code{\link{updateGGV}}}

\keyword{classes}
