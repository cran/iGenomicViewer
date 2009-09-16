\name{makeGGV}
\alias{makeGGV}

\title{Creates a Series of Interactive Plots}

\description{
  This functions acts on a GGV 'generic genomic viewer' object to create
  a series of interactive plots.
}

\usage{
makeGGV(GGV,
        goodDX=NA,
        smplDX=NA,
        smp.color=NA,
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
        cleanDir = TRUE)
}

\arguments{
  \item{GGV}{a GGVobj, 'generic genomic viewer' object }
  \item{goodDX}{range of acceptable y values. This is used as a
    way to removed spots of known bad quality. should correspond to
    locations in the mapObj\$mapping.info object}
  \item{smplDX}{index or ordering for x axis}
  \item{smp.color}{colors for x axis}
  \item{tileNum}{For tiled image heatmaps, the number of tracks to split
    spot.IDs}
  \item{buffer}{An additional number of spots to plot surrounding know
    regions.}
  \item{makeWinArms}{Should within Arm sub plots be generated. This
    creates all chromosome sub plots for chromosomes listed}
  \item{tiledMat}{matrix indicating layout. This argument will be passed into
    the graphics package layout call as mat.Each value in the matrix
    must be '0' or a positive integer.  If N is the largest positive
    integer in the matrix, then the integers \{1,...,N-1\} must also
    appear at least once in the matrix. '0' indicates region of no
    plotting. This may be left as NA, and a default will be used. This
    matrix will be used for tiled images only. mat for chromosome arms
    and sub arms are specified in initGGV}
  \item{tiledMai.mat}{n x 4 matrix of values to be passed in for each plots
    par mai. n will be 3 if plot.call is NA, and 4 if plot.calls is
    specified. This will be used for tiled images only. mai.mat for
    chromosome arms and sub arms are specified in initGGV}
  \item{tiledMai.prc}{logical indicating if mai mat values are percentages
    or hard coded values. If mai.prc is T, indicates percentage. This
    will be used for tiled images only. mai.prc for chromosome arms and
    sub arms are specified in initGGV}
  \item{fname.root}{Base name to use for index file and genomic plot if
    applicable.}
  \item{dir}{directory path to where files and sub directories should be created}
  \item{overwriteSourcePlot}{character, should static image generated be
    a postscript, png, jpeg, or tiff. see \code{makeSplot} for more details}
  \item{header}{May either be v1,v2, or v3. Determines which
    tooltip header will be in the html file. see \code{makeSplot} for more details}
  \item{window.size}{size of the html window for chromosome arm and sub
    chromosome arm plots. see \code{makeSplot} for more details}
  \item{image.size}{character indicating resize value of
    image,'width'x'height'  for chromosome arms and sub chromosome arm
    plots. see \code{initSplot} for more details}
  \item{tiled.window.size}{size of the html window for tiled image
    plots.see \code{makeSplot} for more details}
  \item{tiled.image.size}{character indicating resize value of
    image,'width'x'height' tiled image plots. see \code{initSplot} for
    more details}
  \item{cleanDir}{logical indicating if temporary files should be
    removed}
}

\details{

  see vignette for better details

  This function creates a series of plots based on user settings. The
  top layer is a chromosome arm index page with listed arms in
  chrArms. If a genomic plot is specified in the
  GGVobj\$values\$plot.call, a genomic plot is also generated for the
  given subregion (based on GGVobj\$values\$plot.dx,
  GGVobj\$values\$plot.vec). While the initial plot will be plotted over
  the given plot.dx, it is important that the plot.vec and plot.call
  arguments are across the genome so the plot may be included on any
  chromosome arm plots specified. The next set of plots are the
  chromosome arm plots with limit interactivity. These plots will be
  generated for all chromosome arms listed in chrArm or utilized in the
  genomic plot plot.dx. The main heatmap is not interactive in these
  plots, only the annotation tracks and, if applicable, the genomic plot
  are interactive. Additional annotation tracks for known regions of
  interest as well as sub chromosome arm regions are added. The known
  regions are tiled image displays; these regions were specified in
  GGV\$values\$trackRegions. see makeTrack. If the user clicks on one
  of these regions the assocaited plot will appear. The sub chromosome
  arm regions are fully interactive plots. These plots are only
  generated if makeWinArms is TRUE. If makeWinArms is FALSE the links in
  these regions will not work properly.  

}

\value{A series of interactive plots are generated within the given
  directory, dir}

\note{see vignette for more details}

\author{Lori A. Shepherd}

\seealso{\code{\link{GGVobj}},\code{\link{initGGV}}}

\examples{

# see vignette for example


}

\keyword{methods}
