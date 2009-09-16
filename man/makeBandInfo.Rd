\name{makeBandInfo}
\alias{makeBandInfo}

\title{Function for creating bandinfo object}

\description{
This function creates a bandinfo object that contains genomic location
information for chromosomes, arms, broad bands, and fine bands. 
}

\usage{
makeBandInfo(file,
             chrom.levels, 
             file.sep="\t",
             autosomes=1:22,
             X.chrom = 23, 
             Y.chrom = 24,
             chr.dx = 1,
             band.dx = 4,
             start.dx = 2,
             stop.dx = 3,                         
             returnVl=TRUE,
             saveFile=FALSE,
             saveName = "BandInfo.RData",
             ...)

}

\arguments{

  \item{file}{path name of file}
  \item{chrom.levels}{vector indicating how the chrom column in the file
    is represented. i.e chr1, chrom1, 1, etc.}
  \item{file.sep}{seperation character for file}
  \item{autosomes}{numeric vector listing autosomes. defaults to human 1:22}
  \item{X.chrom}{numeric indication of x chromosome. defaults to human 23}
  \item{Y.chrom}{numeric indication of y chromosome. defaults to human 24}
  \item{chr.dx}{column in file that represents the chromosome}
  \item{band.dx}{column in file that represents the band. assumed to be
    without chromosome i.e. q11, p13.1, etc.}
  \item{start.dx}{column in file that represents the starting genomic
    location of band region}
  \item{stop.dx}{column in file that represents the stopping genomic
    location of band region}
  \item{returnVl}{logical indicating if object created should be returned}
  \item{saveFile}{logical indicating if object created should be saved to a file}
  \item{saveName}{if saveVl, path name of file}
  \item{...}{additional arguments to the read.table function for reading
  the file}

}

\details{

  This function is designed to make useful data.frames of starting and
  ending genomic locations for chromosomes, arms, broad.bands, and
  fine.bands. Each data.frame has five columns: the region (chrom, arm),
  label, lower, center, and upper. lower, center, and upper refer to
  genomic locations. label is the label that will be used for the region
  when plotting. The bandinfo object will also contain an offset
  object. This is a numeric vector of length equal to the number of
  chromsomes. It contains the numeric buffer that should be added to
  genomic locations given within chromosome to get the location with
  respect to the entire genome.  These data.frames and indicies aid in
  graphing when utilizing the \code{iGGV} function.
  
}

\value{
  
a bandinfo object that contains location information for chormosomes,
arms, broad bands, and fine bands. 
see bandinfo for more details
}

\author{Lori A. Shepherd}

\seealso{\code{\link{bandinfo}}, \code{\link{Band.Info}}, \code{\link{iGGV}}}

\examples{


#
# The following code will build the same object as data(Band.Info)
#

library("iGenomicViewer")

# writes out text files to use
writeExFiles()

# make object
band.info = makeBandInfo(file="cytoBand.txt", chrom.levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"), file.sep="\t", returnVl=TRUE, header=TRUE)


}

\keyword{methods}
