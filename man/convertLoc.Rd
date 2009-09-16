\name{convertLoc}
\alias{convertGloc}
\alias{convertCloc}
\alias{convertLoc}

\title{Converts Genomic Locations}

\description{
  This function will convert genomic locations with respect to
  chromosome to genomic locations with respect to the entire genome or
  vice versa. 
}

\usage{
convertCloc(x,
            chr,
            row=TRUE,
            bandobj=NA)

convertGloc(x,
            chr,
            row=TRUE,
            bandobj=NA)
}

\arguments{
\item{x}{numeric vector or matrix of values}
\item{chr}{indication of chromosome}
\item{row}{when x is a matrix, should the chr value be applied by row or
by column}
\item{bandobj}{an object of the class bandinfo. if NA, useds
  iGenomicViewer's default band.info object}

}

\details{

  Some applications require genomic locations to be with respect to the
  genome instead of by chromosome and vice versa. These functions
  provide a way to convert a numeric vector or matrix of one type to the
  other.

  chr indicates which chromosome the value is on. If x is a vector, chr
  should be a vector of equal length, or a single value indicating all
  values of x are on the same chromosome. If x is a matrix, the values
  may be applied by row or by column. If by row, chr should be equal to
  dim(x)[1] and if by column, chr should be equal to dim(x)[2]. A single
  value may be used if all values are on the same chromosome.

  The offset is calculated through the band.info object. See band.info
  for more information. If a bandobj is not specified the package's
  default bandinfo object is used. 

}

\value{vector or matrix, of same form as input, of converted values}

\author{Lori A. Shepherd}

\seealso{\code{\link{Band.Info}}}

\examples{

#
# this is a mock example
#

x = matrix(0,24,3)
x[,2] = 1
x[,3] = 2

# this is x before
x

chr = 1:24

newx = convertGloc(x=x, chr=chr, row=TRUE, bandobj=NA)

# this is x after
newx

}

\keyword{methods}
