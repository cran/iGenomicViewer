\name{iGGVex}
\docType{data}
\alias{iGGVex}
\alias{mat}
\alias{x.lbls}
\alias{y.lbls}

\title{Example Dataset for iGGV Function}

\description{
  This dataset contains an example data matrix of values that
  corresponds to the default mapping object. It also contains mock
  description data for the spot.ID and for samples, y.lbls and x.lbls
  respectively. 
}

\usage{
data(iGGVex)
}

\format{
  mat is a matrix of values. The number of rows should correspond to the
  number of spots in the mapping object. It this case the default
  mapping object is used. see \code{mapping.info}

  y.lbls is a matrix with descriptive attributes for the spots. The
  number of rows should correspond to the number of spots in the mapping
  object. It this case the default mapping object is used. see
  \code{mapping.info}. Do not subset based on the desired range; the
  program assumes a complete table and will automatically subset.

  x.lbls is a matrix with descriptive attributes for the samples. The
  number of rows should be equal to the number of samples, or the number
  of columns in mat.  Do not subset based on the desired range; the
  program assumes a complete table and will automatically subset.
}

\note{mat matrix corresponds to default mapping object see
  \code{mapping.info}. mock data}

\keyword{datasets}
