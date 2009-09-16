\name{CancerGenes}
\docType{data}
\alias{CancerGenes}

\title{Cancer Gene Annotation File}

\description{
  This file is an annotation file for known cancer genes. It has been
  altered slightly to include hyperlinks to the UCSC Genome Browser.
}

\usage{
data(CancerGenes)
}

\format{
  tab delimited text file with columns for gene name, chromosome,
  start location, end location, and weblink to UCSC Genome Browser.  
}

\details{
  annotation files must minimally have columns for gene, chromosome, start
  location, and end location.
  This file was altered to also include weblinks to the human genome
  browser at UCSC.
  see also \code{annotation} and \code{anninfo}

}

\references{
  Futreal PA, Coin L, Marshall M, Down T, Hubbard T, Wooster R, Rahman
  N, Stratton MR. A census of human cancer genes. Nat Rev Cancer. 2004
  Mar;4(3):177-83

  Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM,
  Haussler D. The human genome browser at UCSC. Genome Res. 2002
  Jun;12(6):996-1006
}

\note{see also \code{annotation}, \code{anninfo}, \code{annObj}}

\keyword{datasets}
