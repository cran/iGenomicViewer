\name{DiseaseGenes}
\docType{data}
\alias{DiseaseGenes}

\title{Disease Gene Annotation File}

\description{
  This file is an annotation file for known disease genes. It has been
  altered slightly to include hyperlinks to the UCSC Genome Browser.
}

\usage{
data(DiseaseGenes)
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
  see also \code{annotation}, \code{anninfo}, and \code{annObj}

}

\references{
  Online Mendelian Inheritance in Man (OMIM), a knowledgebase of human
  genes and genetic disorders. Hamosh A, Scott AF, Amberger J, Bocchini
  C, Valle D, McKusick VA. Nucleic Acids Res. 2002 Jan 1;30(1):52-5

  Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM,
  Haussler D. The human genome browser at UCSC. Genome Res. 2002
  Jun;12(6):996-1006
}

\note{see also \code{annotation}, \code{anninfo}, and \code{annObj}}

\keyword{datasets}
