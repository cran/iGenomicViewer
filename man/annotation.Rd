\name{annotation}

\alias{annotation}
\alias{annotationObj}
\alias{makeAnnotation}

\title{Functions for Creating the Annotation Object}

\description{These functions are designed to create an annotation
  object.}

\usage{

makeAnnotation(file,
               label,
               chrom,
               chrom.levels,
               band.info=NA,
               loc=NA,
               loc.start=NA,
               loc.stop=NA,
               file.sep="\t",
               additional=NA,
               names.additional = NA,
               links=NA,
               names.links=NA,
               images=NA,
               names.images=NA,
               returnVl = TRUE,
               saveVl = FALSE,
               saveName="Annotation.RData",
               ...)


annotationObj(annotation,  
              annotationObj = NA, 
              obj.name = NA,
              returnVl = TRUE,
              saveVl = FALSE,
              saveName="AnnotationObj.RData")

}

\arguments{
  
\item{file}{path name of file}
\item{label}{column in file that represents the label for annotation region}
\item{chrom}{column in file that represents the chromosome location}
\item{chrom.levels}{vector indicating how the chrom column is
  represented. i.e chr1, chrom1, 1, etc.}
\item{band.info}{a band.info object. see makeBandInfo}
\item{loc}{column in file representing genomic location. See details for
more information}
\item{loc.start}{column in file representing starting genomic location. See details for
more information}
\item{loc.stop}{column in file representing ending genomic location. See details for
more information}
\item{file.sep}{seperation character for file}
\item{additional}{indication of additional columns in file that should
  be included in annotation object}
\item{names.additional}{optional header names for additional columns to be included}
\item{links}{indication of columns in file that represent hyperlinks,
  character or numeric, or a data.frame of hyperlinks, with rows equal
  to the number of rows in the file, and in the same order as they
  appear in the file}
\item{names.links}{optional header names for links}
\item{images}{indication of columns in file that represent images,
  character or numeric, or a data.frame, with rows equal to the number of
  rows of the file, of images that must be in the same order as they
  appear in the file}
\item{names.images}{optional header names for images}
\item{returnVl}{logical indicating if object created should be returned}
\item{saveVl}{logical indicating if object created should be saved to a file}
\item{saveName}{if saveVl, path name of file}
\item{annotation}{annotation object created by the makeAnnotation function}
\item{annotationObj}{annotation object created by annotationObj
  function. If NA, creates a new object.}
\item{obj.name}{name of the annotation set in the annotationObj}
\item{...}{additional arguments to the read.table function for reading
  the file}

}



\details{
  
  The annotationObj is a larger object that contains all individual
  annotation information desired for a genomic mapping. These
  annotations will be represented in a track along side of the main
  heatmap of the \code{iGGV} function. Each individual annotation
  information object must be added separatly to the main annotationObj.    
  
  The annotation file must minimally contain columns for name,
  chromosome, and genomic location. The genomic location may be given by
  a single value assumed to be a central location and should be
  specified using the loc argument. This value should be the genomic
  location with respect to the entire genome not the location within a
  chromosome. The other, recommended, way of specifing the genomic
  location is with both a starting and stopping genomic location, and
  should be specified using loc.start and loc.stop arguments. These
  values are the location within the chromosome not over the entire
  genome. If loc.start and loc.stop are used, loc should be NA.   
     
  The band.info argument is an object of the class band.info. This object contains
  information about chromosome, arm, broad.band, and fine.band
  regions. See band.info help files for details. 

  The \code{iGGV} function has the capabilities of including hyperlinks
  in displayed tool-tips. In order for the links to display properly
  they need to be specified in the annotation object. If there are links
  in the file, the column number or head name may be specified as
  links. The user may add their own hyperlinks by assigning the links
  argument to a table of the dimension n by m where n is equivalent to
  the number of items in the original file. The links must be in proper
  format of a valid web address, i.e http://completewebpathhere.com .  
  Similarly, images may also be displayed and set through the images
  argument. 

  
}



\value{
annotationObj returns an object containing annotation information
objects. 

makeAnnotation returns an object of the class anninfo. This object
contains annotation information such as labels, chromosome locations,
genomic locations, and web links. 
  
}

\note{For an example object see \code{annObj}}

\author{Lori A. Shepherd}

\seealso{\code{\link{anninfo}}, \code{\link{annObj}},\code{\link{bandinfo}}, \code{\link{iGGV}}}

\examples{

library("iGenomicViewer")

#
# The following code will build the same object as data(annObj)
#


# writes out text files to use
writeExFiles()

# loads a band.info object
data(Band.Info)

# makes anninfo object for cancerGenes
annotation1 = makeAnnotation(file="CancerGenes.txt", file.sep="\t", label=2, chrom=3,chrom.levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),band.info=band.info,  loc=NA, loc.start=4, loc.stop=5, additional=0, header=TRUE, links=6)

# makes anninfo object for DiseaseGenes
annotation2 = makeAnnotation(file="DiseaseGenes.txt", file.sep="\t", label=2, chrom=3,chrom.levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),band.info=band.info, loc=NA, loc.start=4, loc.stop=5, additional=0, header=TRUE, links=7)


# makes anninfo object for DNArepairGenes
annotation3 = makeAnnotation(file="DNArepairgenes.txt", file.sep="\t", label=1, chrom=2, chrom.levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"), band.info=band.info, loc=NA, loc.start=3, loc.stop=4, additional=0, header=TRUE, links=5)

# take individual anninfo objects and create a master annotation object 
annObj = annotationObj(annotation1, obj.name="CancerGenes")
annObj = annotationObj(annotation2, annotationObj=annObj, obj.name="DiseaseGenes")
annObj = annotationObj(annotation3, annotationObj=annObj, obj.name="DNArepairGenes")



}

\keyword{methods}
