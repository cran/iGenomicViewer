\name{mappingObj}
\alias{mappingObj}
\alias{mappingObjMarray}
\alias{mappingObjADF}
\alias{mappingObjDF}

\title{Creates Mapping Object}

\description{This function creates a mapping object}

\usage{
mappingObj(file,
          spot.ID,
          chrom,
          chrom.levels,
          loc=NA,
          loc.start=NA,
          loc.stop=NA,
          file.sep="\t",
          additional=NA,
          names.additional=NA,
          links=NA,
          names.links=NA,
          images=NA,
          names.images=NA,
          band.info = NA,
          returnVl = TRUE,
          saveFile = FALSE,
          saveName="MapObj.RData",
          ...)

mappingObjMarray(obj,
                 spot.ID=NA,
                 chrom,
                 locBy,
                 base.chrm=NA,
                 reg.exp=NA,           
                 loc=NA,
                 loc.start=NA,
                 loc.stop=NA,
                 additional=NA,
                 names.additional=NA,
                 links=NA,
                 names.links=NA,
                 images=NA,
                 names.images=NA,
                 band.info = NA,
                 returnVl = TRUE,
                 saveFile = FALSE,
                 saveName="MapObj.RData")
 
mappingObjADF(adf,
              spot.ID=NA,
              chrom,
              locBy,
              base.chrm=NA,
              reg.exp=NA,           
              loc=NA,
              loc.start=NA,
              loc.stop=NA,
              additional=NA,
              names.additional=NA,
              links=NA,
              names.links=NA,
              images=NA,
              names.images=NA,
              band.info = NA,
              returnVl = TRUE,
              saveFile = FALSE,
              saveName="MapObj.RData")



mappingObjDF(df,
             spot.ID=NA,
             chrom,
             locBy,
             base.chrm=NA,
             reg.exp=NA,           
             loc=NA,
             loc.start=NA,
             loc.stop=NA,
             additional=NA,
             names.additional=NA,
             links=NA,
             names.links=NA,
             images=NA,
             names.images=NA,
             band.info = NA,
             returnVl = TRUE,
             saveFile = FALSE,
             saveName="MapObj.RData")




}

\arguments{
  \item{file}{path name of file}
  \item{obj}{object of the class marrayInfo, marrayRaw or marrayNorm}
  \item{adf}{annotated data frame object}
  \item{df}{data frame object}
  \item{spot.ID}{column in file or object that represents spot.IDs i.e. clone, BAC }
  \item{chrom}{column in file or object that represents the chromosome location}
  \item{chrom.levels}{vector indicating how the chrom column is
    represented. i.e chr1, chrom1, 1, etc.}
  \item{base.chrm}{vector of characters to be replaced in chromosome to
  achieve numeric vector}
  \item{reg.exp}{logical indication if base.chrm entry is exact match or
  should be used in perl regular expression mapping}
  \item{loc}{column in file or object representing genomic location. See details for
    more information}
  \item{loc.start}{column in file or object representing starting genomic location. See details for
    more information}
  \item{loc.stop}{column in file or object representing ending genomic location. See details for
    more information}
  \item{locBy}{either within or across to represent if the genomic
  locations are within, with respect to, chromosome or across the entire genomic.}
  \item{file.sep}{seperation character for file}
  \item{additional}{indication of additional columns in file or object that should
    be included in mapping object}
  \item{names.additional}{optional header names for additional columns to be included}
  \item{links}{indication of columns in file or object that represent hyperlinks,
    character or numeric, or a data.frame, with rows equal to the number of
    rows of the file, of hyperlinks that must be in the same order as they
    appear in the file}
  \item{names.links}{optional header names for links}
  \item{images}{indication of columns in file or object that represent images,
    character or numeric, or a data.frame, with rows equal to the number of
    rows of the file, of images that must be in the same order as they
    appear in the file}
  \item{names.images}{optional header names for images}
  \item{band.info}{a band.info object. see makeBandInfo}
  \item{returnVl}{logical indicating if object created should be returned}
  \item{saveFile}{logical indicating if object created should be saved to a file}
  \item{saveName}{if saveVl, path name of file}
  \item{...}{additional arguments to the read.table function for reading
    the file}
}

\details{

  The mapping file must minimally contain columns for name,
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
  they need to be specified in the mapping object. If there are links
  in the file, the column number or head name may be specified as
  links. The user may add their own hyperlinks by assigning the links
  argument to a table of the dimension n by m where n is equivalent to
  the number of items in the original file. The links must be in proper
  format of a valid web address, i.e http://completewebpathhere.com .  
  Similarly, images may be displayed using the images argument. 

  
}

\value{
an object of the class mapobj which contains mapping and band information
}

\author{Lori A. Shepherd}

\seealso{\code{\link{mapobj}}, \code{\link{mapping.info}}, \code{\link{bandinfo}}, \code{\link{iGGV}}}

\examples{

#
# The following code will build the same object as data(mapping.info)
#

library("iGenomicViewer")


# writes out text files to use
writeExFiles()

# loads a band.info object
data(Band.Info)

# make mapping information object
mapping.info = mappingObj(file="HB19Kv2.HG18.txt", spot.ID="Clone", chrom="Chromosome",  chrom.levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),loc.start="start", loc.stop="Stop", file.sep="\t", additional=c("Mapped.by", "Flag"), band.info=band.info, header=TRUE, returnVl=TRUE, links=10, names.links="UCSC")


}

\keyword{methods}
