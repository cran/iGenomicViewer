#
# writes out data text files for use with examples
#   assumes data has same name as fileList 
#

writeExFiles <- function(direct="./"){

  fileList = c("CancerGenes", "cytoBand", "DiseaseGenes", "DNArepairgenes", "HB19Kv2.HG18")

  for(i in fileList){
    eval.js(paste("data(",i,")", sep=""))
    eval.js(paste("write.table(",i,",sep='\t', row.names=FALSE, file='",direct,i,".txt')", sep=""))
  }


}
