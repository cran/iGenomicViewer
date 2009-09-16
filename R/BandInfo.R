#########################################################
# function to make band.info object [old cytoband object]
# gives locations for chromosome, arm, broad and fine band
###########################################################
#
# function does assume that chrom.levels - all chromosome
#    entries in file share a common prefix if applicable
#    i.e. chr, chrom, ch, etc. 
#

makeBandInfo <- function(file,
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
                         ...
                         ){

  # read file
  cb = read.table(file, sep=file.sep, ...)
  # set up objects
  offset = 0
  band.info = list()

 cbnm = names(cb)
 if(class(chr.dx)=="character") chr.dx = which(cbnm == chr.dx)
 if(class(start.dx)=="character") start.dx = which(cbnm == start.dx)
 if(class(stop.dx)=="character") stop.dx = which(cbnm == stop.dx)
 if(class(band.dx)=="character") band.dx = which(cbnm == band.dx)

  #
  # start with Chrom 
  #
  Chrom = list()
  Chrom$Chrom = chrom.levels
  Chrom$label = c(autosomes, "X", "Y")
  lower= rep(NA, length(chrom.levels))
  upper= rep(NA, length(chrom.levels))
  for(i in 1:length(chrom.levels)){
    idx = which(as.character(cb[,chr.dx]) == chrom.levels[i])
    if(i == 1){
      lower[i] =  (cb[idx[1], start.dx] + 0)
      upper[i] =  (cb[idx[length(idx)], stop.dx] + 0)
      offset[i+1] =  (cb[idx[length(idx)], stop.dx])
    }
    if(i > 1){
      lower[i] =  (cb[idx[1], start.dx] + offset[i])
      upper[i] =  (cb[idx[length(idx)], stop.dx] + offset[i])
#      if(i < length(chrom.levels)) offset[i+1] =  (offset[i]+ upper[i])
      if(i < length(chrom.levels)) offset[i+1] =  (offset[i]+cb[idx[length(idx)], stop.dx])

    }
  }
  Chrom$lower = lower
  Chrom$center = ((as.numeric(lower) + as.numeric(upper))/2)
  Chrom$upper = upper
  Chrom = as.data.frame(Chrom)



  #
  # Arm 
  #
  
  arm.b = paste(gsub(as.character(cb[,chr.dx]), pattern=(strsplit(chrom.levels[1], split="1")[[1]][1]), replacement=""), substring(as.character(cb[,band.dx]), first=1, last=1), sep="")
  narm = length(unique(arm.b))
  uarm = unique(arm.b)

  Arm = list()
  lower = rep(NA, narm)
  upper = rep(NA, narm)
  for(i in 1:length(uarm)){
    idx = which(arm.b == uarm[i])
    offidx = gsub(uarm[i], pattern="p", replacement="")
    offidx = gsub(offidx, pattern="q", replacement="")
    if(offidx == "X") offidx = X.chrom
    if(offidx == "Y") offidx = Y.chrom
    offidx = as.numeric(offidx) 
    lower[i] = cb[idx[1], start.dx] + offset[offidx]
    upper[i] = cb[idx[length(idx)], stop.dx] + offset[offidx]
  }
  center = (as.numeric(lower)+as.numeric(upper))/2
  odx = order(center)
  Arm$Arm = uarm[odx]
  Arm$label = rep(NA, length(uarm))
  Arm$label[grep(Arm$Arm, pattern="p")] = "p"
  Arm$label[grep(Arm$Arm, pattern="q")] = "q"
  Arm$lower = lower[odx]
  Arm$center =center[odx]
  Arm$upper = upper[odx]
  Arm = as.data.frame(Arm)




  #
  # Broad.Band 
  #
  band = as.character(cb[,band.dx])
  for(i in 1:length(band)){
    band[i] = strsplit(band[i], split=".", fixed=TRUE)[[1]][1]
  }
  c.band = gsub(as.character(cb[,chr.dx]), pattern=(strsplit(chrom.levels[1], split="1")[[1]][1]), replacement="")
  b.band = paste(gsub(as.character(cb[,chr.dx]), pattern=(strsplit(chrom.levels[1], split="1")[[1]][1]), replacement=""),band, sep="")
  ubband = unique(b.band)
  nbband = length(ubband)

  Broad.Band = list()
  lower = rep(NA, nbband)
  upper = rep(NA, nbband)
  for(i in 1:length(ubband)){
    idx = which(b.band == ubband[i])
    
    offidx = unique(c.band[idx])
    if(offidx == "X") offidx = X.chrom
    if(offidx == "Y") offidx = Y.chrom
    offidx = as.numeric(offidx) 
    
    lower[i] = cb[idx[1], start.dx] + offset[offidx]
    upper[i] = cb[idx[length(idx)], stop.dx] + offset[offidx]
  }
  center = (as.numeric(lower)+as.numeric(upper))/2
  odx = order(center)
  Broad.Band$Broad.Band = ubband[odx]
  Broad.Band$label = rep(NA, nbband)
  for(j in 1:nbband){
    Broad.Band$label[j] = band[which(b.band == Broad.Band$Broad.Band[j])[1]]
  }
  Broad.Band$lower = lower[odx]
  Broad.Band$center =center[odx]
  Broad.Band$upper = upper[odx]
  Broad.Band = as.data.frame(Broad.Band)
  

  #
  # Fine.Band
  #
  band = as.character(cb[,band.dx])
  c.band = gsub(as.character(cb[,chr.dx]), pattern=(strsplit(chrom.levels[1], split="1")[[1]][1]), replacement="")
  f.band = paste(gsub(as.character(cb[,chr.dx]), pattern=(strsplit(chrom.levels[1], split="1")[[1]][1]), replacement=""),band, sep="")
  ufband = unique(f.band)
  nfband = length(ufband)

  Fine.Band = list()
  lower = rep(NA, nfband)
  upper = rep(NA, nfband)
  for(i in 1:length(ufband)){
    idx = which(f.band == ufband[i])
    
    offidx = unique(c.band[idx])
    if(offidx == "X") offidx = X.chrom
    if(offidx == "Y") offidx = Y.chrom
    offidx = as.numeric(offidx) 
    
    lower[i] = cb[idx[1], start.dx] + offset[offidx]
    upper[i] = cb[idx[length(idx)], stop.dx] + offset[offidx]
  }
  center = (as.numeric(lower)+as.numeric(upper))/2
  odx = order(center)
  Fine.Band$Fine.Band = ufband[odx]
  Fine.Band$label = rep(NA, nfband)
  for(j in 1:nfband){
    Fine.Band$label[j] = band[which(f.band == Fine.Band$Fine.Band[j])[1]]
  }
  Fine.Band$lower = lower[odx]
  Fine.Band$center =center[odx]
  Fine.Band$upper = upper[odx]
  Fine.Band = as.data.frame(Fine.Band)
  Fine.Band = as.data.frame(Fine.Band)
  
  band.info$offset = offset
  band.info$Chrom = Chrom
  band.info$Arm = Arm
  band.info$Broad.Band = Broad.Band
  band.info$Fine.Band = Fine.Band
  
  class(band.info) <- "bandinfo"
  
  if(saveFile) save(band.info, file=saveName, compress=TRUE)
  if(returnVl) return(band.info)
  
}
