#require("FDb.InfiniumMethylation.hg18")
require("illuminaHumanv3.db")
require("IlluminaHumanMethylation450k.db")

## Expression probe --> Gene Symbol
x <- illuminaHumanv3SYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
exprprobe2symbol <- as.list(x[mapped_probes])

## Expression probe --> Chromsome
x <- illuminaHumanv3CHR
# Get the probe identifiers that are mapped to a chromosome
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
exprprobe2chr = unlist(xx)

## Expression probe --> Start Position
x <- illuminaHumanv3CHRLOC
# Get the probe identifiers that are mapped to chromosome locations
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
exprprobe2start = sapply(xx, function(x){
  if(x<0) x = org.Hs.egCHRLENGTHS[names(x)]+x
  return(x)
})
exprprobe2start = unlist(exprprobe2start)
names(exprprobe2start) = substr(names(exprprobe2start), 1, 12)

## Expression probe --> End position
x <- illuminaHumanv3CHRLOCEND
# Get the probe identifiers that are mapped to chromosome locations
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
exprprobe2end = sapply(xx, function(x){
  if(x<0) x = org.Hs.egCHRLENGTHS[names(x)]+x
  return(x)
})
exprprobe2end = unlist(exprprobe2end)
names(exprprobe2end) = substr(names(exprprobe2end), 1, 12)


## CpG site --> Gene symbol
x <- IlluminaHumanMethylation450kSYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
cpg2symbol <- as.list(x[mapped_probes])

## CpG site --> Chromsome
x <- IlluminaHumanMethylation450kCHR37
# Get the probe identifiers that are mapped to a chromosome
mapped_probes <- mappedkeys(x)
# Convert to a list
cpg2chr <- as.list(x[mapped_probes])
cpg2chr = unlist(cpg2chr)

## CpG site --> Position
cpg2position = as.list(IlluminaHumanMethylation450kCPG37)
cpg2position = unlist(cpg2position)



## Probes of gene expression within the 1M region of a cpg site
exprin1M = function(cpg){
  chr = cpg2chr[cpg]
  position = cpg2position[cpg]
  
  exprprobes = names(exprprobe2chr)[which(exprprobe2chr==as.character(chr))]
  exprprobes.start = exprprobe2start[exprprobes]
  probes = names(exprprobes.start)[which(exprprobes.start<position+1000000 & exprprobes.start>position-1000000)]
  return(probes)
}

exprin = function(cpg, region = 1000000){
  chr = cpg2chr[cpg]
  position = cpg2position[cpg]
  exprprobes = names(exprprobe2chr)[which(exprprobe2chr==as.character(chr))]
  exprprobes.start = exprprobe2start[exprprobes]
  probes = names(exprprobes.start)[which(exprprobes.start<position+region & exprprobes.start>position-region)]
  return(probes)
}

## CpG sites within the 1M region of a expression probe
cpgin1M = function(expr, region = 1000000){
  chr = exprprobe2chr[expr]
  start = exprprobe2start[expr]
  end = exprprobe2start[expr]
  
  cpgs = names(cpg2chr)[which(cpg2chr==as.character(chr))]
  cpgs.pos = cpg2position[cpgs]
  probes = names(cpgs.pos)[which(cpgs.pos>start-region & cpgs.pos<start+region)]
  return(probes)
}


cpgin = function(expr, region = 1000000){
  chr = exprprobe2chr[expr]
  start = exprprobe2start[expr]
  end = exprprobe2start[expr]
  
  cpgs = names(cpg2chr)[which(cpg2chr==as.character(chr))]
  cpgs.pos = cpg2position[cpgs]
  probes = names(cpgs.pos)[which(cpgs.pos>start-region & cpgs.pos<start+region)]
  return(probes)
}
