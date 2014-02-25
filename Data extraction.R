require(doMC)
registerDoMC(cores = 8)
files = dir("methylation/F3", pattern="*txt.gz", full.names = T)
files = dir("methylation/F4", pattern="*txt.gz", full.names = T)
files = files[1:22]
#################################################################
## Identification of candidates
#################################################################
#anno.expr = read.csv("expression/Annotation HumanHT-12v3 final.csv")

source("Genomic_infor.R")

exprprobe.1M = sapply(cpg.candidates, exprin1M)


cpg.1M = sapply(rownames(F4.expression), cpgin1M)
cpg.68K = sapply(rownames(F4.expression), cpgin, region = 68000)
#hm450.hg18 <-getPlatform(platform="HM450", genome = 'hg18')
#show(hm450.hg18)

load("cpg and expression around 1Mb.RData")
#################################################################
## run the code in mehtylation analysis to obtain cpg.confirmed
#################################################################
methy.data = foreach(i = 1:length(files), .combine = rbind ) %dopar% {
	f <- gzfile(files[i], "r")
	line=readLines(f,n=1)
	samples = unlist(strsplit(line, split = " "))[-1]
	samples = substr(samples, 2, 10)
	line=readLines(f,n=1)
	data = NULL
	while( length(line) != 0 ) {
		 tmp = unlist(strsplit(line, split = " "))
		 if(tmp[1] %in% cpg_candidate){
			data = rbind(data, tmp)
		 }
		 line=readLines(f,n=1)
	}
	#colnames(data) = samples
	close(f)
	return(data)
}
write.csv(methy.data, file = "methylation data of candidates_F4.csv")

rst = NULL
for(i in files){
	f<-gzfile(i, "r")
	line = readLines(f, n=1)
	samples = unlist(strsplit(line, split = " "))[-1]
	samples = substr(samples, 2, 10)
	close(f)
	rst = rbind(rst, samples)
}

### extract associations with smoking
#files = dir("Results/F3_adjusted for WBC", pattern="*.csv", full.names = T)
files = dir("Results/F4_adjusted for WBC", pattern="*.csv", full.names = T)
files = files[1:22]
#cpgs = unlist(cpgs_of_gene)
association.F4 = foreach(i = 1:length(files), .combine = rbind ) %dopar% {
  tmp = read.csv(files[i], row.names =1)
#  tmp.cpgs = intersect(rownames(tmp), cpgs)
  tmp.cpgs = which(tmp[,8]<1e-7)
  return(tmp[tmp.cpgs, ])
}

## non-paralellized version
association.F4 = NULL
for(i in 1:length(files)){
  tmp = read.csv(files[i], row.names =1)
  #  tmp.cpgs = intersect(rownames(tmp), cpgs)
  tmp.cpgs = which(tmp[,8]<1e-7)
  association.F4 = rbind(association.F4, tmp[tmp.cpgs, ])
}


dim(association.F4)
write.csv(association.F4, file = "Smoking associated methylation sites in F4.csv")
cpgs=rownames(association.F4)
#association.F4 = association

files = dir("Results/F3_adjusted for WBC/", pattern="*.csv", full.names = T)
files = files[1:22]
association.F3 = foreach(i = 1:length(files), .combine = rbind ) %dopar% {
    tmp = read.csv(files[i], row.names =1)
    tmp.cpgs = which(tmp[,4]<0.05/nrow(association.F4))
    return(tmp[tmp.cpgs, ])
}

## non-paralellized version
association.F3 = NULL
for(i in 1:length(files)){
  tmp = read.csv(files[i], row.names =1)
  tmp.cpgs = intersect(rownames(tmp), cpgs)
  #tmp.cpgs = which(tmp[,4]<1e-5)
  association.F3 = rbind(association.F3, tmp[tmp.cpgs, ])
}

dim(association.F3)

rownames(fullannot) = fullannot[,1]
association.F3 = cbind(association.F3, fullannot[rownames(association.F3), ])

write.csv(association.F3, file = "Smoking associated methylation sites in F4 replication in F3.csv")
cpg_candidate = intersect(rownames(association.F3), rownames(association.F4))


#################################################################
## Candidates of differentially expressed genes: CLDND1, LRRN3
#################################################################
probes = c("ILMN_1773650","ILMN_2048591" ##LRRN3
           ,"ILMN_1710326", "ILMN_2352563"  ##CLDND1
)

"ILMN_1710326","ILMN_1773650"

F3.expression = read.csv("expression/KORA_F3_379_samples_normalized_expression_data_zzupdated.csv", row.names = 1)
F4.expression = load("Expression_s4_adjusted_technical_variables")

F3.expression = F3.expression[probes,]
F4.expression = F4.expression[probes,]

#################################################################
## Candidate metabolites according to publication in BMC medicine 
#################################################################




