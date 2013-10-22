files = dir("methylation/F3", pattern="*txt.gz", full.names = T)
files = dir("methylation/F4", pattern="*txt.gz", full.names = T)
files = files[1:22]
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
		 #if(tmp[1] %in% cpg.confirmed){
     if(sum(unlist(lapply(cpgs_of_gene, function(x) tmp[1] %in% x )))>=1){
			data = rbind(data, tmp)
		 }
		 line=readLines(f,n=1)
	}
	#colnames(data) = samples
	close(f)
	return(data)
}

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
files = dir("Results/F3", pattern="*.csv", full.names = T)
files = dir("Results/F4", pattern="*.csv", full.names = T)
files = files[1:22]
cpgs = unlist(cpgs_of_gene)
association = foreach(i = 1:length(files), .combine = rbind ) %dopar% {
  tmp = read.csv(files[i], row.names =1)
  tmp.cpgs = intersect(rownames(tmp), cpgs)
  return(tmp[tmp.cpgs, ])
}

#################################################################
## Candidates of differentially expressed genes: CLDND1, LRRN3
#################################################################



#################################################################
## Candidate metabolites according to publication in BMC medicine 
#################################################################
probes = c("ILMN_1773650","ILMN_2048591" ##LRRN3
		,"ILMN_1710326", "ILMN_2352563"  ##CLDND1
		)

F3.expression = read.csv("expression/KORA_F3_379_samples_normalized_expression_data_zzupdated.csv", row.names = 1)
F4.expression = load("Expression_s4_adjusted_technical_variables")

F3.expression = F3.expression[probes,]
F4.expression = F4.expression[probes,]



