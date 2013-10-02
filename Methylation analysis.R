files = dir("methylation/F3", pattern="*txt.gz", full.names = T)
F3 = subset(F3, subset = !is.na(zz_nr_f3_meth))
rownames(F3) = F3$zz_nr_f3_meth
confounder = c( "rtalter", "rcsex", "rtalkkon", "rtbmi")
data = F3

files = dir("methylation/F4", pattern="*txt.gz", full.names = T)
F4 = subset(F4, subset = !is.na(zz_nr_f4_meth))
rownames(F4) = F4$zz_nr_f4_meth
confounder = c( "utalteru", "ucsex", "utalkkon", "utbmi")
data = F4

for(i in files){
	chrnum = unlist(strsplit(i, split = "_"))[6]
	f <- gzfile(i, "r")
	line=readLines(f,n=1)
	samples = unlist(strsplit(line, split = " "))[-1]
	samples = substr(samples, 2, 10)
	data = data[samples, ]
	line=readLines(f,n=1)
	rst.chr= NULL; CpG=NULL
	while( length(line) != 0 ) {
		 tmp = unlist(strsplit(line, split = " "))
		 CpG = c(CpG,tmp[1])
		 mvalue = as.numeric(tmp[-1])
		 data$mvalue = mvalue
		 
		 model = glm(mvalue ~ as.factor(my.cigreg)
                + rtalter + as.factor(rcsex)
                + as.factor(my.alkkon) + rtbmi #+ as.factor(rtdiabet)
                , data = data
          #, family = binomial
               )
		rst.chr = rbind(rst.chr, c(summary(model)$coefficients[2,], summary(model)$coefficients[3,]))
		line=readLines(f,n=1)
	}
	close(f)
	rownames(rst.chr) = CpG
	write.csv(rst.chr, paste("association in chr",chrnum,".csv"))
}


require(doMC)
registerDoMC(cores = 8)

associations = foreach(i = 1:length(files), .combine = rbind ) %dopar% {
	chrnum = unlist(strsplit(files[i], split = "_"))[6]
	f <- gzfile(files[i], "r")
	line=readLines(f,n=1)
	samples = unlist(strsplit(line, split = " "))[-1]
	samples = substr(samples, 2, 10)
	data = data[samples, ]
	line=readLines(f,n=1)
	rst.chr= NULL; CpG=NULL
	while( length(line) != 0 ) {
		 tmp = unlist(strsplit(line, split = " "))
		 CpG = c(CpG,tmp[1])
		 mvalue = as.numeric(tmp[-1])
		 data$mvalue = mvalue
		 
		 model = glm(mvalue ~ as.factor(my.cigreg)
                + utalteru + as.factor(ucsex)
                + as.factor(my.alkkon) + utbmi #+ as.factor(rtdiabet)
                , data = data
          #, family = binomial
               )
		rst.chr = rbind(rst.chr, c(summary(model)$coefficients[2,], summary(model)$coefficients[3,]))
		line=readLines(f,n=1)
	}
	close(f)
	rownames(rst.chr) = CpG
	write.csv(rst.chr, paste("association in chr",chrnum,".csv"))
	return(rst.chr)
}

##
for(i in files){
	chrnum = unlist(strsplit(i, split = "_"))[6]
	rst.chr = read.csv(paste("association in chr",chrnum,".csv"), row.names = 1)
	f <- gzfile(i, "r")
	line=readLines(f,n=1)
	line=readLines(f,n=1)
	CpG = NULL
	while( length(line) != 0 ) {
		 tmp = unlist(strsplit(line, split = " "))
		 CpG = c(CpG, tmp[1])
		line=readLines(f,n=1)
	}
	close(f)
	row.names(rst.chr) = CpG
	write.csv(rst.chr, paste("association in chr",chrnum,".csv"))
}

##