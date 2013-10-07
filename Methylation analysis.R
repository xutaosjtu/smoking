files = dir("methylation/F3", pattern="*txt.gz", full.names = T)
F3 = subset(F3, subset = !is.na(zz_nr_f3_meth))
rownames(F3) = F3$zz_nr_f3_meth
confounder = c( "rtalter", "rcsex", "rtalkkon", "rtbmi")
## Sensitivity analysis with WBC
WBC.F3 = read.csv("WBC_estimation/Data/KF3_estimated_cell_proportions.csv", sep = ";")
WBC.F3[,1] = substr(WBC.F3[,1], 2, 10)
colnames(WBC.F3)[1] = "zz_nr_f3_meth"
F3 = merge(WBC.F3, F3)
##
data = F3

files = dir("methylation/F4", pattern="*txt.gz", full.names = T)
F4 = subset(F4, subset = !is.na(zz_nr_f4_meth))
rownames(F4) = F4$zz_nr_f4_meth
confounder = c( "utalteru", "ucsex", "utalkkon", "utbmi")
## Sensitivity analysis with WBC
WBC.F4 = read.csv("WBC_estimation/Data/KF4_estimated_cell_proportions.csv", sep = ";")
WBC.F4[,1] = substr(WBC.F4[,1], 2, 10)
colnames(WBC.F4)[1] = "zz_nr_f4_meth"
F4 = merge(WBC.F4, F4)
##
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
                + CD8T + CD4T + NK + Bcell + Mono + Gran # adjust for WBC
				, data = data
          #, family = binomial
               )
		rst.chr = rbind(rst.chr, c(summary(model)$coefficients[2,], summary(model)$coefficients[3,]))
		line=readLines(f,n=1)
	}
	close(f)
	rownames(rst.chr) = CpG
	write.csv(rst.chr, paste("association in chr",chrnum,"_adjust WBC.csv"))
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

## merge the results
files = dir("Results/F3", pattern="*.csv", full.names = T)
files = files[1:22]
association.cpg.f3 = NULL
association.cpg.f3 = foreach(i = 1:length(files), .combine = rbind) %dopar%{
	tmp = read.csv(files[i], row.names =1)
	 return(tmp)
}


files = dir("Results/F4", pattern="*.csv", full.names = T)
files = files[1:22]
association.cpg.f4 = NULL
association.cpg.f4 = foreach(i = 1:length(files), .combine = rbind) %dopar%{
	tmp = read.csv(files[i], row.names =1)
	 return(tmp)
}

cpg.confirmed = intersect(association.cpg.f3, association.cpg.f4)

## meta-analysis
files = dir("Results/F3", pattern = "*.csv")

rst = foreach(i = 1:length(files), .combine = rbind) %dopar%{
	f4.rst = read.csv(paste("Results/F4/",files[i], sep = ""),row.names = 1)
	f3.rst = read.csv(paste("Results/F3/",files[i], sep = ""), row.names = 1)

	f4.rst = f4.rst[, 5:8]
	f3.rst = f3.rst[, 1:4]

	validcpg = intersect(rownames(f4.rst), rownames(f3.rst))

	analysis.meta<-function(x, rst1, rst2,...){
		#tmp = rbind(rst1[x,], rst2[x,])
		rma.test = rma(yi = c(rst1[x,1], rst2[x,1]), sei = c(rst1[x,2], rst2[x,2]),
					method = "FE",
					measure = "GEN")
		summary = c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp)
		return(summary)
	}

	meta.rst = sapply(validcpg, analysis.meta, rst1 = f4.rst, rst2 = f3.rst)
	meta.rst = t(meta.rst)
	rownames(meta.rst) = validcpg
	colnames(meta.rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
	write.csv(meta.rst, file = files[i])
	#return(meta.rst)
}