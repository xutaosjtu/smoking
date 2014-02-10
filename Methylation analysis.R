files = dir("methylation/F3", pattern="*txt.gz", full.names = T)
F3 = subset(F3, subset = !is.na(zz_nr_f3_meth))
rownames(F3) = F3$zz_nr_f3_meth
confounder = c( "rtalteru", "rcsex", "rtalkkon", "rtbmi")
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
registerDoMC(cores = 16)

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
                + rtalteru + as.factor(rcsex)
                + as.factor(my.alkkon) + rtbmi + #as.factor(rtdiabet)
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

## give row names to the result
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



## heatmap plot
F4.methy = F4.methy[, as.character(F4.sub$zz_nr_f4_meth)]
rst = apply(F4.methy, 1, function(x)  tapply(x, INDEX = F4.sub$my.cigreg, mean, na.rm = T))
pdf("methylation in F4.pdf")
heatmap.2(scale(t(rst)), trace = "none", col = greenred, scale = "none")
dev.off()


outlier = function(x){
  x[x > (mean(x, na.rm = T) + 4*sd(x, na.rm = T))|x < (mean(x, na.rm = T) - 4*sd(x, na.rm = T))] = NA
  return(x)
}

F3.methy = F3.methy[, as.character(F3.sub$zz_nr_f3_meth)]
rst = apply(F3.methy, 1, function(x)  tapply(x, INDEX = F3.sub$my.cigreg, mean, na.rm = T))
heatmap.2(rst, Colv=F, dendrogram="row", trace = "none", scale = "none", col = greenred )


F4.methy = read.csv("Results/WBC/Smoking associated methylation/Smoking associated methylation sites in F4.csv", row.names = 1)
F3.methy = read.csv("Results/WBC/Smoking associated methylation/Smoking associated methylation sites in F4 replication in F3.csv", row.names = 1)
heatmap.2(cbind(F4.methy$Estimate.1, F3.methy$Estimate.1), dendrogram="none", Rowv=NA, scale = "none")

F3.methy = F3.methy[order(F3.methy$CHR, F3.methy$COORDINATE_36),]
F4.methy = F4.methy[rownames(F3.methy),]

color = rainbow(22)
plot((F4.methy$Estimate.1), type="h", col = F3.methy$CHR)
color = rainbow(22,alpha= 0.5)
points((F3.methy$Estimate), type = "p", col = F3.methy$CHR)

range(F3.methy$Pr...t...1)
