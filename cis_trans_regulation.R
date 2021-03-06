setwd("../")
# load("expression/Expression_f4_adjusted_technical_variables.Rdata")
# F4.expr = adj.expr.f4
# remove(adj.expr.f4)

load("intermediate results/expression candidates.RData")
F4.expr = F4.expression

source("source/data import and preprocess.R")
files = dir("methylation/F4", pattern="*txt.gz", full.names = T)
files = files[1:22]
F4 = subset(F4, !is.na(zz_nr_f4_meth) & !is.na(zz_nr_s4f4_genexp) & expr_in_F4!="")
confounder = c( "utalteru", "ucsex", "utalkkon", "utbmi")
## Sensitivity analysis with WBC
WBC.F4 = read.csv("WBC_estimation/Data/KF4_estimated_cell_proportions.csv", sep = ";")
WBC.F4[,1] = substr(WBC.F4[,1], 2, 10)
colnames(WBC.F4)[1] = "zz_nr_f4_meth"
F4$zz_nr_f4_meth = as.character(F4$zz_nr_f4_meth)
F4 = merge(WBC.F4, F4)
rownames(F4) = F4$zz_nr_f4_meth

data = F4

require(doMC)
registerDoMC(cores = 22)

exprmethyAsso.2 = function(x, data, gene.expr){
  data$expr = gene.expr[x,]
  model = lm(mvalue ~ expr
    #+ as.factor(my.cigreg)
    #          + utalteru + as.factor(ucsex)
     #         + as.factor(my.alkkon) + utbmi #+ as.factor(rtdiabet)
              + CD8T + CD4T + NK + Bcell + Mono + Gran # adjust for WBC
              , data = data
              #, family = binomial
  )
  return(summary(model)$coef[2,4])
}


exprmethyAsso = function(x, F4.expr){
  data$expr = F4.expr[x,]
  model = cor.test(data$mvalue, data$expr, method = "spearman")
  return(model$p.value)
}


# for(i in 1:length(files)){
#   chrnum = unlist(strsplit(files[i], split = "_"))[6]
#   f <- gzfile(files[i], "r")
#   
#   out = file(paste("association in chr",chrnum,"_adjust WBC.csv"),"w")
#     
#   line=readLines(f,n=1)
#   samples = unlist(strsplit(line, split = " "))[-1]
#   samples = substr(samples, 2, 10)
#   
#   indx = samples %in% data$zz_nr_f4_meth
#   samples = samples[indx]
#   
#   data = data[samples, ]
#   F4.expr = F4.expr[,as.character(data$zz_nr_s4f4_genexp)]
#   line=readLines(f,n=1)
#   CpG = NULL;
#   while( length(line) != 0 ) {
#     tmp = unlist(strsplit(line, split = " "))
#     CpG = c(CpG,tmp[1])
#     mvalue = as.numeric(tmp[-1])
#     data$mvalue = mvalue[indx]
#     
#     associations = foreach(i = 1:nrow(F4.expr), .combine = rbind ) %dopar% {
#       model = cor.test(data$mvalue, F4.expr[i,])
#       return(model$p.value)
#     }
#     
#     #writeLines(paste(as.character(associations), collapse = ","), out)
#     
#     line=readLines(f,n=1)
#   }
#   close(f)
#   close(out)
# }

associations = foreach(i = 1:length(files), .combine = rbind ) %dopar% {
  chrnum = unlist(strsplit(files[i], split = "_"))[6]
  f <- gzfile(files[i], "r")
  
  line=readLines(f,n=1)
  samples = unlist(strsplit(line, split = " "))[-1]
  samples = substr(samples, 2, 10)

  indx = samples %in% data$zz_nr_f4_meth
  samples = samples[indx]
  
  data = data[samples, ]
  F4.expr = F4.expr[,as.character(data$zz_nr_s4f4_genexp)]
  line=readLines(f,n=1)
  rst.chr= NULL; CpG=NULL
  while( length(line) != 0 ) {
    tmp = unlist(strsplit(line, split = " "))
    CpG = c(CpG,tmp[1])
    mvalue = as.numeric(tmp[-1])
    data$mvalue = mvalue[indx]
        
    asso = sapply(rownames(F4.expr), exprmethyAsso.2, data=data, gene.expr = F4.expr)
    rst.chr = rbind(rst.chr, asso)
    
    line=readLines(f,n=1)
  }
  close(f)
  rownames(rst.chr) = CpG
  colnames(rst.chr) = rownames(F4.expr)
  write.csv(rst.chr, file = paste("Adj WBC/association in chr",chrnum,"_WBC only.csv"))
  return(rst.chr)
}


qqplot.2 = function(pvalues, ...){
  x = -log10(runif(length(pvalues)))
  y = -log10(pvalues)
  #names(y) = rownames(F3.asso.expr)
  x = sort(x)
  y = sort(y)
  plot(x=x, y =y, pch = 19, col = (y>7)+1,
       xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", ...)
  abline(0,1, lty =2)
}

#files = dir("Results/WBC/methylation expression association/", full.names=T)
files = dir("No adj/", full.names=T)
associations = NULL
for(i in 1:length(files)){
  tmp = read.csv(files[i], row.names=1)
  associations = rbind(associations, tmp)
}

png("Association of cpg sites with gene expression_no adj_EWAS.png", width = 40, height = 60, unit = "in", res = 180)
par(mfrow = c(6, 4))
for(i in 1:ncol(associations)){
  qqplot.2(associations[,i], main = rownames(F4.expr)[i] )
}
dev.off()



