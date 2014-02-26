require(gplots)
F4.expr = read.csv("Results/WBC/Smoking associated expression/F4 smoking associated genes_lm.csv", row.names = 1)
#F3.expr = read.csv("Results/WBC/Smoking associated expression/F3 smoking associated genes_lm.csv")
source("Genomic_infor.R")
load("cpg and expression around 1Mb.RData")

extract = function(cpgs){
  alpha = 0.05/length(cpgs)
  require(doMC)
  registerDoMC(cores = 8)
  files = dir("Results/WBC/Smoking associated methylation/F4_adjusted for WBC/", full.names = T)
  association.F4 = foreach(i = 1:length(files), .combine = rbind ) %dopar% {
    tmp = read.csv(files[i], row.names = 1)
    tmp = tmp[rownames(tmp) %in% cpgs,]
    return(tmp)
  }
  return(association.F4)
}

exprSig = function(expr){
  alpha = 0.05/length(unique(unlist(expr)))
  expr.sig = sapply(expr, 
                    function(x){
                      x = x[F4.expr[x, 8]<alpha]
                      return(x)
                    }
  )
  hits = sapply(expr.sig, length)
  expr.sig = expr.sig[(hits!=0)]
  return(expr.sig)
}

cpgSig = function(cpg){
  alpha = 0.05/length(unique(unlist(cpg)))
  cpg.sig = sapply(cpg, 
                   function(x){
                     x = x[F4.methy[x, 8]<alpha]
                     return(x)
                   }
  )
  cpg.sig = sapply(cpg.sig, function(x) return(x[which(!is.na(x))]))
  hits = sapply(cpg.sig, length)
  cpg.sig = cpg.sig[(hits!=0)]
  return(cpg.sig)
}

#exprprobe2chr[exprprobe.1M[[1]]]
cpg.candidates = names(exprprobe.1M)
## Analyze the gene expression around 1M of cpg sites
cpgs.select = cpg.candidates
expr.select = NULL

while(length(cpg.candidates)!=0){
  ## Searching differentially expressed genes around differentially mehtylated sites
  ##1. Extract genes 1M around the differntially methylated cpg sites
  expr.1M=sapply(cpg.candidates, exprin1M)
  
  ##2. Find the significant differential expression
  expr.sig = exprSig(expr.1M)
  expr.candidates = unique(unlist(expr.sig))
  
  expr.select = c(expr.select, expr.candidates)
  
  ##3. Extract the cpg sites 1M around the differentially expressed genes
  cpg.1M = sapply(expr.candidates, cpgin1M)
  cpgs.all = unique(unlist(cpg.1M))
  F4.methy= extract(cpgs.all)
  
  ##4. Find the significant differential methylation
  cpg.sig = cpgSig(cpg.1M)
  
  cpg.candidates = setdiff(unlist(cpg.sig), cpgs.select)
  print(length(cpg.candidates))
  cpgs.select = c(cpgs.select, cpg.candidates)
}

## Analyze the methylation sites expression around 1M of gene
load("cpg and expression around 1Mb.RData")
expr.candidates = names(cpg.1M)
expr.select = expr.candidates
cpgs.select = NULL

while(length(expr.candidates)!=0){
  ## Searching differentially expressed genes around differentially mehtylated sites
    
  ##1. Extract the cpg sites 1M around the differentially expressed genes
  cpg.1M = sapply(expr.candidates, cpgin1M)
  cpgs.all = unique(unlist(cpg.1M))
  F4.methy= extract(cpgs.all)
  
  ##2. Find the significant differential methylation
  cpg.sig = cpgSig(cpg.1M)
  cpg.candidates = unlist(unlist(cpg.sig))
  
  cpgs.select = c(cpg.candidates, cpgs.select)
  
  ##3. Extract genes 1M around the differntially methylated cpg sites
  expr.1M=sapply(cpg.candidates, exprin1M)
  
  ##4. Find the significant differential expression
  expr.sig = exprSig(expr.1M)
  expr.candidates = unique(unlist(expr.sig))
   
  expr.candidates = setdiff(expr.candidates, expr.select)
  print(length(expr.candidates))
  expr.select = c(expr.select,expr.candidates)
}

########
pdf("CpG sites with significant cis effects on expression.pdf")
for(i in 1:length(cpg.cis)){
  genes = cpg.cis[[i]]
  plotCI(x = exprprobe2start[genes], 
         y = F4.expr[genes, 5], 
         uiw = 1.96* F4.expr[genes, 6], 
         col = (F4.expr[genes, 8]<9.049774e-06)+1, 
         label = unlist(exprprobe2symbol[genes]), 
         xlab = paste("chr", exprprobe2chr[genes[1]]),
         main = names(cpg.cis[i])
  )
  abline(h = 0, lty = 2)
  abline(v = cpg2position[names(cpg.cis[i])], col = "red")
}
dev.off()
       

## extract the cpg sites around the smoking related genes
cpgs.all = unique(unlist(cpg.1M))
cpgs.all[1]

require(doMC)
registerDoMC(cores = 8)
files = dir("Results/WBC/Smoking associated methylation/F4_adjusted for WBC/", full.names = T)
association.F4 = foreach(i = 1:length(files), .combine = rbind ) %dopar% {
  tmp = read.csv(files[i], row.names = 1)
  tmp = tmp[rownames(tmp) %in% cpgs.all,]
  return(tmp)
}

F4.methy = read.csv("F4 methylation around 1M of cpgs.csv", row.names = 1)

expr.cis = sapply(cpg.1M, 
                 function(x){
                   x = x[F4.methy[x, 8]<8.487523e-06]
                   return(x)
                 }
)
hits = sapply(expr.cis, length)
expr.cis = cpg.1M[(hits!=0)]

pdf("CpGs 1M around smoking related genes_2.pdf")
for(i in 1:length(expr.cis)){
  tmp = expr.cis[[i]]
  plotCI(x = cpg2position[tmp], 
         y = F4.methy[tmp, 5], 
         uiw = 1.96* F4.methy[tmp, 6], 
         col = (F4.methy[tmp, 8]<8.487523e-06)+1, 
         #label = unlist(cpg2symbol[tmp]), 
         xlab = paste("chr", cpg2chr[tmp[1]]),
         main = names(expr.cis[i])
  )
  abline(h = 0, lty = 2)
  segments(x0 = exprprobe2start[names(expr.cis[i])],  
          y0 = 0, 
          x1 = exprprobe2end[names(expr.cis[i])],
          y1 = 0,
          col = "red", lwd = 10)
}
dev.off()

pdf("CpGs 1M around smoking related genes_2.pdf")
for(i in 1:length(expr.cis)){
  tmp = expr.cis[[i]]
  y = -log10(F4.methy[tmp, 8])
  subset = which(F4.methy[tmp, 5]<0)
  y[subset] = 0-y[subset]
  plot(x = cpg2position[tmp], 
       y = y, 
       col = (F4.methy[tmp, 8]<8.487523e-06)+1, 
       pch = 19,
       #label = unlist(cpg2symbol[tmp]), 
       ylab = "-log10(P)",
       xlab = paste("chr", cpg2chr[tmp[1]]),
       main = names(expr.cis[i]),
       yaxt = 'n'
  )
  y.range = round(range(y, na.rm = T),0)
  pos = seq(from = y.range[1], to = y.range[2], by = 1)
  axis(2, at = pos, labels = abs(pos))
  abline(h = 0, lty = 2)
  
  expr.asso=F4.expr[names(expr.cis)[i],5]
  segments(x0 = exprprobe2start[names(expr.cis[i])],  
           y0 = 0, 
           x1 = exprprobe2end[names(expr.cis[i])],
           y1 = 0,
           col = c("lightblue","pink")[(expr.asso>0)+1], lwd = 15)
}
dev.off()


