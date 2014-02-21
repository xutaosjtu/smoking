require(gplots)
F4.expr = read.csv("Results/WBC/Smoking associated expression/F4 smoking associated genes_lm.csv", row.names = 1)
#F3.expr = read.csv("Results/WBC/Smoking associated expression/F3 smoking associated genes_lm.csv")

exprprobe2chr[exprprobe.1M[[1]]]

cpg.cis = sapply(exprprobe.1M, 
  function(x){
    x = x[F4.expr[x, 8]<8.4488e-06]
    return(x)
  }
)
hits = sapply(cpg.cis, length)
cpg.cis = exprprobe.1M[(hits!=0)]

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

expr.cis = sapply(cpg.68K, 
                 function(x){
                   x = x[F4.methy[x, 8]<8.487523e-06]
                   return(x)
                 }
)
hits = sapply(expr.cis, length)
expr.cis = cpg.68K[(hits!=0)]

pdf("CpGs 68K around smoking related genes.pdf")
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

pdf("CpGs 68K around smoking related genes.pdf")
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


