setwd("../Dropbox/Smoking Systems biology/")
F4.mediation = read.csv("Results/WBC/Mediator analysis from the candidates_F4.csv")
F3.mediation = read.csv("Results/WBC/Mediator analysis from the candidates_F3.csv")

F4.mediation = subset(F4.mediation, F4.mediation$expression %in% c("ILMN_1710326","ILMN_1773650"))

F4.significant = subset(F4.mediation, F4.mediation$P<0.05 & F4.mediation$expr.Pvalue<0.05)
F3.replicate = merge(F4.significant, F3.mediation, by = c("cpg", "metabolites","expression"))
F3.replicate = F3.replicate[, c(1:3, 19:33)]

calculateinterval = function(estimation, SE){
  tmp = paste('(',
  round(estimation-1.96*SE,2), 
  ',', 
  round(estimation+1.96*SE,2),
  ')',
              sep= "")
}

F3.replicate$Indirect.Effect.interval = paste(
  '(',
  round(F3.replicate$Indirect.Effect.y-1.96*F3.replicate$SE.y,2), 
  ',', 
  round(F3.replicate$Indirect.Effect.y+1.96*F3.replicate$SE.y,2),
  ')')
F3.replicate$cpg.estimate.interval = paste(
  '(',
  round(F3.replicate$cpg.estimate.y-1.96*F3.replicate$cpg.SE.y,2), 
  ',',
  round(F3.replicate$cpg.estimate.y+1.96*F3.replicate$cpg.SE.y,2),
  ')')
F3.replicate$expr.estimate.interval = paste(
  '(',
  round(F3.replicate$expr.estimate.y-1.96*F3.replicate$expr.SE.y,2), 
  ',', 
  round(F3.replicate$expr.estimate.y+1.96*F3.replicate$expr.SE.y,2),
  ')')
write.csv(F3.replicate, file = "F3 replicate of meidator analysis.csv", row.names = FALSE)


F4.mediation = read.csv("Results/WBC/Methylation mediated smoking_expression association_F4.csv")
F3.mediation = read.csv("Results/WBC/Methylation mediated smoking_expression association_F3.csv")

F4.mediation = read.csv("Results/WBC/Expression mediated smoking_methylation association_F4.csv")
F3.mediation = read.csv("Results/WBC/Expression mediated smoking_methylation association_F3.csv")

F4.mediation = subset(F4.mediation, F4.mediation$expression %in% c("ILMN_1710326","ILMN_1773650"))

F4.significant = subset(F4.mediation, F4.mediation$p<0.05/nrow(F4.mediation) & F4.mediation$mediator.Pvalue<0.05/nrow(F4.mediation))
F4.significant$Indirect.Effect.interval = calculateinterval(F4.significant$Indirect.Effect, F4.significant$SE)
F4.significant$S.interval = calculateinterval(F4.significant$S.estimate, F4.significant$S.SE)
F4.significant$mediator.interval = calculateinterval(F4.significant$mediator.estimate, F4.significant$mediator.SE)
write.csv(F4.significant, file = "expression mediated smoking methylation association_F4.csv", row.names = F)

F3.replicate = merge(F4.significant[,c("methylation","expression")], F3.mediation, by = c("methylation","expression"))

F3.replicate$Indirect.Effect.interval = calculateinterval(F3.replicate$Indirect.Effect, F3.replicate$SE)
F3.replicate$S.interval = calculateinterval(F3.replicate$S.estimate, F3.replicate$S.SE)
F3.replicate$mediator.interval = calculateinterval(F3.replicate$mediator.estimate, F3.replicate$mediator.SE)
write.csv(F3.replicate, file = "expression mediated smoking methylation association_F3.csv", row.names = F)

#F3.replicate = F3.replicate[, c(1:4, 18:32)]


tmp = merge(F4.significant[,c("methylation", "expression")], F4.mediation)
tmp2 = subset(tmp, tmp$p<0.05/8303 & tmp$mediator.Pvalue<0.05/8303)

###
files = dir("Results/WBC/Mediation analysis_smoking_metabolites_F4/", full.names=T)

rst = list();names = NULL
for(i in 1:length(files)){
  F4.mediation = read.csv(files[i], row.names = 1)
  gene_proxy = grep(pattern = "IL", rownames(F4.mediation))
  F4.mediation = rbind(F4.mediation[c("ILMN_1710326","ILMN_1773650"),], F4.mediation[-gene_proxy,])
  names = c(names, unlist(strsplit(unlist(strsplit(files[i], split = "/"))[4], split = "_"))[1])
  rst[[i]] = F4.mediation[which(F4.mediation$p<0.05 & F4.mediation$mediator.Pvalue<0.05),]
}
names(rst) = names


## Circos plot
require(OmicCircos)
require(IlluminaHumanMethylation450k.db)
require(illuminaHumanv3.db)

## methylation positions
cpgs = as.character(unique(F4.mediation$methylation))
xx = as.list(IlluminaHumanMethylation450kCHR37)
chr = unlist(xx[cpgs])
xx = as.list(IlluminaHumanMethylation450kCPG37)
start = unlist(xx[cpgs])
end = start + 1
desc = NA
methylation= data.frame(chr, start, end, name = cpgs, desc)
methylation$chr = paste("chr", methylation$chr, sep = "")
#methylation$chr = as.factor(methylation$chr)
#methylation$start = as.character(methylation$start)
#methylation$end = as.character(methylation$end)

## expression position 
exprs = as.character(unique(F4.mediation$expression))
xx = as.list(illuminaHumanv3CHR[mappedkeys(illuminaHumanv3CHR)])
chr = unlist(xx[exprs])
chr = paste("chr", chr, sep = "")
xx = as.list(illuminaHumanv3CHRLOC[mappedkeys(illuminaHumanv3CHRLOC)])
start = abs(unlist(xx[exprs]))
start = start[unique(names(start))]
xx = as.list(illuminaHumanv3CHRLOCEND[mappedkeys(illuminaHumanv3CHRLOCEND)])
end = abs(unlist(xx[exprs]))
end = end[unique(names(end))]
desc = NA
expression = data.frame(chr, start, end, name = exprs, desc)

## create mapping file
mapping = rbind(expression, methylation)
seg.num = 22
seg.name = paste("chr", 1:22, sep = "")
db = segAnglePo(mapping, seg = seg.name)

## read in methylation value
F4.methy = read.csv("Results/WBC/Smoking associated methylation/Smoking associated methylation sites in F4.csv", row.names = 1)
F4.methy = F4.methy[cpgs,]
methylation.v = cbind(methylation[,c(1,3)], F4.methy[, 5:8])
names(methylation.v)[1:2] = c("seg.name", "seg.po")

## read in expression value
F4.exprs = read.csv("Results/WBC/Smoking associated expression/F4 smoking associated genes_lm.csv", row.names = 1)
F4.exprs = F4.exprs[exprs,]
expression.v = cbind(expression[,c(1,2)], F4.exprs[, 5:8])
names(expression.v)[1:2] = c("seg.name", "seg.po")

## read in association infor
tmp1 = expression[as.character(F4.significant$expression), 1:3]
tmp2 = methylation[as.character(F4.significant$methylation), 1:3]
link = cbind(tmp1, tmp2)

## read in expression value
colors   <- rainbow(seg.num, alpha=0.5);
pdf("Circos plot2.pdf", height=10, width = 10)
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos(R=180, mapping=methylation.v, type="chr", cir=db, col=colors, print.chr.lab=TRUE, W=10, scale=F);
circos(R=160, cir=db, W=20, mapping=methylation.v, col.v=3, type="s",  B=TRUE, col=colors[2],  scale=T);
circos(R=140, cir=db, W=20, mapping=expression.v, col.v=3, type="b",  B=TRUE, col="blue",  scale=T);
circos(R=120, cir=db, W=40, mapping=link, type="link.pg", col=colors, lwd=2);
dev.off()



F4.methy = read.csv("Results/WBC/Smoking associated methylation/Smoking associated methylation sites in F4.csv", row.names = 1)
F3.methy = read.csv("Results/WBC/Smoking associated methylation/Smoking associated methylation sites in F4 replication in F3.csv", row.names = 1)

start = as.numeric(methylation$start) + as.numeric(db[as.numeric(substr(methylation$chr,4, 5)),"seg.sum.start"])
plot(F4.methy$Estimate.1~ start, cex = 0.5, col = as.numeric(substr(methylation$chr,4, 5)))
start = expression$start + as.numeric(db[as.numeric(substr(expression$chr,4, 5)),"seg.sum.start"])
points(F4.exprs$Estimate.1~ start, type = "h")
abline(h = 0, lty = 2)
text(x = start, y = F4.exprs$Estimate.1, labels = rownames(F4.exprs), cex = 0.5)
