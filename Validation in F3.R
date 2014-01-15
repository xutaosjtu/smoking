## Mediation of gene expression for the association between methylation and metabolites.
## F4 and F3 used different platform for expression measurement.  Not all the candidates discovered in F4 are measured in F3.
rst = NULL; motifs.tmp = NULL
for(i in 1:nrow(motifs)){
  e = as.character(motifs$expression[i])
  if(e %in% rownames(F3.expression)){
    F3.sub$expression = t(F3.expression[e,as.character(F3.sub$zz_nr_f3_genexp)])
    m = as.character(motifs$metabolites[i])
    F3.sub$metabo = scale(log(F3.sub[ ,m]))
    cpg = as.character(motifs$cpg[i])
    F3.sub$cpg = t(F3.methy[cpg, as.character(F3.sub$zz_nr_f3_meth)])
    test = sobel.lm(pred = F3.sub$cpg, med = F3.sub$expression, out = F3.sub$metabo, covariates= data.frame(as.factor(F3.sub$my.cigreg), F3.sub$rtalteru, as.factor(F3.sub$rcsex), F3.sub$rtbmi, F3.sub$rtalkkon))
    rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[2,], test$'Mod2: Y~X+M'[3,]))
  }
  else{
    rst = rbind(rst, rep(NA,13))
  }
}
rst = data.frame(rst[,1:5], fdr = p.adjust(rst[,5], method = "BH"), bonf = p.adjust(rst[,5], method = "bonf"), rst[,6:13])
colnames(rst)[1:5] = c("Indirect.Effect", "SE", "z.value", "N", "P")
colnames(rst)[8:11] = c("cpg.estimate", "cpg.SE", "cpg.tvalue", "cpg.Pvalue")
colnames(rst)[12:15] = c("expr.estimate", "expr.SE", "expr.tvalue", "expr.Pvalue")
write.csv(cbind(motifs, rst), file = "Mediator analysis from the candidates_F3.csv", row.names = FALSE)


## Mediation of cpg or gene for the association between smoking and metabolites.
F3.sub2 = subset(F3, !is.na(F3$zz_nr_f3_genexp))
for(m in metabo.asso){
  rst = NULL
  F3.sub2$metabo = scale(log(F3.sub2[, m]))
    
  for(e in rownames(F3.expression)){
    if(e %in% rownames(F4.expression)){
      F3.sub2$expression = t(F3.expression[e,as.character(F3.sub2$zz_nr_f3_genexp)])
      test = sobel.lm(pred = as.factor(F3.sub2$my.cigreg), med = F3.sub2$expression, out = F3.sub2$metabo, covariates=data.frame(F3.sub2$rtalteru, F3.sub2$rcsex, F3.sub2$rtbmi, F3.sub2$rtalkkon))
      rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[2,], test$'Mod2: Y~X+M'[3,]))
    }
    else(rst = rbind(rst, rep(NA, 13)))
  }
  
  F3.sub$metabo = scale(log(F3.sub[, m]))
  for(cpg in rownames(F3.methy)){
    F3.sub$cpg = t(F3.methy[cpg, as.character(F3.sub$zz_nr_f3_meth)])
    test = sobel.lm(pred = as.factor(F3.sub$my.cigreg), med = F3.sub$cpg, out = F3.sub$metabo, covariates=data.frame(F3.sub$rtalteru, F3.sub$rcsex, F3.sub$rtbmi, F3.sub$rtalkkon))
    rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[2,], test$'Mod2: Y~X+M'[3,]))
  }
  
  colnames(rst)[1:4] = c("Indirect.Effect", "SE", "z.value", "N")
  colnames(rst)[6:9] = c("S.estimate", "S.SE", "S.tvalue", "S.Pvalue")
  colnames(rst)[10:13] = c("mediator.estimate", "mediator.SE", "mediator.tvalue", "mediator.Pvalue")
  
  rst = data.frame(mediator = c(rownames(F3.expression),rownames(F3.methy)), rst)
  
  write.csv(rst, paste(m, "_expression and methylation mediation_F3.csv", sep = ""),row.names=FALSE)
}

## Mediation of methylation for the association between smoking and gene expression.
rst = NULL
for(e in rownames(F3.expression)){
  if(e %in% rownames(F4.expression)){
    F3.sub$expression = t(F3.expression[e,as.character(F3.sub$zz_nr_f3_genexp)])
    
    for(cpg in rownames(F3.methy)){
      F3.sub$cpg = t(F3.methy[cpg, as.character(F3.sub$zz_nr_f3_meth)])
      test = sobel.lm(pred = as.factor(F3.sub$my.cigreg), med = F3.sub$cpg, out = F3.sub$expression, covariates=data.frame(F3.sub$rtalteru, as.factor(F3.sub$rcsex), F3.sub$rtbmi, F3.sub$rtalkkon))
      rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[3,], test$'Mod2: Y~X+M'[4,]))
    }
  }
  #else (rst = rbind(rst, rep(NA, 13)))
}
colnames(rst)[1:4] = c("Indirect.Effect", "SE", "z.value", "N")
colnames(rst)[6:9] = c("S.estimate", "S.SE", "S.tvalue", "S.Pvalue")
colnames(rst)[10:13] = c("mediator.estimate", "mediator.SE", "mediator.tvalue", "mediator.Pvalue")

expression = rep(intersect(rownames(F3.expression), rownames(F4.expression)),each = 361)
expression = cbind(as.character(annotation[expression,2]), expression)
methylation = rep(rownames(F3.methy), 12)
methylation = cbind(methylation, fullannot[methylation,22])
rst = cbind(expression, methylation, rst)
write.csv(rst, file = "Methylation mediated smoking_expression association_F3.csv", row.names = FALSE)  
