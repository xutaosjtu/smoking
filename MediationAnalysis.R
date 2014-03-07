setwd("H:/Smoking and methylation")

## import metabolomics data
source("source/data import and preprocess.R")

## import gene expression data
load("intermediate results/expression candidates.RData")
colnames(F3.expression) = gsub("X", "", colnames(F3.expression))
F3.expression = as.matrix(F3.expression)
F3.expression = F3.expression[c("ILMN_1773650","ILMN_1710326"),]

## import methylation data
F4.methy = read.csv("intermediate results/Methylation_WBC adjust/methylation data of candidates_F4.csv", row.names = 1)
F3.methy = read.csv("intermediate results/Methylation_WBC adjust/methylation data of candidates_F3.csv", row.names = 1)
colnames(F4.methy) = substr(colnames(F4.methy), 2, 10)
colnames(F3.methy) = substr(colnames(F3.methy), 2, 10)
load("metabolites/smoking associated metabolites.RData")

## Subset with all three types of data
F4.sub = subset(F4, expr_in_F4!="" & !is.na(zz_nr_f4_meth) )#& !is.na(zz_nr_f4_bio)
#F4.sub = subset(F4, !is.na(zz_nr_f4_meth) & !is.na(zz_nr_f4_bio))

colnames(F3)[27:189] = gsub("_PTC", "", colnames(F3)[27:189])
colnames(F3)[27:189] = gsub("_", ".", colnames(F3)[27:189])
colnames(F3)[86] = "SM..OH..C22.2"
F3.sub = subset(F3, !is.na(zz_nr_f3_meth) & !is.na(zz_nr_f3_genexp))

## import mediation methods: sobel's test
source("source/SobelTest.R")

## import annotation data
load("methylation/fullannotInd.rda")
annotation = read.csv2("expression/Annotation HumanHT-12v3 final.csv")

require(mediation)


##############
##    F4    ##
##          ##
##############
## Mediation of methylation for the association between smoking and gene expression.
rst = NULL
for(e in rownames(F4.expression)){
  F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
  F4.sub$expression = scale(F4.sub$expression)
  
  for(cpg in rownames(F4.methy)){
    F4.sub$cpg = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
    F4.sub$cpg = scale(F4.sub$cpg)
    test = sobel.lm2(pred = F4.sub$my.cigreg, med = F4.sub$cpg, out = F4.sub$expression, covariates=data.frame(F4.sub$utalter, F4.sub$ucsex, F4.sub$utbmi, F4.sub$utalkkon))
    rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[3,], test$'Mod2: Y~X+M'[4,]))
  }
}
colnames(rst)[1:4] = c("Indirect.Effect", "SE", "z.value", "N")
colnames(rst)[6:9] = c("S.estimate", "S.SE", "S.tvalue", "S.Pvalue")
colnames(rst)[10:13] = c("mediator.estimate", "mediator.SE", "mediator.tvalue", "mediator.Pvalue")

expression = rep(rownames(F4.expression),each = 361)
expression = cbind(expression, as.character(annotation[expression,2]))
methylation = rep(rownames(F4.methy), 23)
methylation = cbind(methylation, fullannot[methylation,22])
rst = cbind(expression, methylation, rst)
write.csv(rst, file = "Methylation mediated smoking_expression association_scaled_F4.csv", row.names = FALSE)  


## Mediation of expression for the association between smoking and methylation.
rst = NULL
for(e in rownames(F4.expression)){
  F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
  F4.sub$expression = scale(F4.sub$expression) # scale the expression level to make all coefficients comparable
  
  for(cpg in rownames(F4.methy)){
    F4.sub$cpg = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
    F4.sub$cpg = scale(F4.sub$cpg) # scale the methlation level to make all coefficients comparable
    test = sobel.lm2(pred = F4.sub$my.cigreg, med = F4.sub$expression, out = F4.sub$cpg, covariates=data.frame(F4.sub$utalter, as.factor(F4.sub$ucsex), F4.sub$utbmi, F4.sub$utalkkon))
    rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[3,], test$'Mod2: Y~X+M'[4,]))
  }
}
colnames(rst)[1:4] = c("Indirect.Effect", "SE", "z.value", "N")
colnames(rst)[6:9] = c("S.estimate", "S.SE", "S.tvalue", "S.Pvalue")
colnames(rst)[10:13] = c("mediator.estimate", "mediator.SE", "mediator.tvalue", "mediator.Pvalue")

expression = rep(rownames(F4.expression),each = 361)
expression = cbind(expression, as.character(annotation[expression,2]))
methylation = rep(rownames(F4.methy), 23)
methylation = cbind(methylation, fullannot[methylation,22])
rst = cbind(expression, methylation, rst)
write.csv(rst, file = "Expression mediated smoking_methylation association_F4.csv", row.names = FALSE)  

## Mediation of cpg or gene for the association between smoking and metabolites.
F4.sub = subset(F4, expr_in_F4!="" & !is.na(zz_nr_f4_meth) & !is.na(zz_nr_f4_bio))
for(m in metabo.asso){
  rst = NULL
  F4.sub$metabo = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio), m]))
  
  
  F4.sub = subset(F4, expr_in_F4!="" & !is.na(zz_nr_f4_bio))
  for(e in rownames(F4.expression)){
    F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
    test = sobel.lm2(pred = F4.sub$my.cigreg, med = F4.sub$metabo, out = F4.sub$expression, covariates=data.frame(F4.sub$utalter, F4.sub$ucsex, F4.sub$utbmi, F4.sub$utalkkon))
    rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[3,], test$'Mod2: Y~X+M'[4,]))
  }
  
  
  F4.sub = subset(F4, !is.na(zz_nr_f4_meth) &!is.na(zz_nr_f4_bio))
  for(cpg in rownames(F4.methy)){
    F4.sub$cpg = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
    test = sobel.lm2(pred = F4.sub$my.cigreg, med = F4.sub$metabo, out = F4.sub$cpg, covariates=data.frame(F4.sub$utalter, F4.sub$ucsex, F4.sub$utbmi, F4.sub$utalkkon))
    rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[3,], test$'Mod2: Y~X+M'[4,]))
  }
  
  colnames(rst)[1:4] = c("Indirect.Effect", "SE", "z.value", "N")
  colnames(rst)[6:9] = c("S.estimate", "S.SE", "S.tvalue", "S.Pvalue")
  colnames(rst)[10:13] = c("mediator.estimate", "mediator.SE", "mediator.tvalue", "mediator.Pvalue")
  
  rst = data.frame(mediator = c(rownames(F4.expression),rownames(F4.methy)), rst)
  
  write.csv(rst, paste(m, "_metabolite mediation.csv", sep = ""),row.names=FALSE)
}


##############
##    F3    ##
##          ##
##############
##
## Mediation of methylation/expression for the association between smoking and gene expression.
rst = NULL
for(e in rownames(F3.expression)){
  F3.sub$expression = F3.expression[e,as.character(F3.sub$zz_nr_f3_genexp)]
  F3.sub$expression = scale(F3.sub$expression)
  
  for(cpg in rownames(F3.methy)){
    F3.sub$cpg = t(F3.methy[cpg, as.character(F3.sub$zz_nr_f3_meth)])
    F3.sub$cpg = scale(F3.sub$cpg)
    test = sobel.lm(pred = as.factor(F3.sub$my.cigreg), med = F3.sub$expression, out = F3.sub$cpg, covariates=data.frame(F3.sub$rtalteru, as.factor(F3.sub$rcsex), F3.sub$rtbmi, F3.sub$rtalkkon))
    rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[2,], test$'Mod2: Y~X+M'[3,]))
  }
}
colnames(rst)[1:4] = c("Indirect.Effect", "SE", "z.value", "N")
colnames(rst)[6:9] = c("S.estimate", "S.SE", "S.tvalue", "S.Pvalue")
colnames(rst)[10:13] = c("mediator.estimate", "mediator.SE", "mediator.tvalue", "mediator.Pvalue")

expression = rep(rownames(F3.expression),each = 361)
expression = cbind(expression, as.character(annotation[expression,2]))
methylation = rep(rownames(F3.methy), nrow(F3.expression))
methylation = cbind(methylation, fullannot[methylation,22])
rst = cbind(expression, methylation, rst)
#write.csv(rst, file = "Methylation mediated smoking_expression association_scaled_F3.csv", row.names = FALSE)  

write.csv(rst, file = "Expression mediated smoking_methylation association_scaled_F3.csv", row.names = FALSE) 


