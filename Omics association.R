setwd("H:/Smoking and methylation")

## import metabolomics data
F3 = read.csv("metabolites/K2912_Xu_F3_trans20130304.csv", sep=";")
rownames(F3)=F3$ZZ_nr
S4.F4 = read.csv("metabolites/K2912_Xu_S4F4_trans120213.csv", sep=";")
S4 = S4.F4[, c(1:33,64:70)]
F4 = S4.F4[, c(34:70)]
#S4.metab = read.csv(file="metabolites/_20120716_imputed_Biocr_S4_ZZ_NR.csv")
F4.metab = read.csv(file="metabolites/KoraF4-metabolomics-quality.controlled-mice.imputed-20100107.csv",sep=";", row.names = 1)
rownames(F4)=F4$sample.id
colnames(F4)
length(which(F4$sample.id %in% S4.F4$zz_nr_f4_bio))
colnames(S4);dim(S4)
length(which(S4$S4Metabo_ZZ %in% S4.F4$zz_nr_s4_bio))

## import gene expression data
load("intermediate results/expression candidates.RData")
colnames(F3.expression) = gsub("X", "", colnames(F3.expression))

## import methylation data
F4.methy = read.csv("intermediate results/Methylation_WBC adjust/methylation data of candidates_F4.csv", row.names = 1)
F3.methy = read.csv("intermediate results/Methylation_WBC adjust/methylation data of candidates_F3.csv", row.names = 1)
colnames(F4.methy) = substr(colnames(F4.methy), 2, 10)
colnames(F3.methy) = substr(colnames(F3.methy), 2, 10)
load("metabolites/smoking associated metabolites.RData")

## Subset with all three types of data
F4.sub = subset(F4, expr_in_F4!="" & !is.na(zz_nr_f4_meth) & !is.na(zz_nr_f4_bio))

#F4.sub = subset(F4, !is.na(zz_nr_f4_meth) & !is.na(zz_nr_f4_bio))

colnames(F3)[27:189] = gsub("_PTC", "", colnames(F3)[27:189])
colnames(F3)[27:189] = gsub("_", ".", colnames(F3)[27:189])
colnames(F3)[86] = "SM..OH..C22.2"
F3.sub = subset(F3, !is.na(zz_nr_f3_meth) & !is.na(zz_nr_f3_genexp))


## expression ----- metabolites
metabo.valid = colnames(F4.metab)
rst = NULL;
for(m in metabo.valid){
  F4.sub$metabolite = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio),m]))
  tmp =NULL
  for(e in 1:nrow(F4.expression)){
    F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]  
    model = lm( expression ~  metabolite
        + as.factor(my.cigreg ==2 ) + as.factor(my.cigreg == 1)
        + utalteru + as.factor(ucsex)  + utbmi + utalkkon
        #+ as.factor(utdiabet)
        , data = F4.sub
        )
    tmp = c(tmp,summary(model)$coef[2,])
  }
  rst = rbind(rst,tmp)
}
rownames(rst) = metabo.valid
metab.expresion = apply(rst[, 4*(1:nrow(F4.expression))], 2, function(x) which(x<0.05))
names(metab.expresion) = rownames(F4.expression)

metab.expresion = sapply(metab.expresion, function(x) return(intersect(names(x), metabo.asso)))
names(metab.expresion) = annotation[match(names(metab.expresion), annotation[,1]),2]

## methylation -----  expression
cpg.valid = rownames(F4.methy)
rst = NULL;
for(cpg in cpg.valid){
  F4.sub$methy = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
  tmp = NULL
  for(e in 1:nrow(F4.expression)){
    F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
    model = lm(expression ~ methy
               + as.factor(my.cigreg) 
               + utalteru + as.factor(ucsex)  + utbmi + utalkkon
               #+ as.factor(utdiabet)
               , data = F4.sub
               )
    tmp = c(tmp, summary(model)$coef[2,])
  }
  rst = rbind(rst, tmp)
}
rownames(rst) = cpg.valid
expresion.methy = apply(rst[,4*(1:nrow(F4.expression))], 2, function(x) which(x<0.05/nrow(F4.methy)))
names(expresion.methy) = rownames(F4.expression)

## methylation -----  metabolites
cpg.valid = rownames(F4.methy)
metabo.valid = colnames(F4.metab)
rst = NULL;
for(cpg in cpg.valid){
  F4.sub$methy = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
  tmp = NULL
  for(m in metabo.asso){
    F4.sub$metabolite = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio),m]))
    model = lm(methy ~ metabolite
               #+ as.factor(my.cigreg) 
               #+ utalteru + as.factor(ucsex)  + utbmi + utalkkon
               #+ as.factor(utdiabet)
               , data = F4.sub
    )
    tmp = c(tmp, summary(model)$coef[2,])
  }
  rst = rbind(rst, tmp)
}
rownames(rst) = cpg.valid
methy.metab = apply(rst[,4*(1:length(metabo.asso))], 2, function(x) which(x<0.05))
names(methy.metab) = metabo.asso

#length(metabo.asso)

##?? where does the methy.asso come from #####################################
methy.metab = sapply(methy.metab, function(x) return(intersect(names(x), methy.asso)))
##############################################################################

## Test: smoking --x-- expression adjusting the methylation sites
## Result:not possible: smoking is like a gene with pleiotropy;
## Can only see reduced effect size after adjusting the methylation sites.
rst = NULL
for(e in nrow(F4.expression)){
  F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
  model = lm(expression ~ #cg02909929+cg03652339+cg05983405+cg13775629+cg18352710+cg22900360+cg26971585+cg10126923+cg12916723+
               cg22900360 + 
              as.factor(my.cigreg) 
             + utalteru + as.factor(ucsex)  + utbmi + utalkkon
             + as.factor(utdiabet)
             , data = cbind(F4.sub,t(F4.methy[cpg.valid, as.character(F4.sub$zz_nr_f4_meth)]))
  )
  rst = rbind(rst, summary(model)$coef[4,])
}


F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
F4.sub$metabo = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio),m]))
model = lm(metabo ~ #cg02909929+cg03652339+cg05983405+cg13775629+cg18352710+cg22900360+cg26971585+cg10126923+cg12916723+
             #cg22900360 + 
             #metabo + 
             #expression
           + as.factor(my.cigreg) 
           + utalteru + as.factor(ucsex) + utbmi + utalkkon
           #+ as.factor(utdiabet)
           , data = cbind(F4.sub,t(F4.methy[cpg.valid, as.character(F4.sub$zz_nr_f4_meth)]))
)


## convert to adjacency list
toList = function(x){
  return(cbind(
    rep(names(x),sapply(x, length)), 
    unlist(sapply(x, function(y) names(y)))
  )
  )
}

toList2 = function(x){
  return(cbind(
    rep(names(x),sapply(x, length)), 
    unlist(sapply(x, function(y) return(y)))
  )
  )
}

methy.metab.list = toList(methy.metab)
colnames(methy.metab.list) = c("metabolites", "cpg")

metab.expresion.list = toList2(metab.expresion)
colnames(metab.expresion.list) = c("expression","metabolites")

expresion.methy.list = toList(expresion.methy)
colnames(expresion.methy.list) = c("expression","cpg")

tmp = rbind(metab.expresion.list, methy.metab.list, expresion.methy.list)
write.csv(unique(c(tmp[,1], tmp[,2])), "property_no WBC.csv", quote=FALSE)

write.csv(rbind(metab.expresion.list, methy.metab.list, expresion.methy.list), file = "List_associations_no WBC.csv",quote=FALSE)


## Find the candidate association with the pattern methylation--expression--metabolite 
motifs = merge(methy.metab.list, expresion.methy.list)

rst = NULL
for(i in 1:nrow(motifs)){
  e = as.character(motifs$expression[i])
  F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
  m = as.character(motifs$metabolites[i])
  F4.sub$metabo = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio),m]))
  cpg = as.character(motifs$cpg[i])
  F4.sub$cpg = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
  
  model = lm(metabo ~ cpg 
      + expression
      + as.factor(my.cigreg) 
     + utalteru + as.factor(ucsex) + utbmi + utalkkon,
     data = F4.sub)
  
  rst = rbind(rst, c(summary(model)$coef[2,], summary(model)$coef[3,]))
}

motifs_selected = motifs[which(rst[,4]>0.05&rst[,8]<0.05),]
motifs_selected = merge(motifs_selected, annotation, by.x = "expression", by.y = "Probe_Id")
motifs_selected = cbind(motifs_selected,  cpg_gene= apply(motifs_selected, 1, function(x) return(unlist(xx[x[2]]))))


## Mediation analysis
#require(multilevel)
# N <- 100
# X <- rnorm(N, 175, 7)
# M <- 0.7*X + rnorm(N, 0, 5)
# Y <- 0.4*M + rnorm(N, 0, 5)
# dfMed <- data.frame(X, M, Y)

F4.sub$my.cigreg = as.factor(F4.sub$my.cigreg)
F4.sub$ucsex = as.factor(F4.sub$ucsex)

sobel.lm = function(pred, med, out, covariates){
  NEWDAT <- data.frame(pred = pred, med = med, out = out, covariates)
  NEWDAT <- na.exclude(NEWDAT)
  covar = names(covariates)
  colnames(NEWDAT)[1:3] = c("pred", "med", "out")
  model1 <- lm(out ~ . - med, data = NEWDAT)
  model2 <- lm(out ~ ., data = NEWDAT)
  model3 <- lm(med ~ . - out, data = NEWDAT)
  mod1.out <- summary(model1)$coef
  mod2.out <- summary(model2)$coef
  mod3.out <- summary(model3)$coef
  indir <- mod3.out[2, 1] * mod2.out[3, 1]
  effvar <- (mod3.out[2, 1])^2 * (mod2.out[3, 2])^2 + (mod2.out[3, 1])^2 * (mod3.out[2, 2])^2
  serr <- sqrt(effvar)
  zvalue = indir/serr
  out <- list(`Mod1: Y~X` = mod1.out, `Mod2: Y~X+M` = mod2.out, 
              `Mod3: M~X` = mod3.out, Indirect.Effect = indir, SE = serr, 
              z.value = zvalue, N = nrow(NEWDAT))
  return(out)
}

## Mediation of gene expression for the association between methylation and metabolites.
rst = NULL
for(i in 1:nrow(motifs)){
  e = as.character(motifs$expression[i])
  F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
  m = as.character(motifs$metabolites[i])
  F4.sub$metabo = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio),m]))
  cpg = as.character(motifs$cpg[i])
  F4.sub$cpg = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
  
  test = sobel.lm(pred = F4.sub$cpg, med = F4.sub$expression, out = F4.sub$metabo, covariates= data.frame(F4.sub$my.cigreg, F4.sub$utalter, F4.sub$ucsex, F4.sub$utbmi, F4.sub$utalkkon))
  
  rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[2,], test$'Mod2: Y~X+M'[3,]))
}
rst = data.frame(rst[,1:5], fdr = p.adjust(rst[,5], method = "BH"), bonf = p.adjust(rst[,5], method = "bonf"), rst[,6:13])

colnames(rst)[1:5] = c("Indirect.Effect", "SE", "z.value", "N", "P")
colnames(rst)[8:11] = c("cpg.estimate", "cpg.SE", "cpg.tvalue", "cpg.Pvalue")
colnames(rst)[12:15] = c("expr.estimate", "expr.SE", "expr.tvalue", "expr.Pvalue")

write.csv(cbind(motifs, rst), file = "Mediator analysis from the candidates_F4.csv", row.names = FALSE)

## Mediation of cpg or gene for the association between smoking and metabolites.
sobel.lm2 = function(pred, med, out, covariates){
  NEWDAT <- data.frame(pred = pred, med = med, out = out, covariates)
  NEWDAT <- na.exclude(NEWDAT)
  covar = names(covariates)
  colnames(NEWDAT)[1:3] = c("pred", "med", "out")
  model1 <- lm(out ~ . - med, data = NEWDAT)
  model2 <- lm(out ~ ., data = NEWDAT)
  model3 <- lm(med ~ . - out, data = NEWDAT)
  mod1.out <- summary(model1)$coef
  mod2.out <- summary(model2)$coef
  mod3.out <- summary(model3)$coef
  indir <- mod3.out[3, 1] * mod2.out[4, 1]
  effvar <- (mod3.out[3, 1])^2 * (mod2.out[4, 2])^2 + (mod2.out[4, 1])^2 * (mod3.out[3, 2])^2
  serr <- sqrt(effvar)
  zvalue = indir/serr
  out <- list(`Mod1: Y~X` = mod1.out, `Mod2: Y~X+M` = mod2.out, 
              `Mod3: M~X` = mod3.out, Indirect.Effect = indir, SE = serr, 
              z.value = zvalue, N = nrow(NEWDAT))
  return(out)
}

for(m in metabo.asso){
  rst = NULL
  F4.sub$metabo = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio), m]))
  
  for(e in rownames(F4.expression)){
    F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
    test = sobel.lm2(pred = F4.sub$my.cigreg, med = F4.sub$metabo, out = F4.sub$expression, covariates=data.frame(F4.sub$utalter, F4.sub$ucsex, F4.sub$utbmi, F4.sub$utalkkon))
    rst = rbind(rst, c(test$Indirect.Effect, test$SE, test$z.value, test$N, p = 1- pnorm(abs(test$z.value)), test$'Mod2: Y~X+M'[3,], test$'Mod2: Y~X+M'[4,]))
  }
  
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

## Mediation of methylation for the association between smoking and gene expression.
rst = NULL
for(e in rownames(F4.expression)){
  F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
  
  for(cpg in rownames(F4.methy)){
    F4.sub$cpg = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
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
write.csv(rst, file = "Methylation mediated smoking_expression association_F4.csv", row.names = FALSE)  


rst = NULL
for(e in rownames(F4.expression)){
  F4.sub$expression = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
  
  for(cpg in rownames(F4.methy)){
    F4.sub$cpg = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
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

# ## association between the metabolites with smoking in the subset
# rst = NULL
# for(m in metabo.asso){
#   F4.sub$metabo = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio), m]))
#   
#   model = lm(metabo ~ my.cigreg + utalter + ucsex + utbmi + utalkkon, F4.sub)
#   rst = rbind(rst, summary(model)$coef[3,])
# }
# rst = data.frame(rst, fdr = p.adjust(rst[,4], method = "BH"), bonf = p.adjust(rst[,4], method = "bonf"))
# rownames(rst) = metabo.asso
# write.csv(rst, "association between smoking and the metabolite in the subset.csv")
# 
# ## association between the metabolite with smoking in the whole set
# rst = NULL
# for(m in metabo.asso){
#   F4$metabo = scale(log(F4.metab[as.character(F4$zz_nr_f4_bio), m]))
#   
#   model = lm(metabo ~ my.cigreg + utalter + as.factor(ucsex) + utbmi + utalkkon, F4)
#   rst = rbind(rst, summary(model)$coef[3,])
# }
# rst = data.frame(rst, fdr = p.adjust(rst[,4], method = "BH"), bonf = p.adjust(rst[,4], method = "bonf"))
# rownames(rst) = metabo.asso

F4.sub$cpg = t(F4.methy["cg09837977", as.character(F4.sub$zz_nr_f4_meth)])
F4.sub$expression = F4.expression["ILMN_1773650",as.character(F4.sub$zz_nr_s4f4_genexp)]
