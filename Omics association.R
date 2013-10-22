setwd("H:/Smoking and methylation")

## Importing metabolomics data
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

F4.sub = subset(F4, expr_in_F4!="" & !is.na(zz_nr_f4_meth) & !is.na(zz_nr_f4_bio))


## load gene expression data
load("intermediate results/expression candidates.RData")

## import methylation data
F4.methy = read.csv("intermediate results/Methylation/methylation data of candidates_F4.csv", row.names = 1)
F3.methy = read.csv("intermediate results/Methylation/methylation data of candidates_F3.csv", row.names = 1)
colnames(F4.methy) = substr(colnames(F4.methy), 2, 10)
colnames(F3.methy) = substr(colnames(F3.methy), 2, 10)

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
methy.metab = apply(rst[,4*(1:length(metabo.asso))], 2, function(x) which(x<0.05))
names(methy.metab) = metabo.asso

#/length(metabo.asso)

methy.metab = sapply(methy.metab, function(x) return(intersect(names(x), methy.asso)))


## prove smoking --x-- expression adjusting the methylation sites
## not possible: smoking like a gene with pleiotropy;
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
           + utalteru + as.factor(ucsex)  + utbmi + utalkkon
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

methy.metab.list = toList2(methy.metab)

metab.expresion.list = toList2(metab.expresion)

expresion.methy.list = toList(expresion.methy)

tmp = rbind(metab.expresion.list, methy.metab.list, expresion.methy.list)
write.csv(unique(c(tmp[,1], tmp[,2])), "property.csv", quote=FALSE)

write.csv(rbind(metab.expresion.list, methy.metab.list, expresion.methy.list), file = "List_associations.csv",quote=FALSE)

