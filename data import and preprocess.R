require(lumi)
require(doMC)
## Importing metabolomics data
F3 = read.csv("metabolites/K2912_Xu_F3_trans20130304.csv", sep=";")
rownames(F3)=F3$ZZ_nr
S4.F4 = read.csv("metabolites/K2912_Xu_S4F4_trans120213.csv", sep=";")
S4 = S4.F4[, c(1:33,64:70)]
F4 = S4.F4[, c(34:70)]
S4.meta = read.csv(file="metabolites/_20120716_imputed_Biocr_S4_ZZ_NR.csv")
F4.meta = read.csv(file="metabolites/KoraF4-metabolomics-quality.controlled-mice.imputed-20100107.csv",sep=";")
rownames(F4)=F4$sample.id
colnames(F4)
length(which(F4$sample.id %in% S4.F4$zz_nr_f4_bio))
colnames(S4);dim(S4)
length(which(S4$S4Metabo_ZZ %in% S4.F4$zz_nr_s4_bio))

##importing expression data
#S4_expression = read.csv(file="expression/KORA_S4_657_samples_Expressionsdaten.csv")
#colnames(S4_expression)=sapply(colnames(S4_expression),function(x) unlist(strsplit(x,split="X", fixed=T))[2])
#colnames(S4_expression)[1]="probeid"
#sum(colnames(S4_expression) %in% S4.F4$zz_nr_s4f4_genexp)
#dim(S4_expression)
#F4_expression = read.csv(file = "expression/KORA_F4_993_samples_Expressionsdaten.csv")
#colnames(F4_expression)[1:5]
#colnames(F4_expression)=sapply(colnames(F4_expression),function(x) unlist(strsplit(x,split="X",fixed = T))[2])
#colnames(F4_expression)[1]="probeid"
#sum(colnames(F4_expression) %in% S4.F4$zz_nr_s4f4_genexp)

F3.expression = read.csv(file = "expression/KORA_F3_379_samples_normalized_expression_data_zzupdated.csv", row.names=1)
colnames(F3.expression)=sapply(colnames(F3.expression),function(x) unlist(strsplit(x,split="X",fixed = T))[2])
sum(colnames(F3.expression) %in% F3$zz_nr_f3_genexp)

#data = read.table("methylation/F3/KORA.F3_methylation_m.value_callrate80_chr_1_.txt",sep=" ")
#colnames(data)=sapply(colnames(data),function(x) unlist(strsplit(x,split="X",fixed = T))[2])


## importing expression data
load("expression/Expression_s4_adjusted_technical_variables.Rdata")
load("expression/Expression_f4_adjusted_technical_variables.Rdata")
F4.expression = adj.expr.f4
S4.expression = adj.expr.s4
rm(adj.expr.s4)
rm(adj.expr.f4)

## alcohol consumption
F3$my.alkkon = rep(0, nrow(F3))
F3$my.alkkon[which(F3$rtalkkon >=40 & F3$rcsex==1)] = 1
F3$my.alkkon[which(F3$rtalkkon >=20 & F3$rcsex==2)] = 1

S4$my.alkkon = rep(0, dim(S4)[1])
S4$my.alkkon[which(S4$ltalkkon >=40 & S4$lcsex==1 )] = 1
S4$my.alkkon[which(S4$ltalkkon >=20 & S4$lcsex==2 )] = 1

F4$my.alkkon = rep(0, dim(F4)[1])
F4$my.alkkon[which(F4$utalkkon >=40 & F4$ucsex==1 )] = 1
F4$my.alkkon[which(F4$utalkkon >=20 & F4$ucsex==2 )] = 1

## smoking
F3$my.cigreg = F3$rtcigreg
F3$my.cigreg[which(F3$rtcigreg == 1)]=2
F3$my.cigreg = 4-F3$my.cigreg

S4$my.cigreg = S4$ltcigreg
S4$my.cigreg[which(S4$ltcigreg==1)]=2
S4$my.cigreg = 4-S4$my.cigreg

F4$my.cigreg = F4$utcigreg
F4$my.cigreg[which(F4$utcigreg==1)]=2
F4$my.cigreg = 4-F4$my.cigreg
#F4$my.cigreg = as.factor(F4$my.cigreg)

## diabetes
F3$my.diab = F3$rtdiabet

S4$my.diab = S4$utdiabet

F4$my.diab=F4$utdiabet

require(doMC)
registerDoMC(cores = 4)

data = F3[!is.na(F3$zz_nr_f3_genexp),]
rownames(data) = data$zz_nr_f3_genexp
data = data[colnames(F3.expression),]
association = foreach(i = 1:nrow(F3.expression), .combine = rbind ) %dopar% {
              data$expr = t(F3.expression[i,]) 
              model = glm(expr ~ as.factor(my.cigreg == 2)+ as.factor(my.cigreg == 1)
                + rtalter + as.factor(rcsex)
                + as.factor(my.alkkon) + rtbmi #+ as.factor(rtdiabet)
                , data = data
          #, family = binomial
               )
        return(c(summary(model)$coefficients[2,], summary(model)$coefficients[3,]))
}
rownames(association)=rownames(F3.expression)
write.csv(association, file = "F3 smoking associated genes_lm_combine FS and NS_diab unadj.csv")


data = S4[S4$expr_in_S4!="",]
rownames(data) = data$zz_nr_s4f4_genexp
data = data[colnames(S4.expression),]
association = foreach(i = 1:nrow(S4.expression), .combine = rbind ) %dopar% {
        data$expr = S4.expression[i,]
        model = glm(expr ~ as.factor(my.cigreg == 2) #+ as.factor(my.cigreg == 1)
          + ltalter + as.factor(lcsex)
          + as.factor(my.alkkon) + ltbmi #+ as.factor(ltdiabet)
          , data = data
          #, family = binomial
          )
        return(c(summary(model)$coefficients[2,], summary(model)$coefficients[3,]))
      }
rownames(association)=rownames(S4.expression)
write.csv(association, file = "S4 smoking associated genes_lm_combine FS and NS_diab unadj.csv")

data = F4[(F4$expr_in_F4!=""), ]
rownames(data) = data$zz_nr_s4f4_genexp
data = data[colnames(F4.expression),]
association = foreach(i = 1:nrow(F4.expression), .combine = rbind ) %dopar% {
        data$expr = F4.expression[i,]
        model = glm(expr ~ as.factor(my.cigreg ==2 ) + as.factor(my.cigreg == 1)
          + utalter + as.factor(ucsex)
          + as.factor(my.alkkon) + utbmi #+ as.factor(utdiabet)
          , data = data
         # , family = binomial
          )
        return(c(summary(model)$coefficients[2,], summary(model)$coefficients[3,]))
      }
rownames(association) = rownames(F4.expression)
write.csv(association, file = "F4 smoking associated genes_lm_combine FS and NS_diab unadj.csv")


