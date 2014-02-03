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



require(doMC)
registerDoMC(cores = 8)

data = F3[!is.na(F3$zz_nr_f3_genexp),]
rownames(data) = data$zz_nr_f3_genexp
data = data[colnames(F3.expression),]
association = foreach(i = 1:nrow(F3.expression), .combine = rbind ) %dopar% {
  data$expr = t(F3.expression[i,]) 
  model = glm(expr ~ as.factor(my.cigreg)
              + rtalter + as.factor(rcsex)
              + as.factor(my.alkkon) + rtbmi #+ as.factor(rtdiabet)
              , data = data
              #, family = binomial
  )
  return(c(summary(model)$coefficients[2,], summary(model)$coefficients[3,]))
}
rownames(association)=rownames(F3.expression)
write.csv(association, file = "F3 smoking associated genes_lm.csv")


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
  model = glm(expr ~ as.factor(my.cigreg)
              + utalter + as.factor(ucsex)
              + as.factor(my.alkkon) + utbmi #+ as.factor(utdiabet)
              , data = data
              # , family = binomial
  )
  return(c(summary(model)$coefficients[2,], summary(model)$coefficients[3,]))
}
rownames(association) = rownames(F4.expression)
write.csv(association, file = "F4 smoking associated genes_lm.csv")


############################################
## annotation
############################################
annotation.S4F4 = read.csv("expression/Annotation HumanHT-12v3 final.csv", sep = ";", row.names = 1)
annotation.F3 = read.csv("expression/annotation_Illumina_WG_V2_KORA_F3.csv", sep = ";", row.names = 1)