require(gplots)
require(ggplot2)

cpg.candidate = scan(what = character())
cg23771366
cg00073090
cg02532700
cg00501876
cg14753356
cg26729380
cg09837977

cg21280392
cg12593793


gene.candidate = scan(what = character())
ILMN_1773650
ILMN_1710326

F4.sub = subset(F4, !is.na(zz_nr_f4_meth) & !is.na(zz_nr_s4f4_genexp) & expr_in_F4!="")

## 
png("correlation of cpg sites with gene expression_2.png", width = 20, height = 20, unit = "in", res = 180)
par(mfrow = c(3, 3))
for (i in 1:length(cpg.candidate)){
  methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
  expr = as.vector(F4.expression["ILMN_1773650", as.character(F4.sub$zz_nr_s4f4_genexp)])
  plot(x = methy, 
       y = expr, 
       xlab = "LRRN3",
       ylab = cpg.candidate[i]
       )
  correlation = cor.test(methy, expr)
  legend("topright", 
         legend = c(paste("corr = ", format(correlation$estimate, digits = 2), "(", format(correlation$conf.int[1], digits = 2), ",", format(correlation$conf.int[2], digits = 2),")", sep=""),
           paste("P =", format(correlation$p.value,digits = 4, scientific = TRUE))
           ),
         cex = 2
  )
}
dev.off()

## correlation among all genes and methylation sites
png("Scatter plot matrix of the cpg sites and expression_2.png", width = 20, height = 20, unit = "in", res = 180)
data = cbind(t(F4.methy[cpg.candidate, as.character(F4.sub$zz_nr_f4_meth)]), t(F4.expression[gene.candidate,as.character(F4.sub$zz_nr_s4f4_genexp)]))
pairs(~., data)
dev.off()
require(corrplot)
corrplot(cor(data, use = "pair"),
         method = "square",
         type = "lower", diag = T, tl.pos = "l",
         col = bluered(200))
corrplot(cor(data, use = "pair"), add = T,
         method = "number",
         type = "upper", diag = F, tl.pos = "n", cl.pos = "n",
         col = bluered(200), addCoef.col = F)


## association in the linear models with adjustment of covariates
for(i in 1:length(cpg.candidate)){
  
}

## check wheter the association between the two cpg sites were due to smoking
F4.sub$methy1 = t(F4.methy["cg00501876", as.character(F4.sub$zz_nr_f4_meth)])
F4.sub$methy2 = t(F4.methy["cg09837977", as.character(F4.sub$zz_nr_f4_meth)])
## a model with adjustment of smoking
model.1 = lm(methy1 ~ methy2 + utalteru + ucsex + utbmi + utalkkon + as.factor(utcigreg), data = F4.sub)
## a model without adjustment of smoking
model.2 = update(model.1, .~. - as.factor(utcigreg))
summary(model.1)
summary(model.2)



