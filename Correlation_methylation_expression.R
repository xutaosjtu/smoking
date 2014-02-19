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
png("correlation of cpg sites with gene expression_only smoker.png", width = 20, height = 20, unit = "in", res = 180)
par(mfrow = c(3, 3))
F4.sub = subset(F4.sub, my.cigreg==2)
for (i in 1:length(cpg.candidate)){
  methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
  expr = as.vector(F4.expression["ILMN_1773650", as.character(F4.sub$zz_nr_s4f4_genexp)])
  plot(y = methy, 
       x = expr, 
       xlab = "LRRN3 expression",
       ylab = cpg.candidate[i], 
       col = c("green", "orange","red")[F4.sub$my.cigreg+1], pch = 19
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
rst = NULL
for(i in 1:length(cpg.candidate)){
  F4.sub$methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
  F4.sub$expr = as.vector(F4.expression["ILMN_1773650", as.character(F4.sub$zz_nr_s4f4_genexp)])
  model = lm(methy ~ expr + utalteru + ucsex + utbmi + utalkkon + as.factor(utcigreg), data = F4.sub)
  rst = rbind(rst,c(coef(model)[2], confint(model)[2,], summary(model)$coef[1,4]))
}
rownames(rst) = cpg.candidate

## check wheter the association between the two cpg sites were due to smoking
F4.sub$methy1 = t(F4.methy["cg00501876", as.character(F4.sub$zz_nr_f4_meth)])
F4.sub$methy2 = t(F4.methy["cg09837977", as.character(F4.sub$zz_nr_f4_meth)])
## a model with adjustment of smoking
model.1 = lm(methy1 ~ methy2 + utalteru + ucsex + utbmi + utalkkon + as.factor(utcigreg), data = F4.sub)
## a model without adjustment of smoking
model.2 = update(model.1, .~. - as.factor(utcigreg))
summary(model.1)
summary(model.2)


## at different category of pack years
tapply(F4.sub$utpyrs, INDEX = F4.sub$my.cigreg, mean, na.rm = T)
tapply(F4.sub$utpyrs_ai, INDEX = F4.sub$my.cigreg, mean, na.rm = T)
hist(F4.sub$utpyrs)

boxplot(expr ~ utcigreg_sf, F4.sub, ylim = c(-2, 6), main = "LRRN3 expression and cg09837977 methylation")
boxplot(methy ~ utcigreg_sf, F4.sub, col = "red", add = T)
legend("bottomleft", legend = c("expression", "methylation"), pch = 19, col = c("black", "red"))

##

boxplot(methy~packy_categ, F4.sub, col = "red", ylim = c(-2,6))
boxplot(expr~packy_categ, F4.sub, add = T)
legend("bottomleft", legend = c("expression", "methylation"), pch = 19, col = c("black", "red"))

boxplot(methy~uc075_1, F4.sub,  col = "red",ylim = c(-2,6))
boxplot(expr~uc075_1, F4.sub, add=T)

plot(methy~utpyrs_ai, F4.sub, subset = (F4.sub$my.cigreg==2), ylim = c(-2,6))
lines(lowess(methy~utpyrs_ai, F4.sub, subset = (F4.sub$my.cigreg==2)))
points(expr~utpyrs_ai, F4.sub, subset = (F4.sub$my.cigreg==2), col = "red")
lines(lowess(expr~utpyrs_ai, F4.sub, subset = (F4.sub$my.cigreg==2)), col = "red")

##
rst = NULL
for(i in 1:length(cpg.candidate)){
  F4.sub$methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
  F4.sub$expr = as.vector(F4.expression["ILMN_1773650", as.character(F4.sub$zz_nr_s4f4_genexp)])
  model = lm(methy ~ expr + utalteru + ucsex + utbmi + utalkkon + as.factor(utcigreg), data = F4.sub)
  rst = rbind(rst,c(coef(model)[2], confint(model)[2,], summary(model)$coef[1,4]))
}
rownames(rst) = cpg.candidate

## Interaction of different category of packyears on expression association with methylation
rst = NULL
for(i in 1:length(cpg.candidate)){
  F4.sub$methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
  F4.sub$expr = as.vector(F4.expression["ILMN_1773650", as.character(F4.sub$zz_nr_s4f4_genexp)])
  model = lm(methy ~ expr*packy_categ + utalteru + ucsex + utbmi + utalkkon + as.factor(utcigreg), data = F4.sub)
  rst = rbind(rst,c(coef(model)[2], confint(model)[2,], summary(model)$coef[1,4]))
}
rownames(rst) = cpg.candidate


pdf("Smoking effects(estimate)_quit years_7 methy LRRN3_both data.pdf")
#F4.sub = subset(F4, !is.na(zz_nr_s4f4_genexp) & expr_in_F4!="")
F4.sub$expr = as.vector(F4.expression["ILMN_1773650", as.character(F4.sub$zz_nr_s4f4_genexp)])
model = lm(expr ~ packy_categ + utalteru + ucsex + utbmi + utalkkon, data = F4.sub, subset = F4.sub$my.cigreg!=1)
plotCI(x=c(0, coef(model)[2:6]),uiw =c(0,  confint(model)[2:6,2]-coef(model)[2:6]), 
       pch = 19,
       ylim = c(-2,2), xlim = c(0.7,6.5),
       main = #"LRRN3 and 2 methylation sites with trans regulation effect",
       "LRRN3 and 7 methylation sites with bidirectional mediation effect",
       #labels = levels(F4.sub$packy_categ)
       xaxt = "n",
       xlab = "Smoking packyears", ylab = "Association coefficients")

#F4.sub = subset(F4, !is.na(zz_nr_f4_meth))
for(i in 1:length(cpg.candidate)){
  F4.sub$methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
  model = lm(methy ~ quity_categ + utalteru + ucsex + utbmi + utalkkon, data = F4.sub, subset = F4.sub$my.cigreg!=2)
  tmp = c(0, coef(model)[2:6])
  #tmp = 2^tmp/(2^tmp + 1)
  plotCI(x =tmp,uiw =c(0,  confint(model)[2:6,2]-coef(model)[2:6]), 
         col = 1+i, pch = 19,
         add = T
  )
  axis(1, 1:6, levels(F4.sub$quity_categ), col.axis = "blue")
}
abline(h = 0, lty = 2)
legend("topright", legend = c("ILMN_1773650", cpg.candidate), pch = 19, col = c("black", 1+1:length(cpg.candidate)))
dev.off()


pdf("Adjusted methylation and expression level_smoking years_7 methy LRRN3_whole population.pdf")
F4.sub = subset(F4, !is.na(zz_nr_s4f4_genexp) & expr_in_F4!="")
F4.sub$expr = as.vector(F4.expression["ILMN_1773650", as.character(F4.sub$zz_nr_s4f4_genexp)])
methy * model = lm(expr ~ utalteru + ucsex + utbmi + utalkkon, data = F4.sub)
F4.sub$expr.resid=NA
F4.sub$expr.resid[-model$na.action]= model$residuals
plotmeans(expr.resid ~ quity_categ, data = F4.sub, subset = F4.sub$my.cigreg!=2,col = "black",barcol="black", pch = 19, connect=FALSE)

F4.sub = subset(F4, !is.na(zz_nr_f4_meth))
for(i in 1:length(cpg.candidate)){
  F4.sub$methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
  model = lm(scale(methy) ~ utalteru + ucsex + utbmi + utalkkon, data = F4.sub)
  F4.sub$methy.resid=NA
  F4.sub$methy.resid[-model$na.action]= model$residuals
  plotmeans(methy.resid ~ quity_categ, data = F4.sub, subset = F4.sub$my.cigreg!=2,col = 1+i,barcol=1+i, pch = 19, connect=FALSE, add = T)
}

abline(h = 0, lty = 2)
legend("topright", legend = c("ILMN_1773650", cpg.candidate), pch = 19, col = c("black", 1:length(cpg.candidate)+1))
dev.off()


rst = list()
F4.sub$expr = as.vector(F4.expression["ILMN_1773650", as.character(F4.sub$zz_nr_s4f4_genexp)])
for(i in 1:7){
  F4.sub$methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
  model = lm(expr ~ methy * as.factor(my.cigreg) + utalter + ucsex + utbmi + utalkkon, F4.sub)
  estimate = summary(model)$coef
  rst[[i]] = estimate
}
names(rst) = cpg.candidate
rst.interaction = rst
rst.adj = rst

sapply(rst.interaction, 
       function(x){
         return(x[10,4]<0.05)
       }
)

plot(expr ~ methy, F4.sub, col = c("green", "orange","red")[F4.sub$my.cigreg+1], pch = 19)

summary(lm(methy~ as.factor(my.cigreg) + utalter + ucsex + utbmi + utalkkon, F4.sub))


F4.sub$methy = t(F4.methy[cpg.candidate[i],as.character(F4.sub$zz_nr_f4_meth)])
model = lm(methy ~ as.factor(my.cigreg) + utalter + ucsex + utbmi + utalkkon, F4.sub)


## F3
F3.sub = subset(F3, !is.na(zz_nr_f3_meth) & !is.na(zz_nr_f3_genexp))
dim(F3.sub)
F3.sub$expr = t(F3.expression["ILMN_1773650", as.character(F3.sub$zz_nr_f3_genexp)])
F3.sub$methy = t(F3.methy[cpg.candidate[i], as.character(F3.sub$zz_nr_f3_meth)])

png("correlation of cpg sites with gene expression_smoker_F3.png", width = 20, height = 20, unit = "in", res = 180)
par(mfrow = c(3, 3))
#F3.sub = subset(F3.sub, my.cigreg==2)
for (i in 1:length(cpg.candidate)){
  methy = t(F3.methy[cpg.candidate[i],as.character(F3.sub$zz_nr_f3_meth)])
  expr = t(F3.expression["ILMN_1773650", as.character(F3.sub$zz_nr_f3_genexp)])
  plot(y = methy, 
       x = expr, 
       xlab = "LRRN3 expression",
       ylab = cpg.candidate[i], 
       col = c("green", "orange","red")[F3.sub$my.cigreg+1], pch = 19
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

F3.sub
