F4.sub2 = F4[which(!is.na(F4$zz_nr_f4_bio)),]

rst = NULL
for(m in colnames(F4.metab)[-1]){
  F4.sub2$metabolite = scale(log(F4.metab[as.character(F4.sub2$zz_nr_f4_bio), m]))
  model = lm(metabolite ~ as.factor(my.cigreg)
             + utalteru + as.factor(ucsex) + utbmi + as.factor(my.alkkon) + as.factor(my.diab)
             ,data = F4.sub2
             )
  rst = rbind(rst, summary(model)$coef[3, ])
}
rownames(rst) = metabo.valid
rst = data.frame(rst, fdr = p.adjust(rst[,4], method = "BH"), bonferroni = p.adjust(rst[,4], method = "bonf"))

metabo.asso2 = subset(metabo.valid, subset = rst$fdr<0.05)
