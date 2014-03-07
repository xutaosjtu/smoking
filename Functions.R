metabo_methy_asso = function(cpg, m){
  F4.sub = subset(F4, !is.na(zz_nr_f4_meth) &!is.na(zz_nr_f4_bio))
  
  F4.sub$metabo = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio), m]))
  
  F4.sub$methy = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
  
  model = lm(methy ~ metabo  + as.factor(my.cigreg)
             + utalteru + as.factor(ucsex)
             + as.factor(my.alkkon) + utbmi #+ as.factor(rtdiabet)
             #+ CD8T + CD4T + NK + Bcell + Mono + Gran # adjust for WBC
             , data = F4.sub
  )
  
  return(model)
  
}


metabo_methy_expr_asso = function(cpg, m, e){
  F4.sub = subset(F4, expr_in_F4!="" & !is.na(zz_nr_f4_meth) & !is.na(zz_nr_f4_bio))
  
  F4.sub$metabo = scale(log(F4.metab[as.character(F4.sub$zz_nr_f4_bio), m]))
  
  F4.sub$methy = t(F4.methy[cpg, as.character(F4.sub$zz_nr_f4_meth)])
  
  F4.sub$expr = F4.expression[e,as.character(F4.sub$zz_nr_s4f4_genexp)]
  
  model = lm(metabo ~  expr + methy  #+ as.factor(my.cigreg)
             + utalteru + as.factor(ucsex)
             + as.factor(my.alkkon) + utbmi #+ as.factor(rtdiabet)
             #+ CD8T + CD4T + NK + Bcell + Mono + Gran # adjust for WBC
             , data = F4.sub
  )
  
  return(model)
}
