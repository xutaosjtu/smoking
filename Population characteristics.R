#Population characteristics

population = function(data, subset, feature.continue, feature.discret){
  ## create working data set
  data = data[subset, ]
  
  characteristics = NULL
  for(i in feature.continue){
    M = tapply(data[,i], INDEX = data[, "my.cigreg"], mean, na.rm = T)
    S = tapply(data[,i], INDEX = data[, "my.cigreg"], sd, na.rm=T)
    tmp.test = pairwise.wilcox.test(data[, i], data[, "my.cigreg"])
    p = tmp.test$p.value[,1]
    characteristics = rbind(characteristics, c(paste(round(M,2),"(", round(S,2), ")", sep = ""), p))
  }
  for(i in feature.discret){
    tmp.table = table(data[,i], data[,"my.cigreg"])
    p = fisher.test(tmp.table[,1:2])$p.value
    if(ncol(tmp.table)==3){
      p2 = fisher.test(tmp.table[,c(1,3)])$p.value
      p = c(p, p2)
    }
    characteristics = rbind(characteristics, c(tmp.table[2,], p))
  }
  rownames(characteristics) = c(feature.continue, feature.discret)
  return(characteristics)
}

## F3
F3.subset = which(!is.na(F3$Arg_PTC))#!is.na(F3$zz_nr_f3_meth) & 

F3.feature.continue = scan(what = character())
rtalteru
rtbmi
rtalkkon

F3.feature.discret = scan(what = character())
rcsex
my.diab
my.alkkon

tmp = population(data=F3, subset = F3.subset, F3.feature.continue, F3.feature.discret)
write.csv(tmp, file = "F3 population characteristics_metabolite.csv")

##F4
F4.subset = which(!is.na(F4$zz_nr_f4_bio))#& !is.na(F4$zz_nr_f4_meth) & !is.na(F4$zz_nr_f4_bio)

F4.feature.continue = scan(what = character())
utalteru
utbmi
utalkkon

F4.feature.discret = scan(what = character())
ucsex
my.diab
my.alkkon

tmp = population(data=F4, subset = F4.subset, F4.feature.continue, F4.feature.discret)
write.csv(tmp, file = "F4 population characteristics_metabolite.csv")

