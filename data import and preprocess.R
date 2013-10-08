require(lumi)
require(doMC)
## Importing metabolomics data
F3 = read.csv("metabolites/K2912_Xu_F3_trans20130304.csv", sep=";")
rownames(F3)=F3$ZZ_nr
S4.F4 = read.csv("metabolites/K2912_Xu_S4F4_trans120213.csv", sep=";")
S4 = S4.F4[, c(1:33,64:70)]
F4 = S4.F4[, c(34:70)]
S4.metab = read.csv(file="metabolites/_20120716_imputed_Biocr_S4_ZZ_NR.csv")
F4.metab = read.csv(file="metabolites/KoraF4-metabolomics-quality.controlled-mice.imputed-20100107.csv",sep=";")
rownames(F4)=F4$sample.id
colnames(F4)
length(which(F4$sample.id %in% S4.F4$zz_nr_f4_bio))
colnames(S4);dim(S4)
length(which(S4$S4Metabo_ZZ %in% S4.F4$zz_nr_s4_bio))

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

