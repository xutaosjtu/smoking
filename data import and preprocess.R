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
rownames(F4.metab)=F4.metab$sample.id
colnames(F4)
length(which(F4$sample.id %in% S4.F4$zz_nr_f4_bio))
colnames(S4);dim(S4)
length(which(S4$S4Metabo_ZZ %in% S4.F4$zz_nr_s4_bio))

## alcohol consumption
F3$my.alkkon = rep(0, nrow(F3))
F3$my.alkkon[which(F3$rtalkkon >=40 & F3$rcsex==1)] = 1
F3$my.alkkon[which(F3$rtalkkon >=20 & F3$rcsex==2)] = 1
#F3$my.alkkon = as.factor(F3$my.alkkon)

S4$my.alkkon = rep(0, dim(S4)[1])
S4$my.alkkon[which(S4$ltalkkon >=40 & S4$lcsex==1 )] = 1
S4$my.alkkon[which(S4$ltalkkon >=20 & S4$lcsex==2 )] = 1
#S4$my.alkkon = as.factor(S4$my.alkkon)

F4$my.alkkon = rep(0, dim(F4)[1])
F4$my.alkkon[which(F4$utalkkon >=40 & F4$ucsex==1 )] = 1
F4$my.alkkon[which(F4$utalkkon >=20 & F4$ucsex==2 )] = 1
#F4$my.alkkon = as.factor(F4$my.alkkon)

## smoking
F3$my.cigreg = F3$rtcigreg
F3$my.cigreg[which(F3$rtcigreg == 1)]=2
F3$my.cigreg = 4-F3$my.cigreg
#F3$my.cigreg = as.factor(F3$my.cigreg)

S4$my.cigreg = S4$ltcigreg
S4$my.cigreg[which(S4$ltcigreg==1)]=2
S4$my.cigreg = 4-S4$my.cigreg
#S4$my.cigreg = as.factor(S4$my.cigreg)

F4$my.cigreg = F4$utcigreg_sf
F4$my.cigreg[which(F4$utcigreg_sf==1)]=2
F4$my.cigreg = 4-F4$my.cigreg
#F4$my.cigreg = as.factor(F4$my.cigreg)

F4$quit_year = F4$uc075_3
F4$quit_year[which(F4$quit_year==-999)] = NA
F4$quit_year = F4$utalteru - F4$quit_year
F4$quit_year[which(F4$quit_year<1)]=1
F4$quit_year[which(F4$uc075_1 %in% c(1:4))] = 0

F4$smoke_year = F4$uc069_2
F4$smoke_year[which(F4$smoke_year== -999)] = NA
F4$smoke_year = F4$utalteru - F4$smoke_year

F4$packy_categ = cut(F4$utpyrs_ai, quantile(F4$utpyrs_ai[which(F4$my.cigreg!=0)], probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
#F4$packy_categ = cut(F4$smoke_year, quantile(F4$smoke_year[which(F4$my.cigreg!=0)], probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
tmp = levels(F4$packy_categ)
F4$packy_categ = as.character(F4$packy_categ)
F4$packy_categ[which(F4$my.cigreg==0)] = "NS"
F4$packy_categ = factor(F4$packy_categ, levels = c("NS", tmp))

F4$quity_categ = cut(F4$quit_year, quantile(F4$quit_year[which(F4$my.cigreg!=0)], probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
tmp = levels(F4$quity_categ)
F4$quity_categ = as.character(F4$quity_categ)
F4$quity_categ[which(F4$my.cigreg==0)] = "NS"
F4$quity_categ = factor(F4$quity_categ, levels = c("NS", tmp))


table(F4$packy_categ, useNA="always")
table(F4$my.cigreg, useNA = "alwa")

## diabetes
F3$my.diab = F3$rtdiabet

S4$my.diab = S4$utdiabet

F4$my.diab=F4$utdiabet

## physical activity
F4$my.physical = F4$utphact
F4$my.physical[which(F4$utphact<=2)]=1
F4$my.physical[which(F4$utphact>2)]=0


