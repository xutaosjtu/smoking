require(lumi)
## Importing metabolomics data
F3 = read.csv("metabolites/K2912_Xu_F3_trans20130304.csv", sep=";")
S4.F4 = read.csv("metabolites/K2912_Xu_S4F4_trans120213.csv", sep=";")
S4 = read.csv(file="metabolites/_20120716_imputed_Biocr_S4_ZZ_NR.csv")
F4 = read.csv(file="metabolites/KoraF4-metabolomics-quality.controlled-mice.imputed-20100107.csv",sep=";")
colnames(F4)
length(which(F4$sample.id %in% S4.F4$zz_nr_f4_bio))
colnames(S4);dim(S4)
length(which(S4$S4Metabo_ZZ %in% S4.F4$zz_nr_s4_bio))

##importing expression data
S4_expression = read.csv(file="expression/KORA_S4_657_samples_Expressionsdaten.csv")
colnames(S4_expression)=sapply(colnames(S4_expression),function(x) unlist(strsplit(x,split="X", fixed=T))[2])
colnames(S4_expression)[1]="probeid"
#colnames(S4_expression)[2]
sum(colnames(S4_expression) %in% S4.F4$zz_nr_s4f4_genexp)
dim(S4_expression)
F4_expression = read.csv(file = "expression/KORA_F4_993_samples_Expressionsdaten.csv")
#colnames(F4_expression)[1:5]
colnames(F4_expression)=sapply(colnames(F4_expression),function(x) unlist(strsplit(x,split="X",fixed = T))[2])
colnames(F4_expression)[1]="probeid"
sum(colnames(F4_expression) %in% S4.F4$zz_nr_s4f4_genexp)
F3_expression = read.csv(file = "expression/KORA_F3_381_samples_normalized_expression_data.csv")
#dim(F4_expression)
#colnames(F3_expression)[1:5]
colnames(F3_expression)=sapply(colnames(F3_expression),function(x) unlist(strsplit(x,split="X",fixed = T))[2])
colnames(F3_expression)[1]="probeid"
sum(colnames(F3_expression) %in% F3$zz_nr_f3_genexp)

data = read.table("methylation/F3/KORA.F3_methylation_m.value_callrate80_chr_1_.txt",sep=" ")
colnames(data)=sapply(colnames(data),function(x) unlist(strsplit(x,split="X",fixed = T))[2])
