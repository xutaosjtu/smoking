## annotation of the expression data
annotation.S4F4 = read.csv("expression/Annotation HumanHT-12v3 final.csv", sep = ";")
association.S4 = read.csv("S4 smoking associated genes_lm_combine FS and NS.csv")
association.F4 = read.csv("F4 smoking associated genes_lm_combine FS and NS.csv")

association.S4 = merge(association.S4, annotation.S4F4[, c(1, 3)], by.x = "X", by.y = "Probe_Id")
association.F4 = merge(association.F4, annotation.S4F4[, c(1, 3)], by.x = "X", by.y = "Probe_Id")
write.csv(association.S4, file = "S4 smoking associated genes_lm.csv")
write.csv(association.F4, file = "F4 smoking associated genes_lm.csv") 


annotation.F3 = read.csv("expression/annotation_Illumina_WG_V2_KORA_F3.csv", sep = ";")
association.F3 = read.csv("F3 smoking associated genes_lm_combine FS and NS.csv")
association.F3 = merge(association.F3, annotation.F3[, c(1, 2)], by.x = "X", by.y = "Probe_Id")

write.csv(association.F3, file = "F3 smoking associated genes_lm.csv") 


genes.F3 = association.F3$Illumina_Gene[which(p.adjust(association.F3[,5],method = "BH")<0.05)]
genes.F4 = association.F4$ILMN_Gene[which(p.adjust(association.F4[,5], method="BH")<0.05)]
genes.S4 = association.S4$ILMN_Gene[which(p.adjust(association.S4[,5], method="BH")<0.05)]
genes.common = Reduce(intersect, list(genes.F3, genes.S4, genes.F4))
intersect(genes.F4, genes.S4)
