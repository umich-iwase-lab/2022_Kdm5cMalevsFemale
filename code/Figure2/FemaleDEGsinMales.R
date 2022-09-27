# expression of female DEGs in males and females 22.03.14

#male DEGs M5cKO_MWT.DEG

#all female DEGs F5cHet_FWT.DEG
p6femaleonly = subset(F5cHet_FWT.DEG, !(gene %in% M5cKO_MWT.DEG$gene))

fDEGinM <- subset(M5cKO_MWT.tb, gene %in% p6femaleonly$gene)
fDEGinM_alt <- subset(fDEGinM, select = c(gene, log2FoldChange, padj))
colnames(fDEGinM_alt) <- c("gene", "XY_log2FC","XY_padj")


p6femaleonly_alt <- subset(p6femaleonly, select = c(gene, baseMean, log2FoldChange, padj))
colnames(p6femaleonly_alt) <- c("gene","baseMean", "XX_log2FC","XX_padj")

fDEGinall <- merge(p6femaleonly_alt, fDEGinM_alt, by = "gene")
write.table(fDEGinall, file = "fDEGinall.csv", sep = ",")