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


### supplemental table 1: list of all male and female DEGs
#get list of DEGs in male or female
allDEGs <- c(F5cHet_FWT.DEG$gene, M5cKO_MWT.DEG$gene)
#no repeats
allDEGs <-unique(allDEGs)

#make dataframe with everything together
female_padjbm <- subset(F5cHet_FWT.tb, F5cHet_FWT.tb$gene %in% allDEGs)
female_padjbm <- subset(female_padjbm, select = c("gene", "baseMean", "log2FoldChange", "padj"))
colnames(female_padjbm) <- c("gene", "baseMean", "Female_L2FC", "Female_padj")

male_padjbm <- subset(M5cKO_MWT.tb, M5cKO_MWT.tb$gene %in% allDEGs)
male_padjbm <- subset(male_padjbm, select = c("gene", "log2FoldChange", "padj"))
colnames(male_padjbm) <- c("gene", "Male_L2FC", "Male_padj")

MvsF_DEGs <- merge(male_padjbm, female_padjbm, by = "gene")
nrow(MvsF_DEGs)

#get the gene symbols
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
MvsF_DEGs_mgi <- bitr(MvsF_DEGs$gene, fromType = "ENSEMBL",
                          toType = "SYMBOL",
                          OrgDb = org.Mm.eg.db, drop = FALSE)


colnames(MvsF_DEGs_mgi) <- c("gene", "SYMBOL")
comb = merge(MvsF_DEGs_mgi, MvsF_DEGs, by = "gene" )

#determine which sex the gene is a DEG
maleDEG <- comb$Male_padj < 0.1
femaleDEG <- comb$Female_padj < 0.1
comb$DEG_sex <- ifelse(maleDEG & femaleDEG, "Both", ifelse(maleDEG & !femaleDEG, "Male Only", "Female Only"))

MvsF_DEGs <- subset(comb, select = c("gene", "SYMBOL", "DEG_sex", "baseMean", "Male_L2FC", "Male_padj", "Female_L2FC", "Female_padj"))
colnames(MvsF_DEGs) <- c("Ensembl", "MGI Symbol", "DEG_sex", "baseMean",  "Male_L2FC", "Male_padj", "Female_L2FC", "Female_padj")

write.table(MvsF_DEGs, file = "SupTable1_MalevsFemaleDEGs.csv", sep = ",", row.names = FALSE, col.names = TRUE)
