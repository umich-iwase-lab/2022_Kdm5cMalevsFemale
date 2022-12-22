### 21.09.29 Analysis of male and female WT, Kdm5c-KO, and Kdm5c-Het P6 hippocampus and cortex 
#Load the libraries
library(VennDiagram)
#install.packages("DESeq2", dependencies=TRUE)
library(magrittr)
library(tibble)
library(tidyverse)
library(plyr)
library(dplyr)
#BiocManager::install("bit64")
library("DESeq2")
#library(tximport)
#library("ggpubr")
library("readxl")
library(ggplot2)
library("pheatmap")
library(VennDiagram)
library("EnhancedVolcano")


library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)


##################################################################
#                       Setting Up DESeq2                        #
##################################################################

#read in the data
#before reading in the data make sure the counts files are all integers and there are no isoform gene names
Kdm5c_P6_hipctx <- read.csv("C:/Users/kmbon/Documents/Projects/21.09.28_P6MaleandFemaleBrains/ms_P6brains_readcount_210929_FOR_DESEQ_withoutliers.txt", header=TRUE,sep="\t", row.names=1)
head(Kdm5c_P6_hipctx)

#sample information, made in excel as a separate sheet
SampleInfo <- read.csv("SampleInfo.csv", header=TRUE, sep=",", row.names=NULL)
SampleInfo <- data.frame(SampleInfo)
SampleInfo

coldata = data.frame(row.names = SampleInfo$�..Sample, sex = SampleInfo$Sex, genotype = SampleInfo$Genotype, type = c(rep("paired.end",16)))

#Which samples you want to compare, in this case we are comparing everything 
SAMPLES = coldata$type == "paired.end"

#keep the counts columns that you want to compare
countTable = Kdm5c_P6_hipctx[ , SAMPLES]

#Assign the genotype to each column 
genotype = coldata$genotype[ SAMPLES]
#check what genotypes you have
genotype

#make a DESeqDataSet (dds) from the counts table, starting with comparing just genotype (design)
dds <- DESeqDataSetFromMatrix(countData = Kdm5c_P6_hipctx,
                              colData = coldata,
                              design= ~ genotype)

#set WT as the default level for comparison 
dds$genotype <- relevel(dds$genotype, ref = "WT")

#add in sex so that we can account for the interaction between genotype and sex 
#essentially makes a new category (called "group") that has the conditions of sex and genotype in it
dds$group <- factor(paste0(dds$sex, dds$genotype))
#new design using the group variable
design(dds) <- ~ group

#run deseq on the DESeqDataSet to generate all the results tables
dds <- DESeq(dds)

#get the names of results, this is useful to seeing how their names are formatted
resultsNames(dds)
#pick which you want to compare
#in this case, we would want to compare the category "group" and then "F5cHet vs "FWT or M5cKO vs MWT


##################################################################
#                        Summary Analyses                        #
##################################################################
#these analyses are to see how the data looks overall (which samples are similar vs different)


ntd <- normTransform(dds)
#variance stabilizing transformation
#this is essentially accounts for sampling variablity of low counts while also normalizing variability across all samples 
#thus if a sample had all of its counts decreased, all of those genes wouldn't be automatically called as significant 
vsd <- vst(dds, blind=TRUE)
#rlog shriking then run into MA plot (vst and rlog)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

#simple PCA plot
#can set the intgroup equal to any category (genotype, group, or sex)
plotPCA(vsd, intgroup="group")


#change the color and symbol numbers to make the categories stand out

pdf(file = "PCA_P6HIPCTX_group.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches

data <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#adding sex as a separate variable back in so that the shape of the point is dependent upon sex
data <- data.frame(data, Sex = SampleInfo$Sex)

ggplot(data, aes(PC1, PC2, color=group, shape = Sex)) +
  geom_point(size=4.5) +
  ggtitle("PCA Plot of Postnatal Day 6 Forebrain")+
  scale_color_manual(values=c("#E36752", "#732113", "#556FAF","#2D3C61")) +
  scale_shape_manual(values=c(15, 17)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

dev.off()



pdf(file = "PCA_P6HIPCTX_maleonly.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches

data <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#adding sex as a separate variable back in so that the shape of the point is dependent upon sex
data2 <- subset(data, Sex == "M")


ggplot(data2, aes(PC1, PC2, color=group, shape = Sex)) +
  geom_point(size=4.5) +
  ggtitle("PCA Plot of Postnatal Day 6 Forebrain")+
  scale_color_manual(values=c("#E36752", "#732113", "#556FAF","#2D3C61")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

dev.off()




######another way of representing similarities between samples - how do they cluster with each other?
####how distant each sample is from each other
sampleDists <- dist(t(assay(vsd)))
pdf(file = "distmatrix_P6HIPCTX_group.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


######another way of representing similarities between samples - how do they cluster with each other?
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("genotype","sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

##################################################################
#                            MA Plots                            #
##################################################################

alph = 0.1
malim = c(-4,4)


########## Female Het vs Female WT #######################
res_F5cHet_FWT_unsh <- results(dds, contrast=c("group", "F5cHet", "FWT"), alpha = alph)
plotMA(res_F5cHet_FWT_unsh, xlim=c(0.1,100000), ylim=malim, xlab="Mean of Normalized Counts",ylab="log2 Fold Change",main=paste("MA Plot, XX WT vs. XX Kdm5c-KO,", "a =", alph))

#used the default shrinkage
F5cHet_FWT <- lfcShrink(dds, contrast=c("group", "F5cHet", "FWT"), res=res_F5cHet_FWT_unsh)
#make an MA plot of the results (unshrunken)
pdf(file = "MA_F5cHet_FWT_a0.1.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches
plotMA(F5cHet_FWT, xlim=c(0.1,100000), ylim=malim, xlab="Mean of Normalized Counts",ylab="log2 Fold Change",main=paste("MA Plot, XX WT vs. XX Kdm5c-KO,", "a =", alph))
dev.off()
write.table(F5cHet_FWT, file = "P6Forebrain_res_F5cHet_FWT_0.1.txt", sep = "\t")

summary(F5cHet_FWT)


########## Male KO vs WT #######################
res_M5cKO_MWT_unsh <- results(dds, contrast=c("group", "M5cKO", "MWT"), alpha = alph)
#used the default shrinkage
M5cKO_MWT <- lfcShrink(dds, contrast=c("group", "M5cKO", "MWT"), res=res_M5cKO_MWT_unsh)
#make an MA plot of the results (unshrunken)
pdf(file = "MA_M5cKO_MWT_a0.1.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches
plotMA(M5cKO_MWT, xlim=c(0.1,100000), ylim=malim, xlab="Mean of Normalized Counts",ylab="log2 Fold Change",main=paste("MA Plot, XY WT vs. XY Kdm5c-KO,", "a =", alph))
dev.off()
write.table(M5cKO_MWT, file = "P6Forebrain_res_M5cKO_MWT_0.1.txt", sep = "\t")

summary(M5cKO_MWT)


##################################################################
#                           Biomart                              #
##################################################################
#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html

#genelist would be a dataframe where the first column is the list of ensembl
#this returns a list of genes with the mgi symbol 
ensembl2mgi <- function(genelist){
  library("biomaRt")
  ensembl <- useMart("ensembl")
  mart <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
  genie <- getBM(
    filters="ensembl_gene_id",
    attributes=c("ensembl_gene_id", "mgi_symbol"),
    values=genelist[,1],
    mart=mart)
  colnames(genie) <- c("gene","mgi")
  changedlist = merge(genelist, genie, by="gene")
  
  return(changedlist)
}

library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)


gene.df <- bitr(gene, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)


##################################################################
#                          Sort out DEGs                         #
##################################################################

#change into a tibble
F5cHet_FWT.tb <- F5cHet_FWT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
#get the DEGs
F5cHet_FWT.DEG = subset(F5cHet_FWT.tb, padj < alph)
nrow(F5cHet_FWT.DEG)
#upregulated DEGs
F5cHet_FWT.UPDEG = subset(F5cHet_FWT.DEG, log2FoldChange > 0)
#downregulated DEGs
F5cHet_FWT.DOWNDEG = subset(F5cHet_FWT.DEG, log2FoldChange < 0)


#convert into MGI symbols
F5cHet_FWT.DEG_mgi = ensembl2mgi(F5cHet_FWT.DEG)
#write the DEGs as a text file for opening with Excel
write.table(F5cHet_FWT.DEG, file="F5cHet_FWT.DEG_a0.1.txt", sep="\t")
write.table(F5cHet_FWT.DEG_mgi, file="F5cHet_FWT.DEG_a0.1_mgi.txt", sep="\t")


#change into tibble
M5cKO_MWT.tb <- M5cKO_MWT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
#get the DEGs
M5cKO_MWT.DEG = subset(M5cKO_MWT.tb, padj < alph)
#up DEGs
M5cKO_MWT.UPDEG = subset(M5cKO_MWT.DEG, log2FoldChange > 0)
#down DEGs
M5cKO_MWT.DOWNDEG = subset(M5cKO_MWT.DEG, log2FoldChange < 0)



M5cKO_MWT.DEG_mgi <- bitr(M5cKO_MWT.DEG$gene, fromType = "ENSEMBL",
                         toType = c("ENTREZID", "SYMBOL"),
                         OrgDb = org.Mm.eg.db)


colnames(M5cKO_MWT.DEG_mgi) <- c("gene", "ENTREZID", "SYMBOL")
comb = merge(M5cKO_MWT.DEG_mgi, M5cKO_MWT.DEG, by = "gene" )


write.table(M5cKO_MWT.DEG, file="M5cKO_MWT.DEG_a0.1.txt", sep="\t")
write.table(comb, file="M5cKO_MWT.DEG_a0.1_mgi.txt", sep="\t")

write.table(M5cKO_MWT.tb, file="M5cKO_MWT.allgenes.a0.1.txt", sep="\t")


# #Number of DEGs in each category
# sex.chromo <- c("XY", "XX", "XX")
# geno.5c <- c("Kdm5c-KO", "Kdm5c-Het", "Kdm5c-KO")
# DEG.count <- c(nrow(XY_KO_WT.DEG), nrow(XX_Het_WT.DEG), nrow(XX_KO_WT.DEG))
# 
# DEG.counts.all <- data.frame(sex.chromo, geno.5c, DEG.count)
# DEG.counts.all


##################################################################
#                          Shared DEGs                           #
##################################################################

############## shared between male and female p6 ##########################
#upregulated genes
mandfDEGs = subset(M5cKO_MWT.DEG, gene %in% F5cHet_FWT.DEG$gene)
nrow(mandfDEGs)

## all of the shared DEGs with the log2fc of male and female
mandfDEGs <- subset(mandfDEGs, select = c(gene, log2FoldChange, padj))
colnames(mandfDEGs) <- c("gene","Male_L2FC","Male_padj")

F_short <- subset(F5cHet_FWT.DEG, select = c(gene, log2FoldChange, padj))
colnames(F_short) <- c("gene","Female_L2FC","Female_padj")

mandfDEGs <- merge(mandfDEGs, F_short, by = "gene" )
write.table(mandfDEGs, file = "mandfsharedDEGs.txt", sep ="\t")


mandfDEGs_up_mgi <- bitr(mandfDEGs_up$gene, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)


### downregulated genes
mandfDEGs_down = subset(M5cKO_MWT.DOWNDEG, gene %in% F5cHet_FWT.DOWNDEG$gene)
nrow(mandfDEGs_down)


#only in female
p6femaleonly_up = subset(F5cHet_FWT.UPDEG, !(gene %in% M5cKO_MWT.UPDEG$gene))
p6maleonly_up = subset(M5cKO_MWT.UPDEG, !(gene %in% F5cHet_FWT.UPDEG$gene))


#all female only DEGs (up or down)
p6femaleonly = subset(F5cHet_FWT.DEG, !(gene %in% M5cKO_MWT.DEG$gene))
write.table(p6femaleonly, file="femaleonly_DEG_a0.1.txt", sep="\t")

#all male only DEGs (up or down)
p6maleonly = subset(M5cKO_MWT.DEG, !(gene %in% F5cHet_FWT.DEG$gene))
write.table(p6maleonly, file="maleonly_DEG_a0.1.txt", sep="\t")


########################################################################################
################ expression of all up DEGs in males vs females #########################
################# to show sex-specific uniqueness and het loss #########################
########################################################################################

##to see differences in dysregulation
#make a dataframe with all the germline DEGs and their l2fc in the male and female p6 brain

#Male up genes #M5cKO_MWT.UPDEG$gene
#male up degs with log2fc
M5cKO_MWT.UPDEG_l2fc <- data.frame(gene = M5cKO_MWT.UPDEG$gene, log2FoldChange= M5cKO_MWT.UPDEG$log2FoldChange)
  #all male genes #M5cKO_MWT.tb

#Femle up genes #F5cHet_FWT.UPDEG$gene
#female up degs with log2fc
F5cHet_FWT.UPDEG_l2fc <- data.frame(gene = F5cHet_FWT.UPDEG$gene, log2FoldChange= F5cHet_FWT.UPDEG$log2FoldChange)
  #all female genes #F5cHet_FWT.tb


#list of male and female up DEGs with their log2fc

#brain and liver testis DEGs
p6malefemale.upgene <- c(M5cKO_MWT.UPDEG$gene, F5cHet_FWT.UPDEG$gene)
length(p6malefemale.upgene)

#label which sex they're DEGs in 
sex <- c()
count = 1
for (t in p6malefemale.upgene){
  
  if (t %in% M5cKO_MWT.UPDEG$gene & t %in% F5cHet_FWT.UPDEG$gene){
    a = "Both"
  } else if (t %in% M5cKO_MWT.UPDEG$gene) {
    a = "Male Up DEG"
  } else {
    a = "Female Up DEG"
  }
  
  sex[count] <- a
  count = count+1
}


p6malefemale.up.df <- data.frame(gene = p6malefemale.upgene, sex = sex) 
#remove duplicates
p6malefemale.up.df <- unique(p6malefemale.up.df)
nrow(p6malefemale.up.df)
head(p6malefemale.up.df)
write.table(p6malefemale.up.df, file = "Male_Female_UPdegs_L2fc.txt", sep = "\t")

#add the log2fold change
  #the log2fc for the Female samples
Female_L2FC <- subset(F5cHet_FWT.tb, F5cHet_FWT.tb$gene %in% p6malefemale.up.df$gene)
nrow(Female_L2FC)
Female_L2FC <- data.frame(gene = Female_L2FC$gene, Female_L2FC = Female_L2FC$log2FoldChange)

  #the log2fc for the Male samples
Male_L2FC <- subset(M5cKO_MWT.tb, M5cKO_MWT.tb$gene %in% p6malefemale.up.df$gene)
nrow(Male_L2FC) #lost one gene
Male_L2FC <- data.frame(gene = Male_L2FC$gene, Male_L2FC = Male_L2FC$log2FoldChange)


#merge males and females into one dataframe with labels
p6malefemale.up_l2fc_symb <- read.csv("Male_Female_UPdegs_L2fc_sym.txt", sep='\t', header=T)

p6malefemale.up_l2fc <- merge(Male_L2FC, Female_L2FC, by = "gene")
p6malefemale.up_l2fc_symb <- merge(p6malefemale.up_l2fc, p6malefemale.up_l2fc_symb, by = "gene")

#write this dataframe as a table
write.table(p6malefemale.up_l2fc_symb, )

pdf(file = "P6MalevsFemale_Brain_L2FC_scatter_empty.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 3.5) # The height of the plot in inches

ggplot(p6malefemale.up_l2fc_symb, aes(x=Male_L2FC, y=Female_L2FC, shape=sex, color =sex)) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed") +
  geom_point(size=3, alpha=0.7) +
  xlim(-0.5,2) +
  ylim(-0.5,2) +
  scale_color_manual(values=c("#CC8C4C", "#E36752", "#556FAF"))+
  scale_shape_manual(values=c(1, 0, 2))

dev.off()


p6malefemale.up_l2fc_symb_shared <- subset(p6malefemale.up_l2fc_symb, sex == "Both")

##calculate the correlation
library(tidyverse)
library(ggpubr)


model <- lm(Female_L2FC ~ Male_L2FC, data = p6malefemale.up_l2fc_symb_shared)
model
summary(model)
rsq <- function (x, y) cor(x, y) ^ 2
rsq(p6malefemale.up_l2fc_symb_shared$Male_L2FC, p6malefemale.up_l2fc_symb_shared$Female_L2FC)


pdf(file = "P6MalevsFemale_Brain_L2FC_shared_scatter_trend.pdf",   # The directory you want to save the file in
    width = 4.5, # The width of the plot in inches
    height = 3.5) # The height of the plot in inches

ggplot(p6malefemale.up_l2fc_symb_shared, aes(x=Male_L2FC, y=Female_L2FC, shape=sex, color =sex)) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed") +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=c("#CC8C4C"))+
  scale_shape_manual(values=c(19)) +
  geom_smooth(method=lm, se=FALSE, col='#91500f', size=1) +
  xlim(-0.5,2) +
  ylim(-0.5,2)

dev.off()



##############################################################
# Write a table with all the male vs female upregulated DEGs #

#add the adjusted p-values and basemean
female_padjbm <- subset(F5cHet_FWT.tb, F5cHet_FWT.tb$gene %in% p6malefemale.up_l2fc$gene)
female_padjbm <- subset(female_padjbm, select = c("gene", "baseMean", "padj"))
colnames(female_padjbm) <- c("gene", "baseMean", "Female_padj")

male_padjbm <- subset(M5cKO_MWT.tb, M5cKO_MWT.tb$gene %in% p6malefemale.up_l2fc$gene)
male_padjbm <- subset(male_padjbm, select = c("gene", "padj"))
colnames(male_padjbm) <- c("gene", "Male_padj")

#merge to final df
MvsF_DEGs <- merge(p6malefemale.up_l2fc_symb, female_padjbm, by = "gene")
MvsF_DEGs <- merge(MvsF_DEGs, male_padjbm, by = "gene")
MvsF_DEGs <- subset(MvsF_DEGs, select = c("gene", "MGI", "baseMean", "sex", "Male_L2FC", "Male_padj", "Female_L2FC", "Female_padj"))
colnames(MvsF_DEGs) <- c("Ensembl", "Symbol", "baseMean", "DEG_sex", "Male_L2FC", "Male_padj", "Female_L2FC", "Female_padj")

write.table(MvsF_DEGs, file = "../../results/SupTable1_MalevsFemaleUPDEGs.csv", sep = ",", row.names = FALSE, col.names = TRUE)






##################################################################################
################ male p6 shared with adult male brain ############################
#all of the hippocampus genes in the results table

#results of deseq2 in hippocampus
wt_5c_hip <- read.csv("C:/Users/kmbon/Documents/Projects/KDM5C_DM_RNAseq/Analysis/wt_5c_hip_a0.1.txt", header=TRUE,sep="\t", row.names=1)
head(wt_5c_hip)


wt_5c_hip.tb <- wt_5c_hip %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

HIP_KO.DEG = subset(wt_5c_hip.tb, padj < alph)
HIP_KO.UPDEG = subset(HIP_KO.DEG, log2FoldChange > 0)


adultHIP_mp6_shared = subset(M5cKO_MWT.UPDEG, gene %in% HIP_KO.UPDEG$gene)
nrow(adultHIP_mp6_shared)

#unique to p6
nrow(M5cKO_MWT.UPDEG) - nrow(adultHIP_mp6_shared)

#unique to adult
nrow(HIP_KO.UPDEG) - nrow(adultHIP_mp6_shared)


#get the mgi names for shared ones
adultHIP_mp6_shared_mgi <- bitr(adultHIP_mp6_shared$gene, fromType = "ENSEMBL",
                         toType = c("ENTREZID", "SYMBOL"),
                         OrgDb = org.Mm.eg.db)



########### for the adult amygdala
#results of deseq2 in amygdala
wt_5c_amy <- read.csv("C:/Users/kmbon/Documents/Projects/KDM5C_DM_RNAseq/Analysis/wt_5c_amy_a0.1.txt", header=TRUE,sep="\t", row.names=1)
head(wt_5c_amy)


wt_5c_amy.tb <- wt_5c_amy %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

AMY_KO.DEG = subset(wt_5c_amy.tb, padj < alph)
AMY_KO.UPDEG = subset(AMY_KO.DEG, log2FoldChange > 0)


adultAMY_mp6_shared = subset(M5cKO_MWT.UPDEG, gene %in% AMY_KO.UPDEG$gene)
nrow(adultAMY_mp6_shared)

#unique to p6
nrow(M5cKO_MWT.UPDEG) - nrow(adultAMY_mp6_shared)

#unique to adult
nrow(AMY_KO.UPDEG) - nrow(adultAMY_mp6_shared)



##### comparisons between all 3 (adult hippocampus, adult amydgala, and p6 forebrain)

#shared between all
alladult_mp6_shared = subset(M5cKO_MWT.UPDEG, gene %in% AMY_KO.UPDEG$gene & gene %in% HIP_KO.UPDEG$gene)
nrow(alladult_mp6_shared)

#shared between p6 forebrain and amygdala only
adultAMY_mp6_ONLY = subset(M5cKO_MWT.UPDEG, gene %in% AMY_KO.UPDEG$gene & !(gene %in% HIP_KO.UPDEG$gene))
nrow(adultAMY_mp6_ONLY)

#shared between p6 forebrain and hippocampus only
adultHIP_mp6_ONLY = subset(M5cKO_MWT.UPDEG, gene %in% HIP_KO.UPDEG$gene & !(gene %in% AMY_KO.UPDEG$gene))
nrow(adultHIP_mp6_ONLY)


#shared between amygdala and hippocampus only
adultHIP_AMY_ONLY = subset(AMY_KO.UPDEG, gene %in% HIP_KO.UPDEG$gene & !(gene %in% M5cKO_MWT.UPDEG$gene))
nrow(adultHIP_AMY_ONLY)


#only in hippocampus
adultHIP_ONLY = subset(HIP_KO.UPDEG, !(gene %in% AMY_KO.UPDEG$gene) & !(gene %in% M5cKO_MWT.UPDEG$gene))
nrow(adultHIP_ONLY)

#only in amygdala
adultAMY_ONLY = subset(AMY_KO.UPDEG, !(gene %in% HIP_KO.UPDEG$gene) & !(gene %in% M5cKO_MWT.UPDEG$gene))
nrow(adultAMY_ONLY)


#only in p6 forebrain
p6_ONLY = subset(M5cKO_MWT.UPDEG, !(gene %in% HIP_KO.UPDEG$gene) & !(gene %in% AMY_KO.UPDEG$gene))
nrow(p6_ONLY)











##################################################################
#                     Testis Specific DEGs                       #
#                           Page Lab                             #
##################################################################

#list of testis genes based on page paper
pagetest <- read.csv("C:/Users/kmbon/Documents/Projects/KDM5C_GermlineGenes/KDM5C_EpiLC_RNAseq/PageTestisGenesAnalysis/page2020_testisgenes_ensembl.txt", header=TRUE,sep="\t")
head(pagetest) 


#how many up DEGs are testis genes? Males
pagetest.5cKO = subset(M5cKO_MWT.UPDEG, gene %in% pagetest$gene)

pagetest.5cKO = merge(pagetest, pagetest.5cKO, by = "gene")
nrow(pagetest.5cKO)

write.table(pagetest.5cKO, file = "pagetest.5cKOP6BrainDEG.txt", sep ="\t")



#how many up DEGs are testis genes? Females
pagetest.5cHet = subset(F5cHet_FWT.UPDEG, gene %in% pagetest$gene)

pagetest.5cHet = merge(pagetest, pagetest.5cHet, by = "gene")
nrow(pagetest.5cHet)
write.table(pagetest.5cHet, file = "pagetest.5cHetP6BrainDEG.txt", sep ="\t")




############################# ovary genes ################################
#import testis-specific gene data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5482823/

library("readxl")
tissue_genes <- read_excel("C:/Users/kmbon/Documents/Projects/KDM5C_DM_RNAseq/Analysis/DM_Plots/Testis_Plots/Tissue_specific.xlsx")
Li_ovary_genes = subset(tissue_genes, tissue=="Ov")
#write.table(testis_genes, file="testis_genes.txt", sep="\t")


#look at testis genes in Kdm5c KO hippocampus
LiOvary.5cKO = subset(M5cKO_MWT.UPDEG, gene %in% Li_ovary_genes$gene_id)
LiOvary.5cKO_sym <- bitr(LiOvary.5cKO$gene, fromType = "ENSEMBL",
                          toType = c("SYMBOL"),
                          OrgDb = org.Mm.eg.db)
colnames(LiOvary.5cKO_sym) <- c("gene", "symbol")
LiOvary.5cKO = merge(LiOvary.5cKO_sym, LiOvary.5cKO, by = "gene")
nrow(LiOvary.5cKO)




LiOvary.5cHet = subset(F5cHet_FWT.UPDEG, gene %in% Li_ovary_genes$gene_id)

LiOvary.5cHet_sym <- bitr(LiOvary.5cHet$gene, fromType = "ENSEMBL",
                          toType = c("SYMBOL"),
                          OrgDb = org.Mm.eg.db)
colnames(LiOvary.5cHet_sym) <- c("gene", "symbol")
LiOvary.5cHet = merge(LiOvary.5cHet_sym, LiOvary.5cHet, by = "gene")
nrow(LiOvary.5cHet)





LiOvary.AMY = subset(AMY_KO.UPDEG, gene %in% Li_ovary_genes$gene_id)

LiOvary.AMY_sym <- bitr(LiOvary.5cHet$gene, fromType = "ENSEMBL",
                          toType = c("SYMBOL"),
                          OrgDb = org.Mm.eg.db)
colnames(LiOvary.AMY_sym) <- c("gene", "symbol")
LiOvary.AMY = merge(LiOvary.AMY_sym, LiOvary.AMY, by = "gene")
nrow(LiOvary.AMY)


LiOvary.HIP = subset(HIP_KO.UPDEG, gene %in% Li_ovary_genes$gene_id)

LiOvary.HIP_sym <- bitr(LiOvary.5cHet$gene, fromType = "ENSEMBL",
                        toType = c("SYMBOL"),
                        OrgDb = org.Mm.eg.db)
colnames(LiOvary.HIP_sym) <- c("gene", "symbol")
LiOvary.HIP = merge(LiOvary.HIP_sym, LiOvary.HIP, by = "gene")
nrow(LiOvary.HIP)

############ pie charts
library(ggplot2)
library(scales)

# write.table(pagetest.5cKOXY, file="C:/Users/kmbon/Documents/Projects/KDM5C_GermlineGenes/KDM5C_EpiLC_RNAseq/PageTestisGenesAnalysis/XY_5cKO_pagetestDEGS.txt", sep="\t")

#### 5c Ko Male
p65cKO_page_pie <- data.frame(
  group=c("TestisDEG","Non-TestisDEG"),
  value=c(46,662)
)

# Basic piechart
p.p65cKO_page_pie <- ggplot(p65cKO_page_pie, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=3, color="black") +
  coord_polar("y", start=0) +
  
  theme_void()+ # remove background, grid, numeric labels
  scale_fill_manual(name = " ", labels = c("Non-Testis UP DEG","Testis DEG"), values=c("#FFBD04", "#6BA96B"))+
  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value/sum(p65cKO_page_pie[,2]))), size=5)

p.p65cKO_page_pie

pdf(file = "PieChart_p65cKO_page.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 3) # The height of the plot in inches
p.p65cKO_page_pie
dev.off()

#statistics
P65cKO.page.dat <- data.frame(
  "gcg_no" = c(662, 51787),
  "gcg_yes" = c(46, 1724),
  row.names = c("Up DEG", "Not Up DEG"),
  stringsAsFactors = FALSE
)
colnames(P65cKO.page.dat) <- c("Non-germline", "Germline")
P65cKO.page.dat

chisq.test(P65cKO.page.dat)$expected
P65cKO.page.test <- fisher.test(P65cKO.page.dat)
P65cKO.page.test




#### 5cHet female
p65cHet_page_pie <- data.frame(
  group=c("TestisDEG","Non-TestisDEG"),
  value=c(10,112)
)

# Basic piechart
p.p65cHet_page_pie <- ggplot(p65cHet_page_pie, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=3, color="black") +
  coord_polar("y", start=0) +
  
  theme_void()+ # remove background, grid, numeric labels
  scale_fill_manual(name = " ", labels = c("Non-Testis UP DEG","Testis DEG"), values=c("#FFBD04", "#6BA96B"))+
  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value/sum(p65cHet_page_pie[,2]))), size=5)

p.p65cHet_page_pie

pdf(file = "PieChart_p65cHet_page.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 3) # The height of the plot in inches
p.p65cHet_page_pie
dev.off()

#statistics
P65cHet.page.dat <- data.frame(
  "gcg_no" = c(112, 52337),
  "gcg_yes" = c(10, 1760),
  row.names = c("Up DEG", "Not Up DEG"),
  stringsAsFactors = FALSE
)
colnames(P65cHet.page.dat) <- c("Non-germline", "Germline")
P65cHet.page.dat

chisq.test(P65cHet.page.dat)$expected
P65cHet.page.test <- fisher.test(P65cHet.page.dat)
P65cHet.page.test






#21.10.30
##################################################################
#                       WWv Germline Genes                       #
##################################################################

allgermgenes_WWv <- read.csv("C:/Users/kmbon/Documents/Projects/KDM5C_GermlineGenes/21.10.28_WhatisaGermCellGene/allgermgenes_WWv.txt", header=TRUE,sep="\t")
head(allgermgenes_WWv)


#### expression of germline genes in P6 male brain
#P6 Male DEGs that are germline genes
M5cKO_MWT.UPDEG_WWvgerm <- subset(allgermgenes_WWv, allgermgenes_WWv$ENSEMBL %in% M5cKO_MWT.UPDEG$gene)


#### expression of germline genes in P6 female brain
#P6 female DEGs that are germline genes
F5cHet_FWT.UPDEG_WWvgerm <- subset(allgermgenes_WWv, allgermgenes_WWv$ENSEMBL %in% F5cHet_FWT.UPDEG$gene)


## adult male brains wwv germline genes
AMY_KO.UPDEG_WWvgerm <- subset(allgermgenes_WWv, allgermgenes_WWv$ENSEMBL %in% AMY_KO.UPDEG$gene)
HIP_KO.UPDEG_WWvgerm <- subset(allgermgenes_WWv, allgermgenes_WWv$ENSEMBL %in% HIP_KO.UPDEG$gene)


##################################################################
#                     Li 2017 Testis or Ovary Genes              #
##################################################################

Li_TeOvBr_1 <- read.csv("C:/Users/kmbon/Documents/Projects/KDM5C_GermlineGenes/21.10.28_WhatisaGermCellGene/Li_TeOvBr_1.txt", header=TRUE,sep="\t")
head(Li_TeOvBr_1)

Li_TeOvBr_Test <- subset(Li_TeOvBr_1, Li_TeOvBr_1$Tissue == "Testis")
Li_TeOvBr_TestandOvary <- subset(Li_TeOvBr_1, Li_TeOvBr_1$Tissue == "Testis and Ovary")
Li_TeOvBr_Ovary <- subset(Li_TeOvBr_1, Li_TeOvBr_1$Tissue == "Ovary")
Li_TeOvBr_Brain<- subset(Li_TeOvBr_1, Li_TeOvBr_1$Tissue == "Brain")


##################################################################
#                 Volcano with opposite sex DEGs                 #
##################################################################

library('DESeq2')
library("EnhancedVolcano")


#all of the p6 male genes in the results table
head(M5cKO_MWT)


#genes that are DEGs in females but not males
head(p6femaleonly_up)
p6femaleonly_down = subset(F5cHet_FWT.DOWNDEG, !(gene %in% M5cKO_MWT.DOWNDEG$gene))
head(p6femaleonly_down)

p6femaleonly_all <- rbind(p6femaleonly_up,p6femaleonly_down)

#genes that are DEGs in males but not females
head(p6maleonly_up)
p6maleonly_down = subset(M5cKO_MWT.DOWNDEG, !(gene %in% F5cHet_FWT.DOWNDEG$gene))
head(p6maleonly_down)

p6maleonly_all <- rbind(p6maleonly_up,p6maleonly_down)

#genes that are DEGs in both
mandfDEGs_up = subset(M5cKO_MWT.UPDEG, gene %in% F5cHet_FWT.UPDEG$gene)
nrow(mandfDEGs_up)

mandfDEGs_down = subset(M5cKO_MWT.DOWNDEG, gene %in% F5cHet_FWT.DOWNDEG$gene)
nrow(mandfDEGs_down)

mandfDEGs_all <- rbind(mandfDEGs_up, mandfDEGs_down)


#genes regulated in opposite directions (none)
mandfDEGs_updown = subset(M5cKO_MWT.UPDEG, gene %in% F5cHet_FWT.DOWNDEG$gene)
nrow(mandfDEGs_updown)
mandfDEGs_downup = subset(M5cKO_MWT.DOWNDEG, gene %in% F5cHet_FWT.UPDEG$gene)
nrow(mandfDEGs_downup)



pcut = 1e-30

keyvals <- ifelse(
  row.names(M5cKO_MWT) %in% mandfDEGs_all$gene, '#f2b722ff',
  ifelse(row.names(M5cKO_MWT) %in% p6maleonly_all$gene, '#556FAF',
         ifelse(row.names(M5cKO_MWT) %in% p6femaleonly_all$gene, '#e36752ff',
                ifelse(M5cKO_MWT$log2FoldChange >= 0.1 & M5cKO_MWT$padj <= 0.1, '#F6B000',
                       ifelse(M5cKO_MWT$log2FoldChange <= -0.1 & M5cKO_MWT$padj <= 0.1, '#6991D1',
                              'gray')))))


keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#F6B000'] <- 'Up DEGs'
names(keyvals)[keyvals == '#6991D1'] <- 'Down DEGs'
names(keyvals)[keyvals == '#f2b722ff'] <- 'Shared DEGs'
names(keyvals)[keyvals == '#556FAF'] <- 'Male Only DEGs'
names(keyvals)[keyvals == '#e36752ff'] <- 'Female Only DEGs'
names(keyvals)[keyvals == 'gray'] <- 'NS'

#making shapes different outside the axis range
keyvals.shape <- ifelse(
  M5cKO_MWT$log2FoldChange > 1.5, 9,
  ifelse(M5cKO_MWT$padj < pcut, 9,
         1))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 9] <- 'Outside Axis'
names(keyvals.shape)[keyvals.shape == 1] <- 'Normal'

M5cKO_MWT_plot <- M5cKO_MWT
M5cKO_MWT_plot$log2FoldChange[M5cKO_MWT_plot$log2FoldChange>1.5] <- 1.5
M5cKO_MWT_plot$padj[M5cKO_MWT_plot$padj<pcut] <- pcut


# M5cKO_MWT.sort = M5cKO_MWT.UPDEG[order(M5cKO_MWT.UPDEG[,3]), ]
# top10_M5cKO <- tail(M5cKO_MWT.sort, 10)

VolcM5cKO_shared <- EnhancedVolcano(M5cKO_MWT_plot,
                                 lab = row.names(M5cKO_MWT_plot),
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 selectLab = NA,
                                 xlim=c(-0.5,1.5),
                                 ylim=c(0,30),
                                 xlab = bquote(~Log[2]~ 'fold change'),
                                 ylab = bquote(~-Log[10]~adjusted~italic(P)),
                                 pCutoff = 0.1,
                                 FCcutoff = 0,
                                 title = "Kdm5c-KO Male P6 Forebrain RNA-seq",
                                 labSize = 4.0,
                                 colAlpha = 5/5,
                                 colCustom = keyvals,
                                 pointSize = 4.0,
                                 shapeCustom = keyvals.shape,
                                 #col=c('lightsteelblue4', 'lightsteelblue4', 'lightsteelblue4', 'red3'),
                                 legend=c('NS','Log2 FC','Adjusted p-value',
                                          'Adjusted p-value & Log2 FC'),
                                 legendPosition = 'bottom',
                                 gridlines.major = FALSE,
                                 gridlines.minor = FALSE,
                                 legendLabSize = 10,
                                 legendIconSize = 3.0)

VolcM5cKO_shared




pdf(file = "Volcano_P6Male_SharedDEGs.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 8) # The height of the plot in inches
VolcM5cKO_shared
dev.off()


### female

#all of the p6 Female genes in the results table


keyvals <- ifelse(
  row.names(F5cHet_FWT) %in% mandfDEGs_all$gene, '#f2b722ff',
  ifelse(row.names(F5cHet_FWT) %in% p6maleonly_all$gene, '#556FAF',
         ifelse(row.names(F5cHet_FWT) %in% p6femaleonly_all$gene, '#e36752ff',
                ifelse(F5cHet_FWT$log2FoldChange >= 0.1 & F5cHet_FWT$padj <= 0.1, '#F6B000',
                       ifelse(F5cHet_FWT$log2FoldChange <= -0.1 & F5cHet_FWT$padj <= 0.1, '#6991D1',
                              'gray')))))


keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#F6B000'] <- 'Up DEGs'
names(keyvals)[keyvals == '#6991D1'] <- 'Down DEGs'
names(keyvals)[keyvals == '#f2b722ff'] <- 'Shared DEGs'
names(keyvals)[keyvals == '#556FAF'] <- 'Male Only DEGs'
names(keyvals)[keyvals == '#e36752ff'] <- 'Female Only DEGs'
names(keyvals)[keyvals == 'gray'] <- 'NS'

#making shapes different outside the axis range
keyvals.shape <- ifelse(
  F5cHet_FWT$log2FoldChange > 1.5, 9,
  ifelse(F5cHet_FWT$padj < pcut, 9,
         1))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 9] <- 'Outside Axis'
names(keyvals.shape)[keyvals.shape == 1] <- 'Normal'

F5cHet_FWT_plot <- F5cHet_FWT
F5cHet_FWT_plot$log2FoldChange[F5cHet_FWT_plot$log2FoldChange>1.5] <- 1.5
F5cHet_FWT_plot$padj[F5cHet_FWT_plot$padj<pcut] <- pcut


# F5cHet_FWT.sort = F5cHet_FWT.UPDEG[order(F5cHet_FWT.UPDEG[,3]), ]
# top10_F5cHet <- tail(F5cHet_FWT.sort, 10)

VolcF5cHet_shared <- EnhancedVolcano(F5cHet_FWT_plot,
                                    lab = row.names(F5cHet_FWT_plot),
                                    x = 'log2FoldChange',
                                    y = 'padj',
                                    selectLab = NA,
                                    xlim=c(-0.5,1.5),
                                    ylim=c(0,30),
                                    xlab = bquote(~Log[2]~ 'fold change'),
                                    ylab = bquote(~-Log[10]~adjusted~italic(P)),
                                    pCutoff = 0.1,
                                    FCcutoff = 0,
                                    title = "Kdm5c-Het Female P6 Forebrain RNA-seq",
                                    labSize = 4.0,
                                    colAlpha = 5/5,
                                    pointSize = 4.0,
                                    colCustom = keyvals,
                                    shapeCustom = keyvals.shape,
                                    #col=c('lightsteelblue4', 'lightsteelblue4', 'lightsteelblue4', 'red3'),
                                    legend=c('NS','Log2 FC','Adjusted p-value',
                                             'Adjusted p-value & Log2 FC'),
                                    legendPosition = 'bottom',
                                    gridlines.major = FALSE,
                                    gridlines.minor = FALSE,
                                    legendLabSize = 10,
                                    legendIconSize = 3.0)

VolcF5cHet_shared





pdf(file = "Volcano_P6Female_SharedDEGs.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 8) # The height of the plot in inches
VolcF5cHet_shared
dev.off()



##################################################################
#                      Volcano with Li Genes                     #
##################################################################

library('DESeq2')
library("EnhancedVolcano")


#all of the p6 male genes in the results table
head(M5cKO_MWT)

#get the mgi names
# M5cKO_MWT_mgi <- bitr(row.names(M5cKO_MWT), fromType = "ENSEMBL",
#                 toType = c("SYMBOL"),
#                 OrgDb = org.Mm.eg.db)





keyvals <- ifelse(
  row.names(M5cKO_MWT) %in% Li_TeOvBr_Test$gene_id, '#6BA96B',
  ifelse(row.names(M5cKO_MWT) %in% Li_TeOvBr_Ovary$gene_id, 'red',
         ifelse(row.names(M5cKO_MWT) %in% Li_TeOvBr_TestandOvary$gene_id, 'black',
                ifelse(row.names(M5cKO_MWT) %in% Li_TeOvBr_Brain$gene_id, 'blue',
                       ifelse(M5cKO_MWT$log2FoldChange >= 0.1 & M5cKO_MWT$padj <= 0.1, '#F6B000',
                              ifelse(M5cKO_MWT$log2FoldChange <= -0.1 & M5cKO_MWT$padj <= 0.1, '#6991D1',
                                     'gray'))))))


keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#F6B000'] <- 'Up DEGs'
names(keyvals)[keyvals == '#6991D1'] <- 'Down DEGs'
names(keyvals)[keyvals == '#6BA96B'] <- 'Testis DEGs'
names(keyvals)[keyvals == 'red'] <- 'Ovary DEGs'
names(keyvals)[keyvals == 'black'] <- 'Testis and Ovary DEGs'
names(keyvals)[keyvals == 'blue'] <- 'Brain DEGs'
names(keyvals)[keyvals == 'gray'] <- 'NS'





# M5cKO_MWT.sort = M5cKO_MWT.UPDEG[order(M5cKO_MWT.UPDEG[,3]), ]
# top10_M5cKO <- tail(M5cKO_MWT.sort, 10)

VolcM5cKO_all <- EnhancedVolcano(M5cKO_MWT,
                                lab = row.names(M5cKO_MWT),
                                x = 'log2FoldChange',
                                y = 'padj',
                                selectLab = top10_M5cKO$gene,
                                xlim=c(-1,1),
                                ylim=c(0,10),
                                xlab = bquote(~Log[2]~ 'fold change'),
                                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                                pCutoff = 0.1,
                                FCcutoff = 0.1,
                                title = "Kdm5c-KO Male P6 Forebrain RNA-seq",
                                labSize = 4.0,
                                colAlpha = 4/5,
                                colCustom = keyvals,
                                #col=c('lightsteelblue4', 'lightsteelblue4', 'lightsteelblue4', 'red3'),
                                legend=c('NS','Log2 FC','Adjusted p-value',
                                         'Adjusted p-value & Log2 FC'),
                                legendPosition = 'bottom',
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE,
                                legendLabSize = 10,
                                legendIconSize = 3.0)


VolcM5cKO_all




# pdf(file = "C:/Users/kmbon/Documents/Projects/KDM5C_GermlineGenes/KDM5C_EpiLC_RNAseq/PageTestisGenesAnalysis/GermVolcanoPlot_HIP_page_zoom.pdf",   # The directory you want to save the file in
#     width = 7, # The width of the plot in inches
#     height = 8) # The height of the plot in inches
# VolcHIP_page
# dev.off()


### female

#all of the p6 Female genes in the results table
head(F5cHet_FWT)

#get the mgi names
# F5cHet_FWT_mgi <- bitr(row.names(F5cHet_FWT), fromType = "ENSEMBL",
#                 toType = c("SYMBOL"),
#                 OrgDb = org.Mm.eg.db)





keyvals <- ifelse(
  row.names(F5cHet_FWT) %in% Li_TeOvBr_Test$gene_id, '#6BA96B',
  ifelse(row.names(F5cHet_FWT) %in% Li_TeOvBr_Ovary$gene_id, 'red',
         ifelse(row.names(F5cHet_FWT) %in% Li_TeOvBr_TestandOvary$gene_id, 'black',
                ifelse(row.names(F5cHet_FWT) %in% Li_TeOvBr_Brain$gene_id, 'blue',
                       ifelse(F5cHet_FWT$log2FoldChange >= 0.1 & F5cHet_FWT$padj <= 0.1, '#F6B000',
                              ifelse(F5cHet_FWT$log2FoldChange <= -0.1 & F5cHet_FWT$padj <= 0.1, '#6991D1',
                                     'gray'))))))


keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#F6B000'] <- 'Up DEGs'
names(keyvals)[keyvals == '#6991D1'] <- 'Down DEGs'
names(keyvals)[keyvals == '#6BA96B'] <- 'Testis DEGs'
names(keyvals)[keyvals == 'red'] <- 'Ovary DEGs'
names(keyvals)[keyvals == 'black'] <- 'Testis and Ovary DEGs'
names(keyvals)[keyvals == 'blue'] <- 'Brain DEGs'
names(keyvals)[keyvals == 'gray'] <- 'NS'


#making shapes outside the axis range
df$shape <- ifelse(abs(df$values)>5, "triangle", "circle")
df$values[df$values>5] <- 5
df$values[df$values< -5] <- -5

# F5cHet_FWT.sort = F5cHet_FWT.UPDEG[order(F5cHet_FWT.UPDEG[,3]), ]
# top10_F5cHet <- tail(F5cHet_FWT.sort, 10)

VolcF5cHet_all <- EnhancedVolcano(F5cHet_FWT,
                                 lab = row.names(F5cHet_FWT),
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 selectLab = top10_F5cHet$gene,
                                 xlim=c(-1,2),
                                 ylim=c(0,10),
                                 xlab = bquote(~Log[2]~ 'fold change'),
                                 ylab = bquote(~-Log[10]~adjusted~italic(P)),
                                 pCutoff = 0.1,
                                 FCcutoff = 0,
                                 title = "Kdm5c-Het Female P6 Forebrain RNA-seq",
                                 labSize = 4.0,
                                 colAlpha = 1/5,
                                 colCustom = keyvals,
                                 #col=c('lightsteelblue4', 'lightsteelblue4', 'lightsteelblue4', 'red3'),
                                 legend=c('NS','Log2 FC','Adjusted p-value',
                                          'Adjusted p-value & Log2 FC'),
                                 legendPosition = 'bottom',
                                 gridlines.major = FALSE,
                                 gridlines.minor = FALSE,
                                 legendLabSize = 10,
                                 legendIconSize = 3.0)


VolcF5cHet_all


