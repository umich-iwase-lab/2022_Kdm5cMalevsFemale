#Generates bargraph of KDM5C tpm compared across different groups
# KDM5C expression in P6 Males and females
library(ggplot2)

P6_Brain_TPM <- read.csv(file = "~/Projects/21.09.28_P6MaleandFemaleBrains/ms_P6brains_readcount_210929_tpm_DESEQ2.txt", header = T, sep = "\t")
head(P6_Brain_TPM)
nrow(P6_Brain_TPM)
ncol(P6_Brain_TPM)

head(P6_Brain_TPM)
P6_Brain_TPM$Geneid
View(P6_Brain_TPM)

#pull out Kdm5c
KDM5C_tpm <- subset(P6_Brain_TPM, P6_Brain_TPM$Symbol == "Kdm5c")
#transpose
KDM5C_tpm.t <- t(KDM5C_tpm)
KDM5C_tpm <- as.double(KDM5C_tpm.t[3:18,])

#read in the genotypes of the samples from DESeq2 analysis
samplegeno <- SampleInfo
colnames(samplegeno) <- c("Sample", "Genotype", "Sex")
#make a new group called genosex that includes both the genotype and the sex info
samplegeno <- data.frame(samplegeno, GenoSex = paste(samplegeno$Genotype, samplegeno$Sex))

KDM5C_tpm <- data.frame(samplegeno, KDM5C_tpm)
#set the order
KDM5C_tpm$GenoSex <- factor(KDM5C_tpm$GenoSex, levels = c("WT M", "5cKO M", "WT F","5cHet F"))

kdm5ctpm_all <- ggplot(KDM5C_tpm, aes(x=GenoSex, y=KDM5C_tpm, fill=GenoSex)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), dotsize=1) +
  ggtitle("KDM5C TPM P6 Forebrain") +
  xlab(" ") + ylab("TPM") + scale_fill_manual("Genotype", values = c("WT M" = "#21315bff", "5cKO M" = "#6179b5ff","WT F" = "#892717ff", "5cHet F" = "#e36752ff" ))+
  theme(text = element_text(size = 20))   

pdf(file = "TPM_KDM5C_all.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches
kdm5ctpm_all
dev.off()

##Wild type only - comparing XX vs XY expression of Kdm5c

KDM5C_tpm_WT <- subset(KDM5C_tpm, Genotype == "WT")
write.csv(KDM5C_tpm_WT, "KDM5C_tpm_WT.csv")

kdm5ctpm_WT <- ggplot(KDM5C_tpm_WT, aes(x=GenoSex, y=KDM5C_tpm, fill=GenoSex)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), dotsize=1) +
  ggtitle("KDM5C TPM P6 Forebrain") +
  ylim(c(0,85)) +
  xlab(" ") + ylab("TPM") + scale_fill_manual("Genotype", values = c("WT M" = "#21315bff", "WT F" = "#892717ff"))+
  theme(text = element_text(size = 20))   

pdf(file = "TPM_KDM5C_WT.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
kdm5ctpm_WT
dev.off()



#test for statistical difference using welch's t-test
KDM5C_tpm_WT_M <- subset(KDM5C_tpm_WT, Sex == "M")
KDM5C_tpm_WT_F <- subset(KDM5C_tpm_WT, Sex == "F")
#statistiscs
Kdm5c_WTFvsM <- t.test(KDM5C_tpm_WT_M$KDM5C_tpm, KDM5C_tpm_WT_F$KDM5C_tpm)


