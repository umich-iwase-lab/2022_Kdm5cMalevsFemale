#MA Plots 22.03.08

# Basic scatter plot
#male MA plot

head(M5cKO_MWT.tb)

ggplot(M5cKO_MWT.tb, aes(x=baseMean, y=log2FoldChange)) +
  xlim(c(0.1,100000)) + ylim(c(-3,3)) +
  geom_point(size=2)



summary(F5cHet_FWT)


#try with ggpubr
library(ggpubr)

# write.table(M5cKO_MWT.tb, file = "M5cKO_MWT_gene.txt", sep = "\t")
M5cKO_MWT_gene <- read.csv(file = "M5cKO_MWT_gene.txt", sep = "\t")

# write.table(F5cHet_FWT.tb, file = "F5cHet_FWT_gene.txt", sep = "\t")
F5cHet_FWT_gene <- read.csv(file = "F5cHet_FWT_gene.txt", sep = "\t")


genmaplot <- function(DESEQ){
  ggmaplot(DESEQ, 
           fdr = 0.1, fc = 0, size = 0.6,
           palette = c("#f2b722ff", "#1465AC", "darkgray"),
           genenames = as.vector(DESEQ$gene),
           legend = "top", top = 5,
           font.label = c("bold", 11), label.rectangle = TRUE,
           font.legend = "bold",
           ylim = c(-1, 2.5),
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())
}


pdf(file = "MAPlot_P6MaleBrain.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches
genmaplot(M5cKO_MWT_gene)
dev.off()


pdf(file = "MAPlot_P6FemaleBrain.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches
genmaplot(F5cHet_FWT_gene)
dev.off()



