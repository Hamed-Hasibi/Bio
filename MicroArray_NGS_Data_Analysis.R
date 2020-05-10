#setRepositories()
#add BioC software => 1 2
#install.packages(c("GEOquery", "limma", "pheatmap", "ggplot2", "gplots", "reshape2", "plyr"))

setwd("~/Desktop/Bio/R/")

library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)

series <- "GSE9476"
platform <- "GPL96"
####### Section1: Quality check

#### load data
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gr <- c("CD34",rep("BM",10),rep("CD34",7),rep("AML",26),rep("PB",10),rep("CD34",10))

ex <- exprs(gset)

##### log2 sacle if req.
# ex <- log2(ex+1)
# exprs(gset) <- ex

pdf("Results/boxplot.pdf", width = 64, height = 64)
boxplot(ex)
dev.off()

#### normalize if req.
# temp <- normalizeQuantiles(ex)
# boxplot(temp)

### correlation heatmap
pdf("Results/CoHeatmap.pdf", width = 15, height = 15)
pheatmap(cor(ex), labels_row = gr, labels_col = gr, color = bluered(256), border_color = NA)
dev.off()

### principal component analysis (PCA)
pc <- prcomp(ex)
pdf("Results/pc.pdf")
plot(pc)
plot(pc$x[,1:2])
#names(pc)
#dim(px$x)
dev.off()

ex.scale <- t(scale(t(ex), scale = F))

### scale makes every column to have mean = 0 (actually new value = old value - mean)
### it also for scale = 0, divides each for SD (Standard divation)

pc <- prcomp(ex.scale)
pdf("Results/Pc-scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

pcr <- data.frame(pc$rotation[,1:3], Group=gr)
pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2 , color=Group)) + geom_point(size=3) + theme_bw()
dev.off()


#### section2: Differential Expression analysis

gr <- factor(gr)
gset$description <- gr
design <- model.matrix(~ description + 0, gset)
#head(design)
colnames(design) <- levels(gr)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(AML-CD34, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
write.table(tT, "Results/AML_CD34.txt", row.names=F, sep="\t", quote = F)

#### secion3: analysis
aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file = "Results/AML_CD34_up.txt", quote = F, row.names=F, col.names=F)

aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file = "Results/AML_CD34_down.txt", quote = F, row.names=F, col.names=F)






 