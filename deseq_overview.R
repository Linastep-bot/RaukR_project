library(tximport)
library(GenomicFeatures)
library(readr)
library(readxl)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)

setwd("/Documents and Settings/llpo0001/Documents/slubi/99_rnaseq/")
txdb <- makeTxDbFromGFF("DM8.1_gene.gff3.gz")
k <- keys(txdb, keytype = "GENEID")
tx3gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
head(tx3gene)


samples <- read_xlsx("metadata2.xlsx")
#files <- file.path("slubi/rnaseq/3_salmon", samples$name, ".sf")
files <- list.files("salmon/3_salmon/",pattern=".sf",recursive = TRUE,full.names = TRUE)
names(files) <- paste0(samples$name)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx3gene, ignoreTxVersion = TRUE)


dds_sample_type <- DESeqDataSetFromTximport(txi.salmon, samples, ~sample_type)
dds_sample_type <- DESeq(dds_sample_type)
res <- results(dds_sample_type)

plotDispEsts(dds_sample_type, main="Dispersion plot")

rld_sample_type <- rlogTransformation(dds_sample_type)
head(assay(rld_sample_type))

(mycols <- brewer.pal(8, "Dark2")[1:length(unique(samples$sample_type))])
sampleDists <- as.matrix(dist(t(assay(rld_sample_type))))
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[samples$line],
          RowSideColors=mycols[samples$line],
          margin=c(10, 10), main="Sample Distance Matrix Sample Type")

pcaData = plotPCA(rld_sample_type, intgroup=c("line","sample_type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=line, shape=sample_type, label = name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA of DESeq2 from Salmon") +geom_text(hjust=0, vjust=0)
#facet_wrap(~ group) 

table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds_sample_type, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

hist(res$pvalue, breaks=50, col="grey")
DESeq2::plotMA(dds_sample_type, ylim=c(-1,1), cex=1)

# Volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
