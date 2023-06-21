# Install and load clusterProfiler
#install.packages("clusterProfiler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(tidyverse)
library(readr)
library(enrichplot)
library(data.table)

#Import dataset
deseq2_results <- read_csv("deseq2_results.csv")
View(deseq2_results)

colnames(deseq2_results)[colnames(deseq2_results) == "...1"] <- "ID"
deseq2_results_2 <- deseq2_results %>% drop_na(log2FoldChange)

# Define the thresholds for overexpression/underexpression based on DE genes with a minimum 2-fold change of expression at a maximum FDR of 5%
up <- deseq2_results$padj < 0.05 &
  abs(deseq2_results$log2FoldChange) > log2(2)
deGenesUp <- deseq2_results[up, ]
head(deGenesUp)

down <- deseq2_results$padj < 0.05 &
  abs(deseq2_results$log2FoldChange) < log2(2)
deGenesDown <- deseq2_results[down, ]
head(deGenesDown)

#Remove NAs from deGenes
deGenesUp_2 <- deGenesUp %>% drop_na(ID)
deGenesUp_3 <- deGenesUp_2$ID

deGenesDown_2 <- deGenesDown %>% drop_na(ID)
deGenesDown_3 <- deGenesDown_2$ID

#Create geneUniverse
geneUniverse <- deseq2_results$ID


#KEEG enrichment
ans.kegg <- enrichKEGG(gene = deGenes_3,
                       organism = 'sce', 
                       universe = geneUniverse,
                       pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(ans.kegg)
tab.kegg <- subset(tab.kegg, Count>5)
tab.kegg[1:5, 1:6]

#unlist(split(ans.kegg@result$Description, "-" )) 



#Plot KEGG results
p1 <- dotplot(ans.kegg, showCategory=10) + ggtitle("KEGG")
p1


x2 <- pairwise_termsim(ans.kegg) 
p2 <- emapplot(x2, showCategory = 10, cex.params = list(category_node = 1, category_label = 0.6, line = 0.5))
p2

# Perform gene set enrichment analysis using cluster profiler package
#enrich_result <- enricher(gene = deGenes_2["entrez_gene_id"],
#                          TERM2GENE = NULL,
#                          universe = NULL,
#                          pvalueCutoff = 0.05,
#                          pAdjustMethod = "BH",
#                          qvalueCutoff = 0.05,
#                          )

# Visualize the enriched categories
#barplot(enrich_result, showCategory = 10)


