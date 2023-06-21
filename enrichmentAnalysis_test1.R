library(tidyverse)
#Import dataset
library(readr)
deseq2_results <- read_csv("deseq2_results.csv")
View(deseq2_results)

colnames(deseq2_results)[colnames(deseq2_results) == "...1"] <- "ID"

#Load Annotation file for yeast
#install.packages("ape")
library(ape)
annotationFile<-read.gff("GCF_000146045.2_R64_genomic.gff")


tidy_attributes <- annotationFile %>%
  filter(type == "gene") %>%
  separate(col = attributes, into = c("ID", "GeneID", "Name", "endRage", "gbkey", "gene_biotype", "locus_tag", "partial", "start_range"), sep = ";") %>%
  separate(col = GeneID, into = c("trash", "entrez_gene_id"), sep = ":") %>% separate(col = ID, into = c("trash", "ID"), sep = "ID=gene-") 


#Annotate results
annotated_results <- merge(deseq2_results, tidy_attributes,
                           by.x= "ID", all.x = TRUE) 

#
annotated_results_2 <- annotated_results %>% drop_na(entrez_gene_id) #remove NAs
geneUniverse <- annotated_results_2$entrez_gene_id
#
annotated_results_test <- annotated_results %>% drop_na(ID) #remove NAs
geneUniverse_test <- annotated_results_2$ID

# Define the thresholds for overexpression/underexpression based on DE genes with a minimum 2-fold change of expression at a maximum FDR of 5%
mask <- annotated_results$padj < 0.05 &
  abs(annotated_results$log2FoldChange) > log2(2)
deGenes <- annotated_results[mask, ]
head(deGenes)

#Remove NA from entrezID
deGenes_2 <- deGenes %>% drop_na(entrez_gene_id)

deGenes_3 <- deGenes_2[,"entrez_gene_id"]


#
deGenes_2_test <- deGenes %>% drop_na(ID)

deGenes_3_test <- deGenes_2_test[,"ID"]

# Install and load clusterProfiler
#install.packages("clusterProfiler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

library(clusterProfiler)


#KEEG enrichment
ans.kegg <- enrichKEGG(gene = deGenes_3_test,
                       organism = 'sce', 
                       universe = geneUniverse_test,
                       pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(ans.kegg)
tab.kegg <- subset(tab.kegg, Count>5)
tab.kegg[1:5, 1:6]


#Plot KEGG results
p1 <- dotplot(ans.kegg, showCategory=10) + ggtitle("KEGG")
p1

library(enrichplot)
x2 <- pairwise_termsim(ans.kegg) 
emapplot(x2, showCategory = 10, cex.params = list(category_node = 1, category_label = 0.6, line = 0.5))

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


