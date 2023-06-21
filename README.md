# RaukR_project

![](Untitled.png){width="321"}

**NO** a shiny and an associated package to analyse RNASeq data that has been mapped to a transcriptome with [Salmon](https://salmon.readthedocs.io/en/latest/#)

**NO** will accept salmon (.sf) as well as count matrices as input for the [DESeq2 analysis](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

From the DESeq2 output, genes that are differentially expressed are extracted and: 
- Compared to NDEx Integrated Query to find genes related to differentially expressed genes 
- Pathways impacted by differentially expressed genes are annotated from the [KEGG database](https://www.genome.jp/kegg/) 
- Genes of interests can be searched in cytoscape window without leaving the shiny app


