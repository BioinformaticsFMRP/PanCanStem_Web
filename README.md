## Introduction

This page contains additional information to the TCGA PanCanAtlas Stemness project. Here we provide two components: 1) the Stemness Index Workflow to reproduce the generation of our Stemness Indices (mRNAsi and mDNAsi). This workflow can also be applied to additional dataset (e.g., users own samples) to derive corresponding stemness indices. 2) the PanCanStemness Web portal, which contains the enrichment analysis of clinical and molecular features in relation to the Stemness indices. 


### Stemness Index Workflow

This is a step by step  explanation on how to download stem/progenitor cells from the Progenitor Cell Biology Consortium(PCBC) and TCGA PanCanAtlas datasets, train a stemness signature using normal stem cells, and apply it to score TCGA tumor samples in order to measure the stemness status. The mRNAsi section contains the RNA expression-based stemness index workflow; while the mDNAsi contains the DNA methylation-based stemness index workflow. After stepping through the workflow, we encourage the user to replace TCGA PanCanAtlas dataset with their own samples to derive the corresponding Stemness Indices.

Link: [http://tcgabiolinks.fmrp.usp.br/PanCanStem/](http://tcgabiolinks.fmrp.usp.br/PanCanStem/) 

### PanCanStemness Web portal

To evaluate the performance of our stemness indices across the entire TCGA cohort, we performed an enrichment analysis by sorting TCGA samples by stemness index for each tumor type and looked for associations with all available genomic features, molecular features, and clinical features. We used the fgsea R/Bioconductor package to compute the enrichment scores (Sergushichev 2016). Briefly, for each tumor type we ranked the TCGA samples according to their stemness index (from -low to -high stemness index) and tested if any particular genomic/molecular/clinical feature was associated with either -low or -high stemness index in a non-random behavior. We performed 10,000 permutations for each parameter analyzed to calculated our enrichment score. We then normalized the enrichment scores to mean enrichment of random samples of the same size (NES - normalized enrichment score).The results of this evaluation were summarized and made available through a shiny app to help the results data mining.

To evaluate our stemness score, we performed an enrichement analysis for clinical or molecular features.
The results of this evaluation were summarized and made available through a shiny app
to help the results data mining.

This shiny app can be accessed either by a [hosted version](http://143.107.143.246:3838/PanCanStem_Web/)  
or run a local version of PanCanStem Web with the following code in R
```r
shiny::runGitHub("BioinformaticsFMRP/PanCanStem_Web/",subdir = "PanCanStem_Web")
```

