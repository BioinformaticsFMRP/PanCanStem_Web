## Introduction

This page contains additional information to support the TCGA PanCanAtlas Stemness project. We provide two components: 

1. **Stemness Index Workflow** provides the steps and methods to regenerate our Stemness Indices (mRNAsi and mDNAsi). This workflow also provides steps to generate corresponding stemness indices for non-TCGA samples (e.g., users own samples). 
2. **PanCanStemness Web portal** contains the enrichment analysis reported in our study.  


### Stemness Index Workflow

This workflow describes the steps we developed to generate our Stemness Indices (mRNAsi and mDNAsi).The step-by-step explanation includes a) how to download stem/progenitor cells from the Progenitor Cell Biology Consortium(PCBC), b)download TCGA PanCanAtlas datasets, c) train a stemness signature using normal stem cells, and d) apply the one-class algorithm to define a stemness index for each  tumor sample. The mRNAsi section contains the RNA expression-based stemness index workflow; while the mDNAsi contains the DNA methylation-based stemness index workflow. Users can replace TCGA PanCanAtlas dataset with their own samples to derive the corresponding Stemness Indices.

Link: [http://tcgabiolinks.fmrp.usp.br/PanCanStem/](http://tcgabiolinks.fmrp.usp.br/PanCanStem/) 

### PanCanStemness Web portal

To evaluate the association of our stemness indices to known molecular sub-classfication, mutation events and/or clinical features (such as age, survival, treatment, etc.), we performed a statistical enrichment analysis by harnessing the fgsea R/Bioconductor package (Sergushichev 2016). By sorting TCGA samples by each stemness index within each tumor type, we were able to find associations with all available molecular, and clinical features. Briefly, for each tumor type we ranked the TCGA samples according to their stemness index (from -low to -high stemness index) and tested if any particular molecular/clinical feature was associated with either -low or -high stemness index in a non-random behavior. We performed 10,000 permutations for each parameter analyzed to calculate our enrichment score. We then normalized the enrichment scores to mean enrichment of random samples of the same size (NES - normalized enrichment score).



To access the results of this enrichment analysis, one can either visit a [hosted version](http://143.107.143.246:3838/PanCanStem_Web/)  
page or run a local version of the results by executing  the following R code.
```r
shiny::runGitHub("BioinformaticsFMRP/PanCanStem_Web/",subdir = "PanCanStem_Web")
```
We recommend running the R code on your local machine for faster access.

