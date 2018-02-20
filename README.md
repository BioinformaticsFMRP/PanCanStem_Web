## Introduction

This page contains...  

### Workflows

[http://tcgabiolinks.fmrp.usp.br/PanCanStem/](http://tcgabiolinks.fmrp.usp.br/PanCanStem/) contains workflows, explaining how to download stem/progenitor cells from the Progenitor Cell Biology Consortium(PCBC) and TCGA PanCan33 datasets, train a stemness signature using normal stem cells, and apply it to score TCGA tumor samples. The mRNAsi section contains the RNA expression-based stemness index workflow; while the mDNAsi contains the DNA methylation-based stemness index workflow.

### PanCanStem Web

To evaluate our stemness score, we performed an enrichement analysis for clinical or molecular features.
The results of this evaluation were summarized and made available through a shiny app
to help the results data mining.

This shiny app can be accessed either by a [hosted version](http://143.107.143.246:3838/PanCanStem_Web/)  
or run a local version of PanCanStem Web with the following code in R
```r
shiny::runGitHub("BioinformaticsFMRP/PanCanStem_Web")
```

