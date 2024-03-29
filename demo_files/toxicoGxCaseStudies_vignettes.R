# ---
#   title: 'ToxicoGx: Case Studies'
# author:
#   - Sisira Kadambat Nair
# - Esther Yoo
# - Christopher Eeles
# - Nehme El-Hachem
# - Petr Smirnov
# - Benjamin Haibe-Kains, benjamin.haibe.kains@utoronto.ca
# date: "15-10-2019"

## Plotting gene expression changes from TG-GATEs dataset
  
# Our first case study is based on a plot from a paper by Rangel-Escareño et al., in which gene expression changes of CYP1A1 (gene associated with xenobiotic metabolism) has been plotted at all concentrations and time points. The plot shows clear differential expression at time 8(hr) suggesting that higher the dose, larger the impact of CCL4 on this gene.
# For plotting the gene expression under same conditions using the package, the first step is to load the datasets from disk or download them using the downloadTSet function above. In the following example, we use the toy dataset provided with the package to illustrate the process. 
# To plot, the function drugGeneResponseCurve has been used wherein mandatory inputs such as dataset, drug name, cell-line, molecular type, gene name, dose and time points should be specified.

require(devtools)
devtools::install_github("bhklab/ToxicoGx", ref = "CRAN_submission")

#**Plot time dependent dose response of Carbon tetra chloride on CYP1A1 gene.**
library(ToxicoGx)
# Load the tset 
TGGATESsmall <- readRDS("TGGATESsmall.rds")
drugGeneResponseCurve(tSets = TGGATESsmall, duration = c("2", "8", "24"), 
                      cellline = "Hepatocyte", mDataTypes = "rna", 
                      features = "ENSG00000140465_at", 
                      dose = c("Control", "Low", "Middle","High"),
                      drug = "carbon tetrachloride", 
                      plot.type = "Actual",
                      title = "Effect of Carbon tetra chloride on CYP1A1.",
                      cex = 0.5, cex.main = 1, legend.loc = "topright",
                      mycol = c("red", "green", "dark blue", "turquoise"),
                      verbose = T)



## Connectivity map analysis on TG-GATEs and human hepatocarcinoma signatures.

# For the second case study, we will recreate an analysis from the paper by Jos Kleinjans et al., wherein connectivity mapping has been used to predict compound carcinogenicity
# by linking in vivo human hepatocarcinoma (HCC) signature geneset with in vitro TG-GATEs primary human hepatocyte data. In this example, we are using the toy dataset. The full dataset has to be downloaded to carry out the whole analysis done in the paper.
# The HCC signature, already mapped to the gene level, has been included in this package and it can be loaded by calling data(HCC_sig). Once the dataset is loaded, recreate drug signatures for each drug using the function drugPerturbationSig to perform statistical modelling of the transcriptomic response to the application of each drug. We then compare the observed up-regulated and down-regulated genes to HCC signature published in the paper. The output will be the GSEA connectivity score with FDR values that can be used to determine the correlation between the two signatures.

require(xtable)
# PharmacoGx package is required for this vignette
library(PharmacoGx)

TGGATESsmall <- readRDS("TGGATESsmall.rds")
# To compute the effect of drug concentration on the molecular profile of the cell
# drug.perturbation <- drugPerturbationSig(TGGATESsmall,
#                                          mDataType = "rna", cells = "Hepatocyte",
#                                          duration = "24",
#                                          dose = c("Control", "Low"),
#                                          drugs = drugNames(TGGATESsmall),
#                                          verbose = FALSE)
drug.perturbation <- readRDS("drug.perturbation.rds")
HCC_sig <- readRDS("HCC_sig.rds")
res <- apply(drug.perturbation[,,c("tstat", "fdr")],
             2, function(x, HCC){
               return(PharmacoGx::connectivityScore(x=x,
                                                    y=HCC[,2,drop=FALSE],
                                                    method="fgsea", nperm=100))
             }, HCC=HCC_sig)
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <-  cbind(res,"FDR" = p.adjust(res[,2], method="fdr"))
res <- res[order(res[,3]),]
View(res)
xtable::xtable(res,
               caption='Connectivity Score results for HCC and TG-GATEs PHH gene
       signature.')
```

#In the table, omeprazole showed a positive connectivity score in contrast to isoniazid. This observation aligns with the trends reported in the paper. The above example is to demonstrate the ease with which drug perturbation analysis can be done using ToxicoGx package.
