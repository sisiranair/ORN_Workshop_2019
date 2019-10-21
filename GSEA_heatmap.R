library(ToxicoGx)
library(BiocGenerics)
library(Biobase)
library(piano)
library(ComplexHeatmap)

tggates_humanLdh <- readRDS("data/TGGATES_humanLDH_tSet.rds")

drug.pert <- drugPerturbationSig(tggates_humanLdh,mDataType = "rna", drugs =  drugNames(tggates_humanLdh)
                                 , cells = "Hepatocyte", features = fNames(tggates_humanLdh, mDataType = "rna"), dose = c("Control","Low", "Middle", "High")
                                 , duration = c("2","8","24"), verbose = T)
#extract estimate values
est_drug.pert <- drug.pert[,,1]
est_drug.pert <- as.data.frame(est_drug.pert[apply(est_drug.pert, 1, function(x)!any(is.na(x))),, drop=F])

#to map to symbol name

ff <- fData(tggates_humanLdh@molecularProfiles$rna)


est_drug.pert$rnames <- rownames(est_drug.pert)

est_drug.pert_merged <- merge(est_drug.pert, ff, by.x = "rnames", by.y = "gene_id", all.x = T)
#subset till symbol
est_drug.pert_sub <- subset(est_drug.pert_merged, select = c(1:148), drop=F)

est_drug.pert_sub <- subset(est_drug.pert_sub, subset = !is.na(est_drug.pert_sub$Symbol), drop=F)

rownames(est_drug.pert_sub) <- est_drug.pert_sub$Symbol

est_drug.pert_allDrugs <- as.matrix(est_drug.pert_sub[,c(2:147)])

#function for GSEA - run_GSA 
source("TGx_function_testing/GSEA.R")

#gmt file - Reactome GS V7
SigDb_gene_symbols_Reactome <- "data/c2.cp.reactome.v7.0.symbols.gmt"

op_list <- list()

for(dr in 1:ncol(est_drug.pert_allDrugs)){
  
  gseaOutput_Reactome <- run_GSEA(est_drug.pert_allDrugs[,dr], gmt.file = SigDb_gene_symbols_Reactome,return.piano=FALSE)
  
  op_list[[colnames(est_drug.pert_allDrugs)[dr]]] <- gseaOutput_Reactome
}

saveRDS(op_list, "data/op_list.Reactome.rds")

###################################################----FDR----###################################################################
#output of reactome GSEA on all drugs
op_list.Reactome <- readRDS("data/op_list.Reactome.rds")

#extracting Stat (dist.dir) and re arranging GSEA output.
gsea_op <- do.call(cbind,lapply(op_list.Reactome, function(x){
  
  rownames(x) <- x[,"Name"]
  return(x[sort(x[,"Name"]),"Stat (dist.dir)"])
}))


colnames(gsea_op) <- names(op_list.Reactome)
rownames(gsea_op) <- sort(op_list.Reactome$`naphthyl isothiocyanate`[,"Name"])

#extracting FDR.

gsea_op_fdr <- do.call(cbind,lapply(op_list.Reactome, function(x){
  
  rownames(x) <- x[,"Name"]
  return(x[sort(x[,"Name"]),"fdr"])
}))

colnames(gsea_op_fdr) <- names(op_list.Reactome)
rownames(gsea_op_fdr) <- sort(op_list.Reactome$`naphthyl isothiocyanate`[,"Name"])


#replace all values above 0.05 as 0
gsea_op_zero <- gsea_op

for(cl in 1:ncol(gsea_op_fdr)){
  for(rn in 1:nrow(gsea_op_fdr)){
    if (gsea_op_fdr[rn,cl] > 0.05){
      gsea_op_zero[rn,cl] <- 0
    }
  } 
}

#to filter out rows with no fdr significant patwhays

rf <- apply(gsea_op_fdr, 1, function(x){sum(x<0.05)})

gsea_op_zero_rowfilter <- gsea_op_zero[which(rf > 50),]
gsea_op_zero_rowfilter <- gsea_op_zero_rowfilter[,colSums(gsea_op_zero_rowfilter) != 0]


pdf("FDR selected pathways filter > 50est.pdf", width = 20, height = 20)
ComplexHeatmap::Heatmap(gsea_op_zero_rowfilter, name = "FDR selected pathways", row_names_gp = gpar(fontsize = 2),cluster_rows = T,cluster_columns = T
                        , column_names_gp = gpar(fontsize = 5), na_col = "grey", cell_fun = function(j, i, x, y, width, height, fill) {
                          if(gsea_op_zero_rowfilter[i, j] == 0)
                            grid.rect(x = x, y = y, width = width, height = height, 
                                      gp = gpar(col = "grey", fill = "grey"))
                        })
dev.off()





