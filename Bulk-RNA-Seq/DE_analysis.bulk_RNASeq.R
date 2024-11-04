```{r load.packages, echo=FALSE, message=FALSE, results='hide'}
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(colortools))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(gage))
suppressPackageStartupMessages(library(HTSanalyzeR))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(pheatmap))

options(xtable.type="html")

knitr::opts_chunk$set(
  echo=FALSE,
  warning=FALSE,
  message=TRUE,
  error = FALSE,
  tidy = FALSE,
  cache = TRUE,
  dev = c("png", "pdf"),
  fig.width=7, fig.height=7)

options(markdown.HTML.stylesheet = 'css/custom.css')
```
#geneCounts = as.matrix(read.csv("geneCounts.csv", header=TRUE, row.names=1, check.names=FALSE))
setwd(".")
data <- read.table("all_reps.hg38_unique.featureCounts.txt2", header=F,sep="\t")

gene=data[,1]
length=data[,55]
count <- data[,2:54]

info1 = read.table("NRXN1_Bulk_RNA_info.csv", header=T, sep=",")
info1$Donor = as.factor(info1$Donor)
rownames(info1) = info1$ID

cpm <- cpm(count)
keep <- rowSums(cpm(count) > 1) >= 6
geneCounts <- count[keep,]
length2 <- length[keep]
gene2 <- gene[keep]
cpm2 <- cpm(geneCounts, log=T)

colnames(cpm) = info1$ID
rownames(cpm) = gene

colnames(geneCounts) = info1$ID
rownames(geneCounts) = gene2

##
geneInfo = read.table("ENSEMBLv70_gene_info.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t")
#
#```{r voom}
genes = DGEList(counts=geneCounts)
genes = calcNormFactors(genes)
design = model.matrix(~Time_point, info1)

#
vobj_init = voom(genes, design, plot=FALSE )

dupcor_init = duplicateCorrelation( vobj_init, design, block=info1$Donor)
# dupcor_init$consensus

# # voom / dupcor 2
# #################
# # re-estimate weights using correlation values
vobj <- voom(genes, design, block = info1$Donor, correlation=dupcor_init$consensus)

dupcor = duplicateCorrelation( vobj, design, block=info1$Donor)
# dupcor$consensus
```
form = ~ Group + Time_point + Cell_type + MEF + FB

design = model.matrix( form, info1)
fit = lmFit( vobj, design, block=info1$Donor, correlation=dupcor$consensus )
fit = eBayes(fit)
```
## Differential expression within each cell type
```{r DEsubsets}
info1$Interaction = with(info1, gsub(" ", ".", paste( Group, Time_point, sep='_')))

form = ~ Interaction - 1 + MEF + FB + hiP
design = model.matrix( form , info1)

#SZinteraction = (Interaction3_Prime_Del_DIV14 - InteractionCT_NPC) - (InteractionSZ_6_wk_FB_neuron - InteractionCT_6_wk_FB_neuron),
cons_separate_cell_type = makeContrasts(
    NRXN1_3PrimeDel_DIV14 = Interaction3_Prime_Del_DIV14 - InteractionControl_DIV14,
    NRXN1_3PrimeDel_DIV21 = Interaction3_Prime_Del_DIV21 - InteractionControl_DIV21,
    NRXN1_3PrimeDel_DIV35 = Interaction3_Prime_Del_DIV35 - InteractionControl_DIV35,
    NRXN1_5PrimeDel_DIV14 = Interaction5_Prime_Del_DIV14 - InteractionControl_DIV14,
    NRXN1_5PrimeDel_DIV21 = Interaction5_Prime_Del_DIV21 - InteractionControl_DIV21,
    NRXN1_5PrimeDel_DIV35 = Interaction5_Prime_Del_DIV35 - InteractionControl_DIV35,
    NRXN1_5vs3_DIV14 = Interaction5_Prime_Del_DIV14 - Interaction3_Prime_Del_DIV14,
    NRXN1_5vs3_DIV21 = Interaction5_Prime_Del_DIV21 - Interaction3_Prime_Del_DIV21,
    NRXN1_5vs3_DIV35 = Interaction5_Prime_Del_DIV35 - Interaction3_Prime_Del_DIV35,
    levels = colnames(design))

fit = lmFit( vobj, design, block=info1$Donor, correlation=dupcor$consensus )

res0 <- residuals( fit, vobj)

fit2 <- contrasts.fit(fit, cons_separate_cell_type)
fit2 = eBayes(fit2)

resList = list()
for(key in colnames(cons_separate_cell_type)){
  #res = topTable(fit2, coef=key, number=Inf, sort.by='P')
  res = topTable(fit2, coef=key, number=Inf, sort.by='none')
  res$qvalue = qvalue( res$P.Value )$qvalue
  res$gene = geneInfo$geneName[match(rownames(res), geneInfo$Geneid)]
  resList[[key]] = res
}

for( key in names(resList)){
  ##file = paste0("results/differential_expression/Schizophrenia/COS_DE_", key, ".tsv")
  file = paste0("NRXN1_bulk_RNA_DE_regressHiPSC_", key, ".csv")
  write.table( data.frame(geneName = rownames(resList[[key]]), resList[[key]]), file, quote=FALSE, sep="\t", row.names=FALSE)
}
```
save.image("NRXN1.RNA_seq.DEG_analysis.RData")
load("NRXN1.RNA_seq.DEG_analysis.RData")

##Call DEGs within each pair of donors.
info1$Interaction2 = with(info1, gsub(" ", ".", paste( Group, Time_point, Donor, sep='_')))

form2 = ~ Interaction2 - 1 + MEF + FB + hiP

design = model.matrix( form2 , info1)

cons_separate_donor_celltype = makeContrasts(
    NRXN1_3PrimeDel_DIV14_581_2607 = Interaction23_Prime_Del_DIV14_581 - Interaction2Control_DIV14_2607,
    NRXN1_3PrimeDel_DIV14_641_2607 = Interaction23_Prime_Del_DIV14_641 - Interaction2Control_DIV14_2607,
    NRXN1_3PrimeDel_DIV14_581_553 = Interaction23_Prime_Del_DIV14_581 - Interaction2Control_DIV14_553,
    NRXN1_3PrimeDel_DIV14_641_553 = Interaction23_Prime_Del_DIV14_641 - Interaction2Control_DIV14_553,
    NRXN1_3PrimeDel_DIV21_581_2607 = Interaction23_Prime_Del_DIV21_581 - Interaction2Control_DIV21_2607,
    NRXN1_3PrimeDel_DIV21_641_2607 = Interaction23_Prime_Del_DIV21_641 - Interaction2Control_DIV21_2607,
    NRXN1_3PrimeDel_DIV21_581_553 = Interaction23_Prime_Del_DIV21_581 - Interaction2Control_DIV21_553,
    NRXN1_3PrimeDel_DIV21_641_553 = Interaction23_Prime_Del_DIV21_641 - Interaction2Control_DIV21_553,
    NRXN1_3PrimeDel_DIV35_581_2607 = Interaction23_Prime_Del_DIV35_581 - Interaction2Control_DIV35_2607,
    NRXN1_3PrimeDel_DIV35_641_2607 = Interaction23_Prime_Del_DIV35_641 - Interaction2Control_DIV35_2607,
    NRXN1_3PrimeDel_DIV35_581_553 = Interaction23_Prime_Del_DIV35_581 - Interaction2Control_DIV35_553,
    NRXN1_3PrimeDel_DIV35_641_553 = Interaction23_Prime_Del_DIV35_641 - Interaction2Control_DIV35_553,
    NRXN1_5PrimeDel_DIV14_972_2607 = Interaction25_Prime_Del_DIV14_972 - Interaction2Control_DIV14_2607,
    NRXN1_5PrimeDel_DIV14_973_2607 = Interaction25_Prime_Del_DIV14_973 - Interaction2Control_DIV14_2607,
    NRXN1_5PrimeDel_DIV14_972_553 = Interaction25_Prime_Del_DIV14_972 - Interaction2Control_DIV14_553,
    NRXN1_5PrimeDel_DIV14_973_553 = Interaction25_Prime_Del_DIV14_973 - Interaction2Control_DIV14_553,
    NRXN1_5PrimeDel_DIV21_972_2607 = Interaction25_Prime_Del_DIV21_972 - Interaction2Control_DIV21_2607,
    NRXN1_5PrimeDel_DIV21_973_2607 = Interaction25_Prime_Del_DIV21_973 - Interaction2Control_DIV21_2607,
    NRXN1_5PrimeDel_DIV21_972_553 = Interaction25_Prime_Del_DIV21_972 - Interaction2Control_DIV21_553,
    NRXN1_5PrimeDel_DIV21_973_553 = Interaction25_Prime_Del_DIV21_973 - Interaction2Control_DIV21_553,
    NRXN1_5PrimeDel_DIV35_972_2607 = Interaction25_Prime_Del_DIV35_972 - Interaction2Control_DIV35_2607,
    NRXN1_5PrimeDel_DIV35_973_2607 = Interaction25_Prime_Del_DIV35_973 - Interaction2Control_DIV35_2607,
    NRXN1_5PrimeDel_DIV35_972_553 = Interaction25_Prime_Del_DIV35_972 - Interaction2Control_DIV35_553,
    NRXN1_5PrimeDel_DIV35_973_553 = Interaction25_Prime_Del_DIV35_973 - Interaction2Control_DIV35_553,
    NRXN1_5vs3_DIV14_972_581 = Interaction25_Prime_Del_DIV14_972 - Interaction23_Prime_Del_DIV14_581,
    NRXN1_5vs3_DIV14_972_641 = Interaction25_Prime_Del_DIV14_972 - Interaction23_Prime_Del_DIV14_641,
    NRXN1_5vs3_DIV14_973_581 = Interaction25_Prime_Del_DIV14_973 - Interaction23_Prime_Del_DIV14_581,
    NRXN1_5vs3_DIV14_973_641 = Interaction25_Prime_Del_DIV14_973 - Interaction23_Prime_Del_DIV14_641,
    NRXN1_5vs3_DIV21_972_581 = Interaction25_Prime_Del_DIV21_972 - Interaction23_Prime_Del_DIV21_581,
    NRXN1_5vs3_DIV21_972_641 = Interaction25_Prime_Del_DIV21_972 - Interaction23_Prime_Del_DIV21_641,
    NRXN1_5vs3_DIV21_973_581 = Interaction25_Prime_Del_DIV21_973 - Interaction23_Prime_Del_DIV21_581,
    NRXN1_5vs3_DIV21_973_641 = Interaction25_Prime_Del_DIV21_973 - Interaction23_Prime_Del_DIV21_641,
    NRXN1_5vs3_DIV35_972_581 = Interaction25_Prime_Del_DIV35_972 - Interaction23_Prime_Del_DIV35_581,
    NRXN1_5vs3_DIV35_972_641 = Interaction25_Prime_Del_DIV35_972 - Interaction23_Prime_Del_DIV35_641,
    NRXN1_5vs3_DIV35_973_581 = Interaction25_Prime_Del_DIV35_973 - Interaction23_Prime_Del_DIV35_581,
    NRXN1_5vs3_DIV35_973_641 = Interaction25_Prime_Del_DIV35_973 - Interaction23_Prime_Del_DIV35_641,
    levels = colnames(design))

fit = lmFit( vobj, design, block=info1$Donor, correlation=dupcor$consensus )

res <- residuals( fit, vobj)

fit2 <- contrasts.fit(fit, cons_separate_donor_celltype)
fit2 = eBayes(fit2)

resList2 = list()
for(key in colnames(cons_separate_donor_celltype)){
  #res = topTable(fit2, coef=key, number=Inf, sort.by='P')
  res = topTable(fit2, coef=key, number=Inf, sort.by='none')
  res$qvalue = qvalue( res$P.Value )$qvalue
  res$gene = geneInfo$geneName[match(rownames(res), geneInfo$Geneid)]
  resList2[[key]] = res
}

for( key in names(resList2)){
  file = paste0("NRXN1_bulk_RNA_DE_regressHiPSC_", key, ".csv")
  write.table( data.frame(geneName = rownames(resList2[[key]]), resList2[[key]]), file, quote=FALSE, sep="\t", row.names=FALSE)
}

##Filter DEGs by requiring them also DEGs between each pair of samples from each donor.
fc=1.5
qvalue=0.1
for (key in names(resList)) {
  pairlist=grep(key, names(resList2))
  uplist = resList[[key]]$logFC>=log2(fc) & resList[[key]]$qvalue<qvalue
  #uplist = uplist & resList2[[pairlist[1]]]$logFC>=log2(fc) & resList2[[pairlist[2]]]$logFC>=log2(fc) & resList2[[pairlist[3]]]$logFC>=log2(fc) & resList2[[pairlist[4]]]$logFC>=log2(fc)
  uplist = rownames(resList[[key]])[uplist]
  #uplist

  downlist = resList[[key]]$logFC<=-log2(fc) & resList[[key]]$qvalue<qvalue
  #downlist = downlist & resList2[[pairlist[1]]]$logFC<=-log2(fc) & resList2[[pairlist[2]]]$logFC<=-log2(fc) & resList2[[pairlist[3]]]$logFC<=-log2(fc) & resList2[[pairlist[4]]]$logFC<=-log2(fc)
  downlist = rownames(resList[[key]])[downlist]
  #downlist

  oplist = append(uplist, downlist)
  print(c(key, length(uplist), length(downlist), length(oplist)))
  outlist = resList[[key]][oplist, ]
  file = paste0("NRXN1_bulk_RNA_DE_regressHiPSC_", key, ".fc_op_", fc, ".q_",qvalue,".csv")
  write.table( outlist, file, quote=FALSE, sep=",", row.names=T)
}
