# code to calculate gene set enrichment scores from a gene expression matrix
library(GSVA)
library(GSEABase)
library(org.Hs.eg.db)
library(dplyr)
library(RColorBrewer)

# Put your path to geneset RDS object here and load it in:
gs <- readRDS("bindea_genesets.rds")

# Next, load in your gene expression object.
# "tpm.expset.hgnc" object is basically the same data as your previous "tpm_expset", but 
# has HGNC gene symbols rather than Ensembl ID's as the gene identifier (ie. rownames)
tpm.expset.hgnc <- readRDS("tpm_expset_hgnc.Rds")
head(rownames(tpm.expset.hgnc))
# It's because our geneset object is also written in gene symbols.
# In an ideal scenario, our gene expression matrix stays the same and we just convert the ID's in the geneset object.

# Run ssgsea to calculate single sample GSEA scores
param <- ssgseaParam(tpm.expset.hgnc, gs, minSize=10, maxSize=500)
sgsea <- gsva(param, verbose=T)
# The output looks like an expression matrix, but instead of genes, it has gene scores as rows.

indications <- factor(phenoData(sgsea)$indication)



palette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(palette)
colors <- c('red','orange','yellow','limegreen','green','lightblue','blue','purple','magenta','pink')
col.cell <- colors[indications]
heatmap.2(exprs(sgsea), col=rev(morecols(50)), trace="none", main="Gene Enrichment Scores", ColSideColors=col.cell, scale="row")


pca <- prcomp(t(exprs(sgsea)))
ggbiplot(pca, var.axes=F, groups=indications) +
  labs(color="Indications")
