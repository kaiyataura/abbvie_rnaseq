## 2024-03-05
## Tumor explants (PMED) RNAseq
## Analysis on TPM (protein-coding genes only)

############ SET UP ############

library(edgeR)
library(limma)
library(genefilter)
library(ggbiplot)
library(pheatmap)
library(RColorBrewer)
library(xCell)
library(dplyr)
library(tidyr)

# input.dir<-("/Users/lohf/Documents/R_scripts/2024_intern_explant_project/data_files_for_Kaiya")
# setwd("/Users/lohf/Documents/R_scripts/2024_intern_explant_project/data_files_for_Kaiya")

# Load data:
tpm.expset<-readRDS("tpm_expset.Rds")
tpm.expset.hgnc <- readRDS("tpm_expset_hgnc.Rds")
rna.meta<-read.csv("sample_metadata.csv",row.names=1)


###################
#### PCA using TPM object ####
##################
rna.exprs<-log2(exprs(tpm.expset)+1)
rna.exprs<-rna.exprs[fData(tpm.expset)$gene_biotype == "protein_coding",] # Focus only on protein-coding genes
rna.exprs<-rna.exprs[which(rowVars(rna.exprs)>0),] # remove invariable genes (variance = 0)
vars.rna<-rowVars(rna.exprs)
rna.exprs<-rna.exprs[vars.rna >= quantile(vars.rna,0.95),] # top 5% variable genes

Indication = tpm.expset$indication

prcomp.pca<-prcomp(t(rna.exprs), center = TRUE, scale. = TRUE)
pca.plot<-ggbiplot(prcomp.pca, var.axes=FALSE, choices= c(1,2),
                   obs.scale = 1, var.scale = 1,alpha = 0) +
  geom_point(aes(color=Indication),size = 2.5)+
  theme_classic() # I ended up combining ggbiplot and ggplot options for more flexibility
pca.plot

##############
#### Target expression ####
###############
genes <- c("_____","MET","_____","_____")

genes.int<-tpm.expset.hgnc[match(genes,rownames(tpm.expset.hgnc)),]
genes.int<-exprs(genes.int) %>% t() %>% cbind(rna.meta$indication) %>% data.frame()
colnames(genes.int)<-c(genes,"indication")
genes.int <- genes.int %>% mutate_at(genes, ~log2(as.numeric(.)+1)) # convert gene columns to numeric and to Log2 TPM+1


genes.int %>%
  tibble::rownames_to_column("sample_id") %>%
  pivot_longer(all_of(genes)) %>%
  ggplot(.,aes(x=indication,y=value)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1, aes(color=indication))+
  facet_wrap(~name, scales = "free") +
  ylab("log2 TPM")+
  theme_bw(base_size=15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


###############
#### xCELL ####
###############
exprs.Tmp<-exprs(tpm.expset.hgnc)
xcell.score<-xCellAnalysis(data.frame(exprs.Tmp))
cells.int<-c("CD4+ T-cells","CD8+ T-cells","Tregs","B-cells","NK cells","cDC","Macrophages","Monocytes","Neutrophils",
             "Epithelial cells","Fibroblasts","Endothelial cells","Smooth muscle","Pericytes","Platelets","Erythrocytes",
             "Hepatocytes","Mesangial cells","ImmuneScore","StromaScore","MicroenvironmentScore")
xcell.score.int<-xcell.score[cells.int,] # ordered

# Subset colside annotations
annotsSubset <- rna.meta$indication[match(rownames(rna.meta),colnames(xcell.score))] %>% data.frame(indication=.)
rownames(annotsSubset) <- colnames(xcell.score)

IndPalette <- brewer.pal(10, "Paired")
names(IndPalette) <- unique(rna.meta$indication)
IndPalette<-list("indication"=IndPalette)

pheatmap(xcell.score.int,annotation_col = annotsSubset,
         annotation_colors = IndPalette, scale="none",cluster_rows =T, cluster_cols=T,
         clustering_method = "ward.D2")

