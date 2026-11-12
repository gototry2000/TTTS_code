set.seed(1234)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
#library(ggstatsplot)
library(patchwork)
library(tidyverse)
library(Seurat)
library(RSpectra)
library(harmony)
library(qs)
library(clustree)
library(readxl)
library(data.table)

#####自测数据----------------------------
base_dir <- "/home/iyun1806/MLdata/rawdata"
sample_folders <- c("D1","D2","R1","R2")
sample_dirs <- file.path(base_dir, sample_folders)

seurat_list<-list()
for (i in 1:length(sample_folders)){
  seurat_data<-Read10X(data.dir = sample_dirs[i])
  seurat_list[[i]]<-CreateSeuratObject(counts = seurat_data,
                                       project=sample_folders[[i]],
                                       min.features=200,
                                       min.cells=3)
  seurat_list[[i]]$ID<-sample_folders[[i]]
}
names(seurat_list)<-sample_folders
TTTS_list<-seurat_list
for (i in 1:length(TTTS_list)) {
  TTTS_list[[i]]$Source<-"own"
  TTTS_list[[i]]$AssayType<-"sn"
}
qsave(TTTS_list,"/home/iyun1806/MLdata/analysis/results/1_dataread_list.qs")
TTTS_list<-qread("/home/iyun1806/MLdata/analysis/results/1_dataread_list.qs")

sum(sapply(TTTS_list,function(x){
  ncol(x)
}))
# [1] 58521


######GSE168408的前额叶部分（3个）-------------------------
ga22 <- Read10X_h5(file = "/home/iyun1806/MLdata/GSE168408_RAW/GSM5138509_RL2103_ga22_snRNAseq_filtered_feature_bc_matrix.h5")
ga22 <- CreateSeuratObject(counts = ga22,
                                 project = "ga22",
                                 min.features = 200,
                                 min.cells = 3)
ga24 <- Read10X_h5(file = "/home/iyun1806/MLdata/GSE168408_RAW/GSM5138511_RL2107_ga24_snRNAseq_filtered_feature_bc_matrix.h5")
ga24 <- CreateSeuratObject(counts = ga24,
                           project = "ga24",
                           min.features = 200,
                           min.cells = 3)
ga34 <- Read10X_h5(file = "/home/iyun1806/MLdata/GSE168408_RAW/GSM5138513_RL2121_ga34_snRNAseq_filtered_feature_bc_matrix.h5")
ga34 <- CreateSeuratObject(counts = ga34,
                           project = "ga34",
                           min.features = 200,
                           min.cells = 3)
ga_list<-list(ga22,ga24,ga34)
names(ga_list)<-c("ga22","ga24","ga34")
for (i in 1:length(ga_list)) {
  ga_list[[i]]$Source<-"GSE168408"
  ga_list[[i]]$AssayType<-"sn"
}
qsave(ga_list,"/home/iyun1806/MLdata/analysis/results/1_dataread_galist.qs")
ga_list<-qread("/home/iyun1806/MLdata/analysis/results/1_dataread_galist.qs")
sum(sapply(ga_list,function(x){
  ncol(x)
}))
#[1] 30900

######GSE162170的部分-------------------------
expr_data_sn <- as.data.frame(fread("/home/iyun1806/MLdata/GSE162170/GSE162170_rna_counts.tsv.gz"))
library(clusterProfiler)
library(org.Hs.eg.db)
result <- bitr(expr_data_sn$V1,
               fromType = "ENSEMBL",
               toType = "SYMBOL",
               OrgDb = org.Hs.eg.db)

expr_data <- merge(result, expr_data_sn,
                     by.x = "ENSEMBL",
                     by.y = "V1")

# 按行均值保留表达量最高的
expr_data$mean_expr <- rowMeans(expr_data[, -c(1, 2)], na.rm = TRUE)
expr_data <- expr_data[order(expr_data$SYMBOL, -expr_data$mean_expr), ]
expr_data <- expr_data[!duplicated(expr_data$SYMBOL), ]
expr_data$mean_expr <- NULL

rownames(expr_data)<-expr_data$SYMBOL
expr_data<-expr_data[,-c(1,2)]

cell_matadata<-fread("/home/iyun1806/MLdata/GSE162170/GSE162170_rna_cell_metadata.txt.gz")
rownames(cell_matadata)<-cell_matadata$Cell.ID
cell_matadata<-cell_matadata[,-1]
# hft_w16_p7_r1 hft_w16_p7_r2 hft_w20_p3_r1 hft_w20_p3_r2 hft_w21_p5_r1 hft_w21_p5_r2 hft_w24_p6_r1 hft_w24_p6_r2 
# 9294         10079          6393          6174          6169          6388          6131          7240 
GSE162170 <- CreateSeuratObject(counts = expr_data,
                                meta.data = cell_matadata,
                                project = "GSE162170",
                                min.features = 200,
                                min.cells = 3)
qsave(GSE162170,"/home/iyun1806/MLdata/analysis/results/1_dataread_GSE162170.qs")
GSE162170<-qread("/home/iyun1806/MLdata/analysis/results/1_dataread_GSE162170.qs")

GSE162170_list <- SplitObject(GSE162170, split.by = "Sample.ID")
qsave(GSE162170_list,"/home/iyun1806/MLdata/analysis/results/1_dataread_GSE162170_list.qs")
GSE162170_list<-qread("/home/iyun1806/MLdata/analysis/results/1_dataread_GSE162170_list.qs")

GSE162170_list<-GSE162170_list[c(1,3,5,7)]
for (i in 1:length(GSE162170_list)) {
  GSE162170_list[[i]]$Source<-"GSE162170"
  GSE162170_list[[i]]$AssayType<-"sc"
}
qsave(GSE162170_list,"/home/iyun1806/MLdata/analysis/results/1_dataread_GSE162170_list2.qs")
GSE162170_list<-qread("/home/iyun1806/MLdata/analysis/results/1_dataread_GSE162170_list2.qs")
sum(sapply(GSE162170_list,function(x){
  ncol(x)
}))
# [1] 27987





