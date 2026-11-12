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

#####读取数据--------------
TTTS_list<-qread("/home/iyun1806/MLdata/analysis/results/1_dataread_list.qs")
ga_list<-qread("/home/iyun1806/MLdata/analysis/results/1_dataread_galist.qs")
GSE162170_list<-qread("/home/iyun1806/MLdata/analysis/results/1_dataread_GSE162170_list2.qs")

TTTS_useList<-c(TTTS_list,ga_list,GSE162170_list)
sum(sapply(TTTS_useList,function(x){
  ncol(x)
}))
#[1] 117408
#####检测双细胞--------------
library(scDblFinder)
library(decontX)
library(SingleCellExperiment)
library(BiocParallel)

TTTS_scDblfinderList<-list()
for (i in 1:length(TTTS_useList)) {
  filter_tmp<-as.SingleCellExperiment(TTTS_useList[[i]])
  
  TTTS_scDblfinderList[[i]]<-scDblFinder(filter_tmp,BPPARAM = MulticoreParam(10))
  TTTS_useList[[i]]$scDbfinder<-TTTS_scDblfinderList[[i]]$scDblFinder.class
}

for (i in 1:length(TTTS_useList)) {
  print(table(TTTS_useList[[i]]$scDbfinder))
}
for (i in 1:length(TTTS_useList)) {
  TTTS_useList[[i]]<-subset(TTTS_useList[[i]],subset=scDbfinder=="singlet")
}
sum(sapply(TTTS_useList,function(x){
  ncol(x)
}))
#[1] 100576
qsave(TTTS_useList,"/home/iyun1806/MLdata/analysis/results/2_TTTS_scDbDecontX.qs")
TTTS_useList<-qread("/home/iyun1806/MLdata/analysis/results/2_TTTS_scDbDecontX.qs")

#####MAD过滤---------------
#过滤低质量细胞(同时sn计算线粒体百分比会出现100和NaN，所以记为0)
for (i in 1:length(TTTS_useList)){
  TTTS_useList[[i]][["percent.mt"]]<-PercentageFeatureSet(TTTS_useList[[i]], pattern = "^MT-")
  TTTS_useList[[i]]$percent.mt[is.nan(TTTS_useList[[i]]$percent.mt)] <- 0
}

plot_list<-list()
for (i in 1:length(TTTS_useList)){
  plot_list[[i]]<-VlnPlot(TTTS_useList[[i]], 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                          pt.size = 0, 
                          ncol = 3) +
    plot_annotation(title = names(TTTS_useList)[[i]])
}

#MAD自适应函数计算
filter_seurat2 <- function(obj, nmads = 3) {
  # nFeature上下限
  nFeature_med <- median(obj$nFeature_RNA)
  nFeature_mad <- mad(obj$nFeature_RNA)
  nFeature_low <- max(200, nFeature_med - nmads * nFeature_mad)
  nFeature_high <- nFeature_med + nmads * nFeature_mad
  
  # nCount上限
  nCount_med <- median(obj$nCount_RNA)
  nCount_mad <- mad(obj$nCount_RNA)
  nCount_high <- nCount_med + nmads * nCount_mad
  
  
  obj <- subset(obj, subset = nFeature_RNA > nFeature_low & 
                  nFeature_RNA < nFeature_high & 
                  nCount_RNA < nCount_high &
                  percent.mt <= 5)
  return(obj)
}

# 记录过滤前后细胞数
filter_log <- data.frame(sample = names(TTTS_useList), 
                         before = NA, after = NA)
all_filter_list<-list()
for (i in 1:length(TTTS_useList)) {
  filter_log$before[i] <- ncol(TTTS_useList[[i]])
  all_filter_list[[i]] <- filter_seurat2(TTTS_useList[[i]])
  filter_log$after[i] <- ncol(all_filter_list[[i]])
}

filter_log$removed <- filter_log$before - filter_log$after
filter_log$pct_removed <- round(filter_log$removed / filter_log$before * 100, 1)
qsave(filter_log,"/home/iyun1806/MLdata/analysis/results/2_filter_log.txt")
filter_log<-qread("/home/iyun1806/MLdata/analysis/results/2_filter_log.txt")
filter_log$sample<-c("D1","D2","R1","R2","ga1","ga2","ga3","hft1","hft2","hft3","hft4")

sum(sapply(all_filter_list, ncol))
#[1] 95682
TTTS_filterList<-all_filter_list

for (i in 1:length(TTTS_filterList)){
  TTTS_filterList[[i]]$dataset<-paste0(TTTS_filterList[[i]]$Source,"-",TTTS_filterList[[i]]$AssayType)
}
names(TTTS_filterList)<-names(TTTS_useList)

qsave(TTTS_filterList,"/home/iyun1806/MLdata/analysis/results/2_TTTSAllFilterList.qs")
TTTS_filterList<-qread("/home/iyun1806/MLdata/analysis/results/2_TTTSAllFilterList.qs")

########删除基因------------------
remove_patterns <- c(
  "^MT-",          # 线粒体基因
  "^RP[SL]\\d",    # 核糖体蛋白基因 (RPS, RPL)
  "^MRP[SL]",      # 线粒体核糖体基因
  "^LINC",         # 长链非编码RNA (LINC)
  "^AC\\d{6}",     # Ensembl自动注释的非编码基因 (AC开头)
  "^AL\\d{6}",     # Ensembl自动注释的非编码基因 (AL开头)
  "^AP\\d{6}",     # Ensembl自动注释的非编码基因 (AP开头)
  "^BX\\d{6}",     # Ensembl自动注释的非编码基因 (BX开头)
  "^CR\\d{6}",     # Ensembl自动注释的非编码基因 (CR开头)
  "^CT\\d{6}",     # Ensembl自动注释的非编码基因 (CT开头)
  "^FP\\d{6}",     # Ensembl自动注释的非编码基因 (FP开头)
  "^AF\\d{6}",     # Ensembl自动注释的非编码基因 (AF开头)
  "^KB\\d{6}",     # 非编码基因
  "-AS\\d+$",      # 反义转录本 (如OTX2-AS1, GMDS-AS1)
  "-DT$",          # 发散转录本
  "-IT\\d*$",      # 内含子转录本
  "^SNORD",        # 小核仁RNA
  "^SNHG",         # 小核仁RNA宿主基因
  "^RNU\\d",       # 小核RNA
  "^SCARNA",       # 小Cajal体特异性RNA
  "^snoU",         # snoRNA
  "^MIR\\d",       # microRNA
  "^MTRNR",        # 线粒体rRNA
  "^MT-",          # 线粒体基因（重复确保）
  "^MALAT1$",      # 高表达非编码RNA，常为技术噪音
  "^NEAT1$",       # 核仁富集非编码RNA
  "^XIST$"         # X染色体失活基因（根据需要保留或去除）
)

# 合并所有pattern
remove_regex <- paste(remove_patterns, collapse = "|")

for (i in 1:length(TTTS_filterList)) {
  gene_all<-rownames(TTTS_filterList[[i]])
  genes_remove <- grep(remove_regex, gene_all, value = TRUE)
  genes_keep <- setdiff(gene_all, genes_remove)
  TTTS_filterList[[i]] <- subset(TTTS_filterList[[i]], features = genes_keep)
}
for (i in 1:length(TTTS_filterList)) {
  cat("过滤后每个对象的基因数:", nrow(TTTS_filterList[[i]]), "\n")
}
# 过滤后每个对象的基因数: 25676 
# 过滤后每个对象的基因数: 25984 
# 过滤后每个对象的基因数: 26511 
# 过滤后每个对象的基因数: 26308 
# 过滤后每个对象的基因数: 22948 
# 过滤后每个对象的基因数: 23034 
# 过滤后每个对象的基因数: 21207 
# 过滤后每个对象的基因数: 17990 
# 过滤后每个对象的基因数: 17990 
# 过滤后每个对象的基因数: 17990 
# 过滤后每个对象的基因数: 17990 

qsave(TTTS_filterList,"/home/iyun1806/MLdata/analysis/results/2_TTTSAllFilterList_geneClean.qs")
TTTS_filterList<-qread("/home/iyun1806/MLdata/analysis/results/2_TTTSAllFilterList_geneClean.qs")

################merge数据------------
TTTS_filter_merge <- merge(
  TTTS_filterList[[1]],
  y = TTTS_filterList[-1],
  add.cell.ids = names(TTTS_filterList)
) %>%
  JoinLayers()
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/3_TTTS_filter_merge.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_TTTS_filter_merge.qs")



