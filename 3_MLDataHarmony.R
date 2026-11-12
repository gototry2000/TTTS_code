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

TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_TTTS_filter_merge.qs")

TTTS_filter_merge <- NormalizeData(TTTS_filter_merge, verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:50, verbose = T)
TTTS_filter_merge@meta.data<-TTTS_filter_merge@meta.data[,1:9]

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_before.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_before.qs")

DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "dataset",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  labs(title = "Before integration-dataset")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_before.pdf",width = 6,height = 5)

####完善metadata的部分
atad_cells <- rownames(TTTS_filter_merge@meta.data[TTTS_filter_merge@meta.data$Source == "GSE162170", ])
new_ident <- sub("^(hft_w\\d+)_.*", "\\1", atad_cells)
head(data.frame(cell = atad_cells, new_ident = new_ident))
TTTS_filter_merge@meta.data[atad_cells, "orig.ident"] <- new_ident
####添加临床信息
Table_S1<-read_csv("/home/iyun1806/MLdata/TableS1.csv")
TTTS_filter_merge$Sample_ID <- TTTS_filter_merge$orig.ident

###先存好metadata的行名，取交集后会失去行名导致找差异基因的时候报错
TTTS_filter_merge$barcode <- rownames(TTTS_filter_merge@meta.data)

TTTS_filter_merge@meta.data <- TTTS_filter_merge@meta.data %>%
  left_join(Table_S1, by = "Sample_ID")
rownames(TTTS_filter_merge@meta.data)<-TTTS_filter_merge$barcode

#####重新运行PCA
TTTS_filter_merge <- NormalizeData(TTTS_filter_merge, verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:50, verbose = T)

####删除所有ENSG开头，基因中带.和-的基因
library(stringr)
keep_genes <- rownames(TTTS_filter_merge)[
  !str_detect(rownames(TTTS_filter_merge), "^ENSG|\\.|-")
]
TTTS_filter_merge <- TTTS_filter_merge[keep_genes, ]

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
#####################Harmony-------------------------------
###0.用默认参数进行合并
TTTS_filter_merge <- TTTS_filter_merge %>% 
  RunHarmony("orig.ident", plot_convergence = T,max_iter=100) %>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)

DimPlot(TTTS_filter_merge,reduction = "umap",raster = F,label = F,group.by = "Age/PCW") + 
  scale_color_viridis_d()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/2_merge_AfterClin.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/2_merge_AfterClin.qs")

#######用这个做后续分析

###########用来进行比较的------------
###1.按照样本用RPCA进行整合----------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
#拆分counts
TTTS_filter_merge[["RNA"]] <- split(TTTS_filter_merge[["RNA"]], f = TTTS_filter_merge$orig.ident)
TTTS_filter_merge <- IntegrateLayers(TTTS_filter_merge, method = RPCAIntegration,
                                     orig.reduction = "pca",
                                     new.reduction = "RPCA",
                                     dims = 1:50, k.anchor = 20)
TTTS_filter_merge <- JoinLayers(TTTS_filter_merge)
TTTS_filter_merge<-RunUMAP(TTTS_filter_merge,reduction = "RPCA",dims = 1:50,verbose = T)

DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "Age/PCW",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
#  labs(title = "After integration-dataset")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_RPCA.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_RPCA.qs")

########2.R版bbknnR（国人自己开发的）--------
#devtools::install_github("ycli1995/bbknnR")
#remotes::install_github("ycli1995/bbknnR")
#install.packages("bbknnR")
library(bbknnR)  #2.0.2
options(future.globals.maxSize = 1e9)

TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
TTTS_filter_merge <- RunBBKNN(TTTS_filter_merge, "orig.ident", method = "nndescent",run_TSNE = F,run_UMAP = T,n_pcs = 50) 
# 自己跑 UMAP，注意输入是 BBKNN 图（默认存为 assay + "_bbknn"，即 "RNA_bbknn"，但是要显式指定“uwot”）
# TTTS_filter_merge <- RunUMAP(TTTS_filter_merge, graph = "RNA_bbknn",umap.method = "uwot",verbose = T)

DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "Age/PCW",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_bbknnR.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_bbknnR.qs")

###或者用R中思路的
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
TTTS_filter_merge <- RunBBKNN(TTTS_filter_merge, "orig.ident", method = "nndescent",run_TSNE = F,run_UMAP = F,n_pcs = 50) 
# 自己跑 UMAP，注意输入是 BBKNN 图（默认存为 assay + "_bbknn"，即 "RNA_bbknn"，但是要显式指定“uwot”）
#但是也会报错，就不用了：Error in RunUMAP.Graph(object = data.use, reduction.model = reduction.model,  : 
#Cannot find UMAP, please install through pip (e.g. pip install umap-learn).

TTTS_filter_merge <- RunUMAP(TTTS_filter_merge, graph = "RNA_bbknn",umap.method = "uwot",verbose = T)

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_bbknnR2.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_bbknnR2.qs")

###################python版本bbknn---------------



###################3.scVI整合---------------
############先转成python版本
##linux中
conda activate scbridge

##R中（在linux中运行吧）
R
#remotes::install_github("xiaoqqjun/scBridge")
library(reticulate)
library(Matrix)
library(scBridge)
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
use_python("/home/iyun1806/miniconda3/envs/scbridge/bin/python")
seurat2py(TTTS_filter_merge, "/home/iyun1806/MLdata/analysis/results/2_harmony/3_TTTS_filter_merge.h5ad")

conda activate scVI
python
#####python中运行
import anndata as ad
import scvi
import numpy as np
import pandas as pd

# ---------- 1. 设置种子（保证可复现）
scvi.settings.seed = 1234          # scVI 全局种子
# 注意：GPU 上由于并行计算，严格逐位复现不保证，但基本趋势可复现

# ---------- 2. 读数据 
adata = ad.read_h5ad("/home/iyun1806/MLdata/analysis/results/2_harmony/3_TTTS_filter_merge.h5ad")
print("obs 列:", adata.obs.columns.tolist())   # 先确认列名
adata.X.max()
# ---------- 3. 关键：scVI 必须用原始计数 
# 如果 X 是 log 归一化数据、counts 在 layers 里，覆盖 X：
if "counts" in adata.layers:
 adata.X = adata.layers["counts"].copy()
print("已用 layers['counts'] 作为 scVI 输入")
# 如果两者都没有（只有 log 数据），只能 expm1 近似，结果可靠性下降

# ---------- 4. 数据登记
scvi.model.SCVI.setup_anndata(
 adata,
 batch_key="orig.ident",                 # 批次 = 测序文库（技术差异）
 categorical_covariate_keys=["Age/PCW"] # 周龄，保留不抹掉
)

# ---------- 5. 建模 + 训练（train会很耗时）
model = scvi.model.SCVI(adata, n_latent=50)
model.train(max_epochs=200, early_stopping=True)

# ---------- 7. 保存 
#因为Age/PCW这个斜杠导致保存都会报错，以后还是要提前改过来
# model.save("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scVI/scvi_model", save_anndata=F, overwrite=True)
# adata.obs = adata.obs.rename(columns={"Age/PCW": "Age_PCW"})
# adata.obsm["X_scvi"] = latent
# adata.write_h5ad("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scVI/1_TTTS_scvi.h5ad")    # 输出 1：h5ad（内含 X_scvi）

# 推荐：再单独导出一份 latent CSV，回 R 最稳
latent = model.get_latent_representation()
print("latent 形状:", latent.shape)          # 应为 细胞数 x 50
#latent 形状: (95682, 50)

pd.DataFrame(
  latent,
  index=adata.obs_names,
  columns=[f"scVI_{i+1}" for i in range(latent.shape[1])]
).to_csv(
  "/home/iyun1806/MLdata/analysis/results/2_harmony/4_scVI/1_TTTS_scvi_latent.csv.gz",
  compression="gzip",
  index_label="barcode"
)
##R 中直接加入单细胞对象中
latent <- read.csv("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scVI/1_TTTS_scvi_latent.csv.gz", row.names = 1, check.names = FALSE)
TTTS_filter_merge[["scvi"]] <- CreateDimReducObject(as.matrix(latent), key = "scvi_", assay = DefaultAssay(TTTS_filter_merge))
TTTS_filter_merge <- RunUMAP(TTTS_filter_merge, reduction = "scvi", dims = 1:ncol(latent),
                             reduction.name = "umap_scvi", reduction.key = "scviUMAP_")
# seu <- FindNeighbors(seu, reduction = "scanvi", dims = 1:ncol(lat))
# seu <- FindClusters(seu, resolution = 0.5)
DimPlot(TTTS_filter_merge,reduction = "scvi",
        group.by = "Age/PCW",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/4_scVI/2_TTTS_scvi_UMAP.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/2_TTTS_scvi_UMAP.qs")














##################4.scANVI-----------------
# 用 scANVI 且 LABELS_KEY=None（全 unknown），功能上等价于"无监督"，此时换成 scVI 会更快；一旦填了真实的注释列，scANVI 的半监督优势才会真正发挥出来。
# 那就用之前的4_MKDataCluster.R中注释的,但是只用小胶和内皮细胞当作确定好的标签
TTTS_filter_merge2<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster/2_TTTS__clusterType.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
identical(colnames(TTTS_filter_merge), colnames(TTTS_filter_merge2))
TTTS_filter_merge$Type <- TTTS_filter_merge2$Type
TTTS_filter_merge@meta.data <- TTTS_filter_merge@meta.data %>%
  mutate(Type = case_when(
    Type %in% c("Microglia", "Vasc.") ~ Type,
    TRUE ~ "unknown"
  ))
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/1_TTTS_forAnvi.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/1_TTTS_forAnvi.qs")

##linux中
conda activate scbridge

##R中（在linux中运行吧）
R
#remotes::install_github("xiaoqqjun/scBridge")
library(reticulate)
library(Matrix)
library(scBridge)
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/1_TTTS_forAnvi.qs")
use_python("/home/iyun1806/miniconda3/envs/scbridge/bin/python")
seurat2py(TTTS_filter_merge, "/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/2_TTTS_forAnvi.h5ad")

conda deactivate
conda activate scVI
#####python中运行
import anndata as ad
import scvi
import numpy as np
import pandas as pd

# ---------- 1. 设置种子（保证可复现）
scvi.settings.seed = 1234          # scVI 全局种子
# 注意：GPU 上由于并行计算，严格逐位复现不保证，但基本趋势可复现

# ---------- 2. 读数据 
adataANVI = ad.read_h5ad("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/2_TTTS_forAnvi.h5ad")
print("obs 列:", adataANVI.obs.columns.tolist())   # 先确认列名

# ---------- 3. 关键：scVI 必须用原始计数 
# 如果 X 是 log 归一化数据、counts 在 layers 里，覆盖 X：
if "counts" in adataANVI.layers:
  adataANVI.X = adataANVI.layers["counts"].copy()
print("已用 layers['counts'] 作为 scVI 输入")
adataANVI.X.max()
# 如果两者都没有（只有 log 数据），只能 expm1 近似，结果可靠性下降

scvi.model.SCANVI.setup_anndata(
  adataANVI, batch_key="orig.ident", labels_key="Type", unlabeled_category="unknown",categorical_covariate_keys=["Age/PCW"]
)
model = scvi.model.SCANVI(adataANVI,n_latent=50)
##这一步耗时，与max_epoches相关
model.train(max_epochs=200,n_samples_per_label=500)

##锁定训练状态
model.eval()
##存储latent space
adataANVI.obsm["X_scANVI"] = model.get_latent_representation()

model.write_h5ad("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/3_TTTS_scANVImodel.h5ad")
adataANVI.write_h5ad("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/3_TTTS_scANVIresults.h5ad")

####之后就是转成seurat格式
## R中,在linux中运行好些？
library(reticulate)
library(Matrix)
library(scBridge)
use_python("/home/iyun1806/miniconda3/envs/scbridge/bin/python")
TTTS_filter_merge <- py2seurat("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/3_TTTS_scANVIresults.h5ad")

# 把潜空间设为降维对象，后续 FindNeighbors/UMAP/聚类都用它
TTTS_filter_merge[["scanvi"]] <- CreateDimReducObject(
  t(reducedDim(sce, "X_scANVI")),        # 注意转置：细胞 × 特征
  key = "scanvi_", assay = DefaultAssay(TTTS_filter_merge))
DefaultAssay(TTTS_filter_merge)
DefaultAssay(TTTS_filter_merge) <- "RNA"   # 确认 assay 名称，可能是 "X"

DefaultAssay(TTTS_filter_merge) <- RunUMAP(DefaultAssay(TTTS_filter_merge), reduction = "scanvi", dims = 1:ncol(seu[["scanvi"]]),
               reduction.name = "umap_scanvi", reduction.key = "scanviUMAP_")
# seu <- FindNeighbors(seu, reduction = "scanvi", dims = 1:ncol(seu[["scanvi"]]))
# seu <- FindClusters(seu, resolution = 0.5)


###这些是在python中用到的
# sc.pp.neighbors(adataANVI, use_rep="X_scANVI")
# sc.tl.leiden(adataANVI, resolution=RES, key_added="leiden_scANVI")
# sc.tl.umap(adataANVI)
# 
# adataANVI.write_h5ad(os.path.join(OUT_DIR, "scanvi_integrated.h5ad"))
# model.save(os.path.join(OUT_DIR, "scanvi_model"), overwrite=True)
# pd.DataFrame(
#   adataANVI.obsm["X_scANVI"], index=adataANVI.obs_names,
#   columns=[f"scANVI_{i + 1}" for i in range(adataANVI.obsm["X_scANVI"].shape[1])]
# ).to_csv(os.path.join(OUT_DIR, "scanvi_latent.csv"), index_label="barcode")

##########或者只保存潜空间成csv，直接加到原seurat对象中，这样更方便些
import pandas as pd
latent = model.get_latent_representation(adataANVI)

pd.DataFrame(
  latent,
  index=adataANVI.obs_names,
  columns=[f"scANVI_{i+1}" for i in range(latent.shape[1])]
).to_csv(
  "/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/3_TTTS_scANVIresults_latent.csv.gz",
  compression="gzip",
  index_label="barcode"
)

##R 中
lat <- read.csv("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/3_TTTS_scANVIresults_latent.csv.gz", row.names = 1, check.names = FALSE)
TTTS_filter_merge[["scanvi"]] <- CreateDimReducObject(as.matrix(lat), key = "scanvi_", assay = DefaultAssay(TTTS_filter_merge))
TTTS_filter_merge <- RunUMAP(TTTS_filter_merge, reduction = "scanvi", dims = 1:ncol(lat),
               reduction.name = "umap_scanvi", reduction.key = "scanviUMAP_")
# seu <- FindNeighbors(seu, reduction = "scanvi", dims = 1:ncol(lat))
# seu <- FindClusters(seu, resolution = 0.5)

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/4_TTTS_scanvi_UMAP.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/4_TTTS_scanvi_UMAP.qs")


##################5.LIGER-----------------
#library(devtools)
##临时把路径换到4.4.2版本，避免又因为版本依赖安装失败
#Sys.setenv(PATH = paste0("/home/iyun1806/tools/cmake-4.4.2-linux-x86_64/bin:", Sys.getenv("PATH")))
#Sys.which("cmake") 
# cmake 
# "/home/iyun1806/tools/cmake-4.4.2-linux-x86_64/bin/cmake" 
#install.packages("RcppPlanc")
#devtools::install_local("/home/iyun1806/Rpackages2/liger-2.2.1.tar.gz")
library(rliger)
library(cowplot)
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")

lig <- rliger::createLiger(Seurat::SplitObject(TTTS_filter_merge, split.by = "orig.ident"),organism = "human")
lig <- rliger::normalize(lig)    %>% 
  selectGenes(nGenes = 5000)  %>% 
  scaleNotCenter() %>% 
  runIntegration(k = 50, lambda = 5)  %>% 
  alignFactors(method = "centroidAlign")   %>% 
  runUMAP(seed = 1234)    
qsave(lig,"/home/iyun1806/MLdata/analysis/results/2_harmony/4_Liger/1_TTTS_liger.qs")

##写下细胞x因子的载荷矩阵,后续定量比较要用（也可以存到原来的数据重新画UMAP）
Hnorm <- tryCatch(rliger::getMatrix(lig, "H.norm"), error = function(e) lig@H.norm)
write.csv(as.matrix(Hnorm), file.path("/home/iyun1806/MLdata/analysis/results/2_harmony/4_Liger", "1_H_norm.csv"), row.names = TRUE)

####在R中读取Hn加到原seurat对象中
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
hn <- read.csv("/home/iyun1806/MLdata/analysis/results/2_harmony/4_Liger/1_H_norm.csv", row.names = 1, check.names = FALSE)
identical(rownames(hn), colnames(TTTS_filter_merge))
#cells <- intersect(rownames(hn), colnames(TTTS_filter_merge))
#hn <- as.matrix(hn[cells, ])
TTTS_filter_merge[["liger"]] <- CreateDimReducObject(as.matrix(hn), key = "liger_", assay = DefaultAssay(TTTS_filter_merge))
TTTS_filter_merge <- RunUMAP(TTTS_filter_merge, reduction = "liger", dims = 1:ncol(hn),
               reduction.name = "umap_liger", reduction.key = "ligerUMAP_",
               verbose = T)
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/4_Liger/2_TTTS_ligerUmap.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_Liger/2_TTTS_ligerUmap.qs")
DimPlot(TTTS_filter_merge,reduction = "umap_liger",
        group.by = "Age/PCW",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')


##################6.fastMNN-----------------
library(SeuratWrappers)
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")

TTTS_filter_merge <- RunFastMNN(object.list = SplitObject(TTTS_filter_merge ,split.by = "orig.ident"),features = 5000,d=50,k=20) %>% 
  RunUMAP(reduction = "mnn", dims = 1:50)
DimPlot(TTTS_filter_merge,reduction = "mnn",
        group.by = "Age.PCW",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/4_FastMNN/1_TTTS_FastMNN.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_FastMNN/1_TTTS_FastMNN.qs")


##################定量比较聚类效果-----------------















###################后续的不再考虑------------------

###1.用组合进行整合
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_before.qs")
TTTS_filter_merge <- TTTS_filter_merge %>% 
  RunHarmony(c("Source","AssayType"), plot_convergence = F,theta=c(4,4),max_iter=100) %>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_after.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_after.qs")
###按照样本进行整合
atad_cells <- rownames(TTTS_filter_merge@meta.data[TTTS_filter_merge@meta.data$Source == "GSE162170", ])
new_ident <- sub("^(hft_w\\d+)_.*", "\\1", atad_cells)
head(data.frame(cell = atad_cells, new_ident = new_ident))
TTTS_filter_merge@meta.data[atad_cells, "orig.ident"] <- new_ident

TTTS_filter_merge <- TTTS_filter_merge %>% 
  RunHarmony(c("orig.ident"), plot_convergence = F) %>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)

DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "dataset",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  labs(title = "After integration-dataset")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_after.pdf",width = 6,height = 5)



###按照数据来源用RPCA进行整合
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_before.qs")
#修改orig.ident
atad_cells <- rownames(TTTS_filter_merge@meta.data[TTTS_filter_merge@meta.data$Source == "GSE162170", ])
new_ident <- sub("^(hft_w\\d+)_.*", "\\1", atad_cells)
head(data.frame(cell = atad_cells, new_ident = new_ident))
TTTS_filter_merge@meta.data[atad_cells, "orig.ident"] <- new_ident
#拆分counts
TTTS_filter_merge[["RNA"]] <- split(TTTS_filter_merge[["RNA"]], f = TTTS_filter_merge$dataset)
TTTS_filter_merge <- IntegrateLayers(TTTS_filter_merge, method = RPCAIntegration,
                                     orig.reduction = "pca",
                                     new.reduction = "RPCA",
                                     dims = 1:50, k.anchor = 5)
TTTS_filter_merge <- JoinLayers(TTTS_filter_merge)
TTTS_filter_merge<-RunUMAP(TTTS_filter_merge,reduction = "RPCA",dims = 1:50,verbose = T)
DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "dataset",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  labs(title = "After integration-dataset")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_RPCA.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_RPCA.qs")

###用SCT和harmony进行整合
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_before.qs")
#修改orig.ident
atad_cells <- rownames(TTTS_filter_merge@meta.data[TTTS_filter_merge@meta.data$Source == "GSE162170", ])
new_ident <- sub("^(hft_w\\d+)_.*", "\\1", atad_cells)
head(data.frame(cell = atad_cells, new_ident = new_ident))
TTTS_filter_merge@meta.data[atad_cells, "orig.ident"] <- new_ident

DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "orig.ident",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  labs(title = "Before integration-sample")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony/1_SCTmergebefore.pdf",width = 8,height = 5)

TTTS_filter_merge <- SCTransform(TTTS_filter_merge, vars.to.regress = "percent.mt", verbose = T)

TTTS_filter_merge <- TTTS_filter_merge %>% 
  RunHarmony(c("orig.ident","Source","AssayType"), plot_convergence = F,theta=c(2,4,4),max_iter=100) %>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)
DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "orig.ident",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  labs(title = "After integration-sample")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony/1_SCTmergeafter.pdf",width = 8,height = 5)

DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "dataset",raster = F,label = F) + 
  scale_color_d3(palette = "category20")+
  labs(title = "After integration-dataset")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_SCT.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_SCT.qs")
