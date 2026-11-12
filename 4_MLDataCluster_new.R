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

TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/2_merge_AfterClin.qs")

TTTS_filter_merge <- TTTS_filter_merge %>%  
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution=seq(0,1,0.05))
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/1_TTTS_merge_cluster.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/1_TTTS_merge_cluster.qs")

clustree(TTTS_filter_merge,prefix='RNA_snn_res.')
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/1_clustree.pdf",width = 15,height = 15)

DimPlot(TTTS_filter_merge,reduction = "umap",raster = F,label = T) + 
  #  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')+NoLegend()

####看marker的表达
##大类的表达
FeaturePlot(TTTS_filter_merge,features = c("RBFOX3","SNAP25","SYT1","STMN2","ENO2",#泛神经元表达
                                           "SLC17A7","SLC17A6",#兴奋性神经元表达
                                           "GAD1","GAD2","SLC32A1",#抑制性神经元表达
                                           "PTPRC","OLIG1","PECAM1","PDGFRB","COL1A1"))
VlnPlot(TTTS_filter_merge,features = c("RBFOX3","SNAP25","SYT1","STMN2","ENO2",#泛神经元表达
                                           "SLC17A7","SLC17A6",#兴奋性神经元表达
                                           "GAD1","GAD2","SLC32A1",#抑制性神经元表达
                                           "PTPRC","OLIG1","PECAM1","PDGFRB","COL1A1"),pt.size = 0)
FeaturePlot(TTTS_filter_merge,features = c("GFAP","SLC1A3","SOX9","HEXB","FGFR3",#星胶前体
                                           "AQP4","ALDH1L1","S100B",#成熟星胶
                                           "OLIG1","MBP","PLP1","MOG",#少突
                                           "PDGFRA","CSPG4","OLIG2","PTPRZ1"))#OPC
###30是内皮细胞+壁细胞；27是免疫细胞；19，14，4，29，7，21是抑制性；23是OPC；13，15是星形胶质细胞；

VlnPlot(TTTS_filter_merge,features = c("”SOX9","NFIA","NFIB","GFAP"),pt.size = 0)
###23群的细胞有OPC的表达
subset_23<-subset(TTTS_filter_merge,idents = "23")
subset_23 <- subset_23 %>% NormalizeData(verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>%
  RunHarmony("orig.ident",plot_convergence = F) %>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50)%>% 
  FindClusters(resolution=seq(0,1,0.1))

FeaturePlot(subset_23,features = c("MOG","MBP","PLP1",#少突胶质细胞
                                           "PDGFRA","CSPG4","OLIG2","PTPRZ1"))#OPC
VlnPlot(subset_23,features = c("MOG","MBP","PLP1",#少突胶质细胞
                               "PDGFRA","CSPG4","OLIG2","PTPRZ1"))#OPC)


subset_23<-RenameIdents(subset_23,"4"="37","7"="37","9"="37",
                        "0"="38","1"="38","2"="38","3"="38","5"="38","6"="38","8"="38")
Idents(TTTS_filter_merge, cells = colnames(subset_23)) <- Idents(subset_23)

###30是内皮细胞+壁细胞；27是免疫细胞；19，14，4，29，7，21是抑制性；23是OPC；13是星形胶质前体细胞，15是RG细胞；37是ODC，38是OPC；20是IPC_pro（增殖型）

##求差异基因
DEGs_TTTS <-FindAllMarkers(TTTS_filter_merge,logfc.threshold=0.25,min.pct=0.25,only.pos = T)
DEGs_TTTS <- DEGs_TTTS %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(DEGs_TTTS,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/2_DEGs_TTTS.csv",row.names = T)
DEGs_TTTS<-read.csv("/home/iyun1806/MLdata/analysis/results/3_cluster_new/2_DEGs_TTTS.csv")

TTTS_filter_merge$num_cluster<-Idents(TTTS_filter_merge)
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/1_TTTS_merge_cluster2.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/1_TTTS_merge_cluster2.qs")


#########################大类命名----------------------------
################1.命名-----
TTTS_filter_merge<-RenameIdents(TTTS_filter_merge,
                                "27"="Microglia","38"="OPC","37"="NOL","13"="Astrocyte","15"="NPC","20"="NPC",
                                "17"="ExN","1"="ExN","18"="ExN","33"="ExN","36"="ExN","2"="ExN","8"="ExN","0"="ExN","5"="ExN","6"="ExN","32"="ExN","11"="ExN","26"="ExN","12"="ExN","24"="ExN","25"="ExN","31"="ExN","35"="ExN","16"="ExN","10"="ExN","3"="ExN","22"="ExN",
                                "14"="InN","19"="InN","7"="InN","21"="InN","4"="InN","29"="InN",
                                "30"="Vascular",
                                "9"="Poor quality","34"="Poor quality","28"="Poor quality")
TTTS_filter_merge$broad_cluster<-Idents(TTTS_filter_merge)
levels(TTTS_filter_merge)<-c("ExN","InN","NPC","Microglia","OPC","NOL","Astrocyte","Vascular","Poor quality")
col_broad<-c("#324376","#6b8b23","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
names(col_broad)<-levels(TTTS_filter_merge)
DimPlot(TTTS_filter_merge,reduction = "umap",raster = F,label = F,cols = col_broad) + 
  #  scale_color_d3(palette = "category20")+
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
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_broad_cluster_Dimplot.pdf",width = 6.3,height = 4.5)

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_TTTS_merge_named.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_TTTS_merge_named.qs")


###########2.画harmony之前的umap图----------
tmp_TTTS <- RunUMAP(TTTS_filter_merge, reduction = "pca", dims = 1:50,
               reduction.name = "umap_before_harmony",
               reduction.key = "UMAPb_")
DimPlot(tmp_TTTS,reduction = "umap_before_harmony",raster = F,label = F,cols = col_broad) + 
  #  scale_color_d3(palette = "category20")+
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
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_broad_cluster_Dimplot_beforeHarmony.pdf",width = 6.3,height = 4.5)

qsave(tmp_TTTS,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_TTTS_merge_namedtmp.qs")
tmp_TTTS<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_TTTS_merge_namedtmp.qs")

########2.展示基因----------

############2.1 Dotplot---------
TTTS_filt <- subset(TTTS_filter_merge, idents = "Poor quality",invert=T)
genes <- c(
  "SLC17A7","NEUROD6","CUX2","RORB","TLE4","BCL11B",
  "GAD1","GAD2","DLX2","ADARB2","LHX6","SST",
  "HOPX","FABP7","EOMES","MKI67",
  "GFAP","SLC1A3","TNC","SLCO1C1","LRIG1",
  "PDGFRA","OLIG1","OLIG2","SOX10",
  "BCAS1","GPR17","ENPP6","TCF7L2","MBP",
  "P2RY12","C1QA","C1QB","CX3CR1","FOLR2",
  "CLDN5","FLT1","PDGFRB","RGS5","COL1A2","COL3A1"
)

DotPlot(
  TTTS_filt,
  features = genes,
  dot.scale = 10
) +
  theme_bw() +
  scale_color_gradient2(
    low = "#32577e",
    mid = "#F6D6B8",
    high = "#9E0015",
    midpoint = 0
  ) +
  scale_size_continuous(
    range = c(1, 10),          # 控制点大小范围
    name = "% expr",           # legend title
    guide = guide_legend(
      title.size = 20,
      label.size = 18
    )
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 20,
      color = "black"
    ),
    axis.text.y = element_text(size = 20, colour = "black"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  labs(x = "", y = "") +
  guides(
    color = guide_colorbar(
      title = "scaled expr",
      title.size = 20,
      labels.size = 18
    )
  )
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_broad_cluster_Dotplot.pdf",width = 18,height = 5)

####################################3.写回总的单细胞对象-精细注释----------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_TTTS_merge_named.qs")
TTTS_ExN<-qread("/home/iyun1806/MLdata/analysis/results/4_ExN/4_TTTS_ExN_named.qs")
TTTS_InN<-qread("/home/iyun1806/MLdata/analysis/results/4_InN/4_TTTS_InN_named.qs")

Idents(TTTS_filter_merge, cells = colnames(TTTS_ExN)) <- Idents(TTTS_ExN)
Idents(TTTS_filter_merge, cells = colnames(TTTS_InN)) <- Idents(TTTS_InN)


levels(TTTS_filter_merge)<-c("ExN-IPC","ExN-dev","ExN-L2/3","ExN-L3/4","ExN-L4","ExN-L4/5","ExN-L5PT","ExN-L5/6","ExN-L6CT","ExN-L6","InN-CGE-dev","InN-CGE","InN-MGE-dev","InN-MGE","NPC","Microglia","OPC","NOL","Astrocyte","Vascular","Poor quality")
col_all<-c("#d4dbff","#b8c4f5","#9caef0","#8098e8","#6583dd",
             "#4f6dcc","#4158b5","#3a449e","#343087","#2e1b70","#c8d96f", "#9cbf3e", "#6b8b23", "#3d6b2e","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
names(col_all)<-levels(TTTS_filter_merge)
DimPlot(TTTS_filter_merge,reduction = "umap",raster = F,label = F,cols = col_all) + 
  #  scale_color_d3(palette = "category20")+
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
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/4_TTTS_all_Dimplot.pdf",width = 8.5,height = 4.5)

DimPlot(TTTS_filter_merge,reduction = "umap",raster = F,label = F,cols = col_all) + 
  #  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')+NoLegend()
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/4_TTTS_all_Dimplot2.pdf",width = 4.5,height = 4.5)

TTTS_filter_merge$last_cluster<-Idents(TTTS_filter_merge)

###########3.1.画harmony之前的umap图----------
TTTS_filter_merge <- RunUMAP(TTTS_filter_merge, reduction = "pca", dims = 1:50,
                    reduction.name = "umap_before_harmony",
                    reduction.key = "UMAPb_")
DimPlot(TTTS_filter_merge,reduction = "umap_before_harmony",raster = F,label = F,cols = col_all) + 
  #  scale_color_d3(palette = "category20")+
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
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/4_all_cluster_Dimplot_beforeHarmony.pdf",width = 8.5,height = 4.5)


####################################4.修改之前的单细胞对象-粗略注释----------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/4_TTTS_merge_all.qs")
TTTS_filter_merge$broad_cluster[TTTS_filter_merge$last_cluster == "Poor quality"] <- "Poor quality"

qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/5_TTTS_merge_all.qs")
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/5_TTTS_merge_all.qs")

Idents(TTTS_filter_merge)<-TTTS_filter_merge$broad_cluster

levels(TTTS_filter_merge)<-c("ExN","InN","NPC","Microglia","OPC","NOL","Astrocyte","Vascular","Poor quality")
col_broad<-c("#324376","#6b8b23","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
names(col_broad)<-levels(TTTS_filter_merge)
DimPlot(TTTS_filter_merge,reduction = "umap",raster = F,label = F,cols = col_broad) + 
  #  scale_color_d3(palette = "category20")+
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
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/5_broad_cluster_Dimplot.pdf",width = 6.3,height = 4.5)

###########4.1.画harmony之前的umap图----------
TTTS_filter_merge <- RunUMAP(TTTS_filter_merge, reduction = "pca", dims = 1:50,
                             reduction.name = "umap_before_harmony",
                             reduction.key = "UMAPb_")
DimPlot(TTTS_filter_merge,reduction = "umap_before_harmony",raster = F,label = F,cols = col_broad) + 
  #  scale_color_d3(palette = "category20")+
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
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/5_broad_cluster_Dimplot_beforeHarmony.pdf",width = 6.3,height = 4.5)




#########################5.精细分群的细胞计数-------------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/5_TTTS_merge_all.qs")
dat<-TTTS_filter_merge@meta.data[,c(1:12,37,38)]
dat$`Age/PCW`<-as.factor(dat$`Age/PCW`)

##############5.1 精细分群统计计数--------------
count_by_cluster <- dat %>% group_by(last_cluster) %>% summarise(n = n())


# 图1：细胞数柱状图
col_all<-c("#d4dbff","#b8c4f5","#9caef0","#8098e8","#6583dd",
           "#4f6dcc","#4158b5","#3a449e","#343087","#2e1b70","#c8d96f", "#9cbf3e", "#6b8b23", "#3d6b2e","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")

ggplot(count_by_cluster, aes(last_cluster, n, fill = last_cluster)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = col_all) +
  labs(x = "", y = "Cell count") + theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "none"
  )
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/6_cellCount.pdf",width = 15,height = 6.5)

################5.2百分比图(按照年龄统计大的分群)------------------
col_broad<-c("#324376","#6b8b23","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
count_by_Age  <- dat %>% group_by(`Age/PCW`, broad_cluster) %>% summarise(n = n())

# 占比（样本内归一化）
pct_by_age <- count_by_Age %>%
  group_by(`Age/PCW`) %>% mutate(percent = n / sum(n) * 100) %>% ungroup()
ggplot(pct_by_age, aes(x = `Age/PCW`, y = percent, fill = broad_cluster)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  # geom_text(aes(label = sprintf("%.1f%%", percent)),
  #           position = position_stack(vjust = 0.5),
  #           color = "white", size = 3) +
  scale_fill_manual(values = col_broad) +
  labs(x = "Age (PCW)", y = "Percentage (%)", fill = "Cluster") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "right"
  )+theme_bw(base_size = 14)+theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
  )
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/6_cellpercent.pdf",width = 8,height = 6.5)


##################5.3柱状图(按照年龄统计大的分群)--------------
count_by_Age <- dat %>% group_by(`Age/PCW`) %>% summarise(n = n())
ggplot(count_by_Age, aes(`Age/PCW`, n, fill = `Age/PCW`)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = col_broad) +
  labs(x = "Age (PCW)", y = "Cell count") + theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "none"
  )
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/6_Age_cellCount.pdf",width = 6,height = 5.5)


##################5.3柱状图(按照数据集来源统计)--------------
col_dataset<-c("#7A142C","#5D90BA","#431A3D")
count_by_dataset <- dat %>% group_by(dataset) %>% summarise(n = n())
ggplot(count_by_dataset, aes(dataset, n, fill = dataset)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = col_dataset) +
  labs(x = "Dataset", y = "Cell count") + theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "none"
  )
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/6_dataset_cellCount.pdf",width = 5.5,height = 6.5)

##################5.4百分比图(按照数据集来源)--------------
col_broad<-c("#324376","#6b8b23","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
count_by_dataset  <- dat %>% group_by(dataset, broad_cluster) %>% summarise(n = n())

# 占比（样本内归一化）
pct_by_dataset <- count_by_dataset %>%
  group_by(dataset) %>% mutate(percent = n / sum(n) * 100) %>% ungroup()
ggplot(pct_by_dataset, aes(x = dataset, y = percent, fill = broad_cluster)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  # geom_text(aes(label = sprintf("%.1f%%", percent)),
  #           position = position_stack(vjust = 0.5),
  #           color = "white", size = 3) +
  scale_fill_manual(values = col_broad) +
  labs(x = "Dataset", y = "Percentage (%)", fill = "Cluster") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "right"
  )+theme_bw(base_size = 14)+theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
  )
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/6_datasetpercent.pdf",width = 7.5,height = 6.5)

##################6,提取自己的4个样本重新分群--------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/5_TTTS_merge_all.qs")


####################6.1 D1---------------

##############################6.1.1 Dimplot图-------------------
D1<-subset(TTTS_filter_merge,orig.ident=="D1")

D1 <- D1 %>% NormalizeData(verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)

Idents(D1)<-D1$last_cluster

levels(TTTS_filter_merge)<-c("ExN-IPC","ExN-dev","ExN-L2/3","ExN-L3/4","ExN-L4","ExN-L4/5","ExN-L5PT","ExN-L5/6","ExN-L6CT","ExN-L6","InN-CGE-dev","InN-CGE","InN-MGE-dev","InN-MGE","NPC","Microglia","OPC","NOL","Astrocyte","Vascular","Poor quality")
col_all<-c("#d4dbff","#b8c4f5","#9caef0","#8098e8","#6583dd",
           "#4f6dcc","#4158b5","#3a449e","#343087","#2e1b70","#c8d96f", "#9cbf3e", "#6b8b23", "#3d6b2e","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
names(col_all)<-levels(TTTS_filter_merge)

library(mascarade)
df_D1 <- data.frame(Embeddings(D1, "umap"), celltype = D1$last_cluster)
write.csv(df,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_df_D1.csv",row.names = T)
df<-read.csv("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_df_D1.csv")

# 2) 生成每个群的遮罩（轮廓线坐标表）
mask <- generateMask(dims    = Embeddings(D1, "umap"),
                     clusters = D1$last_cluster,
                     expand   = 0.01)   # 圈往外扩一点，别贴着点
# DimPlot(D1,reduction = "umap",raster = F,label = T,cols = col_all,repel = T) + 
#   #  scale_color_d3(palette = "category20")+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_text(size = 16),
#         legend.text = element_text(size = 16),
#         legend.title = element_blank(),
#         plot.title = element_text(size = 16))+
#   xlab('UMAP_1')+
#   ylab('UMAP_2')+geom_path(data = mask,
#                            aes(x = umap_1, y = umap_2, group = group),
#                            linetype = "dashed", linewidth = 0.6) +NoLegend()


DimPlot(D1,reduction = "umap",raster = F,label = F,cols = col_all,repel = T) + 
  #  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')+geom_path(data = mask,
                           aes(x = umap_1, y = umap_2, group = group),
                           linetype = "dashed", linewidth = 0.6) +NoLegend()

ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_D1_Dimplot.pdf",width = 4.5,height = 4.5)


##########################6.1.2 柱状图+百分比图---------------------
dat_D1<-D1@meta.data[,c(1:12,37,38)]
dat_D1$`Age/PCW`<-as.factor(dat$`Age/PCW`)
# 判断有细分的大类
refined_D1 <- dat %>%
  group_by(broad_cluster) %>%
  summarise(n_sub = n_distinct(last_cluster), .groups = "drop") %>%
  filter(n_sub > 1) %>%
  pull(broad_cluster) %>%
  as.character()

col_broad<-c("#324376","#6b8b23","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
names(col_broad)<-c("ExN","InN","NPC","Microglia","OPC","NOL","Astrocyte","Vascular","Poor quality")

col_all<-c("#d4dbff","#b8c4f5","#9caef0","#8098e8","#6583dd",
           "#4f6dcc","#4158b5","#3a449e","#343087","#2e1b70","#c8d96f", "#9cbf3e", "#6b8b23", "#3d6b2e","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
names(col_all)<-c("ExN-IPC","ExN-dev","ExN-L2/3","ExN-L3/4","ExN-L4","ExN-L4/5","ExN-L5PT","ExN-L5/6","ExN-L6CT","ExN-L6","InN-CGE-dev","InN-CGE","InN-MGE-dev","InN-MGE","NPC","Microglia","OPC","NOL","Astrocyte","Vascular","Poor quality")

plot_df_D1 <- dat_D1 %>%
  mutate(
    broad_cluster = as.character(broad_cluster),
    last_cluster  = as.character(last_cluster)
  ) %>%
  group_by(broad_cluster, last_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    fill_id = ifelse(broad_cluster %in% refined_D1, last_cluster, broad_cluster)
  ) %>%
  mutate(
    broad_cluster = factor(broad_cluster, levels = names(col_broad)),
    fill_id = factor(fill_id, levels = names(col_all))
  )

####百分比标签数
pct_labels_D1 <- plot_df_D1 %>%
  filter(broad_cluster %in% refined_D1) %>%
  group_by(broad_cluster, .drop = FALSE) %>%
  arrange(fill_id) %>%                 # ← 这一行是核心，不能省
  mutate(
    ypos = cumsum(n) - n / 2,          # 每段真实几何中心
    percent = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  filter(percent >= 1.5)               # 太小的不标，避免溢出
write.csv(pct_labels_D1,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_D1_pecent.csv",row.names = T)
pct_labels_D1<-read.csv("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_D1.csv")

# ---- 5. 准备柱顶总数标签 
count_summary_D1 <- plot_df_D1 %>%
  group_by(broad_cluster) %>%
  summarise(n = sum(n), .groups = "drop")

# ---- 6. 绘图
ggplot(plot_df_D1, aes(broad_cluster, n, fill = fill_id)) +
  geom_col(width = 0.62, color = "white") +
  geom_text(data = count_summary_D1,
            aes(broad_cluster, n, label = n),
            inherit.aes = FALSE, vjust = -0.35, size = 3.5) +
  # geom_text(data = pct_labels,
  #           aes(broad_cluster, ypos, label = sprintf("%.1f%%", percent)),
  #           inherit.aes = FALSE, color = "white", size = 2.5) +   # ← 改这里：2.5 比原来小一号
  scale_fill_manual(values = col_all, name = "Subtype (last_cluster)") +
  labs(x = "", y = "Cell count") +
  guides(fill = "none") +                                          # ← 新增：隐藏图例
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_D1_countPlot.pdf",width = 5.5,height = 5)

###################6.2 R1-----------------------
R1<-subset(TTTS_filter_merge,orig.ident=="R1")

R1 <- R1 %>% NormalizeData(verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)

Idents(R1)<-R1$last_cluster
df_R1 <- data.frame(Embeddings(R1, "umap"), celltype = R1$last_cluster)
# 2) 生成每个群的遮罩（轮廓线坐标表）
mask <- generateMask(dims    = Embeddings(R1, "umap"),
                     clusters = R1$last_cluster,
                     expand   = 0.01)   # 圈往外扩一点，别贴着点

DimPlot(R1,reduction = "umap",raster = F,label = F,cols = col_all,repel = T) + 
  #  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')+geom_path(data = mask,
                           aes(x = umap_1, y = umap_2, group = group),
                           linetype = "dashed", linewidth = 0.6) +NoLegend()

ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_R1_Dimplot.pdf",width = 4.5,height = 4.5)


#########################6.2.2 柱状图+百分比图-----------------
dat_R1<-R1@meta.data[,c(1:12,37,38)]
# 判断有细分的大类
refined_R1 <- dat %>%
  group_by(broad_cluster) %>%
  summarise(n_sub = n_distinct(last_cluster), .groups = "drop") %>%
  filter(n_sub > 1) %>%
  pull(broad_cluster) %>%
  as.character()


plot_df_R1 <- dat_R1 %>%
  mutate(
    broad_cluster = as.character(broad_cluster),
    last_cluster  = as.character(last_cluster)
  ) %>%
  group_by(broad_cluster, last_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    fill_id = ifelse(broad_cluster %in% refined_R1, last_cluster, broad_cluster)
  ) %>%
  mutate(
    broad_cluster = factor(broad_cluster, levels = names(col_broad)),
    fill_id = factor(fill_id, levels = names(col_all))
  )

####百分比标签数
pct_labels_R1 <- plot_df_R1 %>%
  filter(broad_cluster %in% refined_R1) %>%
  group_by(broad_cluster, .drop = FALSE) %>%
  arrange(fill_id) %>%                 # ← 这一行是核心，不能省
  mutate(
    ypos = cumsum(n) - n / 2,          # 每段真实几何中心
    percent = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  filter(percent >= 1.5)               # 太小的不标，避免溢出
write.csv(pct_labels_R1,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_R1_pecent.csv",row.names = T)
pct_labels_R1<-read.csv("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_R1.csv")

# ---- 5. 准备柱顶总数标签 
count_summary_R1 <- plot_df_R1 %>%
  group_by(broad_cluster) %>%
  summarise(n = sum(n), .groups = "drop")

# ---- 6. 绘图
ggplot(plot_df_R1, aes(broad_cluster, n, fill = fill_id)) +
  geom_col(width = 0.62, color = "white") +
  geom_text(data = count_summary_R1,
            aes(broad_cluster, n, label = n),
            inherit.aes = FALSE, vjust = -0.35, size = 3.5) +
  # geom_text(data = pct_labels,
  #           aes(broad_cluster, ypos, label = sprintf("%.1f%%", percent)),
  #           inherit.aes = FALSE, color = "white", size = 2.5) +   # ← 改这里：2.5 比原来小一号
  scale_fill_manual(values = col_all, name = "Subtype (last_cluster)") +
  labs(x = "", y = "Cell count") +
  guides(fill = "none") +                                          # ← 新增：隐藏图例
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_R1_countPlot.pdf",width = 5.5,height = 5)
###########################6.3 D2----------------------
D2<-subset(TTTS_filter_merge,orig.ident=="D2")

D2 <- D2 %>% NormalizeData(verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)

Idents(D2)<-D2$last_cluster
# 2) 生成每个群的遮罩（轮廓线坐标表）
mask <- generateMask(dims    = Embeddings(D2, "umap"),
                     clusters = D2$last_cluster,
                     expand   = 0.01)   # 圈往外扩一点，别贴着点

DimPlot(D2,reduction = "umap",raster = F,label = F,cols = col_all,repel = T) + 
  #  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')+geom_path(data = mask,
                           aes(x = umap_1, y = umap_2, group = group),
                           linetype = "dashed", linewidth = 0.6) +NoLegend()

ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_D2_Dimplot.pdf",width = 4.5,height = 4.5)

##########################6.3.2 柱状图+百分比图---------------------
dat_D2<-D2@meta.data[,c(1:12,37,38)]
# 判断有细分的大类
refined_D2 <- dat %>%
  group_by(broad_cluster) %>%
  summarise(n_sub = n_distinct(last_cluster), .groups = "drop") %>%
  filter(n_sub > 1) %>%
  pull(broad_cluster) %>%
  as.character()


plot_df_D2 <- dat_D2 %>%
  mutate(
    broad_cluster = as.character(broad_cluster),
    last_cluster  = as.character(last_cluster)
  ) %>%
  group_by(broad_cluster, last_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    fill_id = ifelse(broad_cluster %in% refined_D2, last_cluster, broad_cluster)
  ) %>%
  mutate(
    broad_cluster = factor(broad_cluster, levels = names(col_broad)),
    fill_id = factor(fill_id, levels = names(col_all))
  )

####百分比标签数
pct_labels_D2 <- plot_df_D2 %>%
  filter(broad_cluster %in% refined_D2) %>%
  group_by(broad_cluster, .drop = FALSE) %>%
  arrange(fill_id) %>%                 # ← 这一行是核心，不能省
  mutate(
    ypos = cumsum(n) - n / 2,          # 每段真实几何中心
    percent = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  filter(percent >= 1.5)               # 太小的不标，避免溢出
write.csv(pct_labels_D2,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_D2_pecent.csv",row.names = T)
pct_labels_D2<-read.csv("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_D2.csv")

# ---- 5. 准备柱顶总数标签 
count_summary_D2 <- plot_df_D2 %>%
  group_by(broad_cluster) %>%
  summarise(n = sum(n), .groups = "drop")

# ---- 6. 绘图
ggplot(plot_df_D2, aes(broad_cluster, n, fill = fill_id)) +
  geom_col(width = 0.62, color = "white") +
  geom_text(data = count_summary_D2,
            aes(broad_cluster, n, label = n),
            inherit.aes = FALSE, vjust = -0.35, size = 3.5) +
  # geom_text(data = pct_labels,
  #           aes(broad_cluster, ypos, label = sprintf("%.1f%%", percent)),
  #           inherit.aes = FALSE, color = "white", size = 2.5) +   # ← 改这里：2.5 比原来小一号
  scale_fill_manual(values = col_all, name = "Subtype (last_cluster)") +
  labs(x = "", y = "Cell count") +
  guides(fill = "none") +                                          # ← 新增：隐藏图例
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_D2_countPlot.pdf",width = 5.5,height = 5)


###########################6.3 R2----------------------
R2<-subset(TTTS_filter_merge,orig.ident=="R2")

R2 <- R2 %>% NormalizeData(verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>%
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)

Idents(R2)<-R2$last_cluster
# 2) 生成每个群的遮罩（轮廓线坐标表）
mask <- generateMask(dims    = Embeddings(R2, "umap"),
                     clusters = R2$last_cluster,
                     expand   = 0.01)   # 圈往外扩一点，别贴着点

DimPlot(R2,reduction = "umap",raster = F,label = F,cols = col_all,repel = T) + 
  #  scale_color_d3(palette = "category20")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))+
  xlab('UMAP_1')+
  ylab('UMAP_2')+geom_path(data = mask,
                           aes(x = umap_1, y = umap_2, group = group),
                           linetype = "dashed", linewidth = 0.6) +NoLegend()

ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_R2_Dimplot.pdf",width = 4.5,height = 4.5)

##########################6.3.2 柱状图+百分比图---------------------
dat_R2<-R2@meta.data[,c(1:12,37,38)]
# 判断有细分的大类
refined_R2 <- dat %>%
  group_by(broad_cluster) %>%
  summarise(n_sub = n_distinct(last_cluster), .groups = "drop") %>%
  filter(n_sub > 1) %>%
  pull(broad_cluster) %>%
  as.character()


plot_df_R2 <- dat_R2 %>%
  mutate(
    broad_cluster = as.character(broad_cluster),
    last_cluster  = as.character(last_cluster)
  ) %>%
  group_by(broad_cluster, last_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    fill_id = ifelse(broad_cluster %in% refined_R2, last_cluster, broad_cluster)
  ) %>%
  mutate(
    broad_cluster = factor(broad_cluster, levels = names(col_broad)),
    fill_id = factor(fill_id, levels = names(col_all))
  )

####百分比标签数
pct_labels_R2 <- plot_df_R2 %>%
  filter(broad_cluster %in% refined_R2) %>%
  group_by(broad_cluster, .drop = FALSE) %>%
  arrange(fill_id) %>%                 # ← 这一行是核心，不能省
  mutate(
    ypos = cumsum(n) - n / 2,          # 每段真实几何中心
    percent = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  filter(percent >= 1.5)               # 太小的不标，避免溢出
write.csv(pct_labels_R2,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_R2_pecent.csv",row.names = T)
pct_labels_R2<-read.csv("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_R2.csv")

# ---- 5. 准备柱顶总数标签 
count_summary_R2 <- plot_df_R2 %>%
  group_by(broad_cluster) %>%
  summarise(n = sum(n), .groups = "drop")

# ---- 6. 绘图
ggplot(plot_df_R2, aes(broad_cluster, n, fill = fill_id)) +
  geom_col(width = 0.62, color = "white") +
  geom_text(data = count_summary_R2,
            aes(broad_cluster, n, label = n),
            inherit.aes = FALSE, vjust = -0.35, size = 3.5) +
  # geom_text(data = pct_labels,
  #           aes(broad_cluster, ypos, label = sprintf("%.1f%%", percent)),
  #           inherit.aes = FALSE, color = "white", size = 2.5) +   # ← 改这里：2.5 比原来小一号
  scale_fill_manual(values = col_all, name = "Subtype (last_cluster)") +
  labs(x = "", y = "Cell count") +
  guides(fill = "none") +                                          # ← 新增：隐藏图例
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/7_R2_countPlot.pdf",width = 5.5,height = 5)

#########################7.比较Dornor和Receptor--------------------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/5_TTTS_merge_all.qs")
dat<-TTTS_filter_merge@meta.data[,c(1:12,37,38)]
#####################7.1 D1和R1  ExN----------------
exn <- dat %>%
  filter(broad_cluster == "ExN", orig.ident %in% c("D1", "R1"))

pct <- exn %>%
  group_by(orig.ident, last_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

plot_df <- pct %>%
  select(orig.ident, last_cluster, percent) %>%
  pivot_wider(names_from = orig.ident, values_from = percent) %>%
  arrange(D1) %>%                                           # 按 D1 占比排序
  mutate(last_cluster = factor(last_cluster, levels = last_cluster))
write.csv(plot_df,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/8_D1R1_pecent.csv",row.names = T)

ggplot(plot_df, aes(x = last_cluster)) +
  geom_segment(aes(y = D1, yend = R1), color = "grey60", linewidth = 1) +  # 连接杆
  geom_point(aes(y = D1), color = "#AD494A", size = 4) +                       # D1
  geom_point(aes(y = R1), color = "#E7CB94", size = 4) +                       # R1
  geom_text(aes(y = D1, label = sprintf("%.1f%%", D1)), vjust = 1.3, size = 4) +
  geom_text(aes(y = R1, label = sprintf("%.1f%%", R1)), vjust = -1.3, size = 4) +
#  coord_flip() +                                                          # 让亚群名横向排布（可选）
  labs(x = "", y = "Percentage within ExN (%)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 15,
          color = "black"
        ),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/8_D1R1_Plot.pdf",width = 4,height = 4.5)

#####################7.2 D1和R1  InN----------------
InN <- dat %>%
  filter(broad_cluster == "InN", orig.ident %in% c("D1", "R1"))

pct <- InN %>%
  group_by(orig.ident, last_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

plot_df <- pct %>%
  select(orig.ident, last_cluster, percent) %>%
  pivot_wider(names_from = orig.ident, values_from = percent) %>%
  arrange(D1) %>%                                           # 按 D1 占比排序
  mutate(last_cluster = factor(last_cluster, levels = last_cluster))
write.csv(plot_df,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/8_D1R1_InN_pecent.csv",row.names = T)
ggplot(plot_df, aes(x = last_cluster)) +
  geom_segment(aes(y = D1, yend = R1), color = "grey60", linewidth = 1) +  # 连接杆
  geom_point(aes(y = D1), color = "#AD494A", size = 4) +                       # D1
  geom_point(aes(y = R1), color = "#E7CB94", size = 4) +                       # R1
  geom_text(aes(y = D1, label = sprintf("%.1f%%", D1)), vjust = 2, size = 4) +
  geom_text(aes(y = R1, label = sprintf("%.1f%%", R1)), vjust = -1.5, size = 4) +
  #  coord_flip() +                                                          # 让亚群名横向排布（可选）
  labs(x = "", y = "Percentage within InN (%)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 15,
          color = "black"
        ),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/8_D1R1_InN_Plot.pdf",width = 3,height = 4.5)

#####################7.2 D2和R2  ExN----------------
ExN <- dat %>%
  filter(broad_cluster == "ExN", orig.ident %in% c("D2", "R2"))

pct <- ExN %>%
  group_by(orig.ident, last_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()


plot_df <- pct %>%
  select(orig.ident, last_cluster, percent) %>%
  pivot_wider(names_from = orig.ident, values_from = percent) %>%
  arrange(D2) %>%                                           # 按 D2 占比排序
  mutate(last_cluster = factor(last_cluster, levels = last_cluster))
write.csv(plot_df,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/8_D2R2_ExN_pecent.csv",row.names = T)
ggplot(plot_df, aes(x = last_cluster)) +
  geom_segment(aes(y = D2, yend = R2), color = "grey60", linewidth = 1) +  # 连接杆
  geom_point(aes(y = D2), color = "#AD494A", size = 4) +                       # D2
  geom_point(aes(y = R2), color = "#E7CB94", size = 4) +                       # R2
  geom_text(aes(y = D2, label = sprintf("%.1f%%", D2)), vjust = 2, size = 4) +
  geom_text(aes(y = R2, label = sprintf("%.1f%%", R2)), vjust = -1.5, size = 4) +
  #  coord_flip() +                                                          # 让亚群名横向排布（可选）
  labs(x = "", y = "Percentage within ExN (%)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 15,
          color = "black"
        ),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/8_D2R2_ExN_Plot.pdf",width = 4,height = 4.5)

#####################7.4 D2和R2  InN----------------
InN <- dat %>%
  filter(broad_cluster == "InN", orig.ident %in% c("D2", "R2"))

pct <- InN %>%
  group_by(orig.ident, last_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

plot_df <- pct %>%
  select(orig.ident, last_cluster, percent) %>%
  pivot_wider(names_from = orig.ident, values_from = percent) %>%
  arrange(D2) %>%                                           # 按 D2 占比排序
  mutate(last_cluster = factor(last_cluster, levels = last_cluster))
write.csv(plot_df,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/8_D2R2_InN_pecent.csv",row.names = T)

ggplot(plot_df, aes(x = last_cluster)) +
  geom_segment(aes(y = D2, yend = R2), color = "grey60", linewidth = 1) +  # 连接杆
  geom_point(aes(y = D2), color = "#AD494A", size = 4) +                       # D2
  geom_point(aes(y = R2), color = "#E7CB94", size = 4) +                       # R2
  geom_text(aes(y = D2, label = sprintf("%.1f%%", D2)), vjust = 2, size = 4) +
  geom_text(aes(y = R2, label = sprintf("%.1f%%", R2)), vjust = -1.5, size = 4) +
  #  coord_flip() +                                                          # 让亚群名横向排布（可选）
  labs(x = "", y = "Percentage within InN (%)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 15,
          color = "black"
        ),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/8_D2R2_InN_Plot.pdf",width = 3,height = 4.5)

#######################8.按照样本计数--------------------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/5_TTTS_merge_all.qs")
dat<-TTTS_filter_merge@meta.data[,c(1:12,37,38)]
dat$`Age/PCW`<-as.factor(dat$`Age/PCW`)
col_broad <- c("#324376","#6b8b23","#bbbbbb","#c63312","#e8c92a","#c7a86c","#8b2043","#99934a","black")
names(col_broad) <- c("ExN","InN","NPC","Microglia","OPC","NOL","Astrocyte","Vascular","Poor quality")

plot_df <- dat %>%
  group_by(orig.ident, broad_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percent = n / sum(n) * 100) %>%      # ← 样本内归一化
  ungroup() %>%
  mutate(
    broad_cluster = factor(broad_cluster, levels = names(col_broad)),
    orig.ident        = factor(orig.ident, levels = unique(TTTS_filter_merge$orig.ident))   # 横轴样本顺序
  )

# 柱顶标总细胞数
count_orig.ident <- plot_df %>%
  group_by(orig.ident) %>% summarise(n = sum(n), .groups = "drop")

# 自检：每样本百分比总和 = 100
plot_df %>% group_by(orig.ident) %>% summarise(total = sum(percent))
write.csv(plot_df,"/home/iyun1806/MLdata/analysis/results/3_cluster_new/9_sample_pecent.csv",row.names = T)
ggplot(plot_df, aes(x = orig.ident, y = n, fill = broad_cluster)) +
  geom_col(width = 0.55, color = "white") +
  # 柱内部：百分比（堆叠居中）
  geom_text(aes(label = ifelse(percent >= 3, sprintf("%.1f%%", percent), "")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3) +
  # 柱顶：总细胞数
  geom_text(data = count_orig.ident, aes(x = orig.ident, y = n, label = n),
            inherit.aes = FALSE, vjust = -0.35, size = 4) +
  scale_fill_manual(values = col_broad, name = "Major class") +
  labs(x = "Sample", y = "Cell count") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 14))
ggsave("/home/iyun1806/MLdata/analysis/results/3_cluster_new/9_sample_pecent.pdf",width =15,height = 6.5 )











