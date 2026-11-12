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

TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/3_cluster_new/3_TTTS_merge_named.qs")
TTTS_InN<-subset(TTTS_filter_merge,idents = "InN")
TTTS_InN@meta.data<-TTTS_InN@meta.data[,c(1:13,36,37)]
TTTS_InN <- TTTS_InN %>% NormalizeData(verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>%
  RunHarmony("orig.ident",plot_convergence = F) %>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50)%>% 
  FindClusters(resolution=seq(0,1,0.1))
TTTS_InN$num_cluster2<-Idents(TTTS_InN)

qsave(TTTS_InN,"/home/iyun1806/MLdata/analysis/results/4_InN/1_TTTSInN_cluster.qs")
TTTS_InN<-qread("/home/iyun1806/MLdata/analysis/results/4_InN/1_TTTSInN_cluster.qs")

p1<-clustree(TTTS_InN,prefix='RNA_snn_res.')
ggsave(filename = "/home/iyun1806/MLdata/analysis/results/4_InN/1_InN_clustree.pdf",plot = p1,width = 15,height = 15)

Idents(TTTS_InN)<-TTTS_InN$RNA_snn_res.0.3
DimPlot(TTTS_InN,reduction = "umap",raster = F,label = T) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/1_InN_diimplot.png",width = 4.5,height = 4.5)

##求差异基因
DEGs_InN <-FindAllMarkers(TTTS_InN,logfc.threshold=0.25,min.pct=0.25,only.pos = T)
DEGs_InN <- DEGs_InN %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(DEGs_InN,"/home/iyun1806/MLdata/analysis/results/4_InN/2_DEGs_InN.csv",row.names = T)
DEGs_InN<-read.csv("/home/iyun1806/MLdata/analysis/results/4_InN/2_DEGs_TTTS.csv")

subset_del<-subset(TTTS_InN,idents = c("9","8","10"))
table(subset_del$num_cluster)


##########################注释----------
TTTS_InN <- RenameIdents(TTTS_InN,
                           "0"  = "InN-CGE-dev",
                           "3"  = "InN-MGE-dev",
                           "2"  = "InN-CGE", "5" = "InN-CGE", "7" = "InN-CGE", "12" = "InN-CGE",
                           "1"  = "InN-MGE", "4" = "InN-MGE", "6" = "InN-MGE", "11" = "InN-MGE",
                           "8"  = "Poor quality", "9" = "Poor quality", "10" = "Poor quality"
                         )
TTTS_InN$InN_cluster<-Idents(TTTS_InN)

qsave(TTTS_InN,"/home/iyun1806/MLdata/analysis/results/4_InN/4_TTTS_InN_named.qs")
TTTS_InN<-qread("/home/iyun1806/MLdata/analysis/results/4_InN/4_TTTS_InN_named.qs")

#######################重新绘图-----------------
TTTS_InN2<-subset(TTTS_InN,idents = "Poor quality",invert=T)
# TTTS_InN2 <- TTTS_InN2 %>% NormalizeData(verbose = T) %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = 5000, verbose = T) %>% 
#   ScaleData(verbose = T) %>% 
#   RunPCA(npcs = 50, verbose = T) %>%
#   RunHarmony("orig.ident",plot_convergence = F) %>% 
#   RunUMAP(reduction = "harmony", dims = 1:50, verbose = T)

Idents(TTTS_InN2)<-TTTS_InN2$InN_cluster

levels(TTTS_InN2)<-c("InN-CGE-dev","InN-CGE","InN-MGE-dev","InN-MGE")
col_InN <- c("#c8d96f", "#9cbf3e", "#6b8b23", "#3d6b2e")
names(col_InN)<-levels(TTTS_InN2)

DimPlot(TTTS_InN2,reduction = "umap",raster = F,label = T,cols = col_InN) + 
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

ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/4_InN_cluster_Dimplot.pdf",width = 4.5,height = 4.5)

###########3.画harmony之前的umap图----------
TTTS_InN2 <- RunUMAP(TTTS_InN2, reduction = "pca", dims = 1:50,
                    reduction.name = "umap_before_harmony",
                    reduction.key = "UMAPb_")
DimPlot(TTTS_InN2,reduction = "umap_before_harmony",raster = F,label = F,cols = col_InN) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/4_InN_cluster_beforeharmony_Dimplot.pdf",width = 6.5,height = 4.5)

qsave(TTTS_InN2,"/home/iyun1806/MLdata/analysis/results/4_InN/4_TTTS_InN2_named.qs")
TTTS_InN2<-qread("/home/iyun1806/MLdata/analysis/results/4_InN/4_TTTS_InN2_named.qs")

###########4.画marker基因的Dotplot图----------
# InN_genes <- c(
#   "SCGN", "SP8", "NR2F2", "PROX1", "DLX2", "DCX", "STMN2",  # InN-CGE_dev
#   "LHX6", "MAF", "DLX2", "DCX", "STMN2",                     # InN-MGE_dev
#   "ADARB2", "NR2F2", "PROX1", "KIRREL3", "KCNIP1", "SYT10",  # InN-CGE
#   "SOX6", "LHX6", "ERBB4", "GRIA4", "NPY", "NOS1"            # InN-MGE
# )
InN_genes <-c(
  # InN-CGE_dev (c0)：CGE 来源特异 + 未成熟
  "ADARB2", "SCGN", "SP8", "NR2F2", "PROX1", "CALB2", "NR2F1",
  "DLX1", "DLX2", "DLX5", "DCX", "STMN2", "SOX11", "SOX4",
  "SP9", "ASCL1", "INSM1", "ID4", "CXCR4",
  # InN-MGE_dev (c3)：MGE 来源特异 + 未成熟
  "LHX6", "SOX6", "MAF", "ERBB4", "GRIA4", "SST",
  "DLX2", "DLX5", "DCX", "STMN2", "SOX11", "SOX4",
  "SP9", "ASCL1", "CXCR4",
  # InN-CGE (c2,5,7,12)：CGE 来源特异 + 成熟特征
  "ADARB2", "NR2F2", "PROX1", "SP8", "SCGN",
  "KIRREL3", "KCNIP1", "CCK", "CRH", "RGS6", "SYT10", "EYA2",
  # InN-MGE (c1,4,6,11)：MGE 来源特异 + 成熟特征
  "LHX6", "SOX6", "MAF", "ERBB4", "GRIA4",
  "TRHDE", "GRIK3", "CHRM2", "DPP6", "RELN", "PDE1A",
  "NPY", "NOS1", "RIT2"
)
InN_genes<-unique(InN_genes)
DotPlot(
  TTTS_InN2,
  features = InN_genes,
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
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/4_InN_Dotplot.pdf",width = 21,height = 5.5)

##################density-----------------
library(ggpointdensity)
dat1 <- FetchData(object=TTTS_InN2, vars=c("umap_1","umap_2","orig.ident"))
colnames(dat1)<-c("UMAP_1","UMAP_2","sample")
head(dat1)

ggplot(data = dat1, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_pointdensity(size=1) +
  scale_color_viridis_c(option="D", alpha = 1,begin = 0.35,end = 0.75) + 
  theme_bw(base_size = 15) +
  facet_wrap(~sample,ncol=6)+
  theme(
    strip.text = element_text(size = 25),       
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    plot.title = element_text(size = 16)
  )+
  stat_density_2d(
    linewidth = 0.8,
    colour = "#2E8B57"
  ) +
  labs(color="counts")+NoLegend()
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/5_InN2_density.pdf",width = 18.5,height = 8)

dat2 <- FetchData(object=TTTS_InN2, vars=c("umap_1","umap_2","Age/PCW"))
colnames(dat2)<-c("UMAP_1","UMAP_2","Age")
head(dat2)

ggplot(data = dat2, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_pointdensity(size=1) +
  scale_color_viridis_c(option="D", alpha = 1,begin = 0.35,end = 0.75) + 
  theme_bw(base_size = 15) +
  facet_wrap(~Age,ncol=6)+
  theme(
    strip.text = element_text(size = 25),       
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    plot.title = element_text(size = 16)
  )+
  stat_density_2d(
    linewidth = 0.8,
    colour = "#2E8B57"
  ) +
  labs(color="counts")+NoLegend()
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/5_InN2_density_age.pdf",width = 14.5,height = 4)


###########5.Ro/e----------
######Age=20
library(ComplexHeatmap)
library(circlize)
meta20 <- TTTS_InN2@meta.data %>% 
  filter(`Age/PCW` == "20")
tab_samp <- table(meta20$InN_cluster, meta20$orig.ident)
tab_samp <- tab_samp[, colSums(tab_samp) > 0]      # 去掉无细胞的样本（保险起见）

fit_samp <- chisq.test(tab_samp)
roe_samp <- fit_samp$observed / fit_samp$expected  # 行 = 亚群，列 = 样本

col_fun <- colorRamp2(c(0, 1, 3), c("#2166AC", "white", "#B2182B"))
Heatmap(roe_samp,
        show_heatmap_legend = T, 
        cluster_rows = F, 
        cluster_columns = T,
        row_names_side = 'right', 
        show_column_names = T,
        show_row_names = T,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        # heatmap_legend_param = list(
        #   title = "R o/e",
        #   at =  seq(0, 3, by = 0.5),
        #   labels =  seq(0, 3, by = 0.5),
        #   legend_gp = gpar(fill = col_fun(at))
        # ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", roe_samp[i, j]), x, y, gp = gpar(fontsize = 10))
        }
)
##只有20.（5，4）
######Age=16
library(ComplexHeatmap)
library(circlize)
meta16 <- TTTS_InN2@meta.data %>% 
  filter(`Age/PCW` == "16") 
tab_samp <- table(meta16$InN_cluster, meta16$orig.ident)
tab_samp <- tab_samp[, colSums(tab_samp) > 0]      # 去掉无细胞的样本（保险起见）
tab_samp <- tab_samp[rowSums(tab_samp) > 0,]      # 去掉无细胞的亚群（保险起见）

fit_samp <- chisq.test(tab_samp)
roe_samp <- fit_samp$observed / fit_samp$expected  # 行 = 亚群，列 = 样本

col_fun <- colorRamp2(c(0, 1, 3), c("#2166AC", "white", "#B2182B"))
Heatmap(roe_samp,
        show_heatmap_legend = T, 
        cluster_rows = F, 
        cluster_columns = T,
        row_names_side = 'right', 
        show_column_names = T,
        show_row_names = T,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        # heatmap_legend_param = list(
        #   title = "R o/e",
        #   at =  seq(0, 3, by = 0.5),
        #   labels =  seq(0, 3, by = 0.5),
        #   legend_gp = gpar(fill = col_fun(at))
        # ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", roe_samp[i, j]), x, y, gp = gpar(fontsize = 10))
        }
)

################7.画16，20的umap图----------------
TTTS_InN2<-qread("/home/iyun1806/MLdata/analysis/results/4_InN/4_TTTS_InN2_named.qs")
library(mascarade)


TTTS_age_16<-subset(TTTS_InN2,subset=`Age/PCW` ==16)
col_InN <- c("#c8d96f", "#9cbf3e", "#6b8b23", "#3d6b2e")

DimPlot(TTTS_age_16,reduction = "umap",raster = F,label = T,cols = col_InN,repel = T) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/7_InN2_Age16_Dimplot.pdf",width = 4.5,height = 4.5)

TTTS_age_20<-subset(TTTS_InN2,subset=`Age/PCW` ==20)
DimPlot(TTTS_age_20,reduction = "umap",raster = F,label = T,cols = col_InN,repel = T) + 
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
  ylab('UMAP_2') +NoLegend()
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/7_InN2_Age20_Dimplot.pdf",width = 4.5,height = 4.5)

################8.画16，20的细胞计数----------------
dat<-TTTS_age_16@meta.data[,c(1:12,29)]
count_by_cluster <- dat %>% group_by(InN_cluster) %>% summarise(n = n())

col_InN <- c("#c8d96f", "#9cbf3e", "#6b8b23", "#3d6b2e")
ggplot(count_by_cluster, aes(InN_cluster, n, fill = InN_cluster)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = col_InN) +
  labs(x = "", y = "Cell count") + theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "none"
  )
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/8_Age16_cellCount.pdf",width = 6,height = 6)

dat<-TTTS_age_20@meta.data[,c(1:12,29)]
count_by_cluster <- dat %>% group_by(InN_cluster) %>% summarise(n = n())

ggplot(count_by_cluster, aes(InN_cluster, n, fill = InN_cluster)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = col_InN) +
  labs(x = "", y = "Cell count") + theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "none"
  )
ggsave("/home/iyun1806/MLdata/analysis/results/4_InN/8_Age20_cellCount.pdf",width = 6,height = 6.5)
