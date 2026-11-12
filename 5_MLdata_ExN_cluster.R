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
TTTS_ExN<-subset(TTTS_filter_merge,idents = "ExN")
TTTS_ExN@meta.data<-TTTS_ExN@meta.data[,c(1:13,36,37)]
TTTS_ExN <- TTTS_ExN %>% NormalizeData(verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>%
  RunHarmony("orig.ident",plot_convergence = F) %>% 
  RunUMAP(reduction = "harmony", dims = 1:50, verbose = T) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50)%>% 
  FindClusters(resolution=seq(0,1,0.1))

qsave(TTTS_ExN,"/home/iyun1806/MLdata/analysis/results/4_ExN/1_TTTSExN_cluster.qs")
TTTS_ExN<-qread("/home/iyun1806/MLdata/analysis/results/4_ExN/1_TTTSExN_cluster.qs")

clustree(TTTS_ExN,prefix='RNA_snn_res.')
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/1_ExN_clustree.pdf",width = 10,height = 10)

Idents(TTTS_ExN)<-TTTS_ExN$RNA_snn_res.0.6
DimPlot(TTTS_ExN,reduction = "umap",raster = F,label = T) + 
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

ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/3_ExNDimplot.png",width = 5,height = 5)

##求差异基因
DEGs_ExN <-FindAllMarkers(TTTS_ExN,logfc.threshold=0.25,min.pct=0.25,only.pos = T)
DEGs_ExN <- DEGs_ExN %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(DEGs_ExN,"/home/iyun1806/MLdata/analysis/results/4_ExN/2_DEGs_ExN.csv",row.names = T)
DEGs_ExN<-read.csv("/home/iyun1806/MLdata/analysis/results/4_ExN/2_DEGs_TTTS.csv")


##########################注释精细分群------------------
#######################1.注释-----------------
TTTS_ExN$num_cluster2<-Idents(TTTS_ExN)
TTTS_ExN <- RenameIdents(TTTS_ExN,
                         "19"="ExN-IPC",
                         "0"="ExN-dev","12"="ExN-dev","13"="ExN-dev",
                         "1"="ExN-L2/3","7"="ExN-L2/3","2"="ExN-L2/3","20"="ExN-L2/3",
                         "4"="ExN-L3/4","6"="ExN-L3/4",
                         "3"="ExN-L4","14"="ExN-L4","16"="ExN-L4","10"="ExN-L4","17"="ExN-L4",
                         "5"="ExN-L4/5","11"="ExN-L4/5",
                         "15"="ExN-L5PT","21"="ExN-L5PT",
                         "8"="ExN-L5/6","23"="ExN-L5/6",
                         "9"="ExN-L6CT","18"="ExN-L6CT",
                         "22"="ExN-L6")
TTTS_ExN$ExN_cluster<-Idents(TTTS_ExN)

###################2.画umap图--------------------------
levels(TTTS_ExN)<-c("ExN-IPC","ExN-dev","ExN-L2/3","ExN-L3/4","ExN-L4","ExN-L4/5","ExN-L5PT","ExN-L5/6","ExN-L6CT","ExN-L6")
col_ExN <- c("#d4dbff","#b8c4f5","#9caef0","#8098e8","#6583dd",
                  "#4f6dcc","#4158b5","#3a449e","#343087","#2e1b70")
names(col_ExN)<-levels(TTTS_ExN)


###画包络线
library(mascarade)
df <- data.frame(Embeddings(TTTS_ExN, "umap"), celltype = TTTS_ExN$ExN_cluster)
write.csv(df,"/home/iyun1806/MLdata/analysis/results/4_ExN/4_df_ExN.csv",row.names = T)
df<-read.csv("/home/iyun1806/MLdata/analysis/results/4_ExN/4_df_ExN.csv")

# 2) 生成每个群的遮罩（轮廓线坐标表）
mask <- generateMask(dims    = Embeddings(TTTS_ExN, "umap"),
                     clusters = TTTS_ExN$ExN_cluster,
                     expand   = 0.01)   # 圈往外扩一点，别贴着点
DimPlot(TTTS_ExN,reduction = "umap",raster = F,label = T,cols = col_ExN) + 
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

ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/4_ExN_cluster_Dimplot.pdf",width = 4.5,height = 4.5)

###########3.画harmony之前的umap图----------
TTTS_ExN <- RunUMAP(TTTS_ExN, reduction = "pca", dims = 1:50,
                    reduction.name = "umap_before_harmony",
                    reduction.key = "UMAPb_")
DimPlot(TTTS_ExN,reduction = "umap_before_harmony",raster = F,label = F,cols = col_ExN) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/4_ExN_cluster_beforeharmony_Dimplot.pdf",width = 6,height = 4.5)

qsave(TTTS_ExN,"/home/iyun1806/MLdata/analysis/results/4_ExN/4_TTTS_ExN_named.qs")
TTTS_ExN<-qread("/home/iyun1806/MLdata/analysis/results/4_ExN/4_TTTS_ExN_named.qs")

###########4.基因Dotplot----------
ExN_genes <- c(
  "EOMES","NEUROG1","NEUROD4","NHLH1",            # ExN-IPC
  "DCX","STMN2","NEUROD6","SOX11",                # ExN-dev
  "CUX2","SATB2","KIF26B",                        # ExN-L2/3
  "CHN2","PRSS12","STK32B",                       # ExN-L3/4
  "RORB","NTNG1","LRRK1","CA10",                  # ExN-L4
  "NTNG2","FOXP1","RORB",                         # ExN-L4/5
  "ETV1","FEZF2","BCL11B","SULF2","CRYM",         # ExN-L5PT
  "TLE4","FOXP2","CRYM","SULF2","PCP4",           # ExN-L5/6
  "TLE4","PDE8B","ABI3BP",                        # ExN-L6CT
  "TLE4","FOXP2","PDE8B","KCTD8"                  # ExN-L6
)
ExN_genes<-unique(ExN_genes)
DotPlot(
  TTTS_ExN,
  features = ExN_genes,
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
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/4_ExN_Dotplot.pdf",width = 18,height = 5.5)


###########4.银河图----------
TTTS_ExN<-qread("/home/iyun1806/MLdata/analysis/results/4_ExN/4_TTTS_ExN_named.qs")

library(ggpointdensity)
dat1 <- FetchData(object=TTTS_ExN, vars=c("umap_1","umap_2","orig.ident"))
colnames(dat1)<-c("UMAP_1","UMAP_2","sample")
head(dat1)

ggplot(data = dat1, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_pointdensity(size=1) +
  scale_color_viridis_c(option="F", alpha = 1) + 
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
    linewidth = 0.8
  ) +
  labs(color="counts")+NoLegend()
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/5_ExN_density.pdf",width = 18.5,height = 8)

dat2 <- FetchData(object=TTTS_ExN, vars=c("umap_1","umap_2","Age/PCW"))
colnames(dat2)<-c("UMAP_1","UMAP_2","Age")
head(dat2)

ggplot(data = dat2, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_pointdensity(size=1) +
  scale_color_viridis_c(option="F", alpha = 1) + 
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
    linewidth = 0.8
  ) +
  labs(color="counts")+NoLegend()
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/5_ExN_density_age.pdf",width = 14.5,height = 4)

###########5.Ro/e----------
######Age=16,20
library(ComplexHeatmap)
library(circlize)
meta20 <- TTTS_ExN@meta.data %>% 
  filter(`Age/PCW` %in% c("16","20")) 
tab_samp <- table(meta20$ExN_cluster, meta20$orig.ident)
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
##只有20.（4，4）
######Age=16
library(ComplexHeatmap)
library(circlize)
meta16 <- TTTS_ExN@meta.data %>% 
  filter(`Age/PCW` == "16") 
tab_samp <- table(meta16$ExN_cluster, meta16$orig.ident)
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

##############################6.统计细胞数量和百分比--------------
library(ggplot2)
library(dplyr)

TTTS_ExN<-qread("/home/iyun1806/MLdata/analysis/results/4_ExN/4_TTTS_ExN_named.qs")
dat<-TTTS_ExN@meta.data[,c(1:12,29)]
# 统计计数
count_by_cluster <- dat %>% group_by(ExN_cluster) %>% summarise(n = n())
count_by_sample  <- dat %>% group_by(orig.ident, ExN_cluster) %>% summarise(n = n())

# 占比（样本内归一化）
pct_by_sample <- count_by_sample %>%
  group_by(orig.ident) %>% mutate(percent = n / sum(n) * 100) %>% ungroup()

# 图1：细胞数柱状图
levels(TTTS_ExN)<-c("ExN-IPC","ExN-dev","ExN-L2/3","ExN-L3/4","ExN-L4","ExN-L4/5","ExN-L5PT","ExN-L5/6","ExN-L6CT","ExN-L6")
col_ExN <- c("#d4dbff","#b8c4f5","#9caef0","#8098e8","#6583dd",
             "#4f6dcc","#4158b5","#3a449e","#343087","#2e1b70")
names(col_ExN)<-levels(TTTS_ExN)
ggplot(count_by_cluster, aes(ExN_cluster, n, fill = ExN_cluster)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = col_ExN) +
  labs(x = "Cell cluster", y = "Cell count") + theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "none"
  )
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/6_cellCount.pdf",width = 4.5,height = 6)

# 图2：百分比堆叠图
ggplot(pct_by_sample, aes(sample, percent, fill = ExN_cluster)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", percent)),
            position = position_stack(vjust = 0.5), color = "white", size = 3) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.9) +
  labs(x = "Sample", y = "Percentage (%)", fill = "Cluster") + theme_bw()

################7.画16，20的umap图----------------
TTTS_ExN<-qread("/home/iyun1806/MLdata/analysis/results/4_ExN/4_TTTS_ExN_named.qs")
library(mascarade)
##用总细胞的包络线
mask <- generateMask(dims    = Embeddings(TTTS_ExN, "umap"),
                     clusters = TTTS_ExN$ExN_cluster,
                     expand   = 0.01)   # 圈往外扩一点，别贴着点


TTTS_age_16<-subset(TTTS_ExN,subset=`Age/PCW` ==16)
col_ExN <- c("#d4dbff","#b8c4f5","#9caef0","#8098e8","#6583dd",
             "#4f6dcc","#4158b5","#3a449e","#343087","#2e1b70")
names(col_ExN)<-levels(TTTS_ExN)



DimPlot(TTTS_age_16,reduction = "umap",raster = F,label = T,cols = col_ExN,repel = T) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/7_ExN_Age16_Dimplot.pdf",width = 4.5,height = 4.5)

TTTS_age_20<-subset(TTTS_ExN,subset=`Age/PCW` ==20)
DimPlot(TTTS_age_20,reduction = "umap",raster = F,label = T,cols = col_ExN,repel = T) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/7_ExN_Age20_Dimplot.pdf",width = 4.5,height = 4.5)

################8.画16，20的细胞计数----------------
dat<-TTTS_age_16@meta.data[,c(1:12,29)]
count_by_cluster <- dat %>% group_by(ExN_cluster) %>% summarise(n = n())

col_ExN <- c("#d4dbff","#b8c4f5","#9caef0","#8098e8","#6583dd",
             "#4f6dcc","#4158b5","#3a449e","#343087","#2e1b70")
ggplot(count_by_cluster, aes(ExN_cluster, n, fill = ExN_cluster)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = col_ExN) +
  labs(x = "", y = "Cell count") + theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "none"
  )
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/8_Age16_cellCount.pdf",width = 5.5,height = 5.5)

dat<-TTTS_age_20@meta.data[,c(1:12,29)]
count_by_cluster <- dat %>% group_by(ExN_cluster) %>% summarise(n = n())

ggplot(count_by_cluster, aes(ExN_cluster, n, fill = ExN_cluster)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = n), vjust = -0.4) +
  scale_fill_manual(values = col_ExN) +
  labs(x = "", y = "Cell count") + theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),                # 去掉网格
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # 横轴斜45度
    legend.position = "none"
  )
ggsave("/home/iyun1806/MLdata/analysis/results/4_ExN/8_Age20_cellCount.pdf",width = 5.5,height = 5.5)




