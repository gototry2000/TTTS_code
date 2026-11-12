














###############1没有整合的----------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin.qs")
TTTS_filter_merge <- NormalizeData(TTTS_filter_merge, verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = T) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(npcs = 50, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:50, verbose = T)

# 修改 orig.ident
TTTS_filter_merge$orig.ident <- recode(TTTS_filter_merge$orig.ident,
                                       "ga22" = "ga1",
                                       "ga24" = "ga2",
                                       "ga34" = "ga3",
                                       "hft_w20" = "hft1",
                                       "hft_w21" = "hft2",
                                       "hft_w24" = "hft3",
                                       "hft_w16" = "hft4"
)
qsave(TTTS_filter_merge,"/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_beforeClin_umap.qs")
DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "orig.ident",raster = F,label = F) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/1_origident.pdf",width = 5.3,height = 4.5)

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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/1_AgePCW.pdf",width = 5,height = 4.5)

###############2.harmony的----------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/2_merge_AfterClin.qs")
# 修改 orig.ident
TTTS_filter_merge$orig.ident <- recode(TTTS_filter_merge$orig.ident,
                                       "ga22" = "ga1",
                                       "ga24" = "ga2",
                                       "ga34" = "ga3",
                                       "hft_w20" = "hft1",
                                       "hft_w21" = "hft2",
                                       "hft_w24" = "hft3",
                                       "hft_w16" = "hft4"
)
DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "orig.ident",raster = F,label = F) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/2_Harmony_origident.pdf",width = 5.3,height = 4.5)

###############3.RPCA----------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_RPCA.qs")
# 修改 orig.ident
TTTS_filter_merge$orig.ident <- recode(TTTS_filter_merge$orig.ident,
                                       "ga22" = "ga1",
                                       "ga24" = "ga2",
                                       "ga34" = "ga3",
                                       "hft_w20" = "hft1",
                                       "hft_w21" = "hft2",
                                       "hft_w24" = "hft3",
                                       "hft_w16" = "hft4"
)
DimPlot(TTTS_filter_merge,reduction = "RPCA",
        group.by = "orig.ident",raster = F,label = F) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/3_RPCA_origident.pdf",width = 5.3,height = 4.5)



###############################4.bbknn-----------------------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/1_merge_bbknnR.qs")

# 修改 orig.ident
TTTS_filter_merge$orig.ident <- recode(TTTS_filter_merge$orig.ident,
                                       "ga22" = "ga1",
                                       "ga24" = "ga2",
                                       "ga34" = "ga3",
                                       "hft_w20" = "hft1",
                                       "hft_w21" = "hft2",
                                       "hft_w24" = "hft3",
                                       "hft_w16" = "hft4"
)
DimPlot(TTTS_filter_merge,reduction = "umap",
        group.by = "orig.ident",raster = F,label = F) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/4_bbknn_origident.pdf",width = 5.3,height = 4.5)


############################5.scVI---------------------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/2_TTTS_scvi_UMAP.qs")
# 修改 orig.ident
TTTS_filter_merge$orig.ident <- recode(TTTS_filter_merge$orig.ident,
                                       "ga22" = "ga1",
                                       "ga24" = "ga2",
                                       "ga34" = "ga3",
                                       "hft_w20" = "hft1",
                                       "hft_w21" = "hft2",
                                       "hft_w24" = "hft3",
                                       "hft_w16" = "hft4"
)
DimPlot(TTTS_filter_merge,reduction = "scvi",
        group.by = "orig.ident",raster = F,label = F) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/5_scvi_origident.pdf",width = 5.3,height = 4.5)

############################6.scanvi-------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_scANVI/4_TTTS_scanvi_UMAP.qs")
# 修改 orig.ident
TTTS_filter_merge$orig.ident <- recode(TTTS_filter_merge$orig.ident,
                                       "ga22" = "ga1",
                                       "ga24" = "ga2",
                                       "ga34" = "ga3",
                                       "hft_w20" = "hft1",
                                       "hft_w21" = "hft2",
                                       "hft_w24" = "hft3",
                                       "hft_w16" = "hft4"
)
DimPlot(TTTS_filter_merge,reduction = "scanvi",
        group.by = "orig.ident",raster = F,label = F) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/6_scanvi_origident.pdf",width = 5.3,height = 4.5)

############################7,Liger-------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_Liger/2_TTTS_ligerUmap.qs")
# 修改 orig.ident
TTTS_filter_merge$orig.ident <- recode(TTTS_filter_merge$orig.ident,
                                       "ga22" = "ga1",
                                       "ga24" = "ga2",
                                       "ga34" = "ga3",
                                       "hft_w20" = "hft1",
                                       "hft_w21" = "hft2",
                                       "hft_w24" = "hft3",
                                       "hft_w16" = "hft4"
)
DimPlot(TTTS_filter_merge,reduction = "umap_liger",
        group.by = "orig.ident",raster = F,label = F) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/7_liger_origident.pdf",width = 5.3,height = 4.5)

############################8.fastMNN-------------
TTTS_filter_merge<-qread("/home/iyun1806/MLdata/analysis/results/2_harmony/4_FastMNN/1_TTTS_FastMNN.qs")
TTTS_filter_merge$orig.ident <- recode(TTTS_filter_merge$orig.ident,
                                       "ga22" = "ga1",
                                       "ga24" = "ga2",
                                       "ga34" = "ga3",
                                       "hft_w20" = "hft1",
                                       "hft_w21" = "hft2",
                                       "hft_w24" = "hft3",
                                       "hft_w16" = "hft4"
)
DimPlot(TTTS_filter_merge,reduction = "mnn",
        group.by = "orig.ident",raster = F,label = F) + 
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
ggsave("/home/iyun1806/MLdata/analysis/results/2_harmony_plot/8_fastmnn_origident.pdf",width = 5.3,height = 4.5)
























