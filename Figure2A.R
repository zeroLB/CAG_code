##Figure2A
p1 <- DimPlot(Epi,split.by = "sample")
ggsave(filename="/Users/lin/Desktop/胃炎/epi_UMAP_cell_number_bysample.tiff",plot = p1,width = 8,height = 3,dpi = 300)
##Figure2B
Ratio <- Epi@meta.data %>%group_by(sample,celltype) %>%
  count()%>%
  group_by(sample) %>%
  mutate(Freq = n/sum(n)*100)
p <- ggplot(Ratio, aes(x = sample, y = Freq, fill = celltype))+
  geom_col()+theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                               #"#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                               "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "/Users/lin/Desktop/胃炎/tcell.ratio.tiff", width = 4, height = 5, plot = p)
##Figure2C
p1 <- FeaturePlot(Epi,features = c("MUC5AC","TFF1","MUC6","PGA5","PGA4","ATP4A","ATP4B","CHGA","PLP1","S100B"),cols = c("grey", "#CC0033"),ncol = 4)
ggsave(filename="/Users/lin/Desktop/胃炎/epi-markers-gene-umap.tiff",plot = p1,width = 14,height = 9,dpi = 300)
##Figure2D
library(COSG)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(stringr)
library(COSG)
library(GseaVis)
marker_cosg.smc.sample <- cosg(
  smc,
  groups=c('all'),
  assay='SCT',
  slot='data',
  mu=100,n_genes_user=200)
x_BP.smc = compareCluster(marker_cosg.epi.celltype[[1]], fun='enrichGO', 
                          OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont="BP")
res <- x_BP.smc@compareClusterResult
dt <- enrich
dt <- dt[order(dt$Cluster), ]
dt$Description <- factor(dt$Description, levels = dt$Description)
enrich <- res %>% 
  group_by(Cluster) %>% 
  top_n(n = 10, wt = -pvalue)%>%filter(Cluster %in% c("Control","chronic atrophic gastritis"))
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)
p <- ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description), fill = Cluster)) +
  scale_fill_manual(values =c('#6bb9d2')) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "Chronic atrophic gastritis", title = " Pathway enrichment") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  #geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust = 2.5, color=rep(c('#6bb9d2'),each=5)) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme +
  NoLegend()
ggsave(filename = "/Users/lin/Desktop/胃炎/smcgo.tiff", width = 5.2, height = 5.1, plot = p)
##Figure2E
gene <- c("LCN2","DUOX2")
smg.deg <- FindMarkers(smc, ident.1 = "chronic atrophic gastritis",ident.2 = "control")
smgdeg.antimibrobial <- smg.deg[smg.deg$gene%in%gene,]
smgdeg.antimibrobial %>% 
  filter(avg_log2FC>0) %>% 
  ggplot(aes(x = reorder(gene, avg_log2FC),
             y = avg_log2FC,
             fill = cluster)) +
  geom_bar(stat = "identity",
           color = NA)+scale_fill_manual(values =c("#B17BA6"))+
  coord_flip()+
  #geom_text(aes(label =  gene, color = "white", hjust = 0.1, vjust = 0.5), size = 2.5)+
  theme_classic()->p.up2 
ggsave(filename="/Users/lin/Desktop/胃炎/smc_santimibrobial.tiff",plot = p.up2,
       width = 6,height = 4,dpi = 300)
genes <-c("SMARCC2","BCL7C","CDKN1A")
mncdeg.cellcycle <- mncdeg[mncdeg$gene%in%genes,]
mncdeg.cellcycle$cluster <- "CAG"
mncdeg.cellcycle %>% 
  filter(avg_log2FC>0) %>% 
  ggplot(aes(x = reorder(gene, avg_log2FC),
             y = avg_log2FC,
             fill = cluster)) +
  geom_bar(stat = "identity",
           color = NA)+scale_fill_manual(values =c("#FF9900"))+
  coord_flip()+
  #geom_text(aes(label =  gene, color = "white", hjust = 0.1, vjust = 0.5), size = 2.5)+
  theme_classic()->p.up2 
ggsave(filename="/Users/lin/Desktop/胃炎/mnc_cellcycle.tiff",plot = p.up2,
       width = 4,height = 1.5,dpi = 300)
##Figure2F 
p1 <-VlnPlot(cc,"MT1M",pt.size = 0,  stack = T)+stat_compare_means()
p2 <- VlnPlot(cc,"MT1H",pt.size = 0)+stat_compare_means()
p3 <- VlnPlot(cc,"MT1F",pt.size = 0)+stat_compare_means()
ggsave(filename="/Users/lin/Desktop/胃炎/vlnMT1M.tiff",plot = p1,width = 6,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/vlnMT1H.tiff",plot = p2,width = 6,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/vlnMT1F.tiff",plot = p3,width = 6,height = 5,dpi = 300)
p1 <- VlnPlot(pc,"GPX4",pt.size = 0)+stat_compare_means()
p2 <- VlnPlot(pc,"PRDX3",pt.size = 0)+stat_compare_means()
p3 <- VlnPlot(pc,"PRDX1",pt.size = 0)+stat_compare_means()
ggsave(filename="/Users/lin/Desktop/胃炎/vlnGPX4.tiff",plot = p1,width = 6,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/vlnPRDX3.tiff",plot = p2,width = 6,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/vlnPRDX1.tiff",plot = p3,width = 6,height = 5,dpi = 300)

