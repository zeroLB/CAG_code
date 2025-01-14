##Figure4A
mycolor =c("#89558D","#435B95","#79B99D")
p <- DimPlot(myeloid,cols =mycolor,split.by = "sample")
ggsave(filename="/Users/lin/Desktop/胃炎/macrophage-cluster.tiff",plot = p,width = 8,height = 4,dpi = 300)
##Figure4B
Idents(myeloid) <- "sample"
marker_cosg.myeloid.sample <- cosg(
  myeloid,
  groups=c('all'),
  assay='SCT',
  slot='data',
  mu=100,
  n_genes_user=200)
x_BP.myeloid.sample = compareCluster(marker_cosg.myeloid.sample[[1]], fun='enrichGO', 
                                     OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont="BP")
go <- x_BP.myeloid.sample@compareClusterResult$Description[1:15]
p<-dotplot(x_BP.myeloid.sample, showCategory=go)+theme(axis.text.x = element_text(angle=45, hjust=1),
                                                       axis.text.y = element_text(size=12),
                                                       panel.spacing = unit(5, "mm"))+
  scale_colour_gradientn(colours =c("#E54924","#498EA4"))
##Figure4C
p <- FeaturePlot(myeloid,c("C1QA","C1QB","C1QC","CD1C","FCN1","S100A8","S100A9"),cols = c("grey", "#CC0033"),ncol = 4)
ggsave(filename="/Users/lin/Desktop/胃炎/macrophage-related-gene.tiff",plot = p,width = 12,height = 5,dpi = 300)
##Figure4D
p1 <- FeaturePlot(myeloid.C1Q,c("MRC1","MSR1"),
                  split.by = "sample",cols = c("grey", "#CC0033"),ncol = 4)
p2 <- FeaturePlot(myeloid.C1Q,c("CD36","FPR1"),
                  split.by = "sample",cols = c("grey", "#CC0033"),ncol = 4)
ggsave(filename="/Users/lin/Desktop/胃炎/macrophage-related-gene-1.tiff",plot = p1,width = 6,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/macrophage-related-gene-2.tiff",plot = p2,width = 6,height = 5,dpi = 300)
##Figure4E
myeloid.C1Q <- subset(myeloid,idents = "C1Q+ Macrophage")
p <- DotPlot(myeloid.C1Q,c("MRC1","CD36","MSR1","FPR1"),cols = c("grey", "#CC0033"),split.by = "sample")
df <- p$data
colnames(df)[4] <- "dataset"
df$dataset <-c(rep("control",4),rep("CAG",4))

p <-ggplot(df,aes(x=features.plot,y=pct.exp,fill=dataset))+
  geom_bar(stat = "identity",position = position_dodge(0.55),color="black",width = 0.5,linewidth=0.3)+
  theme_bw()+
  xlab("gene")+ylab("cell proportion in C1Q+ Macrophage cells")+scale_fill_brewer(palette = 'Set2')+theme(axis.line = element_blank())+coord_flip()
ggsave(filename="/Users/lin/Desktop/胃炎/MRC1macrophagecellratio1.tiff",plot = p,width = 4,height = 4,dpi = 300)
##Figure4F
p1 <- DimPlot(sce.all.wy.cell.annotation,cells.highlight = WhichCells(sce.all.wy.cell.annotation,idents = "Myeloid cell"),
              cols.highlight = "orange",cols = "grey",split.by = "sample")
ggsave(filename="/Users/lin/Desktop/胃炎/mastcell.tiff",plot = p1,width = 8,height = 4,dpi = 300)
##Figure4G
marker_cosg <- cosg(
  mastcell,
  groups=c('all'),
  assay='SCT',
  slot='data',
  mu=100,
  n_genes_user=200)
x_BP = compareCluster(marker_cosg[[1]], fun='enrichGO', 
                      OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont="BP")
res <-x_BP@compareClusterResult
enrich <- res[res$Description%in%c("cellular response to extracellular stimulus",
                                   "cellular response to external stimulus","regulation of TOR signaling",
                                   "positive regulation of pattern recognition receptor signaling pathway",
                                   "myeloid cell differentiation","immune response-regulating signaling pathway"),]
dt <- enrich
dt <- dt[order(dt$Cluster),]
dt$Description <- factor(dt$Description, levels = dt$Description)
colnames(dt)
# 先自定义主题：
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
ggsave(filename = "/Users/lin/Desktop/胃炎/mastgo.tiff", width = 5.2, height = 5.1, plot = p)
##Figure4H
p1 <- VlnPlot(mastcell,features = c("IL1RL1"),split.by = "sample",cols = c("#0099CC", "#CC0033"))+stat_compare_means()
p2 <- VlnPlot(mastcell,features = c("MAPK1"),split.by = "sample",cols = c("#0099CC", "#CC0033"))+stat_compare_means()
p3 <- VlnPlot(mastcell,features = c("POU2F1"),split.by = "sample",cols = c("#0099CC", "#CC0033"))+stat_compare_means()
p4 <- VlnPlot(mastcell,features = c("LYN"),split.by = "sample",cols = c("#0099CC", "#CC0033"))+stat_compare_means()
ggsave(filename="/Users/lin/Desktop/胃炎/mast-diffgene1.tiff",plot = p1|p2,width = 12,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/mast-diffgene2.tiff",plot = p3|p4,width = 12,height = 5,dpi = 300)



