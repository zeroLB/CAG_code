##Figure3A
p1  <-  DimPlot(Tcell_remove_mito.1,split.by = "sample")
p2 <-  DimPlot(Tcell_remove_mito.1)
ggsave(filename="/Users/lin/Desktop/胃炎/UMAP_cell_number.tiff",plot = p2,width = 6,height = 3,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/UMAP_cell_number_bysample.tiff",plot = p1,width = 8,height = 3,dpi = 300)
##Figure3B
p1  <-  DotPlot(Tcell_remove_mito.1,features = c("CD4","CD8A"),cols = c("grey","red"))
ggsave(filename="/Users/lin/Desktop/胃炎/CD4-CD8A-gene.tiff",plot = p1,width = 8,height = 3,dpi = 300)
##Figure3C
Ratio <- Tcell_remove_mito.1@meta.data %>%group_by(sample,celltype) %>%
  count() %>%
  group_by(sample) %>%
  mutate(Freq = n/sum(n)*100)

p <- ggplot(Ratio, aes(x = sample, y = Freq, fill = celltype))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                               #"#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                               "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "/Users/lin/Desktop/胃炎/Tcell.ratio.tiff", width = 4, height = 5, plot = p)
##Figure3D
Tcelldeg_cluster <-FindAllMarkers(Tcell_remove_mito.1,only.pos = T)
Tcelldeg_cluster %>% group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> Tcelldeg_cluster3.top5
mean_gene_exp <- AverageExpression(Tcell_remove_mito.1,
                                   features = Tcelldeg_cluster3.top5$gene,
                                   group.by = 'celltype',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()
mean_gene_exp <- mean_gene_exp[,c(4:6)]
# add colnames
colnames(mean_gene_exp) <- c("CD4+ Tex","CD8+ Tef","CD8+ Tsr")
# Z-score
htdf <- t(scale(t(mean_gene_exp),scale = T,center = T))
colnames(htdf) <- c("CD4+ Tex","CD8+ Tef","CD8+ Tsr")
# color
col_fun = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
# anno_col <- pal_npg()(9)
anno_col <- brewer.pal(3, "Paired")
# top annotation
column_ha = HeatmapAnnotation(cluster = colnames(htdf),
                              col = list(cluster = anno_col))
# plot
p1 <- Heatmap(htdf,
              name = "Z-score",
              cluster_columns = F,cluster_rows = F,
              row_title = "Cluster top 5 Marker genes",
              column_title = "Clusters",
              row_names_gp = gpar(fontface = 'italic',fontsize = 10),
              row_names_side = 'left',
              border = T,
              rect_gp = gpar(col = "white", lwd = 1),
              column_names_side = 'top',
              column_names_rot = 45,
              #top_annotation = column_ha,
              # column_split = paste('clsuter ',0:8,sep = ''),
              col = col_fun)
tiff("/Users/lin/Desktop/胃炎/heatmap.tiff",width = 6.5, height = 14, units = "cm",res = 300)
print(p1)
dev.off()
##Figure3E
gene <- c("CTLA4","PDCD1","GZMA","NKG7")
p1 <- FeaturePlot(Tcell_remove_mito.1,features = gene,split.by = "sample",cols = c("grey", "#CC0033"),label.size = 2)
ggsave(filename="/Users/lin/Desktop/胃炎/dfg.tiff",plot = p1,width = 6,height = 10,dpi = 300)
##Figure3F
p <- DotPlot(Tcell,c("CTLA4","GZMA","NKG7","PDCD1"),cols = c("grey", "#CC0033"),split.by = "sample")
df <- p$data
df <- df[c(1:8),]
df$dataset <-c(rep("control",4),rep("CAG",4))
p <-ggplot(df,aes(x=features.plot,y=pct.exp,fill=dataset))+
  geom_bar(stat = "identity",position = position_dodge(0.55),color="black",width = 0.5,linewidth=0.3)+
  theme_bw()+
  xlab("gene")+ylab("cell propotion in T cell")+scale_fill_brewer(palette = 'Set2')+theme(axis.line = element_blank())
ggsave(filename="/Users/lin/Desktop/胃炎/makerTcellratio.tiff",plot = p,width = 6,height = 4,dpi = 300)
##Figure3G
plot_expr <- function(data, group_col, expression_col, marker_col, 
                      plot_type = "box_violin", facet_ncol = 5,
                      comparisons = NULL, diff_method = "wilcox.test", 
                      theme_choice = theme_clean(), legend_position = "none") {
  # packages
  if (!require(ggsci)) {
    install.packages("ggsci")
    library(ggsci)
  }
  
  if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  
  if (!require(ggdist)) {
    install.packages("ggdist")
    library(ggdist)
  }
  
  # Base ggplot
  p <- ggplot(data, aes_string(x = group_col, y = expression_col, fill = group_col)) +
    
    # Plot type options
    switch(plot_type,
           "box_violin" = list(
             geom_violin(trim = TRUE, adjust = 0.95),
             geom_boxplot(outlier.shape = NA, color = "white", width = 0.18)
           ),
           "box_jitter" = list(
             geom_boxplot(outlier.shape = NA),
             geom_jitter(size = 1.5, alpha = 0.55, position = position_jitter(0.22))
           ),
           "halfeye" = ggdist::stat_halfeye(),
           "eye" = ggdist::stat_eye()
    ) +
    
    # Statistical comparison
    stat_compare_means(
      comparisons = comparisons,
      label = "p.signif", 
      method = diff_method, 
      label.x = 1.5, 
      col = "red", 
      size = 3
    ) +
    
    # Facet and labels
    facet_wrap(as.formula(paste("~", marker_col)), scales = "free", ncol = facet_ncol) +
    
    # Theme and colors
    theme_choice +
    ggsci::scale_fill_nejm() +
    ggsci::scale_color_nejm() +
    theme(
      legend.position = legend_position,
      axis.title.y = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}
exhaustion <-list(c("PDCD1","LAYN","HAVCR2","LAG3","CTLA4","TIGIT","TOX","VSIR","BTLA","ENTPD1"))
cytotoxicity <- list(c("GZMA","GZMB","GZMH","GZMK","GZMH","GNLY","PRF1","IFNG","TNF","SERPINB1",
                       "SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7","CSTB",
                       "LAMP1","LAMP3","CAPN2"))
Inscore <- AddModuleScore(Inscore,
                          features = exhaustion, nbin = 12,
                          ctrl = 100,
                          name = "exhaustion")
Inscore <- AddModuleScore(Inscore,
                          features = cytotoxicity, nbin = 12,
                          ctrl = 100,
                          name = "cytotoxicity")
Inscore1 <- subset(Inscore,idents = "CD8+ effector T cell")
df1 <- Inscore@meta.data[,c(20,21)]
df1$symbol <- "exhaustion"
colnames(df1)[2] <- "score"
Idents(Inscore1) <- "sample"
Inscore2 <- subset(Inscore,idents = "CD4+ exhausted T cell")
df2 <- Inscore@meta.data[,c(20,22)]
df2$symbol <- "exhaustion"
colnames(df2)[2] <- "score"
Idents(Inscore2) <- "sample"
plot1 <- plot_expr(data = df1, group_col = "sample", expression_col = "score", marker_col = "symbol", 
                   plot_type = "box_violin",theme_choice = theme_pubclean())+ scale_fill_manual(values = c("#00468BFF","#925E9FFF"))
plot2 <- plot_expr(data = df2, group_col = "sample", expression_col = "score", marker_col = "symbol", 
                   plot_type = "box_violin",theme_choice = theme_pubclean())+ scale_fill_manual(values = c("#00468BFF","#925E9FFF"))
ggsave(filename="/Users/lin/Desktop/胃炎/cytotoxicity-1.tiff",plot = plot1|plot2,width = 4,height = 4,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/exhaution-1.tiff",plot = plot1,width = 4,height = 3,dpi = 300)
