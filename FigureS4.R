##FigureS4A
p1 <- DoHeatmap(fibro,degfibro.top10$gene,size = 3.5,angle = 45) +
  scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
ggsave(filename="/Users/lin/Desktop/胃炎/fibro_cell_cluster.tiff",plot = p1,width = 12,height = 10,dpi = 300)
##FigureS4B
Ratio <- fibro@meta.data %>%group_by(sample,celltype) %>%
  count() %>%
  group_by(sample) %>%
  mutate(Freq = n/sum(n)*100)
Ratio$sample <- c(rep("CAG",4),rep("N",4))
p <- ggplot(Ratio, aes(x = sample, y = Freq, fill = celltype))+
  geom_col()+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                               #"#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                               "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "/Users/lin/Desktop/胃炎/fibrocell.ratio.tiff", width = 4, height = 5, plot = p)
##FigureS4C
Idents(fibro) <- "sample"
fibro.CAG <- subset(fibro, idents = 
                      "chronic atrophic gastritis")
mean_gene_exp <- AverageExpression(fibro,
                                   features = c(chemokine,cytokine),
                                   group.by = 'celltype',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()
mean_gene_exp <- mean_gene_exp[,c(5:8)]
# add colnames
colnames(mean_gene_exp) <- c("CCL11+APOE+ fibroblasts","CXCL14+SOX6+ fibroblasts",
                             "HIPP+ myofibroblasts","RGS5+ pericytes")
# Z-score
htdf <- t(scale(t(mean_gene_exp),scale = T,center = T))
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
tiff("/Users/lin/Desktop/胃炎/fibroheatmap.tiff",width = 10, height = 14, units = "cm",res = 300)
print(p1)
dev.off()
##FigureS4G
p <-FeaturePlot(Tcell,"CCR5",split.by = "sample",cols = c("grey", "#CC0033"))
ggsave(filename = "/Users/lin/Desktop/胃炎/CCR5.tiff", width = 5, height = 3, plot = p)
