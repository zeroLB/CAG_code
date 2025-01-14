##FigureS1A&B
gene <- read.delim("/Users/lin/Downloads/易感基因geneset 2.csv",sep = ",")
genelist <-  list(gene$Symbol)
Inscore1 <- AddModuleScore(epi,
                           features = genelist,nbin = 12,
                           ctrl = 100,
                           name = "diseasegene")
p <- VlnPlot(Inscore,"diseasegene1")
p <- FeaturePlot(Inscore,"susceptible gene",cols = c("grey", "#CC0033"))
ggsave(filename="/Users/lin/Desktop/胃炎/yigangenescore-3.tiff",plot = p,width = 5,height = 4,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/yigangenescore-4.tiff",plot = p,width = 8,height = 4.5,dpi = 300)
##FigureS1C
df <- Inscore1@meta.data[,c(19,20,21)]
df$Symbol <- "Disease related gene"
colnames(df)[3] <-"score"
df1 <- df[df$celltype%in%"Surface mucous cell",]
df2 <- df[df$celltype%in%"Parietal cell",]
df3 <- df[df$celltype%in%"Mucous neck cell",]
df4 <- df[df$celltype%in%"Gial cell",]
df5 <- df[df$celltype%in%"Enteroendocrine cell",]
df6 <- df[df$celltype%in%"Chief cell",]

plot1 <- plot_expr(data = df1, group_col = "sample", expression_col = "score", marker_col = "celltype", 
                   plot_type = "box_violin",theme_choice = theme_pubclean())+ scale_fill_manual(values = c("#00468BFF","#925E9FFF"))

plot2 <- plot_expr(data = df2, group_col = "sample", expression_col = "score", marker_col = "celltype", 
                   plot_type = "box_violin",theme_choice = theme_pubclean())+ scale_fill_manual(values = c("#00468BFF","#925E9FFF"))

plot3 <- plot_expr(data = df3, group_col = "sample", expression_col = "score", marker_col = "celltype", 
                   plot_type = "box_violin",theme_choice = theme_pubclean())+ scale_fill_manual(values = c("#00468BFF","#925E9FFF"))

plot4 <- plot_expr(data = df4, group_col = "sample", expression_col = "score", marker_col = "celltype", 
                   plot_type = "box_violin",theme_choice = theme_pubclean())+ scale_fill_manual(values = c("#00468BFF","#925E9FFF"))

plot5 <- plot_expr(data = df5, group_col = "sample", expression_col = "score", marker_col = "celltype", 
                   plot_type = "box_violin",theme_choice = theme_pubclean())+ scale_fill_manual(values = c("#00468BFF","#925E9FFF"))

plot6 <- plot_expr(data = df6, group_col = "sample", expression_col = "score", marker_col = "celltype", 
                   plot_type = "box_violin",theme_choice = theme_pubclean())+ scale_fill_manual(values = c("#00468BFF","#925E9FFF"))
ggsave(filename="/Users/lin/Desktop/胃炎/yigangenescore-1.tiff",plot = plot1|plot2|plot3,width = 12,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/yigangenescore-2.tiff",plot = plot4|plot5|plot6,width = 12,height = 5,dpi = 300)
##FigureS1D
genelist <- read_excel("/Users/lin/Downloads/41467_2024_52615_MOESM10_ESM.xlsx")
genelist <- genelist[-1,]
colnames(genelist) <- genelist[1,]
cycle <- list(genelist$Cycle)
epi<- AddModuleScore(epi, features = cycle,
                     ctrl = 100,
                     name = "Cycle")
library(ggpubr)
p <-ggboxplot(epi@meta.data, 
              x="celltype", y="Cycle1", 
              width = 0.4,color = "black",#轮廓颜色
              fill="sample",#填充 
              xlab = F, #不显示x轴的标签 
              bxp.errorbar=T,#显示误差条         
              bxp.errorbar.width=0.05, #误差条大小          
              size=0.5, #箱型图边线的粗细
              palette = "npg", 
              legend = "right")+stat_compare_means(aes(group=sample),
                                                   method="wilcox.test",
                                                   symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                                                   label = "p.signif")
ggsave(filename="/Users/lin/Desktop/胃炎Fig/Cycle.tiff",plot=p,width =12,height=3,dpi=600)
pEMT <- list(genelist$pEMT)
epi<- AddModuleScore(epi, features = pEMT,
                     ctrl = 100,
                     name = "pEMT")
p <- ggviolin(epi@meta.data, 
              x="celltype", y="pEMT1", 
              width = 0.8,color = "black",#轮廓颜色
              fill="sample",#填充 
              xlab = F, #不显示x轴的标签 
              add = 'mean_sd', 
              bxp.errorbar=T,#显示误差条         
              bxp.errorbar.width=0.05, #误差条大小          
              size=0.5, #箱型图边线的粗细
              palette = "npg", 
              legend = "right")
ggsave(filename="/Users/lin/Desktop/胃炎Fig/pEMT.tiff",plot=p,width =12,height=3,dpi=600)
