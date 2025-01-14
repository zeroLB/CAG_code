##Figure1C
p <- DimPlot(sce.all)
ggsave(filename = "/Users/lin/Desktop/胃炎/celltype.tiff", width = 7, height = 5, plot = p)
##Figure1D
Idents(sce.all) <- "sample"
p <- DimPlot(sce.all,split.by = "sample")
ggsave(filename="/Users/lin/Desktop/胃炎/celltype_sample.tiff",plot =p,width = 10,height = 4,dpi = 300)
##Figure1E
marker <- c("CDH1","EPCAM","CD3D","CD3E","CD79A","BANK1","MS4A1","PLVAP","PECAM1","DCN","POSTN","CD74","CD14","CD68","KIT","TPSAB1","RGS5","JCHAIN","MZB1")
p <- DotPlot(sce.all,marker,cols=c("#ffffff","#448444"))+RotatedAxis()
##Figure1F
#提取样本分类信息
cellinfo <- FetchData(sce.all, vars = c("orig.ident", "sample"))
Ratio <- sce.all@meta.data %>% 
  group_by(orig.ident, celltype1) %>% 
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
cli_neu=cli[match(Ratio$Sample,cli$sample),]
Ratio_cli=cbind(Ratio,cli_neu)
head(Ratio_cli)
Ratio$group <- c(rep("N",27),rep("CAG",27))
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
cellname <- unique(Ratio$celltype1)
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[1],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
col <- c("#5CB85C","#D9534F", "#F0AD4E")
p1 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[2],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
p2 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[3],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
p3 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[4],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
p4 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[5],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
p5 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[6],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
p6 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[7],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
p7 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[8],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
p8 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
bar_tmp <- Ratio[Ratio$celltype1%in%cellname[9],]
bar_tmp$group <- c(rep("N",3),rep("CAG",3))
p9 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                x = "group", # X-axis is for groups.
                y = "relative_freq", # Y-axis is for expression levels.
                color = "group", # Fill by sample group.
                fill = NULL,
                #add = "jitter", # Add jitter points.
                bxp.errorbar.width = 0.8,
                width = 0.5,
                size = 0.1,
                font.label = list(size = 20), 
                palette = col)+
  theme(panel.background = element_blank())
ggsave(filename="/Users/lin/Desktop/胃炎/cell_ration_group1.tiff",plot =p1|p2|p3|p4|p5|p6|p7|p8|p9,width = 14,height = 3,dpi = 600)
