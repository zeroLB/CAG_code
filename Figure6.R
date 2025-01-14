##Figure6A
sam1 <- subset(INF,idents = c("XHI001I","XHWS001"))
degsam1 <- FindAllMarkers(sam1)
degsam1 <- degsam1[degsam1$p_val_adj < 0.05,]
degsam1 <- degsam1[degsam1$avg_log2FC > 0,]
degsam1 <- degsam1[degsam1$cluster == "XHWS001",]
sam2 <- subset(INF,idents = c("XHI002I-1","XHWS002"))
degsam2 <- FindAllMarkers(sam2)
degsam2 <- degsam2[degsam2$p_val_adj < 0.05,]
degsam2 <- degsam2[degsam2$avg_log2FC > 0,]
degsam2 <- degsam2[degsam2$cluster == "XHWS002",]
sam3 <- subset(INF,idents = c("XHI003I-1","XHWS003-1"))
degsam3 <- FindAllMarkers(sam3)
degsam3 <- degsam3[degsam3$p_val_adj < 0.05,]
degsam3 <- degsam3[degsam3$avg_log2FC > 0,]
degsam3 <- degsam3[degsam3$cluster == "XHWS003-1",]
patient1 <- degsam1$gene
patient2 <- degsam2$gene
patient3 <- degsam3$gene
venn_list <- list(patient1=S1, patient2=S2, patient3=S3)
p <- ggvenn(venn_list, # 数据列表
            columns = c("patient1","patient2","patient3"), # 指定需要比较的数据集
            show_percentage = F, # 显示每一组的百分比
            digits = 1, # 百分比的小数点位数
            fill_color = c("#E41A1C", "#1E90FF", "#FF8C00"))
ggsave("/Users/lin/Desktop/胃炎Fig/venn.tiff",plot = p,width = 8, height = 4,dpi = 300)
##Figure6B
p<-FeaturePlot(object = IAF, features = "ASPN",cols = c("grey","#CC0033"),split.by = "sample")
ggsave(filename = "/Users/lin/Desktop/胃炎Fig/ASPNfeature.tiff",plot = p,width = 8,height = 4,dpi = 300)
##Figure6C
gastritis3 <- read.delim2("/Users/lin/Desktop/胃癌课题/胃癌公开数据/GSVA_score/已算-GSE130823_GCnormalized_data.txt")
ASPN1 <- gastritis3[gastritis3$X..Notes...All.Entities%in%c("ProbeName","A_23_P216429"),]
ASPN1 <- t(ASPN1)
colnames(ASPN1) <- ASPN1[1,]
ASPN1 <- ASPN1[-1,]
ASPN1 <- as.data.frame(ASPN1)
group <- Gsva$stage
group <- rep(c("Gastritis","LGIN"),47)
ASPN1$group <- group
colnames(ASPN1)[2] <- "value"
ASPN1$value <- as.numeric(ASPN1$value)
my_comparisons <- list(c("gastritis","LGIN"))
my_comparisons1 <- list(c("gastritis","HGIN"))
my_comparisons2 <- list(c("gastritis","GC"))

p <- ggplot(ASPN1, aes(x=group, y=value, fill=group)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  ggtitle("ASPN") +
  labs(fill="Category")+
  xlab("")+
  ylab("")+
  theme_bw()+stat_compare_means(method ="wilcox.test")
ggsave("/Users/lin/Desktop/胃炎Fig/ASPNboxplot-1.tiff",plot = p,width = 4, height = 4,dpi = 300)
p <- ggplot(ASPN1, aes(group, value)) + 
  geom_boxplot(aes(fill = group)) + 
  geom_signif(comparisons = my_comparisons)+geom_signif(comparisons = my_comparisons1,
                                                        y_position = 4.2)+geom_signif(comparisons = my_comparisons2,
                                                                                      y_position = 5.2)+theme_bw()+xlab("")+ ggtitle("ASPN")
ylab("")
ggsave("/Users/lin/Desktop/胃炎Fig/ASPNboxplot-1.tiff",plot = p,width = 4, height = 4,dpi = 300)

ASPN2 <- gastritis1[gastritis1$gene_name%in%"ASPN",]
ASPN2 <- ASPN2[,c(6,14:43)]
ASPN2 <- as.data.frame(t(ASPN2))
group <- c(rep("GC",10),rep("IM",10),rep("NAG",10))
ASPN2 <- as.data.frame(ASPN2)
ASPN2$group <- group
ASPN2<- ASPN2[-1,]
ASPN2$ASPN2 <- as.numeric(ASPN2$ASPN2)
ASPNdat <- ASPN2[ASPN2$group%in%c("GC","NAG"),]
p<- ggplot(ASPNdat, aes(x=group, y=ASPN2, fill=group)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  ggtitle("ASPN") +
  labs(fill="Category")+
  xlab("")+
  ylab("")+
  theme_bw()+stat_compare_means(method ="wilcox.test")
ggsave(filename = "/Users/lin/Desktop/胃炎Fig/ASPNboxplot-2.tiff",plot = p,width = 4, height = 4,dpi = 300)

TCGA <- readRDS("/Users/lin/Desktop/TCGA/pancancer_drawdata.rds")
STAD <- drawdata[drawdata$Cancer=="STAD",]
p2 <- ggplot(STAD, aes(x=Type, y=ASPN, fill=Type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  ggtitle("ASPN") +
  labs(fill="Category")+
  xlab("")+
  ylab("")+
  theme_bw()+stat_compare_means(method ="wilcox.test")
ggsave("/Users/lin/Desktop/胃炎Fig/ASPNboxplotTCGA.tiff",plot = p2,width = 4.5, height = 4,dpi = 300)
##Figure6D
gastritis <- read_excel("/Users/lin/Downloads/GSE153224_mRNA_Expression_Profiling.xlsx")
ASPN <- gastritis[gastritis$gene_short_name%in%"ASPN",]
ASPN <- ASPN[,c(7:16)]
ASPN <- as.data.frame(t(ASPN))
ASPN$group <- c(rep("control",5),rep("chronic atrophic gastritis",5))
colnames(ASPN) <- c("value","name")
p <- ggplot(ASPN, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  ggtitle("ASPN") +
  labs(fill="Category")+
  xlab("")+
  ylab("")+
  theme_bw()+stat_compare_means(method ="wilcox.test")
##Figure6F
Idents(fibro) <- "celltype"
IF <- subset(fibro,ident = "CCL11+APOE+ fibroblasts")
Idents(IF) <- "sample"
IF.CAG  <- subset(IF, ident = "chronic atrophic gastritis")
exprss <-as.data.frame(IF.CAG@assays$RNA@data)
ASPN <- exprss[rownames(exprss)%in%("ASPN"),]
ASPN <- as.data.frame(t(ASPN))
gene_expression_sel <- subset(ASPN,ASPN > 0)
positive <- rownames(gene_expression_sel) 
IF.CAG$ASPN <- ifelse(rownames(IF.CAG@meta.data)%in% positive, yes = "positive", no = "negative")
Idents(IF.CAG) <- "ASPN"
dif <- FindMarkers(IF.CAG,ident.1 = "positive",ident.2 ="negative")
dif=data.frame(
  symbol=rownames(deg_all),
  log2FoldChange=deg_all$avg_log2FC,
  padj=deg_all$p_val_adj
)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
dif$logpvalue <- -log10(dif$padj)
# 标记上下调的差异基因
dif$group ="ns"
dif$group[which((dif$avg_log2FC > 0.5)&(dif$p_val_adj < 0.05))] = "up"   
#标记上调基因
dif$group[which((dif$avg_log2FC < -0.5) & (dif$p_val_adj < 0.05))] = "down"  
#标记下调基因
dif$gene <- rownames(dif)
dif$logpvalue <- -log10(dif$p_val_adj)
# 设置颜色并排序
dif$color <- ifelse(dif$group == "none" & dif$gene == "", "color1", 
                    ifelse(dif$group == "up" & dif$gene == "", "color2",  
                           ifelse(dif$group == "down" & dif$gene == "", "color3", 
                                  ifelse(dif$group == "up" & dif$gene != "", "color4", "color5")))) 
dif <- arrange(dif, color)  #dplyr里的函数  #这里排序是为了稍后颜色color和palette对应上
tb2=table(dif$threshold); print(tb2)
cols=c("#497aa2", "#ae3137")
g1 = ggplot(data=dif, aes(x=avg_log2FC, y=logpvalue, color=group))+scale_color_manual(values=c("#497aa2", "grey","#ae3137"))+
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype=2, color="blue")+
  geom_hline(yintercept = -log10(0.05), linetype=2, color="blue")+
  labs(paste("DEG:", "DSS vs WT"))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*P.adj))+
  theme_classic(base_size = 14) +
  theme(legend.box = "horizontal",
        legend.position="top",
        legend.spacing.x = unit(0, 'pt'),
        legend.text = element_text( margin = margin(r = 20) ),
        legend.margin=margin(b= -10, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size=10))
dif.sig=dif[which(dif$group != "ns" ), ]
label.symbols=rownames(dif.sig)
dd_text = dif[label.symbols, ]
g2 <- g1 + geom_text_repel(data=dd_text, aes(x=avg_log2FC, y=logpvalue, label=row.names(dd_text)),
                           colour="black",alpha=1)+labs(title=paste("DEG:", "ASPN gain"))
ggsave(filename="/Users/lin/Desktop/胃炎Fig/ASPNdeg.tiff",plot=g2,width = 4,height= 4,dpi= 600)
