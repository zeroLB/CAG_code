###FigureS3A
PPR <- list(c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR10"))
genelist <- list(c("CXCL8","CXCL12","CXCL3","CXCL12","CXCL16","CXCL1","CXCL9","CXCL10", "CCL3","CCL2","CCL4","CCL8","CCL5","CCL20","CCL19","CCL14","CCL18","CXCL12","CXCL3","CXCL12","CXCL16","CXCL1","CXCL9","CXCL10"))
myeloid.C1Q <- subset(myeloid,idents = "C1Q+ macrophages")
Inscore <- AddModuleScore(myeloid.C1Q,
                          features = PPR,nbin = 12,
                          ctrl = 100,
                          name = "Toll like receptors")
Inscore <- AddModuleScore(Inscore,
                          features = genelist,nbin = 12,
                          ctrl = 100,
                          name = "chemokines")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[25] <- 'Toll like receptors' 
colnames(Inscore@meta.data)[26] <- 'chemokines' 
p1 <- VlnPlot(Inscore,features = c('Toll like receptors'),
              pt.size = 0.05, adjust = 2,group.by = "sample")+stat_compare_means()
p2 <- VlnPlot(Inscore,features =  "chemokines", 
              pt.size = 0.05, adjust = 2,group.by = "sample")+stat_compare_means()
gsave(filename="/Users/lin/Desktop/胃炎/Toll like receptors.tiff",plot = p1,width = 6,height = 6,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/chemokines_score.tiff",plot = p2,width = 6,height = 6,dpi = 300)
###FigureS3B
p <- DotPlot(myeloid.C1Q,c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR10"),cols = c("grey", "#CC0033"),split.by = "sample")
df <- p$data
colnames(df)[4] <- "dataset"
df$dataset <-c(rep("control",9),rep("CAG",9))

p <-ggplot(df,aes(x=features.plot,y=pct.exp,fill=dataset))+
  geom_bar(stat = "identity",position = position_dodge(0.55),color="black",width = 0.5,linewidth=0.3)+
  theme_bw()+
  xlab("gene")+ylab("cell propotion in myeloid cells")+scale_fill_brewer(palette = 'Set2')+theme(axis.line = element_blank())
ggsave(filename="/Users/lin/Desktop/胃炎/TLRmacrophagecellratio1.tiff",plot = p,width = 4,height = 6,dpi = 300)
###FigureS3C
p1<-FeaturePlot(myeloid.C1Q,c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR10"),cols = c("grey", "#CC0033"),ncol = 5)
ggsave(filename="/Users/lin/Desktop/胃炎/TLR-related-gene-1.tiff",plot = p1,width = 14,height = 5,dpi = 300)
###FigureS3D
p1<-FeaturePlot(myeloid.C1Q.CAG,c("MMP2","MMP9","MMP12","MMP14","MMP19","MMP25","TIMP1","TIMP2"),cols = c("grey", "#CC0033"),ncol = 5)
ggsave(filename="/Users/lin/Desktop/胃炎/MMP-related-gene-1.tiff",plot = p1,width = 14,height = 5,dpi = 300)
###FigureS3E&F
countexp.Seurat <- scMetabolism::sc.metabolism.Seurat(obj = myeloid.C1Q,method = "VISION",imputation = F, ncores = 1,
                                                      metabolism.type = "KEGG")
saveRDS(countexp.Seurat,"/Users/lin/Downloads/myeloid.metabolism.rds")
input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:30]
p <- scMetabolism::DotPlot.metabolism(obj= countexp.Seurat,pathway = input.pathway,phenotype= "sample",norm = "y")
exprSet <- as.data.frame(as.matrix(myeloid.C1Q@assays[["SCT"]]@data))
exprSet<- as.data.frame(t(exprSet))
exprSet1 <- exprSet[exprSet$FPR1 >0,]
exprSet2 <- exprSet[exprSet$MRC1 >0,]
exprSet3 <- exprSet[exprSet$CD36 >0,]
exprSet4 <- exprSet[exprSet$MSR1 >0,]
countexp.Seurat@meta.data$CD36  <-rownames(countexp.Seurat@meta.data)%in%rownames(exprSet3)
countexp.Seurat@meta.data$CD36[which(countexp.Seurat@meta.data$CD36 %in% c('TRUE'))] = 'CD36+ macrophage'
countexp.Seurat@meta.data$CD36[which(countexp.Seurat@meta.data$CD36 %in% c('FALSE'))] = 'CD36- macrophage'
countexp.Seurat@meta.data$FPR1  <-rownames(countexp.Seurat@meta.data)%in%rownames(exprSet1)
countexp.Seurat@meta.data$FPR1[which(countexp.Seurat@meta.data$FPR1 %in% c('TRUE'))] = 'FPR1+ macrophage'
countexp.Seurat@meta.data$FPR1[which(countexp.Seurat@meta.data$FPR1 %in% c('FALSE'))] = 'FPR1- macrophage'
countexp.Seurat@meta.data$MRC1  <-rownames(countexp.Seurat@meta.data)%in%rownames(exprSet2)
countexp.Seurat@meta.data$MRC1[which(countexp.Seurat@meta.data$MRC1 %in% c('TRUE'))] = 'MRC1+ macrophage'
countexp.Seurat@meta.data$MRC1[which(countexp.Seurat@meta.data$MRC1 %in% c('FALSE'))] = 'MRC1- macrophage'
countexp.Seurat@meta.data$MSR1  <-rownames(countexp.Seurat@meta.data)%in%rownames(exprSet4)
countexp.Seurat@meta.data$MSR1[which(countexp.Seurat@meta.data$MSR1 %in% c('TRUE'))] = 'MSR1+ macrophage'
countexp.Seurat@meta.data$MSR1[which(countexp.Seurat@meta.data$MSR1 %in% c('FALSE'))] = 'MSR1- macrophage'
p4 <- BoxPlot.metabolism(obj = countexp.Seurat.CAG, 
                         pathway = "Nitrogen metabolism", 
                         phenotype = "MRC1", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
p5 <- BoxPlot.metabolism(obj = countexp.Seurat.CAG, 
                         pathway = "Nitrogen metabolism", 
                         phenotype = "FPR1", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
p6 <- BoxPlot.metabolism(obj = countexp.Seurat.CAG, 
                         pathway = "Nitrogen metabolism", 
                         phenotype = "MSR1", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
p7 <- BoxPlot.metabolism(obj = countexp.Seurat.CAG, 
                         pathway = "Nitrogen metabolism", 
                         phenotype = "CD36", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
gsave(filename="/Users/lin/Desktop/胃炎/macrophageNitrogenmetabolism.tiff",plot = p4|p5,width = 10,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/macrophageNitrogenmetabolism-2.tiff",plot = p6|p7,width = 10,height = 5,dpi = 300)

p8 <- BoxPlot.metabolism(obj = countexp.Seurat.CAG, 
                         pathway = "Glycerophospholipid metabolism", 
                         phenotype = "MRC1", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
p9 <- BoxPlot.metabolism(obj = countexp.Seurat.CAG, 
                         pathway = "Glycerophospholipid metabolism", 
                         phenotype = "FPR1", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
p10 <- BoxPlot.metabolism(obj = countexp.Seurat.CAG, 
                          pathway = "Glycerophospholipid metabolism", 
                          phenotype = "MSR1", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
p11 <- BoxPlot.metabolism(obj = countexp.Seurat.CAG, 
                          pathway = "Glycerophospholipid metabolism", 
                          phenotype = "CD36", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
ggsave(filename="/Users/lin/Desktop/胃炎/macrophageGlycerophospholipid-1.tiff",plot = p8|p9,width = 10,height = 5,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/macrophageGlycerophospholipid-2.tiff",plot = p10|p11,width = 10,height = 5,dpi = 300)
p12 <- BoxPlot.metabolism(obj = countexp.Seurat, 
                          pathway = "Glycerophospholipid metabolism", 
                          phenotype = "sample", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
p13 <- BoxPlot.metabolism(obj = countexp.Seurat, 
                          pathway = "Nitrogen metabolism", 
                          phenotype = "sample", ncol = 5)+stat_compare_means(method = "wilcox.test", hide.ns = FALSE,label = "p.signif",vjust=0.02,bracket.size=0.6)
ggsave(filename="/Users/lin/Desktop/胃炎/macrophageGlycerophospholipid-3.tiff",plot = p12|p13,width = 10,height = 5,dpi = 300)

###FigureS3G
p1 <- FeaturePlot(myeloid.C1Q,c("LDLR"),split.by = "sample",cols = c("grey", "#CC0033"))
ggsave(filename="/Users/lin/Desktop/胃炎/LDLR.tiff",plot = p1,width = 8,height = 4,dpi = 300)
