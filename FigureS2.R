###Figure2A,B&C
p1<-DotPlot(Tcell_remove_mito.1,features = c("NKTR","KLRG1","FGFBP2","FCGR3A","IFNG","NKG7","PRF1","CX3CR1","GZMA","GZMB","GZMK","GZMH","GNLY","CD8A"),
            cols =  c("blue", "yellow"))
pggsave(filename="/Users/lin/Desktop/胃炎/CD8+related-gene.tiff",plot = p1,width = 12,height = 4,dpi = 300)

p1<-DotPlot(Tcell_remove_mito.1,features = c("HSP90AA1","HSPH1","HSPD1","HSPA1B","HSP90AB1","CD8A"),
            cols =  c("blue", "yellow"))

ggsave(filename="/Users/lin/Desktop/胃炎/CD8+stressresponse-gene.tiff",plot = p1,width = 16,height = 4,dpi = 300)

p1<-DotPlot(Tcell_remove_mito.1,features = c("CD4","FOXP3","PDCD1","CTLA4","ICOS","TIGIT","LAG3","IL2RA"),
            cols =  c("blue", "yellow"))
ggsave(filename="/Users/lin/Desktop/胃炎/CD4+exhausted-gene.tiff",plot = p1,width = 12,height = 4,dpi = 300)

##FigureS2D&E
## We load the required packages
library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## We compute Viper Scores 
Tcell <- run_viper(Tcell, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = FALSE))
## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = Tcell) <- "dorothea"
## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(Tcell, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

## We create a data frame containing the cells and their clusters
Idents(epi) <- "celltype"
CellsClusters <- data.frame(cell = names(Idents(Tcell)), 
                            cell_type = as.character(Idents(Tcell)),
                            check.names = F)
##
Tcell@meta.data$celltype1 <-"NA"
for (i in 1:length(rownames(Tcell@meta.data))) {
  Tcell@meta.data$celltype1[i] <- paste(Tcell@meta.data$sample[i],Tcell@meta.data$celltype[i])
}
##
Idents(epi) <- "celltype1"
CellsClusters <- data.frame(cell = names(Idents(Tcell)), 
                            cell_type = as.character(Idents(Tcell)),
                            check.names = F)
## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(180, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 
ggsave(filename="/Users/lin/Desktop/胃炎Fig/viper.map.tiff",plot=viper_hmap,width =10,height=10,dpi=300)
ggsave(filename="/Users/lin/Desktop/胃炎Fig/viper.map1.tiff",plot=viper_hmap,width =10,height=10,dpi=300)

