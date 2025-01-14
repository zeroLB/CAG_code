##Figure5A
p1 <- DimPlot(fibro,split.by = "sample")
ggsave(filename="/Users/lin/Desktop/胃炎/fibro_UMAP_cell_number_bysample.tiff",plot = p1,width = 8,height = 3,dpi = 300)
##Figure5B
p <- FeaturePlot(fibro,c("CCL11","APOE","CXCL14","SOX6","RGS5","ACTG2","MYH11","HHIP"),cols = c("grey", "#CC0033"),ncol = 4)
ggsave(filename="/Users/lin/Desktop/胃炎/fibro-markers-gene-umap.tiff",plot = p,width = 14,height = 6,dpi = 300)
##Figure5C
DefaultAssay(fibro) <- "SCT"
#granulocyte_chemotaxis_score
granulocyte_chemotaxis <- list(c("CCL11","TNFAIP6","CCL8","CXCL1","CCL2",
                                 "DPP4","SYK","IL34","MDK","SLIT2"))
#cytokine activity
cytokine_activity <- list(c("CCL11","CCL8","CXCL1","TNFSF13B","CCL2",
                            "IL33","KITLG","IL7","IL34"))
##BMP_signaling
BMP_signaling <- list(c("PCSK6","BMP5","WNT5A","BMPR1B","TMEM100",
                        "HIPK2","BMP2","SFRP1"))
Inscore <- AddModuleScore(fibro,
                          features = cytokine_activity,nbin = 12,
                          ctrl = 100,
                          name = "cytokine_activity")
Inscore <- AddModuleScore(Inscore,
                          features = granulocyte_chemotaxis,nbin = 12,
                          ctrl = 100,
                          name = "granulocyte_chemotaxis")
Inscore <- AddModuleScore(Inscore,
                          features = BMP_signaling,nbin = 12,
                          ctrl = 100,
                          name = "BMP_signaling")
colnames(Inscore@meta.data)[23] <- 'granulocyte_chemotaxis_score' 
colnames(Inscore@meta.data)[22] <- 'cytokine_activity_score' 
colnames(Inscore@meta.data)[24] <- 'BMP_signaling_score' 
p1 <- VlnPlot(Inscore,features = 'granulocyte_chemotaxis_score', 
              pt.size = 0, adjust = 2,group.by = "celltype")
p2 <- VlnPlot(Inscore,features = 'cytokine_activity_score', 
              pt.size = 0, adjust = 2,group.by = "celltype")
p3 <- VlnPlot(Inscore,features = 'BMP_signaling_score', 
              pt.size = 0, adjust = 2,group.by = "celltype")
ggsave(filename="/Users/lin/Desktop/胃炎/granulocyte_chemotaxis_score.tiff",plot = p1,width = 10,height = 6,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/cytokine_activity_score.tiff",plot = p2,width = 10,height = 6,dpi = 300)
ggsave(filename="/Users/lin/Desktop/胃炎/BMP_signaling_score.tiff",plot = p3,width = 10,height = 6,dpi = 300)
##Figure5D
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
IAF <- subset(fibro,idents = "CCL11+APOE+ fibroblasts")
Idents(IAF) <- "sample"
deg <-FindAllMarkers(IAF,only.pos = T)
#1. Define a set of potential ligands for both the sender-agnostic and sender-focused approach
receiver = "C1Q+ Macrophage"
expressed_genes_receiver <- get_expressed_genes(receiver, myeloid.C1Q, pct = 0.05)
lr_network <- lr_network %>% distinct(from, to)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
sender_celltypes <- c("CCL11+APOE+ fibroblasts")
# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes,IAF, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
#2. Define the gene set of interest
condition_oi <-  "chronic atrophic gastritis"
condition_reference <- "control"
#seurat_obj_receiver <- subset(seuratObj, idents = receiver)
DE_table_receiver <-  FindMarkers(object = myeloid.C1Q,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "sample",
                                  min.pct = 0.05) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
best_upstream_ligands <- c("IL7","TNFSF13B","CCL2","CCL8","IL34","IL33","KITLG","CSF1")
best_upstream_ligands <-unique(chemokine)
best_upstream_ligands <- c(cytokine)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

make_heatmap_ggplot(vis_ligand_aupr,
                    "Prioritized ligands", "Ligand activity", 
                    legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()
nrow(active_ligand_target_links_df)
## [1] 637
head(active_ligand_target_links_df)
## # A tibble: 6 × 3
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 
nrow(active_ligand_target_links)
## [1] 86
head(active_ligand_target_links)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
p_ligand_target<- make_heatmap_ggplot(vis_ligand_target, "Ligands of interest", "Predicted target genes",
                                      color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high ="purple",breaks = c(0,0.005,0.009))
p_ligand_target
ggsave(filename="/Users/lin/Desktop/胃炎/fibro_macrophage_ligand_target.tiff",plot = p_ligand_target,
       width = 12,height = 3,dpi = 300)
##Figure5E
receiver = "CD8+ effector T cell"
expressed_genes_receiver <- get_expressed_genes(receiver, Tcell, pct = 0.05)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
sender_celltypes <- c("CCL11+APOE+ fibroblasts")
# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, cell, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
# Also check 
length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)
condition_oi <-  "chronic atrophic gastritis"
condition_reference <- "control"
seurat_obj_receiver <- subset(Tcell, idents = receiver)
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "sample",
                                  min.pct = 0.05) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
make_heatmap_ggplot(vis_ligand_aupr,
                    "Prioritized ligands", "Ligand activity", 
                    legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()
nrow(active_ligand_target_links_df)
## [1] 637
head(active_ligand_target_links_df)
## # A tibble: 6 × 3
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 
nrow(active_ligand_target_links)
## [1] 86
head(active_ligand_target_links)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
p_ligand_target<- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                      color = "royalblue", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "royalblue",breaks = c(0,0.0025,0.005))
##Figure5F
p <- DotPlot(IAF,c("IL33","CCL2","CSF1","IL15"),cols = c("grey", "#CC0033"),split.by = "sample")
df <- p$data
colnames(df)[4] <- "dataset"
df$dataset <-c(rep("control",4),rep("CAG",4))

p <-ggplot(df,aes(x=features.plot,y=pct.exp,fill=dataset))+
  geom_bar(stat = "identity",position = position_dodge(0.55),color="black",width = 0.5,linewidth=0.3)+
  theme_bw()+
  xlab("gene")+ylab("cell propotion in CCL11+APOE+ fibroblasts")+scale_fill_brewer(palette = 'Set2')+theme(axis.line = element_blank())
ggsave(filename="/Users/lin/Desktop/胃炎/CCL11+APOE+ fibroblastsratio.tiff",plot = p,width = 8,height = 6,dpi = 300)




