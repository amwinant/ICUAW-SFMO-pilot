
# Integration proteins to genes --------------------------------------------

#Full dataset integration

# common_cells <- intersect(colnames(prot_matrix), rownames(seurat_obj@meta.data))
# 
# cluster_cells <- intersect(
#   rownames(seurat_obj@meta.data)[seurat_obj$in_cluster == "cluster"],
#   common_cells
# )
# other_cells <- intersect(
#   rownames(seurat_obj@meta.data)[seurat_obj$in_cluster == "other"],
#   common_cells
# )
# 
# rna_logFC <- rowMeans(rna_matrix[, cluster_cells, drop = FALSE], na.rm = TRUE) -
#   rowMeans(rna_matrix[, other_cells, drop = FALSE], na.rm = TRUE)
# 
# prot_logFC <- rowMeans(prot_matrix[, cluster_cells, drop = FALSE], na.rm = TRUE) -
#   rowMeans(prot_matrix[, other_cells, drop = FALSE], na.rm = TRUE)
# 
# common_genes <- intersect(rownames(rna_matrix), rownames(prot_matrix))
# 
# all_df <- data.frame(
#   Gene = common_genes,
#   logFC_RNA = rna_logFC[common_genes],
#   logFC_prot = prot_logFC[common_genes]
# )
# 
# rna_sig <- read.csv("objects/significant_genes_cluster.csv")
# prot_sig <- read.csv("objects/significant_proteins_cluster.csv")
# 
# rna_genes_sig <- toupper(rna_sig$gene)  
# prot_genes_sig <- toupper(prot_sig$X)   
# 
# all_df <- all_df %>%
#   mutate(
#     sig = toupper(Gene) %in% rna_genes_sig & toupper(Gene) %in% prot_genes_sig
#   )
# 
# ggplot(all_df, aes(x = logFC_RNA, y = logFC_prot)) +
#   geom_point(aes(color = sig), alpha = 0.6) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
#   scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
#   geom_text_repel(
#     data = subset(all_df, sig),
#     aes(label = Gene),
#     max.overlaps = 20,
#     size = 3
#   ) +
#   theme_minimal(base_size = 14) +
#   labs(
#     x = "log2FC RNA (cluster vs other)",
#     y = "log2FC Protein (cluster vs other)",
#     color = "Significant",
#     title = "RNA vs Protein log2FC in cluster vs other fibers"
#   )

#Significant integration
rna <- read.csv("objects/significant_genes_cluster.csv")
prot <- read.csv("objects/significant_proteins_cluster.csv")
colnames(rna)[1] <- "Gene"
colnames(prot)[1] <- "Protein"

overlap_genes <- intersect(rna$Gene, prot$Protein)
length(overlap_genes) 
overlap_genes

overlap_df <- merge(
  rna %>% filter(Gene %in% overlap_genes) %>% dplyr::select(Gene, logFC_RNA = logFC),
  prot %>% filter(Protein %in% overlap_genes) %>% dplyr::select(Protein, logFC_prot = logFC),
  by.x = "Gene", by.y = "Protein"
)

overlap_df$same_direction <- sign(overlap_df$logFC_RNA) == sign(overlap_df$logFC_prot)

table(overlap_df$same_direction)
xlim <- range(overlap_df$logFC_RNA) * 1.1
ylim <- range(overlap_df$logFC_prot) * 1.1

quadrants <- data.frame(
  xmin = c(0, -Inf, -Inf, 0),
  xmax = c(Inf, 0, 0, Inf),
  ymin = c(0, 0, -Inf, -Inf),
  ymax = c(Inf, Inf, 0, 0),
  fill = c("lightgreen", "orange", "lightgreen", "red")
)

ggplot(overlap_df, aes(x = logFC_RNA, y = logFC_prot, label = Gene)) +
  
  geom_rect(data = quadrants, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE, alpha = 0.2) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_text_repel(max.overlaps = 20) +
  scale_fill_manual(values = c("lightgreen" = "lightgreen", "orange" = "orange", "red" = "red")) +
  theme_minimal() +
  labs(x = "log2FC RNA", y = "log2FC Protein", title = "RNA vs Protein Fold Change for Overlapping Genes") +
  theme(legend.position = "none")


# Scoring and integration with functional data ---------------------------------
cluster_fibers <- c(
  "c22f17", "s14_2f10", "s14_2f8", "s15_2f12", "s15_2f13",
  "s15_2f2", "s15_2f3", "s15_2f7", "s15_2f9", "s20f1",
  "s20f3", "s21f15", "s23f1", "s23f10", "s23f12", "s23f2",
  "s23f3", "s23f4", "s23f5", "s23f6", "s23f7", "s23f8",
  "s23f9", "s26f1", "s26f12", "s26f3", "s5_2f7", "s5_2f8",
  "s8_2f16", "s8_2f2", "s8_2f3", "s8_2f4"
)
rna_norm <- readRDS("../pilot-rna/objects/seurat_normalized.rds")
DefaultAssay(rna_norm) <- "RNA"
rna_matrix <- GetAssayData(rna_norm, slot = "data")
prot_matrix <- normalized_data

rna_overlap <- rna_matrix[rownames(rna_matrix) %in% overlap_genes, , drop = FALSE]
prot_overlap <- prot_matrix[rownames(prot_matrix) %in% overlap_genes, , drop = FALSE]

common_cells <- intersect(colnames(rna_overlap), colnames(prot_overlap))
rna_overlap_sub  <- rna_overlap[, common_cells, drop = FALSE]
prot_overlap_sub <- prot_overlap[, common_cells, drop = FALSE]

rna_score  <- colMeans(rna_overlap_sub,  na.rm = TRUE)
prot_score <- colMeans(prot_overlap_sub, na.rm = TRUE)

rna_score_z  <- as.numeric(scale(rna_score))
prot_score_z <- as.numeric(scale(prot_score))

combined_score <- (rna_score_z + prot_score_z) / 2
names(combined_score) <- common_cells

rna_norm@meta.data <- rna_norm@meta.data[1:length(colnames(rna_overlap)), ]
rownames(rna_norm@meta.data) <- colnames(rna_overlap)
rna_norm@meta.data$in_cluster <- ifelse(
  rownames(rna_norm@meta.data) %in% cluster_fibers,
  "cluster",
  "other"
)

meta_data_subset <- rna_norm@meta.data[common_cells, ]


meta_data_subset$Combined_Score <- combined_score
meta_data_subset$in_cluster <- factor(
  ifelse(rownames(meta_data_subset) %in% cluster_fibers, "cluster", "other"),
  levels = c("other", "cluster")
)

func_metrics <- c("srx","drx","t1","t2")
meta_data_subset[, func_metrics] <- lapply(meta_data_subset[, func_metrics], as.numeric)
meta_data_subset <- meta_data_subset[complete.cases(meta_data_subset[, func_metrics]), ]
meta_data_subset$ATP_turnover <- 220 * ( (meta_data_subset$drx * 60) / (100 * meta_data_subset$t1) + (meta_data_subset$srx * 60) / (100 * meta_data_subset$t2) )
p_srx <- wilcox.test(meta_data_subset$srx ~ meta_data_subset$in_cluster)$p.value
p_drx <- wilcox.test(meta_data_subset$drx ~ meta_data_subset$in_cluster)$p.value
p_t1  <- wilcox.test(meta_data_subset$t1  ~ meta_data_subset$in_cluster)$p.value
p_t2  <- wilcox.test(meta_data_subset$t2  ~ meta_data_subset$in_cluster)$p.value
atp  <- wilcox.test(meta_data_subset$ATP_turnover  ~ meta_data_subset$in_cluster)$p.value

plot_srx <- ggplot(meta_data_subset, aes(x = srx, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_srx, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "SRX",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "SRX vs Combined Score"
  )

plot_drx <- ggplot(meta_data_subset, aes(x = drx, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_drx, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "DRX",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "DRX vs Combined Score"
  )

plot_t1 <- ggplot(meta_data_subset, aes(x = t1, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_t1, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "T1",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "T1 vs Combined Score"
  )

plot_t2 <- ggplot(meta_data_subset, aes(x = t2, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_t2, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "T2",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "T2 vs Combined Score"
  )

combined_plot <- (plot_drx + plot_srx) / (plot_t1 + plot_t2)
combined_plot

atp_plot <- ggplot(meta_data_subset, aes(x = ATP_turnover, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    aes(group = in_cluster, color = in_cluster),
    method = "lm",
    se = FALSE
  ) +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(atp, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "Theoretical myosin ATP consumption (a.u.)",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "Theoretical ATP turnover time"
  )
combined_plot + atp_plot
