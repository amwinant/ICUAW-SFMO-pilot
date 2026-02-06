# Load packages and data -------------------------------------------------------
library(vroom)
library(here)
library(PhosR)
library(msigdbr)
library(readr)
library(gt)
library(edgeR)

data_raw <- vroom(here("objects/sepsis_pilot_diann_report.unique_genes_matrix.tsv"))
metadata <- read.csv(here("objects/metadata_protein.csv")) |>
  mutate(fiber_id = Subject.ID, subject = gsub(pattern = "f.*", replacement = "", Subject.ID)
  ) |>
  dplyr::select(!Subject.ID) |>
  mutate(condition = case_when(stringr::str_starts(subject, "s") ~ "s", stringr::str_starts(subject, "c") ~ "c", TRUE ~ "error")
  )
# Data normalizing and filtering -----------------------------------------------

data <- data_raw |>
  tibble::column_to_rownames("Genes") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Analytical.Sample.ID") |>
  dplyr::mutate(
    Analytical.Sample.ID = gsub("_S2-.*", "", Analytical.Sample.ID),
    Analytical.Sample.ID = gsub(".*_", "", Analytical.Sample.ID)
  ) |>
  dplyr::inner_join(
    metadata |> dplyr::select(Analytical.Sample.ID, fiber_id),
    by = "Analytical.Sample.ID"
  ) |>
  dplyr::select(!Analytical.Sample.ID) |>
  tibble::column_to_rownames("fiber_id") |>
  t() |>
  as.data.frame()

transformed_data <- log2(data)

condition_vec <- metadata$condition[match(colnames(transformed_data), metadata$fiber_id)]

stopifnot(!any(is.na(condition_vec)))  

transformed_data |>
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "intensities") |>
  ggplot(aes(x = sample_id, y = intensities)) +
  geom_boxplot(outlier.size = 0.25) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(face = "bold", colour = "black", size = 6)
  )

filtered_data <- PhosR::selectGrps(
  mat = transformed_data,
  percent = 0.7,
  grps = condition_vec
)

normalized_data <- limma::normalizeBetweenArrays(filtered_data, method = "scale") |> as.data.frame()

normalized_data |>
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "intensities") |>
  ggplot(aes(x = sample_id, y = intensities)) +
  geom_boxplot(outlier.size = 0.25) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(face = "bold", colour = "black", size = 6)
  )

filtering_columns_Na <- function(.data, percentage_accepted_missing) {
  keep_vector <- .data |>
    is.na() |>
    colSums()
  filtering_vector <- keep_vector / nrow(.data)
  filtering_vector <- filtering_vector <= percentage_accepted_missing
  vector_id <- data.frame(
    ID = colnames(.data),
    keep = filtering_vector) |>
    dplyr::filter(keep == TRUE)
  data_filtered <- .data |>
    dplyr::select(vector_id$ID)
  return(data_filtered)
}

filtered_data <- transformed_data |>
  filtering_columns_Na(percentage_accepted_missing = 0.5)
filtered_metadata <- metadata %>%
  filter(fiber_id %in% colnames(filtered_data)) %>%
  arrange(match(fiber_id, colnames(filtered_data)))

stopifnot(identical(filtered_metadata$fiber_id, colnames(filtered_data)))

filtered_data <- filtered_data |>
  PhosR::selectGrps(grps = filtered_metadata$condition, percent = 0.7)

normalized_data <- filtered_data |>
  limma::normalizeBetweenArrays(method = "scale") |>
  as.data.frame()

# Dynamic range
log10_protein_data <- log10(normalized_data + 1) 

protein_means <- rowMeans(log10_protein_data, na.rm = TRUE)

protein_df <- data.frame(
  protein = names(protein_means),
  mean_log10_intensity = as.numeric(protein_means)  
)
krt_proteins <- grep("^KRT", protein_df$protein, value = TRUE, ignore.case = TRUE)
protein_df <- protein_df[!protein_df$protein %in% krt_proteins, ]

protein_df <- protein_df[order(protein_df$mean_log10_intensity, decreasing = TRUE), ]
protein_df$rank <- seq_len(nrow(protein_df))

top_proteins <- protein_df[1:10, ]
dyn_prot <- ggplot(protein_df, aes(x = rank, y = mean_log10_intensity)) +
  geom_line(color = "steelblue") +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_text_repel(
    data = top_proteins,
    aes(label = protein),
    size = 3,
    max.overlaps = 20
  ) +
  theme_classic() +
  labs(
    x = "Protein rank",
    y = "Mean log10 intensity",
    title = "Dynamic range of protein abundance"
  )

dyn_rna + dyn_prot

# Seurat object ----------------------------------------------------------------
seurat_object <- normalized_data |>
  tImpute(m = 1.8, s = 0.3) |>
  as.data.frame()

seurat_object <- Seurat::CreateSeuratObject(counts = seurat_object,
                                            project = "sepsis",
                                            meta.data = filtered_metadata)
# Seurat processing of filtered data --------------------------------------

seurat_object <- Seurat::FindVariableFeatures(seurat_object,
                                              selection.method = "vst")

seurat_object[["RNA"]]$data <- seurat_object[["RNA"]]$counts

seurat_object <- Seurat::ScaleData(seurat_object)

seurat_object <- Seurat::RunPCA(seurat_object, features = Seurat::VariableFeatures(object = seurat_object))

seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:10)

pca_df <- data.frame(
  PC1 = Embeddings(seurat_object[["pca"]])[,1],
  PC2 = Embeddings(seurat_object[["pca"]])[,2],
  seurat_object@meta.data
)

var_expl <- (seurat_object[["pca"]]@stdev^2) / sum(seurat_object[["pca"]]@stdev^2) * 100

# Fiber typing -----------------------------------------------------------------
# myh_genes <- c("MYH7", "MYH2", "MYH1")
# myh_mat_log <- seurat_object@assays$RNA$counts[myh_genes, , drop = FALSE]
# myh_mat <- 2^myh_mat_log
# myh_df <- as.data.frame(t(myh_mat))
# myh_df$fiber_ID <- rownames(myh_df)
# perc_MYHs <- myh_df %>%
#   mutate(
#     sum_MYH = MYH7 + MYH2 + MYH1,
#     MYH7 = ifelse(sum_MYH > 0, MYH7 / sum_MYH * 100, 0),
#     MYH2 = ifelse(sum_MYH > 0, MYH2 / sum_MYH * 100, 0),
#     MYH1 = ifelse(sum_MYH > 0, MYH1 / sum_MYH * 100, 0)
#   )
# find_bottom_knee <- function(x) {
#   x_sorted <- sort(x, decreasing = TRUE)
#   y_inv <- 100 - x_sorted
#   curvature <- abs(diff(diff(y_inv)))
#   knee_index <- which.max(curvature) + 1
#   100 - y_inv[knee_index]
# }
# 
# bottom_knee <- list(
#   MYH7 = find_bottom_knee(perc_MYHs$MYH7),
#   MYH2 = find_bottom_knee(perc_MYHs$MYH2),
#   MYH1 = find_bottom_knee(perc_MYHs$MYH1)
# )
# perc_MYHs <- perc_MYHs %>%
#   mutate(
#     fiber_type = case_when(
#       MYH7 >= bottom_knee$MYH7 & MYH2 < bottom_knee$MYH2 & MYH1 < bottom_knee$MYH1 ~ "Type I",
#       MYH2 >= bottom_knee$MYH2 & MYH7 < bottom_knee$MYH7 & MYH1 < bottom_knee$MYH1 ~ "Type IIA",
#       MYH1 >= bottom_knee$MYH1 & MYH7 < bottom_knee$MYH7 & MYH2 < bottom_knee$MYH2 ~ "Type IIX",
#       MYH7 >= bottom_knee$MYH7 & MYH2 >= bottom_knee$MYH2 ~ "Hybrid I/IIA",
#       MYH2 >= bottom_knee$MYH2 & MYH1 >= bottom_knee$MYH1 ~ "Hybrid IIA/IIX",
#       TRUE ~ "Hybrid (other)"
#     )
#   )
# fiber_type_long <- perc_MYHs %>%
#   pivot_longer(
#     cols = c(MYH7, MYH2, MYH1),
#     names_to = "MYH",
#     values_to = "percent"
#   )
# # Cluster fiber type composition -----------------------------------------------
# fiber_df <- fiber_type_long %>%
#   mutate(
#     in_cluster = ifelse(fiber_ID %in% cluster_fibers, "Cluster", "Other")
#   ) %>%
#   filter(!is.na(fiber_type))
# fiber_df <- fiber_type_long %>%
#   mutate(
#     in_cluster = ifelse(fiber_ID %in% cluster_fibers, "Cluster", "Other")
#   ) %>%
#   filter(!is.na(fiber_type))
# plot_df <- fiber_df %>%
#   group_by(in_cluster, fiber_type) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(in_cluster) %>%
#   mutate(percent = 100 * n / sum(n))
# 
# ggplot(plot_df, aes(x = in_cluster, y = percent, fill = fiber_type)) +
#   geom_col(color = "black", width = 0.6) +
#   theme_minimal(base_size = 14) +
#   labs(
#     x = NULL,
#     y = "Percentage of fibers",
#     fill = "Fiber type",
#     title = "Fiber type composition of cluster vs non-cluster fibers"
#   ) +
#   scale_y_continuous(expand = c(0, 0))

# DE CvS -----------------------------------------------------------------------

mat_counts <- as.matrix(GetAssayData(seurat_object, slot = "counts")) 

lin_mat <- 2^as.matrix(GetAssayData(seurat_object, slot = "counts")) - 1 

donor_vec <- seurat_object$subject[colnames(lin_mat)] 
donor_levels <- unique(donor_vec) 
pseudobulk_donor <- sapply(donor_levels, function(d) { rowSums(lin_mat[, donor_vec == d, drop = FALSE], na.rm = TRUE) }) 
rownames(pseudobulk_donor) <- rownames(lin_mat) 
colnames(pseudobulk_donor) <- donor_levels 

metadata$donor_id <- sub("f.*", "", metadata$fiber_id) 

donor_meta <- metadata |> dplyr::select(donor_id, condition) |> distinct() 
donor_meta$condition <- factor(donor_meta$condition, levels = c("c", "s")) 
donor_meta <- donor_meta[match(colnames(pseudobulk_donor), donor_meta$donor_id), ] 
stopifnot(all(colnames(pseudobulk_donor) == donor_meta$donor_id)) 

dge <- DGEList(counts = pseudobulk_donor) 
dge <- calcNormFactors(dge) 

design <- model.matrix(~ 0 + donor_meta$condition) 
colnames(design) <- levels(factor(donor_meta$condition)) 

v <- voom(dge, design) 
fit <- lmFit(v, design) 

contr <- makeContrasts(SvsC = s - c, levels = design) 
fit2 <- contrasts.fit(fit, contr) 
fit2 <- eBayes(fit2) 

DE_results <- topTable(fit2, number = Inf)

DE_results <- DE_results %>% 
  mutate(
    negLogP = -log10(P.Value),
    Signif_FDR = adj.P.Val < 0.05 & abs(logFC) > 1
  )
p_fdr_cutoff <- max(DE_results$P.Value[DE_results$adj.P.Val < 0.05], na.rm = TRUE)

ggplot(DE_results, aes(x = logFC, y = negLogP, color = Signif_FDR)) + 
  geom_point(size = 1.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + 
  geom_hline(yintercept = -log10(p_fdr_cutoff), linetype = "dashed") + 
  scale_color_manual(values = c("grey", "red")) + 
  labs(
    title = "Volcano — ICU-AW vs Control",
    y = "-log10(p-value)"
  ) + 
  theme_minimal()

# DE CvS filtering donors ------------------------------------------------------

fiber_counts <- table(donor_vec) 
donors_keep <- names(fiber_counts)[fiber_counts > 3]
pseudobulk_donor <- pseudobulk_donor[, donors_keep, drop = FALSE]
donor_meta <- donor_meta[donor_meta$donor_id %in% donors_keep, ]
donor_meta <- donor_meta[match(colnames(pseudobulk_donor), donor_meta$donor_id), ]
stopifnot(all(colnames(pseudobulk_donor) == donor_meta$donor_id))

dge <- DGEList(counts = pseudobulk_donor)
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + donor_meta$condition)
colnames(design) <- levels(factor(donor_meta$condition))

v <- voom(dge, design)
fit <- lmFit(v, design)

contr <- makeContrasts(SvsC = s - c, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

DE_results <- topTable(fit2, number = Inf)

DE_results <- DE_results %>% 
  mutate(
    negLogP = -log10(P.Value),
    Signif_FDR = adj.P.Val < 0.05 & abs(logFC) > 1
  )

p_fdr_cutoff <- max(DE_results$P.Value[DE_results$adj.P.Val < 0.05], na.rm = TRUE)

ggplot(DE_results, aes(x = logFC, y = negLogP, color = Signif_FDR)) + 
  geom_point(size = 1.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + 
  geom_hline(yintercept = -log10(p_fdr_cutoff), linetype = "dashed") + 
  scale_color_manual(values = c("grey", "red")) + 
  labs(
    title = "Volcano — ICU-AW vs Control filtered",
    y = "-log10(p-value)"
  ) + 
  theme_minimal()

fiber_counts <- table(donor_vec) %>% as.data.frame()
colnames(fiber_counts) <- c("donor_id", "num_fibers")
fiber_counts$donor_id <- as.character(fiber_counts$donor_id)


fiber_counts <- data.frame(
  donor_raw = donor_vec,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    condition = ifelse(grepl("^c", donor_raw), "Control", "ICU-AW"),
    donor = donor_raw %>% gsub("^[csCS]", "", .),
    donor = recode(donor, !!!fiber_map),
    donor_id = paste(condition, donor, sep = "_")
  ) %>%
  dplyr::group_by(donor_id, condition) %>%
  dplyr::summarise(
    num_fibers = n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    status = ifelse(num_fibers > 3, "kept", "filtered")
  )

ggplot(fiber_counts, aes(
  x = donor_id,
  y = num_fibers,
  fill = status
)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("kept" = "steelblue", "filtered" = "tomato")) +
  scale_x_discrete(labels = function(x) gsub(".*_", "", x)) +
  facet_wrap(~condition, scales = "free_x") +
  labs(
    x = "Donor",
    y = "Number of fibers",
    fill = "Status",
    title = "Number of fibers per donor"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# DE cluster v other -----------------------------------------------------------

cluster_fibers <- c(
  "c22f17", "s14_2f10", "s14_2f8", "s15_2f12", "s15_2f13",
  "s15_2f2", "s15_2f3", "s15_2f7", "s15_2f9", "s20f1",
  "s20f3", "s21f15", "s23f1", "s23f10", "s23f12", "s23f2",
  "s23f3", "s23f4", "s23f5", "s23f6", "s23f7", "s23f8",
  "s23f9", "s26f1", "s26f12", "s26f3", "s5_2f7", "s5_2f8",
  "s8_2f16", "s8_2f2", "s8_2f3", "s8_2f4"
)

filtered_metadata <- filtered_metadata %>%
  mutate(cluster_status = ifelse(fiber_id %in% cluster_fibers, "cluster", "other"))

mat <- as.matrix(normalized_data)
stopifnot(all(colnames(mat) == filtered_metadata$fiber_id))

donor_cluster <- paste(filtered_metadata$subject, filtered_metadata$cluster_status, sep = "_")
split_idx <- split(seq_along(donor_cluster), donor_cluster)

pseudobulk_cluster <- sapply(names(split_idx), function(d) {
  Matrix::rowSums(mat[, split_idx[[d]], drop = FALSE])
})
pseudobulk_cluster[is.na(pseudobulk_cluster)] <- 0
keep <- rowSums(edgeR::cpm(pseudobulk_cluster) > 1) >= 2  
pseudobulk_cluster <- pseudobulk_cluster[keep, ]
colnames(pseudobulk_cluster) <- names(split_idx)

donor_meta_donor <- donor_meta

donor_meta <- filtered_metadata %>%
  dplyr::select(subject, cluster_status) %>%
  distinct() %>%
  dplyr::slice(match(colnames(pseudobulk_cluster), paste(subject, cluster_status, sep = "_")))

stopifnot(identical(paste(donor_meta$subject, donor_meta$cluster_status, sep = "_"),
                    colnames(pseudobulk_cluster)))

design <- model.matrix(~ factor(cluster_status, levels = c("other", "cluster")),
                       data = donor_meta)

dge <- DGEList(counts = pseudobulk_cluster)
dge <- calcNormFactors(dge, method = "TMM")

v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

coef_name <- colnames(design)[2]
results_pseudobulk <- topTable(fit, coef = coef_name, number = Inf)

volcano_df <- results_pseudobulk %>%
  rownames_to_column("protein") %>%
  mutate(
    Signif_FDR = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "sig", "ns"),
    negLogP = -log10(P.Value)
  )

p_fdr_cutoff <- max(
  volcano_df$P.Value[volcano_df$adj.P.Val < 0.05],
  na.rm = TRUE
)

ggplot(volcano_df, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = Signif_FDR), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("sig" = "red", "ns" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_fdr_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = subset(volcano_df, Signif_FDR == "sig"),
    aes(label = protein),
    size = 3,
    max.overlaps = 30
  ) +
  labs(
    title = "Volcano — Cluster vs Other (pseudobulk by donor)",
    x = "log2 Fold Change (cluster / other)",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

# Feature count
feature_count_protein <- colSums(!is.na(filtered_data))

filtered_metadata$feature_count_protein <- feature_count_protein[filtered_metadata$fiber_id]
plot_data_protein <- filtered_metadata %>%
  dplyr::mutate(
    cluster_status = ifelse(fiber_id %in% cluster_fibers, "cluster", "other"),
    cluster_status = factor(cluster_status, levels = c("other", "cluster"))
  ) %>%
  dplyr::select(fiber_id, cluster_status, feature_count_protein)
protein_features <- ggplot(plot_data_protein, aes(x = cluster_status, y = feature_count_protein, fill = cluster_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("other" = "grey70", "cluster" = "steelblue")) +
  theme_classic(base_size = 14) +
  labs(
    x = "Fiber group",
    y = "Number of proteins detected",
    title = "Protein counts in cluster vs other"
  ) +
  theme(legend.position = "none")
rna_features + protein_features

# DEP annotation ---------------------------------------------------------------

sig_df <- results_pseudobulk %>%
  rownames_to_column("protein") %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

up_proteins <- sig_df %>% filter(logFC > 0) %>% pull(protein)
down_proteins <- sig_df %>% filter(logFC < 0) %>% pull(protein)

ego_up <- enrichGO(
  gene = up_proteins,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "fdr",
  readable = TRUE
)

ego_up_sim <- pairwise_termsim(ego_up)

ego_down <- enrichGO(
  gene = down_proteins,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "fdr",
  readable = TRUE
)

ego_down_sim <- pairwise_termsim(ego_down)
tree_up <- treeplot(ego_up_sim, showCategory = 30) +
  ggtitle("GO Treeplot — Upregulated Proteins")

ego_down_df <- as.data.frame(ego_down)

down <- ggplot(ego_down_df, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "C", direction = -1) +
  labs(
    title = "GO Terms — Downregulated Proteins",
    x = NULL,
    y = "Gene Count",
    fill = "FDR"
  ) +
  theme_minimal()

tree_up / down

tree_up
down