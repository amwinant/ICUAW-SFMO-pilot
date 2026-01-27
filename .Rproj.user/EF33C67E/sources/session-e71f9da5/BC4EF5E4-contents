# Count data -------------------------------------------------------------------

# Feature counts
rna_counts <- GetAssayData(seurat_obj, slot = "counts")
genes_per_cell <- colSums(rna_counts > 0)

rna_df <- data.frame(
  cell = names(genes_per_cell),
  genes = as.numeric(genes_per_cell),
  sample = seurat_obj$sample,
  condition = recode(seurat_obj$group, "C" = "Control", "S" = "ICU-AW")
)

prot_counts <- colSums(!is.na(normalized_data))

prot_df <- data.frame(
  cell = names(prot_counts),
  proteins = prot_counts
) %>%
  left_join(filtered_metadata, by = c("cell" = "fiber_id")) %>%  
  mutate(
    condition = recode(condition, "c" = "Control", "s" = "ICU-AW")
  )

prot_df$sample <- gsub("^[csCS]", "", prot_df$subject)
rna_long <- rna_df %>%
  dplyr::select(sample, condition, count = genes) %>%
  mutate(omic = "RNA: Detected Genes")

prot_long <- prot_df %>%
  dplyr::select(sample, condition, count = proteins) %>%
  mutate(omic = "Proteomics: Detected Proteins")

combined_long <- bind_rows(rna_long, prot_long)

rna <- ggplot(rna_df, aes(x = sample, y = genes, fill = sample)) +
  geom_violin(alpha = 0.65) +
  geom_point(size = 1.5, alpha = 0.75) +
  facet_grid(~ condition, scales = "free_x") +
  scale_fill_viridis_d(option = "plasma") +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Detected Genes per Fiber (RNA-seq)",
    x = "Donor",
    y = "Detected Genes"
  )
prot <- ggplot(prot_df, aes(x = sample, y = proteins, fill = sample)) +
  geom_violin(alpha = 0.65) +
  geom_point(size = 1.5, alpha = 0.75) +
  facet_grid(~ condition, scales = "free_x") +
  scale_fill_viridis_d(option = "plasma") +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Detected Proteins per Fiber (Proteomics)",
    x = "Donor",
    y = "Detected Proteins"
  )
rna / prot
