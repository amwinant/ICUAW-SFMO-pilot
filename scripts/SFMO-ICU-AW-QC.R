# Count data -------------------------------------------------------------------

# Feature counts
rna_counts <- GetAssayData(seurat_obj, slot = "counts")
genes_per_cell <- colSums(rna_counts > 0)

fiber_map <- c(
  "14"   = "1",
  "15"   = "2",
  "17"   = "3",
  "18"   = "4",
  "19"   = "5",
  "22"   = "6",
  "24"   = "7",
  "25"   = "8",
  "14_2" = "1",
  "15_2" = "2",
  "20"   = "3",
  "21"   = "4",
  "23"   = "5",
  "26"   = "6",
  "5_2"  = "7",
  "8_2"  = "8"
)


rna_df <- data.frame(
  cell = names(genes_per_cell),
  genes = as.numeric(genes_per_cell),
  sample = seurat_obj$sample,
  condition = recode(seurat_obj$group, "C" = "Control", "S" = "ICU-AW")
)
rna_df <- rna_df %>%
  mutate(sample = recode(as.character(sample), !!!fiber_map))

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
prot_df <- prot_df %>%
  mutate(subject = recode(as.character(sample), !!!fiber_map))
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
prot <- ggplot(prot_df, aes(x = subject, y = proteins, fill = sample)) +
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
