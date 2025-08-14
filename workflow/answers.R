## Map GeneID to Symbol for DE results
r_df <- as.data.frame(r)
r_df$GeneID <- rownames(r_df)
r_df$Symbol <- annot[r_df$GeneID, "Symbol"]
r_df$Symbol <- ifelse(is.na(r_df$Symbol) | r_df$Symbol == "",
                      r_df$GeneID,  # fallback to GeneID if no symbol
                      r_df$Symbol)

## q1 – Approximate % variance explained by UMAP1
umap_coords <- as.data.frame(ump$layout)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_var <- apply(umap_coords, 2, var)
q1_umap1_var <- round(100 * umap_var[1] / sum(umap_var), 2)

## q2 – Most upregulated gene by absolute log2FC (Symbol)
q2_up_gene <- r_df$Symbol[which.max(abs(r_df$log2FoldChange))]

## q3 – Most statistically significant gene by smallest padj (Symbol)
q3_top_gene_by_padj <- r_df$Symbol[which.min(r_df$padj)]

## q4 – Number of significantly DE genes (FDR < 0.05)
q4_n_de <- sum(!is.na(r_df$padj) & r_df$padj < 0.05)

## q5 – Proportion of genes with p-value < 0.001
q5_prop_p_below <- mean(r_df$pvalue < 0.001, na.rm = TRUE)
q5_prop_str <- paste0(round(100 * q5_prop_p_below, 2), "%")

## q6 – Top enriched Biological Process pathway from go_df
if (exists("go_df") && !is.null(go_df) && "Biological Process" %in% go_df$Ontology) {
  q6_top_bp <- go_df$Description[go_df$Ontology == "Biological Process"][1]
} else {
  q6_top_bp <- NA
}

## Output results
cat(
  "q1 UMAP1 variance explained:", q1_umap1_var, "%\n",
  "q2 most upregulated gene (abs log2FC):", q2_up_gene, "\n",
  "q3 most significant gene (lowest padj):", q3_top_gene_by_padj, "\n",
  "q4 n_DE (FDR<0.05):", q4_n_de, "\n",
  "q5 prop p<0.001:", q5_prop_str, "\n",
  "q6 top enriched BP pathway:", q6_top_bp, "\n"
)
