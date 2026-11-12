setwd("~/HBV_HCC")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table)
  library(SingleR); library(celldex); library(BiocParallel)
  library(svglite); library(ggrastr)
})

ann_dir <- "results/annotation"; fig_dir <- "results/figures"
dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

hbv <- readRDS("results/HBV_HCC_integrated.rds")
Idents(hbv) <- "seurat_clusters"
n_cores <- max(1, parallel::detectCores() - 2)

save_fig <- function(p, name, w = 7, h = 6, raster_points = FALSE) {
  if (raster_points) p <- ggrastr::rasterise(p, layers = "Point", dpi = 300, dev = "ragg")
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p, width = w, height = h,
         device = svglite::svglite)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  invisible(p)
}
umap_plot <- function(obj, group, label = FALSE, title = group, cols = NULL) {
  DimPlot(obj, reduction = "umap", group.by = group, label = label, repel = TRUE,
          raster = FALSE, pt.size = 0.1, cols = cols, combine = FALSE)[[1]] +
    ggtitle(title) + coord_fixed()
}

markers <- list(
  Hepatocyte      = c("ALB", "APOA1", "APOC3", "TTR", "SERPINA1", "HP", "FGB"),
  Cholangiocyte   = c("KRT19", "KRT7", "EPCAM", "SOX9", "SPP1"),
  Endothelial     = c("PECAM1", "VWF", "PLVAP", "CDH5", "ENG", "CLEC4G", "FCN3"),
  Fibroblast_HSC  = c("COL1A1", "COL1A2", "COL3A1", "DCN", "ACTA2", "RGS5", "TAGLN"),
  T_cell          = c("CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD8A"),
  NK_cell         = c("NKG7", "GNLY", "KLRD1", "KLRF1", "FCGR3A"),
  B_cell          = c("MS4A1", "CD79A", "CD19", "BANK1"),
  Plasma_cell     = c("JCHAIN", "MZB1", "IGHG1", "SDC1"),
  Myeloid         = c("CD68", "CD14", "LYZ", "CSF1R", "C1QA", "C1QB"),
  Dendritic_cell  = c("CD1C", "FCER1A", "CLEC10A", "CLEC9A", "LILRA4"),
  Neutrophil      = c("FCGR3B", "CSF3R", "CXCR2"),
  Mast_cell       = c("TPSAB1", "CPA3", "KIT"),
  Proliferating   = c("MKI67", "TOP2A", "CENPF"),
  Erythrocyte     = c("HBB", "HBA1", "HBA2")
)
markers <- lapply(markers, intersect, rownames(hbv))
print(lengths(markers))


ref  <- celldex::HumanPrimaryCellAtlasData()          
expr <- GetAssayData(hbv, assay = "RNA", layer = "data")


pred_cl <- SingleR(test = expr, ref = ref, labels = ref$label.main,
                   clusters = hbv$seurat_clusters,
                   BPPARAM = MulticoreParam(n_cores))

pred_cell <- SingleR(test = expr, ref = ref, labels = ref$label.main,
                     BPPARAM = MulticoreParam(n_cores))
hbv$SingleR_cell    <- pred_cell$pruned.labels
hbv$SingleR_cluster <- pred_cl$labels[match(hbv$seurat_clusters, rownames(pred_cl))]
rm(expr); gc()

svglite(file.path(fig_dir, "SingleR_cluster_score_heatmap.svg"), width = 9, height = 7)
plotScoreHeatmap(pred_cl, show.labels = TRUE)
dev.off()


hbv <- AddModuleScore(hbv, features = markers, name = "MS_", seed = 1)
ms_cols <- paste0("MS_", seq_along(markers))
colnames(hbv@meta.data)[match(ms_cols, colnames(hbv@meta.data))] <- paste0("MS_", names(markers))
ms_cols <- paste0("MS_", names(markers))

ms <- as.data.table(hbv@meta.data[, c("seurat_clusters", ms_cols)])
ms_mean <- ms[, lapply(.SD, mean), by = seurat_clusters, .SDcols = ms_cols]
setorder(ms_mean, seurat_clusters)
sc_mat  <- as.matrix(ms_mean[, -1]); rownames(sc_mat) <- ms_mean$seurat_clusters
colnames(sc_mat) <- names(markers)

top1 <- apply(sc_mat, 1, function(x) names(sort(x, decreasing = TRUE))[1])
top2 <- apply(sc_mat, 1, function(x) names(sort(x, decreasing = TRUE))[2])
gap  <- apply(sc_mat, 1, function(x) diff(sort(x, decreasing = TRUE)[2:1]))


hm <- melt(as.data.table(scale(sc_mat), keep.rownames = "cluster"),
           id.vars = "cluster", variable.name = "panel", value.name = "z")
hm$cluster <- factor(hm$cluster, levels = rownames(sc_mat))
p <- ggplot(hm, aes(panel, cluster, fill = z)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  theme_minimal(base_size = 11) + labs(x = NULL, y = "Cluster", fill = "Scaled\nscore") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_fig(p, "MarkerScore_heatmap_by_cluster", w = 8, h = 8)


final_file <- file.path(ann_dir, "cluster_annotation_final.csv")
lab <- if (file.exists(final_file)) fread(final_file) else agree
cat("Using labels from:", if (file.exists(final_file)) final_file else "automatic markers", "\n")
hbv$celltype <- lab$final[match(as.character(hbv$seurat_clusters),
                                as.character(lab$seurat_clusters))]

ct_levels <- intersect(c(names(markers), sort(unique(hbv$celltype))), unique(hbv$celltype))
hbv$celltype <- factor(hbv$celltype, levels = ct_levels)
ct_cols <- setNames(scales::hue_pal()(length(ct_levels)), ct_levels)


save_fig(umap_plot(hbv, "seurat_clusters", label = TRUE, title = "Clusters"),
         "UMAP_clusters", w = 8, h = 7, raster_points = TRUE)
save_fig(umap_plot(hbv, "celltype", label = TRUE, title = "Cell type", cols = ct_cols),
         "UMAP_celltype", w = 8.5, h = 7, raster_points = TRUE)
save_fig(umap_plot(hbv, "celltype", title = "Cell type", cols = ct_cols) + NoLegend(),
         "UMAP_celltype_nolegend", w = 7, h = 7, raster_points = TRUE)
save_fig(umap_plot(hbv, "SingleR_cluster", label = TRUE, title = "SingleR (HPCA)"),
         "UMAP_SingleR", w = 8.5, h = 7, raster_points = TRUE)
for (g in c("dataset", "patient", "tissue")) {
  save_fig(umap_plot(hbv, g), paste0("UMAP_", g), w = 8, h = 7, raster_points = TRUE)
}

p <- DimPlot(hbv, reduction = "umap", group.by = "celltype", split.by = "tissue",
             raster = FALSE, pt.size = 0.1, cols = ct_cols, combine = FALSE)[[1]]
save_fig(p, "UMAP_celltype_by_tissue", w = 13, h = 6.5, raster_points = TRUE)


feat_dir <- file.path(fig_dir, "FeaturePlots"); dir.create(feat_dir, showWarnings = FALSE)
key_genes <- intersect(c("ALB", "APOA1", "KRT19", "EPCAM", "PECAM1", "PLVAP", "COL1A1",
                         "ACTA2", "CD3D", "CD8A", "IL7R", "FOXP3", "NKG7", "MS4A1",
                         "JCHAIN", "CD68", "C1QA", "CD1C", "TPSAB1", "MKI67"), rownames(hbv))
for (gn in key_genes) {
  p <- FeaturePlot(hbv, gn, reduction = "umap", raster = FALSE, pt.size = 0.1,
                   order = TRUE, combine = FALSE)[[1]] + coord_fixed()
  p <- ggrastr::rasterise(p, layers = "Point", dpi = 300, dev = "ragg")
  ggsave(file.path(feat_dir, paste0(gn, ".svg")), p, width = 5, height = 4.5, device = svglite)
}

dp_genes <- unique(unlist(lapply(markers, head, 3)))
p <- DotPlot(hbv, features = dp_genes, group.by = "celltype") + RotatedAxis() +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 8))
save_fig(p, "DotPlot_markers_by_celltype", w = 14, h = 6)
p <- DotPlot(hbv, features = dp_genes, group.by = "seurat_clusters") + RotatedAxis() +
  labs(x = NULL, y = "Cluster") + theme(axis.text.x = element_text(size = 8))
save_fig(p, "DotPlot_markers_by_cluster", w = 14, h = 9)

for (q in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
  p <- VlnPlot(hbv, q, group.by = "sample", pt.size = 0, combine = FALSE)[[1]] +
    NoLegend() + labs(x = NULL)
  save_fig(p, paste0("QC_after_", q), w = 10, h = 4)
}

comp <- as.data.table(hbv@meta.data)[, .N, by = .(sample, patient, tissue, celltype)]
comp[, frac := N / sum(N), by = sample]
p <- ggplot(comp, aes(sample, frac, fill = celltype)) +
  geom_col(width = 0.85) + scale_fill_manual(values = ct_cols) +
  facet_grid(~ tissue, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of cells", fill = "Cell type") +
  theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "Composition_by_sample", w = 11, h = 5)

p <- ggplot(comp[, .(N = sum(N)), by = .(tissue, celltype)][, frac := N / sum(N), by = tissue],
            aes(tissue, frac, fill = celltype)) +
  geom_col(width = 0.7) + scale_fill_manual(values = ct_cols) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of cells", fill = "Cell type") + theme_classic(base_size = 11)
save_fig(p, "Composition_by_tissue", w = 4.5, h = 5)

fwrite(dcast(comp, sample + patient + tissue ~ celltype, value.var = "N", fill = 0),
       file.path(ann_dir, "celltype_counts_by_sample.csv"))


saveRDS(hbv, "results/HBV_HCC_annotated.rds")
cat("\nDone. Figures in", fig_dir, "| tables in", ann_dir, "\n")


