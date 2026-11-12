setwd("~/HBV_HCC")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(DESeq2)
  library(ggplot2); library(ggrepel); library(svglite)
})
out_dir <- "results/DEG"; fig_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS   <- 20      
MIN_PAIRS   <- 3       
PADJ        <- 0.05
LFC         <- 1

hbv    <- readRDS("results/HBV_HCC_annotated_cnv.rds")
counts <- LayerData(hbv, assay = "RNA", layer = "counts")
meta   <- as.data.table(hbv@meta.data, keep.rownames = "cell")


pseudobulk <- function(cells, group) {
  g  <- factor(group)
  mm <- sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1,
                     dims = c(length(g), nlevels(g)), dimnames = list(NULL, levels(g)))
  as.matrix(counts[, cells, drop = FALSE] %*% mm)
}

volcano <- function(res, title, name) {
  d <- as.data.table(res, keep.rownames = "gene")[!is.na(padj)]
  d[, cls := fifelse(padj < PADJ & log2FoldChange >  LFC, "Up",
                     fifelse(padj < PADJ & log2FoldChange < -LFC, "Down", "n.s."))]
  top <- rbind(d[cls == "Up"][order(padj)][1:min(.N, 12)],
               d[cls == "Down"][order(padj)][1:min(.N, 12)])
  p <- ggplot(d, aes(log2FoldChange, -log10(padj), colour = cls)) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_text_repel(data = top, aes(label = gene), size = 2.2, max.overlaps = 30,
                    segment.size = 0.2, colour = "black") +
    geom_vline(xintercept = c(-LFC, LFC), linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = -log10(PADJ), linetype = 2, linewidth = 0.3) +
    scale_colour_manual(values = c(Up = "#E64B35", Down = "#3C5488", n.s. = "grey80")) +
    labs(title = title, x = "log2 fold change", y = "-log10 adj. P", colour = NULL) +
    theme_classic(base_size = 8) + theme(plot.title = element_text(size = 9, face = "bold"))
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p, width = 4, height = 3.6, device = svglite)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = 4, height = 3.6, dpi = 200)
}


meta[, ct_deg := as.character(celltype)]
meta[celltype_cnv == "Malignant", ct_deg := NA]
summary_A <- list()
for (ct in sort(unique(na.omit(meta$ct_deg)))) {
  m <- meta[ct_deg == ct]
  n <- m[, .N, by = .(patient, tissue)][N >= MIN_CELLS]
  pts <- n[, .N, by = patient][N == 2, patient]
  if (length(pts) < MIN_PAIRS) {
    cat(ct, ": only", length(pts), "paired patients -> skipped\n"); next
  }
  m  <- m[patient %in% pts]
  pb <- pseudobulk(m$cell, paste(m$patient, m$tissue, sep = "|"))
  cd <- data.frame(row.names = colnames(pb),
                   patient = factor(sub("\\|.*", "", colnames(pb))),
                   tissue  = factor(sub(".*\\|", "", colnames(pb)), levels = c("Normal", "Tumor")))
  pb <- pb[rowSums(pb >= 10) >= MIN_PAIRS, ]
  dds <- DESeqDataSetFromMatrix(round(pb), cd, design = ~ patient + tissue)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("tissue", "Tumor", "Normal"))
  res <- res[order(res$padj), ]
  fwrite(as.data.table(as.data.frame(res), keep.rownames = "gene"),
         file.path(out_dir, paste0("DEG_TumorVsNormal_", ct, ".csv")))
  volcano(res, paste0(ct, ": Tumor vs Normal (n = ", length(pts), " pairs)"),
          paste0("DEG_volcano_TumorVsNormal_", ct))
  up <- sum(res$padj < PADJ & res$log2FoldChange >  LFC, na.rm = TRUE)
  dn <- sum(res$padj < PADJ & res$log2FoldChange < -LFC, na.rm = TRUE)
  summary_A[[ct]] <- data.table(celltype = ct, n_pairs = length(pts), up = up, down = dn)
  cat(sprintf("%-16s pairs=%d  up=%d  down=%d\n", ct, length(pts), up, dn))
}
sA <- rbindlist(summary_A)
fwrite(sA, file.path(out_dir, "DEG_TumorVsNormal_summary.csv"))
pl <- melt(sA, id.vars = c("celltype", "n_pairs"), variable.name = "dir", value.name = "n")
pl[dir == "down", n := -n]
p <- ggplot(pl, aes(reorder(celltype, abs(n)), n, fill = dir)) + geom_col(width = 0.7) +
  coord_flip() + geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(up = "#E64B35", down = "#3C5488"),
                    labels = c(up = "Up in tumour", down = "Down in tumour"), name = NULL) +
  labs(x = NULL, y = sprintf("Number of DEGs (padj < %.2g, |log2FC| > %g)", PADJ, LFC)) +
  theme_classic(base_size = 8)
ggsave(file.path(fig_dir, "DEG_TumorVsNormal_counts.svg"), p, width = 4.5, height = 3, device = svglite)
ggsave(file.path(fig_dir, "DEG_TumorVsNormal_counts.png"), p, width = 4.5, height = 3, dpi = 200)


mB <- rbind(
  meta[celltype_cnv == "Malignant" & tissue == "Tumor"][, grp := "Malignant"],
  meta[celltype == "Hepatocyte" & tissue == "Normal" & celltype_cnv != "Malignant"][, grp := "NormalHep"])
nB <- mB[, .N, by = .(patient, grp)][N >= MIN_CELLS]
print(nB)
mB <- merge(mB, nB[, .(patient, grp)], by = c("patient", "grp"))
pb <- pseudobulk(mB$cell, paste(mB$patient, mB$grp, mB$dataset, sep = "|"))
cd <- data.frame(row.names = colnames(pb),
                 patient = sub("\\|.*", "", colnames(pb)),
                 group   = factor(sapply(strsplit(colnames(pb), "\\|"), `[`, 2),
                                  levels = c("NormalHep", "Malignant")),
                 dataset = factor(sapply(strsplit(colnames(pb), "\\|"), `[`, 3)))
print(cd)
pb <- pb[rowSums(pb >= 10) >= 3, ]
form <- if (nlevels(droplevels(cd$dataset)) > 1 &&
            all(table(cd$dataset, cd$group) > 0)) ~ dataset + group else ~ group
cat("Design for malignant vs normal hepatocytes:", deparse(form), "\n")
dds <- DESeq(DESeqDataSetFromMatrix(round(pb), cd, design = form), quiet = TRUE)
res <- results(dds, contrast = c("group", "Malignant", "NormalHep"))
res <- res[order(res$padj), ]
fwrite(as.data.table(as.data.frame(res), keep.rownames = "gene"),
       file.path(out_dir, "DEG_Malignant_vs_NormalHepatocyte.csv"))
volcano(res, "Malignant vs normal hepatocytes", "DEG_volcano_Malignant_vs_NormalHep")
cat("Malignant vs normal hepatocytes: up =",
    sum(res$padj < PADJ & res$log2FoldChange > LFC, na.rm = TRUE), " down =",
    sum(res$padj < PADJ & res$log2FoldChange < -LFC, na.rm = TRUE), "\n")
cat("Done. Tables in", out_dir, "\n")

