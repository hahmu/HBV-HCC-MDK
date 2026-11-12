setwd("~/HBV_HCC")
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(data.table); library(ggplot2)
  library(svglite); library(ggrastr); library(ggsci)
})
set.seed(1)
options(future.globals.maxSize = 50 * 1024^3)

COMPARTMENTS <- list(
  TNK     = c("T_cell", "NK_cell"),
  Myeloid = c("Myeloid", "Dendritic_cell"),
  Stromal = c("Endothelial", "Fibroblast_HSC"))

COMP_LINEAGE <- list(TNK = "TNK", Myeloid = "Myeloid",
                     Stromal = c("Endothelial", "Fibroblast"))


RES      <- c(TNK = 0.8, Myeloid = 0.6, Stromal = 1.0)
RES_GRID <- c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0)
RES_GRID <- sort(unique(c(RES_GRID, RES)))      
NPCS     <- c(TNK = 30, Myeloid = 25, Stromal = 30)
N_HVG    <- 2000
N_MARK   <- 5            
ASSIGN_PROLIFERATING <- TRUE

out_root <- "results/subcluster"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)


AMBIENT <- c("ALB", "APOA1", "APOA2", "APOC3", "APOH", "APOB", "TTR", "HP", "HPX", "FGA", "FGB",
             "FGG", "SERPINA1", "TF", "AHSG", "RBP4", "ORM1", "ORM2", "SAA1", "SAA2", "AMBP",
             "VTN", "FABP1", "CYP2E1")
IEG     <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "IER2", "IER3", "ATF3", "DUSP1",
             "ZFP36", "HSPB1")
hvg_exclude <- function(g) {
  grepl(paste0("^MT-|^RP[SL][0-9]|^HSPA|^HSPH|^HSP90|^DNAJ|^IG[HKL][VDJ]|^IGH[GAMDE]|",
               "^IGKC$|^IGLC|^TR[ABDG][VDJ]|^HB[AB]|^MALAT1$|^NEAT1$|^XIST$|^RPS4Y|",
               "^DDX3Y$|^UTY$|^KDM5D$|^EIF1AY$"), g) | g %in% c(AMBIENT, IEG)
}


LINEAGE <- list(
  TNK         = c("CD3D", "CD3E", "CD2", "CD247", "NKG7", "GNLY", "KLRD1"),
  Myeloid     = c("CD68", "CD14", "LYZ", "CSF1R", "C1QA", "AIF1", "FCER1G"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP", "CLEC4G", "EGFL7"),
  Fibroblast  = c("COL1A1", "COL1A2", "COL3A1", "DCN", "RGS5", "ACTA2", "PDGFRB"),
  B_plasma    = c("MS4A1", "CD79A", "JCHAIN", "MZB1"),
  Epithelial  = c("ALB", "APOA1", "TTR", "KRT19", "EPCAM", "KRT8"))

STATE <- list(
  Cycling = c("MKI67", "TOP2A", "CENPF", "NUSAP1", "UBE2C"),
  Stress  = c("HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "FOS", "JUN"),
  ISG     = c("ISG15", "IFI6", "IFIT3", "MX1", "IFI44L"),
  Ambient = c("ALB", "APOA1", "APOC3", "TTR", "HP", "SERPINA1"))

PANELS <- list(
  TNK = list(
    CD4_naive_CM     = c("CCR7", "SELL", "LEF1", "TCF7", "IL7R"),
    CD4_memory       = c("CD4", "CD40LG", "ANXA1", "IL7R", "KLRB1"),
    Tfh_CXCL13       = c("CXCL13", "PDCD1", "ICOS", "TOX2", "CD4"),
    Treg             = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF4", "IKZF2"),
    CD8_GZMK_EM      = c("CD8A", "GZMK", "EOMES", "CCL5", "CD8B"),
    CD8_exhausted    = c("PDCD1", "HAVCR2", "LAG3", "TOX", "ENTPD1", "TIGIT", "CXCL13"),
    CD8_TRM          = c("ITGAE", "ZNF683", "CXCR6", "ITGA1", "CD69"),
    CD8_cytotoxic    = c("GZMB", "GZMH", "PRF1", "CX3CR1", "FGFBP2"),
    MAIT             = c("SLC4A10", "KLRB1", "ZBTB16", "RORC", "TRAV1-2"),
    gdT              = c("TRDC", "TRGC1", "TRGC2", "TRDV2"),
    NK_CD56bright_LR = c("NCAM1", "XCL1", "XCL2", "CD160", "KLRC1"),
    NK_CD56dim       = c("FCGR3A", "FGFBP2", "CX3CR1", "SPON2", "PRF1"),
    ILC              = c("KIT", "IL7R", "TCF7", "TOX2", "IL1R1", "IL4I1")),
  Myeloid = list(
    Kupffer          = c("MARCO", "CD5L", "VSIG4", "TIMD4", "C1QC"),
    FOLR2_TAM        = c("FOLR2", "LYVE1", "MRC1", "SELENOP", "F13A1"),
    SPP1_TAM         = c("SPP1", "TREM2", "GPNMB", "CD9", "FABP5", "APOE"),
    Inflammatory_mac = c("IL1B", "CXCL8", "CCL3", "CCL4", "EREG"),
    Classical_mono   = c("CD14", "FCN1", "S100A8", "S100A9", "VCAN"),
    Nonclassical_mono= c("FCGR3A", "CDKN1C", "LST1", "MS4A7"),
    cDC1             = c("CLEC9A", "XCR1", "CADM1", "IRF8"),
    cDC2             = c("CD1C", "FCER1A", "CLEC10A", "CD1E"),
    mregDC_LAMP3     = c("LAMP3", "CCR7", "FSCN1", "CCL19"),
    pDC              = c("LILRA4", "CLEC4C", "IL3RA", "IRF7"),
    Neutrophil       = c("FCGR3B", "CSF3R", "CXCR2")),
  
  Fibroblast = list(
    cyHSC            = c("HGF", "LRAT", "RELN", "COLEC11", "COLEC10", "NGFR", "RSPO3"),
    myHSC            = c("COL1A1", "COL1A2", "COL3A1", "LUM", "TIMP1", "POSTN"),
    Perivascular_HSC = c("RGS5", "KCNJ8", "HIGD1B", "NDUFA4L2", "PDGFRB"),
    CAF_matrix       = c("POSTN", "FAP", "COL11A1", "CTHRC1", "MMP11"),
    iCAF             = c("CCL19", "CCL21", "CCL2", "CXCL2", "ICAM1"),
    Portal_fibroblast= c("PI16", "CD34", "MFAP5", "CD248", "SFRP2"),
    vSMC             = c("MYH11", "CNN1", "DES", "ACTG2", "RERGL"),
    Mesothelial      = c("MSLN", "UPK3B", "KRT19")),
  Endothelial = list(
    LSEC             = c("CLEC4G", "CLEC4M", "FCN3", "STAB2", "CRHBP"),
    Capillarized_PLVAP = c("PLVAP", "VWA1", "COL4A1", "COL4A2", "SPARC"),
    Tip_cell         = c("ESM1", "APLN", "ANGPT2", "PGF", "KCNE3"),
    Arterial         = c("GJA5", "HEY1", "SEMA3G", "EFNB2", "SOX17"),
    Venous_ACKR1     = c("ACKR1", "SELP", "CCL14", "AQP1"),
    Lymphatic        = c("PROX1", "LYVE1", "CCL21", "TFF3", "PDPN")))
PANELS$Stromal <- c(PANELS$Endothelial, PANELS$Fibroblast)


MARKERS_SHOW <- list(
  TNK = list(
    "CD4 T"             = c("CD4", "CD40LG", "IL7R", "CCR7", "ANXA1"),
    "Tfh-like CD4 T"    = c("CXCL13", "ICOS", "PDCD1", "TOX2", "CD4"),
    "Treg"              = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF4", "IKZF2"),
    "CD8 Tem"           = c("CD8A", "CD8B", "GZMK", "EOMES", "CCL5"),
    "CD8 Trm"           = c("CD8A", "CXCR6", "ITGAE", "ZNF683", "ITGA1"),
    "CD8 Temra"         = c("CD8A", "GZMH", "GZMB", "FGFBP2", "CX3CR1"),
    "MAIT"              = c("SLC4A10", "KLRB1", "RORC", "ZBTB16", "CCR6"),
    "gdT"               = c("TRDC", "TRGC1", "TRGC2", "TRDV2", "KLRD1"),
    "NK"                = c("FCGR3A", "FGFBP2", "SPON2", "PRF1", "GNLY"),
    "Liver-resident NK" = c("XCL1", "XCL2", "CD160", "NCAM1", "KLRC1"),
    "ILC"               = c("KIT", "IL7R", "TCF7", "TOX2", "IL1R1"),
    "Proliferating"     = c("MKI67", "TOP2A", "CENPF", "NUSAP1", "UBE2C")),
  Myeloid = list(
    "Kupffer"                = c("MARCO", "CD5L", "VSIG4", "TIMD4", "C1QC"),
    "FOLR2+ TAM"             = c("FOLR2", "LYVE1", "MRC1", "SELENOP", "F13A1"),
    "SPP1+ TAM"              = c("SPP1", "TREM2", "GPNMB", "CD9", "FABP5"),
    "Classical monocyte"     = c("CD14", "FCN1", "S100A8", "S100A9", "VCAN"),
    "Non-classical monocyte" = c("FCGR3A", "CDKN1C", "LST1", "MS4A7", "CX3CR1"),
    "cDC1"                   = c("CLEC9A", "XCR1", "CADM1", "IRF8", "IDO1"),
    "cDC2"                   = c("CD1C", "FCER1A", "CLEC10A", "CD1E", "HLA-DQA1"),
    "mregDC"                 = c("LAMP3", "CCR7", "FSCN1", "CCL19", "CD40"),
    "pDC"                    = c("LILRA4", "CLEC4C", "IL3RA", "IRF7", "TCF4"),
    "Proliferating"          = c("MKI67", "TOP2A", "CENPF", "NUSAP1", "UBE2C")),
  
  Stromal = list(
    "LSEC"              = c("CLEC4G", "CLEC4M", "FCN3", "STAB2", "CRHBP"),
    "Venous EC"         = c("ACKR1", "SELP", "SELE", "VCAM1", "CCL14"),
    "Arterial EC"       = c("GJA5", "HEY1", "SEMA3G", "FBLN5", "SOX17"),
    "Tumour EC"         = c("PLVAP", "VWA1", "COL4A1", "ESM1", "SPARC"),
    "Lymphatic EC"      = c("PROX1", "PDPN", "CCL21", "TFF3", "MMRN1"),
    "Proliferating EC"  = c("PECAM1", "MKI67", "TOP2A", "CENPF", "NUSAP1"),
    "cyHSC"             = c("HGF", "RELN", "LRAT", "COLEC10", "RSPO3"),
    "Mesenchymal-like HSC"  = c("RGS5", "KCNJ8", "HIGD1B", "PLXDC1", "COLEC11"),
    "myHSC"             = c("COL1A1", "COL3A1", "POSTN", "FAP", "CTHRC1"),
    "iCAF"              = c("CCL19", "CCL21", "CCL2", "CXCL2", "ICAM1"),
    "vSMC"              = c("MYH11", "CNN1", "DES", "ACTG2", "RERGL"),
    "Portal fibroblast" = c("PI16", "CD34", "MFAP5", "CD248", "SFRP2"),
    "Mesothelial"       = c("MSLN", "UPK3B", "KRT19", "CALB2", "WT1"),
    "Proliferating HSC" = c("PDGFRB", "MKI67", "TOP2A", "CENPF", "NUSAP1")))


npg16 <- c(pal_npg("nrc")(10), "#FFDC91", "#6F99AD", "#EE4C97", "#7876B1", "#20854E", "#E18727")
npg_cols  <- function(n) if (n <= length(npg16)) npg16[seq_len(n)] else colorRampPalette(npg16)(n)
name_cols <- function(x) { lv <- levels(factor(x)); setNames(npg_cols(length(lv)), lv) }
tissue_cols <- c(Normal = "#4DBBD5", Tumor = "#E64B35")

save_fig <- function(p, path, w = 7, h = 6, raster_points = FALSE) {
  if (raster_points) p <- ggrastr::rasterise(p, layers = "Point", dpi = 300, dev = "ragg")
  ggsave(paste0(path, ".svg"), p, width = w, height = h, device = svglite::svglite)
  ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 150, bg = "white")
}
umap_sub <- function(obj, group, label = FALSE, title = group) {
  cols <- if (group == "tissue") tissue_cols else name_cols(obj@meta.data[[group]])
  DimPlot(obj, reduction = "umap_sub", group.by = group, label = label, repel = TRUE,
          raster = FALSE, pt.size = 0.2, cols = cols, combine = FALSE)[[1]] +
    ggtitle(title) + coord_fixed()
}
res_col <- function(r) paste0("RNA_snn_res.", r)

add_scores <- function(obj, panels, prefix) {
  panels <- lapply(panels, intersect, rownames(obj))
  panels <- panels[lengths(panels) >= 2]
  obj <- AddModuleScore(obj, features = panels, name = paste0(prefix, "tmp"), seed = 1)
  old <- paste0(prefix, "tmp", seq_along(panels))
  colnames(obj@meta.data)[match(old, colnames(obj@meta.data))] <- paste0(prefix, names(panels))
  obj
}


hbv <- readRDS("results/HBV_HCC_annotated_cnv.rds")
ct  <- as.character(hbv$celltype)
compartment <- rep(NA_character_, ncol(hbv))
for (cp in names(COMPARTMENTS)) compartment[ct %in% COMPARTMENTS[[cp]]] <- cp

if (ASSIGN_PROLIFERATING) {
  pr  <- which(ct == "Proliferating" & hbv$cnv_call != "malignant")
  pro <- subset(hbv, cells = colnames(hbv)[pr])
  pro <- add_scores(pro, LINEAGE, "LIN_")
  sc  <- as.matrix(pro@meta.data[, paste0("LIN_", names(LINEAGE))])
  colnames(sc) <- names(LINEAGE)
  srt  <- t(apply(sc, 1, sort, decreasing = TRUE))
  best <- colnames(sc)[max.col(sc, ties.method = "first")]
  lin2comp <- setNames(rep(names(COMP_LINEAGE), lengths(COMP_LINEAGE)), unlist(COMP_LINEAGE))
  ok   <- best %in% names(lin2comp) & srt[, 1] > 0.2 & (srt[, 1] - srt[, 2]) > 0.1
  compartment[pr[ok]] <- lin2comp[best[ok]]
  cat("\nProliferating cells (non-malignant):", length(pr), "\n")
  print(table(best_lineage = best, moved = ok))
  fwrite(data.table(cell = colnames(pro), best_lineage = best, moved = ok, sc),
         file.path(out_root, "proliferating_assignment.csv"))
  rm(pro); gc()
}
hbv$compartment <- compartment
print(table(hbv$compartment, hbv$tissue, useNA = "ifany"))


build_compartment <- function(cp) {
  dir   <- file.path(out_root, cp); fd <- file.path(dir, "figures")
  dir.create(fd, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(dir, paste0(cp, ".rds"))
  
  if (file.exists(cache)) {
    cat("\n==", cp, ": cached object loaded\n")
    obj <- readRDS(cache)
  } else {
    cat("\n==", cp, ": subclustering", sum(hbv$compartment == cp, na.rm = TRUE), "cells\n")
    obj <- subset(hbv, cells = colnames(hbv)[which(hbv$compartment == cp)])
    obj <- DietSeurat(obj, assays = "RNA", layers = c("counts", "data"),
                      dimreducs = NULL, graphs = NULL)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = N_HVG + 1000, verbose = FALSE)
    hv  <- head(VariableFeatures(obj)[!hvg_exclude(VariableFeatures(obj))], N_HVG)
    VariableFeatures(obj) <- hv
    obj <- ScaleData(obj, features = hv, verbose = FALSE)
    obj <- RunPCA(obj, features = hv, npcs = 50, verbose = FALSE)
    obj <- RunHarmony(obj, group.by.vars = "patient", reduction.use = "pca",
                      dims.use = 1:NPCS[[cp]], verbose = FALSE)
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:NPCS[[cp]], verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:NPCS[[cp]],
                   reduction.name = "umap_sub", seed.use = 1, verbose = FALSE)
  }
  
  
  if (!isTRUE(obj@misc$inherited_fixed)) {
    inh <- intersect(grep("^RNA_snn_res\\.", colnames(hbv@meta.data), value = TRUE),
                     colnames(obj@meta.data))
    inh <- c(inh, intersect("seurat_clusters", colnames(obj@meta.data)))
    if (length(inh)) {
      obj@meta.data <- obj@meta.data[, setdiff(colnames(obj@meta.data), inh)]
      for (r in sub("RNA_snn_res.", "", grep("^RNA_snn_res", inh, value = TRUE), fixed = TRUE)) {
        rx    <- gsub(".", "\\.", r, fixed = TRUE)
        stale <- c(list.files(fd,  paste0("_res", rx, "\\.(png|svg)$"), full.names = TRUE),
                   list.files(dir, paste0("_res", rx, "\\.csv$"), full.names = TRUE))
        file.remove(stale)
      }
      file.remove(list.files(fd, "^clustree\\.", full.names = TRUE))
      cat(cp, ": dropped inherited global clustering(s):", paste(inh, collapse = ", "), "\n")
    }
    obj@misc$inherited_fixed <- TRUE
  }
  
  for (g in c("tissue", "patient", "dataset", "subtype")) {
    f <- file.path(fd, paste0("UMAP_", g))
    if (!file.exists(paste0(f, ".png")))
      save_fig(umap_sub(obj, g, title = paste(cp, "-", g)), f, 8.5, 6.5, raster_points = TRUE)
  }
  
  
  obj@meta.data <- obj@meta.data[, !grepl("^SC_|^LIN_", colnames(obj@meta.data))]
  obj <- add_scores(obj, c(PANELS[[cp]], STATE), "SC_")
  obj <- add_scores(obj, LINEAGE, "LIN_")
  
  
  miss <- RES_GRID[!res_col(RES_GRID) %in% colnames(obj@meta.data)]
  if (length(miss))
    obj <- FindClusters(obj, graph.name = "RNA_snn", resolution = miss,
                        random.seed = 1, verbose = FALSE)
  saveRDS(obj, cache)
  
  
  genes <- intersect(unique(unlist(lapply(c(PANELS[[cp]], STATE["Cycling"]), head, 4))),
                     rownames(obj))
  for (r in RES_GRID) {
    f <- file.path(fd, paste0("UMAP_res", r))
    if (!file.exists(paste0(f, ".png")))
      save_fig(umap_sub(obj, res_col(r), label = TRUE, title = paste(cp, "res", r)),
               f, 7.5, 6.5, raster_points = TRUE)
    f <- file.path(fd, paste0("DotPlot_panels_res", r))
    if (!file.exists(paste0(f, ".png"))) {
      Idents(obj) <- res_col(r)
      p <- DotPlot(obj, features = genes, cols = c("grey88", "#DC0000")) + RotatedAxis() +
        labs(x = NULL, y = "Subcluster") + theme(axis.text.x = element_text(size = 7))
      save_fig(p, f, 0.22 * length(genes) + 3, 0.3 * nlevels(Idents(obj)) + 2)
    }
  }
  
  f <- file.path(fd, "clustree")
  if (requireNamespace("clustree", quietly = TRUE) && !file.exists(paste0(f, ".png"))) {
    ct_plot <- clustree::clustree(obj@meta.data, prefix = "RNA_snn_res.")
    ok <- tryCatch({ save_fig(ct_plot, f, 9, 11); TRUE }, error = function(e) FALSE)
    if (!ok)          
      ok <- tryCatch({
        save_fig(ct_plot + ggplot2::guides(edge_colour = "none", edge_alpha = "none"), f, 9, 11)
        TRUE
      }, error = function(e) { message("clustree skipped: ", conditionMessage(e)); FALSE })
  }
  obj
}


suggest_table <- function(obj, cp) {
  dir <- file.path(out_root, cp); fd <- file.path(dir, "figures")
  r   <- RES[[cp]]
  Idents(obj) <- res_col(r)
  md  <- as.data.table(obj@meta.data, keep.rownames = "cell")
  md[, cluster := as.character(get(res_col(r)))]
  sc_cols  <- grep("^SC_",  names(md), value = TRUE)
  lin_cols <- grep("^LIN_", names(md), value = TRUE)
  
  tab <- md[, c(list(n_cells         = .N,
                     pct_tumor       = round(100 * mean(tissue == "Tumor"), 1),
                     n_patients_ge10 = sum(table(patient) >= 10),
                     top_patient_pct = round(100 * max(table(patient)) / .N, 1),
                     median_genes    = as.numeric(median(nFeature_RNA)),
                     median_mt       = round(median(percent.mt), 1),
                     old_subtype     = names(sort(table(subtype), decreasing = TRUE))[1]),
                lapply(.SD, function(x) round(mean(x), 3))),
            by = cluster, .SDcols = c(sc_cols, lin_cols)]
  tab <- tab[order(as.numeric(cluster))]
  
  id_cols <- paste0("SC_", intersect(names(PANELS[[cp]]), sub("^SC_", "", sc_cols)))
  m   <- as.matrix(tab[, ..id_cols]); colnames(m) <- sub("^SC_", "", id_cols)
  srt <- t(apply(m, 1, sort, decreasing = TRUE))
  tab[, suggested := colnames(m)[max.col(m, ties.method = "first")]]
  tab[, id_gap := round(srt[, 1] - srt[, 2], 3)]
  
  z <- function(x) if (sd(x) > 0) (x - mean(x)) / sd(x) else x * 0
  flag <- vector("list", nrow(tab))
  add  <- function(i, txt) flag[[i]] <<- c(flag[[i]], txt)
  for (st in c("Cycling", "Stress", "ISG", "Ambient")) {
    v <- tab[[paste0("SC_", st)]]
    if (!is.null(v)) for (i in which(z(v) > 1.5 & v > 0.1)) add(i, st)
  }
  own <- COMP_LINEAGE[[cp]]
  for (o in setdiff(names(LINEAGE), own))
    for (i in which(tab[[paste0("LIN_", o)]] > 0.15)) add(i, paste0("doublet?:", o))
  if (length(own) > 1) {               # e.g. Stromal: both EC and mesenchymal programmes
    both <- Reduce(`&`, lapply(own, function(o) tab[[paste0("LIN_", o)]] > 0.15))
    for (i in which(both)) add(i, paste0("doublet?:", paste(own, collapse = "+")))
  }
  tab[, lineage := own[max.col(as.matrix(tab[, paste0("LIN_", own), with = FALSE]))]]
  for (i in which(tab$top_patient_pct > 80)) add(i, "patient-specific")
  for (i in which(tab$median_genes < 0.6 * median(tab$median_genes))) add(i, "low genes")
  for (i in which(tab$id_gap < 0.05)) add(i, "low confidence")
  tab[, flags := vapply(flag, function(f) paste(f, collapse = "; "), "")]
  
  deg <- as.data.table(FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                                      logfc.threshold = 0.5, max.cells.per.ident = 1000,
                                      random.seed = 1, verbose = FALSE))
  fwrite(deg, file.path(dir, paste0("auto_markers_res", r, ".csv")))
  top <- deg[!hvg_exclude(gene)][order(-avg_log2FC),
                                 .(top_genes = paste(head(gene, 12), collapse = ", ")), by = cluster]
  tab <- merge(tab, top[, .(cluster = as.character(cluster), top_genes)], by = "cluster",
               all.x = TRUE)
  tab <- tab[order(as.numeric(cluster))]
  tab[, `:=`(resolution = r, final = suggested, note = "")]
  setcolorder(tab, c("resolution", "cluster", "n_cells", "lineage", "suggested", "flags", "final",
                     "note", "top_genes", "old_subtype", "pct_tumor", "n_patients_ge10",
                     "top_patient_pct", "median_genes", "median_mt", "id_gap"))
  fwrite(tab, file.path(dir, paste0("subcluster_annotation_auto_res", r, ".csv")))
  
  hm <- melt(data.table(cluster = tab$cluster, scale(as.matrix(tab[, c(sc_cols, lin_cols),
                                                                   with = FALSE]))), id.vars = "cluster", variable.name = "panel", value.name = "z")
  hm[, cluster := factor(cluster, levels = tab$cluster)]
  p <- ggplot(hm, aes(panel, cluster, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35") +
    theme_minimal(base_size = 9) + labs(x = NULL, y = "Subcluster", fill = "Scaled\nscore") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  save_fig(p, file.path(fd, paste0("ScoreHeatmap_res", r)),
           0.3 * length(c(sc_cols, lin_cols)) + 3, 0.3 * nrow(tab) + 2)
  cat("Written: subcluster_annotation_auto_res", r, ".csv for ", cp, "\n", sep = "")
  print(tab[, .(cluster, n_cells, suggested, flags, old_subtype, pct_tumor)])
  invisible(tab)
}


apply_labels <- function(obj, cp) {
  dir <- file.path(out_root, cp); fd <- file.path(dir, "figures")
  lab <- fread(file.path(dir, "subcluster_annotation_final.csv"), colClasses = "character")
  r   <- as.numeric(unique(lab$resolution))
  if (length(r) != 1) stop(cp, ": final table must have one resolution")
  if (!res_col(r) %in% colnames(obj@meta.data))
    stop(cp, ": resolution ", r, " not computed; add it to RES_GRID")
  cl   <- as.character(obj@meta.data[[res_col(r)]])
  miss <- setdiff(unique(cl), lab$cluster)
  if (length(miss)) stop(cp, ": clusters missing in final table: ", paste(miss, collapse = ", "))
  obj$subtype_fine <- trimws(lab$final[match(cl, lab$cluster)])
  
  removed <- colnames(obj)[grepl("^Remove", obj$subtype_fine)]
  cat("\n==", cp, ": labels applied (res", r, ");", length(removed), "cells removed\n")
  obj <- subset(obj, cells = setdiff(colnames(obj), removed))
  
  
  show <- MARKERS_SHOW[[cp]]
  labs_ <- unique(obj$subtype_fine)
  ord  <- c(intersect(names(show), labs_), sort(setdiff(labs_, names(show))))
  obj$subtype_fine <- factor(obj$subtype_fine, levels = ord)
  Idents(obj) <- "subtype_fine"
  print(table(obj$subtype_fine, obj$tissue))
  
  
  deg <- as.data.table(FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                                      logfc.threshold = 0.5, max.cells.per.ident = 1000,
                                      random.seed = 1, verbose = FALSE))
  fwrite(deg, file.path(dir, "final_markers_by_subtype.csv"))
  
  
  pick <- lapply(ord, function(l) {
    g     <- intersect(show[[l]], rownames(obj))
    extra <- deg[as.character(cluster) == l & !hvg_exclude(gene)][order(-avg_log2FC)]$gene
    head(unique(c(g, extra)), max(N_MARK, length(g)))
  })
  names(pick) <- ord
  fwrite(data.table(subtype = rep(ord, lengths(pick)), gene = unlist(pick)),
         file.path(dir, "final_dotplot_genes.csv"))
  genes <- unique(unlist(pick))
  Idents(obj) <- factor(obj$subtype_fine, levels = rev(ord))     # first subtype on top
  p <- DotPlot(obj, features = genes, cols = c("grey88", "#DC0000")) + RotatedAxis() +
    labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 7, face = "italic"))
  save_fig(p, file.path(fd, "Final_DotPlot_markers"), 0.2 * length(genes) + 4,
           0.32 * length(ord) + 2)
  Idents(obj) <- "subtype_fine"
  
  
  cnt <- as.data.table(obj@meta.data)[, .N, by = .(sample, patient, tissue, subtype_fine)]
  fwrite(dcast(cnt, sample + patient + tissue ~ subtype_fine, value.var = "N", fill = 0),
         file.path(dir, "final_counts_by_sample.csv"))
  
  
  save_fig(umap_sub(obj, "subtype_fine", label = TRUE, title = cp),
           file.path(fd, "Final_UMAP_subtype"), 9, 7, raster_points = TRUE)
  save_fig(umap_sub(obj, "subtype_fine", title = cp) + NoLegend(),
           file.path(fd, "Final_UMAP_subtype_nolegend"), 7, 7, raster_points = TRUE)
  p <- DimPlot(obj, reduction = "umap_sub", group.by = "subtype_fine", split.by = "tissue",
               cols = name_cols(obj$subtype_fine), raster = FALSE, pt.size = 0.2,
               combine = FALSE)[[1]]
  save_fig(p, file.path(fd, "Final_UMAP_subtype_by_tissue"), 13, 6.5, raster_points = TRUE)
  
  
  cnt[, tot := sum(N), by = sample]
  fr <- CJ(sample = unique(cnt$sample), subtype_fine = unique(as.character(cnt$subtype_fine)))
  fr <- merge(fr, unique(cnt[, .(sample, patient, tissue, tot)]), by = "sample")
  fr <- merge(fr, cnt[, .(sample, subtype_fine = as.character(subtype_fine), N)],
              by = c("sample", "subtype_fine"), all.x = TRUE)
  fr[is.na(N), N := 0][, frac := N / tot][, subtype_fine := factor(subtype_fine, levels = ord)]
  p <- ggplot(fr[tot >= 20], aes(subtype_fine, frac, fill = tissue)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 0.8) +
    scale_fill_manual(values = tissue_cols) +
    labs(x = NULL, y = paste("Fraction of", cp, "cells per sample")) +
    theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p, file.path(fd, "Final_fraction_by_tissue"), 0.5 * length(ord) + 3, 5)
  
  list(obj = obj, removed = removed)
}


subtype_fine <- setNames(sub(" \\(.*", "", as.character(hbv$subtype)), colnames(hbv))
removed_all  <- character(0)
finished     <- character(0)

for (cp in names(COMPARTMENTS)) {
  obj <- build_compartment(cp)
  dir <- file.path(out_root, cp)
  if (!file.exists(file.path(dir, paste0("subcluster_annotation_auto_res", RES[[cp]], ".csv"))))
    suggest_table(obj, cp)
  if (file.exists(file.path(dir, "subcluster_annotation_final.csv"))) {
    res <- apply_labels(obj, cp)
    subtype_fine[colnames(res$obj)] <- as.character(res$obj$subtype_fine)
    removed_all <- c(removed_all, res$removed)
    finished    <- c(finished, cp)
    saveRDS(res$obj, file.path(dir, paste0(cp, "_annotated.rds")))
  } else {
    cat(cp, ": no subcluster_annotation_final.csv yet -> old subtype kept\n")
  }
  rm(obj); gc()
}

if (length(finished)) {
  hbv$subtype_fine <- unname(subtype_fine[colnames(hbv)])
  fwrite(data.table(cell = removed_all), file.path(out_root, "removed_cells.csv"))
  hbv <- subset(hbv, cells = setdiff(colnames(hbv), removed_all))
  fwrite(as.data.table(table(subtype_fine = hbv$subtype_fine, tissue = hbv$tissue)),
         file.path(out_root, "subtype_fine_counts_by_tissue.csv"))
  saveRDS(hbv, "results/HBV_HCC_subclustered.rds")
  cat("\nFinished compartments:", paste(finished, collapse = ", "),
      "\nRemoved cells:", length(removed_all),
      "\nSaved results/HBV_HCC_subclustered.rds (column subtype_fine)\n")
} else {
  cat("\nPass 1 done. Review results/subcluster/<compartment>/ and create",
      "subcluster_annotation_final.csv for each compartment, then source again.\n")
}