setwd("~/HBV_HCC")
dir.create("results", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(harmony)
})
options(future.globals.maxSize = 50 * 1024^3)

keep_sharma <- paste0("P", c(1, 2, 3, 4, 7, 9, 10, 11, 14))
keep_lu     <- c("HCC06", "HCC09", "HCC10")


check_raw <- function(m, label) {
  x  <- m@x[seq_len(min(length(m@x), 1e7))]
  is_int <- all(x == round(x))
  cs <- Matrix::colSums(m)
  cat("\n==", label, "==\n",
      " genes x cells       :", dim(m), "\n",
      " all whole numbers   :", is_int, "\n",
      " min / max value     :", min(m@x), "/", max(m@x), "\n",
      " median UMIs per cell:", median(cs), "\n")
  if (is_int && max(m@x) > 50) cat("  -> RAW counts: OK to use\n")
  else cat("  -> NOT raw counts (normalised/log): do not use as counts\n")
  invisible(is_int)
}

m1  <- readMM("GSE156625/GSE156625_HCCmatrix.mtx.gz")
g1  <- fread("GSE156625/GSE156625_HCCgenes.tsv.gz",    header = FALSE)
bc1 <- fread("GSE156625/GSE156625_HCCbarcodes.tsv.gz", header = FALSE)[[1]]
stopifnot(nrow(m1) == nrow(g1), ncol(m1) == length(bc1))

m1 <- as(m1, "CsparseMatrix")
rownames(m1) <- make.unique(g1[[ncol(g1)]])   
colnames(m1) <- bc1
stopifnot(!anyDuplicated(bc1))                


nUMI <- Matrix::colSums(m1)
m1   <- m1[, nUMI >= 500]
rm(nUMI, bc1); gc()
cat("GSE156625 barcodes with >= 500 UMIs:", ncol(m1), "\n")


print(head(colnames(m1)))
lib <- sub("^[ACGTN]+", "", colnames(m1))
print(table(lib))

map <- fread("GSE156625/sample_map.csv", colClasses = "character")
stopifnot(nrow(map) == 58, all(unique(lib) %in% map$lib))

idx   <- match(lib, map$lib)
meta1 <- data.frame(
  row.names = colnames(m1),
  lib       = lib,
  gsm       = map$gsm[idx],
  patient   = map$patient[idx],
  tissue    = map$tissue[idx],      
  sector    = map$sector[idx]
)


s1a <- data.frame(
  patient = c("P1","P2","P3","P4","P7","P9","P10","P11","P14"),
  Normal_paper = c(400, 865, 396, 1226, 1266, 636, 4234, 479, 26),
  Tumor_paper  = c(2098, 3224, 10946, 518, 3933, 5591, 556, 223, 2007))
obs <- as.data.frame.matrix(table(meta1$patient, meta1$tissue))
obs$patient <- rownames(obs)
print(merge(s1a, obs[, c("patient", "Normal", "Tumor")], by = "patient"))

sel1  <- meta1$patient %in% keep_sharma & meta1$tissue %in% c("Tumor", "Normal")
m1    <- m1[, sel1]
meta1 <- meta1[sel1, ]
print(table(meta1$patient, meta1$tissue))

check_raw(m1, "GSE156625 - HBV patients")

s1 <- CreateSeuratObject(m1, meta.data = meta1, project = "GSE156625")
s1$dataset <- "GSE156625"
rm(m1, meta1, lib, sel1); gc()


f_cnt  <- list.files("GSE149614", pattern = "count",    full.names = TRUE)[1]
f_meta <- list.files("GSE149614", pattern = "metadata", full.names = TRUE)[1]
cat("\nCounts file  :", f_cnt, "\nMetadata file:", f_meta, "\n")

meta2 <- fread(f_meta)
print(head(meta2)); print(colnames(meta2))


cell_col   <- colnames(meta2)[1]   
sample_col <- "sample"             

meta2[, patient := sub("^(HCC\\d+).*$", "\\1", get(sample_col))]
meta2[, site    := sub("^HCC\\d+", "",       get(sample_col))]   
meta2[, tissue  := fifelse(site == "T", "Tumor",
                           fifelse(site == "N", "Normal", NA_character_))]
print(table(meta2$patient, meta2$site))

sel2   <- meta2$patient %in% keep_lu & !is.na(meta2$tissue)      
cells2 <- meta2[[cell_col]][sel2]
stopifnot(length(cells2) > 0)


hdr <- names(fread(f_cnt, nrows = 0))
stopifnot(all(cells2 %in% hdr))
cnt <- fread(f_cnt, select = c(hdr[1], cells2))
m2  <- as(as.matrix(cnt[, -1]), "CsparseMatrix")
rownames(m2) <- make.unique(cnt[[1]])
rm(cnt, hdr); gc()

check_raw(m2, "GSE149614 - HBV patients")

meta2_df <- data.frame(meta2[sel2, ], row.names = cells2)
s2 <- CreateSeuratObject(m2, meta.data = meta2_df[colnames(m2), c("patient", "tissue")],
                         project = "GSE149614")
s2$dataset <- "GSE149614"
rm(m2, meta2, meta2_df); gc()


genes <- intersect(rownames(s1), rownames(s2))
cat("\nShared genes:", length(genes), "(GSE156625:", nrow(s1),
    "| GSE149614:", nrow(s2), ")\n")

hbv <- merge(subset(s1, features = genes), subset(s2, features = genes),
             add.cell.ids = c("SH", "LU"))
hbv <- JoinLayers(hbv)
hbv$sample <- paste(hbv$patient, hbv$tissue, sep = "_")
hbv[["percent.mt"]] <- PercentageFeatureSet(hbv, pattern = "^MT-")
rm(s1, s2); gc()

pdf("results/QC_before_filter.pdf", width = 14, height = 5)
print(VlnPlot(hbv, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "sample", pt.size = 0))
dev.off()


hbv <- subset(hbv, nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt < 25)
print(table(hbv$patient, hbv$tissue))


hbv[["RNA"]] <- split(hbv[["RNA"]], f = hbv$patient)

hbv <- NormalizeData(hbv)
hbv <- FindVariableFeatures(hbv, nfeatures = 2000)
hbv <- ScaleData(hbv)
hbv <- RunPCA(hbv, npcs = 50)
hbv <- RunUMAP(hbv, reduction = "pca", dims = 1:30,
               reduction.name = "umap.unintegrated")

hbv <- IntegrateLayers(hbv, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = "harmony")
hbv <- FindNeighbors(hbv, reduction = "harmony", dims = 1:30)
hbv <- FindClusters(hbv, resolution = 0.5)
hbv <- RunUMAP(hbv, reduction = "harmony", dims = 1:30, reduction.name = "umap")
hbv <- JoinLayers(hbv)

pdf("results/UMAP_before_after_integration.pdf", width = 16, height = 6)
for (g in c("dataset", "patient", "tissue")) {
  print(DimPlot(hbv, reduction = "umap.unintegrated", group.by = g) +
          DimPlot(hbv, reduction = "umap",              group.by = g))
}
print(DimPlot(hbv, reduction = "umap", label = TRUE))
dev.off()

saveRDS(hbv, "results/HBV_HCC_integrated.rds")
