suppressPackageStartupMessages({
  library(survival); library(ggplot2); library(dplyr)
})

WORK_DIR   <- "./Bulk_HBV"
setwd(WORK_DIR)
GDC_DIR    <- WORK_DIR                                                  
VIRUS_XLSX <- file.path(WORK_DIR, "tcga_virus_status.xlsx")             
SC_RDS     <- "./HBV_HCC/results/HBV_HCC_subclustered.rds"   
CDR_XLSX   <- file.path(WORK_DIR, "TCGA-CDR-SupplementalTableS1.xlsx")  
OUT_DIR    <- file.path(WORK_DIR, "TCGA_deconv_LR")
N_CORES    <- max(1, parallel::detectCores() - 2)   
MAX_CELLS_PER_STATE <- 400        
MIN_CELLS_PER_STATE <- 20         
CENSOR_DAYS <- 5 * 365            
HBV_INCLUDE_COINFECTED <- FALSE   
options(timeout = 3600)           
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(file.path(OUT_DIR, "01_tcga_lihc.rds")) && !dir.exists(file.path(GDC_DIR, "TCGA-LIHC")))
  stop("TCGA download not found: expected ", file.path(GDC_DIR, "TCGA-LIHC"),
       "\n  Set GDC_DIR to the folder you copied (the one that CONTAINS the TCGA-LIHC folder)")
if (!file.exists(VIRUS_XLSX)) stop("Copy tcga_virus_status.xlsx into ", WORK_DIR)
if (!file.exists(file.path(OUT_DIR, "02_bayesprism_fractions.csv")) && !file.exists(SC_RDS))
  stop("scRNA reference not found: ", SC_RDS, "\n  Set SC_RDS to your HBV_HCC_subclustered.rds")


CD8  <- c("CD8 Tem", "CD8 Trm", "CD8 Temra")
TNK  <- c(CD8, "CD4 T", "Treg", "MAIT", "gdT", "NK", "Liver-resident NK")
TAM  <- c("SPP1+ TAM", "FOLR2+ TAM", "Kupffer")
MONO <- c("Classical monocyte", "Non-classical monocyte")
HSC  <- c("cyHSC", "Mesenchymal-like HSC", "myHSC")
EC   <- c("Tumour EC", "Arterial EC", "Venous EC", "LSEC", "Proliferating EC")
AXES <- list(
  list("GPC3-CD81",    "Malignant",                    c(TAM, MONO),        "GPC3",              "CD81"),
  list("MDK-SDC2",     "Malignant",                    c(HSC, "vSMC"),      "MDK",               "SDC2"),
  list("SPP1-CD44",    "SPP1+ TAM",                    c(TNK, MONO),        "SPP1",              "CD44"),
  list("C1QA-CD93",    c("FOLR2+ TAM", "Kupffer"),     EC,                  "C1QA",              "CD93"),
  list("FN1-ITGA6",    c(HSC, "vSMC"),                 EC,                  "FN1",               "ITGA6"),
  list("COL4A1-CD93",  c(HSC, EC),                     EC,                  "COL4A1",            "CD93"),
  list("CXCL12-CXCR4", c(EC, HSC),                     c(TNK, TAM, MONO),   "CXCL12",            "CXCR4"),
  list("HLA-I-CD8",    c(CD8, "Malignant", "Tumour EC"), CD8,               "HLA-A|HLA-B|HLA-C", "CD8A+CD8B"),
  list("HLA-E-KLRC1",  "Malignant",                    c("NK", "Liver-resident NK", CD8), "HLA-E", "KLRC1+KLRD1"),
  list("CCL21-CCR7",   "Lymphatic EC",                 c("mregDC", "Treg", "B cell", "CD4 T", "Tfh-like CD4 T"), "CCL21", "CCR7")
)


MAJOR <- c(
  setNames(rep("Malignant", 1), "Malignant"),
  setNames("Hepatocyte", "Hepatocyte"), setNames("Cholangiocyte", "Cholangiocyte"),
  setNames(rep("T/NK", 14), c(TNK, "Tfh-like CD4 T", "ILC", "Proliferating T/NK", "Stressed T")),
  setNames(rep("B/Plasma", 2), c("B cell", "Plasma cell")),
  setNames(rep("Myeloid", 12), c(TAM, MONO, "cDC1", "cDC2", "mregDC", "pDC", "Mast cell", "Proliferating myeloid", "Neutrophil")),
  setNames(rep("Endothelial", 6), c(EC, "Lymphatic EC")),
  setNames(rep("Mesenchymal", 4), c(HSC, "vSMC"))
)

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")


f1 <- file.path(OUT_DIR, "01_tcga_lihc.rds")
if (!file.exists(f1)) {
  suppressPackageStartupMessages({ library(TCGAbiolinks); library(SummarizedExperiment); library(readxl) })
  msg("Stage 1: TCGA-LIHC")
  query <- GDCquery(project = "TCGA-LIHC", data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
  GDCdownload(query, directory = GDC_DIR, files.per.chunk = 10)   # finds the copied files; downloads only missing ones
  se <- GDCprepare(query, directory = GDC_DIR)
  cd <- as.data.frame(colData(se))

  ## primary tumours, one aliquot per patient
  keep <- cd$sample_type == "Primary Tumor"
  se <- se[, keep]; cd <- cd[keep, ]
  cd <- cd[order(cd$barcode), ]
  cd <- cd[!duplicated(cd$patient), ]
  se <- se[, cd$barcode]

  ## gene symbols (sum duplicated symbols)
  sym <- rowData(se)$gene_name
  collapse <- function(m) { m <- rowsum(m, sym[match(rownames(m), rownames(se))]); m }
  counts <- collapse(assay(se, "unstranded"))
  tpm    <- collapse(assay(se, "tpm_unstrand"))
  colnames(counts) <- colnames(tpm) <- cd$patient

  ## virus status (core set, consensus calls) - same definition as your HCV script
  virus <- read_excel(VIRUS_XLSX, sheet = "Sheet1")
  virus$etiology <- "NBNC"
  virus$etiology[virus$HBV_consensus == "pos" & virus$HCV_consensus == "pos"] <- "HBV_HCV"
  virus$etiology[virus$HBV_consensus == "pos" & virus$HCV_consensus != "pos"] <- "HBV"
  virus$etiology[virus$HCV_consensus == "pos" & virus$HBV_consensus != "pos"] <- "HCV"
  virus$patient <- substr(virus$Barcode, 1, 12)

  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  ## clinical column names differ between TCGAbiolinks versions -> take the first one that exists
  pick <- function(...) {
    for (nm in c(...)) if (nm %in% colnames(cd) && !all(is.na(cd[[nm]]))) return(cd[[nm]])
    message("  clinical column not found: ", paste(c(...), collapse = " / "), " -> NA")
    rep(NA, nrow(cd))
  }
  vital <- tolower(as.character(pick("vital_status", "paper_vital_status")))
  stage <- gsub("STAGE ", "", toupper(as.character(pick("ajcc_pathologic_stage", "tumor_stage",
                                                         "paper_pathologic_stage", "pathologic_stage"))))
  age   <- num(pick("age_at_index", "paper_age_at_initial_pathologic_diagnosis", "age_at_diagnosis"))
  if (median(age, na.rm = TRUE) > 200) age <- age / 365.25          
  d_death <- num(pick("days_to_death", "paper_days_to_death"))
  d_fu    <- num(pick("days_to_last_follow_up", "paper_days_to_last_followup"))
  clin <- data.frame(
    patient   = cd$patient,
    age       = age,
    gender    = as.character(pick("gender", "paper_gender")),
    stage     = ifelse(grepl("^III|^IV", stage), "III-IV", ifelse(grepl("^I", stage), "I-II", NA)),
    os_event  = as.integer(vital == "dead"),
    os_time   = ifelse(vital == "dead", d_death, d_fu),
    stringsAsFactors = FALSE)
  msg("  clinical: OS available for", sum(!is.na(clin$os_time)), "tumours; stage for", sum(!is.na(clin$stage)))
  clin$etiology <- virus$etiology[match(clin$patient, virus$patient)]
  clin$batch    <- virus$Batch[match(clin$patient, virus$patient)]

  ## PFI from the TCGA Pan-Cancer Clinical Data Resource (recommended endpoint for LIHC)
  if (!file.exists(CDR_XLSX))
    try(download.file("https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
                      CDR_XLSX, mode = "wb"), silent = TRUE)
  if (file.exists(CDR_XLSX)) {
    cdr <- as.data.frame(read_excel(CDR_XLSX, sheet = "TCGA-CDR"))
    cdr <- cdr[cdr$type == "LIHC", ]
    m <- match(clin$patient, cdr$bcr_patient_barcode)
    clin$pfi_event <- num(cdr$PFI[m]); clin$pfi_time <- num(cdr$PFI.time[m])
    clin$os_event_cdr <- num(cdr$OS[m]); clin$os_time_cdr <- num(cdr$OS.time[m])
    use <- !is.na(clin$os_time_cdr)                      # curated TCGA-CDR OS preferred over GDC fields
    clin$os_event[use] <- clin$os_event_cdr[use]; clin$os_time[use] <- clin$os_time_cdr[use]
    if (all(is.na(clin$age))) clin$age <- num(cdr$age_at_initial_pathologic_diagnosis[m])
    if (all(is.na(clin$stage))) {
      st <- gsub("STAGE ", "", toupper(as.character(cdr$ajcc_pathologic_tumor_stage[m])))
      clin$stage <- ifelse(grepl("^III|^IV", st), "III-IV", ifelse(grepl("^I", st), "I-II", NA))
    }
  } else msg("  TCGA-CDR not available -> OS from GDC only")

  saveRDS(list(counts = counts, tpm = tpm, clin = clin), f1)
  msg("  ", ncol(counts), "primary tumours;", sum(!is.na(clin$etiology)), "with virus calls")
}
D <- readRDS(f1)
print(table(D$clin$etiology, useNA = "ifany"))


f2 <- file.path(OUT_DIR, "02_bayesprism_fractions.csv")
if (!file.exists(f2)) {
  suppressPackageStartupMessages({ library(Seurat); library(BayesPrism) })
  msg("Stage 2: reference")
  sc <- readRDS(SC_RDS)
  md <- sc@meta.data
  lab <- as.character(md$subtype_fine)
  ## same relabelling as the cell2location reference (config.py)
  lab[lab == "Proliferating" & md$compartment %in% "TNK"]     <- "Proliferating T/NK"
  lab[lab == "Proliferating" & md$compartment %in% "Myeloid"] <- "Proliferating myeloid"
  epi <- lab %in% c("Hepatocyte", "Cholangiocyte", "Proliferating hepatocyte", "Proliferating")
  mal <- epi & md$cnv_call %in% "malignant"
  lab[mal] <- "Malignant"
  drop <- epi & !mal & (md$tissue %in% "Tumor" | md$cnv_call %in% c("unresolved", "undetermined"))
  lab[lab == "Proliferating hepatocyte"] <- "Hepatocyte"
  drop <- drop | lab == "Proliferating" | is.na(lab)
  cells <- rownames(md)[!drop]; lab <- lab[!drop]

  n <- table(lab)
  keep_states <- names(n)[n >= MIN_CELLS_PER_STATE]
  set.seed(1)
  pick <- unlist(lapply(keep_states, function(s) {
    i <- which(lab == s); if (length(i) > MAX_CELLS_PER_STATE) sample(i, MAX_CELLS_PER_STATE) else i }))
  cells <- cells[pick]; state <- lab[pick]
  unmapped <- setdiff(unique(state), names(MAJOR))
  if (length(unmapped)) stop("Add these labels to MAJOR: ", paste(unmapped, collapse = ", "))
  type <- unname(MAJOR[state])
  msg("  reference:", length(cells), "cells,", length(unique(state)), "states,", length(unique(type)), "types")
  print(table(type))

  cnt <- tryCatch(GetAssayData(sc, assay = "RNA", layer = "counts"),
                  error = function(e) GetAssayData(sc, assay = "RNA", slot = "counts"))
  sc.dat <- t(as.matrix(cnt[, cells]))                      # cells x genes
  rm(sc, cnt); gc()

  sc.dat.f <- cleanup.genes(input = sc.dat, input.type = "count.matrix", species = "hs",
                            gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
                            exp.cells = 5)
  sc.dat.f <- select.gene.type(sc.dat.f, gene.type = "protein_coding")
  stat <- get.exp.stat(sc.dat = sc.dat[, colSums(sc.dat > 0) > 3], cell.type.labels = type,
                       cell.state.labels = state, pseudo.count = 0.1, cell.count.cutoff = 50,
                       n.cores = N_CORES)
  sc.dat.sig <- select.marker(sc.dat = sc.dat.f, stat = stat, pval.max = 0.01, lfc.min = 0.1)
  msg("  marker genes:", ncol(sc.dat.sig))

  bk.dat <- t(D$counts)                                     # samples x genes (raw counts)
  msg("Stage 2: BayesPrism on", nrow(bk.dat), "tumours (slow)")
  prism <- new.prism(reference = sc.dat.sig, mixture = bk.dat, input.type = "count.matrix",
                     cell.type.labels = type, cell.state.labels = state, key = "Malignant",
                     outlier.cut = 0.01, outlier.fraction = 0.1)
  bp <- run.prism(prism = prism, n.cores = N_CORES)
  saveRDS(bp, file.path(OUT_DIR, "02_bayesprism_result.rds"))
  th_state <- get.fraction(bp = bp, which.theta = "first", state.or.type = "state")
  th_type  <- get.fraction(bp = bp, which.theta = "final", state.or.type = "type")
  colnames(th_type) <- paste0("TYPE_", colnames(th_type))
  out <- data.frame(patient = rownames(th_state), th_state, th_type, check.names = FALSE)
  write.csv(out, f2, row.names = FALSE)
}
FR <- read.csv(f2, check.names = FALSE)
rownames(FR) <- FR$patient


msg("Stage 3: LR scores")
expr <- log2(D$tpm + 1)
gene_vec <- function(spec) {
  parts <- lapply(strsplit(spec, "+", fixed = TRUE)[[1]], function(p) {
    g <- intersect(strsplit(p, "|", fixed = TRUE)[[1]], rownames(expr))
    if (!length(g)) return(NULL)
    colMeans(expr[g, , drop = FALSE])
  })
  if (any(vapply(parts, is.null, TRUE))) return(NULL)
  M <- do.call(cbind, parts)
  if (ncol(M) == 1) return(M[, 1])
  exp(rowMeans(log(M + 1e-3))) - 1e-3                       # geometric mean (complex)
}
pats <- intersect(colnames(expr), FR$patient)
S <- data.frame(patient = pats)
for (ax in AXES) {
  nm <- ax[[1]]
  L <- gene_vec(ax[[4]]); R <- gene_vec(ax[[5]])
  if (is.null(L) || is.null(R)) { msg("  ", nm, ": gene missing, skipped"); next }
  snd <- intersect(ax[[2]], colnames(FR)); rcv <- intersect(ax[[3]], colnames(FR))
  lr <- sqrt(pmax(L[pats], 0) * pmax(R[pats], 0))
  S[[paste0(nm, "|LR")]] <- lr
  if (length(snd) && length(rcv)) {
    fs <- rowSums(FR[pats, snd, drop = FALSE]); fr <- rowSums(FR[pats, rcv, drop = FALSE])
    S[[paste0(nm, "|LRxCT")]]    <- lr * sqrt(fs * fr)      # same weighting as spatial s03
    S[[paste0(nm, "|sender")]]   <- fs
    S[[paste0(nm, "|receiver")]] <- fr
  } else msg("  ", nm, ": sender/receiver states not in deconvolution -> only LR expression score")
}
S <- merge(S, D$clin, by = "patient")
S$cohort_HBV <- if (HBV_INCLUDE_COINFECTED) S$etiology %in% c("HBV", "HBV_HCV") else S$etiology %in% "HBV"
write.csv(S, file.path(OUT_DIR, "03_LR_scores_per_tumour.csv"), row.names = FALSE)


msg("Stage 4: survival")
endpoints <- list(OS = c("os_time", "os_event"))
if ("pfi_time" %in% names(S) && any(!is.na(S$pfi_time))) endpoints$PFI <- c("pfi_time", "pfi_event")
censor <- function(t, e) {
  if (is.na(CENSOR_DAYS)) return(list(t = t, e = e))
  list(t = pmin(t, CENSOR_DAYS), e = ifelse(t > CENSOR_DAYS, 0L, e))
}
cohorts <- list(HBV = S$cohort_HBV, NBNC = S$etiology %in% "NBNC", HCV = S$etiology %in% "HCV",
                All_core = !is.na(S$etiology), All_TCGA = rep(TRUE, nrow(S)))
score_cols <- grep("\\|(LR|LRxCT|sender|receiver)$", names(S), value = TRUE)

cox_one <- function(d, sc, adjust) {
  d$x <- as.numeric(scale(d[[sc]]))
  f <- if (adjust) Surv(t, e) ~ x + age + stage else Surv(t, e) ~ x
  dd <- if (adjust) d[complete.cases(d[, c("t", "e", "x", "age", "stage")]), ] else d[complete.cases(d[, c("t", "e", "x")]), ]
  if (sum(dd$e) < 5 || length(unique(dd$x)) < 3) return(NULL)
  fit <- tryCatch(coxph(f, data = dd), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(n = nrow(dd), events = sum(dd$e), HR_per_SD = s$conf.int["x", 1],
             lo = s$conf.int["x", 3], hi = s$conf.int["x", 4], p = s$coefficients["x", 5])
}

res <- list()
for (ep in names(endpoints)) for (co in names(cohorts)) {
  d0 <- S[cohorts[[co]] & !is.na(S[[endpoints[[ep]][1]]]) & S[[endpoints[[ep]][1]]] > 0, ]
  cs <- censor(d0[[endpoints[[ep]][1]]], d0[[endpoints[[ep]][2]]]); d0$t <- cs$t; d0$e <- cs$e
  for (sc in score_cols) for (adj in c(FALSE, TRUE)) {
    r <- cox_one(d0, sc, adj)
    if (!is.null(r)) res[[length(res) + 1]] <- data.frame(endpoint = ep, cohort = co,
      axis = sub("\\|.*", "", sc), score = sub(".*\\|", "", sc),
      model = ifelse(adj, "adj_age_stage", "univariate"), r)
  }
}
res <- bind_rows(res) %>% group_by(endpoint, cohort, score, model) %>%
  mutate(q = p.adjust(p, "BH")) %>% ungroup() %>% arrange(endpoint, cohort, score, model, p)
write.csv(res, file.path(OUT_DIR, "04_cox_per_axis.csv"), row.names = FALSE)

## HBV-specific effect: score x HBV interaction in the core set 
inter <- list()
for (ep in names(endpoints)) {
  d0 <- S[!is.na(S$etiology) & !is.na(S[[endpoints[[ep]][1]]]) & S[[endpoints[[ep]][1]]] > 0, ]
  cs <- censor(d0[[endpoints[[ep]][1]]], d0[[endpoints[[ep]][2]]]); d0$t <- cs$t; d0$e <- cs$e
  d0$HBV <- as.integer(d0$cohort_HBV)
  for (sc in grep("\\|(LR|LRxCT)$", score_cols, value = TRUE)) {
    d0$x <- as.numeric(scale(d0[[sc]]))
    fit <- tryCatch(coxph(Surv(t, e) ~ x * HBV, data = d0), error = function(e) NULL)
    if (is.null(fit)) next
    cf <- summary(fit)$coefficients
    inter[[length(inter) + 1]] <- data.frame(endpoint = ep, axis = sub("\\|.*", "", sc),
      score = sub(".*\\|", "", sc), HR_interaction = exp(cf["x:HBV", 1]), p_interaction = cf["x:HBV", 5])
  }
}
inter <- bind_rows(inter) %>% group_by(endpoint, score) %>% mutate(q = p.adjust(p_interaction, "BH")) %>% ungroup()
write.csv(inter, file.path(OUT_DIR, "04_cox_HBV_interaction.csv"), row.names = FALSE)

## ---- figures
cohort_cols <- c(HBV = "#C0392B", NBNC = "#2E86C1", HCV = "#8E44AD", All_core = "grey40", All_TCGA = "black")
pdf(file.path(OUT_DIR, "04_forest_plots.pdf"), width = 9, height = 6)
for (ep in names(endpoints)) for (scn in c("LRxCT", "LR")) for (mdl in c("univariate", "adj_age_stage")) {
  d <- res %>% filter(endpoint == ep, score == scn, model == mdl, cohort %in% c("HBV", "NBNC", "HCV", "All_TCGA"))
  if (!nrow(d)) next
  d$axis <- factor(d$axis, levels = rev(unique(sapply(AXES, `[[`, 1))))
  d$lab <- ifelse(d$q < 0.05, "**", ifelse(d$p < 0.05, "*", ""))
  p <- ggplot(d, aes(HR_per_SD, axis, colour = cohort)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.3) +
    geom_text(aes(x = hi, label = lab), position = position_dodge(width = 0.6), hjust = -0.3, size = 3,
              show.legend = FALSE) +
    scale_x_log10() + scale_colour_manual(values = cohort_cols) +
    labs(x = "Hazard ratio per 1 SD (log scale)", y = NULL,
         title = sprintf("%s - %s score (%s)", ep, ifelse(scn == "LRxCT", "cell-type-weighted LR", "LR expression"), mdl),
         subtitle = "* p < 0.05   ** BH q < 0.05 within cohort") +
    theme_classic(base_size = 10)
  print(p)
}
dev.off()

km_df <- function(fit, strata_names) {
  s <- summary(fit, censored = TRUE)
  st <- if (is.null(s$strata)) rep(strata_names[1], length(s$time)) else sub(".*=", "", as.character(s$strata))
  first <- do.call(rbind, lapply(unique(st), function(g) data.frame(time = 0, surv = 1, strata = g)))
  rbind(first, data.frame(time = s$time, surv = s$surv, strata = st))
}
pdf(file.path(OUT_DIR, "04_KM_median_split.pdf"), width = 4.2, height = 3.8)
for (ep in names(endpoints)) for (co in c("HBV", "NBNC", "All_TCGA")) for (sc in grep("\\|LRxCT$", score_cols, value = TRUE)) {
  d <- S[cohorts[[co]] & !is.na(S[[endpoints[[ep]][1]]]) & S[[endpoints[[ep]][1]]] > 0, ]
  cs <- censor(d[[endpoints[[ep]][1]]], d[[endpoints[[ep]][2]]]); d$t <- cs$t / 30.44; d$e <- cs$e
  d$grp <- factor(ifelse(d[[sc]] > median(d[[sc]], na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
  if (sum(d$e) < 5) next
  fit <- survfit(Surv(t, e) ~ grp, data = d)
  lr <- survdiff(Surv(t, e) ~ grp, data = d)
  p <- 1 - pchisq(lr$chisq, 1)
  kd <- km_df(fit, levels(d$grp))
  kd$strata <- factor(kd$strata, levels = c("Low", "High"))
  nn <- table(d$grp)
  g <- ggplot(kd, aes(time, surv, colour = strata)) + geom_step(linewidth = 0.7) +
    scale_colour_manual(values = c(Low = "#2E86C1", High = "#C0392B"), breaks = c("Low", "High"),
                        labels = c(Low = sprintf("Low (n=%d)", nn[["Low"]]), High = sprintf("High (n=%d)", nn[["High"]]))) +
    scale_y_continuous(limits = c(0, 1)) +
    annotate("text", x = 0, y = 0.05, hjust = 0, size = 3.2, label = sprintf("log-rank p = %.3g", p)) +
    labs(x = "Months", y = paste(ep, "probability"), colour = NULL,
         title = sub("\\|LRxCT", "", sc), subtitle = sprintf("%s, median split of cell-type-weighted LR score", co)) +
    theme_classic(base_size = 9) + theme(legend.position = c(0.75, 0.9))
  print(g)
}
dev.off()

## deconvolution overview by etiology
tcols <- grep("^TYPE_", colnames(FR), value = TRUE)
long <- data.frame(patient = rep(FR$patient, length(tcols)),
                   type = rep(sub("TYPE_", "", tcols), each = nrow(FR)),
                   fraction = unlist(FR[, tcols], use.names = FALSE))
long$etiology <- D$clin$etiology[match(long$patient, D$clin$patient)]
long <- long[!is.na(long$etiology), ]
g <- ggplot(long, aes(etiology, fraction, fill = etiology)) + geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~type, scales = "free_y") + theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "BayesPrism fraction", title = "TCGA-LIHC deconvolution by etiology")
ggsave(file.path(OUT_DIR, "02_fractions_by_etiology.pdf"), g, width = 8, height = 6)

## console summary
cat("\nEvents per cohort (", names(endpoints)[1], "):\n")
print(res %>% filter(endpoint == names(endpoints)[1], model == "univariate") %>% group_by(cohort) %>%
        summarise(n = max(n), events = max(events)))
cat("\nTop associations in the HBV cohort (cell-type-weighted score, univariate):\n")
print(res %>% filter(cohort == "HBV", score == "LRxCT", model == "univariate") %>%
        select(endpoint, axis, n, events, HR_per_SD, lo, hi, p, q) %>% arrange(p), n = 30)
msg("Done. Outputs in", OUT_DIR)
