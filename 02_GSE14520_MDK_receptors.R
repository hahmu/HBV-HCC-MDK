suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2); library(survival); library(stringr)
  library(affy); library(AnnotationDbi)
})
select <- dplyr::select   # AnnotationDbi masks dplyr::select

OUT  <- path.expand("~/GSE14520/analysis")
OUT2 <- file.path(OUT, "MDK"); dir.create(OUT2, showWarnings = FALSE, recursive = TRUE)
GENES <- c("MDK", "SDC1", "SDC2", "SDC3", "SDC4", "LRP1", "NCL", "ITGA4", "ITGA6",
           "ITGB1", "NOTCH2", "PTPRZ1", "SORL1")
COLS <- c("Non-tumour" = "#4DBBD5", "Tumour" = "#E64B35")
HBV_COLS <- c("AVR-CC" = "#E64B35", "CC" = "#F39B7F", "N" = "#8491B4")
theme_fig <- theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

info <- read_csv(file.path(OUT, "GSE14520_sample_info.csv"), show_col_types = FALSE) %>%
  filter(has_CEL, tissue %in% c("Tumour", "Non-tumour"))
platform_db <- c("GPL3921_HT-HG-U133A" = "hthgu133a.db", "GPL571_HG-U133A-2" = "hgu133a2.db")
MAIN <- names(which.max(table(info$platform)))
message("Main platform: ", MAIN)


rma_platform <- function(plat) {
  cache <- file.path(OUT, paste0("rma_", plat, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  s <- info %>% filter(platform == plat)
  message("RMA on ", nrow(s), " arrays of ", plat, " (the HT set takes a while and ~8 GB RAM)")
  dirs <- unique(dirname(s$cel_file))
  if (length(dirs) != 1) stop("CEL files are in several folders; move them into one: ", paste(dirs, collapse = ", "))
  missing <- s$cel_file[!file.exists(s$cel_file)]
  if (length(missing)) stop(length(missing), " CEL files not found, e.g. ", missing[1])
  # justRMA joins celfile.path + filenames, so pass the folder and bare file names separately
  m <- exprs(justRMA(filenames = basename(s$cel_file), celfile.path = dirs))
  colnames(m) <- str_extract(colnames(m), "GSM[0-9]+")
  saveRDS(m, cache); m
}

load_annot <- function(db) {
  for (pkg in unique(c(db, "hgu133a.db", "hgu133a2.db"))) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing ", pkg, " ...")
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      try(BiocManager::install(pkg, update = FALSE, ask = FALSE), silent = TRUE)
    }
    if (requireNamespace(pkg, quietly = TRUE)) {
      if (pkg != db) message("Using ", pkg, " for probe annotation (", db, " unavailable)")
      return(getExportedValue(pkg, pkg))
    }
  }
  stop("No annotation package could be installed: ", db)
}

gene_matrix <- function(m, db) {
  annot <- load_annot(db)
  ks  <- intersect(rownames(m), AnnotationDbi::keys(annot, keytype = "PROBEID"))
  message(length(ks), " of ", nrow(m), " probes found in the annotation")
  ann <- suppressMessages(AnnotationDbi::select(annot, keys = ks, columns = "SYMBOL", keytype = "PROBEID")) %>%
    filter(SYMBOL %in% GENES) %>% distinct(PROBEID, .keep_all = TRUE)
  probe <- as_tibble(m[ann$PROBEID, , drop = FALSE], rownames = "PROBEID") %>%
    left_join(ann, by = "PROBEID") %>%
    pivot_longer(starts_with("GSM"), names_to = "gsm", values_to = "expr")
  probe_stats <- probe %>% group_by(SYMBOL, PROBEID) %>%
    summarise(mean_expr = mean(expr), sd_expr = sd(expr), .groups = "drop") %>%
    group_by(SYMBOL) %>% mutate(selected = mean_expr == max(mean_expr)) %>% ungroup()
  sel <- probe_stats %>% filter(selected) %>% pull(PROBEID)      # one probe per gene: highest mean
  list(gene = probe %>% filter(PROBEID %in% sel) %>% dplyr::select(gene = SYMBOL, gsm, expr),
       probes = probe_stats)
}

L <- list(); P <- list()
for (plat in intersect(names(platform_db), unique(info$platform))) {
  m  <- rma_platform(plat)                     # run (or load cached) RMA first
  gm <- gene_matrix(m, platform_db[[plat]])
  L[[plat]] <- gm$gene %>% mutate(platform = plat)
  P[[plat]] <- gm$probes %>% mutate(platform = plat)
}
expr_long <- bind_rows(L) %>%
  left_join(info %>% dplyr::select(gsm, tissue, patient_id, paired_TN, HBV_status), by = "gsm") %>%
  mutate(gene = factor(gene, levels = GENES))
write_csv(bind_rows(P), file.path(OUT2, "probes_used.csv"))
write_csv(expr_long, file.path(OUT2, "expr_MDK_receptors_long.csv"))
missing_genes <- setdiff(GENES, unique(as.character(expr_long$gene)))
if (length(missing_genes)) message("No probe for: ", paste(missing_genes, collapse = ", "))


paired <- expr_long %>% filter(paired_TN) %>%
  dplyr::select(platform, gene, patient_id, tissue, expr) %>%
  pivot_wider(names_from = tissue, values_from = expr) %>%
  filter(!is.na(Tumour), !is.na(`Non-tumour`))

statA <- paired %>% group_by(platform, gene) %>%
  summarise(n_pairs = n(), median_log2FC = median(Tumour - `Non-tumour`),
            n_up = sum(Tumour > `Non-tumour`),
            p = wilcox.test(Tumour, `Non-tumour`, paired = TRUE, exact = FALSE)$p.value, .groups = "drop") %>%
  group_by(platform) %>% mutate(padj = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(sig = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "n.s."))
write_csv(statA, file.path(OUT2, "A_paired_TvsN_stats.csv"))

for (plat in unique(paired$platform)) {
  d <- paired %>% filter(platform == plat) %>%
    pivot_longer(c(Tumour, `Non-tumour`), names_to = "tissue", values_to = "expr") %>%
    mutate(tissue = factor(tissue, c("Non-tumour", "Tumour")))
  st <- statA %>% filter(platform == plat) %>% mutate(lab = sprintf("%s\n%d/%d up", sig, n_up, n_pairs))
  pA <- ggplot(d, aes(tissue, expr)) +
    geom_line(aes(group = patient_id), colour = "grey80", linewidth = 0.1) +
    geom_boxplot(aes(fill = tissue), width = 0.5, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
    geom_text(data = st, aes(x = 1.5, y = Inf, label = lab), vjust = 1.1, size = 2.2, inherit.aes = FALSE) +
    facet_wrap(~gene, scales = "free_y", nrow = 1) +
    scale_x_discrete(labels = c("N", "T")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    scale_fill_manual(values = COLS, guide = "none") +
    labs(x = NULL, y = "log2 RMA expression",
         title = sprintf("GSE14520 (%s): MDK and receptors, %d paired HBV-HCC tumour / non-tumour; paired Wilcoxon, BH",
                         plat, max(st$n_pairs))) + theme_fig
  ggsave(file.path(OUT2, paste0("A_paired_TvsN_", plat, ".pdf")), pA, width = 11, height = 2.8)
}


tum <- expr_long %>% filter(platform == MAIN, tissue == "Tumour")
grp_n <- tum %>% distinct(gsm, HBV_status) %>% count(HBV_status)
write_csv(grp_n, file.path(OUT2, "B_HBV_group_sizes_tumour.csv")); print(grp_n)

statB <- tum %>% filter(HBV_status %in% c("AVR-CC", "CC")) %>% group_by(gene) %>%
  summarise(n_AVR_CC = sum(HBV_status == "AVR-CC"), n_CC = sum(HBV_status == "CC"),
            median_AVR_CC = median(expr[HBV_status == "AVR-CC"]),
            median_CC = median(expr[HBV_status == "CC"]),
            diff_AVR_minus_CC = median_AVR_CC - median_CC,
            p = wilcox.test(expr[HBV_status == "AVR-CC"], expr[HBV_status == "CC"])$p.value, .groups = "drop") %>%
  mutate(padj = p.adjust(p, "BH"))
write_csv(statB, file.path(OUT2, "B_tumour_AVRCC_vs_CC_stats.csv"))

tumB <- tum %>% filter(HBV_status %in% names(HBV_COLS)) %>%
  mutate(HBV_status = factor(HBV_status, names(HBV_COLS)))
labB <- statB %>% mutate(lab = sprintf("AVR-CC vs CC\np = %.2g", p))
pB <- ggplot(tumB, aes(HBV_status, expr, fill = HBV_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.5) +
  geom_text(data = labB, aes(x = 1.5, y = Inf, label = lab), vjust = 1.1, size = 2, inherit.aes = FALSE) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  scale_fill_manual(values = HBV_COLS, guide = "none") +
  labs(x = "HBV viral status", y = "log2 RMA expression (tumour)",
       title = sprintf("GSE14520 tumours by HBV status (AVR-CC n=%d, CC n=%d, N n=%d)",
                       sum(grp_n$n[grp_n$HBV_status %in% "AVR-CC"]), sum(grp_n$n[grp_n$HBV_status %in% "CC"]),
                       sum(grp_n$n[grp_n$HBV_status %in% "N"]))) + theme_fig
ggsave(file.path(OUT2, "B_tumour_by_HBV_status.pdf"), pB, width = 9, height = 4.5)


wide <- tum %>% dplyr::select(gsm, gene, expr) %>% pivot_wider(names_from = gene, values_from = expr)
statC <- bind_rows(lapply(setdiff(intersect(GENES, names(wide)), "MDK"), function(g) {
  ct <- cor.test(wide$MDK, wide[[g]], method = "spearman", exact = FALSE)
  tibble(receptor = g, rho = unname(ct$estimate), p = ct$p.value)
})) %>% mutate(padj = p.adjust(p, "BH"))
write_csv(statC, file.path(OUT2, "C_MDK_receptor_correlation_tumour.csv"))


z <- function(x) as.numeric(scale(x))
clin <- info %>% filter(platform == MAIN, tissue == "Tumour") %>%
  dplyr::select(gsm, OS_event, OS_months, RFS_event, RFS_months, TNM_group, HBV_status)
sv <- wide %>% left_join(clin, by = "gsm") %>%
  mutate(`MDK-SDC2 score` = (z(MDK) + z(SDC2)) / 2,   # mean z-score; not the TCGA LR-score formula
         `MDK-SDC3 score` = (z(MDK) + z(SDC3)) / 2,
         `MDK-LRP1 score` = (z(MDK) + z(LRP1)) / 2,
         TNM_group = factor(TNM_group, c("I", "II", "III")))
FEATS <- intersect(c("MDK", "SDC2", "SDC3", "LRP1", "SDC4", "ITGA6", "NCL",
                     "MDK-SDC2 score", "MDK-SDC3 score", "MDK-LRP1 score"), names(sv))
ENDP <- list(OS = c("OS_months", "OS_event"), RFS = c("RFS_months", "RFS_event"))

cox_rows <- list()
pdf(file.path(OUT2, "D_KM_median_split.pdf"), width = 9, height = 6)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), cex = 0.7)
for (ep in names(ENDP)) {
  time <- sv[[ENDP[[ep]][1]]]; ev <- sv[[ENDP[[ep]][2]]]
  for (f in FEATS) {
    x <- z(sv[[f]]); ok <- !is.na(time) & !is.na(ev) & !is.na(x)
    fit_row <- function(fm, label) {
      s <- summary(fm)
      tibble(endpoint = ep, feature = f, model = label, n = s$n, events = s$nevent,
             HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5])
    }
    d <- data.frame(time, ev, x, TNM = sv$TNM_group)[ok, ]
    cox_rows[[paste(ep, f, "u")]] <- fit_row(coxph(Surv(time, ev) ~ x, d), "univariate")
    d2 <- d[!is.na(d$TNM), ]
    cox_rows[[paste(ep, f, "a")]] <- fit_row(coxph(Surv(time, ev) ~ x + TNM, d2), "+ TNM stage")
    if (f %in% c("MDK", "SDC2", "MDK-SDC2 score")) {
      d$grp <- factor(ifelse(d$x > median(d$x), "High", "Low"), c("Low", "High"))
      km <- survfit(Surv(time, ev) ~ grp, d)
      plr <- 1 - pchisq(survdiff(Surv(time, ev) ~ grp, d)$chisq, 1)
      plot(km, col = c("#4DBBD5", "#E64B35"), lwd = 1.5, xlab = "Months",
           ylab = paste(ep, "probability"), main = sprintf("%s: %s (log-rank p = %.3g)", ep, f, plr))
      legend("bottomleft", sprintf("%s (n=%d)", c("Low", "High"), as.integer(table(d$grp))),
             col = c("#4DBBD5", "#E64B35"), lwd = 1.5, bty = "n")
    }
  }
}
dev.off()
cox <- bind_rows(cox_rows) %>% group_by(endpoint, model) %>% mutate(padj = p.adjust(p, "BH")) %>% ungroup()
write_csv(cox, file.path(OUT2, "D_cox_per_SD.csv"))

cox$feature <- factor(cox$feature, rev(FEATS))
pD <- ggplot(cox, aes(HR, feature, colour = model)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), size = 1.5) +
  facet_wrap(~endpoint) + scale_x_log10() +
  scale_colour_manual(values = c("univariate" = "grey50", "+ TNM stage" = "#3C5488")) +
  labs(x = "Hazard ratio per 1 SD (log scale)", y = NULL, colour = NULL,
       title = sprintf("GSE14520 HBV-HCC tumours (%s): Cox per SD", MAIN)) + theme_fig
ggsave(file.path(OUT2, "D_cox_forest.pdf"), pD, width = 7, height = 3.5)

message("\n=== A: paired T vs N (", MAIN, ") ===")
print(as.data.frame(statA %>% filter(platform == MAIN) %>% dplyr::select(gene, n_pairs, median_log2FC, n_up, padj, sig)))
message("\n=== B: AVR-CC vs CC ===");  print(as.data.frame(statB %>% dplyr::select(gene, diff_AVR_minus_CC, p, padj)))
message("\n=== D: Cox (OS, + TNM) ==="); print(as.data.frame(cox %>% filter(endpoint == "OS", model == "+ TNM stage")))
message("\nAll outputs in ", OUT2)
