suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2); library(survival)
})

BASE <- path.expand("~/LIRI-JP/analysis")
OUT  <- file.path(BASE, "MDK"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
GENES <- c("MDK", "SDC1", "SDC2", "SDC3", "SDC4", "LRP1", "NCL", "ITGA4", "ITGA6",
           "ITGB1", "NOTCH2", "PTPRZ1", "SORL1")
VCOL <- c(HBV = "#E64B35", HCV = "#3C5488", NBNC = "#00A087")
TCOL <- c("Non-tumour" = "#4DBBD5", "Tumour" = "#E64B35")
theme_fig <- theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

info <- read_csv(file.path(BASE, "LIRIJP_sample_info.csv"), show_col_types = FALSE, guess_max = 1e5)
if (all(is.na(info$histology))) stop("histology is empty: put virus_status.csv in ~/Desktop/LIRI-JP and rerun 01_LIRIJP_build.R")
fpkm <- readRDS(file.path(BASE, "LIRIJP_fpkm.rds"))

tpm <- sweep(fpkm, 2, colSums(fpkm), "/") * 1e6
miss <- setdiff(GENES, rownames(tpm)); if (length(miss)) message("Not in matrix: ", paste(miss, collapse = ", "))
G <- intersect(GENES, rownames(tpm))
lx <- log2(tpm[G, , drop = FALSE] + 1)

s <- info %>% filter(histology == "HCC", tissue %in% c("Tumour", "Non-tumour"), icgc_sample_id %in% colnames(lx))
long <- as_tibble(t(lx[, s$icgc_sample_id, drop = FALSE]), rownames = "icgc_sample_id") %>%
  pivot_longer(-icgc_sample_id, names_to = "gene", values_to = "expr") %>%
  left_join(s %>% dplyr::select(icgc_sample_id, icgc_donor_id, submitted_donor_id, tissue, virus_group,
                                stage_group, OS_event, OS_months), by = "icgc_sample_id") %>%
  group_by(icgc_donor_id, submitted_donor_id, tissue, virus_group, stage_group, OS_event, OS_months, gene) %>%
  summarise(expr = mean(expr), n_samples = n(), .groups = "drop") %>%
  mutate(gene = factor(gene, GENES))
write_csv(long, file.path(OUT, "expr_MDK_receptors_donor_level.csv"))

n_tab <- long %>% filter(gene == "MDK") %>% count(tissue, virus_group)
write_csv(n_tab, file.path(OUT, "n_by_tissue_virus.csv")); print(as.data.frame(n_tab))

V3 <- c("HBV", "HCV", "NBNC")
sig_of <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "n.s.")


paired <- long %>% dplyr::select(icgc_donor_id, virus_group, gene, tissue, expr) %>%
  pivot_wider(names_from = tissue, values_from = expr) %>%
  filter(!is.na(Tumour), !is.na(`Non-tumour`)) %>% mutate(d = Tumour - `Non-tumour`)

pstat <- function(df) df %>% summarise(n_pairs = n(), median_log2FC = median(d), n_up = sum(d > 0),
                                       p = if (n() >= 5) wilcox.test(Tumour, `Non-tumour`, paired = TRUE, exact = FALSE)$p.value else NA_real_,
                                       .groups = "drop")
statA <- bind_rows(
  paired %>% filter(virus_group %in% c(V3, "HBV+HCV")) %>% group_by(gene) %>% pstat() %>% mutate(group = "All HCC"),
  paired %>% filter(virus_group %in% V3) %>% group_by(virus_group, gene) %>% pstat() %>% rename(group = virus_group)) %>%
  group_by(group) %>% mutate(padj = p.adjust(p, "BH"), sig = sig_of(padj)) %>% ungroup() %>%
  dplyr::select(group, gene, everything())
# is the tumour-vs-normal change different between viruses?
statA_virus <- paired %>% filter(virus_group %in% V3) %>% group_by(gene) %>%
  summarise(p_kruskal_FC_by_virus = kruskal.test(d ~ factor(virus_group))$p.value,
            median_FC_HBV = median(d[virus_group == "HBV"]), median_FC_HCV = median(d[virus_group == "HCV"]),
            median_FC_NBNC = median(d[virus_group == "NBNC"]),
            p_FC_HBV_vs_nonHBV = wilcox.test(d[virus_group == "HBV"], d[virus_group != "HBV"])$p.value,
            .groups = "drop") %>%
  mutate(padj_kruskal = p.adjust(p_kruskal_FC_by_virus, "BH"), padj_HBV_vs_nonHBV = p.adjust(p_FC_HBV_vs_nonHBV, "BH"))
write_csv(statA, file.path(OUT, "A_paired_TvsN_by_virus.csv"))
write_csv(statA_virus, file.path(OUT, "A_TvsN_change_differs_by_virus.csv"))

SHOW <- intersect(c("MDK", "SDC2", "SDC3", "LRP1", "NCL", "ITGA6"), G)
dA <- paired %>% filter(virus_group %in% V3, gene %in% SHOW) %>%
  pivot_longer(c(`Non-tumour`, Tumour), names_to = "tissue", values_to = "expr") %>%
  mutate(tissue = factor(tissue, c("Non-tumour", "Tumour")), virus_group = factor(virus_group, V3),
         gene = factor(gene, SHOW))
labA <- statA %>% filter(group %in% V3, gene %in% SHOW) %>%
  mutate(virus_group = factor(group, V3), gene = factor(gene, SHOW), lab = sprintf("%s %d/%d", sig, n_up, n_pairs))
pA <- ggplot(dA, aes(tissue, expr)) +
  geom_line(aes(group = icgc_donor_id), colour = "grey80", linewidth = 0.15) +
  geom_boxplot(aes(fill = tissue), width = 0.5, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
  geom_text(data = labA, aes(x = 1.5, y = Inf, label = lab), vjust = 1.2, size = 2, inherit.aes = FALSE) +
  facet_grid(gene ~ virus_group, scales = "free_y") +
  scale_x_discrete(labels = c("N", "T")) + scale_fill_manual(values = TCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "log2(TPM + 1)", title = "LIRI-JP HCC: paired tumour vs non-tumour by virus (paired Wilcoxon, BH)") + theme_fig
ggsave(file.path(OUT, "A_paired_TvsN_by_virus.pdf"), pA, width = 5, height = 1.6 * length(SHOW) + 0.6)


bstat <- function(df) df %>% filter(virus_group %in% V3) %>% group_by(gene) %>%
  summarise(n_HBV = sum(virus_group == "HBV"), n_HCV = sum(virus_group == "HCV"), n_NBNC = sum(virus_group == "NBNC"),
            median_HBV = median(expr[virus_group == "HBV"]), median_HCV = median(expr[virus_group == "HCV"]),
            median_NBNC = median(expr[virus_group == "NBNC"]),
            p_kruskal = kruskal.test(expr ~ factor(virus_group))$p.value,
            p_HBV_vs_HCV  = wilcox.test(expr[virus_group == "HBV"], expr[virus_group == "HCV"])$p.value,
            p_HBV_vs_NBNC = wilcox.test(expr[virus_group == "HBV"], expr[virus_group == "NBNC"])$p.value,
            p_HBV_vs_nonHBV = wilcox.test(expr[virus_group == "HBV"], expr[virus_group != "HBV"])$p.value,
            .groups = "drop") %>%
  mutate(across(starts_with("p_"), ~ p.adjust(.x, "BH"), .names = "padj_{sub('p_', '', .col)}"))
statB_T <- bstat(long %>% filter(tissue == "Tumour"))
statB_N <- bstat(long %>% filter(tissue == "Non-tumour"))
write_csv(statB_T, file.path(OUT, "B_tumour_by_virus.csv"))
write_csv(statB_N, file.path(OUT, "B_nontumour_by_virus.csv"))


mdk_t <- long %>% filter(tissue == "Tumour", gene == "MDK", virus_group %in% V3) %>%
  mutate(virus_group = relevel(factor(virus_group), "HBV"), stage_group = factor(stage_group))
lm_adj <- summary(lm(expr ~ virus_group + stage_group, mdk_t))$coefficients
write_csv(as_tibble(lm_adj, rownames = "term"), file.path(OUT, "B_MDK_tumour_virus_lm_stage_adjusted.csv"))

dB <- long %>% filter(virus_group %in% V3, gene %in% SHOW) %>%
  mutate(virus_group = factor(virus_group, V3), gene = factor(gene, SHOW),
         tissue = factor(tissue, c("Non-tumour", "Tumour")))
labB <- bind_rows(statB_T %>% mutate(tissue = "Tumour"), statB_N %>% mutate(tissue = "Non-tumour")) %>%
  filter(gene %in% SHOW) %>%
  mutate(gene = factor(gene, SHOW), tissue = factor(tissue, c("Non-tumour", "Tumour")),
         lab = sprintf("HBV vs HCV p=%.2g\nHBV vs NBNC p=%.2g", p_HBV_vs_HCV, p_HBV_vs_NBNC))
pB <- ggplot(dB, aes(virus_group, expr, fill = virus_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.5) +
  geom_text(data = labB, aes(x = 2, y = Inf, label = lab), vjust = 1.1, size = 1.8, inherit.aes = FALSE) +
  facet_grid(gene ~ tissue, scales = "free_y") +
  scale_fill_manual(values = VCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = NULL, y = "log2(TPM + 1)", title = "LIRI-JP HCC: expression by virus (Wilcoxon, unadjusted p)") + theme_fig
ggsave(file.path(OUT, "B_by_virus.pdf"), pB, width = 4.5, height = 1.6 * length(SHOW) + 0.6)


wide <- long %>% filter(tissue == "Tumour") %>% dplyr::select(icgc_donor_id, virus_group, stage_group, OS_event, OS_months, gene, expr) %>%
  pivot_wider(names_from = gene, values_from = expr)
statC <- bind_rows(lapply(setdiff(G, "MDK"), function(g) {
  ct <- cor.test(wide$MDK, wide[[g]], method = "spearman", exact = FALSE)
  tibble(receptor = g, rho = unname(ct$estimate), p = ct$p.value)
})) %>% mutate(padj = p.adjust(p, "BH"))
write_csv(statC, file.path(OUT, "C_MDK_receptor_correlation_tumour.csv"))


z <- function(x) as.numeric(scale(x))
sv <- wide %>% filter(!is.na(OS_months), !is.na(OS_event), OS_months > 0) %>%
  mutate(stage_group = factor(stage_group), HBV = factor(ifelse(virus_group == "HBV", "HBV", "non-HBV"), c("non-HBV", "HBV")),
         `MDK-SDC2 score` = (z(MDK) + z(SDC2)) / 2, `MDK-LRP1 score` = (z(MDK) + z(LRP1)) / 2,
         `MDK-SDC3 score` = (z(MDK) + z(SDC3)) / 2)
FEATS <- intersect(c("MDK", "SDC2", "SDC3", "LRP1", "NCL", "ITGA6", "MDK-SDC2 score", "MDK-LRP1 score", "MDK-SDC3 score"), names(sv))
cox1 <- function(df, f, adj, grp) {
  df <- df %>% mutate(x = z(.data[[f]]))
  fm <- if (adj) Surv(OS_months, OS_event) ~ x + stage_group else Surv(OS_months, OS_event) ~ x
  m <- tryCatch(coxph(fm, df), error = function(e) NULL); if (is.null(m)) return(NULL)
  s <- summary(m)
  tibble(group = grp, feature = f, model = ifelse(adj, "+ stage", "univariate"), n = s$n, events = s$nevent,
         HR = s$conf.int["x", 1], lo = s$conf.int["x", 3], hi = s$conf.int["x", 4], p = s$coefficients["x", 5])
}
groups <- list(`All HCC` = sv, HBV = sv %>% filter(virus_group == "HBV"),
               HCV = sv %>% filter(virus_group == "HCV"), NBNC = sv %>% filter(virus_group == "NBNC"))
coxD <- bind_rows(lapply(names(groups), function(g) bind_rows(lapply(FEATS, function(f)
  bind_rows(cox1(groups[[g]], f, FALSE, g), cox1(groups[[g]], f, TRUE, g))))))
coxD <- coxD %>% group_by(group, model) %>% mutate(padj = p.adjust(p, "BH")) %>% ungroup()
write_csv(coxD, file.path(OUT, "D_cox_OS.csv"))

inter <- summary(coxph(Surv(OS_months, OS_event) ~ z(MDK) * HBV + stage_group, sv))$coefficients
write_csv(as_tibble(inter, rownames = "term"), file.path(OUT, "D_cox_MDK_x_HBV_interaction.csv"))

message("\n=== A: paired T vs N, MDK and key receptors ===")
print(as.data.frame(statA %>% filter(gene %in% SHOW) %>% dplyr::select(group, gene, n_pairs, median_log2FC, n_up, padj, sig)))
message("\n=== A: does the T-N change differ by virus? ===")
print(as.data.frame(statA_virus %>% dplyr::select(gene, median_FC_HBV, median_FC_HCV, median_FC_NBNC, p_kruskal_FC_by_virus, padj_kruskal, p_FC_HBV_vs_nonHBV)))
message("\n=== B: tumours by virus ===")
print(as.data.frame(statB_T %>% dplyr::select(gene, n_HBV, n_HCV, n_NBNC, median_HBV, median_HCV, median_NBNC, p_HBV_vs_HCV, p_HBV_vs_NBNC, p_HBV_vs_nonHBV, padj_HBV_vs_nonHBV)))
message("\n=== B: MDK tumour ~ virus + stage (reference = HBV) ==="); print(lm_adj)
message("\n=== D: OS Cox, MDK and MDK-SDC2 ===")
print(as.data.frame(coxD %>% filter(feature %in% c("MDK", "MDK-SDC2 score"))))
message("\n=== D: MDK x HBV interaction ==="); print(inter)
message("\nOutputs in ", OUT)
nt <- long %>% filter(tissue == "Non-tumour", gene == "MDK", virus_group %in% c("HBV","HCV","NBNC")) %>%
  left_join(info %>% distinct(icgc_donor_id, paper_fibrosis_inuyama), by = "icgc_donor_id") %>%
  mutate(fibrosis = as.numeric(paper_fibrosis_inuyama))
print(table(nt$virus_group, nt$fibrosis))
summary(lm(expr ~ virus_group + fibrosis, nt))
