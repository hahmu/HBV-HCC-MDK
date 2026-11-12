suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2); library(survival)
})

OUT  <- path.expand("~/Desktop/GSE14520/analysis")
OUT3 <- file.path(OUT, "MDK", "clinical"); dir.create(OUT3, showWarnings = FALSE, recursive = TRUE)
GENE <- "MDK"
theme_fig <- theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

info <- read_csv(file.path(OUT, "GSE14520_sample_info.csv"), show_col_types = FALSE)
ex   <- read_csv(file.path(OUT, "MDK", "expr_MDK_receptors_long.csv"), show_col_types = FALSE)
MAIN <- names(which.max(table(ex$platform)))

d <- ex %>% filter(platform == MAIN, tissue == "Tumour", gene == GENE) %>%
  dplyr::select(gsm, MDK = expr) %>%
  left_join(info %>% filter(tissue == "Tumour"), by = "gsm") %>%
  mutate(
    MDK_z     = as.numeric(scale(MDK)),
    MDK_group = factor(ifelse(MDK > median(MDK), "MDK-high", "MDK-low"), c("MDK-low", "MDK-high")),
    MDK_tert  = cut(MDK, quantile(MDK, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("Low", "Mid", "High")),
    age_group = factor(ifelse(age > 50, ">50", "<=50"), c("<=50", ">50")),
    TNM_group = factor(TNM_group, c("I", "II", "III")),
    BCLC      = factor(BCLC, c("0", "A", "B", "C")),
    CLIP_group = factor(case_when(CLIP %in% c(0, "0") ~ "0", CLIP %in% c(1, "1") ~ "1",
                                  !is.na(CLIP) ~ ">=2"), c("0", "1", ">=2")),
    HBV_status = factor(ifelse(HBV_status %in% c("AVR-CC", "CC"), HBV_status, NA), c("CC", "AVR-CC")),
    across(c(sex, ALT, tumour_size, multinodular, cirrhosis, AFP, metastasis_signature), factor))
message(nrow(d), " tumours on ", MAIN)


VARS <- c(sex = FALSE, age_group = FALSE, HBV_status = FALSE, ALT = FALSE, cirrhosis = FALSE,
          tumour_size = FALSE, multinodular = FALSE, AFP = FALSE, metastasis_signature = FALSE,
          TNM_group = TRUE, BCLC = TRUE, CLIP_group = TRUE)
VARS <- VARS[names(VARS) %in% names(d)]


e1 <- bind_rows(lapply(names(VARS), function(v) {
  x <- d %>% filter(!is.na(.data[[v]])) %>% mutate(g = droplevels(.data[[v]]))
  k <- nlevels(x$g); if (k < 2) return(NULL)
  med <- x %>% group_by(g) %>% summarise(n = n(), med = median(MDK), .groups = "drop")
  p <- if (k == 2) wilcox.test(MDK ~ g, x)$p.value else kruskal.test(MDK ~ g, x)$p.value
  trend <- if (VARS[[v]]) cor.test(as.numeric(x$g), x$MDK, method = "spearman", exact = FALSE) else NULL
  tibble(variable = v, test = ifelse(k == 2, "Wilcoxon", "Kruskal-Wallis"),
         groups = paste(sprintf("%s (n=%d, median %.2f)", med$g, med$n, med$med), collapse = " | "),
         diff_last_minus_first = med$med[k] - med$med[1], p = p,
         trend_rho = if (!is.null(trend)) unname(trend$estimate) else NA_real_,
         trend_p = if (!is.null(trend)) trend$p.value else NA_real_)
})) %>% mutate(padj = p.adjust(p, "BH")) %>% arrange(p)
age_cor <- cor.test(d$age, d$MDK, method = "spearman", exact = FALSE)
e1 <- bind_rows(e1, tibble(variable = "age (continuous)", test = "Spearman",
                           trend_rho = unname(age_cor$estimate), p = age_cor$p.value))
write_csv(e1, file.path(OUT3, "E1_MDK_vs_clinical_tests.csv"))

long <- bind_rows(lapply(names(VARS), function(v)
  d %>% filter(!is.na(.data[[v]])) %>% transmute(variable = v, level = as.character(.data[[v]]), MDK)))
long$variable <- factor(long$variable, names(VARS))
lab <- e1 %>% filter(variable %in% names(VARS)) %>%
  mutate(variable = factor(variable, names(VARS)),
         lab = ifelse(is.na(trend_p), sprintf("p = %.2g", p), sprintf("p = %.2g\ntrend p = %.2g", p, trend_p)))
pE1 <- ggplot(long, aes(level, MDK)) +
  geom_boxplot(fill = "#E64B35", alpha = 0.6, outlier.shape = NA, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.5) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = lab), hjust = 1.05, vjust = 1.1, size = 2, inherit.aes = FALSE) +
  facet_wrap(~variable, scales = "free_x", nrow = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = NULL, y = "MDK log2 RMA (tumour)", title = sprintf("GSE14520 (%s): tumour MDK by clinical variable", MAIN)) +
  theme_fig
ggsave(file.path(OUT3, "E1_MDK_by_clinical_boxplots.pdf"), pE1, width = 10, height = 4.5)


e2 <- bind_rows(lapply(names(VARS), function(v) {
  tb <- table(d[[v]], d$MDK_group)
  if (nrow(tb) < 2) return(NULL)
  p <- tryCatch(fisher.test(tb, workspace = 2e7)$p.value, error = function(e) chisq.test(tb)$p.value)
  tibble(variable = v, level = rownames(tb),
         `MDK-low n (%)`  = sprintf("%d (%.1f%%)", tb[, 1], 100 * tb[, 1] / sum(tb[, 1])),
         `MDK-high n (%)` = sprintf("%d (%.1f%%)", tb[, 2], 100 * tb[, 2] / sum(tb[, 2])),
         p_Fisher = c(p, rep(NA, nrow(tb) - 1)))
}))
age_p <- wilcox.test(age ~ MDK_group, d)$p.value
e2 <- bind_rows(tibble(variable = "age", level = "median (IQR)",
                       `MDK-low n (%)`  = with(d[d$MDK_group == "MDK-low", ], sprintf("%.0f (%.0f-%.0f)", median(age, na.rm = TRUE), quantile(age, .25, na.rm = TRUE), quantile(age, .75, na.rm = TRUE))),
                       `MDK-high n (%)` = with(d[d$MDK_group == "MDK-high", ], sprintf("%.0f (%.0f-%.0f)", median(age, na.rm = TRUE), quantile(age, .25, na.rm = TRUE), quantile(age, .75, na.rm = TRUE))),
                       p_Fisher = age_p), e2)
write_csv(e2, file.path(OUT3, "E2_MDK_high_low_table.csv"))


ENDP <- list(OS = c("OS_months", "OS_event"), RFS = c("RFS_months", "RFS_event"))
zph_rows <- list()
pdf(file.path(OUT3, "E3_cox_zph_schoenfeld.pdf"), width = 8, height = 6)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), cex = 0.7)
for (ep in names(ENDP)) {
  dd <- d %>% mutate(time = .data[[ENDP[[ep]][1]]], ev = .data[[ENDP[[ep]][2]]]) %>%
    filter(!is.na(time), !is.na(ev))
  fits <- list(univariate = coxph(Surv(time, ev) ~ MDK_z, dd),
               `+ TNM stage` = coxph(Surv(time, ev) ~ MDK_z + TNM_group, dd %>% filter(!is.na(TNM_group))))
  for (m in names(fits)) {
    z <- cox.zph(fits[[m]])
    zt <- as.data.frame(z$table)
    zph_rows[[paste(ep, m)]] <- tibble(endpoint = ep, model = m, term = rownames(zt),
                                       chisq = zt$chisq, df = zt$df, p = zt$p)
    plot(z[1], resid = TRUE, se = TRUE, col = "#E64B35",
         main = sprintf("%s, %s: MDK Schoenfeld residuals (p = %.2g)", ep, m, zt["MDK_z", "p"]))
    abline(h = coef(fits[[m]])[1], lty = 2, col = "grey40")
  }
}
dev.off()
zph <- bind_rows(zph_rows)
write_csv(zph, file.path(OUT3, "E3_cox_zph.csv"))


covs_full <- c("TNM_group", "AFP", "tumour_size", "multinodular", "cirrhosis", "ALT")
covs_full <- covs_full[covs_full %in% names(d)]
e4_rows <- list()
for (ep in names(ENDP)) {
  dd <- d %>% mutate(time = .data[[ENDP[[ep]][1]]], ev = .data[[ENDP[[ep]][2]]]) %>%
    filter(!is.na(time), !is.na(ev))
  models <- list(`MDK + TNM` = c("TNM_group"),
                 `MDK + TNM + AFP` = c("TNM_group", "AFP"),
                 `MDK + all clinical` = covs_full,
                 `MDK + BCLC` = c("BCLC"))
  for (m in names(models)) {
    cv <- models[[m]]
    dm <- dd %>% filter(if_all(all_of(cv), ~ !is.na(.x)))
    f  <- as.formula(paste("Surv(time, ev) ~ MDK_z +", paste(cv, collapse = " + ")))
    s  <- summary(coxph(f, dm))
    ci <- as.data.frame(s$conf.int); co <- as.data.frame(s$coefficients)
    e4_rows[[paste(ep, m)]] <- tibble(endpoint = ep, model = m, n = s$n, events = s$nevent,
                                      term = rownames(ci), HR = ci[, 1], lo = ci[, 3], hi = ci[, 4],
                                      p = co[, "Pr(>|z|)"])
  }
}
e4 <- bind_rows(e4_rows)
write_csv(e4, file.path(OUT3, "E4_cox_multivariable.csv"))

full <- e4 %>% filter(model == "MDK + all clinical") %>%
  mutate(term = factor(term, rev(unique(term))),
         lab = sprintf("%.2f (%.2f-%.2f)  p=%.2g", HR, lo, hi, p))
pE4 <- ggplot(full, aes(HR, term)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "#3C5488") +
  geom_point(aes(colour = term == "MDK_z"), size = 1.8, show.legend = FALSE) +
  geom_text(aes(x = max(hi) * 1.15, label = lab), hjust = 0, size = 2.1) +
  scale_colour_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "#3C5488")) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.9))) +
  facet_wrap(~endpoint) +
  labs(x = "Hazard ratio (log scale); MDK per 1 SD", y = NULL,
       title = sprintf("GSE14520: multivariable Cox, MDK + %s", paste(covs_full, collapse = " + "))) + theme_fig
ggsave(file.path(OUT3, "E4_cox_multivariable_forest.pdf"), pE4, width = 9, height = 3.8)


pdf(file.path(OUT3, "E5_KM_tertiles.pdf"), width = 8, height = 3.6)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), cex = 0.7)
for (ep in names(ENDP)) {
  dd <- d %>% mutate(time = .data[[ENDP[[ep]][1]]], ev = .data[[ENDP[[ep]][2]]]) %>% filter(!is.na(time), !is.na(ev))
  km <- survfit(Surv(time, ev) ~ MDK_tert, dd)
  p_all <- 1 - pchisq(survdiff(Surv(time, ev) ~ MDK_tert, dd)$chisq, 2)
  dhl <- dd %>% filter(MDK_tert != "Mid") %>% mutate(MDK_tert = droplevels(MDK_tert))
  p_hl <- 1 - pchisq(survdiff(Surv(time, ev) ~ MDK_tert, dhl)$chisq, 1)
  cols <- c("#4DBBD5", "grey55", "#E64B35")
  plot(km, col = cols, lwd = 1.5, xlab = "Months", ylab = paste(ep, "probability"),
       main = sprintf("%s by MDK tertile (3-group p = %.2g; high vs low p = %.2g)", ep, p_all, p_hl))
  legend("bottomleft", sprintf("%s (n=%d)", levels(dd$MDK_tert), as.integer(table(dd$MDK_tert))),
         col = cols, lwd = 1.5, bty = "n")
}
dev.off()

message("\n=== E1: MDK vs clinical (sorted by p) ===");   print(as.data.frame(e1 %>% dplyr::select(variable, groups, p, padj, trend_rho, trend_p)))
message("\n=== E3: cox.zph (p < 0.05 = PH violated) ==="); print(as.data.frame(zph))
message("\n=== E4: MDK term in multivariable models ===")
print(as.data.frame(e4 %>% filter(term == "MDK_z") %>% dplyr::select(endpoint, model, n, events, HR, lo, hi, p)))
message("\nOutputs in ", OUT3)
dd <- d %>% filter(!is.na(RFS_months), !is.na(RFS_event)) %>%
  mutate(t24 = pmin(RFS_months, 24), e24 = ifelse(RFS_months <= 24 & RFS_event == 1, 1, 0))
summary(coxph(Surv(t24, e24) ~ MDK_z + TNM_group, dd))
cox.zph(coxph(Surv(t24, e24) ~ MDK_z + TNM_group, dd))
