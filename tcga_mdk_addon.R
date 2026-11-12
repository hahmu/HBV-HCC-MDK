suppressPackageStartupMessages({ library(survival); library(ggplot2); library(dplyr); library(tidyr) })

WORK_DIR <- "./Bulk_HBV"
IN  <- file.path(WORK_DIR, "TCGA_deconv_LR")
OUT <- file.path(IN, "MDK"); dir.create(OUT, showWarnings = FALSE)
CENSOR_DAYS <- 5 * 365

D  <- readRDS(file.path(IN, "01_tcga_lihc.rds"))
FR <- read.csv(file.path(IN, "02_bayesprism_fractions.csv"), check.names = FALSE); rownames(FR) <- FR$patient
expr <- log2(D$tpm + 1)
pats <- intersect(colnames(expr), FR$patient)

HSC <- c("cyHSC", "Mesenchymal-like HSC", "myHSC"); TAM <- c("SPP1+ TAM", "FOLR2+ TAM", "Kupffer")
AX <- list(   # name, receptor, receiver states (sender = Malignant for all)
  list("MDK-SDC2", "SDC2", c(HSC, "vSMC")),
  list("MDK-LRP1", "LRP1", c("myHSC", TAM)),
  list("MDK-SDC3", "SDC3", TAM))
frac <- function(states) { s <- intersect(states, colnames(FR)); rowSums(FR[pats, s, drop = FALSE]) }

S <- data.frame(patient = pats, MDK = expr["MDK", pats], f_Malignant = frac("Malignant"))
for (a in AX) {
  S[[a[[2]]]] <- expr[a[[2]], pats]
  S[[paste0(a[[1]], "|LR")]] <- sqrt(pmax(S$MDK, 0) * pmax(S[[a[[2]]]], 0))   # same LR score as the main script
  S[[paste0("f_", a[[1]])]] <- frac(a[[3]])
}
S <- merge(S, D$clin, by = "patient")
write.csv(S, file.path(OUT, "MDK_scores_per_tumour.csv"), row.names = FALSE)

## ================================================================== A  boxplots by virus
feat <- c("MDK-SDC2|LR", "MDK-LRP1|LR", "MDK-SDC3|LR", "MDK", "SDC2", "LRP1", "SDC3")
core <- S %>% filter(etiology %in% c("HBV", "HCV", "NBNC"))       # HBV_HCV co-infected left out
long <- core %>% select(patient, etiology, all_of(feat)) %>%
  pivot_longer(all_of(feat), names_to = "feature", values_to = "value") %>%
  mutate(feature = factor(sub("\\|LR", " (LR score)", feature), levels = sub("\\|LR", " (LR score)", feat)),
         etiology = factor(etiology, levels = c("HBV", "HCV", "NBNC")))
tests <- long %>% group_by(feature) %>% summarise(
  n_HBV = sum(etiology == "HBV"), n_other = sum(etiology != "HBV"),
  median_HBV = median(value[etiology == "HBV"]), median_other = median(value[etiology != "HBV"]),
  p_HBV_vs_other = wilcox.test(value[etiology == "HBV"], value[etiology != "HBV"])$p.value,
  p_kruskal = kruskal.test(value, etiology)$p.value, .groups = "drop")
write.csv(tests, file.path(OUT, "A_by_virus_tests.csv"), row.names = FALSE)
print(tests, width = 200)

lab <- tests %>% mutate(txt = sprintf("HBV vs others p = %.2g\nKruskal p = %.2g", p_HBV_vs_other, p_kruskal))
g <- ggplot(long, aes(etiology, value, fill = etiology)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.5) +
  geom_text(data = lab, aes(x = 2, y = Inf, label = txt), inherit.aes = FALSE, vjust = 1.2, size = 2.6) +
  facet_wrap(~feature, scales = "free_y", nrow = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  scale_fill_manual(values = c(HBV = "#C0392B", HCV = "#8E44AD", NBNC = "#2E86C1")) +
  labs(x = NULL, y = "log2(TPM + 1) / LR score", title = "MDK and its receptors in TCGA-LIHC by virus status",
       subtitle = sprintf("core set with virus calls: HBV n = %d, HCV n = %d, NBNC n = %d",
                          sum(core$etiology == "HBV"), sum(core$etiology == "HCV"), sum(core$etiology == "NBNC"))) +
  theme_classic(base_size = 9) + theme(legend.position = "none")
ggsave(file.path(OUT, "A_MDK_by_virus.pdf"), g, width = 10, height = 5.5)
ggsave(file.path(OUT, "A_MDK_by_virus.png"), g, width = 10, height = 5.5, dpi = 200)

## ================================================================== B  Cox
censor <- function(t, e) list(t = pmin(t, CENSOR_DAYS), e = ifelse(t > CENSOR_DAYS, 0L, e))
endpoints <- list(OS = c("os_time", "os_event"))
if ("pfi_time" %in% names(S) && any(!is.na(S$pfi_time))) endpoints$PFI <- c("pfi_time", "pfi_event")
cohorts <- list(All_TCGA = rep(TRUE, nrow(S)), HBV = S$etiology %in% "HBV", NBNC = S$etiology %in% "NBNC")
## variable -> fraction covariates for the "+ cell fractions" model
vars <- list(MDK = "f_Malignant",
             SDC2 = "f_MDK-SDC2", LRP1 = "f_MDK-LRP1", SDC3 = "f_MDK-SDC3",
             `MDK-SDC2|LR` = c("f_Malignant", "f_MDK-SDC2"),
             `MDK-LRP1|LR` = c("f_Malignant", "f_MDK-LRP1"),
             `MDK-SDC3|LR` = c("f_Malignant", "f_MDK-SDC3"))
cox <- function(d, v, covs) {
  d$x <- as.numeric(scale(d[[v]]))
  for (c in covs) d[[make.names(c)]] <- as.numeric(scale(d[[c]]))
  rhs <- paste(c("x", make.names(covs)), collapse = " + ")
  dd <- d[complete.cases(d[, unique(c("t", "e", "x", make.names(covs)))]), ]
  if (sum(dd$e) < 5) return(NULL)
  fit <- tryCatch(coxph(as.formula(paste("Surv(t, e) ~", rhs)), data = dd), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(n = nrow(dd), events = sum(dd$e), HR = s$conf.int["x", 1], lo = s$conf.int["x", 3],
             hi = s$conf.int["x", 4], p = s$coefficients["x", 5])
}
res <- list()
for (ep in names(endpoints)) for (co in names(cohorts)) {
  d0 <- S[cohorts[[co]] & !is.na(S[[endpoints[[ep]][1]]]) & S[[endpoints[[ep]][1]]] > 0, ]
  cs <- censor(d0[[endpoints[[ep]][1]]], d0[[endpoints[[ep]][2]]]); d0$t <- cs$t; d0$e <- cs$e
  for (v in names(vars)) {
    models <- list(univariate = character(0), `+ age + stage` = c("age", "stage"),
                   `+ age + stage + cell fractions` = c("age", "stage", vars[[v]]))
    for (m in names(models)) {
      covs <- models[[m]]
      dd <- d0
      if ("stage" %in% covs) dd$stage <- as.numeric(dd$stage == "III-IV")
      r <- cox(dd, v, covs)
      if (!is.null(r)) res[[length(res) + 1]] <- data.frame(endpoint = ep, cohort = co, variable = v, model = m, r)
    }
  }
}
res <- bind_rows(res)
write.csv(res, file.path(OUT, "B_cox_MDK.csv"), row.names = FALSE)
cat("\nAll TCGA, OS:\n"); print(res %>% filter(endpoint == "OS", cohort == "All_TCGA") %>%
                                  mutate(across(c(HR, lo, hi), ~ round(.x, 2)), p = signif(p, 2)), row.names = FALSE)

pdf(file.path(OUT, "B_forest_MDK.pdf"), width = 9, height = 5.5)
for (ep in names(endpoints)) {
  d <- res %>% filter(endpoint == ep) %>%
    mutate(variable = factor(sub("\\|LR", " (LR score)", variable),
                             levels = rev(sub("\\|LR", " (LR score)", names(vars)))),
           star = ifelse(p < 0.05, "*", ""),
           model = factor(model, levels = c("univariate", "+ age + stage", "+ age + stage + cell fractions")))
  p <- ggplot(d, aes(HR, variable, colour = model)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.25) +
    geom_text(aes(x = hi, label = star), position = position_dodge(width = 0.6), hjust = -0.4, size = 3.5,
              show.legend = FALSE) +
    facet_wrap(~cohort, nrow = 1) + scale_x_log10() +
    scale_colour_manual(values = c("grey60", "#2E86C1", "#C0392B")) +
    labs(x = "Hazard ratio per 1 SD (log scale)", y = NULL, colour = NULL,
         title = sprintf("%s: MDK, its receptors and MDK-receptor LR scores (TCGA-LIHC)", ep),
         subtitle = "cell fractions = BayesPrism malignant (MDK) and receiver fractions; * p < 0.05") +
    theme_classic(base_size = 9) + theme(legend.position = "bottom")
  print(p)
}
dev.off()

## ================================================================== C  KM, LR expression score
km_df <- function(fit) {
  s <- summary(fit, censored = TRUE); st <- sub(".*=", "", as.character(s$strata))
  rbind(data.frame(time = 0, surv = 1, strata = unique(st)), data.frame(time = s$time, surv = s$surv, strata = st))
}
pdf(file.path(OUT, "C_KM_MDK_LR_expression.pdf"), width = 4.2, height = 3.8)
for (ep in names(endpoints)) for (co in names(cohorts)) for (v in c("MDK-SDC2|LR", "MDK-LRP1|LR", "MDK-SDC3|LR", "MDK")) {
  d <- S[cohorts[[co]] & !is.na(S[[endpoints[[ep]][1]]]) & S[[endpoints[[ep]][1]]] > 0, ]
  cs <- censor(d[[endpoints[[ep]][1]]], d[[endpoints[[ep]][2]]]); d$t <- cs$t / 30.44; d$e <- cs$e
  d$grp <- factor(ifelse(d[[v]] > median(d[[v]], na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
  if (sum(d$e) < 5) next
  pv <- 1 - pchisq(survdiff(Surv(t, e) ~ grp, data = d)$chisq, 1)
  kd <- km_df(survfit(Surv(t, e) ~ grp, data = d)); kd$strata <- factor(kd$strata, levels = c("Low", "High"))
  nn <- table(d$grp)
  g <- ggplot(kd, aes(time, surv, colour = strata)) + geom_step(linewidth = 0.7) +
    scale_colour_manual(values = c(Low = "#2E86C1", High = "#C0392B"),
                        labels = c(Low = sprintf("Low (n=%d)", nn[["Low"]]), High = sprintf("High (n=%d)", nn[["High"]]))) +
    scale_y_continuous(limits = c(0, 1)) +
    annotate("text", x = 0, y = 0.05, hjust = 0, size = 3.2, label = sprintf("log-rank p = %.3g", pv)) +
    labs(x = "Months", y = paste(ep, "probability"), colour = NULL, title = sub("\\|LR", " (LR expression score)", v),
         subtitle = sprintf("%s, median split", co)) +
    theme_classic(base_size = 9) + theme(legend.position = c(0.75, 0.9))
  print(g)
}
dev.off()
cat("\nDone. Outputs in", OUT, "\n")
