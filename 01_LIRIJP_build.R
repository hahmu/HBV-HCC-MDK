suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(tidyr); library(stringr)
})

BASE <- path.expand("~/Desktop/LIRI-JP")
RAW  <- file.path(BASE, "raw"); HDR <- file.path(BASE, "headers")
OUT  <- file.path(BASE, "analysis"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

header_of <- function(type) strsplit(readLines(gzfile(file.path(HDR, paste0(type, ".tsv.gz"))), n = 1), "\t")[[1]]

read_type <- function(type, select = NULL) {
  files <- list.files(RAW, pattern = "\\.gz$", recursive = TRUE, full.names = TRUE)
  files <- files[basename(dirname(files)) == type]
  h <- header_of(type)
  message(type, ": ", length(files), " files")
  sel <- if (is.null(select)) seq_along(h) else match(select, h)
  rbindlist(lapply(files, function(f) {
    x <- fread(f, header = FALSE, sep = "\t", quote = "", fill = TRUE, colClasses = "character",
               select = sel, na.strings = c("", "NA"))
    setnames(x, h[sel]); x
  }), fill = TRUE)
}

donor    <- unique(read_type("donor"))
specimen <- unique(read_type("specimen"))
sample   <- unique(read_type("sample"))
expr     <- read_type("exp_seq", c("icgc_donor_id", "icgc_specimen_id", "icgc_sample_id", "submitted_sample_id",
                                   "gene_id", "normalized_read_count", "raw_read_count"))
expr[, `:=`(fpkm = as.numeric(normalized_read_count), count = as.numeric(raw_read_count))]


fpkm  <- dcast(expr, gene_id ~ icgc_sample_id, value.var = "fpkm",  fun.aggregate = sum, fill = 0)
count <- dcast(expr, gene_id ~ icgc_sample_id, value.var = "count", fun.aggregate = sum, fill = 0)
fm <- as.matrix(fpkm[, -1]);  rownames(fm) <- fpkm$gene_id
cm <- as.matrix(count[, -1]); rownames(cm) <- count$gene_id
saveRDS(fm, file.path(OUT, "LIRIJP_fpkm.rds")); saveRDS(cm, file.path(OUT, "LIRIJP_counts.rds"))
message("Expression: ", nrow(fm), " genes x ", ncol(fm), " samples")


rna <- unique(expr[, .(icgc_donor_id, icgc_specimen_id, icgc_sample_id, submitted_sample_id)])
lib <- expr[, .(total_raw_count = sum(count, na.rm = TRUE), n_genes_detected = sum(fpkm > 0, na.rm = TRUE)),
            by = icgc_sample_id]

sp <- specimen %>% as_tibble() %>%
  dplyr::select(icgc_specimen_id, submitted_specimen_id, specimen_type, specimen_donor_treatment_type,
                specimen_donor_treatment_type_other, tumour_confirmed, tumour_histological_type,
                tumour_grading_system, tumour_grade, tumour_stage_system, tumour_stage, percentage_cellularity) %>%
  distinct(icgc_specimen_id, .keep_all = TRUE)

dn <- donor %>% as_tibble() %>%
  transmute(icgc_donor_id, submitted_donor_id, sex = donor_sex,
            age = as.numeric(donor_age_at_diagnosis),
            vital_status = donor_vital_status, disease_status_last_followup,
            relapse_type = donor_relapse_type,
            relapse_interval_days = as.numeric(donor_relapse_interval),
            survival_time_days = as.numeric(donor_survival_time),
            followup_interval_days = as.numeric(donor_interval_of_last_followup),
            stage_system = donor_tumour_staging_system_at_diagnosis,
            stage = donor_tumour_stage_at_diagnosis, prior_malignancy) %>%
  distinct(icgc_donor_id, .keep_all = TRUE) %>%
  mutate(
    OS_event  = case_when(vital_status == "deceased" ~ 1L, vital_status == "alive" ~ 0L, TRUE ~ NA_integer_),
    OS_months = coalesce(survival_time_days, followup_interval_days) / 30.44,
    stage_group = case_when(stage %in% c("1", "I") ~ "I", stage %in% c("2", "II") ~ "II",
                            stage %in% c("3", "III") ~ "III", stage %in% c("4", "IV", "4A", "4B") ~ "IV",
                            TRUE ~ NA_character_))

info <- as_tibble(rna) %>%
  left_join(sp, by = "icgc_specimen_id") %>%
  left_join(dn, by = "icgc_donor_id") %>%
  left_join(as_tibble(lib), by = "icgc_sample_id") %>%
  mutate(tissue = case_when(
    grepl("^Primary tumour", specimen_type)          ~ "Tumour",
    grepl("^Normal", specimen_type)                  ~ "Non-tumour",
    grepl("Recurrent", specimen_type)                ~ "Recurrent tumour",
    grepl("Metasta", specimen_type)                  ~ "Metastasis",
    TRUE ~ paste0("Other: ", specimen_type)))


vfile <- file.path(BASE, "virus_status.csv")
if (file.exists(vfile)) {
  # virus_status.csv built from Fujimoto 2016 Supp. Table 1: one row per patient, with virus, histology
  # (HCC / ICC / cHCC/CC) and the paper's clinical columns. Extra columns get the prefix "paper_".
  v <- read_csv(vfile, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
    mutate(submitted_donor_id = str_trim(submitted_donor_id), virus = str_trim(virus)) %>%
    distinct(submitted_donor_id, .keep_all = TRUE)
  keep <- setdiff(names(v), c("submitted_donor_id", "virus", "histology"))
  v <- v %>% rename_with(~ paste0("paper_", .x), all_of(keep))
  if (!"histology" %in% names(v)) v$histology <- NA_character_
  info <- info %>% left_join(v, by = "submitted_donor_id")
  unmatched <- setdiff(unique(info$submitted_donor_id), v$submitted_donor_id)
  if (length(unmatched)) message("RNA-seq donors without a row in virus_status.csv: ", paste(unmatched, collapse = ", "))
  message("Virus status merged for ", sum(!is.na(info$virus)), " of ", nrow(info), " samples")
} else {
  info$virus <- NA_character_; info$histology <- NA_character_
  tmpl <- info %>% distinct(submitted_donor_id, icgc_donor_id) %>% arrange(submitted_donor_id) %>% mutate(virus = "")
  write_csv(tmpl, file.path(BASE, "virus_status_TEMPLATE.csv"))
  message("No virus_status.csv yet -> wrote virus_status_TEMPLATE.csv (", nrow(tmpl), " donors). ",
          "Fill 'virus' (HBV / HCV / HBV+HCV / NBNC), save as virus_status.csv, rerun.")
}
info <- info %>% mutate(virus_group = case_when(
  is.na(virus) | virus == "" ~ "Unknown",
  grepl("HBV", virus, ignore.case = TRUE) & grepl("HCV", virus, ignore.case = TRUE) ~ "HBV+HCV",
  grepl("HBV", virus, ignore.case = TRUE) ~ "HBV",
  grepl("HCV", virus, ignore.case = TRUE) ~ "HCV",
  grepl("NBNC|non|none|neg", virus, ignore.case = TRUE) ~ "NBNC",
  TRUE ~ paste0("Other: ", virus)))


info <- info %>% group_by(icgc_donor_id) %>%
  mutate(paired_TN = any(tissue == "Tumour") & any(tissue == "Non-tumour")) %>% ungroup() %>%
  dplyr::select(icgc_sample_id, submitted_sample_id, icgc_donor_id, submitted_donor_id, tissue, specimen_type,
                paired_TN, virus, virus_group, histology, sex, age, stage, stage_group, tumour_grade,
                OS_event, OS_months, vital_status, relapse_type, relapse_interval_days,
                everything()) %>%
  arrange(submitted_donor_id, tissue)
write_csv(info, file.path(OUT, "LIRIJP_sample_info.csv"))


pt <- info %>% filter(tissue %in% c("Tumour", "Non-tumour")) %>%
  group_by(icgc_donor_id, submitted_donor_id, tissue) %>%
  summarise(samples = paste(icgc_sample_id, collapse = ";"), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = samples, names_prefix = "RNA_") %>%
  left_join(info %>% distinct(icgc_donor_id, paired_TN, virus_group, histology, sex, age, stage_group,
                              OS_event, OS_months, relapse_type), by = "icgc_donor_id")
write_csv(pt, file.path(OUT, "LIRIJP_patient_info.csv"))


cnt <- function(df, v) df %>% count(level = as.character(.data[[v]])) %>%
  transmute(section = v, level = ifelse(is.na(level), "missing", level), n)
tum <- info %>% filter(tissue == "Tumour") %>% distinct(icgc_donor_id, .keep_all = TRUE)
summ <- bind_rows(
  info %>% count(level = tissue) %>% transmute(section = "RNA-seq samples by tissue", level, n),
  tibble(section = "Donors", level = c("with RNA-seq", "with tumour RNA", "paired tumour + non-tumour RNA"),
         n = c(n_distinct(info$icgc_donor_id), n_distinct(tum$icgc_donor_id),
               n_distinct(info$icgc_donor_id[info$paired_TN]))),
  info %>% filter(tissue == "Tumour") %>% count(icgc_donor_id) %>% count(level = paste(n, "tumour RNA samples")) %>%
    transmute(section = "Tumour RNA samples per donor", level, n),
  cnt(tum, "virus_group") %>% mutate(section = "Virus (tumour donors)"),
  cnt(tum %>% filter(paired_TN), "virus_group") %>% mutate(section = "Virus (paired donors)"),
  tum %>% count(level = paste(coalesce(histology, "missing"), virus_group, sep = " | ")) %>%
    transmute(section = "Histology | virus (tumour donors)", level, n),
  tum %>% filter(paired_TN) %>% count(level = paste(coalesce(histology, "missing"), virus_group, sep = " | ")) %>%
    transmute(section = "Histology | virus (paired donors)", level, n),
  cnt(tum, "sex"), cnt(tum, "stage"), cnt(tum, "tumour_grade"),
  cnt(tum, "vital_status"), cnt(tum, "relapse_type"), cnt(tum, "specimen_donor_treatment_type"),
  tibble(section = "OS_months (tumour donors)", level = c("median", "missing"),
         n = c(round(median(tum$OS_months, na.rm = TRUE), 1), sum(is.na(tum$OS_months)))))
write_csv(summ, file.path(OUT, "LIRIJP_summary_counts.csv"))
print(as.data.frame(summ), right = FALSE)
message("\nWritten to ", OUT)
