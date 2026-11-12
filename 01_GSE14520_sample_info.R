suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr); library(affyio)
})

BASE  <- path.expand("~/GSE14520")
RAW   <- file.path(BASE, "GSE14520_RAW")
EXTRA <- file.path(BASE, "GSE14520_Extra_Supplement.txt")
OUT   <- file.path(BASE, "analysis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)


for (t in list.files(RAW, pattern = "\\.tar$", full.names = TRUE)) {
  message("Untarring ", basename(t)); untar(t, exdir = RAW)
}
cel <- list.files(RAW, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
if (!length(cel)) stop("No CEL files in ", RAW)
cel_df <- tibble(cel_file = cel, gsm = str_extract(basename(cel), "GSM[0-9]+")) %>%
  distinct(gsm, .keep_all = TRUE)
cel_df$cdf <- vapply(cel_df$cel_file, function(f)
  tryCatch(read.celfile.header(f)$cdfName, error = function(e) NA_character_), character(1))
cel_df$platform <- case_when(cel_df$cdf == "HT_HG-U133A" ~ "GPL3921_HT-HG-U133A",
                             cel_df$cdf == "HG-U133A_2"  ~ "GPL571_HG-U133A-2",
                             TRUE ~ paste0("other_", cel_df$cdf))
message(nrow(cel_df), " CEL files"); print(table(cel_df$platform, useNA = "ifany"))


ex <- read.delim(EXTRA, check.names = FALSE, stringsAsFactors = FALSE,
                 na.strings = c("", ".", "NA"), strip.white = TRUE)
names(ex) <- str_trim(names(ex))
ex <- ex %>% filter(!is.na(Affy_GSM))

ex <- ex %>%
  transmute(
    gsm          = Affy_GSM,
    lcs_id       = `LCS ID`,
    patient_id   = str_remove(`LCS ID`, "[AB]$"),
    hospital_id  = str_trim(ID),
    tissue       = recode(`Tissue Type`, "Tumor" = "Tumour", "Non-Tumor" = "Non-tumour",
                          "Normal" = "Healthy liver pool"),
    sex          = Gender,
    age          = as.numeric(Age),
    HBV_status   = `HBV viral status`,
    ALT          = `ALT(>/<=50U/L)`,
    tumour_size  = `Main Tumor Size (>/<=5 cm)`,
    multinodular = Multinodular,
    cirrhosis    = Cirrhosis,
    TNM          = `TNM staging`,
    BCLC         = `BCLC staging`,
    CLIP         = `CLIP staging`,
    AFP          = `AFP (>/<=300ng/ml)`,
    metastasis_signature = `Predicted risk Metastasis Signature`,
    OS_event     = as.numeric(`Survival status`),   # 1 = death
    OS_months    = as.numeric(`Survival months`),
    RFS_event    = as.numeric(`Recurr status`),     # 1 = recurrence
    RFS_months   = as.numeric(`Recurr months`),
    agilent_cgh_gsm = Agilent_GSM)


clin_vars <- c("sex", "age", "HBV_status", "ALT", "tumour_size", "multinodular", "cirrhosis",
               "TNM", "BCLC", "CLIP", "AFP", "metastasis_signature",
               "OS_event", "OS_months", "RFS_event", "RFS_months")
tum_clin <- ex %>% filter(tissue == "Tumour") %>% select(patient_id, all_of(clin_vars))
ex <- ex %>% select(-all_of(clin_vars)) %>% left_join(tum_clin, by = "patient_id")

ex <- ex %>% mutate(
  HBV_label = case_when(
    HBV_status == "AVR-CC" ~ "AVR-CC (HBV carrier, active viral replication)",
    HBV_status == "CC"     ~ "CC (HBV chronic carrier)",
    HBV_status == "N"      ~ "N (code not defined in GEO; likely no HBV carrier status)",
    tissue == "Healthy liver pool" ~ "Healthy liver pool",
    TRUE                   ~ "Missing"),
  TNM_group = case_when(TNM == "I" ~ "I", TNM == "II" ~ "II", grepl("^III", TNM) ~ "III", TRUE ~ NA_character_))


info <- full_join(cel_df, ex, by = "gsm") %>%
  mutate(has_CEL = !is.na(cel_file), in_supplement = !is.na(tissue))

info <- info %>% group_by(patient_id) %>%
  mutate(paired_TN = has_CEL & tissue %in% c("Tumour", "Non-tumour") &
           sum(tissue == "Tumour" & has_CEL) == 1 & sum(tissue == "Non-tumour" & has_CEL) == 1 &
           n_distinct(platform[tissue %in% c("Tumour", "Non-tumour")]) == 1) %>%
  ungroup() %>%
  select(gsm, platform, tissue, patient_id, lcs_id, paired_TN, HBV_status, HBV_label, TNM, TNM_group,
         BCLC, everything()) %>%
  arrange(platform, patient_id, desc(tissue))
write_csv(info, file.path(OUT, "GSE14520_sample_info.csv"))


pt <- info %>% filter(tissue %in% c("Tumour", "Non-tumour")) %>%
  select(patient_id, tissue, gsm, platform) %>%
  pivot_wider(names_from = tissue, values_from = c(gsm, platform)) %>%
  left_join(info %>% filter(tissue == "Tumour") %>%
              select(patient_id, hospital_id, paired_TN, HBV_status, HBV_label, TNM, TNM_group, BCLC, CLIP,
                     sex, age, ALT, tumour_size, multinodular, cirrhosis, AFP,
                     OS_event, OS_months, RFS_event, RFS_months),
            by = "patient_id")
write_csv(pt, file.path(OUT, "GSE14520_patient_info.csv"))


counts <- info %>% filter(has_CEL | in_supplement) %>%
  count(platform, tissue, HBV_label, name = "n_samples") %>%
  left_join(info %>% filter(paired_TN, tissue == "Tumour") %>%
              count(platform, HBV_label, name = "n_paired_patients") %>% mutate(tissue = "Tumour"),
            by = c("platform", "tissue", "HBV_label"))
write_csv(counts, file.path(OUT, "GSE14520_summary_counts.csv"))

tum <- info %>% filter(tissue == "Tumour")
cat_vars <- c("platform", "HBV_label", "sex", "ALT", "tumour_size", "multinodular", "cirrhosis",
              "TNM", "TNM_group", "BCLC", "CLIP", "AFP", "OS_event", "RFS_event", "paired_TN")
clin_sum <- bind_rows(
  lapply(cat_vars, function(v) tum %>% count(level = as.character(.data[[v]])) %>%
           transmute(variable = v, level = ifelse(is.na(level), "missing", level), value = as.character(n))),
  lapply(c("age", "OS_months", "RFS_months"), function(v) {
    x <- tum[[v]]
    tibble(variable = v, level = c("median", "IQR", "missing"),
           value = c(as.character(median(x, na.rm = TRUE)),
                     paste(quantile(x, c(.25, .75), na.rm = TRUE), collapse = " - "),
                     as.character(sum(is.na(x)))))
  }))
write_csv(clin_sum, file.path(OUT, "GSE14520_clinical_summary.csv"))

message("\nCEL files without supplement row: ", sum(info$has_CEL & !info$in_supplement),
        " | supplement rows without CEL: ", sum(!info$has_CEL & info$in_supplement))
message("Paired patients (T + NT on same platform): ", sum(pt$paired_TN, na.rm = TRUE))
print(as.data.frame(counts))
message("Written to ", OUT)
