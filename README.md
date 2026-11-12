# HBV-HCC-MDK
# HBV-HCC MDK analysis scripts

Code for: "Malignant-cell MDK signaling at the tumor–stroma interface and its association
with poor survival in hepatitis B virus-related hepatocellular carcinoma".

## Run order
1. scRNA-seq (GSE156625, GSE149614): Integration_script.R → Annotation.R → Subclustering.R
   → run_infercnv_whole.py → DESeq2.R → run_sccoda.py
2. Cell–cell communication (LIANA): run_liana_normal_fine.py → run_liana_de.py → run_liana_figs.py
3. Spatial (GSE281759 Visium, cell2location): s00 → s01 → s02 → s03 → s04 → s05 → s06 → s07 → s08
4. Bulk validation: tcga_lihc_deconv_LR_survival.R → tcga_mdk_addon.R (TCGA-LIHC);
   01–03_GSE14520_*.R (GSE14520); 01–02_LIRIJP_*.R (ICGC LIRI-JP)

## Data
All data are public: GEO GSE156625, GSE149614, GSE281759, GSE14520; TCGA-LIHC (GDC); ICGC LIRI-JP.
File paths inside the scripts point to local folders and need to be edited before running.
EOF
