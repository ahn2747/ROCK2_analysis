# =============================================================
# 00_setup.R
# 필요한 패키지 설치 (최초 1회)
#   Rscript 00_setup.R          # CRAN + Bioconductor 전체
#   Rscript 00_setup.R --core   # 01~04에 필요한 CRAN 패키지만
# =============================================================

core_only <- "--core" %in% commandArgs(trailingOnly = TRUE)

cran_pkgs <- c(
  "haven", "dplyr", "tidyr", "readr", "stringr", "purrr", "janitor",   # 01, 02b
  "gtsummary", "flextable", "officer",                                 # 03 Table 1
  "survival", "survminer", "ggplot2", "ggpubr", "forestmodel",         # 04 생존분석
  "survRM2", "jsonlite", "ggrepel",                                    # 04 RMST, 01 manifest, 05 volcano
  "msigdbr", "BiocManager"                                             # 05 GSEA
)
bioc_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment",                              # 02a
  "DESeq2", "apeglm", "BiocParallel", "limma", "fgsea",               # 05
  "clusterProfiler", "enrichplot",                                     # 05
  "GEOquery", "org.Hs.eg.db"                                           # 06
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)

if (!core_only) {
  missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
  if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

pkgs <- if (core_only) cran_pkgs else c(cran_pkgs, bioc_pkgs)
still_missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(still_missing)) {
  stop("설치 실패: ", paste(still_missing, collapse = ", "))
}
message("패키지 준비 완료 (", length(pkgs), "개)")
