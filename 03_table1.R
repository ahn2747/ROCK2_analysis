# =============================================================
# 03_table1.R
# Table 1: 유전자 발현 그룹(Low/High)별 임상 특성 + 전체 코호트 요약
#   - 연속형: median [IQR], Wilcoxon rank-sum
#   - 범주형: chi-square (기대빈도 < 5 셀이 있으면 Fisher's exact로 자동 전환)
#   - 결측은 "Missing" 행으로 표시, p-value 계산에서는 제외
# 입력: data/processed/<CANCER>_merged.rds
# 결과 (파일명 앞에 <CANCER>_): output/tables/<CANCER>/table1_<GENE>.{docx,csv}, table1_overall.{docx,csv}
# =============================================================

source("config.R")
source("R/utils.R")
suppressPackageStartupMessages({
  library(gtsummary)
  library(flextable)
})

# ---- 범주형 검정: chi-square ↔ Fisher 자동 선택 ----
# 결측 제외 후 실제 관측된 수준만으로 분할표를 만들고 기대빈도 확인
# 선택된 검정명은 test_log에 기록 → CSV의 Test 열
test_log <- new.env()

test_chisq_fisher <- function(data, variable, by, ...) {
  d   <- data[!is.na(data[[variable]]) & !is.na(data[[by]]), ]
  tab <- table(droplevels(factor(d[[variable]])), droplevels(factor(d[[by]])))
  res <- if (nrow(tab) < 2 || ncol(tab) < 2) {
    data.frame(p.value = NA_real_, method = "Not tested (single category)")
  } else if (all(suppressWarnings(chisq.test(tab, correct = FALSE)$expected) >= 5)) {
    data.frame(p.value = chisq.test(tab, correct = FALSE)$p.value,
               method = "Pearson's chi-squared test")
  } else {
    # 큰 분할표에서 정확 검정 메모리 초과 시 Monte Carlo로 대체 (seed 고정)
    tryCatch(
      data.frame(p.value = fisher.test(tab, workspace = 2e8)$p.value,
                 method = "Fisher's exact test"),
      error = function(e) {
        set.seed(2026)
        data.frame(p.value = fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value,
                   method = "Fisher's exact test (Monte Carlo, B = 100,000)")
      }
    )
  }
  assign(variable, res$method, envir = test_log)
  res
}

# 표의 각 행에 대응하는 검정명 (변수 라벨 행에만)
test_column <- function(tbl) {
  b <- tbl$table_body
  if (!"p.value" %in% names(b)) return(NULL)
  logged <- unlist(mget(b$variable, envir = test_log, ifnotfound = NA))
  ifelse(b$row_type != "label", "",
         ifelse(b$var_type == "continuous", "Wilcoxon rank-sum test", logged))
}

make_summary <- function(d, vars, by = NULL) {
  labs <- var_labels[intersect(vars, names(var_labels))]
  tbl_summary(
    d, by = all_of(by), include = all_of(vars),
    type      = list(all_dichotomous() ~ "categorical", any_of(cfg$table1_continuous) ~ "continuous"),
    statistic = list(all_continuous() ~ "{median} [{p25}, {p75}]",
                     all_categorical() ~ "{n} ({p}%)"),
    digits    = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    label     = as.list(labs),
    missing = "ifany", missing_text = "Missing"
  ) |>
    bold_labels()
}

save_tbl <- function(tbl, cancer, stem) {
  dir <- out_dir(cancer, "tables")
  # flextable 함수는 네임스페이스로 호출 (ggpubr::font 등이 flextable::font를 가림)
  ft <- gtsummary::as_flex_table(tbl) |>
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = 9, part = "all") |>
    flextable::autofit()
  docx <- file.path(dir, out_file(cancer, paste0(stem, ".docx")))
  flextable::save_as_docx(ft, path = docx)

  # CSV: 마크다운(**, __) 및 헤더 줄바꿈 제거, 검정명 열 추가
  csv_df <- as_tibble(tbl, col_labels = TRUE)
  names(csv_df) <- str_squish(gsub("\\*\\*|__", "", gsub("\\s*\n", ", ", names(csv_df))))
  csv_df[[1]] <- gsub("__", "", csv_df[[1]])
  tests <- test_column(tbl)
  if (!is.null(tests)) csv_df$Test <- tests
  csv <- save_table(csv_df, cancer, paste0(stem, ".csv"))
  c(docx, csv)
}

created <- character()

for (cancer in cfg$cancers) {
  cat("\n==========", cancer, "==========\n")
  df    <- readRDS(processed_path(cancer, "merged.rds"))
  vars  <- intersect(cfg$table1_vars, names(df))
  genes <- detect_genes(df)
  cat("변수", length(vars), "개:", paste(vars, collapse = ", "), "\n")

  # ---- 전체 코호트 ----
  tbl <- make_summary(df, vars) |>
    modify_caption(paste0("Table 1. Clinical characteristics of the TCGA-", cancer,
                          " cohort (N = ", nrow(df), ")"))
  created <- c(created, save_tbl(tbl, cancer, "table1_overall"))

  # ---- 유전자별 ----
  for (g in genes) {
    gene <- toupper(g)
    grp  <- paste0(g, "_group")
    d    <- df[!is.na(df[[grp]]), ]
    cat(sprintf("  %-8s n = %d (%s %d / %s %d)\n", gene, nrow(d),
                group_ref(), sum(d[[grp]] == group_ref()), group_alt(), sum(d[[grp]] == group_alt())))

    rm(list = ls(test_log), envir = test_log)
    tbl <- make_summary(d, vars, by = grp) |>
      add_p(test = list(all_continuous() ~ "wilcox.test",
                        all_categorical() ~ test_chisq_fisher),
            pvalue_fun = label_style_pvalue(digits = 3)) |>
      add_overall(last = FALSE) |>
      modify_spanning_header(c(stat_1, stat_2) ~ paste0("**", gene, " expression**")) |>
      modify_caption(paste0("Table 1. Clinical characteristics of TCGA-", cancer,
                            " patients by ", gene, " expression"))
    created <- c(created, save_tbl(tbl, cancer, paste0("table1_", gene)))

    # ---- 층화 변수(기본 TSS) × 발현 그룹 (배치/기관 효과 점검) ----
    sv  <- cfg$strata_var
    ref <- group_ref()
    alt <- group_alt()
    tab <- table(d[[sv]], d[[grp]])
    test <- test_chisq_fisher(d, sv, grp)
    tss_df <- data.frame(rownames(tab), as.integer(tab[, ref]), as.integer(tab[, alt]),
                         as.integer(rowSums(tab)), sprintf("%.1f", 100 * tab[, alt] / rowSums(tab)))
    names(tss_df) <- c(toupper(sv), ref, alt, "Total", paste0(alt, ", %"))
    tss_df <- tss_df[order(-tss_df$Total), ]
    tss_df$`p-value` <- c(fmt_p(test$p.value), rep("", nrow(tss_df) - 1))
    tss_df$Test      <- c(test$method, rep("", nrow(tss_df) - 1))
    created <- c(created, save_table(tss_df, cancer, paste0(sv, "_by_group_", gene, ".csv")))
    cat(sprintf("           %s × 그룹: %d개 수준, p = %s (%s)\n",
                toupper(sv), nrow(tss_df), fmt_p(test$p.value), test$method))
  }

  # 더 이상 없는 유전자의 이전 결과 → _stale/
  move_stale_outputs(cancer, c("table1_", paste0(cfg$strata_var, "_by_group_")), created)
}

cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")
