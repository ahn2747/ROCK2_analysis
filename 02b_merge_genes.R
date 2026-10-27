# =============================================================
# 02b_merge_genes.R
# database/gene_files/<CANCER>_<split>_<GENE>.csv 를 임상 데이터에 병합
#   - 병합 열: Patient(→ sample_id), Expression, Group 만
#   - 이미 .sav에 있는 유전자 → 병합 없이 QC만 (CSV 값 == .sav 값 확인)
#   - 모든 CSV: Status(Alive/Dead) ↔ .sav status(0/1) 교차 확인, 불일치 시 중단
#   - cfg$genes_from_gdc: 02a 발현 행렬에서 추출 + median split
# 결과:
#   data/processed/<CANCER>_merged.{rds,csv,sav}
#   output/tables/<CANCER>/<CANCER>_merge_QC.csv, _merge_QC_status_crosstab.csv, _merge_QC_no_survival.csv
# =============================================================

source("config.R")
source("R/utils.R")

gene_files <- list.files(cfg$gene_dir, pattern = cfg$gene_pattern, full.names = TRUE)

for (cancer in cfg$cancers) {
  cat("\n==========", cancer, "==========\n")
  # 전체 환자 기준으로 병합/QC (그룹은 전체 코호트 기준으로 정의됨) → 마지막에 제외 적용
  clin <- readRDS(processed_path(cancer, "clinical_all.rds"))
  existing <- detect_genes(clin)
  surv <- clin[, c("sample_id", "surv_time", "status")]
  names(surv) <- c("sample_id", "sav_time", "sav_status")

  files <- Filter(function(f) toupper(str_match(basename(f), cfg$gene_pattern)[2]) == cancer,
                  gene_files)
  cat("유전자 파일", length(files), "개:", paste(basename(files), collapse = ", "), "\n")

  qc_rows <- list(); xtabs <- list()

  for (f in files) {
    gf <- read_gene_file(f)
    d  <- gf$data
    m  <- inner_join(d, surv, by = "sample_id")
    not_in_clin <- setdiff(d$sample_id, clin$sample_id)

    # ---- Status 교차 확인 (Dead = 1, Alive = 0) ----
    if (any(!is.na(d$csv_status))) {
      tab <- as.data.frame(table(csv_status = factor(m$csv_status, 0:1, c("Alive", "Dead")),
                                 sav_status = m$sav_status, useNA = "ifany"))
      xtabs[[gf$file]] <- data.frame(file = gf$file, tab)
      bad <- m$sample_id[!is.na(m$csv_status) & !is.na(m$sav_status) & m$csv_status != m$sav_status]
      if (length(bad)) {
        stop(gf$file, ": CSV Status와 .sav status 불일치 (Dead=1, Alive=0 아님) → ",
             paste(bad, collapse = ", "))
      }
    }

    row <- data.frame(
      gene = gf$gene, file = gf$file,
      action = if (gf$key %in% existing) "QC only (already in .sav)" else "merged",
      n_file = nrow(d), n_matched = nrow(m),
      n_file_not_in_clinical = length(not_in_clin),
      file_not_in_clinical = paste(not_in_clin, collapse = " "),
      status_mismatch = 0
    )

    if (gf$key %in% existing) {
      # ---- 기존 유전자: .sav 값 유지, 일치 여부만 확인 ----
      sav_expr  <- clin[[paste0(gf$key, "_expression")]][match(m$sample_id, clin$sample_id)]
      sav_group <- clin[[paste0(gf$key, "_group")]][match(m$sample_id, clin$sample_id)]
      both <- !is.na(sav_expr) & !is.na(m$expression)
      row$n_compared          <- sum(both)
      row$expr_identical_pct  <- round(100 * mean(abs(sav_expr[both] - m$expression[both]) < 0.005), 2)
      row$expr_correlation    <- round(cor(sav_expr[both], m$expression[both]), 6)
      row$group_identical_pct <- round(100 * mean(as.character(sav_group[both]) == as.character(m$group[both])), 2)
      if (row$group_identical_pct < 100) {
        warning(cancer, " ", gf$gene, ": CSV Group과 .sav 그룹 일치율 ", row$group_identical_pct, "%")
      }
    } else {
      # ---- 새 유전자: 발현값 + 그룹(파일 값 그대로) 병합 ----
      add <- d[, c("sample_id", "expression", "group")]
      names(add) <- c("sample_id", paste0(gf$key, c("_expression", "_group")))
      clin <- left_join(clin, add, by = "sample_id")
      existing <- c(existing, gf$key)
    }
    # 파일 자체 행 기준 median split QC (허용 오차 cfg$median_tie_tolerance)
    row <- cbind(row, check_median_split(d$expression, d$group))
    qc_rows[[gf$file]] <- row
  }

  # ---- 02a 발현 행렬에서 추가 유전자 ----
  gdc_new <- setdiff(gene_key(cfg$genes_from_gdc), existing)
  if (length(gdc_new)) {
    want <- cfg$genes_from_gdc[gene_key(cfg$genes_from_gdc) %in% gdc_new]
    add  <- extract_gdc_genes(cancer, want, clin$sample_id)
    clin <- left_join(clin, add, by = "sample_id")
    for (g in setdiff(detect_genes(add), existing)) {
      x <- clin[[paste0(g, "_expression")]]
      qc_rows[[g]] <- data.frame(gene = toupper(g), file = "GDC matrix (02a)", action = "merged (median split)",
                                 n_file = sum(!is.na(x)), n_matched = sum(!is.na(x)),
                                 check_median_split(x, clin[[paste0(g, "_group")]]))
    }
  }

  # ---- 발현값은 있으나 생존 결측인 환자 (생존분석에서 자동 제외) ----
  genes <- detect_genes(clin)
  has_expr <- rowSums(!is.na(clin[paste0(genes, "_expression")])) > 0
  no_surv  <- has_expr & (is.na(clin$surv_time) | is.na(clin$status))
  no_surv_df <- data.frame(sample_id = clin$sample_id[no_surv],
                           n_genes_with_expression = rowSums(!is.na(clin[no_surv, paste0(genes, "_expression"), drop = FALSE])))
  cat("발현 있음 + 생존 결측:", if (nrow(no_surv_df)) paste(no_surv_df$sample_id, collapse = ", ") else "없음", "\n")

  # ---- log2 복사본, 라벨 ----
  for (g in genes) {
    clin[[paste0(g, "_expression_log2")]] <- log2(clin[[paste0(g, "_expression")]] + 1)
    attr(clin[[paste0(g, "_group")]], "label") <- paste(toupper(g), "expression")
  }

  # ---- 신보조요법 제외 ----
  n_before <- nrow(clin)
  if (cfg$exclude_neoadjuvant) clin <- clin[!clin$neoadjuvant_flag, ]
  cat("환자:", n_before, "→", nrow(clin), "\n")

  # ---- QC 출력 ----
  qc <- bind_rows(qc_rows)
  if (nrow(qc)) print(qc[, intersect(c("gene", "action", "n_file", "n_matched", "n_file_not_in_clinical",
                                       "expr_identical_pct", "group_identical_pct",
                                       "n_diff_from_median_split", "median_split_pass"), names(qc))],
                      row.names = FALSE)
  if (length(xtabs)) {
    xt <- bind_rows(xtabs)
    cat("Status 교차표 (CSV × .sav):\n"); print(xt, row.names = FALSE)
    save_table(xt, cancer, "merge_QC_status_crosstab.csv")
  }
  save_table(qc, cancer, "merge_QC.csv")
  save_table(no_surv_df, cancer, "merge_QC_no_survival.csv")

  # ---- 저장 (database/*.sav는 덮어쓰지 않음) ----
  saveRDS(clin, processed_path(cancer, "merged.rds"))
  write.csv(clin, processed_path(cancer, "merged.csv"), row.names = FALSE, na = "")
  haven::write_sav(clin, processed_path(cancer, "merged.sav"))
  cat("저장:", processed_path(cancer, "merged.{rds,csv,sav}"), "—", nrow(clin), "명,", length(genes), "개 유전자\n")
}
