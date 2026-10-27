# =============================================================
# 01_import_clinical.R
# <CANCER>.sav 스키마 검증 → 읽기 → 표준 열 이름 → 정리/재코딩 → QC → manifest
# 결과:
#   data/processed/<CANCER>_clinical_all.rds  (전체 환자, neoadjuvant_flag 포함)
#   data/processed/<CANCER>_clinical.rds      (분석용: cfg$exclude_neoadjuvant 적용)
#   data/processed/db_manifest.json           (.sav md5, n, 생존 n/사건, 열, 유전자, 신보조요법 ID)
#   output/tables/<CANCER>/<CANCER>_clinical_QC_median_split.csv, _clinical_QC_labels.csv
# manifest 규칙:
#   같은 md5인데 n/생존 n/사건 수가 다름 → 중단 (코드 변경이 결과를 바꿈)
#   md5가 바뀜 → 이전 대비 변경 내용 출력 후 계속 (manifest 갱신)
# =============================================================

source("config.R")
source("R/utils.R")

surv_counts <- function(d) {
  ok <- !is.na(d$surv_time) & !is.na(d$status)
  list(n = sum(ok), events = sum(d$status[ok] == 1))
}

# ---- manifest ----------------------------------------------------------

read_manifest <- function() {
  if (!file.exists(cfg$manifest)) return(list())
  jsonlite::read_json(cfg$manifest, simplifyVector = TRUE)
}

manifest_entry <- function(cancer, df) {
  f <- sav_path(cancer)
  list(
    file            = if (file.exists(f)) f else "GDC",
    md5             = if (file.exists(f)) unname(tools::md5sum(f)) else NA_character_,
    updated         = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n               = nrow(df),
    surv_incl       = surv_counts(df),
    surv_excl       = surv_counts(df[!df$neoadjuvant_flag, ]),
    columns         = names(df),
    genes           = toupper(detect_genes(df)),
    neoadjuvant_ids = df$sample_id[df$neoadjuvant_flag]
  )
}

# 이전 항목 대비 변경 내용 (문자열 벡터)
manifest_diff <- function(old, new) {
  out <- character()
  cnt <- function(x) sprintf("%s (%s)", x$n, x$events)
  if (!identical(as.integer(old$n), as.integer(new$n))) out <- c(out, sprintf("환자 수 %s → %s", old$n, new$n))
  for (k in c("surv_incl", "surv_excl")) {
    if (!identical(cnt(old[[k]]), cnt(new[[k]]))) {
      out <- c(out, sprintf("생존 n (사건) [%s] %s → %s", k, cnt(old[[k]]), cnt(new[[k]])))
    }
  }
  for (k in c("columns", "genes", "neoadjuvant_ids")) {
    add <- setdiff(new[[k]], old[[k]])
    rem <- setdiff(old[[k]], new[[k]])
    if (length(add)) out <- c(out, sprintf("%s 추가: %s", k, paste(add, collapse = ", ")))
    if (length(rem)) out <- c(out, sprintf("%s 삭제: %s", k, paste(rem, collapse = ", ")))
  }
  out
}

counts_changed <- function(old, new) {
  !identical(as.integer(old$n), as.integer(new$n)) ||
    !identical(unlist(lapply(old[c("surv_incl", "surv_excl")], unlist)) |> as.integer(),
               unlist(lapply(new[c("surv_incl", "surv_excl")], unlist)) |> as.integer())
}

# ---- 실행 ---------------------------------------------------------------

manifest <- read_manifest()

for (cancer in cfg$cancers) {
  cat("\n==========", cancer, "==========\n")
  validate_schema(cancer)

  df <- load_clinical(cancer) |> recode_clinical(cancer)
  genes <- detect_genes(df)
  cat("환자", nrow(df), "명 /", ncol(df), "열\n")
  cat("유전자:", paste(toupper(genes), collapse = ", "), "\n")

  # ---- 라벨 일관성 (그룹 변수 vs 원 변수 cutoff) ----
  checks <- attr(df, "label_checks")
  if (!is.null(checks)) {
    cat("라벨 일관성:", paste(sprintf("%s ~ %s > %s: %d명 일치", checks$variable, checks$source,
                                      checks$cutoff, checks$n_compared), collapse = "; "), "\n")
    save_table(checks, cancer, "clinical_QC_labels.csv")
  }

  # ---- 신보조요법 환자 ----
  neo_ids <- df$sample_id[df$neoadjuvant_flag]
  cat("신보조요법(neoadjuvant) 환자", length(neo_ids), "명:",
      if (length(neo_ids)) paste(neo_ids, collapse = ", ") else "-", "\n")

  # ---- 생존 n / 사건 수 ----
  entry <- manifest_entry(cancer, df)
  cat(sprintf("생존분석 n (사건): 제외 전 %d (%d), 신보조요법 제외 후 %d (%d)\n",
              entry$surv_incl$n, entry$surv_incl$events, entry$surv_excl$n, entry$surv_excl$events))

  # ---- manifest 비교 ----
  old <- manifest[[cancer]]
  if (is.null(old)) {
    cat("manifest: 새 항목 생성 (md5", entry$md5, ")\n")
  } else if (identical(old$md5, entry$md5)) {
    if (counts_changed(old, entry)) {
      stop(cancer, ": .sav는 그대로(md5 동일)인데 결과가 달라짐 → 코드 변경이 원인\n  ",
           paste(manifest_diff(old, entry), collapse = "\n  "),
           "\n  의도한 변경이면 ", cfg$manifest, " 에서 ", cancer, " 항목을 삭제하고 다시 실행")
    }
    d <- manifest_diff(old, entry)
    cat("manifest: .sav 변경 없음, n/생존 수 동일", if (length(d)) paste0("(열/유전자 변경: ", paste(d, collapse = "; "), ")") else "", "\n")
  } else {
    d <- manifest_diff(old, entry)
    cat("manifest: .sav 변경됨 (md5 ", old$md5, " → ", entry$md5, ")\n  ",
        if (length(d)) paste(d, collapse = "\n  ") else "내용 차이 없음", "\n", sep = "")
  }
  manifest[[cancer]] <- entry

  # ---- 기존 유전자 그룹 QC (재계산 값으로 덮어쓰지 않음, 전체 코호트 기준) ----
  qc <- bind_rows(lapply(genes, function(g) {
    data.frame(gene = toupper(g),
               check_median_split(df[[paste0(g, "_expression")]], df[[paste0(g, "_group")]]))
  }))
  print(qc, row.names = FALSE)
  if (any(!qc$median_split_pass)) {
    warning(cancer, ": median split과 ", cfg$median_tie_tolerance, "명 넘게 다른 유전자 → ",
            paste(qc$gene[!qc$median_split_pass], collapse = ", "))
  }
  save_table(qc, cancer, "clinical_QC_median_split.csv")

  # ---- 저장 ----
  saveRDS(df, processed_path(cancer, "clinical_all.rds"))
  analytic <- if (cfg$exclude_neoadjuvant) df[!df$neoadjuvant_flag, ] else df
  cat("분석용 환자:", nrow(df), "→", nrow(analytic),
      if (cfg$exclude_neoadjuvant) "(신보조요법 제외)" else "(제외 없음)", "\n")
  saveRDS(analytic, processed_path(cancer, "clinical.rds"))
}

# 모든 암종이 통과한 뒤에만 manifest 저장
dir.create(dirname(cfg$manifest), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(manifest, cfg$manifest, auto_unbox = TRUE, pretty = TRUE)
cat("\nmanifest 저장:", cfg$manifest, "\n")
