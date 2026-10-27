# =============================================================================
# test_geo_inspect.R
# GEO 데이터셋을 받아서 raw data(발현 행렬, phenotype 필드, probe)를 확인하는 테스트 스크립트
# - 분석은 하지 않고, 06_external_geo.R config 작성에 필요한 정보만 출력/저장
# - 사용법: 아래 PARAMETERS만 바꾸고 source("test_geo_inspect.R")
# =============================================================================

# ---- PARAMETERS -------------------------------------------------------------
gse      <- "GSE39582"      # 확인할 GSE ID
platform <- "GPL570"        # 여러 플랫폼이 섞인 GSE일 때 선택용 (NULL이면 첫 번째)
gene     <- "ROCK2"         # probe를 확인할 유전자
geo_dir  <- "data/geo"      # 다운로드 캐시 폴더 (git-ignore)
out_dir  <- file.path("output/tables/GEO", gse, "inspect")

# ---- PACKAGES ---------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("GEOquery", "Biobase")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE)
}
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

options(timeout = 1200)                          # 큰 파일 다운로드 대비
Sys.setenv(VROOM_CONNECTION_SIZE = 131072 * 64)  # 긴 series matrix 줄 읽기 오류 방지

dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. DOWNLOAD (cache) ----------------------------------------------------
rds_file <- file.path(geo_dir, paste0(gse, "_eset.rds"))

if (file.exists(rds_file)) {
  message("[cache] ", rds_file, " 불러옴")
  eset <- readRDS(rds_file)
} else {
  message("[download] ", gse, " -> ", geo_dir, " (수 분 걸릴 수 있음)")
  gse_list <- getGEO(gse, GSEMatrix = TRUE, destdir = geo_dir, getGPL = TRUE)
  message("series matrix 개수: ", length(gse_list), "  (", paste(names(gse_list), collapse = ", "), ")")
  idx <- 1
  if (!is.null(platform) && length(gse_list) > 1) {
    idx <- grep(platform, names(gse_list))
    if (length(idx) != 1) stop("platform '", platform, "'에 해당하는 matrix를 하나로 특정할 수 없음")
  }
  eset <- gse_list[[idx]]
  saveRDS(eset, rds_file)
}

cat("\n=============================================================\n")
cat(gse, " | platform: ", annotation(eset), "\n", sep = "")
cat("=============================================================\n")

# ---- 2. EXPRESSION MATRIX ---------------------------------------------------
ex <- exprs(eset)
cat("\n[Expression] dim (probes x samples): ", nrow(ex), " x ", ncol(ex), "\n", sep = "")
cat("NA 개수: ", sum(is.na(ex)), "\n", sep = "")
q <- quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
print(round(q, 3))
is_log2 <- max(ex, na.rm = TRUE) <= 100
cat("log2 scale 추정: ", is_log2,
    if (!is_log2) "  -> 06에서는 log2 변환 필요" else "", "\n", sep = "")
cat("\n발현값 미리보기 (5 probes x 5 samples):\n")
print(round(ex[1:5, 1:5], 3))

# ---- 3. PHENOTYPE: raw fields -----------------------------------------------
pd <- pData(eset)
cat("\n[Phenotype] samples: ", nrow(pd), " | columns: ", ncol(pd), "\n", sep = "")
cat("\n전체 pData 컬럼명:\n")
print(colnames(pd))

# 정상/종양 구분용 후보 필드 (어떤 필드로 normal을 빼야 할지 확인)
for (f in intersect(c("source_name_ch1", "characteristics_ch1"), colnames(pd))) {
  cat("\n---", f, "(상위 값) ---\n")
  print(head(sort(table(pd[[f]], useNA = "ifany"), decreasing = TRUE), 10))
}
cat("\n--- title 예시 ---\n")
print(head(pd$title, 10))

# 첫 샘플의 raw characteristics 원문 ("key: value" 형태)
char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
cat("\n--- 첫 샘플 raw characteristics (", length(char_cols), "개) ---\n", sep = "")
print(unlist(pd[1, char_cols]), quote = FALSE)

# GEOquery가 "key:ch1" 형태로 파싱한 필드 요약 -> config의 field mapping에 쓸 이름
ch1_cols <- grep(":ch1$", colnames(pd), value = TRUE)
is_missing <- function(x) is.na(x) | trimws(x) %in% c("", "NA", "N/A", "na", "n/a", "NaN")

field_summary <- do.call(rbind, lapply(ch1_cols, function(f) {
  x <- pd[[f]]
  u <- unique(x[!is_missing(x)])
  data.frame(
    field     = f,
    key       = sub(":ch1$", "", f),
    n_unique  = length(u),
    n_missing = sum(is_missing(x)),
    numeric   = length(u) > 0 && all(!is.na(suppressWarnings(as.numeric(u)))),
    examples  = paste(head(u, 6), collapse = " | "),
    stringsAsFactors = FALSE
  )
}))

cat("\n[Parsed fields] (config fields 매핑에 사용할 key)\n")
print(field_summary[, c("key", "n_unique", "n_missing", "numeric", "examples")], right = FALSE)

# 범주형 필드(값 종류 <= 15)는 빈도표 출력
cat("\n[범주형 필드 빈도]\n")
for (f in field_summary$field[!field_summary$numeric & field_summary$n_unique <= 15]) {
  cat("\n--", sub(":ch1$", "", f), "--\n")
  print(table(pd[[f]], useNA = "ifany"))
}

# 수치형 필드(생존시간, 나이 등) 요약
cat("\n[수치형 필드 요약]\n")
for (f in field_summary$field[field_summary$numeric]) {
  cat("\n--", sub(":ch1$", "", f), "--\n")
  print(summary(suppressWarnings(as.numeric(pd[[f]]))))
}

# ---- 4. PROBES for the gene -------------------------------------------------
fd <- fData(eset)
cat("\n[Feature annotation] 컬럼명:\n")
print(colnames(fd))

sym_col <- intersect(c("Gene Symbol", "Gene.Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment"),
                     colnames(fd))[1]
if (is.na(sym_col)) {
  warning("Gene symbol 컬럼을 찾지 못함 -> 위 컬럼명 확인 후 sym_col 직접 지정")
} else {
  cat("사용한 symbol 컬럼: ", sym_col, "\n", sep = "")
  # "ROCK2 /// XXX" 같은 다중 매핑도 정확히 일치하는 경우만 잡음
  sym_split <- strsplit(as.character(fd[[sym_col]]), "\\s*///\\s*")
  hit <- vapply(sym_split, function(s) gene %in% s, logical(1))
  probes <- rownames(fd)[hit]

  cat("\n[", gene, "] 매칭 probe 수: ", length(probes), "\n", sep = "")
  if (length(probes) > 0) {
    sub_ex <- ex[probes, , drop = FALSE]
    probe_tab <- data.frame(
      probe_id = probes,
      symbol   = fd[probes, sym_col],
      mean     = round(rowMeans(sub_ex, na.rm = TRUE), 3),
      median   = round(apply(sub_ex, 1, median, na.rm = TRUE), 3),
      IQR      = round(apply(sub_ex, 1, IQR, na.rm = TRUE), 3),
      sd       = round(apply(sub_ex, 1, sd, na.rm = TRUE), 3),
      row.names = NULL
    )
    probe_tab <- probe_tab[order(-probe_tab$mean), ]
    print(probe_tab, row.names = FALSE)
    if (length(probes) > 1) {
      cat("\nprobe 간 상관 (Pearson):\n")
      print(round(cor(t(sub_ex), use = "pairwise.complete.obs"), 3))
    }
    write.csv(probe_tab, file.path(out_dir, paste0("probes_", gene, ".csv")), row.names = FALSE)
  }
}

# ---- 5. SAVE ----------------------------------------------------------------
write.csv(pd, file.path(out_dir, "pheno_raw.csv"), row.names = TRUE)
write.csv(field_summary, file.path(out_dir, "pheno_field_summary.csv"), row.names = FALSE)

cat("\n저장 위치: ", normalizePath(out_dir), "\n", sep = "")
cat("  - pheno_raw.csv            : pData 전체\n")
cat("  - pheno_field_summary.csv  : 필드별 unique/missing/예시값\n")
cat("  - probes_", gene, ".csv        : 매칭 probe QC\n", sep = "")
cat("캐시: ", rds_file, " (다음 실행부터 다운로드 생략)\n", sep = "")
