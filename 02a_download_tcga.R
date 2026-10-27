# =============================================================
# 02a_download_tcga.R
# GDC에서 TCGA RNA-seq(STAR - Counts) 다운로드 + 전처리 (암종은 config.R의 cfg$cancers)
# 사용:
#   Rscript 02a_download_tcga.R                        # 전체 암종 다운로드 + 전처리 + provenance
#   Rscript 02a_download_tcga.R COAD                   # 지정 암종만
#   Rscript 02a_download_tcga.R --provenance-only [COAD]
#       다운로드/GDCprepare 없이 기존 data/processed/<C>_se.rda + GDC 캐시로 provenance만 생성
#       (GDC API에는 파일 메타데이터 조회만 함; counts/발현 파일은 쓰지 않음)
# 결과:
#   database/TCGA_<COHORT>_RNAseq_Expression.csv  (log2(TPM+1), 유전자 x 환자)
#   data/processed/<COHORT>_counts.rds            (raw counts, GSEA/DESeq2용)
#   data/processed/<COHORT>_se.rda                (GDCprepare 결과 전체)
#   data/processed/<COHORT>_gdc_provenance.json   (release, 쿼리, 샘플 필터, 유전자 모델, 버전)
#   data/processed/<COHORT>_gdc_files.csv         (파일별 file_id, 바코드, md5, updated_datetime, 캐시 대조)
# =============================================================

# ---- 0. 패키지 (최초 1회만) ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("TCGAbiolinks", "SummarizedExperiment")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE)
}
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

source("config.R")
source("R/utils.R")

# ---- 명령줄 인수: --provenance-only, 암종 ----
args <- commandArgs(trailingOnly = TRUE)
provenance_only <- "--provenance-only" %in% args
cancer_args <- toupper(setdiff(args, "--provenance-only"))
if (length(cancer_args)) {
  bad <- setdiff(cancer_args, cfg$cancers)
  if (length(bad)) stop("cfg$cancers에 없는 암종: ", paste(bad, collapse = ", "))
  cfg$cancers <- cancer_args
  cat("명령줄 지정 암종만 실행:", paste(cancer_args, collapse = ", "), "\n")
}
if (provenance_only) cat("--provenance-only: 다운로드/GDCprepare 없이 기존 _se.rda로 provenance만 생성\n")

dir.create(cfg$clinical_dir, showWarnings = FALSE)
dir.create(cfg$processed_dir, recursive = TRUE, showWarnings = FALSE)

# Windows 경로 길이 문제를 피하려고 짧은 폴더에 받음
gdc_dir <- cfg$gdc_dir
gq <- cfg$gdc

expr_csv_path <- function(cohort) file.path(cfg$clinical_dir, paste0("TCGA_", cohort, "_RNAseq_Expression.csv"))
se_path       <- function(cohort) processed_path(cohort, "se.rda")
fmt_time      <- function(t) format(t, "%Y-%m-%d %H:%M:%S %Z")

# GDC 캐시 폴더 (TCGAbiolinks 규칙: <gdc_dir>/<project>/<category>/<type>/<file_id>/<file>)
cache_dir <- function(project) {
  file.path(gdc_dir, project, gsub(" ", "_", gq$data.category), gsub(" ", "_", gq$data.type))
}

# ---- 1. 쿼리 (GDCquery 메시지에서 참조 유전체 문구도 수집) ----
gdc_query <- function(project) {
  msgs <- character()
  q <- withCallingHandlers(
    GDCquery(project = project, data.category = gq$data.category,
             data.type = gq$data.type, workflow.type = gq$workflow.type),
    message = function(m) msgs <<- c(msgs, conditionMessage(m)))
  ref <- regmatches(msgs, regexpr("Genome of reference:\\s*\\S+", msgs))
  attr(q, "genome_message") <- if (length(ref)) trimws(ref[1]) else NA_character_
  q
}

# ---- 2. 종양 샘플 + 환자 단위 중복 정리 + 유전자 기호 (다운로드/provenance 공통) ----
process_se <- function(se) {
  info <- list(n_samples_all = ncol(se))

  # 종양 샘플만
  bc <- colnames(se)
  se <- se[, substr(bc, gq$sample_code_pos[1], gq$sample_code_pos[2]) == gq$tumor_sample_code]
  info$n_tumor_samples <- ncol(se)

  # 환자 단위로 중복 정리: 첫 번째 aliquot 유지
  pid <- substr(colnames(se), 1, gq$patient_id_chars)
  dup <- duplicated(pid)
  if (any(dup)) {
    message("중복 환자 ", length(unique(pid[dup])), "명 (추가 aliquot ", sum(dup), "개) → 첫 aliquot만 유지: ",
            paste(unique(pid[dup]), collapse = ", "))
  }
  dup_pat <- unique(pid[dup])
  info$duplicate_patients <- lapply(dup_pat, function(p) {
    al <- colnames(se)[pid == p]
    list(patient = p, kept_aliquot = al[1], dropped_aliquots = al[-1])
  })
  se <- se[, !dup]
  info$kept_aliquots <- colnames(se)
  colnames(se) <- substr(colnames(se), 1, gq$patient_id_chars)
  info$n_patients <- ncol(se)

  # 유전자 ID → Symbol, 중복 Symbol은 평균 발현이 가장 높은 것만
  sym  <- rowData(se)$gene_name
  tpm  <- assay(se, gq$tpm_assay)
  keep <- !is.na(sym) & sym != ""
  tpm  <- tpm[keep, ]; sym <- sym[keep]
  ord  <- order(rowMeans(tpm), decreasing = TRUE)
  first <- ord[!duplicated(sym[ord])]
  tpm  <- tpm[first, ]
  rownames(tpm) <- sym[first]

  cnt <- assay(se, gq$count_assay)[keep, ][first, ]
  rownames(cnt) <- sym[first]

  info$n_features_all <- nrow(se)
  info$n_genes <- nrow(cnt)
  list(expr = log2(tpm + 1), cnt = cnt, info = info)
}

# ---- 3. provenance ----

# 캐시 파일: file_id, 경로, mtime, md5, 유전자 모델 헤더 ("# gene-model: ...")
scan_cache <- function(project) {
  d <- cache_dir(project)
  if (!dir.exists(d)) {
    warning(project, ": GDC 캐시 폴더 없음 (", d, ")", call. = FALSE)
    return(data.frame(file_id = character(), cache_path = character()))
  }
  ids <- basename(list.dirs(d, recursive = FALSE))
  paths <- vapply(ids, function(i) {
    f <- list.files(file.path(d, i), full.names = TRUE)
    f <- f[!dir.exists(f) & !grepl("^\\.|\\.partial$", basename(f))]
    if (length(f)) f[1] else NA_character_
  }, "")
  hdr <- vapply(paths, function(f) {
    if (is.na(f)) return(NA_character_)
    l <- readLines(f, n = 1, warn = FALSE)
    if (length(l) && startsWith(l, "#")) trimws(sub("^#\\s*", "", l)) else NA_character_
  }, "")
  data.frame(file_id = ids, cache_path = unname(paths),
             cache_mtime = file.mtime(paths),
             cache_md5 = unname(tools::md5sum(paths)),
             cache_header = unname(hdr), stringsAsFactors = FALSE)
}

build_provenance <- function(cohort, se_raw, query, proc, mode, times) {
  project <- paste0(gq$project_prefix, cohort)
  now <- Sys.time()

  # 파일 목록: 쿼리 결과 + 캐시 대조
  res <- if (!is.null(query)) getResults(query) else NULL
  cache <- scan_cache(project)
  if (is.null(res)) {
    files <- cache
    files$in_query <- NA
  } else {
    cols <- intersect(c("file_id", "file_name", "cases", "sample.submitter_id", "sample_type",
                        "md5sum", "file_size", "version", "created_datetime", "updated_datetime",
                        "data_release", "state", "analysis_workflow_type", "analysis_workflow_version"),
                      names(res))
    files <- merge(res[, cols], cache, by = "file_id", all = TRUE)
    files$in_query <- files$file_id %in% res$file_id
    names(files)[names(files) == "cases"] <- "barcode"
  }
  files$in_cache  <- files$file_id %in% cache$file_id
  files$md5_match <- if ("md5sum" %in% names(files)) files$md5sum == files$cache_md5 else NA
  bcol <- if ("barcode" %in% names(files)) files$barcode else rep(NA_character_, nrow(files))
  files$in_se     <- bcol %in% colnames(se_raw)
  files$used      <- bcol %in% proc$info$kept_aliquots
  files <- files[order(!files$used, files$barcode %||% files$file_id), ]

  # release: 다운로드 시점 값 우선 (full 모드의 쿼리 시점 또는 metadata(se))
  current <- tryCatch(getGDCInfo(), error = function(e) NULL)
  meta_release <- metadata(se_raw)$data_release %||% NA_character_
  if (mode == "download") {
    release <- current$data_release %||% meta_release
    release_source <- "getGDCInfo() at download"
  } else if (!is.na(meta_release)) {
    release <- meta_release
    release_source <- "metadata(se)$data_release (recorded by GDCprepare at download)"
  } else {
    release <- current$data_release %||% NA_character_
    release_source <- paste0("queried_at_", format(now, "%Y-%m-%d"), ", not at download")
    warning(cohort, ": _se.rda에 다운로드 당시 release가 없음 → 현재 release(", release,
            ")를 기록. 다운로드 시점과 다를 수 있음", call. = FALSE)
  }
  if (!is.null(current) && !is.na(meta_release) && current$data_release != meta_release) {
    warning(cohort, ": 현재 GDC release(", current$data_release, ")가 다운로드 당시(", meta_release,
            ")와 다름 → gdc_files.csv의 md5_match로 파일 변경 여부 확인", call. = FALSE)
  }

  gm <- unique(na.omit(sub("^gene-model:\\s*", "", files$cache_header[grepl("^gene-model:", files$cache_header)])))
  gnm <- unique(genome(rowRanges(se_raw)))
  outputs <- c(counts = processed_path(cohort, "counts.rds"), expression_csv = expr_csv_path(cohort))

  prov <- list(
    cancer = cohort,
    mode = mode,
    provenance_created = fmt_time(now),
    gdc_data_release = list(
      release = release,
      release_source = release_source,
      metadata_se = meta_release,
      current_at_provenance = current$data_release %||% NA,
      current_commit = current$commit %||% NA
    ),
    timing = times,
    query = list(project = project, data.category = gq$data.category, data.type = gq$data.type,
                 workflow.type = gq$workflow.type,
                 queried_at = if (!is.null(query)) fmt_time(attr(query, "queried_at")) else NA,
                 n_files = if (!is.null(res)) nrow(res) else NA),
    assays = list(counts = gq$count_assay, tpm = gq$tpm_assay, se_assays = assayNames(se_raw)),
    reference = list(
      gene_model = if (length(gm)) gm else NA,
      gene_model_source = "first line of cached GDC STAR-Counts TSV files",
      genome_build = if (length(na.omit(gnm))) na.omit(gnm) else NA,
      genome_build_source = if (length(na.omit(gnm))) "genome(rowRanges(se))" else
        paste0("genome(rowRanges(se)) is NA; TCGAbiolinks GDCquery message: ",
               attr(query, "genome_message") %||% "not available")
    ),
    samples = list(
      n_se_columns = proc$info$n_samples_all,
      tumor_filter = paste0("barcode characters ", gq$sample_code_pos[1], "-", gq$sample_code_pos[2],
                            " == \"", gq$tumor_sample_code, "\""),
      n_tumor_samples = proc$info$n_tumor_samples,
      duplicate_rule = "first aliquot per patient (column order of GDCprepare)",
      duplicate_patients = proc$info$duplicate_patients,
      n_patients = proc$info$n_patients
    ),
    genes = list(n_features_se = proc$info$n_features_all, n_genes_final = proc$info$n_genes,
                 rule = "non-empty gene_name; duplicate symbols keep the highest mean TPM"),
    final_matrix = sprintf("%d genes x %d patients", proc$info$n_genes, proc$info$n_patients),
    cache = list(
      dir = cache_dir(project),
      n_files = nrow(cache),
      mtime_min = if (nrow(cache)) fmt_time(min(cache$cache_mtime)) else NA,
      mtime_max = if (nrow(cache)) fmt_time(max(cache$cache_mtime)) else NA,
      n_query_files_in_cache = sum(files$in_query & files$in_cache, na.rm = TRUE),
      n_md5_match = sum(files$md5_match, na.rm = TRUE),
      n_cache_only = sum(files$in_cache & !files$in_query, na.rm = TRUE),
      n_query_only = sum(files$in_query & !files$in_cache, na.rm = TRUE)
    ),
    se_rda = list(path = se_path(cohort), mtime = fmt_time(file.mtime(se_path(cohort)))),
    outputs_md5 = as.list(setNames(unname(tools::md5sum(outputs[file.exists(outputs)])),
                                   names(outputs)[file.exists(outputs)])),
    software = list(R = R.version.string,
                    TCGAbiolinks = as.character(packageVersion("TCGAbiolinks")),
                    SummarizedExperiment = as.character(packageVersion("SummarizedExperiment")))
  )
  if (!is.null(proc$counts_match_existing)) prov$outputs_counts_identical <- proc$counts_match_existing

  files$cache_mtime <- fmt_time(files$cache_mtime)
  files$cache_path <- NULL
  readr::write_excel_csv(files, processed_path(cohort, "gdc_files.csv"), na = "")
  jsonlite::write_json(prov, processed_path(cohort, "gdc_provenance.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  cat(cohort, ": provenance → ", processed_path(cohort, "gdc_provenance.json"), ", ",
      processed_path(cohort, "gdc_files.csv"), "\n", sep = "")
  if (prov$cache$n_md5_match < sum(files$used)) {
    warning(cohort, ": 사용한 파일 중 캐시 md5가 GDC 현재 md5와 다른 것이 있음 → gdc_files.csv 확인", call. = FALSE)
  }
  invisible(prov)
}

# ---- 4. 암종 1개 ----
process_cohort <- function(cohort) {
  project <- paste0(gq$project_prefix, cohort)
  message("==== ", project, " ====")
  times <- list()

  if (provenance_only) {
    f <- se_path(cohort)
    if (!file.exists(f)) stop(cohort, ": ", f, " 없음 → --provenance-only 없이 02a를 실행해 다운로드/prepare 먼저")
    env <- new.env()
    obj <- load(f, envir = env)
    se_raw <- env[[obj[1]]]
    if (!is(se_raw, "SummarizedExperiment")) stop(f, ": SummarizedExperiment가 아님 (객체 ", obj[1], ")")
    times$download_prepare <- "not re-run (--provenance-only); see se_rda.mtime and cache.mtime_min/max"

    # 메타데이터 조회만 (파일 다운로드 없음). 실패하면 캐시 정보만 기록
    query <- tryCatch({
      q <- gdc_query(project); attr(q, "queried_at") <- Sys.time(); q
    }, error = function(e) {
      warning(cohort, ": GDCquery 실패 (", conditionMessage(e), ") → 파일 목록은 캐시만 기록", call. = FALSE)
      NULL
    })

    proc <- process_se(se_raw)
    # 기존 counts와 같은 결과인지 확인 (파일은 쓰지 않음)
    cf <- processed_path(cohort, "counts.rds")
    proc$counts_match_existing <- if (file.exists(cf)) identical(proc$cnt, readRDS(cf)) else NA
    if (isFALSE(proc$counts_match_existing)) {
      warning(cohort, ": _se.rda에서 다시 만든 counts가 기존 ", cf, "와 다름", call. = FALSE)
    }
    build_provenance(cohort, se_raw, query, proc, "provenance-only", times)
    return(invisible(NULL))
  }

  # ---- 쿼리 ----
  times$query <- fmt_time(Sys.time())
  query <- gdc_query(project)
  attr(query, "queried_at") <- Sys.time()

  # ---- 다운로드 (중간에 끊기면 다시 실행하면 이어서 받음) ----
  times$download_start <- fmt_time(Sys.time())
  GDCdownload(query, method = "api", directory = gdc_dir, files.per.chunk = gq$files_per_chunk)
  times$download_end <- fmt_time(Sys.time())

  # ---- SummarizedExperiment로 합치기 ----
  se_raw <- GDCprepare(query, directory = gdc_dir, save = TRUE, save.filename = se_path(cohort))
  times$prepare_end <- fmt_time(Sys.time())

  # ---- 종양 / 중복 / 유전자 ----
  proc <- process_se(se_raw)

  # ---- 저장 ----
  expr <- proc$expr
  out_csv <- expr_csv_path(cohort)
  write.csv(data.frame(Gene = rownames(expr), expr, check.names = FALSE),
            out_csv, row.names = FALSE)
  saveRDS(proc$cnt, processed_path(cohort, "counts.rds"))
  message(cohort, ": ", nrow(expr), " genes x ", ncol(expr), " patients → ", out_csv)

  build_provenance(cohort, se_raw, query, proc, "download", times)
  invisible(expr)
}

for (cohort in cfg$cancers) process_cohort(cohort)

# ---- 5. 기존 .sav와 겹치는 환자 수 확인 ----
if (!provenance_only) {
  for (cohort in cfg$cancers) {
    sav <- sav_path(cohort)
    if (!file.exists(sav)) next
    clin <- haven::read_sav(sav)
    expr <- read.csv(expr_csv_path(cohort), check.names = FALSE, nrows = 1)
    ids  <- trimws(clin[[resolve_schema(cohort)$id]])
    cat(cohort, ": 임상", length(ids), "명 중 발현 매칭",
        sum(ids %in% colnames(expr)), "명\n")
  }
}
