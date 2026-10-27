# =============================================================
# R/utils.R
# 모든 스크립트가 공유하는 함수 (암종/유전자 무관)
# 사용: source("config.R"); source("R/utils.R")
# =============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ---- 상수 --------------------------------------------------------------

# TCGA 결측 표기
na_strings <- c("", "[Not Available]", "[Unknown]", "[Not Applicable]",
                "[Discrepancy]", "[Not Evaluated]")

# 스키마 항목 → 표준 열 이름 (원본 열 이름은 cfg$schema_default / cfg$schema)
canon <- c(id = "sample_id", time = "surv_time", status = "status",
           stage = "pathologic_stage", age = "age", sex = "gender",
           neoadjuvant = "neoadjuvant", last_contact = "last_contact")
schema_required <- c("id", "time", "status", "stage", "age", "sex")

# 값 라벨 / 라벨 일관성 규칙: 기본값(cfg$*_default) + 암종별 재정의
labels_for <- function(cancer) modifyList(cfg$value_labels_default, cfg$value_labels[[cancer]] %||% list())
rules_for  <- function(cancer) modifyList(cfg$label_rules_default, cfg$label_rules[[cancer]] %||% list())

# ---- 공변량 / 노출 그룹 (설정은 config.R) -------------------------------

group_ref <- function() cfg$group_levels[1]   # 기준 그룹 (Low)
group_alt <- function() cfg$group_levels[2]   # 비교 그룹 (High)
group_contrast_label <- function() paste(group_alt(), "vs", group_ref())

# 연속형 공변량의 모형 항 이름 (age → age_per10) / 라벨
scaled_term <- function(v) {
  s <- cfg$covariate_scale[[v]]
  if (is.null(s)) v else paste0(v, "_per", s$by)
}
scaled_label <- function(term) {
  for (v in names(cfg$covariate_scale)) {
    if (identical(term, scaled_term(v))) return(cfg$covariate_scale[[v]]$label)
  }
  NULL
}
add_scaled_terms <- function(d) {
  for (v in intersect(names(cfg$covariate_scale), names(d))) {
    d[[scaled_term(v)]] <- d[[v]] / cfg$covariate_scale[[v]]$by
  }
  d
}

# cfg$collapse_small_levels에 있는 변수: 소수 수준 → cfg$collapse_other_label
collapse_small <- function(x, var) {
  n_min <- cfg$collapse_small_levels[var]
  if (is.null(n_min) || is.na(n_min)) return(x)
  n <- table(x)
  ifelse(is.na(x), NA, ifelse(x %in% names(n)[n >= n_min], as.character(x), cfg$collapse_other_label))
}

# cfg$reference_levels 적용 (존재하는 factor 변수, 해당 수준이 있을 때만)
apply_reference_levels <- function(df) {
  for (v in intersect(names(cfg$reference_levels), names(df))) {
    ref <- cfg$reference_levels[[v]]
    if (is.factor(df[[v]]) && ref %in% levels(df[[v]])) {
      lab <- attr(df[[v]], "label")
      df[[v]] <- relevel(df[[v]], ref = ref)
      attr(df[[v]], "label") <- lab
    }
  }
  df
}

# Table/그림용 변수 라벨
var_labels <- c(
  age                                 = "Age, years",
  age_g                               = "Age group, years",
  gender                              = "Sex",
  stage                               = "Pathologic stage",
  pathologic_stage                    = "Pathologic stage",
  pathologic_stage_12_34              = "Pathologic stage group",
  pathologic_t                        = "Pathologic T",
  pathologic_n                        = "Pathologic N",
  pathologic_m                        = "Pathologic M",
  cea_g                               = "Pretreatment CEA",
  histologic_diagnosis                = "Histology",
  anatomic_neoplasm_subdivision       = "Tumor site",
  residual_tumor                      = "Residual tumor",
  lymphovascular_invasion_indicator   = "Lymphovascular invasion",
  vascular_invasion_indicator         = "Vascular invasion",
  perineural_invasion                 = "Perineural invasion",
  tumor_status                        = "Tumor status",
  tss                                 = "Tissue source site",
  # 06 GEO 필드 (결과 변수: GEO 데이터셋이 실제 판정 결과를 제공)
  mmr_status                          = "MMR status",
  kras_mutation                       = "KRAS mutation",
  braf_mutation                       = "BRAF mutation",
  tumor_location                      = "Tumor location",
  adjuvant_chemo                      = "Adjuvant chemotherapy",
  cit_subtype                         = "CIT molecular subtype",
  # 검사 "시행 여부" 변수 — 변이/MMR 결과가 아님
  kras_gene_analysis_indicator        = "KRAS testing performed",
  braf_gene_analysis_indicator        = "BRAF testing performed",
  mismatch_rep_proteins_tested_by_ihc = "MMR IHC performed"
)

# ---- 파일/경로 ---------------------------------------------------------

# subdir: output/<type>/<cancer>/ 아래 하위 폴더 (예: "gsea/MS4A1")
out_dir <- function(cancer, type = c("tables", "figures"), subdir = NULL) {
  type <- match.arg(type)
  d <- do.call(file.path, as.list(c(cfg$output_dir, type, cancer, subdir)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

processed_path <- function(cancer, suffix) {
  dir.create(cfg$processed_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(cfg$processed_dir, paste0(cancer, "_", suffix))
}

# 결과 파일 이름: <prefix>_<분석>[_<GENE>] (prefix 기본값 = 암종, 06은 GSE ID)
out_file <- function(prefix, name) paste0(prefix, "_", name)

# UTF-8 (BOM) CSV → Excel에서 ≤, – 등 기호가 깨지지 않음
save_table <- function(df, cancer, name, subdir = NULL, prefix = cancer) {
  path <- file.path(out_dir(cancer, "tables", subdir), out_file(prefix, name))
  readr::write_excel_csv(df, path, na = "")
  invisible(path)
}

# data.frame → <prefix>_<stem>.docx (flextable) + .csv. 생성된 경로 반환
save_df_table <- function(df, cancer, stem, caption = NULL, landscape = FALSE, font_size = 9,
                          subdir = NULL, prefix = cancer) {
  ft <- flextable::flextable(df) |>
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = font_size, part = "all") |>
    flextable::bold(part = "header") |>
    flextable::autofit()
  if (!is.null(caption)) ft <- flextable::set_caption(ft, caption)   # 일반 텍스트 (마크다운 X)
  sect <- if (landscape) {
    officer::prop_section(page_size = officer::page_size(orient = "landscape"))
  }
  docx <- file.path(out_dir(cancer, "tables", subdir), out_file(prefix, paste0(stem, ".docx")))
  flextable::save_as_docx(ft, path = docx, pr_section = sect)
  c(docx, save_table(df, cancer, paste0(stem, ".csv"), subdir, prefix))
}

# ggplot/patchwork → <prefix>_<stem>.pdf (cairo, 유니코드 OK) + .tiff (300 dpi, LZW)
save_fig <- function(plot, cancer, stem, width = 7, height = 6, dpi = 300, subdir = NULL,
                     prefix = cancer) {
  dir <- out_dir(cancer, "figures", subdir)
  pdf <- file.path(dir, out_file(prefix, paste0(stem, ".pdf")))
  tif <- file.path(dir, out_file(prefix, paste0(stem, ".tiff")))
  ggplot2::ggsave(pdf, plot, width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(tif, plot, width = width, height = height, dpi = dpi,
                  device = "tiff", compression = "lzw")
  c(pdf, tif)
}

# 이번 실행에서 만들지 않은 결과 파일(<prefix>_<prefixes>로 시작; NULL = 모든 파일) → <dir>/_stale/
# <prefix>_로 시작하지 않는 이전 이름 형식 파일도 함께 이동
# (삭제하지 않고 이동; Word 잠금 파일 ~$* 와 하위 폴더는 무시)
move_stale_outputs <- function(cancer, prefixes, keep, subdir = NULL, prefix = cancer) {
  norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  keep <- norm(keep)
  moved <- character()
  for (type in c("tables", "figures")) {
    dir <- out_dir(cancer, type, subdir)
    f <- list.files(dir, full.names = TRUE)
    f <- f[!dir.exists(f) & !startsWith(basename(f), "~$")]
    legacy <- !startsWith(basename(f), out_file(prefix, ""))
    has_prefix <- if (is.null(prefixes)) rep(TRUE, length(f)) else
      Reduce(`|`, lapply(prefixes, function(p) startsWith(basename(f), out_file(prefix, p))), legacy)
    stale <- f[has_prefix & !norm(f) %in% keep]
    if (length(stale)) {
      dir.create(file.path(dir, "_stale"), showWarnings = FALSE)
      dest <- file.path(dir, "_stale", basename(stale))
      unlink(dest)
      file.rename(stale, dest)
      moved <- c(moved, stale)
    }
  }
  if (length(moved)) {
    cat("이전 결과 파일", length(moved), "개 → _stale/ 이동:",
        paste(basename(moved), collapse = ", "), "\n")
  }
  invisible(moved)
}

# 이번 실행에 포함되지 않은 하위 폴더(예: 이전에 분석한 유전자) → <parent>/_stale/<폴더>
move_stale_dirs <- function(cancer, parent, keep_dirs) {
  moved <- character()
  for (type in c("tables", "figures")) {
    base <- out_dir(cancer, type, parent)
    d <- list.dirs(base, recursive = FALSE, full.names = TRUE)
    stale <- d[!basename(d) %in% c(keep_dirs, "_stale")]
    for (s in stale) {
      dest <- file.path(base, "_stale", basename(s))
      dir.create(dirname(dest), showWarnings = FALSE)
      unlink(dest, recursive = TRUE)
      file.rename(s, dest)
      moved <- c(moved, s)
    }
  }
  if (length(moved)) cat("이전 결과 폴더 → _stale/ 이동:", paste(unique(basename(moved)), collapse = ", "), "\n")
  invisible(moved)
}

# ---- 숫자 형식 (논문용) ------------------------------------------------

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_hr <- function(hr, lo, hi) {
  ifelse(is.na(hr), NA_character_, sprintf("%.2f (%.2f–%.2f)", hr, lo, hi))
}

# ---- 열 이름 -----------------------------------------------------------

# 유전자 기호 → 열 이름 키 (MS4A1 → ms4a1, HLA-A → hla_a)
gene_key <- function(gene) tolower(gsub("[^A-Za-z0-9]", "_", gene))

# 유전자 열(<GENE>expression / <GENE>Group 등, 대소문자 무관) → <gene>_expression / <gene>_group
# 나머지 열 → janitor::clean_names (ageG/AgeG → age_g 등)
# clean_names만 쓰면 APCexpression → ap_cexpression 이 되므로 유전자 열을 먼저 처리
standardize_names <- function(df) {
  nm     <- names(df)
  suffix <- str_match(nm, regex("(expression|group)$", ignore_case = TRUE))[, 2]
  stem   <- str_remove(nm, regex("(expression|group)$", ignore_case = TRUE))
  # 유전자 기호: 대문자/숫자로 시작, 밑줄 없음, 대문자 1개 이상
  is_gene <- !is.na(suffix) &
    str_detect(stem, "^[A-Z0-9][A-Za-z0-9.-]*$") & str_detect(stem, "[A-Z]")

  new <- nm
  new[is_gene]  <- paste0(gene_key(stem[is_gene]), "_", tolower(suffix[is_gene]))
  new[!is_gene] <- janitor::make_clean_names(nm[!is_gene])
  if (anyDuplicated(new)) {
    stop("열 이름 표준화 후 중복: ", paste(unique(new[duplicated(new)]), collapse = ", "))
  }
  names(df) <- new
  df
}

# <gene>_expression 과 <gene>_group 이 모두 있는 유전자 키
detect_genes <- function(df) {
  g <- str_remove(grep("_group$", names(df), value = TRUE), "_group$")
  g[paste0(g, "_expression") %in% names(df)]
}

# ---- 값 정리 -----------------------------------------------------------

# 발현 그룹 → factor(Low, High). 0/1 숫자는 1 = High
as_group_factor <- function(x) {
  lv <- cfg$group_levels
  if (is.numeric(x)) {
    if (!all(x %in% c(0, 1, NA))) stop("그룹 값이 0/1이 아님: ", paste(unique(x), collapse = ", "))
    x <- lv[x + 1]                           # 0 → 기준(Low), 1 → 비교(High)
  }
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "Na")] <- NA
  m <- lv[match(tolower(x), tolower(lv))]    # 대소문자 무시하고 group_levels에 맞춤
  bad <- unique(x[!is.na(x) & is.na(m)])
  if (length(bad)) stop("알 수 없는 그룹 값: ", paste(bad, collapse = ", "),
                        " (cfg$group_levels = ", paste(lv, collapse = ", "), ")")
  factor(m, levels = lv)
}

# 공백 제거, TCGA 결측 표기 → NA, 숫자형 문자 → numeric,
# 중복 열 제거, 발현 그룹 → factor
clean_tcga <- function(df) {
  df <- df[, setdiff(names(df), cfg$drop_columns), drop = FALSE]

  df[] <- lapply(df, function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      x <- str_squish(x)
      x[x %in% na_strings] <- NA
      num <- suppressWarnings(as.numeric(x))
      if (any(!is.na(x)) && all(is.na(x) | !is.na(num))) x <- num   # "2.00", ".00" → 2, 0
    }
    x
  })

  for (g in detect_genes(df)) {
    df[[paste0(g, "_group")]] <- as_group_factor(df[[paste0(g, "_group")]])
  }
  df
}

# 코드 → 라벨 factor (매핑에 없는 값이 있으면 중단)
apply_value_labels <- function(df, labels) {
  for (v in intersect(names(labels), names(df))) {
    map <- labels[[v]]
    x   <- as.character(df[[v]])
    bad <- setdiff(unique(na.omit(x)), names(map))
    if (length(bad)) stop(v, ": 라벨 매핑에 없는 값 ", paste(bad, collapse = ", "))
    df[[v]] <- factor(unname(map[x]), levels = unname(map))
  }
  df
}

# 숫자 코드 → 접두어 붙인 factor (T 2 → "T2")
code_factor <- function(x, prefix) {
  lv <- sort(unique(na.omit(x)))
  factor(ifelse(is.na(x), NA, paste0(prefix, x)), levels = paste0(prefix, lv))
}

set_var_labels <- function(df, labels = var_labels) {
  for (v in intersect(names(labels), names(df))) attr(df[[v]], "label") <- labels[[v]]
  df
}

# ---- 스키마 -------------------------------------------------------------

# 암종별 스키마 = 기본값 + 재정의 (NA = 해당 열 없음)
resolve_schema <- function(cancer) {
  s <- modifyList(cfg$schema_default, cfg$schema[[cancer]] %||% list())
  if (!s$time_unit %in% c("days", "months")) {
    stop(cancer, ": time_unit은 \"days\" 또는 \"months\" (현재 ", s$time_unit,
         ") → cfg$schema$", cancer, "$time_unit")
  }
  if (!s$stage_format %in% c("numeric", "ajcc_text")) {
    stop(cancer, ": stage_format은 \"numeric\" 또는 \"ajcc_text\" → cfg$schema$", cancer, "$stage_format")
  }
  s
}

# 스키마에 지정된 원본 열 (필드 → 원본 열 이름, NA 필드는 제외)
schema_columns <- function(schema) {
  f <- intersect(names(canon), names(schema))
  cols <- vapply(f, function(k) as.character(schema[[k]] %||% NA), "")
  cols[!is.na(cols) & nzchar(cols)]
}

sav_path <- function(cancer) file.path(cfg$clinical_dir, paste0(cancer, ".sav"))

read_sav_raw <- function(cancer) {
  haven::read_sav(sav_path(cancer)) |>
    haven::as_factor(only_labelled = TRUE) |>
    haven::zap_formats() |>
    haven::zap_widths() |>
    as.data.frame()
}

# 스키마 열이 원본에 있는지 확인 → 없으면 고칠 위치를 알려주고 중단
check_schema_columns <- function(raw, schema, cancer) {
  cols <- schema_columns(schema)
  miss_req <- setdiff(schema_required, names(cols))
  if (length(miss_req)) {
    stop(cancer, ": 필수 스키마 항목이 비어 있음 (", paste(miss_req, collapse = ", "),
         ") → config.R의 cfg$schema$", cancer, " 에 지정")
  }
  for (k in names(cols)) {
    if (!cols[[k]] %in% names(raw)) {
      near <- agrep(cols[[k]], names(raw), max.distance = 0.3, value = TRUE, ignore.case = TRUE)
      stop(cancer, ": 열 '", cols[[k]], "' 없음 (schema$", k, ") → config.R의 cfg$schema$", cancer,
           "$", k, " 를 실제 열 이름으로 수정",
           if (k %in% schema_required) "" else " (없는 항목이면 NA로 지정)",
           if (length(near)) paste0(". 비슷한 열: ", paste(head(near, 5), collapse = ", ")) else "")
    }
  }
  invisible(cols)
}

# 스키마 열 → 표준 이름, 나머지 → standardize_names
# (표준 이름과 겹치는 나머지 열은 _src 접미사)
apply_schema_names <- function(raw, schema) {
  cols  <- schema_columns(schema)
  other <- standardize_names(raw[setdiff(names(raw), cols)])
  clash <- names(other) %in% canon[names(cols)]
  names(other)[clash] <- paste0(names(other)[clash], "_src")
  sch <- raw[unname(cols)]
  names(sch) <- canon[names(cols)]
  cbind(sch, other)
}

# 상태 → 0/1 (event_value = 사건)
code_status <- function(x, event_value) {
  ifelse(is.na(x), NA_real_, as.numeric(as.character(x) == as.character(event_value)))
}

# stage → 1–4 (numeric: 코드 그대로, ajcc_text: "Stage IIIB" → 3)
parse_stage <- function(x, format) {
  if (format == "numeric") {
    v <- suppressWarnings(as.numeric(as.character(x)))
  } else {
    roman <- c(I = 1, II = 2, III = 3, IV = 4)
    m <- str_match(toupper(as.character(x)), "^STAGE\\s+(IV|III|II|I)")[, 2]
    v <- unname(roman[m])
  }
  bad <- unique(as.character(x)[!is.na(x) & (is.na(v) | !v %in% 1:4)])
  if (length(bad)) {
    stop("stage 변환 불가 값 (stage_format = ", format, "): ", paste(bad, collapse = ", "))
  }
  v
}

# 01 시작 시 실행: 스키마 열 존재 + 값 점검. 오류는 중단, 의심 사항은 경고
validate_schema <- function(cancer) {
  schema <- resolve_schema(cancer)
  if (!file.exists(sav_path(cancer))) {
    message(cancer, ": ", sav_path(cancer), " 없음 → 스키마 검증 생략 (GDC 임상 사용)")
    return(invisible(schema))
  }
  raw <- read_sav_raw(cancer)
  check_schema_columns(raw, schema, cancer)
  d <- clean_tcga(apply_schema_names(raw, schema))
  where <- function(k) paste0(" (cfg$schema$", cancer, "$", k, " = '", schema[[k]], "')")
  ids_txt <- function(x) paste(head(x, 10), collapse = ", ")

  # ID
  if (any(is.na(d$sample_id))) stop(cancer, ": ID 결측 ", sum(is.na(d$sample_id)), "건", where("id"))
  dup <- unique(d$sample_id[duplicated(d$sample_id)])
  if (length(dup)) stop(cancer, ": 중복 ID ", ids_txt(dup), where("id"))

  # 생존 상태
  vals <- unique(as.character(na.omit(d$status)))
  ev   <- as.character(schema$event_value)
  if (length(vals) > 2) {
    stop(cancer, ": status 값이 2개 초과 (", paste(vals, collapse = ", "), ")", where("status"))
  }
  if (!ev %in% vals) {
    stop(cancer, ": status에 event_value '", ev, "' 없음 (값: ", paste(vals, collapse = ", "),
         ") → cfg$schema$", cancer, "$event_value 확인")
  }
  status <- code_status(d$status, ev)

  # 생존 시간
  t <- d$surv_time
  if (!is.numeric(t)) stop(cancer, ": 생존 시간이 숫자가 아님", where("time"))
  if (any(t < 0, na.rm = TRUE)) {
    stop(cancer, ": 음수 생존 시간 → ", ids_txt(d$sample_id[which(t < 0)]))
  }
  if (any(t == 0, na.rm = TRUE)) {
    warning(cancer, ": 생존 시간 0인 환자 ", ids_txt(d$sample_id[which(t == 0)]), call. = FALSE)
  }
  if (!identical(is.na(t), is.na(status))) {
    stop(cancer, ": 생존 시간/상태 결측 위치 불일치 → ",
         ids_txt(d$sample_id[is.na(t) != is.na(status)]))
  }
  mx <- max(t, na.rm = TRUE)
  if (schema$time_unit == "days" && mx < 200) {
    warning(cancer, ": time_unit = days인데 최대 ", mx, " → 개월 단위일 수 있음",
            where("time_unit"), call. = FALSE)
  }
  if (schema$time_unit == "months" && mx > 600) {
    warning(cancer, ": time_unit = months인데 최대 ", mx, " → 일 단위일 수 있음",
            where("time_unit"), call. = FALSE)
  }

  # 사건 비율
  rate <- mean(status, na.rm = TRUE)
  if (rate < 0.05 || rate > 0.70) {
    warning(cancer, ": 사건 비율 ", sprintf("%.1f%%", 100 * rate),
            " (5–70% 범위 밖) → event_value 코딩 확인", call. = FALSE)
  }

  # 역코딩 점검: 마지막 접촉일이 없는 환자(보통 사망)는 대부분 사건이어야 함
  if ("last_contact" %in% names(d)) {
    no_lc <- is.na(d$last_contact) & !is.na(status)
    if (sum(no_lc) >= 5 && mean(status[no_lc]) < 0.5) {
      warning(cancer, ": last_contact 결측 ", sum(no_lc), "명 중 사건 ",
              sprintf("%.0f%%", 100 * mean(status[no_lc])),
              " → status 코딩이 반대일 수 있음", where("event_value"), call. = FALSE)
    }
  }

  # stage, age, sex
  parse_stage(d$pathologic_stage, schema$stage_format)
  if (!is.numeric(d$age) || any(d$age < 0 | d$age > 120, na.rm = TRUE)) {
    stop(cancer, ": 나이 값 이상 (숫자 0–120 아님)", where("age"))
  }
  if (length(unique(na.omit(d$gender))) != 2) {
    warning(cancer, ": 성별 값이 2개가 아님 (",
            paste(unique(na.omit(d$gender)), collapse = ", "), ")", call. = FALSE)
  }

  cat(sprintf("스키마 확인: %d명, 생존 %d명 (사건 %d, %.1f%%), 시간 단위 %s, stage %s\n",
              nrow(d), sum(!is.na(status)), sum(status, na.rm = TRUE), 100 * rate,
              schema$time_unit, schema$stage_format))
  invisible(schema)
}

# ---- 임상 데이터 -------------------------------------------------------

# 원본 → 표준 이름 → 정리 → status 0/1, stage 1–4, os_months
load_clinical <- function(cancer) {
  schema <- resolve_schema(cancer)
  if (file.exists(sav_path(cancer))) {
    raw <- read_sav_raw(cancer)
    check_schema_columns(raw, schema, cancer)
    df <- apply_schema_names(raw, schema)
  } else {
    message(cancer, ": ", sav_path(cancer), " 없음 → GDC 임상 데이터 사용")
    df <- load_gdc_clinical(cancer)
    schema$event_value <- 1
    schema$time_unit <- "days"
    schema$stage_format <- "numeric"
  }
  df <- clean_tcga(df)
  df$status <- code_status(df$status, schema$event_value)
  df$pathologic_stage <- parse_stage(df$pathologic_stage, schema$stage_format)
  df$os_months <- if (schema$time_unit == "days") df$surv_time / 30.44 else df$surv_time
  df
}

# .sav가 없는 암종: GDC 임상 → 표준 이름 열 (미검증 경로)
load_gdc_clinical <- function(cancer) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop(cancer, ".sav가 없고 TCGAbiolinks도 설치되지 않음 → Rscript 00_setup.R 실행 후 재시도")
  }
  cl <- TCGAbiolinks::GDCquery_clinic(paste0("TCGA-", cancer), type = "clinical")
  dead <- cl$vital_status == "Dead"
  roman <- c(I = 1, II = 2, III = 3, IV = 4)
  out <- data.frame(
    cl$submitter_id,
    ifelse(dead, cl$days_to_death, cl$days_to_last_follow_up),
    ifelse(is.na(cl$vital_status), NA, as.numeric(dead)),
    unname(roman[str_match(cl$ajcc_pathologic_stage, "Stage (IV|III|II|I)")[, 2]]),
    cl$age_at_index, toupper(cl$gender), cl$days_to_last_follow_up,
    stringsAsFactors = FALSE
  )
  names(out) <- canon[c("id", "time", "status", "stage", "age", "sex", "last_contact")]
  out
}

# 그룹 변수(코드)가 원 변수 > cutoff 와 일치하는지 (라벨 적용 전). 불일치 시 중단
check_label_consistency <- function(df, cancer) {
  rules <- rules_for(cancer)
  res <- list()
  for (v in names(rules)) {
    r <- rules[[v]]
    if (!all(c(v, r$source) %in% names(df))) next
    x   <- df[[v]]
    src <- df[[r$source]]
    ok  <- !is.na(x) & !is.na(src)
    expected <- src[ok] > r$cutoff
    actual   <- as.character(x[ok]) == as.character(r$high)
    bad <- df$sample_id[ok][expected != actual]
    res[[v]] <- data.frame(variable = v, source = r$source, cutoff = r$cutoff,
                           n_compared = sum(ok), n_mismatch = length(bad))
    if (length(bad)) {
      stop(cancer, ": ", v, "이(가) ", r$source, " > ", r$cutoff, " 규칙과 불일치 ",
           length(bad), "명 → ", paste(head(bad, 10), collapse = ", "),
           " (규칙: cfg$label_rules$", cancer, "$", v, ")")
    }
  }
  do.call(rbind, res)
}

# 분석용 변수 생성 (없는 열은 건너뜀 → 다른 암종에도 사용 가능)
recode_clinical <- function(df, cancer) {
  # 평가 불가 등 분석에서 제외할 값 → NA (cfg$analysis_na_values)
  for (v in intersect(names(cfg$analysis_na_values), names(df))) {
    hit <- !is.na(df[[v]]) & as.character(df[[v]]) %in% cfg$analysis_na_values[[v]]
    if (any(hit)) {
      cat(sprintf("  %s: %s → NA %d명\n", v,
                  paste(cfg$analysis_na_values[[v]], collapse = "/"), sum(hit)))
      df[[v]][hit] <- NA
    }
  }

  if ("gender" %in% names(df)) {
    df$gender <- factor(str_to_title(df$gender))   # 기준 수준: cfg$reference_levels
  }
  tnm <- c(pathologic_t = "T", pathologic_n = "N", pathologic_m = "M")
  for (v in intersect(names(tnm), names(df))) df[[v]] <- code_factor(df[[v]], tnm[[v]])

  # stage 그룹 열이 없으면 stage에서 생성 (1–2 → 0, 3–4 → 1)
  if (!cfg$stage_collapsed %in% names(df)) {
    df[[cfg$stage_collapsed]] <- ifelse(is.na(df$pathologic_stage), NA,
                                        as.numeric(df$pathologic_stage >= 3))
  }

  checks <- check_label_consistency(df, cancer)
  df <- apply_value_labels(df, labels_for(cancer))
  df$stage <- df$pathologic_stage

  # YES/NO 열 → factor(No, Yes)
  # (검사 시행 여부 변수도 여기 포함 — 결과 변수로 사용하지 않음)
  for (v in names(df)) {
    x <- df[[v]]
    if (is.character(x) && any(!is.na(x)) && all(toupper(na.omit(x)) %in% c("YES", "NO"))) {
      df[[v]] <- factor(str_to_title(x), levels = c("No", "Yes"))
    }
  }
  if ("tumor_status" %in% names(df)) {
    df$tumor_status <- factor(str_to_sentence(df$tumor_status),
                              levels = c("Tumor free", "With tumor"))
  }

  neo <- if ("neoadjuvant" %in% names(df)) as.character(df$neoadjuvant) else rep(NA, nrow(df))
  neo_yes <- resolve_schema(cancer)$neoadjuvant_yes %||% "Yes"
  df$neoadjuvant_flag <- !is.na(neo) & tolower(neo) == tolower(neo_yes)

  # 검체 제공 기관(tissue source site): TCGA-XX-.... 의 XX (민감도 분석 strata)
  df$tss <- substr(df$sample_id, cfg$tss_barcode_pos[1], cfg$tss_barcode_pos[2])

  df <- apply_reference_levels(df)
  df <- set_var_labels(df)
  attr(df, "label_checks") <- checks
  df
}

# ---- 발현 그룹 ---------------------------------------------------------

# 중앙값 초과 = High (새 유전자에만 사용; 기존 .sav 그룹은 덮어쓰지 않음)
median_split <- function(x) {
  factor(ifelse(is.na(x), NA, ifelse(x > median(x, na.rm = TRUE), group_alt(), group_ref())),
         levels = cfg$group_levels)
}

# 기존 그룹이 median split과 일치하는지 QC (홀수 n의 중앙값 동점 1명 차이 허용)
check_median_split <- function(expr, group, tol = cfg$median_tie_tolerance) {
  ok <- !is.na(expr) & !is.na(group)
  n_diff <- sum(as.character(median_split(expr[ok])) != as.character(group[ok]))
  data.frame(n = sum(ok),
             n_low = sum(group[ok] == group_ref()), n_high = sum(group[ok] == group_alt()),
             median = median(expr[ok]),
             n_diff_from_median_split = n_diff,
             median_split_pass = n_diff <= tol)
}

# <CANCER>_<split>_<GENE>.csv → sample_id, expression, group, csv_status(0/1)
read_gene_file <- function(path) {
  parts <- str_match(basename(path), cfg$gene_pattern)
  if (is.na(parts[1])) stop("파일명 규칙 불일치: ", basename(path))

  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("Patient", "Expression", "Group")
  if (!all(need %in% names(d))) {
    stop(basename(path), ": 필요한 열 없음 (", paste(setdiff(need, names(d)), collapse = ", "), ")")
  }

  status <- rep(NA_real_, nrow(d))
  if ("Status" %in% names(d)) {
    s <- str_to_title(str_squish(d$Status))
    bad <- setdiff(unique(na.omit(s)), c("Alive", "Dead", ""))
    if (length(bad)) stop(basename(path), ": 알 수 없는 Status 값 ", paste(bad, collapse = ", "))
    status <- c(Alive = 0, Dead = 1)[s]
  }

  out <- data.frame(
    sample_id  = str_sub(str_squish(d$Patient), 1, 12),
    expression = as.numeric(d$Expression),
    group      = as_group_factor(d$Group),
    csv_status = unname(status),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(out$sample_id)) {
    stop(basename(path), ": 중복 환자 ", paste(unique(out$sample_id[duplicated(out$sample_id)]), collapse = ", "))
  }
  list(cancer = toupper(parts[2]), split = parts[3], gene = parts[4],
       key = gene_key(parts[4]), file = basename(path), data = out)
}

# 02a 발현 행렬(log2(TPM+1))에서 유전자 추출 → TPM으로 되돌리고 median split
# (다른 발현 열과 같은 선형 척도 유지). ids: median 계산에 쓸 환자
extract_gdc_genes <- function(cancer, genes, ids) {
  f <- file.path(cfg$clinical_dir, paste0("TCGA_", cancer, "_RNAseq_Expression.csv"))
  if (!file.exists(f)) stop(f, " 없음 → 02a_download_tcga.R 먼저 실행")
  m <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
  miss <- setdiff(genes, m$Gene)
  if (length(miss)) warning(cancer, ": 발현 행렬에 없는 유전자 ", paste(miss, collapse = ", "))
  m <- m[m$Gene %in% genes, ]

  out <- data.frame(sample_id = ids)
  for (i in seq_len(nrow(m))) {
    key  <- gene_key(m$Gene[i])
    expr <- 2^unlist(m[i, -1]) - 1
    x    <- unname(expr[ids])
    out[[paste0(key, "_expression")]] <- x
    out[[paste0(key, "_group")]]      <- median_split(x)
  }
  out
}
