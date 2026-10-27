# =============================================================
# 06_external_geo.R
# 주 가설 유전자(cfg$primary_gene)의 GEO 외부 검증 (데이터셋: cfg$geo_datasets)
#   1. GEOquery series matrix + GPL 주석 (cfg$geo_dir 캐시), log2 척도 확인
#   2. probe → 유전자: 매칭 probe 전부 QC, cfg$geo_probe_rule로 하나 선택
#   3. characteristics_ch1 파싱 → config fields 매핑, 비종양 샘플 제외
#   4. 종양 샘플 내 median split (cfg$group_levels, 첫 번째 = 기준)
#   5. 종점별 (첫 번째 = 주 종점): KM + log-rank, 단변량 Cox, 다변량 Cox (EPV),
#      다변량 + extra_covariate (MMR), 하위군 (pMMR), cox.zph
#   6. TCGA 결과와 방향 비교 (validates 코호트의 survival_summary_raw.csv)
# 사용: Rscript 06_external_geo.R [GSE39582 ...]   (인자 없으면 cfg$geo_datasets 전체)
# 결과 (파일명 앞에 <GSE>_): output/tables/GEO/<GSE>/ (probe_QC, pheno_QC, survival_summary, km_summary, cox_<EP>, ph_tests)
#       output/figures/GEO/<GSE>/ (km_<EP>, km_<EP>_<하위군>, forest, forest_multi_<EP>)
# =============================================================

source("config.R")
source("R/utils.R")
source("R/survival.R")
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

options(timeout = max(1200, getOption("timeout")))    # 큰 series matrix 다운로드
Sys.setenv(VROOM_CONNECTION_SIZE = 131072 * 64)      # 긴 series matrix 줄 읽기 오류 방지

geo_where <- function(gse, key) paste0(" → cfg$geo_datasets$", gse, "$", key)
ids_txt <- function(x) paste(head(x, 10), collapse = ", ")

# ---- 1. 다운로드 / 캐시 ------------------------------------------------

load_geo_eset <- function(gse, g) {
  dir.create(cfg$geo_dir, recursive = TRUE, showWarnings = FALSE)
  rds <- file.path(cfg$geo_dir, paste0(gse, "_", g$platform, "_eset.rds"))
  if (file.exists(rds)) {
    cat("[캐시] ", rds, "\n", sep = "")
    eset <- readRDS(rds)
  } else {
    cat("[다운로드] ", gse, " → ", cfg$geo_dir, " (받아 둔 series matrix/GPL 파일이 있으면 재사용)\n", sep = "")
    lst <- getGEO(gse, GSEMatrix = TRUE, destdir = cfg$geo_dir, getGPL = TRUE)
    idx <- if (length(lst) == 1) 1 else grep(g$platform, names(lst), fixed = TRUE)
    if (length(idx) != 1) {
      stop(gse, ": platform '", g$platform, "'의 series matrix를 하나로 특정할 수 없음 (",
           paste(names(lst), collapse = ", "), ")", geo_where(gse, "platform"))
    }
    eset <- lst[[idx]]
    saveRDS(eset, rds)
  }
  if (annotation(eset) != g$platform) {
    stop(gse, ": platform이 ", annotation(eset), " (설정 ", g$platform, ")", geo_where(gse, "platform"))
  }
  eset
}

# log2 척도 확인 (max > 100 → log2 변환)
prepare_expression <- function(eset, gse) {
  ex <- exprs(eset)
  rng <- range(ex, na.rm = TRUE)
  cat(sprintf("발현 행렬: %d probes × %d samples, 범위 %.2f–%.2f, NA %d\n",
              nrow(ex), ncol(ex), rng[1], rng[2], sum(is.na(ex))))
  if (rng[2] > 100) {
    if (rng[1] < 0) stop(gse, ": 선형 척도로 보이는데 음수 값이 있어 log2 변환 불가")
    off <- if (rng[1] < 1) 1 else 0
    ex <- log2(ex + off)
    cat(sprintf("  max > 100 → log2%s 변환, 새 범위 %.2f–%.2f\n",
                if (off) "(x + 1)" else "(x)", min(ex, na.rm = TRUE), max(ex, na.rm = TRUE)))
  } else {
    cat("  이미 log2 척도 (max ≤ 100) → 변환 없음\n")
  }
  ex
}

# ---- 2. probe → 유전자 -------------------------------------------------

# 유전자에 매칭되는 모든 probe의 QC (tumor 샘플 기준) + cfg$geo_probe_rule로 선택
probe_qc <- function(ex, fd, gene, g, gse) {
  sc <- g$symbol_col
  if (!sc %in% names(fd)) {
    stop(gse, ": GPL 주석에 '", sc, "' 열 없음 (열: ", paste(names(fd), collapse = ", "), ")",
         geo_where(gse, "symbol_col"))
  }
  if (!identical(rownames(ex), rownames(fd))) stop(gse, ": 발현 행렬과 GPL 주석의 probe 순서 불일치")
  syms <- lapply(strsplit(as.character(fd[[sc]]), "///", fixed = TRUE), trimws)
  hit  <- vapply(syms, function(s) gene %in% s, logical(1))
  probes <- rownames(fd)[hit]
  if (!length(probes)) stop(gse, ": '", sc, "'에서 ", gene, " probe를 찾지 못함")

  sub <- ex[probes, , drop = FALSE]
  qc <- data.frame(
    probe_id     = probes,
    symbol       = as.character(fd[probes, sc]),
    multi_mapped = lengths(syms[hit]) > 1,
    n            = rowSums(!is.na(sub)),
    mean         = rowMeans(sub, na.rm = TRUE),
    median       = apply(sub, 1, median, na.rm = TRUE),
    IQR          = apply(sub, 1, IQR, na.rm = TRUE),
    sd           = apply(sub, 1, sd, na.rm = TRUE),
    row.names = NULL, stringsAsFactors = FALSE
  )
  rules <- c(max_mean = "mean", max_iqr = "IQR", max_sd = "sd")
  rule  <- cfg$geo_probe_rule
  if (!rule %in% names(rules)) {
    stop("cfg$geo_probe_rule = '", rule, "' 지원 안 함 (", paste(names(rules), collapse = ", "), ")")
  }
  qc$chosen <- seq_len(nrow(qc)) == which.max(qc[[rules[[rule]]]])
  if (length(probes) > 1) {
    r <- cor(t(sub), use = "pairwise.complete.obs")
    for (p in probes) qc[[paste0("r_", p)]] <- r[probes, p]
  }
  qc$rule <- rule
  qc
}

# ---- 3. phenotype ------------------------------------------------------

# characteristics_ch1.* ("key: value") → 샘플 × key 표 (원본 key 이름 그대로)
parse_characteristics <- function(pd, gse) {
  cols <- grep("^characteristics_ch1", names(pd), value = TRUE)
  if (!length(cols)) stop(gse, ": characteristics_ch1 열 없음")
  long <- do.call(rbind, lapply(cols, function(cl) {
    x  <- as.character(pd[[cl]])
    ok <- !is.na(x) & grepl(":", x, fixed = TRUE)
    data.frame(sample_id = rownames(pd)[ok],
               key   = trimws(sub(":.*$", "", x[ok])),
               value = trimws(sub("^[^:]*:", "", x[ok])),
               stringsAsFactors = FALSE)
  }))
  dup <- duplicated(long[, c("sample_id", "key")])
  if (any(dup)) {
    stop(gse, ": 한 샘플에 같은 key가 여러 번 → ", ids_txt(unique(paste(long$sample_id[dup], long$key[dup]))))
  }
  out <- data.frame(sample_id = rownames(pd), stringsAsFactors = FALSE)
  for (k in unique(long$key)) {
    v <- long[long$key == k, ]
    out[[k]] <- v$value[match(out$sample_id, v$sample_id)]
  }
  out
}

is_geo_missing <- function(x) is.na(x) | trimws(x) %in% cfg$geo_na_values

print_raw_keys <- function(raw) {
  keys <- setdiff(names(raw), "sample_id")
  tab <- data.frame(
    key = keys,
    n_present = vapply(keys, function(k) sum(!is_geo_missing(raw[[k]])), 0L),
    n_unique  = vapply(keys, function(k) length(unique(raw[[k]][!is_geo_missing(raw[[k]])])), 0L),
    examples  = vapply(keys, function(k) {
      u <- unique(raw[[k]][!is_geo_missing(raw[[k]])])
      substr(paste(head(u, 5), collapse = " | "), 1, 60)
    }, ""),
    row.names = NULL
  )
  cat("\n[characteristics_ch1 원본 key] (", nrow(raw), " samples; config fields에 이 key를 지정)\n", sep = "")
  print(tab, right = FALSE, row.names = FALSE)
  invisible(tab)
}

# config fields (표준 이름 → 원본 key) 적용. 결측 표기 → NA
map_fields <- function(raw, g, gse) {
  f <- g$fields
  ep_fields <- unlist(lapply(g$endpoints, function(e) c(e$time, e$event)))
  req <- unique(c("sample_type", "age", "gender", "stage", ep_fields,
                  g$extra_covariate, g$subgroup$var))
  miss <- setdiff(req, names(f))
  if (length(miss)) stop(gse, ": fields에 필수 항목 없음 (", paste(miss, collapse = ", "), ")", geo_where(gse, "fields"))
  for (k in names(f)) {
    if (!f[[k]] %in% names(raw)) {
      near <- agrep(f[[k]], names(raw), max.distance = 0.3, value = TRUE, ignore.case = TRUE)
      stop(gse, ": 원본 key '", f[[k]], "' 없음", geo_where(gse, paste0("fields[\"", k, "\"]")),
           if (length(near)) paste0(". 비슷한 key: ", paste(head(near, 5), collapse = ", ")) else "")
    }
  }
  d <- data.frame(sample_id = raw$sample_id, stringsAsFactors = FALSE)
  for (k in names(f)) {
    x <- trimws(raw[[f[[k]]]])
    x[is_geo_missing(x)] <- NA
    d[[k]] <- x
  }
  d
}

num_field <- function(d, k, g, gse) {
  x <- d[[k]]
  v <- suppressWarnings(as.numeric(x))
  bad <- !is.na(x) & is.na(v)
  if (any(bad)) {
    stop(gse, ": ", k, " (원본 key '", g$fields[[k]], "')에 숫자가 아닌 값 ",
         paste(unique(x[bad]), collapse = ", "), " → ", ids_txt(d$sample_id[bad]),
         " (결측 표기면 cfg$geo_na_values에 추가)")
  }
  v
}

# 비종양 제외 + 공변량 / 종점 변수 생성
recode_geo <- function(d, g, gse) {
  if (any(is.na(d$sample_type))) {
    stop(gse, ": sample_type 결측 → ", ids_txt(d$sample_id[is.na(d$sample_type)]),
         geo_where(gse, "fields[\"sample_type\"]"))
  }
  cat("\n샘플 유형:\n")
  print(table(d$sample_type))
  excl <- d$sample_type %in% g$exclude_sample_type
  cat(sprintf("비종양 제외 (%s): %d → 종양 %d\n",
              paste(g$exclude_sample_type, collapse = ", "), sum(excl), sum(!excl)))
  d <- d[!excl, ]
  if (!nrow(d)) stop(gse, ": 종양 샘플 없음", geo_where(gse, "exclude_sample_type"))

  # 나이
  d$age <- num_field(d, "age", g, gse)
  if (any(d$age < 0 | d$age > 120, na.rm = TRUE)) stop(gse, ": 나이 값 이상 (0–120 아님)")

  # 성별 (기준 수준: cfg$reference_levels)
  d$gender <- factor(str_to_title(d$gender))
  if (nlevels(d$gender) != 2) warning(gse, ": 성별 값이 2개가 아님 (", paste(levels(d$gender), collapse = ", "), ")", call. = FALSE)

  # stage: stage_na_values → NA (명시), 나머지는 1–4 코드 (TCGA와 같은 라벨)
  s <- d$stage
  to_na <- !is.na(s) & s %in% g$stage_na_values
  if (any(to_na)) {
    cat(sprintf("stage %s → NA %d명 (cfg$geo_datasets$%s$stage_na_values): %s\n",
                paste(g$stage_na_values, collapse = "/"), sum(to_na), gse, ids_txt(d$sample_id[to_na])))
    s[to_na] <- NA
  }
  stage_labels <- labels_for(g$validates)
  bad <- !is.na(s) & !s %in% names(stage_labels$pathologic_stage)
  if (any(bad)) {
    stop(gse, ": stage 값 ", paste(unique(s[bad]), collapse = ", "), " 해석 불가 → ", ids_txt(d$sample_id[bad]),
         geo_where(gse, "stage_na_values"))
  }
  d$stage_num <- as.numeric(s)
  d$pathologic_stage <- d$stage_num
  d$pathologic_stage_12_34 <- ifelse(is.na(d$stage_num), NA, as.numeric(d$stage_num >= 3))
  d <- apply_value_labels(d, stage_labels[intersect(c("pathologic_stage", "pathologic_stage_12_34"),
                                                    names(stage_labels))])
  d$stage <- d$pathologic_stage

  # 나머지 범주형 필드 → factor (factor_levels가 있으면 그 순서, 첫 번째 = 기준)
  ep_fields <- unlist(lapply(g$endpoints, function(e) c(e$time, e$event)))
  other <- setdiff(names(g$fields), c("sample_type", "age", "gender", "stage", ep_fields))
  for (k in other) {
    lv <- g$factor_levels[[k]]
    if (!is.null(lv)) {
      bad <- !is.na(d[[k]]) & !d[[k]] %in% lv
      if (any(bad)) {
        stop(gse, ": ", k, " 값 ", paste(unique(d[[k]][bad]), collapse = ", "), "이(가) factor_levels에 없음 → ",
             ids_txt(d$sample_id[bad]), geo_where(gse, paste0("factor_levels$", k)))
      }
      d[[k]] <- factor(d[[k]], levels = lv)
    } else {
      d[[k]] <- factor(d[[k]])
    }
  }

  # 종점: <EP>_months, <EP>_event (0/1)
  for (ep in names(g$endpoints)) {
    e <- g$endpoints[[ep]]
    tm <- num_field(d, e$time, g, gse)
    if (!e$time_unit %in% c("days", "months")) stop(gse, " ", ep, ": time_unit은 days/months", geo_where(gse, paste0("endpoints$", ep)))
    if (e$time_unit == "days") tm <- tm / 30.44
    if (any(tm < 0, na.rm = TRUE)) stop(gse, " ", ep, ": 음수 시간 → ", ids_txt(d$sample_id[which(tm < 0)]))
    ev_raw <- d[[e$event]]
    vals <- unique(na.omit(ev_raw))
    if (length(vals) > 2 || !as.character(e$event_value) %in% vals) {
      stop(gse, " ", ep, ": event 값 (", paste(vals, collapse = ", "), ")이 event_value ",
           e$event_value, "와 맞지 않음", geo_where(gse, paste0("endpoints$", ep, "$event_value")))
    }
    d[[paste0(ep, "_months")]] <- tm
    d[[paste0(ep, "_event")]]  <- code_status(ev_raw, e$event_value)
  }

  d <- apply_reference_levels(d)
  set_var_labels(d)
}

# 필드별 결측 / 값 요약 (종양 샘플)
pheno_qc <- function(d, g) {
  rows <- lapply(names(g$fields), function(k) {
    x <- d[[k]]
    if (k %in% c(unlist(lapply(g$endpoints, `[[`, "time")))) {
      ep <- names(g$endpoints)[vapply(g$endpoints, function(e) e$time == k, TRUE)]
      x  <- d[[paste0(ep, "_months")]]
    }
    if (k %in% c(unlist(lapply(g$endpoints, `[[`, "event")))) {
      ep <- names(g$endpoints)[vapply(g$endpoints, function(e) e$event == k, TRUE)]
      x  <- d[[paste0(ep, "_event")]]
    }
    vals <- if (is.numeric(x) && length(unique(na.omit(x))) > 5) {
      sprintf("median %.1f (range %.1f–%.1f)", median(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
    } else {
      t <- table(x)
      paste(paste0(names(t), " = ", as.integer(t)), collapse = "; ")
    }
    data.frame(field = k, raw_key = g$fields[[k]], n = sum(!is.na(x)), n_missing = sum(is.na(x)),
               values = vals, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- 5. 생존분석 -------------------------------------------------------

model_label_extra <- function(g) var_labels[[g$extra_covariate]] %||% g$extra_covariate

# 한 종점 × 한 집단 (All 또는 하위군)
analyze_endpoint <- function(d, ep, g, gse, gene, pop = "All") {
  e     <- g$endpoints[[ep]]
  time  <- paste0(ep, "_months")
  event <- paste0(ep, "_event")
  out   <- list(summary = NULL, km = NULL, cox = NULL, zph = NULL, files = character(), forest = NULL)
  tag   <- if (pop == "All") ep else paste0(ep, "_", pop)

  has_ep <- !is.na(d[[time]]) & !is.na(d[[event]])
  one_side <- is.na(d[[time]]) != is.na(d[[event]])
  if (any(one_side)) {
    cat(sprintf("  [%s] 시간/사건 중 하나만 결측 %d명 → 이 종점에서 제외: %s\n",
                tag, sum(one_side), ids_txt(d$sample_id[one_side])))
  }
  de <- d[has_ep & !is.na(d$group), ]
  if (!is.null(e$exclude_stage)) {
    ex_st <- de$stage_num %in% e$exclude_stage
    cat(sprintf("  [%s] stage %s 제외 %d명 (cfg$geo_datasets$%s$endpoints$%s$exclude_stage)\n",
                tag, paste(e$exclude_stage, collapse = "/"), sum(ex_st), gse, ep))
    de <- de[!ex_st, ]
  }
  de <- add_scaled_terms(de)
  cnt <- table(factor(de$group, levels = cfg$group_levels))
  n_ev <- sum(de[[event]])
  cat(sprintf("  [%s] n = %d (%s %d / %s %d), 사건 %d\n", tag, nrow(de),
              group_ref(), cnt[[group_ref()]], group_alt(), cnt[[group_alt()]], n_ev))
  if (n_ev < cfg$min_events || any(cnt == 0)) {
    warning(gse, " ", tag, ": 사건 ", n_ev, " 또는 빈 그룹 → 건너뜀", call. = FALSE)
    return(out)
  }

  # ---- KM ----
  f  <- as.formula(paste0("Surv(", time, ", ", event, ") ~ group"))
  km <- surv_fit(f, data = de)
  p_lr <- survdiff(f, data = de)$pvalue
  out$km <- data.frame(endpoint = ep, population = pop, km_stats(km, cfg$km_times), logrank_p = p_lr)
  title <- paste0(gene, " expression — ", g$label, if (pop != "All") paste0(", ", pop))
  out$files <- c(out$files, suppressMessages(save_fig(
    km_plot(km, de, gene, NULL, p_lr, title = title, ylab = paste(e$label, "probability")),
    "GEO", paste0("km_", tag), width = 7, height = 6, subdir = gse, prefix = gse)))

  # ---- Cox 모형들 ----
  cov_terms <- intersect(unname(sapply(cfg$cox_covariates, model_term)), names(de))
  labs <- term_labels(c("group", "expr_log2", cov_terms, g$extra_covariate, cfg$stage_collapsed), gene)
  rows <- list(); cox_tabs <- list(); zph <- list()

  add_model <- function(name, m, key_term = "group", keep_table = FALSE) {
    fit <- fit_cox(m$d, m$terms, time = time, event = event)
    tab <- cox_rows(fit, m$d, m$terms, labs[m$terms], event = event)
    z   <- zph_rows(fit, gene, name)
    k   <- tab[tab$term == key_term & !tab$reference, ]
    rows[[length(rows) + 1]] <<- data.frame(
      endpoint = ep, population = pop, model = name, term = key_term,
      n = nrow(m$d), events = m$events, dropped = m$dropped, params = m$k, epv = m$epv,
      stage_levels = m$stage_levels %||% NA, exploratory = m$exploratory %||% FALSE,
      hr = k$hr, lo = k$lo, hi = k$hi, p = k$p,
      ph_p_term = z$p[z$term == key_term], ph_p_global = z$p[z$term == "GLOBAL"],
      stringsAsFactors = FALSE)
    cox_tabs[[name]] <<- data.frame(Model = name, format_cox(tab), check.names = FALSE)
    zph[[name]] <<- data.frame(endpoint = ep, population = pop, z)
    if (keep_table) out$forest <<- list(tab = tab, name = name, m = m)
    fit
  }
  simple <- function(terms) {
    dm <- complete_cases(de, terms)
    list(terms = terms, d = dm, events = sum(dm[[event]]), k = n_params(dm, terms),
         epv = sum(dm[[event]]) / n_params(dm, terms), dropped = nrow(de) - nrow(dm),
         stage_levels = NA, exploratory = FALSE)
  }

  add_model("Univariable", simple("group"))
  add_model("Univariable, continuous (per 1 log2 unit)", simple("expr_log2"), key_term = "expr_log2")

  m <- build_multi(de, c("group", cov_terms), event = event)
  add_model("Multivariable", m, keep_table = is.null(g$extra_covariate) || pop != "All")
  cat(sprintf("  [%s] 다변량: n = %d (결측 제외 %d), 사건 %d, 모수 %d, EPV %.1f, stage = %s%s\n",
              tag, nrow(m$d), m$dropped, m$events, m$k, m$epv, m$stage_levels,
              if (m$exploratory) " → 탐색적(EPV 부족)" else ""))

  # ---- + extra covariate (MMR): 같은 환자에서 보정 전/후 비교 + 상호작용 ----
  x <- g$extra_covariate
  if (!is.null(x) && pop == "All" && x %in% names(de)) {
    xl  <- model_label_extra(g)
    mk  <- build_multi(de[!is.na(de[[x]]), ], c("group", cov_terms), event = event)
    add_model(paste0("Multivariable (", xl, " known)"), mk)
    mx  <- build_multi(de, c("group", cov_terms, x), event = event)
    fx  <- add_model(paste0("Multivariable + ", xl), mx, keep_table = TRUE)
    fi  <- fit_cox(mx$d, c(mx$terms, paste0("group:", x)), time = time, event = event)
    p_int <- anova(fx, fi)[["Pr(>|Chi|)"]][2]
    rows[[length(rows)]]$p_interaction <- p_int
    cat(sprintf("  [%s] + %s: n = %d, 사건 %d, EPV %.1f; group × %s 상호작용 p = %s\n",
                tag, xl, nrow(mx$d), mx$events, mx$epv, x, fmt_p(p_int)))
  }

  # ---- 시간 분할 Cox (04와 같은 규칙: 유전자 그룹 항이 PH 위반일 때만, 다변량 모형 기준) ----
  zz <- bind_rows(zph)
  if (any(zz$flag & zz$term == "group")) {
    cut <- cfg$ph_split_months
    tsc <- time_split_cox(m$d, m$terms, cut, time = time, event = event)
    cat("  [", tag, "] ", tsc$msg, "\n", sep = "")
    for (w in 1:2) {
      sfx <- c("early", "late")[w]
      rows[[length(rows) + 1]] <- data.frame(
        endpoint = ep, population = pop,
        model = if (w == 1) sprintf("Multivariable, 0–%d months", cut) else sprintf("Multivariable, >%d months", cut),
        term = "group",
        n = if (w == 1) nrow(m$d) else sum(m$d[[time]] > cut), events = sum(tsc$events[w, ]),
        dropped = m$dropped, params = m$k, epv = m$epv, stage_levels = m$stage_levels,
        exploratory = m$exploratory,
        hr = tsc$est[[paste0("hr_", sfx)]], lo = tsc$est[[paste0("lo_", sfx)]],
        hi = tsc$est[[paste0("hi_", sfx)]], p = tsc$est[[paste0("p_", sfx)]],
        ph_p_term = NA, ph_p_global = NA, stringsAsFactors = FALSE)
    }
  }

  # ---- RMST 차이 (High − Low, 비례위험 가정 불필요) ----
  rr <- rmst_diff(de, time = time, event = event)
  if (nzchar(rr$note)) warning(gse, " ", tag, ": ", rr$note, call. = FALSE)
  cat(sprintf("  [%s] RMST 차이 (%s − %s, tau %g개월): %.2f (%.2f, %.2f), p %s\n", tag, group_alt(), group_ref(),
              rr$est$rmst_tau, rr$est$rmst_diff, rr$est$rmst_lo, rr$est$rmst_hi, fmt_p(rr$est$rmst_p)))

  out$summary <- cbind(bind_rows(rows), rr$est)
  out$summary$logrank_p <- p_lr
  out$cox <- bind_rows(cox_tabs)
  out$zph <- bind_rows(zph)
  out
}

# TCGA 결과 (04의 survival_summary_raw.csv) — 방향 비교용
tcga_reference <- function(cancer, gene) {
  f <- file.path(cfg$output_dir, "tables", cancer, out_file(cancer, "survival_summary_raw.csv"))
  if (!file.exists(f)) {
    warning(f, " 없음 → TCGA 방향 비교 생략 (04_survival.R 먼저 실행)", call. = FALSE)
    return(NULL)
  }
  s <- read.csv(f, stringsAsFactors = FALSE)
  r <- s[s$gene == gene, ]
  if (!nrow(r) || is.na(r$multi_hr)) return(NULL)
  r
}

# ---- 데이터셋 1개 ------------------------------------------------------

run_geo_validation <- function(gse, cfg_geo) {
  g    <- cfg_geo
  gene <- toupper(cfg$primary_gene)
  created <- character()
  cat("\n==========", gse, "(", g$platform, ", 검증 대상 TCGA-", g$validates, ") ==========\n", sep = "")

  eset <- load_geo_eset(gse, g)
  ex   <- prepare_expression(eset, gse)
  pd   <- pData(eset)
  cat(sprintf("phenotype: %d samples × %d columns\n", nrow(pd), ncol(pd)))

  # ---- phenotype ----
  raw <- parse_characteristics(pd, gse)
  raw_keys <- print_raw_keys(raw)
  d <- recode_geo(map_fields(raw, g, gse), g, gse)

  # ---- probe 선택 (종양 샘플 기준) ----
  pq <- probe_qc(ex[, d$sample_id, drop = FALSE], fData(eset), gene, g, gse)
  cat("\n[", gene, " probe] ", nrow(pq), "개, 선택 규칙 cfg$geo_probe_rule = ", cfg$geo_probe_rule, "\n", sep = "")
  print(transform(pq, mean = round(mean, 3), median = round(median, 3), IQR = round(IQR, 3), sd = round(sd, 3)),
        row.names = FALSE, digits = 3)
  created <- c(created, save_table(pq, "GEO", "probe_QC.csv", subdir = gse, prefix = gse))
  probe <- pq$probe_id[pq$chosen]

  # ---- median split (종양 샘플 전체) ----
  d$expr_log2 <- unname(ex[probe, d$sample_id])
  d$group <- median_split(d$expr_log2)
  cat(sprintf("\nmedian split (%s, 종양 %d명): median %.3f → %s %d / %s %d, 발현 결측 %d\n",
              probe, sum(!is.na(d$expr_log2)), median(d$expr_log2, na.rm = TRUE),
              group_ref(), sum(d$group == group_ref(), na.rm = TRUE),
              group_alt(), sum(d$group == group_alt(), na.rm = TRUE), sum(is.na(d$expr_log2))))

  # ---- pheno QC ----
  pq_tab <- pheno_qc(d, g)
  for (ep in names(g$endpoints)) {
    ok <- !is.na(d[[paste0(ep, "_months")]]) & !is.na(d[[paste0(ep, "_event")]])
    pq_tab <- rbind(pq_tab, data.frame(
      field = paste0(ep, " (usable)"), raw_key = "", n = sum(ok), n_missing = sum(!ok),
      values = paste0("events = ", sum(d[[paste0(ep, "_event")]][ok])), stringsAsFactors = FALSE))
  }
  cat("\n[phenotype QC] 종양", nrow(d), "명\n")
  print(pq_tab, right = FALSE, row.names = FALSE)
  created <- c(created, save_table(pq_tab, "GEO", "pheno_QC.csv", subdir = gse, prefix = gse))

  # ---- 생존분석: 종점 × (전체, 하위군) ----
  pops <- list(All = d)
  sg <- g$subgroup
  if (!is.null(sg)) pops[[sg$level]] <- d[!is.na(d[[sg$var]]) & d[[sg$var]] == sg$level, ]
  res <- list()
  for (ep in names(g$endpoints)) {
    for (pop in names(pops)) {
      cat("\n--", ep, "/", pop, "--\n")
      r <- analyze_endpoint(pops[[pop]], ep, g, gse, gene, pop)
      res[[paste(ep, pop)]] <- r
      created <- c(created, r$files)
    }
  }

  # ---- TCGA 방향 비교 ----
  tcga <- tcga_reference(g$validates, gene)
  s <- bind_rows(lapply(res, `[[`, "summary"))
  if (!"p_interaction" %in% names(s)) s$p_interaction <- NA
  s$tcga_multi_hr <- if (is.null(tcga)) NA else tcga$multi_hr
  s$direction_match <- if (is.null(tcga)) NA else sign(log(s$hr)) == sign(log(tcga$multi_hr))
  created <- c(created, save_table(s, "GEO", "survival_summary_raw.csv", subdir = gse, prefix = gse))

  ep_lab <- vapply(g$endpoints, `[[`, "", "label")
  tab <- data.frame(
    Endpoint = unname(ep_lab[s$endpoint]), Population = s$population, Model = s$model,
    N = s$n, Events = s$events,
    `HR (95% CI)` = fmt_hr(s$hr, s$lo, s$hi), p = fmt_p(s$p),
    EPV = sprintf("%.1f", s$epv), `Stage variable` = s$stage_levels,
    `PH p (gene term)` = fmt_p(s$ph_p_term), `Interaction p` = fmt_p(s$p_interaction),
    `Log-rank p` = ifelse(duplicated(paste(s$endpoint, s$population)), "", fmt_p(s$logrank_p)),
    check.names = FALSE)
  first <- !duplicated(paste(s$endpoint, s$population))
  tab[[paste0("RMST difference at ", cfg$rmst_tau, " mo, months (95% CI)")]] <-
    ifelse(first, sprintf("%.2f (%.2f, %.2f)", s$rmst_diff, s$rmst_lo, s$rmst_hi), "")
  tab[["RMST p"]] <- ifelse(first, fmt_p(s$rmst_p), "")
  tab[[paste0("Direction as TCGA-", g$validates)]] <- ifelse(is.na(s$direction_match), "",
                                                             ifelse(s$direction_match, "Yes", "No"))
  tab[[paste0("Exploratory (EPV < ", cfg$epv_min, ")")]] <- s$exploratory
  created <- c(created, save_df_table(
    tab, "GEO", "survival_summary", landscape = TRUE, font_size = 7, subdir = gse, prefix = gse,
    caption = paste0(
      "External validation of ", gene, " expression (", group_contrast_label(), ", median split of probe ", probe,
      " within tumors) in ", g$label, ". HR from Cox regression; multivariable models adjusted for ", adjust_text(),
      if (!is.null(g$extra_covariate)) paste0("; '+ ", model_label_extra(g), "' additionally adjusts for ",
                                              tolower(model_label_extra(g)), " and '(", model_label_extra(g),
                                              " known)' is the model without it on the same patients") else "",
      if (!is.null(g$extra_covariate)) paste0(". Interaction p = likelihood-ratio test for ", gene, " group × ",
                                              tolower(model_label_extra(g))) else "",
      ". Time-split HRs (0–", cfg$ph_split_months, " / >", cfg$ph_split_months,
      " months, adjusted) are shown only when the gene term violated proportional hazards. RMST difference = ",
      group_alt(), " − ", group_ref(), " restricted mean survival time up to ", cfg$rmst_tau, " months. ",
      paste(vapply(names(g$endpoints), function(ep) {
        e <- g$endpoints[[ep]]
        if (is.null(e$exclude_stage)) "" else paste0(e$label, " excludes stage ", paste(e$exclude_stage, collapse = "/"), ". ")
      }, ""), collapse = ""),
      if (!is.null(tcga)) sprintf("Direction compared with the TCGA-%s multivariable HR %.2f.", g$validates, tcga$multi_hr) else "")))

  # ---- KM 요약 ----
  km <- bind_rows(lapply(res, `[[`, "km"))
  km_out <- data.frame(Endpoint = unname(ep_lab[km$endpoint]), Population = km$population,
                       Group = km$group, N = km$n, Events = km$events,
                       `Median, months (95% CI)` = fmt_median(km$median, km$median_lo, km$median_hi),
                       check.names = FALSE)
  for (tm in cfg$km_times) {
    col <- paste0("surv_", tm / 12, "y")
    km_out[[paste0(tm / 12, "-year (95% CI)")]] <- fmt_surv(km[[col]], km[[paste0(col, "_lo")]], km[[paste0(col, "_hi")]])
  }
  km_out$`Log-rank p` <- ifelse(duplicated(paste(km$endpoint, km$population)), "", fmt_p(km$logrank_p))
  created <- c(created, save_df_table(km_out, "GEO", "km_summary", subdir = gse, prefix = gse,
                                      caption = paste0("Kaplan–Meier estimates by ", gene, " expression, ",
                                                       g$label, ". NR = not reached.")))

  # ---- Cox 표 (종점 × 집단) ----
  for (k in names(res)) {
    r <- res[[k]]
    if (is.null(r$cox)) next
    ep  <- sub(" .*$", "", k); pop <- sub("^\\S+ ", "", k)
    stem <- paste0("cox_", ep, if (pop != "All") paste0("_", pop))
    created <- c(created, save_df_table(
      r$cox, "GEO", stem, subdir = gse, prefix = gse,
      caption = paste0("Cox regression for ", tolower(ep_lab[[ep]]), ": ", gene, " expression, ", g$label,
                       if (pop != "All") paste0(", ", pop, " only") else "")))
  }
  zph <- bind_rows(lapply(res, `[[`, "zph"))
  created <- c(created, save_table(transform(zph, p = signif(p, 4), chisq = round(chisq, 3)),
                                   "GEO", "ph_tests.csv", subdir = gse, prefix = gse))
  if (any(zph$flag)) {
    fl <- zph[zph$flag, ]
    cat("\ncox.zph p < 0.05:", paste(sprintf("%s/%s %s: %s", fl$endpoint, fl$population, fl$model, fl$term),
                                    collapse = "; "), "\n")
  }

  # ---- forest: 모형별 유전자 HR ----
  fs <- s[s$term == "group", ]
  ft <- data.frame(
    label = sprintf("%s, %s: %s (n = %d, events = %d)", fs$endpoint, fs$population, fs$model, fs$n, fs$events),
    hr = fs$hr, lo = fs$lo, hi = fs$hi, p = fs$p, reference = FALSE,
    shape = ifelse(fs$exploratory, "exploratory", "est"))
  created <- c(created, save_fig(
    forest_plot(ft, paste0(gene, " expression (", group_contrast_label(), ") — ", g$label),
                paste0("Probe ", probe, "; multivariable models adjusted for ", adjust_text())),
    "GEO", "forest", width = 13, height = 1.6 + 0.35 * nrow(ft), subdir = gse, prefix = gse))

  # 전체 공변량 forest (종점별, 가장 많이 보정한 모형)
  for (ep in names(g$endpoints)) {
    fo <- res[[paste(ep, "All")]]$forest
    if (is.null(fo)) next
    created <- c(created, save_fig(
      forest_plot(fo$tab, paste0(fo$name, " Cox: ", ep_lab[[ep]], " — ", g$label),
                  sprintf("n = %d, events = %d, EPV = %.1f", nrow(fo$m$d), fo$m$events, fo$m$epv)),
      "GEO", paste0("forest_multi_", ep), width = 9, height = 1.6 + 0.35 * nrow(fo$tab), subdir = gse, prefix = gse))
  }

  move_stale_outputs("GEO", c("km_", "forest", "cox_", "zph_"), created, subdir = gse, prefix = gse)
  list(files = created, summary = s, probe = pq[pq$chosen, ], tcga = tcga, g = g, d = d)
}

# ---- 실행 -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args)) args else names(cfg$geo_datasets)
unknown <- setdiff(datasets, names(cfg$geo_datasets))
if (length(unknown)) stop("cfg$geo_datasets에 없는 GSE: ", paste(unknown, collapse = ", "))

created  <- character()
failures <- character()
console  <- list()
for (gse in datasets) {
  r <- tryCatch(run_geo_validation(gse, cfg$geo_datasets[[gse]]), error = function(e) {
    failures <<- c(failures, paste0(gse, ": ", conditionMessage(e)))
    cat("  !!", gse, "실패:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) {
    created <- c(created, r$files)
    console[[gse]] <- r
  }
}

cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")

cat("\n===== 요약 =====\n")
for (gse in names(console)) {
  r <- console[[gse]]; s <- r$summary; g <- r$g
  cat(sprintf("%s (검증 대상 TCGA-%s), %s probe %s (%s; mean %.2f, IQR %.2f)\n", gse, g$validates,
              toupper(cfg$primary_gene), r$probe$probe_id, r$probe$rule, r$probe$mean, r$probe$IQR))
  if (!is.null(r$tcga)) {
    cat(sprintf("  TCGA-%s 다변량 HR %s, p %s\n", g$validates,
                fmt_hr(r$tcga$multi_hr, r$tcga$multi_lo, r$tcga$multi_hi), fmt_p(r$tcga$multi_p)))
  }
  for (i in which(s$term == "group")) {
    cat(sprintf("  %-4s %-5s %-38s n = %3d, 사건 %3d: HR %s, p %s%s%s\n",
                s$endpoint[i], s$population[i], s$model[i], s$n[i], s$events[i],
                fmt_hr(s$hr[i], s$lo[i], s$hi[i]), fmt_p(s$p[i]),
                if (is.na(s$direction_match[i])) "" else if (s$direction_match[i]) "  [TCGA와 같은 방향]" else "  [TCGA와 반대 방향]",
                if (s$exploratory[i]) " (탐색적)" else ""))
  }
  lr <- unique(s[, c("endpoint", "population", "logrank_p")])
  cat("  log-rank:", paste(sprintf("%s/%s p %s", lr$endpoint, lr$population, fmt_p(lr$logrank_p)), collapse = ", "), "\n")
}
if (length(failures)) {
  cat("\n실패한 분석 (", length(failures), "개):\n", paste0("  ", failures, collapse = "\n"), "\n", sep = "")
} else {
  cat("\n실패한 분석 없음\n")
}
