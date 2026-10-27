# =============================================================
# 04b_survival_mi.R
# 다변량 Cox (04와 같은 모형) — cfg$mi_vars(기본 CEA)만 다중대체(MI)
#   1. 분석 집단: 생존 정보 + 발현 그룹 + mi_vars 외 공변량 완전 (나머지 결측은 제외)
#   2. mice: mi_vars만 대체 (이진 factor → logreg), 유전자별로 따로 대체
#      대체 모형 예측변수: 다른 공변량 + 발현 그룹 + status + Nelson–Aalen H(t) (White & Royston 2009)
#   3. 대체 데이터셋마다 Cox 적합 → Rubin 규칙으로 풀링 (HR, 95% CI, p, FMI)
#   4. 같은 집단의 완전사례 결과와 비교
# 입력: data/processed/<CANCER>_merged.rds
# 결과 (파일명 앞에 <CANCER>_): output/tables/<CANCER>/MI/
#       cox_multi_MI_<GENE>.{docx,csv}, MI_summary.{docx,csv}, MI_diagnostics.csv
# =============================================================

source("config.R")
source("R/utils.R")
source("R/survival.R")   # Cox 공통 함수 (04와 공유)
suppressPackageStartupMessages(library(mice))

mi_sub <- "MI"

# 생존 정보가 있는 환자 + 모형용 변수 (04의 surv_base와 같음)
# haven 라벨 속성 제거 (mice가 일반 factor/numeric만 다루도록)
surv_base <- function(df) {
  d <- add_scaled_terms(df[!is.na(df$os_months) & !is.na(df$status), ])
  for (v in names(d)) attr(d[[v]], "label") <- NULL
  d
}

# ---- Rubin 풀링 결과 → cox_rows와 같은 형식의 표 ----

pooled_rows <- function(pooled, d, terms, labels) {
  ps  <- summary(pooled, conf.int = TRUE)
  fmi <- setNames(pooled$pooled$fmi, as.character(pooled$pooled$term))
  rownames(ps) <- as.character(ps$term)
  rows <- list()
  for (t in terms) {
    x <- d[[t]]
    if (is.factor(x)) {
      lv <- levels(x)
      for (l in lv) {
        nm  <- paste0(t, l)
        ref <- l == lv[1]
        obs <- !is.na(x)
        rows[[length(rows) + 1]] <- data.frame(
          term = t, variable = labels[[t]], level = l, reference = ref,
          n = sum(x[obs] == l), events = sum(d$status[obs][x[obs] == l]),
          hr = if (ref) 1 else exp(ps[nm, "estimate"]),
          lo = if (ref) NA else exp(ps[nm, "2.5 %"]),
          hi = if (ref) NA else exp(ps[nm, "97.5 %"]),
          p  = if (ref) NA else ps[nm, "p.value"],
          p_overall = NA, fmi = if (ref) NA else fmi[[nm]]
        )
      }
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        term = t, variable = labels[[t]], level = "", reference = FALSE,
        n = nrow(d), events = sum(d$status),
        hr = exp(ps[t, "estimate"]), lo = exp(ps[t, "2.5 %"]), hi = exp(ps[t, "97.5 %"]),
        p = ps[t, "p.value"], p_overall = NA, fmi = fmi[[t]]
      )
    }
  }
  do.call(rbind, rows)
}

# ---- 유전자 1개 ---------------------------------------------------------

mi_cox_gene <- function(g, base, cancer) {
  gene  <- toupper(g)
  res   <- list(summary = data.frame(gene = gene), files = character(), diag = NULL)
  terms <- c("group", intersect(unname(sapply(cfg$cox_covariates, model_term)), names(base)))
  mi_terms <- intersect(cfg$mi_vars, terms)
  if (!length(mi_terms)) stop("cfg$mi_vars가 다변량 공변량에 없음")

  d <- base
  d$group <- d[[paste0(g, "_group")]]
  d <- complete_cases(d[, c("os_months", "status", "group", terms[terms != "group"])],
                      setdiff(terms, mi_terms))
  n_miss <- sum(!stats::complete.cases(d[, mi_terms, drop = FALSE]))
  k   <- n_params(d, terms)
  epv <- sum(d$status) / k
  cat(sprintf("  %-8s MI 집단 n = %d, 사건 %d, %s 결측 %d (%.1f%%), EPV %.1f\n",
              gene, nrow(d), sum(d$status), paste(mi_terms, collapse = "/"),
              n_miss, 100 * n_miss / nrow(d), epv))
  if (sum(d$status) < cfg$min_events) {
    warning(cancer, " ", gene, ": 사건 ", sum(d$status), " < ", cfg$min_events, " → 건너뜀", call. = FALSE)
    res$summary$note <- "skipped: too few events"
    return(res)
  }

  # ---- 완전사례 (같은 집단에서 mi_vars 결측 제외 = 04 다변량과 같은 환자) ----
  dcc <- complete_cases(d, terms)
  cc  <- cox_rows(fit_cox(dcc, terms), dcc, terms, term_labels(terms, gene))
  ccg <- cc[cc$term == "group" & !cc$reference, ]

  # ---- 다중대체 ----
  dm <- d
  dm$na_hazard <- nelsonaalen(dm, os_months, status)
  meth <- make.method(dm)
  meth[] <- ""
  meth[mi_terms] <- ifelse(sapply(dm[mi_terms], function(x) is.factor(x) && nlevels(x) == 2), "logreg",
                    ifelse(sapply(dm[mi_terms], is.factor), "polyreg", "pmm"))
  pred <- make.predictorMatrix(dm)
  pred[, "os_months"] <- 0                     # 시간은 H(t)로만 반영
  imp <- mice(dm, m = cfg$mi_m, maxit = cfg$mi_maxit, method = meth, predictorMatrix = pred,
              seed = cfg$seed, printFlag = FALSE)
  if (nrow(imp$loggedEvents %||% data.frame())) {
    warning(cancer, " ", gene, ": mice loggedEvents ", nrow(imp$loggedEvents), "건", call. = FALSE)
  }
  fits <- lapply(seq_len(cfg$mi_m), function(i) fit_cox(complete(imp, i), terms))
  mi <- pooled_rows(pool(as.mira(fits)), d, terms, term_labels(terms, gene))
  mig <- mi[mi$term == "group" & !mi$reference, ]

  # ---- 결과 표 ----
  tab <- format_cox(mi)
  tab$FMI <- ifelse(is.na(mi$fmi), "", sprintf("%.3f", mi$fmi))
  tab$`Overall p` <- NULL                      # 이진 공변량만 → LRT 전체 p 없음
  res$files <- save_df_table(
    tab, cancer, paste0("cox_multi_MI_", gene), subdir = mi_sub,
    caption = paste0("Multivariable Cox regression for overall survival with multiple imputation: ", gene,
                     " expression, TCGA-", cancer, " (n = ", nrow(d), ", events = ", sum(d$status),
                     ", EPV = ", sprintf("%.1f", epv), "). ", mi_caption(mi_terms, n_miss, nrow(d)),
                     " N/Events for imputed variables are observed values only."))

  # ---- 진단: 관측값 vs 대체값 분포 ----
  res$diag <- do.call(rbind, lapply(mi_terms, function(v) {
    x <- d[[v]]
    lv <- levels(x)[nlevels(x)]
    imputed <- unlist(imp$imp[[v]])
    data.frame(gene = gene, variable = v, level = lv, n_observed = sum(!is.na(x)), n_imputed = sum(is.na(x)),
               prop_observed = mean(x[!is.na(x)] == lv), prop_imputed = mean(imputed == lv),
               m = cfg$mi_m, maxit = cfg$mi_maxit, method = meth[[v]],
               predictors = paste(names(which(pred[v, ] == 1)), collapse = ", "),
               logged_events = nrow(imp$loggedEvents %||% data.frame()))
  }))

  res$summary <- data.frame(
    gene = gene, n_mi = nrow(d), events_mi = sum(d$status), mi_missing = n_miss,
    mi_missing_pct = 100 * n_miss / nrow(d), params = k, epv_mi = epv,
    exploratory = epv < cfg$epv_min,
    cc_n = nrow(dcc), cc_events = sum(dcc$status), cc_hr = ccg$hr, cc_lo = ccg$lo, cc_hi = ccg$hi, cc_p = ccg$p,
    mi_hr = mig$hr, mi_lo = mig$lo, mi_hi = mig$hi, mi_p = mig$p, mi_fmi = mig$fmi, note = NA
  )
  res
}

mi_caption <- function(mi_terms, n_miss, n) {
  paste0(paste(tolower(unname(term_labels(mi_terms))), collapse = ", "), " missing in ", n_miss, "/", n,
         " patients was imputed by chained equations (m = ", cfg$mi_m, ", ", cfg$mi_maxit,
         " iterations) using the other covariates, expression group, event status and the Nelson–Aalen ",
         "cumulative hazard; estimates were pooled with Rubin's rules. FMI = fraction of missing information.")
}

# ---- 실행 -------------------------------------------------------------

created  <- character()
failures <- character()
console  <- list()

for (cancer in cfg$cancers) {
  cat("\n==========", cancer, "==========\n")
  df    <- readRDS(processed_path(cancer, "merged.rds"))
  base  <- surv_base(df)
  genes <- detect_genes(df)
  cat("생존 정보 있는 환자", nrow(base), "명, 사건", sum(base$status), "건 / 유전자", length(genes),
      "개 / 대체 변수:", paste(cfg$mi_vars, collapse = ", "), "\n")

  results <- list()
  for (g in genes) {
    r <- tryCatch(mi_cox_gene(g, base, cancer), error = function(e) {
      failures <<- c(failures, paste0(cancer, " ", toupper(g), ": ", conditionMessage(e)))
      cat("  !!", toupper(g), "실패:", conditionMessage(e), "\n")
      list(summary = data.frame(gene = toupper(g), note = paste("failed:", conditionMessage(e))),
           files = character())
    })
    results[[g]] <- r
    created <- c(created, r$files)
  }

  # ---- 요약: 완전사례 vs MI ----
  s <- bind_rows(lapply(results, `[[`, "summary"))
  # 건너뛴/실패한 유전자만 있어도 표가 만들어지도록 모든 열을 NA로 보장
  summary_cols <- c("n_mi", "events_mi", "mi_missing", "mi_missing_pct", "params", "epv_mi", "exploratory",
                    "cc_n", "cc_events", "cc_hr", "cc_lo", "cc_hi", "cc_p",
                    "mi_hr", "mi_lo", "mi_hi", "mi_p", "mi_fmi", "note")
  for (col in setdiff(summary_cols, names(s))) s[[col]] <- NA
  primary <- toupper(cfg$primary_gene %||% "")
  s$role <- ifelse(s$gene == primary, "primary", "exploratory")
  s <- s[order(s$role != "primary"), ]
  out <- data.frame(
    Gene = s$gene, Role = s$role,
    `Complete-case N` = s$cc_n, `Complete-case events` = s$cc_events,
    `Complete-case HR (95% CI)` = fmt_hr(s$cc_hr, s$cc_lo, s$cc_hi), `Complete-case p` = fmt_p(s$cc_p),
    `MI N` = s$n_mi, `MI events` = s$events_mi,
    `Imputed, n (%)` = ifelse(is.na(s$mi_missing), NA, sprintf("%d (%.1f%%)", s$mi_missing, s$mi_missing_pct)),
    `MI HR (95% CI)` = fmt_hr(s$mi_hr, s$mi_lo, s$mi_hi), `MI p` = fmt_p(s$mi_p),
    FMI = ifelse(is.na(s$mi_fmi), NA, sprintf("%.3f", s$mi_fmi)),
    `MI EPV` = ifelse(is.na(s$epv_mi), NA, sprintf("%.1f", s$epv_mi)),
    check.names = FALSE
  )
  out[[paste0("Exploratory (EPV < ", cfg$epv_min, ")")]] <- s$exploratory
  out$Note <- s$note
  ok <- !is.na(s$mi_missing)
  created <- c(created, save_df_table(
    out, cancer, "MI_summary", landscape = TRUE, font_size = 7, subdir = mi_sub,
    caption = paste0(
      "Overall survival by gene expression (", group_contrast_label(), "), TCGA-", cancer,
      ": multivariable Cox regression adjusted for ", adjust_text(), ", complete-case vs multiple imputation. ",
      if (any(ok)) mi_caption(intersect(cfg$mi_vars, cfg$cox_covariates), s$mi_missing[ok][1], s$n_mi[ok][1]) else "",
      " Patients missing other covariates were excluded from both analyses.")))
  created <- c(created, save_table(s, cancer, "MI_summary_raw.csv", subdir = mi_sub))

  dg <- bind_rows(lapply(results, `[[`, "diag"))
  if (nrow(dg)) created <- c(created, save_table(dg, cancer, "MI_diagnostics.csv", subdir = mi_sub))

  move_stale_outputs(cancer, c("cox_multi_MI_", "MI_"), created, subdir = mi_sub)
  console[[cancer]] <- s
}

# ---- 끝: 생성 파일 + 요약 ----
cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")

cat("\n===== 요약 (발현 그룹 HR: 완전사례 → MI) =====\n")
for (cancer in names(console)) {
  s <- console[[cancer]]
  cat(cancer, "\n")
  for (i in seq_len(nrow(s))) {
    if (is.na(s$mi_hr[i])) { cat(sprintf("  %-8s %s\n", s$gene[i], s$note[i])); next }
    cat(sprintf("  %-8s CC n = %d: HR %s, p %s  →  MI n = %d: HR %s, p %s (FMI %.2f)%s\n",
                s$gene[i], s$cc_n[i], fmt_hr(s$cc_hr[i], s$cc_lo[i], s$cc_hi[i]), fmt_p(s$cc_p[i]),
                s$n_mi[i], fmt_hr(s$mi_hr[i], s$mi_lo[i], s$mi_hi[i]), fmt_p(s$mi_p[i]), s$mi_fmi[i],
                if (s$exploratory[i]) " (탐색적)" else ""))
  }
}
if (length(failures)) {
  cat("\n실패한 분석 (", length(failures), "개):\n", paste0("  ", failures, collapse = "\n"), "\n", sep = "")
} else {
  cat("\n실패한 분석 없음\n")
}
