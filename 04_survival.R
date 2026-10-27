# =============================================================
# 04_survival.R
# 전체 생존(OS) 분석: Surv(os_months, status), status 1 = 사망
#   1. Kaplan–Meier (Low vs High) + log-rank, 중앙 생존기간, 3/5년 생존율
#   2. 단변량 Cox: 유전자 그룹, 유전자 연속형(log2), 임상 공변량
#   3. 다변량 Cox: 유전자 그룹 + age(10세 단위) + gender + stage (EPV 부족 시 stage I–II vs III–IV)
#   4. 민감도 분석: 다변량 + strata(cfg$strata_var)
#   5. 비례위험 가정 (cox.zph)
#   6. 유전자 간 BH 보정 q-value
# 입력: data/processed/<CANCER>_merged.rds (신보조요법 제외 적용됨)
# 결과 (파일명 앞에 <CANCER>_): output/tables/<CANCER>/survival_summary, km_summary, cox_uni, cox_multi_<GENE>, ph_tests
#       output/figures/<CANCER>/km_<GENE>, forest_multi_<GENE>, forest_genes, zph_<GENE>
# =============================================================

source("config.R")
source("R/utils.R")
source("R/survival.R")   # Cox/KM/forest 공통 함수 (06과 공유)

# 생존 정보가 있는 환자 + 모형용 변수
surv_base <- function(df) {
  add_scaled_terms(df[!is.na(df$os_months) & !is.na(df$status), ])
}

# ---- 유전자 1개 분석 ---------------------------------------------------

analyze_gene <- function(g, base, cancer, uni_cov_terms) {
  gene <- toupper(g)
  res  <- list(summary = data.frame(gene = gene), files = character(),
               uni = NULL, multi = NULL, zph = NULL, km = NULL)

  d <- base
  d$group     <- d[[paste0(g, "_group")]]
  d$expr_log2 <- d[[paste0(g, "_expression_log2")]]
  d <- d[!is.na(d$group), ]
  n_ev <- sum(d$status)
  cnt  <- table(factor(d$group, levels = cfg$group_levels))
  res$summary <- data.frame(gene = gene, n = nrow(d), events = n_ev,
                            n_low = cnt[[group_ref()]], n_high = cnt[[group_alt()]])
  cat(sprintf("  %-8s n = %d (%s %d / %s %d), 사건 %d\n",
              gene, nrow(d), group_ref(), cnt[[group_ref()]], group_alt(), cnt[[group_alt()]], n_ev))

  if (n_ev < cfg$min_events || any(cnt == 0)) {
    msg <- if (any(cnt == 0)) "빈 발현 그룹" else paste0("사건 ", n_ev, " < ", cfg$min_events)
    warning(cancer, " ", gene, ": ", msg, " → 생존분석 건너뜀", call. = FALSE)
    res$summary$note <- paste("skipped:", msg)
    return(res)
  }

  # ---- 1. KM ----
  km <- surv_fit(Surv(os_months, status) ~ group, data = d)
  p_lr <- survdiff(Surv(os_months, status) ~ group, data = d)$pvalue
  ks <- km_stats(km, cfg$km_times)
  res$km <- data.frame(gene = gene, ks)
  res$files <- c(res$files, suppressMessages(save_fig(km_plot(km, d, gene, cancer, p_lr), cancer,
                                                      paste0("km_", gene), width = 7, height = 6)))
  low  <- ks[ks$group == group_ref(), ]
  high <- ks[ks$group == group_alt(), ]

  # ---- 2. 단변량 Cox (그룹, 연속형) ----
  labs_g <- term_labels(c("group", "expr_log2"), gene)
  fit_u  <- fit_cox(d, "group")
  uni_g  <- cox_rows(fit_u, d, "group", labs_g)
  dc     <- d[!is.na(d$expr_log2), ]
  uni_c  <- cox_rows(fit_cox(dc, "expr_log2"), dc, "expr_log2", labs_g)
  res$uni <- rbind(data.frame(gene = gene, uni_g), data.frame(gene = gene, uni_c))
  zph <- zph_rows(fit_u, gene, "univariable")

  # ---- 3. 다변량 Cox (EPV 확인 후 stage 변수 결정) ----
  cov_terms <- intersect(unname(sapply(cfg$cox_covariates, model_term)), names(d))
  m <- build_multi(d, c("group", cov_terms))
  stage_used   <- m$stage_used
  stage_levels <- m$stage_levels
  epv          <- m$epv
  exploratory  <- m$exploratory
  dropped      <- m$dropped
  cat(sprintf("           다변량: n = %d (결측 제외 %d), 사건 %d, 모수 %d, EPV %.1f, stage = %s%s\n",
              nrow(m$d), dropped, m$events, m$k, epv, stage_used,
              if (exploratory) " → 탐색적(EPV 부족)" else ""))
  if (exploratory) {
    warning(cancer, " ", gene, ": EPV ", sprintf("%.1f", epv), " < ", cfg$epv_min,
            " → 다변량 결과는 탐색적 (exploratory = TRUE)", call. = FALSE)
  }
  fit_m <- fit_cox(m$d, m$terms)
  multi <- cox_rows(fit_m, m$d, m$terms, term_labels(m$terms, gene))
  res$multi <- multi
  zph <- rbind(zph, zph_rows(fit_m, gene, "multivariable"))
  mg <- multi[multi$term == "group" & !multi$reference, ]

  tab_m <- format_cox(multi)
  res$files <- c(res$files, save_df_table(
    tab_m, cancer, paste0("cox_multi_", gene),
    caption = paste0("Multivariable Cox regression for overall survival: ", gene,
                     " expression, TCGA-", cancer,
                     " (n = ", nrow(m$d), ", events = ", m$events, ", EPV = ",
                     sprintf("%.1f", epv), if (exploratory) "; exploratory" else "", ")")))
  res$files <- c(res$files, save_fig(
    forest_plot(multi, paste0("Multivariable Cox: ", gene, " expression — TCGA-", cancer),
                paste0("n = ", nrow(m$d), ", events = ", m$events, ", EPV = ", sprintf("%.1f", epv),
                       if (exploratory) " (exploratory)" else "")),
    cancer, paste0("forest_multi_", gene), width = 9, height = 1.6 + 0.35 * nrow(multi)))

  # ---- 4. 민감도 분석: strata(cfg$strata_var), 소수 수준은 cfg$collapse_small_levels로 병합 ----
  ds <- m$d
  ds$strata_grp <- collapse_small(ds[[cfg$strata_var]], cfg$strata_var)
  fit_s <- fit_cox(ds, m$terms, strata = "strata_grp")
  sg <- cox_rows(fit_s, ds, m$terms, term_labels(m$terms, gene))
  sg <- sg[sg$term == "group" & !sg$reference, ]

  # ---- 5. 비례위험 그림 (위반 모형만) ----
  res$zph <- zph
  if (any(zph$flag)) {
    f <- file.path(out_dir(cancer, "figures"), out_file(cancer, paste0("zph_", gene, ".pdf")))
    grDevices::cairo_pdf(f, width = 7, height = 5, onefile = TRUE)
    for (mod in unique(zph$model[zph$flag])) {
      z <- cox.zph(if (mod == "univariable") fit_u else fit_m)
      for (i in seq_len(nrow(z$table) - 1)) {
        plot(z[i], resid = TRUE, se = TRUE, col = "#C8442F",
             main = paste0(gene, " — ", mod, ": ", rownames(z$table)[i],
                           " (p = ", fmt_p(z$table[i, "p"]), ")"))
        abline(h = coef(if (mod == "univariable") fit_u else fit_m)[i],
               lty = 2, col = "grey40")
      }
    }
    grDevices::dev.off()
    res$files <- c(res$files, f)
  }

  # ---- 6. 시간 분할 Cox: 유전자 그룹 항이 PH 위반(p < 0.05)일 때만 ----
  # 다변량 모형과 같은 공변량/환자, High 효과를 0–cut / >cut 구간으로 분리
  ts <- data.frame(hr_early = NA, lo_early = NA, hi_early = NA, p_early = NA,
                   hr_late = NA, lo_late = NA, hi_late = NA, p_late = NA)
  if (any(zph$flag & zph$term == "group")) {
    tsc <- time_split_cox(m$d, m$terms, cfg$ph_split_months)
    ts  <- tsc$est
    cat("           ", tsc$msg, "\n", sep = "")
  }

  # ---- 7. RMST 차이 (High − Low, 비례위험 가정 불필요) ----
  rr <- rmst_diff(d)
  rmst <- rr$est
  tau_note <- rr$note
  if (nzchar(tau_note)) {
    warning(cancer, " ", gene, ": 추적기간이 짧아 RMST tau를 ", rmst$rmst_tau, "개월로 줄임", call. = FALSE)
  }

  res$summary <- cbind(res$summary, ts, rmst, data.frame(
    median_os_low  = fmt_median(low$median, low$median_lo, low$median_hi),
    median_os_high = fmt_median(high$median, high$median_lo, high$median_hi),
    logrank_p = p_lr,
    uni_hr = uni_g$hr[2], uni_lo = uni_g$lo[2], uni_hi = uni_g$hi[2], uni_p = uni_g$p[2],
    multi_n = nrow(m$d), multi_dropped = dropped, multi_events = m$events,
    multi_hr = mg$hr, multi_lo = mg$lo, multi_hi = mg$hi, multi_p = mg$p,
    params = m$k, epv = epv, stage_var = stage_used, exploratory = exploratory,
    strata_hr = sg$hr, strata_lo = sg$lo, strata_hi = sg$hi, strata_p = sg$p,
    n_strata = length(unique(ds$strata_grp)), stage_levels = stage_levels,
    ph_global_p_uni = zph$p[zph$model == "univariable" & zph$term == "GLOBAL"],
    ph_global_p_multi = zph$p[zph$model == "multivariable" & zph$term == "GLOBAL"],
    ph_flag = any(zph$flag), note = tau_note
  ))
  res
}

# 주 가설 유전자는 제외하고 탐색적 유전자끼리만 BH 보정
bh_exploratory <- function(p, role) {
  q <- rep(NA_real_, length(p))
  i <- role == "exploratory"
  q[i] <- p.adjust(p[i], "BH")
  q
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
  cat("생존 정보 있는 환자", nrow(base), "명, 사건", sum(base$status), "건 / 유전자", length(genes), "개\n")

  # ---- 단변량 Cox: 임상 공변량 ----
  uni_cov_terms <- unique(unname(sapply(cfg$cox_uni_covariates, model_term)))
  missing_cov <- setdiff(uni_cov_terms, names(base))
  if (length(missing_cov)) cat("없는 공변량 (건너뜀):", paste(missing_cov, collapse = ", "), "\n")
  uni_cov_terms <- intersect(uni_cov_terms, names(base))
  uni_cov <- do.call(rbind, lapply(uni_cov_terms, function(t) {
    dt <- complete_cases(base, t)
    cat(sprintf("  단변량 %-24s n = %d, 사건 %d\n", t, nrow(dt), sum(dt$status)))
    cox_rows(fit_cox(dt, t), dt, t, term_labels(t))
  }))

  # ---- 유전자별 ----
  results <- list()
  for (g in genes) {
    r <- tryCatch(analyze_gene(g, base, cancer, uni_cov_terms), error = function(e) {
      failures <<- c(failures, paste0(cancer, " ", toupper(g), ": ", conditionMessage(e)))
      cat("  !!", toupper(g), "실패:", conditionMessage(e), "\n")
      list(summary = data.frame(gene = toupper(g), note = paste("failed:", conditionMessage(e))),
           files = character())
    })
    results[[g]] <- r
    created <- c(created, r$files)
  }

  # ---- 요약 + BH q ----
  s <- bind_rows(lapply(results, `[[`, "summary"))
  # 건너뛴/실패한 유전자만 있어도 표가 만들어지도록 모든 열을 NA로 보장
  summary_cols <- c("n", "events", "n_low", "n_high", "median_os_low", "median_os_high",
                    "logrank_p", "uni_hr", "uni_lo", "uni_hi", "uni_p",
                    "multi_n", "multi_dropped", "multi_events",
                    "multi_hr", "multi_lo", "multi_hi", "multi_p", "params", "epv",
                    "stage_var", "stage_levels", "exploratory", "strata_hr", "strata_lo", "strata_hi", "strata_p",
                    "n_strata", "ph_global_p_uni", "ph_global_p_multi", "ph_flag",
                    "hr_early", "lo_early", "hi_early", "p_early",
                    "hr_late", "lo_late", "hi_late", "p_late",
                    "rmst_tau", "rmst_low", "rmst_high", "rmst_diff", "rmst_lo", "rmst_hi", "rmst_p",
                    "note")
  for (col in setdiff(summary_cols, names(s))) s[[col]] <- NA

  # 주 가설 유전자: raw p 보고, BH는 나머지(탐색적) 유전자끼리
  primary <- toupper(cfg$primary_gene %||% "")
  if (nzchar(primary) && !primary %in% s$gene) {
    warning(cancer, ": primary_gene ", primary, "이(가) 이 암종 유전자 목록에 없음", call. = FALSE)
  }
  s$role <- ifelse(s$gene == primary, "primary", "exploratory")
  s <- s[order(s$role != "primary"), ]
  s$logrank_q <- bh_exploratory(s$logrank_p, s$role)
  s$uni_q     <- bh_exploratory(s$uni_p, s$role)
  s$multi_q   <- bh_exploratory(s$multi_p, s$role)
  created <- c(created, save_table(s, cancer, "survival_summary_raw.csv"))

  q_txt <- function(q) ifelse(s$role == "primary", "–", fmt_p(q))
  sc <- cfg$ph_split_months
  out <- data.frame(
    Gene = s$gene, Role = s$role, N = s$n, Events = s$events,
    `Median OS Low, months (95% CI)`  = s$median_os_low,
    `Median OS High, months (95% CI)` = s$median_os_high,
    `Log-rank p` = fmt_p(s$logrank_p), `Log-rank q` = q_txt(s$logrank_q),
    `Univariable HR (95% CI)` = fmt_hr(s$uni_hr, s$uni_lo, s$uni_hi),
    `Univariable p` = fmt_p(s$uni_p), `Univariable q` = q_txt(s$uni_q),
    `Multivariable HR (95% CI)` = fmt_hr(s$multi_hr, s$multi_lo, s$multi_hi),
    `Multivariable p` = fmt_p(s$multi_p), `Multivariable q` = q_txt(s$multi_q),
    EPV = ifelse(is.na(s$epv), NA, sprintf("%.1f", s$epv)),
    `Stage variable` = s$stage_levels,
    `Strata HR (95% CI)` = fmt_hr(s$strata_hr, s$strata_lo, s$strata_hi),
    `Strata p` = fmt_p(s$strata_p),
    `PH global p` = fmt_p(s$ph_global_p_multi),
    check.names = FALSE
  )
  out[[paste0("HR 0–", sc, " mo (95% CI)")]] <- fmt_hr(s$hr_early, s$lo_early, s$hi_early)
  out[[paste0("p 0–", sc, " mo")]]           <- fmt_p(s$p_early)
  out[[paste0("HR >", sc, " mo (95% CI)")]]       <- fmt_hr(s$hr_late, s$lo_late, s$hi_late)
  out[[paste0("p >", sc, " mo")]]                 <- fmt_p(s$p_late)
  out[[paste0("RMST difference at ", cfg$rmst_tau, " mo, months (95% CI)")]] <-
    ifelse(is.na(s$rmst_diff), NA, sprintf("%.2f (%.2f, %.2f)", s$rmst_diff, s$rmst_lo, s$rmst_hi))
  out[["RMST p"]] <- fmt_p(s$rmst_p)
  names(out) <- sub("^Strata ", paste0("Strata(", toupper(cfg$strata_var), ") "), names(out))
  out[[paste0("Exploratory (EPV < ", cfg$epv_min, ")")]] <- s$exploratory
  out$Note <- s$note
  created <- c(created, save_df_table(
    out, cancer, "survival_summary", landscape = TRUE, font_size = 6,
    caption = paste0(
      "Overall survival by gene expression (", group_contrast_label(), "), TCGA-", cancer,
      ". HR from Cox regression; multivariable model adjusted for ", adjust_text(), ". Sensitivity: stratified by ", cfg$strata_var, ". ",
      "q = Benjamini–Hochberg across exploratory genes only; the primary gene (", primary,
      ") is reported with its unadjusted p. Time-split HRs (0–", sc, " / >", sc,
      " months, adjusted) are shown only for genes whose expression term violated proportional hazards. ",
      "RMST difference = High − Low restricted mean survival time up to ", cfg$rmst_tau, " months.")))

  # ---- KM 요약 ----
  km <- bind_rows(lapply(results, `[[`, "km"))
  if (nrow(km)) {
    km_out <- data.frame(Gene = km$gene, Group = km$group, N = km$n, Events = km$events,
                         `Median OS, months (95% CI)` = fmt_median(km$median, km$median_lo, km$median_hi),
                         check.names = FALSE)
    for (tm in cfg$km_times) {
      col <- paste0("surv_", tm / 12, "y")
      km_out[[paste0(tm / 12, "-year OS (95% CI)")]] <-
        fmt_surv(km[[col]], km[[paste0(col, "_lo")]], km[[paste0(col, "_hi")]])
    }
    created <- c(created, save_df_table(km_out, cancer, "km_summary",
                                        caption = paste0("Kaplan–Meier estimates by gene expression, TCGA-", cancer,
                                                         ". NR = not reached.")))
  }

  # ---- 단변량 Cox 표 ----
  uni_all <- rbind(bind_rows(lapply(results, `[[`, "uni")) |> select(-gene), uni_cov)
  created <- c(created, save_df_table(
    format_cox(uni_all), cancer, "cox_uni",
    caption = paste0("Univariable Cox regression for overall survival, TCGA-", cancer)))

  # ---- 비례위험 검정 ----
  zph <- bind_rows(lapply(results, `[[`, "zph"))
  if (nrow(zph)) {
    created <- c(created, save_table(transform(zph, p = signif(p, 4), chisq = round(chisq, 3)),
                                     cancer, "ph_tests.csv"))
  }

  # ---- 유전자 forest plot (다변량 유전자 HR) ----
  fg <- s[!is.na(s$multi_hr), ]
  if (nrow(fg)) {
    fg_tab <- data.frame(
      label = sprintf("%s%s (n = %d, events = %d)", fg$gene,
                      ifelse(fg$role == "primary", " [primary]", ""), fg$multi_n, fg$multi_events),
      hr = fg$multi_hr, lo = fg$multi_lo, hi = fg$multi_hi, p = fg$multi_p, reference = FALSE,
      shape = ifelse(fg$exploratory, "exploratory", "est"))
    sub <- paste0(group_contrast_label(), " expression; adjusted for ", adjust_text(),
                  if (any(fg$exploratory)) paste0(". Open squares: exploratory (EPV < ", cfg$epv_min, ")") else "")
    created <- c(created, save_fig(
      forest_plot(fg_tab, paste0("Multivariable Cox: gene expression — TCGA-", cancer), sub),
      cancer, "forest_genes", width = 9, height = 1.6 + 0.35 * nrow(fg_tab)))
  }

  # ---- 더 이상 없는 유전자의 이전 결과 → _stale/ ----
  move_stale_outputs(cancer, c("km_", "forest_multi_", "zph_", "cox_multi_"), created)

  console[[cancer]] <- s
}

# ---- 끝: 생성 파일 + 요약 ----
cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")

cat("\n===== 요약 =====\n")
for (cancer in names(console)) {
  s <- console[[cancer]]
  lr <- s$gene[!is.na(s$logrank_p) & s$logrank_p < 0.05]
  mv <- s$gene[!is.na(s$multi_p) & s$multi_p < 0.05]
  pr <- s[s$role == "primary", ]
  cat(cancer, "\n")
  if (nrow(pr)) {
    cat(sprintf("  주 가설 %s: log-rank p %s, 다변량 HR %s p %s, RMST 차이 %s개월 (p %s)\n",
                pr$gene, fmt_p(pr$logrank_p), fmt_hr(pr$multi_hr, pr$multi_lo, pr$multi_hi),
                fmt_p(pr$multi_p), sprintf("%.2f", pr$rmst_diff), fmt_p(pr$rmst_p)))
  }
  cat("  log-rank p < 0.05     :", if (length(lr)) paste(lr, collapse = ", ") else "없음", "\n")
  cat("  다변량 유전자 p < 0.05:", if (length(mv)) paste(mv, collapse = ", ") else "없음",
      if (any(s$exploratory %in% TRUE)) paste0("(일부/전체 탐색적, EPV < ", cfg$epv_min, ")") else "", "\n")
}
if (length(failures)) {
  cat("\n실패한 분석 (", length(failures), "개):\n", paste0("  ", failures, collapse = "\n"), "\n", sep = "")
} else {
  cat("\n실패한 분석 없음\n")
}
