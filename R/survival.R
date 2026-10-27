# =============================================================
# R/survival.R
# 생존분석 공통 함수 (04_survival.R, 06_external_geo.R)
# 사용: source("config.R"); source("R/utils.R"); source("R/survival.R")
# 시간/사건 열 기본값은 04의 os_months / status (06은 종점별로 지정)
# =============================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(patchwork)
})

# ---- 모형 변수 ---------------------------------------------------------

# 설정 변수명 → 모형 항 (cfg$covariate_scale에 있는 연속형은 단위 변환, 예: age → age_per10)
model_term <- function(v) scaled_term(v)

term_labels <- function(terms, gene = NULL) {
  labs <- sapply(terms, function(t) {
    sl <- scaled_label(t)
    if (!is.null(sl)) return(sl)
    if (t == "group") return(paste(gene, "expression"))
    if (t == "expr_log2") return(paste(gene, "expression, log2 (per 1 unit)"))
    if (t %in% names(var_labels)) var_labels[[t]] else t
  })
  setNames(labs, terms)
}

# 캡션용 보정 변수 목록 (cfg$cox_covariates → "age (per 10 years), sex and pathologic stage")
adjust_text <- function(covariates = cfg$cox_covariates) {
  labs <- tolower(unname(term_labels(unname(sapply(covariates, model_term)))))
  if (length(labs) < 2) return(paste(labs, collapse = ""))
  paste0(paste(head(labs, -1), collapse = ", "), " and ", tail(labs, 1))
}

# 모형 항들이 모두 있는 행만 (완전 사례), 빈 factor 수준 제거
complete_cases <- function(d, terms) {
  d <- d[stats::complete.cases(d[, terms, drop = FALSE]), , drop = FALSE]
  for (t in terms) if (is.factor(d[[t]])) d[[t]] <- droplevels(d[[t]])
  d
}

n_params <- function(d, terms) {
  sum(sapply(terms, function(t) if (is.factor(d[[t]])) nlevels(d[[t]]) - 1 else 1))
}

fit_cox <- function(d, terms, strata = NULL, time = "os_months", event = "status") {
  rhs <- c(terms, if (!is.null(strata)) paste0("strata(", strata, ")"))
  coxph(as.formula(paste0("Surv(", time, ", ", event, ") ~ ", paste(rhs, collapse = " + "))), data = d)
}

# 다변량 모형 구성: EPV < cfg$epv_min 이면 cfg$stage_full → cfg$stage_collapsed,
# 그래도 부족하면 exploratory = TRUE
build_multi <- function(d, terms, event = "status") {
  build <- function(tm) {
    dm <- complete_cases(d, tm)
    list(terms = tm, d = dm, events = sum(dm[[event]]), k = n_params(dm, tm))
  }
  m <- build(terms)
  m$stage_used <- if (cfg$stage_full %in% terms) cfg$stage_full else NA
  if (m$events / m$k < cfg$epv_min && cfg$stage_full %in% terms &&
      cfg$stage_collapsed %in% names(d)) {
    m <- build(replace(terms, terms == cfg$stage_full, cfg$stage_collapsed))
    m$stage_used <- cfg$stage_collapsed
  }
  # 표시용 stage 수준 (예: "I–IV", "I–II vs III–IV")
  m$stage_levels <- if (is.na(m$stage_used)) NA else {
    lv <- levels(m$d[[m$stage_used]])
    if (length(lv) == 2) paste(lv, collapse = " vs ") else paste0(lv[1], "–", lv[length(lv)])
  }
  m$epv <- m$events / m$k
  m$exploratory <- m$epv < cfg$epv_min
  m$dropped <- nrow(d) - nrow(m$d)
  m
}

# coxph → 항목별 표 (factor는 기준 수준 행 포함, 다수준 factor는 LRT 전체 p)
cox_rows <- function(fit, d, terms, labels, event = "status") {
  s  <- summary(fit)
  ci <- s$conf.int
  co <- s$coefficients
  overall <- tryCatch(drop1(fit, test = "Chisq"), error = function(e) NULL)
  rows <- list()
  for (t in terms) {
    x <- d[[t]]
    p_all <- if (!is.null(overall) && t %in% rownames(overall)) overall[t, "Pr(>Chi)"] else NA
    if (is.factor(x)) {
      lv <- levels(x)
      for (l in lv) {
        nm  <- paste0(t, l)
        ref <- l == lv[1]
        rows[[length(rows) + 1]] <- data.frame(
          term = t, variable = labels[[t]], level = l, reference = ref,
          n = sum(x == l), events = sum(d[[event]][x == l]),
          hr = if (ref) 1 else ci[nm, "exp(coef)"],
          lo = if (ref) NA else ci[nm, "lower .95"],
          hi = if (ref) NA else ci[nm, "upper .95"],
          p  = if (ref) NA else co[nm, "Pr(>|z|)"],
          p_overall = if (ref && length(lv) > 2) p_all else NA
        )
      }
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        term = t, variable = labels[[t]], level = "", reference = FALSE,
        n = nrow(d), events = sum(d[[event]]),
        hr = ci[t, "exp(coef)"], lo = ci[t, "lower .95"], hi = ci[t, "upper .95"],
        p = co[t, "Pr(>|z|)"], p_overall = NA
      )
    }
  }
  do.call(rbind, rows)
}

# 표 출력용 형식
format_cox <- function(tab) {
  data.frame(
    Variable        = ifelse(duplicated(tab$variable), "", tab$variable),
    Level           = ifelse(tab$reference, paste0(tab$level, " (ref)"), tab$level),
    N               = tab$n,
    Events          = tab$events,
    `HR (95% CI)`   = ifelse(tab$reference, "Reference", fmt_hr(tab$hr, tab$lo, tab$hi)),
    p               = fmt_p(tab$p),
    `Overall p`     = fmt_p(tab$p_overall),
    check.names = FALSE
  )
}

# cox.zph → 항목별 + GLOBAL
zph_rows <- function(fit, gene, model) {
  z <- cox.zph(fit)
  data.frame(gene = gene, model = model, term = rownames(z$table),
             chisq = z$table[, "chisq"], df = z$table[, "df"], p = z$table[, "p"],
             flag = z$table[, "p"] < 0.05, row.names = NULL)
}

# 시간 분할 Cox: 다변량 모형과 같은 공변량/환자, 비교 그룹 효과를 0–cut / >cut 구간으로 분리
# (유전자 그룹 항이 PH 위반일 때만 사용)
time_split_cox <- function(d, terms, cut, time = "os_months", event = "status") {
  # 시간 0인 환자가 있으면 시작점을 그보다 앞으로 (위험집합은 일반 Cox와 동일)
  zero <- if (any(d[[time]] <= 0)) min(d[[time]]) - 1 else 0
  sp <- survSplit(as.formula(paste0("Surv(", time, ", ", event, ") ~ .")), data = d,
                  cut = cut, episode = "period", zero = zero)
  sp$alt_early <- as.numeric(sp$group == group_alt() & sp$period == 1)
  sp$alt_late  <- as.numeric(sp$group == group_alt() & sp$period == 2)
  rhs <- c("alt_early", "alt_late", setdiff(terms, "group"))
  fit <- coxph(as.formula(paste0("Surv(tstart, ", time, ", ", event, ") ~ ", paste(rhs, collapse = " + "))),
               data = sp)
  ci <- summary(fit)$conf.int
  co <- summary(fit)$coefficients
  est <- data.frame(
    hr_early = ci["alt_early", 1], lo_early = ci["alt_early", 3], hi_early = ci["alt_early", 4],
    p_early = co["alt_early", "Pr(>|z|)"],
    hr_late = ci["alt_late", 1], lo_late = ci["alt_late", 3], hi_late = ci["alt_late", 4],
    p_late = co["alt_late", "Pr(>|z|)"])
  ev <- tapply(sp[[event]], list(sp$period, sp$group), sum)
  msg <- sprintf("PH 위반 → 시간 분할 Cox (%d개월): 0–%d HR %s (사건 기준/비교 %d/%d), >%d HR %s (사건 기준/비교 %d/%d)",
                 cut, cut, fmt_hr(est$hr_early, est$lo_early, est$hi_early), ev[1, group_ref()], ev[1, group_alt()],
                 cut, fmt_hr(est$hr_late, est$lo_late, est$hi_late), ev[2, group_ref()], ev[2, group_alt()])
  list(est = est, events = ev, msg = msg)
}

# RMST 차이 (비교 − 기준 그룹, tau = cfg$rmst_tau; 추적이 짧으면 tau를 줄이고 note에 기록)
rmst_diff <- function(d, tau = cfg$rmst_tau, time = "os_months", event = "status") {
  tau_max <- min(tapply(d[[time]], d$group, max))
  note <- ""
  if (tau > tau_max) {
    tau  <- floor(tau_max)
    note <- paste0("RMST tau reduced to ", tau, " months (follow-up)")
  }
  rm <- survRM2::rmst2(d[[time]], d[[event]], as.numeric(d$group == group_alt()), tau = tau)
  ur <- rm$unadjusted.result
  rd <- ur[grep("^RMST \\(arm=1\\)-\\(arm=0\\)", rownames(ur)), ]
  est <- data.frame(rmst_tau = tau,
                    rmst_low = rm$RMST.arm0$rmst[["Est."]], rmst_high = rm$RMST.arm1$rmst[["Est."]],
                    rmst_diff = rd[["Est."]], rmst_lo = rd[["lower .95"]], rmst_hi = rd[["upper .95"]],
                    rmst_p = rd[["p"]])
  list(est = est, note = note)
}

# ---- Kaplan–Meier -----------------------------------------------------

km_stats <- function(fit, times) {
  tb  <- summary(fit)$table
  grp <- sub("^group=", "", rownames(tb))
  out <- data.frame(group = grp, n = tb[, "n.max"], events = tb[, "events"],
                    median = tb[, "median"], median_lo = tb[, "0.95LCL"],
                    median_hi = tb[, "0.95UCL"], row.names = NULL)
  s <- summary(fit, times = times, extend = TRUE)
  for (tm in times) {
    k   <- s$time == tm
    ok  <- s$n.risk[k] > 0            # 추적 범위를 넘으면 추정 불가 → NA
    idx <- match(out$group, sub("^group=", "", as.character(s$strata[k])))
    col <- paste0("surv_", tm / 12, "y")
    out[[col]]          <- ifelse(ok, s$surv[k], NA)[idx]
    out[[paste0(col, "_lo")]] <- ifelse(ok, s$lower[k], NA)[idx]
    out[[paste0(col, "_hi")]] <- ifelse(ok, s$upper[k], NA)[idx]
  }
  out
}

fmt_median <- function(m, lo, hi) {
  f <- function(x) ifelse(is.na(x), "NR", sprintf("%.1f", x))
  paste0(f(m), " (", f(lo), "–", f(hi), ")")
}

fmt_surv <- function(s, lo, hi) {
  ifelse(is.na(s), NA_character_,
         sprintf("%.1f%% (%.1f–%.1f)", 100 * s, 100 * lo, 100 * hi))
}

km_plot <- function(fit, d, gene, cancer, p_lr,
                    title = paste0(gene, " expression — TCGA-", cancer),
                    ylab = "Overall survival probability") {
  # survminer + ggplot2 4.x의 "Ignoring unknown labels" 메시지 억제 (결과에는 영향 없음)
  p <- suppressMessages(ggsurvplot(
    fit, data = d,
    palette = unname(cfg$group_colors[cfg$group_levels]),
    legend.labs = cfg$group_levels, legend.title = paste(gene, "expression"),
    pval = paste("Log-rank p", ifelse(p_lr < 0.001, "< 0.001", paste("=", sprintf("%.3f", p_lr)))),
    pval.size = 4, censor = TRUE, censor.shape = "|", censor.size = 3,
    risk.table = TRUE, risk.table.title = "Number at risk",
    risk.table.y.text = FALSE, fontsize = 3.5,
    break.time.by = 12, xlab = "Time (months)", ylab = ylab,
    title = title,
    ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable()
  ))
  suppressMessages(p$plot / p$table + plot_layout(heights = c(3, 1)))
}

# ---- Forest plot -------------------------------------------------------

forest_plot <- function(tab, title, subtitle = NULL) {
  if (!"label" %in% names(tab)) {
    tab$label <- ifelse(tab$level == "", tab$variable,
                        paste0(tab$variable, ": ", tab$level, ifelse(tab$reference, " (ref)", "")))
  }
  tab$row   <- rev(seq_len(nrow(tab)))
  tab$hr_ci <- ifelse(tab$reference, "Reference", fmt_hr(tab$hr, tab$lo, tab$hi))
  tab$p_txt <- ifelse(tab$reference, "", fmt_p(tab$p))
  if (!"shape" %in% names(tab)) tab$shape <- ifelse(tab$reference, "ref", "est")
  ylim <- c(0.4, nrow(tab) + 0.9)
  head_y <- nrow(tab) + 0.8

  txt <- function(d, x, label, header, hjust = 0) {
    ggplot(d, aes(x = x, y = row)) +
      geom_text(aes(label = .data[[label]]), hjust = hjust, size = 3.3) +
      annotate("text", x = x, y = head_y, label = header, hjust = hjust,
               fontface = "bold", size = 3.3) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = ylim) +
      theme_void()
  }
  left  <- txt(tab, 0, "label", "")
  right <- txt(tab, 0, "hr_ci", "HR (95% CI)") + txt(tab, 0, "p_txt", "p") +
    plot_layout(widths = c(2, 1))

  rng <- range(c(tab$lo, tab$hi, tab$hr), na.rm = TRUE)
  mid <- ggplot(tab, aes(x = hr, y = row)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y", na.rm = TRUE) +
    geom_point(aes(shape = shape), size = 2.6, fill = "white", na.rm = TRUE) +
    scale_shape_manual(values = c(est = 15, ref = 23, exploratory = 22), guide = "none") +
    scale_x_log10(limits = c(min(rng[1], 0.8), max(rng[2], 1.25))) +
    scale_y_continuous(limits = ylim, breaks = NULL) +
    labs(x = "Hazard ratio (log scale)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

  (left | mid | right) + plot_layout(widths = c(2.3, 2.2, 2)) +
    plot_annotation(title = title, subtitle = subtitle)
}
