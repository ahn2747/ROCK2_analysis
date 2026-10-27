# =============================================================
# 05_gsea.R
# 발현 그룹(High vs Low)별 차등발현(DESeq2) + GSEA(fgsea, MSigDB)
#   - 그룹: merged 데이터의 기존 <gene>_group 그대로 사용 (생존분석과 동일, 재분할 없음)
#   - 주 design: ~ <cfg$gsea_design_covariates> + group  (소수 수준은 cfg$collapse_small_levels로 병합)
#   - 민감도: ~ <cfg$gsea_sensitivity_covariates> + group (기본: 공변량 없음) → robust_no_covariate
#   - 순위: DESeq2 Wald stat → MSigDB .chip으로 옛 기호 재매핑 → 중복 평균
# 입력: data/processed/<C>_counts.rds (02a), <C>_merged.rds, database/TCGA_<C>_RNAseq_Expression.csv
# 결과 (파일명 앞에 <C>_): output/tables/<C>/gsea/<GENE>/  de_results.csv, gsea_<collection>.csv, gsea_QC.csv, gsea_summary
#       output/figures/<C>/gsea/<GENE>/ volcano, pca, nes_<collection>, enrichment_<pathway>
#       data/processed/<C>_dds_<GENE>.rds
# =============================================================

source("config.R")
source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
})

# 암종 하나만 실행: Rscript 05_gsea.R COAD  (명령줄 인수가 cfg$cancers보다 우선)
# 메모리 절약을 위해 암종마다 별도 프로세스로 실행 권장
args <- toupper(commandArgs(trailingOnly = TRUE))
if (length(args)) {
  bad <- setdiff(args, cfg$cancers)
  if (length(bad)) stop("알 수 없는 암종: ", paste(bad, collapse = ", "), " (cfg$cancers: ",
                        paste(cfg$cancers, collapse = ", "), ")")
  cfg$cancers <- args
  cat("명령줄 지정 암종만 실행:", paste(args, collapse = ", "), "\n")
}

set.seed(cfg$seed)
bp <-if (.Platform$OS.type == "windows") {
  SnowParam(workers = cfg$n_cores, RNGseed = cfg$seed)
} else {
  MulticoreParam(workers = cfg$n_cores, RNGseed = cfg$seed)
}

# ---- 유전자 세트 --------------------------------------------------------

# 기본 컬렉션 + cfg$gsea_extra_collections. subcollection이 여러 개면 처음 존재하는 것 사용
collection_specs <- c(list(
  Hallmark = list(collection = "H"),
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  KEGG     = list(collection = "C2", subcollection = c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS")),
  GOBP     = list(collection = "C5", subcollection = "GO:BP")
), cfg$gsea_extra_collections)

load_gene_sets <- function(spec) {
  subs <- spec$subcollection %||% list(NULL)
  for (s in subs) {
    df <- msigdbr::msigdbr(species = "Homo sapiens", collection = spec$collection, subcollection = s)
    if (nrow(df)) {
      sets <- lapply(split(df$gene_symbol, df$gs_name), unique)
      attr(sets, "source") <- paste(c(spec$collection, s), collapse = ":")
      return(sets)
    }
  }
  stop("MSigDB 컬렉션 없음: ", spec$collection, " ", paste(unlist(spec$subcollection), collapse = "/"))
}

cat("MSigDB 유전자 세트 불러오는 중...\n")
gene_sets <- lapply(collection_specs, load_gene_sets)
for (k in names(gene_sets)) cat(sprintf("  %-9s %-18s %5d sets\n", k, attr(gene_sets[[k]], "source"), length(gene_sets[[k]])))

# MSigDB .chip: 옛/별칭 기호(Probe Set ID) → 현재 MSigDB 기호(Gene Symbol)
chip_file <- list.files(cfg$clinical_dir, pattern = "\\.chip$", full.names = TRUE)[1]
if (is.na(chip_file)) stop("database/ 에 MSigDB .chip 파일 없음")
chip <- data.table::fread(chip_file, sep = "\t", select = 1:2, col.names = c("from", "to"))
chip <- chip[!is.na(to) & to != "" & !duplicated(from)]
chip_map <- setNames(chip$to, chip$from)

# ---- 함수 ---------------------------------------------------------------

genes_to_run <- function(all_genes) {
  g <- cfg$gsea_genes %||% cfg$primary_gene
  if (identical(tolower(g), "all")) return(all_genes)
  keys <- gene_key(g)
  miss <- setdiff(keys, all_genes)
  if (length(miss)) warning("GSEA 대상 유전자가 데이터에 없음: ", paste(toupper(miss), collapse = ", "), call. = FALSE)
  intersect(keys, all_genes)
}

# 설계 공변량 준비 (cfg$collapse_small_levels에 있는 변수는 소수 수준 병합)
prepare_coldata <- function(d, covs) {
  cd <- data.frame(row.names = d$sample_id, group = factor(d$group, levels = cfg$group_levels))
  for (v in covs) {
    x <- collapse_small(d[[v]], v)
    cd[[v]] <- if (is.numeric(x)) x else factor(x)
  }
  cd
}

run_deseq <- function(cnt, cd, covs) {
  design <- as.formula(paste("~", paste(c(covs, "group"), collapse = " + ")))
  dds <- DESeqDataSetFromMatrix(cnt, cd, design)
  keep <- rowSums(counts(dds) >= 10) >= min(table(cd$group))
  dds <- dds[keep, ]
  suppressMessages(DESeq(dds, parallel = TRUE, BPPARAM = bp))
}

# Wald stat → 결측 제거 → .chip 재매핑 → 중복 평균 → 내림차순
make_ranks <- function(res) {
  st <- res$stat
  names(st) <- rownames(res)
  st <- st[!is.na(st)]
  to <- chip_map[names(st)]
  nm <- ifelse(is.na(to), names(st), to)
  n_remapped <- sum(!is.na(to) & unname(to) != names(st))
  st <- vapply(split(st, nm), mean, numeric(1))
  out <- sort(st, decreasing = TRUE)
  attr(out, "n_remapped") <- n_remapped
  out
}

run_fgsea <- function(sets, ranks) {
  set.seed(cfg$seed)
  r <- suppressWarnings(fgseaMultilevel(sets, ranks, minSize = cfg$gsea_size[1], maxSize = cfg$gsea_size[2],
                                        eps = 0, BPPARAM = bp))
  r <- r[order(r$pval), ]
  data.frame(pathway = r$pathway, size = r$size, ES = r$ES, NES = r$NES, pval = r$pval, padj = r$padj,
             leading_edge = vapply(r$leadingEdge, paste, "", collapse = ";"))
}

pretty_pathway <- function(p, width = 70) {
  p <- sub("^(HALLMARK|REACTOME|KEGG_MEDICUS|KEGG|GOBP|GOCC|GOMF|WP|BIOCARTA|PID)_", "", p)
  p <- gsub("_", " ", p)
  ifelse(nchar(p) > width, paste0(substr(p, 1, width - 1), "…"), p)
}

norm_name <- function(x) gsub("[^A-Z0-9]", "", toupper(x))

# ---- 그림 ---------------------------------------------------------------

volcano_plot <- function(de, gene, cancer) {
  de$y   <- -log10(pmax(de$padj, 1e-300))
  de$dir <- ifelse(de$padj < 0.05 & de$log2FC_shrunk > 0, group_alt(),
                   ifelse(de$padj < 0.05 & de$log2FC_shrunk < 0, group_ref(), "NS"))
  lab <- head(de[order(de$padj), ], 15)
  lab <- unique(rbind(lab, de[de$symbol == gene, ]))
  ggplot(de[!is.na(de$padj), ], aes(log2FC_shrunk, y, colour = dir)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    ggrepel::geom_text_repel(data = lab, aes(label = symbol), size = 3, colour = "black",
                             max.overlaps = 30, min.segment.length = 0) +
    scale_colour_manual(values = c(cfg$group_colors[cfg$group_levels], NS = "grey75"),
                        labels = c(setNames(paste("Up in", cfg$group_levels), cfg$group_levels), NS = "NS"),
                        name = NULL) +
    labs(x = paste0("log2 fold change (", group_contrast_label(), ", apeglm-shrunken)"), y = expression(-log[10]~"adjusted p"),
         title = paste0(gene, " ", group_contrast_label(), " — TCGA-", cancer)) +
    theme_classic(base_size = 11) + theme(legend.position = "top")
}

pca_plot <- function(dds, gene, cancer) {
  vsd <- vst(dds, blind = TRUE)
  extra <- intersect(cfg$gsea_design_covariates, names(colData(dds)))[1]   # 두 번째 패널 색 (기본 tss)
  # intgroup이 2개 이상이면 plotPCA가 group 열을 상호작용("Low:AF")으로 덮어쓰므로 group만 넘기고
  # 두 번째 색 변수는 colData에서 직접 가져옴
  pc  <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
  if (!is.na(extra)) pc[[extra]] <- colData(vsd)[[extra]]
  pv  <- round(100 * attr(pc, "percentVar"))
  base <- ggplot(pc, aes(PC1, PC2)) +
    labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) +
    theme_classic(base_size = 10)
  p1 <- base + geom_point(aes(colour = group), size = 1.2, alpha = 0.8) +
    scale_colour_manual(values = cfg$group_colors, name = paste(gene, "expression"))
  if (is.na(extra)) return(p1 + labs(title = paste0("PCA (vst) — TCGA-", cancer)))
  p2 <- base + geom_point(aes(colour = .data[[extra]]), size = 1.2, alpha = 0.8) +
    labs(colour = if (extra %in% names(var_labels)) var_labels[[extra]] else extra)
  (p1 | p2) + plot_annotation(title = paste0("PCA of variance-stabilized expression — TCGA-", cancer))
}

nes_plot <- function(tab, coll, gene, cancer) {
  sig <- tab[!is.na(tab$padj) & tab$padj < 0.05, ]
  if (nrow(sig)) {
    top <- rbind(head(sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ], 10),
                 head(sig[sig$NES < 0, ][order(sig$NES[sig$NES < 0]), ], 10))
    sub <- paste0(nrow(sig), " pathways with padj < 0.05; top 10 up / down by NES")
  } else {
    top <- rbind(head(tab[order(-tab$NES), ], 10), head(tab[order(tab$NES), ], 10))
    top <- top[!duplicated(top$pathway), ]
    sub <- "No pathway with padj < 0.05; top 10 up / down by NES shown"
  }
  top$label <- pretty_pathway(top$pathway)
  top$dir   <- ifelse(top$NES > 0, group_alt(), group_ref())
  top$padj_txt <- ifelse(top$padj < 0.001, "<0.001", sprintf("%.3f", top$padj))
  ggplot(top, aes(NES, reorder(label, NES), fill = dir)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, colour = "grey30") +
    geom_text(aes(label = padj_txt, hjust = ifelse(NES > 0, -0.15, 1.15)), size = 2.7) +
    scale_fill_manual(values = cfg$group_colors,
                      labels = setNames(paste("Enriched in", gene, cfg$group_levels), cfg$group_levels), name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.18)) +
    labs(x = "Normalized enrichment score (NES); label = adjusted p", y = NULL,
         title = paste0(coll, " — ", gene, " ", group_contrast_label(), ", TCGA-", cancer), subtitle = sub) +
    theme_classic(base_size = 9) + theme(legend.position = "top")
}

enrichment_plot <- function(genes, ranks, row, coll) {
  d <- plotEnrichmentData(genes, ranks)
  col <- if (row$NES > 0) cfg$group_colors[[group_alt()]] else cfg$group_colors[[group_ref()]]
  tick_h <- d$spreadES / 16
  ggplot(d$curve, aes(rank, ES)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(d$posES, d$negES), colour = "grey60", linetype = "dashed") +
    geom_segment(data = d$ticks, aes(x = rank, xend = rank, y = -tick_h, yend = tick_h),
                 inherit.aes = FALSE, linewidth = 0.2) +
    geom_line(colour = col, linewidth = 0.9) +
    labs(x = paste0("Gene rank (left: up in ", group_alt(), ", right: up in ", group_ref(), ")"), y = "Enrichment score",
         title = pretty_pathway(row$pathway, 90),
         subtitle = sprintf("%s | NES = %.2f, padj = %s, size = %d", coll, row$NES, fmt_p(row$padj), row$size)) +
    theme_classic(base_size = 10)
}

# ---- 유전자 1개 ---------------------------------------------------------

analyze_gene <- function(g, cancer, cnt, merged, gdc_expr) {
  gene <- toupper(g)
  sub  <- file.path("gsea", gene)
  files <- character()
  covs <- intersect(cfg$gsea_design_covariates, names(merged))
  on.exit(gc(), add = TRUE)   # 유전자 사이 메모리 반환 (dds 객체가 큼)
  if (length(setdiff(cfg$gsea_design_covariates, covs))) {
    cat("  설계 공변량 없음 (제외):", paste(setdiff(cfg$gsea_design_covariates, covs), collapse = ", "), "\n")
  }

  # ---- 환자 매칭 ----
  d <- merged
  d$group <- d[[paste0(g, "_group")]]
  d <- d[!is.na(d$group), ]
  d <- d[stats::complete.cases(d[, covs, drop = FALSE]), ]
  ids <- intersect(d$sample_id, colnames(cnt))
  only_merged <- setdiff(d$sample_id, colnames(cnt))
  only_counts <- setdiff(colnames(cnt), merged$sample_id)
  cat(sprintf("  %s: 매칭 %d명 | 임상에만 %d명%s | counts에만 %d명 (신보조요법 제외/임상 없음)\n",
              gene, length(ids), length(only_merged),
              if (length(only_merged)) paste0(" (", paste(head(only_merged, 5), collapse = ", "), ")") else "",
              length(only_counts)))
  d  <- d[match(ids, d$sample_id), ]
  cd <- prepare_coldata(d, covs)
  n_grp <- table(cd$group)
  cat(sprintf("           그룹 %s %d / %s %d%s\n", group_ref(), n_grp[[group_ref()]], group_alt(), n_grp[[group_alt()]],
              paste0(vapply(covs, function(v) sprintf(", %s %d개 수준", v, length(unique(cd[[v]]))), ""), collapse = "")))

  # ---- QC: GDC 발현 vs .sav 발현 ----
  qc <- data.frame(gene = gene, cancer = cancer, n_matched = length(ids),
                   n_only_clinical = length(only_merged), n_only_counts = length(only_counts),
                   n_low = n_grp[[group_ref()]], n_high = n_grp[[group_alt()]])
  if (!is.null(gdc_expr) && gene %in% rownames(gdc_expr)) {
    x_gdc <- gdc_expr[gene, ids]
    x_sav <- d[[paste0(g, "_expression")]]
    r <- suppressWarnings(cor(x_gdc, x_sav, method = "spearman", use = "complete.obs"))
    split_gdc <- median_split(x_gdc)
    agree <- 100 * mean(as.character(split_gdc) == as.character(d$group), na.rm = TRUE)
    qc$spearman_gdc_vs_sav <- round(r, 4)
    qc$median_split_agreement_pct <- round(agree, 1)
    cat(sprintf("           QC: GDC log2(TPM+1) vs .sav 발현 Spearman r = %.3f, median split 일치 %.1f%%\n", r, agree))
    if (r < 0.8) warning(cancer, " ", gene, ": GDC 발현과 .sav 발현 상관 낮음 (r = ", round(r, 3),
                         ") → 서로 다른 데이터일 수 있음", call. = FALSE)
  } else {
    qc$spearman_gdc_vs_sav <- NA
    qc$median_split_agreement_pct <- NA
    cat("           QC: GDC 발현 행렬에 유전자 없음 → 상관 확인 생략\n")
  }

  # ---- DESeq2: 주 설계 + 공변량 없는 민감도 ----
  t0 <- Sys.time()
  dds <- run_deseq(cnt[, ids], cd, covs)
  res <- results(dds, contrast = c("group", group_alt(), group_ref()), parallel = TRUE, BPPARAM = bp)
  # apeglm은 일부 유전자에서 최적화 경고를 반복 출력 → 억제하고 개수만 보고 (결과에는 영향 없음)
  n_apeglm_warn <- 0
  shr <- withCallingHandlers(
    lfcShrink(dds, coef = paste0("group_", group_alt(), "_vs_", group_ref()), type = "apeglm",
              parallel = TRUE, BPPARAM = bp, quiet = TRUE),
    warning = function(w) { n_apeglm_warn <<- n_apeglm_warn + 1; invokeRestart("muffleWarning") })
  cat(sprintf("           DESeq2 (~ %s): %d개 유전자, %.1f분%s\n", paste(c(covs, "group"), collapse = " + "),
              nrow(dds), as.numeric(difftime(Sys.time(), t0, units = "mins")),
              if (n_apeglm_warn) sprintf(" (apeglm 경고 %d건 억제)", n_apeglm_warn) else ""))
  saveRDS(dds, processed_path(cancer, paste0("dds_", gene, ".rds")))

  de <- data.frame(symbol = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange,
                   lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                   log2FC_shrunk = shr$log2FoldChange)
  de <- de[order(de$pvalue), ]
  files <- c(files, save_table(de, cancer, "de_results.csv", sub))

  # PCA는 dds가 필요하므로 여기서 그리고 dds/vst를 바로 해제 (메모리)
  files <- c(files, save_fig(pca_plot(dds, gene, cancer), cancer, "pca", 12, 5, subdir = sub))
  rm(dds, shr)
  gc(verbose = FALSE)

  # 기준 유전자 자체는 High에서 강하게 증가해야 함
  rk  <- rank(-de$stat, ties.method = "min", na.last = "keep")
  i   <- match(gene, de$symbol)
  pct <- 100 * rk[i] / sum(!is.na(de$stat))
  cat(sprintf("           확인: %s 순위 %d / %d (상위 %.2f%%), log2FC %.2f (shrunken %.2f), padj %s\n",
              gene, rk[i], sum(!is.na(de$stat)), pct, de$log2FC[i], de$log2FC_shrunk[i], fmt_p(de$padj[i])))
  if (is.na(pct) || pct > 1) warning(cancer, " ", gene, ": 기준 유전자가 상위 1%에 없음 → 그룹/발현 확인 필요",
                                     call. = FALSE)
  qc$gene_rank <- rk[i]
  qc$gene_rank_pct <- round(pct, 3)
  qc$gene_log2FC <- round(de$log2FC[i], 3)
  qc$n_genes_tested <- nrow(de)
  qc$n_de_padj05 <- sum(de$padj < 0.05, na.rm = TRUE)

  t0 <- Sys.time()
  covs0 <- intersect(cfg$gsea_sensitivity_covariates, names(cd))
  dds0 <- run_deseq(cnt[, ids], cd[c(covs0, "group")], covs0)
  res0 <- results(dds0, contrast = c("group", group_alt(), group_ref()), parallel = TRUE, BPPARAM = bp)
  cat(sprintf("           DESeq2 민감도 (~ %s): %.1f분\n", paste(c(covs0, "group"), collapse = " + "),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(dds0)
  gc(verbose = FALSE)

  # ---- 순위 ----
  ranks  <- make_ranks(res)
  ranks0 <- make_ranks(res0)
  qc$n_ranked <- length(ranks)
  qc$n_symbols_remapped <- attr(ranks, "n_remapped")
  cat(sprintf("           순위: %d개 유전자 (MSigDB .chip 재매핑 %d개)\n", length(ranks), attr(ranks, "n_remapped")))

  # ---- GSEA (컬렉션별) ----
  all_res <- list()
  for (coll in names(gene_sets)) {
    t0 <- Sys.time()
    tab  <- run_fgsea(gene_sets[[coll]], ranks)
    tab0 <- run_fgsea(gene_sets[[coll]], ranks0)
    tab$NES_no_cov  <- tab0$NES[match(tab$pathway, tab0$pathway)]
    tab$padj_no_cov <- tab0$padj[match(tab$pathway, tab0$pathway)]
    tab$robust_no_covariate <- ifelse(tab$padj < 0.05,
                                      sign(tab$NES_no_cov) == sign(tab$NES) & tab$padj_no_cov < 0.05, NA)
    tab$collapsed_main <- NA
    if (coll == "GOBP") {
      sig <- head(tab[tab$padj < 0.05, ], 300)   # 계산량 제한: 유의 경로 상위 300개
      if (nrow(sig) > 1) {
        fr <- data.table::as.data.table(sig)
        fr$leadingEdge <- strsplit(sig$leading_edge, ";")
        set.seed(cfg$seed)
        coll_res <- suppressWarnings(collapsePathways(fr, gene_sets[[coll]], ranks))
        tab$collapsed_main[tab$pathway %in% sig$pathway] <- tab$pathway[tab$pathway %in% sig$pathway] %in%
          coll_res$mainPathways
      }
    }
    n_sig <- sum(tab$padj < 0.05, na.rm = TRUE)
    cat(sprintf("           GSEA %-9s padj<0.05 %4d개 (공변량 없이도 유지 %d), %.1f분\n", coll, n_sig,
                sum(tab$robust_no_covariate %in% TRUE),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    files <- c(files, save_table(tab, cancer, paste0("gsea_", coll, ".csv"), sub))
    if (coll == "Hallmark" || n_sig > 0) {
      nt <- min(20, if (n_sig) n_sig else 20)
      files <- c(files, save_fig(nes_plot(tab, coll, gene, cancer), cancer, paste0("nes_", coll),
                                 width = 9, height = 1.8 + 0.27 * nt, subdir = sub))
    }
    all_res[[coll]] <- data.frame(collection = coll, tab)
  }
  allr <- do.call(rbind, all_res)

  # ---- 그림: volcano, enrichment ----
  files <- c(files, save_fig(volcano_plot(de, gene, cancer), cancer, "volcano", 7, 6.5, subdir = sub))

  # enrichment plot 대상: Hallmark 유의 상위 N (padj 순) + gsea_highlight + C8 중 gsea_c8_pattern 유의 상위 N
  n_top <- cfg$gsea_enrichment_top
  top_sig <- function(t) {
    s <- t[!is.na(t$padj) & t$padj < 0.05, ]
    head(s$pathway[order(s$padj, -abs(s$NES))], n_top)
  }
  top_h  <- if (!is.null(all_res$Hallmark)) top_sig(all_res$Hallmark) else character()
  top_c8 <- if (!is.null(all_res$C8) && nzchar(cfg$gsea_c8_pattern %||% "")) {
    top_sig(all_res$C8[grepl(cfg$gsea_c8_pattern, all_res$C8$pathway), ])
  } else character()
  hl <- allr$pathway[norm_name(allr$pathway) %in% norm_name(cfg$gsea_highlight)]
  miss_hl <- cfg$gsea_highlight[!norm_name(cfg$gsea_highlight) %in% norm_name(allr$pathway)]
  if (length(miss_hl)) cat("           강조 경로 중 결과에 없음:", paste(miss_hl, collapse = ", "), "\n")
  cat(sprintf("           enrichment plot: Hallmark %d, 강조 %d, C8 (%s) %d\n",
              length(top_h), length(hl), cfg$gsea_c8_pattern, length(top_c8)))
  for (pw in unique(c(top_h, hl, top_c8))) {
    row  <- allr[allr$pathway == pw, ][1, ]
    genes <- gene_sets[[row$collection]][[pw]]
    files <- c(files, save_fig(enrichment_plot(genes, ranks, row, row$collection), cancer,
                               paste0("enrichment_", pw), width = 7, height = 4.5, subdir = sub))
  }

  # ---- 요약 표 ----
  # Hallmark: 유의 경로 중 NES 상위 5 (High 쪽) + 하위 5 (Low 쪽); 나머지 컬렉션: padj 상위 10
  hallmark_updown <- function(t, n = 5) {
    s  <- t[!is.na(t$padj) & t$padj < 0.05, ]
    up <- s[s$NES > 0, ]
    dn <- s[s$NES < 0, ]
    rbind(head(up[order(-up$NES), ], n), head(dn[order(dn$NES), ], n))
  }
  summ <- do.call(rbind, lapply(names(gene_sets), function(coll) {
    t <- all_res[[coll]]
    t <- if (coll == "Hallmark") hallmark_updown(t) else head(t[t$padj < 0.05, ], 10)
    if (!nrow(t)) return(NULL)
    le <- vapply(strsplit(t$leading_edge, ";"), function(x) paste(head(x, 10), collapse = ", "), "")
    data.frame(Collection = coll, Pathway = t$pathway, Size = t$size, NES = sprintf("%.2f", t$NES),
               padj = fmt_p(t$padj), `Robust (no covariate)` = t$robust_no_covariate,
               `Leading edge (top 10)` = le, check.names = FALSE)
  }))
  if (is.null(summ)) summ <- data.frame(Collection = names(gene_sets), Pathway = "No pathway with padj < 0.05")
  files <- c(files, save_df_table(
    summ, cancer, "gsea_summary", landscape = TRUE, font_size = 7, subdir = sub,
    caption = paste0("GSEA of ", gene, " ", group_contrast_label(), " expression, TCGA-", cancer,
                     " (DESeq2 Wald statistic, design ~ ", paste(c(covs, "group"), collapse = " + "),
                     "). Hallmark: top 5 up and top 5 down by NES; other collections: top 10 by padj; ",
                     "all padj < 0.05. Positive NES = enriched in ", group_alt(), ". ",
                     "Robust = same direction and padj < 0.05 with design ~ ",
                     paste(c(cfg$gsea_sensitivity_covariates, "group"), collapse = " + "), ".")))
  files <- c(files, save_table(qc, cancer, "gsea_QC.csv", sub))
  move_stale_outputs(cancer, NULL, files, subdir = sub)

  list(files = files,
       summary = list(gene = gene,
                      n_sig = vapply(all_res, function(t) sum(t$padj < 0.05, na.rm = TRUE), numeric(1)),
                      n_robust = vapply(all_res, function(t) sum(t$robust_no_covariate %in% TRUE), numeric(1)),
                      rank = qc$gene_rank, n_tested = qc$n_genes_tested,
                      hallmark_up = head(all_res$Hallmark[order(-all_res$Hallmark$NES), ], 5),
                      hallmark_down = head(all_res$Hallmark[order(all_res$Hallmark$NES), ], 5)))
}

# ---- 실행 ---------------------------------------------------------------

created  <- character()
failures <- character()
console  <- list()

for (cancer in cfg$cancers) {
  cat("\n==========", cancer, "==========\n")
  f_cnt <- processed_path(cancer, "counts.rds")
  if (!file.exists(f_cnt)) {
    stop(cancer, ": ", f_cnt, " 없음 → 먼저 Rscript 02a_download_tcga.R 을 실행하세요")
  }
  cnt    <- readRDS(f_cnt)
  merged <- readRDS(processed_path(cancer, "merged.rds"))
  genes  <- genes_to_run(detect_genes(merged))
  cat("counts:", nrow(cnt), "유전자 ×", ncol(cnt), "명 / GSEA 대상:", paste(toupper(genes), collapse = ", "), "\n")

  # QC용 GDC log2(TPM+1) 행렬에서 대상 유전자 행만
  f_tpm <- file.path(cfg$clinical_dir, paste0("TCGA_", cancer, "_RNAseq_Expression.csv"))
  gdc_expr <- NULL
  if (file.exists(f_tpm)) {
    m <- data.table::fread(f_tpm)
    m <- m[m$Gene %in% toupper(genes), ]
    gdc_expr <- as.matrix(m[, -1])
    rownames(gdc_expr) <- m$Gene
  }

  results <- list()
  for (g in genes) {
    r <- tryCatch(analyze_gene(g, cancer, cnt, merged, gdc_expr), error = function(e) {
      failures <<- c(failures, paste0(cancer, " ", toupper(g), ": ", conditionMessage(e)))
      cat("  !!", toupper(g), "실패:", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(r)) {
      created <- c(created, r$files)
      results[[g]] <- r$summary
    }
  }
  # 이번 실행 대상이 아닌 유전자의 이전 GSEA 폴더 → gsea/_stale/
  move_stale_dirs(cancer, "gsea", toupper(genes))
  console[[cancer]] <- results
  rm(cnt, merged, gdc_expr)
  gc(verbose = FALSE)
}

bpstop(bp)

cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")

cat("\n===== 요약 =====\n")
for (cancer in names(console)) {
  for (r in console[[cancer]]) {
    cat(cancer, r$gene, sprintf("(기준 유전자 순위 %d / %d)\n", r$rank, r$n_tested))
    cat("  padj < 0.05 경로 수 (공변량 없이도 유지):",
        paste(sprintf("%s %d (%d)", names(r$n_sig), r$n_sig, r$n_robust), collapse = ", "), "\n")
    for (k in c("up", "down")) {
      h <- r[[paste0("hallmark_", k)]]
      cat(sprintf("  Hallmark NES %s 5 (%s):\n", if (k == "up") "상위" else "하위",
                  if (k == "up") paste(group_alt(), "쪽") else paste(group_ref(), "쪽")))
      cat(sprintf("    %-45s NES %5.2f  padj %s\n", h$pathway, h$NES, fmt_p(h$padj)), sep = "")
    }
  }
}
if (length(failures)) {
  cat("\n실패한 분석 (", length(failures), "개):\n", paste0("  ", failures, collapse = "\n"), "\n", sep = "")
} else {
  cat("\n실패한 분석 없음\n")
}
