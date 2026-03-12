## =========================================================
## 0. 当前项目专用脚本：VTE + PE
## 说明：
##   - 这是一个单文件脚本
##   - 不再 source 其他 core.R
##   - 只跑两个 outcome：
##       1) finn-b-I9_VTE
##       2) finn-b-I9_PULMEMB
##   - 先把“项目参数区”放最上面
##   - 接下来你把原来 core 的全部函数内容贴到下面指定位置
## =========================================================


## =========================================================
## 1. 基础设置
## =========================================================
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)


## =========================================================
## 2. 当前项目输入区（以后主要改这里）
## =========================================================
PROJECT_NAME <- "Venous_thromboembolism_PE_VTE_only"

FIXED_IDS <- c(
  "finn-b-I9_VTE",
  "finn-b-I9_PULMEMB"
)

PIPELINE_ROOT <- "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE"

TRAIT_REGEX <- NULL
ID_REGEX <- NULL
POPULATION_REGEX <- "European|EUR|UKB|UK Biobank|Finnish|FinnGen|Finland"

REBUILD_TOP3 <- TRUE
REBUILD_MR   <- TRUE
REBUILD_AO   <- FALSE

DO_FOREST_PLOT <- TRUE
DO_QC_TOP      <- TRUE
VERBOSE_RUN    <- TRUE

N_SEARCH <- length(FIXED_IDS)


## =========================================================
## 3. 下面请粘贴你原来 core 的全部函数内容
## 粘贴范围：
##   从 slugify <- function(x) {
##   一直到 run_project <- function(...) { ... } 结束
## =========================================================



## ======= 你把原来 core 的内容粘贴在这里下面 =======



## =========================================================
## 4. 先不要在这里写执行代码
## 等你把 core 贴完后，再发给我
## 我再给你最后一段“主执行区”
## =========================================================
slugify <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

calc_sample_size_num <- function(dt) {
  dt <- as.data.table(copy(dt))
  dt[, sample_size_num := suppressWarnings(as.numeric(sample_size))]
  if (all(is.na(dt$sample_size_num)) && all(c("ncase", "ncontrol") %in% names(dt))) {
    dt[, sample_size_num := suppressWarnings(as.numeric(ncase) + as.numeric(ncontrol))]
  }
  dt
}

resolve_outcomes_meta <- function(
    ao,
    keyword = NULL,
    fixed_ids = character(),
    mode = c("mixed", "fixed", "search"),
    n = 3,
    population_regex = "European|EUR|UKB|UK Biobank|Finnish|FinnGen|Finland",
    trait_regex = NULL,
    id_regex = NULL,
    min_year = NA_integer_,
    cache_tsv = NULL,
    rebuild = FALSE,
    verbose = TRUE
) {
  mode <- match.arg(mode)
  
  # 目的：有缓存且不重建 -> 直接读缓存
  if (!is.null(cache_tsv) && file.exists(cache_tsv) && !rebuild) {
    if (verbose) message("[CACHE] loaded selected outcomes: ", cache_tsv)
    out <- fread(cache_tsv)
    return(unique(out, by = "id"))
  }
  
  ao_dt <- as.data.table(ao)
  ao_dt <- ao_dt[!is.na(id)]
  ao_dt <- unique(ao_dt, by = "id")
  
  fixed_ids <- unique(fixed_ids[nzchar(fixed_ids)])
  fixed_hit <- ao_dt[id %in% fixed_ids]
  
  # 目的：固定模式仅使用 fixed_ids
  if (mode == "fixed") {
    if (length(fixed_ids) > 0 && nrow(fixed_hit) == 0) {
      stop("mode='fixed' 但 fixed_ids 在 available_outcomes 中都未命中，请检查 id。")
    }
    fixed_hit <- calc_sample_size_num(fixed_hit)
    if (!is.null(cache_tsv)) fwrite(fixed_hit, cache_tsv, sep = "\t")
    return(fixed_hit)
  }
  
  # 目的：search/mixed 模式需要 keyword；若缺失则 mixed 退化为 fixed
  if (is.null(keyword) || !nzchar(keyword)) {
    if (mode == "search") stop("mode='search' 必须提供 keyword。")
    fixed_hit <- calc_sample_size_num(fixed_hit)
    if (!is.null(cache_tsv)) fwrite(fixed_hit, cache_tsv, sep = "\t")
    return(fixed_hit)
  }
  
  # 目的：按 trait keyword 找候选
  meta <- ao_dt[!is.na(trait) & grepl(keyword, trait, ignore.case = TRUE)]
  
  # 目的：人群过滤（列存在才用）
  if ("population" %in% names(meta) && !is.null(population_regex) && nzchar(population_regex)) {
    meta <- meta[grepl(population_regex, population, ignore.case = TRUE)]
  }
  
  # 目的：trait 二次过滤（更精准）
  if (!is.null(trait_regex) && nzchar(trait_regex)) {
    meta <- meta[grepl(trait_regex, trait, ignore.case = TRUE)]
  }
  
  # 目的：按 id 限定来源（如 ^ukb- / ^finn-）
  if (!is.null(id_regex) && nzchar(id_regex)) {
    meta <- meta[grepl(id_regex, id, ignore.case = TRUE)]
  }
  
  # 目的：限定年份（列存在才用）
  if (!is.na(min_year) && "year" %in% names(meta)) {
    yr <- suppressWarnings(as.integer(meta$year))
    meta <- meta[!is.na(yr) & yr >= min_year]
  }
  
  if (nrow(meta) == 0) {
    if (mode == "search") {
      stop(
        "[STOP] No OpenGWAS outcomes matched keyword = ", keyword,
        "\n建议：换同义词、放宽 population/trait 过滤、或先不限定 EUR。"
      )
    }
    fixed_hit <- calc_sample_size_num(fixed_hit)
    if (!is.null(cache_tsv)) fwrite(fixed_hit, cache_tsv, sep = "\t")
    return(fixed_hit)
  }
  
  # 目的：按样本量排序优先大样本 outcome
  meta <- calc_sample_size_num(meta)
  setorder(meta, -sample_size_num, trait)
  meta <- unique(meta, by = "id")
  
  if (mode == "search") {
    picked <- meta[1:min(n, .N)]
    if (!is.null(cache_tsv)) fwrite(picked, cache_tsv, sep = "\t")
    return(picked)
  }
  
  # 目的：mixed：固定保留 + 搜索补齐到 n
  meta2 <- meta[!id %in% fixed_ids]
  need <- max(0L, n - nrow(fixed_hit))
  picked_search <- if (need > 0L) meta2[1:min(need, .N)] else meta2[0]
  
  out <- rbindlist(list(fixed_hit, picked_search), fill = TRUE)
  out <- unique(out, by = "id")
  out <- calc_sample_size_num(out)
  setorder(out, -sample_size_num, trait)
  
  if (!is.null(cache_tsv)) fwrite(out, cache_tsv, sep = "\t")
  out
}

build_gene_map <- function(vpub_file, vpub_sheet = "S3") {
  stopifnot(file.exists(vpub_file))
  s3 <- as.data.table(readxl::read_xlsx(vpub_file, sheet = vpub_sheet, skip = 1))
  if (nrow(s3) >= 1) s3 <- s3[-.N]
  
  gene_map <- s3[, .(
    GENE_ID = trimws(as.character(geneID)),
    gene_symbol = trimws(as.character(geneSymble))
  )]
  
  gene_map <- gene_map[between(nchar(GENE_ID), 1L, 1e9)]
  gene_map[, GENE_ID := sub("\\.\\d+$", "", GENE_ID)]
  
  gene_map[is.na(gene_symbol), gene_symbol := ""]
  gene_map[, sym_len := nchar(gene_symbol)]
  setorder(gene_map, GENE_ID, -sym_len)
  gene_map <- unique(gene_map, by = "GENE_ID")
  gene_map[, sym_len := NULL]
  gene_map
}

make_forest_plot <- function(key_table, top_table, out_dir, theme_pub) {
  plot_use <- if (!is.null(key_table) && nrow(key_table) > 0) "significant" else "tophits"
  plot_dat <- if (plot_use == "significant") copy(as.data.table(key_table)) else
    copy(as.data.table(top_table))
  
  if (!"gene_symbol" %in% names(plot_dat)) plot_dat[, gene_symbol := NA_character_]
  
  need_cols <- c("CELL_TYPE", "GENE_ID", "method", "nsnp", "pval", "fdr", "OR", "OR_lci95", "OR_uci95")
  miss_cols <- setdiff(need_cols, names(plot_dat))
  if (length(miss_cols) > 0) {
    stop("Forest plot missing columns: ", paste(miss_cols, collapse = ", "),
         "\n请确认主表已生成 OR/CI，且包含这些列。")
  }
  
  num_cols <- intersect(c("OR", "OR_lci95", "OR_uci95", "pval", "fdr", "nsnp"), names(plot_dat))
  plot_dat[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]
  
  if (!"outcome" %in% names(plot_dat)) plot_dat[, outcome := NA_character_]
  if (!"id.outcome" %in% names(plot_dat)) plot_dat[, id.outcome := "outcome"]
  plot_dat[, outcome_show := fifelse(!is.na(outcome) & nzchar(outcome), outcome, id.outcome)]
  
  plot_dat <- plot_dat[
    is.finite(OR) & is.finite(OR_lci95) & is.finite(OR_uci95) &
      OR > 0 & OR_lci95 > 0 & OR_uci95 > 0
  ]
  if (nrow(plot_dat) == 0) stop("No valid rows left for forest plot after filtering OR/CI.")
  
  plot_dat[, label := ifelse(
    !is.na(gene_symbol) & nzchar(gene_symbol),
    paste0(CELL_TYPE, " | ", gene_symbol, " (", GENE_ID, ")"),
    paste0(CELL_TYPE, " | ", GENE_ID)
  )]
  
  plot_dat[, ann := paste0(
    method,
    " | nsnp=", nsnp,
    " | p=", formatC(pval, format = "e", digits = 2),
    " | FDR=", formatC(fdr, format = "e", digits = 2)
  )]
  
  plot_dat[, y_fac := paste0(label, "@@", outcome_show)]
  setorder(plot_dat, outcome_show, pval)
  lev <- plot_dat[, unique(y_fac)]
  plot_dat[, y_fac := factor(y_fac, levels = rev(lev))]
  
  plot_dat[, x_text := pmax(max(OR_uci95, na.rm = TRUE), 1) * 1.5, by = outcome_show]
  
  out_plot_data <- file.path(out_dir, "Fig_MR_forest_data.tsv")
  fwrite(plot_dat, out_plot_data, sep = "\t")
  
  p_forest <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = OR, y = y_fac, xmin = OR_lci95, xmax = OR_uci95, color = method)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    ggplot2::geom_pointrange(fatten = 1.15) +
    ggplot2::geom_text(ggplot2::aes(x = x_text, label = ann),
                       color = "black", size = 2, hjust = 0) +
    ggplot2::scale_x_log10() +
    scico::scale_color_scico_d(palette = "batlow", begin = 0.1, end = 0.9) +
    ggplot2::scale_y_discrete(labels = function(x) sub("@@.*$", "", x)) +
    ggplot2::labs(
      x = "Odds ratio (log10 scale)",
      y = NULL,
      color = "MR method",
      title = "MR forest plot (cell type | gene → outcome)"
    ) +
    theme_pub(base = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 50, 5.5, 5.5),
      axis.text.y = ggplot2::element_text(size = 6)   # <- 加这一行
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::facet_wrap(~ outcome_show, scales = "free_y", ncol = 1)
  
  out_png <- file.path(out_dir, "Fig_MR_forest.png")
  h <- min(26, max(7, 0.25 * nrow(plot_dat)))
  ggplot2::ggsave(out_png, p_forest, width = 14, height = h, dpi = 300)
  
  list(plot = p_forest, data = plot_dat, plot_use = plot_use, out_png = out_png, out_data = out_plot_data)
}

qc_top_table <- function(top_table, out_dir) {
  top_dt <- as.data.table(copy(top_table))
  if (!"gene_symbol" %in% names(top_dt)) top_dt[, gene_symbol := NA_character_]
  
  top_dt[, GENE_ID := sub("\\.\\d+$", "", trimws(GENE_ID))]
  top_dt[, CELL_TYPE := trimws(CELL_TYPE)]
  
  if (!"outcome" %in% names(top_dt)) top_dt[, outcome := NA_character_]
  if (!"id.outcome" %in% names(top_dt)) top_dt[, id.outcome := NA_character_]
  top_dt[, outcome_show := fifelse(!is.na(outcome) & nzchar(outcome), outcome, id.outcome)]
  
  top_dt[, gene_show := fifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, GENE_ID)]
  
  gene_cell <- unique(top_dt[, .(GENE_ID, gene_show, CELL_TYPE)])
  gene_cell_stat <- gene_cell[, .(
    n_celltypes = uniqueN(CELL_TYPE),
    celltypes = paste(sort(unique(CELL_TYPE)), collapse = "; ")
  ), by = .(GENE_ID, gene_show)]
  setorder(gene_cell_stat, -n_celltypes, gene_show)
  
  gene_out <- unique(top_dt[, .(GENE_ID, gene_show, outcome_show)])
  gene_out_stat <- gene_out[, .(
    n_outcomes = uniqueN(outcome_show),
    outcomes = paste(sort(unique(outcome_show)), collapse = "; ")
  ), by = .(GENE_ID, gene_show)]
  setorder(gene_out_stat, -n_outcomes, gene_show)
  
  gene_joint <- merge(gene_cell_stat, gene_out_stat, by = c("GENE_ID", "gene_show"), all = TRUE)
  setorder(gene_joint, -n_outcomes, -n_celltypes, gene_show)
  
  out_joint <- file.path(out_dir, "QC_topTable_gene_across_celltypes_outcomes.tsv")
  fwrite(gene_joint, out_joint, sep = "\t")
  
  ct_stat <- top_dt[, .(n_hits = .N, n_genes = uniqueN(GENE_ID)), by = .(CELL_TYPE)][order(-n_hits)]
  out_ct <- file.path(out_dir, "QC_topTable_hits_by_celltype.tsv")
  fwrite(ct_stat, out_ct, sep = "\t")
  
  out_stat <- top_dt[, .(n_hits = .N, n_genes = uniqueN(GENE_ID)), by = .(outcome_show)][order(-n_hits)]
  out_outcome <- file.path(out_dir, "QC_topTable_hits_by_outcome.tsv")
  fwrite(out_stat, out_outcome, sep = "\t")
  
  ct_out <- top_dt[, .(n_hits = .N, n_genes = uniqueN(GENE_ID)),
                   by = .(outcome_show, CELL_TYPE)][order(outcome_show, -n_hits)]
  out_ct_out <- file.path(out_dir, "QC_topTable_hits_by_celltype_by_outcome.tsv")
  fwrite(ct_out, out_ct_out, sep = "\t")
  
  invisible(list(
    gene_joint = gene_joint,
    ct_stat = ct_stat,
    out_stat = out_stat,
    ct_out = ct_out,
    out_joint = out_joint,
    out_ct = out_ct,
    out_outcome = out_outcome,
    out_ct_out = out_ct_out
  ))
}

# ============================================================
# 主函数：run_project()
# 目的：你外部只改 keyword 或 fixed_ids 就能完整跑完并出表/出图
# ============================================================
run_project <- function(
    keyword,
    fixed_ids = character(),     # 指定 outcome id（不填则走搜索）
    use_search_fill = FALSE,     # TRUE：固定 + 搜索补齐；FALSE：只跑 fixed_ids
    n_search = 3,                # 搜索/补齐时取几个
    trait_regex = NULL,
    id_regex = NULL,
    population_regex = "European|EUR|UKB|UK Biobank|Finnish|FinnGen|Finland",
    min_year = NA_integer_,
    
    # 路径参数（你一般不用改，放这里做默认值）
    pipeline_root = "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE",
    exp_rds = "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE/tmp_onek1k/vpub_s3_p5e8_kb10000_r2_0.001/mr_n982_from_s3_F10.rds",
    theme_file = "W:/fangxiaobin-CTRH-omics/000temple/scico_plot_theme.R",
    ao_cache_rds = "D:/fangxiaobin-CTRH-omics/disease_mr_keygene/Reference/MR/cache/available_outcomes.rds",
    vpub_file  = "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE/literature/V_publish.xlsx",
    vpub_sheet = "S3",
    
    # 分析参数
    fdr_cut = 0.05,
    top_n_per_outcome = 30,
    
    # 重跑开关
    rebuild_top3 = FALSE,
    rebuild_mr   = FALSE,
    rebuild_ao   = FALSE,
    
    # 其它
    sleep = 10,
    do_forest_plot = TRUE,
    do_qc_top = TRUE,
    restore_wd = FALSE,
    verbose = TRUE
) {
  suppressPackageStartupMessages({
    library(data.table)
    library(TwoSampleMR)
    library(ieugwasr)
    library(ggplot2)
    library(scico)
    library(patchwork)
    library(readxl)
  })
  
  # 目的：可选关闭旧图形设备，避免 ggsave/绘图冲突
 # while (length(grDevices::dev.list()) > 0) grDevices::dev.off()
  gc()
  
  fdr_tag <- gsub("\\.", "p", format(fdr_cut, scientific = FALSE, trim = TRUE))
  
  # step1：按 keyword 建目录
  key_slug <- slugify(keyword)
  base_dir <- file.path(pipeline_root, key_slug)
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  
  old_wd <- getwd()
  if (restore_wd) on.exit(setwd(old_wd), add = TRUE)
  setwd(base_dir)
  if (verbose) message("当前项目工作目录：", normalizePath(getwd(), winslash = "/"))
  
  # step2：主题
  if (!is.null(theme_file) && file.exists(theme_file)) source(theme_file)
  if (!exists("theme_pub")) theme_pub <- function(base = 11) ggplot2::theme_minimal(base_size = base)
  
  # step3：暴露
  stopifnot(file.exists(exp_rds))
  mr_exp <- readRDS(exp_rds)
  mr_exp <- as.data.table(mr_exp)
  
  exposure_dat <- mr_exp[, .(
    SNP,
    beta.exposure = beta,
    se.exposure   = se,
    effect_allele.exposure = effect_allele,
    other_allele.exposure  = other_allele,
    eaf.exposure  = eaf,
    pval.exposure = pval,
    exposure      = exposure_id,
    id.exposure   = exposure_id,
    samplesize.exposure = N
  )]
  # ---- DEDUP exposure: (id.exposure, SNP) 必须唯一 ----
  exp_dt <- as.data.table(exposure_dat)
  exp_dt <- unique(exp_dt, by = c("id.exposure", "SNP"))
  exposure_dat <- as.data.frame(exp_dt)
  # -----------------------------------------------
  snps <- unique(exposure_dat$SNP)   # <- 这行保留/放在去重后
  
  
  if (verbose) {
    message("Exposure rows: ", nrow(exposure_dat))
    message("Unique exposures: ", length(unique(exposure_dat$id.exposure)))
    message("Unique SNPs: ", length(snps))
  }
  
  # step5：token
  token <- ieugwasr::get_opengwas_jwt()
  Sys.setenv(OPENGWAS_JWT = token)
  stopifnot(nzchar(token))
  
  # step6：available_outcomes 缓存
  dir.create(dirname(ao_cache_rds), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(ao_cache_rds) && !rebuild_ao) {
    ao <- readRDS(ao_cache_rds)
    if (verbose) message("available_outcomes 从缓存读取：", ao_cache_rds, " 行=", nrow(ao))
  } else {
    ao <- TwoSampleMR::available_outcomes()
    saveRDS(ao, ao_cache_rds, compress = TRUE)
    if (verbose) message("available_outcomes 已缓存：", ao_cache_rds, " 行=", nrow(ao))
  }
  
  # step7：选择 outcomes（只改 fixed_ids/keyword 就能换）
  mode <- if (length(fixed_ids) > 0 && !use_search_fill) "fixed" else
    if (length(fixed_ids) > 0 && use_search_fill) "mixed" else "search"
  
  out_sel <- file.path(getwd(), "opengwas_outcomes_selected.tsv")
  
  top3 <- resolve_outcomes_meta(
    ao = ao,
    keyword = keyword,
    fixed_ids = fixed_ids,
    mode = mode,
    n = n_search,
    population_regex = population_regex,
    trait_regex = trait_regex,
    id_regex = id_regex,
    min_year = min_year,
    cache_tsv = out_sel,
    rebuild = rebuild_top3,
    verbose = verbose
  )
  
  # step8：ss（Steiger 用）
  top3 <- as.data.table(top3)
  if (!"sample_size_num" %in% names(top3)) top3 <- calc_sample_size_num(top3)
  top3[, ss := fifelse(
    !is.na(sample_size_num), as.numeric(sample_size_num),
    suppressWarnings(as.numeric(ncase) + as.numeric(ncontrol))
  )]
  setorder(top3, -ss)
  
  if (verbose) {
    message("\nSelected outcomes (sorted by ss):")
    print(top3[, .(id, trait, population, ss)])
  }
  # --- NEW: export selected outcomes meta to CSV ---
  out_sel_csv <- file.path(getwd(), "opengwas_outcomes_selected_main.csv")
  
  top3_export <- as.data.table(top3)
  keep_cols <- intersect(
    c("id", "trait", "population", "category", "subcategory", "ontology", "year",
      "sample_size", "ncase", "ncontrol", "sample_size_num", "ss"),
    names(top3_export)
  )
  
  data.table::fwrite(top3_export[, ..keep_cols], out_sel_csv)
  if (verbose) message("[OK] wrote selected outcomes (CSV): ", out_sel_csv)
  # ============================================================
  # step9：跑一个 outcome（闭包：自动拿到 exposure_dat/snps/fdr_cut/...）
  # 目的：你的 MR、Steiger、缓存、输出逻辑保持不变
  # ============================================================
  run_mr_one_outcome <- function(outcome_id, outcome_ss, out_prefix) {
    
    # 9.1 输出文件名（MR结果/显著子集/细胞汇总/缓存）
    out_all <- paste0(out_prefix, "_MR_all.tsv.gz")
    out_sig <- paste0(out_prefix, "_MR_FDR", fdr_tag, ".tsv.gz")
    out_ct  <- paste0(out_prefix, "_sig_counts_by_celltype.tsv")
    
    cache_mr_rds      <- paste0(out_prefix, "_MR_res.rds")
    cache_dat_rds     <- paste0(out_prefix, "_harmonised_dat2.rds")
    cache_outcome_rds <- paste0(out_prefix, "_outcome_raw.rds")
    
    # 9.2 MR结果缓存：如果已有 MR_all / MR_res，就直接读并返回（不再跑）
    if (!rebuild_mr) {
      if (file.exists(cache_mr_rds)) {
        mr_res <- as.data.table(readRDS(cache_mr_rds))
        if (verbose) message("[CACHE] MR loaded: ", cache_mr_rds, " rows=", nrow(mr_res))
        return(list(
          mr = mr_res,
          dat = if (file.exists(cache_dat_rds)) readRDS(cache_dat_rds) else NULL
        ))
      }
      if (file.exists(out_all)) {
        mr_res <- fread(out_all)
        if (verbose) message("[CACHE] MR loaded: ", out_all, " rows=", nrow(mr_res))
        return(list(
          mr = mr_res,
          dat = if (file.exists(cache_dat_rds)) readRDS(cache_dat_rds) else NULL
        ))
      }
    }
    
    # 9.3 outcome 提取（优先读 outcome_raw 缓存；否则联网 extract 并缓存）
    if (file.exists(cache_outcome_rds)) {
      outcome_dat <- readRDS(cache_outcome_rds)
    } else {
      args <- list(snps = snps, outcomes = outcome_id, proxies = FALSE)
      if ("splitsize" %in% names(formals(TwoSampleMR::extract_outcome_data))) {
        args$splitsize <- 500L
      }
      
      outcome_dat <- NULL
      for (k in 1:3) {
        outcome_dat <- tryCatch(
          do.call(TwoSampleMR::extract_outcome_data, args),
          error = function(e) e
        )
        
        if (!inherits(outcome_dat, "error")) break
        
        msg <- conditionMessage(outcome_dat)
        message("[WARN] extract_outcome_data failed (try ", k, "/3): ", msg)
        
        # ✅ 如果遇到 401（Unauthorized），自动刷新 token 再重试
        if (grepl("\\b401\\b", msg)) {
          message("[AUTH] Detected 401. Refreshing OpenGWAS JWT...")
          token <- ieugwasr::get_opengwas_jwt()
          Sys.setenv(OPENGWAS_JWT = token)
        }
        
        Sys.sleep(3 * k)
      }
      
      if (inherits(outcome_dat, "error")) stop(conditionMessage(outcome_dat))
      
      
      saveRDS(outcome_dat, cache_outcome_rds, compress = TRUE)
    }
    # ---- DEDUP outcome: (SNP, id.outcome) 必须唯一 ----
    out_dt <- as.data.table(outcome_dat)
    if (all(c("SNP", "id.outcome") %in% names(out_dt))) {
      out_dt <- unique(out_dt, by = c("SNP", "id.outcome"))
      outcome_dat <- as.data.frame(out_dt)
    }
    # --------------------------------------------------
    
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      warning("No outcome data extracted for outcome_id = ", outcome_id)
      return(NULL)
    }
    
    # 9.4 harmonise：等位基因对齐/回文位点处理
    dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2)
    
    # outcome 样本量补齐（Steiger会用）
    if (!"samplesize.outcome" %in% names(dat)) dat$samplesize.outcome <- NA
    dat$samplesize.outcome[is.na(dat$samplesize.outcome)] <- outcome_ss
    
    # ------------------------------------------------------------
    # 关键修复点：dat 是 data.frame，统计时先转 data.table 再用 by=
    # ------------------------------------------------------------
    dat_dt <- as.data.table(dat)
    # ---- DEDUP#1 harmonise: (SNP, exposure, id.outcome) 必须唯一 ----
    if (all(c("SNP","exposure","id.outcome") %in% names(dat_dt))) {
      if ("mr_keep" %in% names(dat_dt)) {
        dat_dt[, mr_keep := as.logical(mr_keep)]
        dat_dt[, mr_keep_int := as.integer(mr_keep)]
        data.table::setorder(dat_dt, exposure, id.outcome, SNP, -mr_keep_int)
        dat_dt[, mr_keep_int := NULL]
      } else {
        data.table::setorder(dat_dt, exposure, id.outcome, SNP)
      }
      dat_dt <- unique(dat_dt, by = c("SNP","exposure","id.outcome"))
      dat <- as.data.frame(dat_dt)  # 关键：让 Steiger 用到去重后的 dat
    }
    # ----------------------------------------------------------------
    dat_dt <- as.data.table(dat)   # 确保 harm_sum 用的是去重后的 dat
    
    # 9.5 统计：harmonise 后每个 exposure 实际保留的 SNP 数（mr_keep==TRUE）
    harm_sum <- NULL
    if ("mr_keep" %in% names(dat_dt)) {
      harm_sum <- dat_dt[mr_keep %in% TRUE,
                         .(nsnp_harmonised = uniqueN(SNP)),
                         by = "exposure"]
    } else {
      harm_sum <- dat_dt[, .(nsnp_harmonised = uniqueN(SNP)), by = "exposure"]
    }
    
    # 9.6 Steiger filtering：方向性检验（失败则跳过）
    steiger_sum <- NULL
    steiger_applied <- FALSE
    
    dat2_dt <- tryCatch({
      tmp <- TwoSampleMR::steiger_filtering(dat)   # 通常返回 data.frame
      tmp_dt <- as.data.table(tmp)
      
      if ("steiger_dir" %in% names(tmp_dt)) {
        steiger_applied <- TRUE
        
        # 统计：Steiger 通过的 SNP 数（按 exposure）
        if ("mr_keep" %in% names(tmp_dt)) {
          steiger_sum <- tmp_dt[mr_keep %in% TRUE,
                                .(nsnp_steiger_keep = uniqueN(SNP[steiger_dir %in% TRUE])),
                                by = "exposure"]
        } else {
          steiger_sum <- tmp_dt[, .(nsnp_steiger_keep = uniqueN(SNP[steiger_dir %in% TRUE])),
                                by = "exposure"]
        }
        
        # 用于 MR 的数据：只保留 steiger_dir==TRUE
        tmp_dt <- tmp_dt[steiger_dir %in% TRUE]
      }
      
      tmp_dt
    }, error = function(e) {
      message("[WARN] Steiger filtering failed, continue without it: ", e$message)
      dat_dt
    })
    # ---- DEDUP#2 steiger: (SNP, exposure, id.outcome) 必须唯一 ----
    if (all(c("SNP","exposure","id.outcome") %in% names(dat2_dt))) {
      if ("mr_keep" %in% names(dat2_dt)) {
        dat2_dt[, mr_keep := as.logical(mr_keep)]
        dat2_dt[, mr_keep_int := as.integer(mr_keep)]
        data.table::setorder(dat2_dt, exposure, id.outcome, SNP, -mr_keep_int)
        dat2_dt[, mr_keep_int := NULL]
      } else {
        data.table::setorder(dat2_dt, exposure, id.outcome, SNP)
      }
      dat2_dt <- unique(dat2_dt, by = c("SNP","exposure","id.outcome"))
    }
    # ---------------------------------------------------------------
    
    if (nrow(dat2_dt) == 0) {
      warning("All instruments removed after harmonise/Steiger for outcome_id = ", outcome_id)
      return(NULL)
    }
    
    # 9.7 MR：单工具 Wald；多工具 IVW
    # 为兼容性，喂给 TwoSampleMR::mr 时用 data.frame
    dat2 <- as.data.frame(dat2_dt)
    
    mr_res <- TwoSampleMR::mr(dat2, method_list = c("mr_wald_ratio", "mr_ivw"))
    mr_res <- as.data.table(mr_res)
    
    # 每个 exposure 只保留“该用的那一个方法”
    mr_res <- mr_res[
      (nsnp == 1 & method == "Wald ratio") |
        (nsnp > 1 & method == "Inverse variance weighted")
    ]
    if (nrow(mr_res) == 0) {
      warning("No MR results left after method selection for outcome_id = ", outcome_id)
      return(NULL)
    }
    
    # 9.8 BH-FDR + 拆 CELL_TYPE/GENE_ID
    mr_res <- unique(mr_res, by = c("exposure", "id.outcome", "method"))
    mr_res[, fdr := p.adjust(pval, method = "BH")]
    
    mr_res[, c("CELL_TYPE", "GENE_ID") := tstrsplit(exposure, "\\|")]
    
    # 9.9 ★把 harmonise/Steiger 的统计 merge 进 MR 结果
    mr_res <- merge(mr_res, harm_sum, by = "exposure", all.x = TRUE)
    if (!is.null(steiger_sum) && nrow(steiger_sum) > 0) {
      mr_res <- merge(mr_res, steiger_sum, by = "exposure", all.x = TRUE)
    }
    
    # steiger 保留比例（如果 steiger 没跑或没有列，就会是 NA）
    if (!"nsnp_steiger_keep" %in% names(mr_res)) mr_res[, nsnp_steiger_keep := NA_real_]
    if (!"nsnp_harmonised" %in% names(mr_res)) mr_res[, nsnp_harmonised := NA_real_]
    mr_res[, steiger_keep_rate := nsnp_steiger_keep / nsnp_harmonised]
    mr_res[, steiger_applied := steiger_applied]
    
    # 9.10 输出 + 缓存（下次命中缓存就不再跑）
    fwrite(mr_res, out_all, sep = "\t", compress = "gzip")
    fwrite(mr_res[mr_res$fdr <= fdr_cut], out_sig, sep = "\t", compress = "gzip")
    
    sig_ct <- mr_res[fdr <= fdr_cut,
                     .(n_sig = uniqueN(exposure)),
                     by = "CELL_TYPE"][order(-n_sig)]
    fwrite(sig_ct, out_ct, sep = "\t")
    
    saveRDS(mr_res, cache_mr_rds, compress = TRUE)
    saveRDS(dat2_dt, cache_dat_rds, compress = TRUE)  # 缓存 SNP-level（data.table 也行）
    
    if (verbose) {
      message("\n[OK] outcome: ", outcome_id)
      message("  harmonised rows: ", nrow(dat2_dt))
      message("  MR tests: ", nrow(mr_res))
      message("  FDR<= ", fdr_cut, ": ", sum(mr_res$fdr <= fdr_cut))
      message("  wrote:\n   ", out_all, "\n   ", out_sig, "\n   ", out_ct)
    }
    
    list(mr = mr_res, dat = dat2_dt)
  }
  
  # step10：循环跑 outcomes
  # step10：循环跑 outcomes（每个 outcome 一个子目录）
  results <- list()
  for (i in seq_len(nrow(top3))) {
    outcome_id <- top3$id[i]
    outcome_ss <- top3$ss[i]
    
    outcome_dir <- file.path(getwd(), paste0("outcome_", outcome_id))
    dir.create(outcome_dir, showWarnings = FALSE, recursive = TRUE)
    
    out_prefix <- file.path(outcome_dir, paste0("MR_", outcome_id))
    
    results[[outcome_id]] <- run_mr_one_outcome(outcome_id, outcome_ss, out_prefix)
    Sys.sleep(sleep)
  }
  
  # step11：合并 MR_all
  # step11：合并 MR_all（递归搜子目录）
  mr_files <- list.files(getwd(), pattern = "^MR_.*_MR_all\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)
  mr_files <- unique(normalizePath(mr_files, winslash = "/", mustWork = TRUE))
  stopifnot(length(mr_files) > 0)
  
  mr_all <- rbindlist(lapply(mr_files, fread), fill = TRUE)
  mr_all <- unique(mr_all, by = c("exposure", "id.outcome", "method", "b", "se", "pval"))
  
  if (verbose) message("[INFO] merged MR_all files: ", length(mr_files), " total rows=", nrow(mr_all))
  
  # step12：主表 + TopHits
  num_cols <- intersect(c("b", "se", "pval", "fdr", "nsnp", "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate"), names(mr_all))
  mr_all[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]
  
  if (!all(c("CELL_TYPE", "GENE_ID") %in% names(mr_all))) {
    src <- if ("exposure" %in% names(mr_all)) "exposure" else "id.exposure"
    mr_all[, c("CELL_TYPE", "GENE_ID") := tstrsplit(get(src), "\\|")]
  }
  
  mr_all[, `:=`(
    OR       = exp(b),
    OR_lci95 = exp(b - 1.96 * se),
    OR_uci95 = exp(b + 1.96 * se),
    OR_CI    = sprintf("%.3f (%.3f–%.3f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
  )]
  
  key_cols <- intersect(c(
    "CELL_TYPE", "GENE_ID",
    "id.outcome", "outcome",
    "method", "nsnp",
    "OR_CI", "OR", "OR_lci95", "OR_uci95",
    "b", "se", "pval", "fdr",
    "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate", "steiger_applied"
  ), names(mr_all))
  
  key_table <- mr_all[is.finite(fdr) & fdr <= fdr_cut, ..key_cols]
  out_key <- file.path(getwd(), paste0("Table_MR_key_FDR", fdr_tag, ".tsv"))
  fwrite(key_table, out_key, sep = "\t")
  if (verbose) message("[OK] wrote (MAIN TABLE): ", out_key, " rows=", nrow(key_table), " (FDR<=", fdr_cut, ")")
  
  top_table <- mr_all[is.finite(pval)][order(id.outcome, pval)][, head(.SD, top_n_per_outcome), by = id.outcome][, ..key_cols]
  out_top <- file.path(getwd(), "Table_MR_key_TopHits.tsv")
  fwrite(top_table, out_top, sep = "\t")
  if (verbose) message("[OK] wrote (TOP HITS): ", out_top, " rows=", nrow(top_table))
  
  # step12.7：加 gene_symbol
  gene_map <- build_gene_map(vpub_file, vpub_sheet)
  
  if (nrow(key_table) > 0) {
    key_table[, GENE_ID := sub("\\.\\d+$", "", GENE_ID)]
    key_table <- merge(key_table, gene_map, by = "GENE_ID", all.x = TRUE)
    front <- intersect(c("CELL_TYPE", "gene_symbol", "GENE_ID"), names(key_table))
    setcolorder(key_table, c(front, setdiff(names(key_table), front)))
    
    out_key2 <- file.path(getwd(), paste0("Table_MR_key_FDR", fdr_tag, "_with_symbol.tsv"))
    fwrite(key_table, out_key2, sep = "\t")
    if (verbose) message("[OK] wrote: ", out_key2, " rows=", nrow(key_table))
  } else {
    out_key2 <- NA_character_
  }
  
  top_table[, GENE_ID := sub("\\.\\d+$", "", GENE_ID)]
  top_table <- merge(top_table, gene_map, by = "GENE_ID", all.x = TRUE)
  front <- intersect(c("CELL_TYPE", "gene_symbol", "GENE_ID"), names(top_table))
  setcolorder(top_table, c(front, setdiff(names(top_table), front)))
  
  out_top2 <- file.path(getwd(), "Table_MR_key_TopHits_with_symbol.tsv")
  fwrite(top_table, out_top2, sep = "\t")
  if (verbose) message("[OK] wrote: ", out_top2, " rows=", nrow(top_table))
  # ============================================================
  # ============================================================
  # step12.7.9：Per-outcome forest plots（每个 outcome 一个森林图）
  # 位置：gene_symbol merge 完成后、决策表/总森林图之前
  # 输出：
  #   - <base_dir>/outcome_<id>/Fig_MR_forest.png
  #   - <base_dir>/outcome_<id>/Fig_MR_forest_data.tsv
  #   - <base_dir>/outcome_<id>/steiger_only/Fig_MR_forest.png （若有）
  # 并在 RStudio Plots 里依次显示
  # ============================================================
  
  forest_by_outcome <- list()
  forest_steiger_by_outcome <- list()
  
  if (isTRUE(do_forest_plot)) {
    
    # 用你选中的 outcomes 顺序画图（更可控）
    outcome_ids <- as.character(top3$id)
    
    for (oid in outcome_ids) {
      
      out_dir_i <- file.path(getwd(), paste0("outcome_", oid))
      dir.create(out_dir_i, showWarnings = FALSE, recursive = TRUE)
      
      # 子集：该 outcome 的 key_table / top_table
      kt_i <- data.table::as.data.table(key_table)[id.outcome %in% oid]
      tt_i <- data.table::as.data.table(top_table)[id.outcome %in% oid]
      
      # 如果该 outcome 没有任何行，跳过
      if (nrow(kt_i) == 0 && nrow(tt_i) == 0) next
      
      # 画该 outcome 的主森林图（优先 kt_i，否则 tt_i）
      forest_i <- make_forest_plot(
        key_table = if (nrow(kt_i) > 0) kt_i else NULL,
        top_table = tt_i,
        out_dir = out_dir_i,
        theme_pub = theme_pub
      )
      forest_by_outcome[[oid]] <- forest_i
      
      if (verbose) message("[OK] saved (per-outcome forest): ", forest_i$out_png)
      if (interactive() && !is.null(forest_i$plot)) print(forest_i$plot)
      
      # ---- Steiger-only（该 outcome）----
      # 只在有 steiger_applied 列时做筛选
      kt_s <- NULL
      if ("steiger_applied" %in% names(kt_i)) {
        kt_i[, steiger_applied := as.logical(steiger_applied)]
        kt_s <- kt_i[steiger_applied %in% TRUE]
      }
      
      tt_s <- tt_i
      if ("steiger_applied" %in% names(tt_s)) {
        tt_s[, steiger_applied := as.logical(steiger_applied)]
        tt_s <- tt_s[steiger_applied %in% TRUE]
      } else {
        tt_s <- tt_s[0]
      }
      
      if (!((is.null(kt_s) || nrow(kt_s) == 0) && nrow(tt_s) == 0)) {
        out_dir_s <- file.path(out_dir_i, "steiger_only")
        dir.create(out_dir_s, showWarnings = FALSE, recursive = TRUE)
        
        forest_s <- make_forest_plot(
          key_table = if (!is.null(kt_s) && nrow(kt_s) > 0) kt_s else NULL,
          top_table = tt_s,
          out_dir = out_dir_s,
          theme_pub = theme_pub
        )
        forest_steiger_by_outcome[[oid]] <- forest_s
        
        if (verbose) message("[OK] saved (per-outcome Steiger-only): ", forest_s$out_png)
        if (interactive() && !is.null(forest_s$plot)) print(forest_s$plot)
      } else {
        if (verbose) message("[INFO] No Steiger-applied rows for per-outcome Steiger-only: ", oid)
      }
    }
  }
  
  # ============================================================
  # step12.8：少而稳 - Go/No-go 决策表（基于 key_table with symbol）
  # 位置：gene_symbol merge 之后、森林图之前
  # ============================================================
  decision_dir <- getwd()
  
  # 1) 只从“显著表”出发（少而稳核心）
  dec_src <- NULL
  if (exists("key_table") && is.data.frame(key_table) && nrow(key_table) > 0) {
    dec_src <- data.table::as.data.table(key_table)
  } else {
    # 少而稳：显著表为空，则 A_go 为空；仍可输出空表做记录
    dec_src <- data.table::data.table()
  }
  
  # 2) 设定保守阈值（你可后面再调）
  min_nsnp <- 3L
  min_steiger_rate <- 0.70
  max_ci_ratio <- 3.0
  min_effect_log_or <- log(1.1)  # 软门槛，用于排序，不一定作为硬过滤
  
  # 3) 生成审计字段与失败原因
  decision_audit <- copy(dec_src)
  if ("steiger_applied" %in% names(decision_audit)) decision_audit[, steiger_applied := as.logical(steiger_applied)]
  if (nrow(decision_audit) > 0) {
    # 确保数值列是数值
    num_cols <- intersect(
      c("fdr", "pval", "OR", "OR_lci95", "OR_uci95", "nsnp",
        "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate"),
      names(decision_audit)
    )
    decision_audit[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]
    
    # 逐条门槛（硬）
    decision_audit[, pass_fdr := is.finite(fdr) & fdr <= fdr_cut]
    decision_audit[, pass_nsnp := is.finite(nsnp) & nsnp >= min_nsnp &
                     is.finite(nsnp_harmonised) & nsnp_harmonised >= min_nsnp]
    
    decision_audit[, pass_steiger := (steiger_applied %in% TRUE) &
                     is.finite(steiger_keep_rate) & steiger_keep_rate >= min_steiger_rate]
    
    decision_audit[, ci_ratio := OR_uci95 / OR_lci95]
    decision_audit[, pass_ci := is.finite(OR_lci95) & is.finite(OR_uci95) &
                     OR_lci95 > 0 & OR_uci95 > 0 & is.finite(ci_ratio) & ci_ratio <= max_ci_ratio]
    
    # 软门槛：效应大小（主要用于排序）
    decision_audit[, effect_log_or := abs(log(OR))]
    decision_audit[, pass_effect := is.finite(effect_log_or) & effect_log_or >= min_effect_log_or]
    
    # 汇总失败原因（便于“少而稳”解释）
    decision_audit[, fail_reason := ""]
    decision_audit[(pass_fdr == FALSE), fail_reason := paste0(fail_reason, "FDR;")]
    decision_audit[(pass_fdr == TRUE)  & (pass_nsnp == FALSE), fail_reason := paste0(fail_reason, "Low_nsnp;")]
    decision_audit[(pass_fdr == TRUE)  & (pass_nsnp == TRUE) & (pass_steiger == FALSE), fail_reason := paste0(fail_reason, "Steiger;")]
    decision_audit[(pass_fdr == TRUE)  & (pass_nsnp == TRUE) & (pass_steiger == TRUE) & (pass_ci == FALSE), fail_reason := paste0(fail_reason, "Wide_CI;")]
    
    
    # A_go：全部硬门槛都过
    decision_audit[, tier := ifelse(pass_fdr & pass_nsnp & pass_steiger & pass_ci, "A_go", "B_hold")]
    
    # 给 A_go 做一个排序：先 fdr，再 steiger_keep_rate，再 nsnp_harmonised，再效应
    setorder(decision_audit, tier, fdr, -steiger_keep_rate, -nsnp_harmonised, -effect_log_or, pval)
  }
  
  # 4) 输出（A_go + 审计）
  out_audit <- file.path(decision_dir, "Table_MR_decision_audit.tsv")
  data.table::fwrite(decision_audit, out_audit, sep = "\t")
  if (verbose) message("[OK] wrote: ", out_audit, " rows=", nrow(decision_audit))
  
  decision_A <- decision_audit[0]
  if (nrow(decision_audit) > 0 && "tier" %in% names(decision_audit)) {
    decision_A <- decision_audit[tier == "A_go"]
  }
  
  out_A <- file.path(decision_dir, "Table_MR_decision_A_go.tsv")
  data.table::fwrite(decision_A, out_A, sep = "\t")
  if (verbose) message("[OK] wrote: ", out_A, " rows=", nrow(decision_A))
  
  # step13：森林图（可选）
  forest <- NULL
  if (isTRUE(do_forest_plot)) {
    forest <- make_forest_plot(
      key_table = if (nrow(key_table) > 0) key_table else NULL,
      top_table = top_table,
      out_dir = getwd(),
      theme_pub = theme_pub
    )
    if (verbose) message("[OK] saved: ", forest$out_png)
    if (interactive() && !is.null(forest)) print(forest$plot)   # ✅ NEW: show in RStudio Plots
  }
  
  # TopHits QC（可选）
  qc <- NULL
  if (isTRUE(do_qc_top)) {
    qc <- qc_top_table(top_table, out_dir = getwd())
  }
  # --- NEW: Forest plot (Steiger-only) ---
  forest_steiger <- NULL
  if (isTRUE(do_forest_plot)) {
    
    # 只保留 steiger_applied == TRUE 的结果
    key_s <- NULL
    if (exists("key_table") && is.data.frame(key_table) && nrow(key_table) > 0 && "steiger_applied" %in% names(key_table)) {
      key_s <- data.table::as.data.table(key_table)
      key_s[, steiger_applied := as.logical(steiger_applied)]
      key_s <- key_s[steiger_applied %in% TRUE]
    }
    
    top_s <- data.table::as.data.table(top_table)
    if ("steiger_applied" %in% names(top_s)) {
      top_s[, steiger_applied := as.logical(steiger_applied)]
      top_s <- top_s[steiger_applied %in% TRUE]
    } else {
      top_s <- top_s[0]
    }
    
    # 若过滤后完全没数据，就不画
    if ((is.null(key_s) || nrow(key_s) == 0) && nrow(top_s) == 0) {
      if (verbose) message("[INFO] No Steiger-applied rows available for Steiger-only forest plot.")
    } else {
      out_dir_steiger <- file.path(getwd(), "steiger_only")
      dir.create(out_dir_steiger, showWarnings = FALSE, recursive = TRUE)
      
      forest_steiger <- make_forest_plot(
        key_table = if (!is.null(key_s) && nrow(key_s) > 0) key_s else NULL,
        top_table = top_s,
        out_dir = out_dir_steiger,
        theme_pub = theme_pub
      )
      if (verbose) message("[OK] saved (Steiger-only): ", forest_steiger$out_png)
      if (interactive() && !is.null(forest_steiger)) print(forest_steiger$plot)  # ✅ NEW: show in Plots
    }
  }
  
  # ============================================================
  # stepX：自动 View 最重要的东西（表格优先）
  # 目的：跑完立刻在 RStudio 看到结果表；Windows 可自动打开目录/图片
  # ============================================================
  if (interactive()) {
    
    # 1) 你最终选了哪些 outcome
    top3_view <- top3[, .(id, trait, population, ss)]
    utils::View(as.data.frame(top3_view), title = "Selected outcomes")
    
    # 2) 主表（显著结果）
    if (exists("key_table") && is.data.frame(key_table) && nrow(key_table) > 0) {
      key_view_cols <- intersect(c(
        "CELL_TYPE", "gene_symbol", "GENE_ID",
        "outcome", "id.outcome",
        "OR_CI", "OR", "OR_lci95", "OR_uci95",
        "b", "se", "pval", "fdr", "nsnp",
        "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate", "steiger_applied"
      ), names(key_table))
      
      key_view <- data.table::as.data.table(key_table)[order(fdr, pval)][, ..key_view_cols]
      utils::View(as.data.frame(key_view), title = paste0("MR KEY (FDR<=", fdr_cut, ")"))
    } else {
      message("[VIEW] key_table 为空（没有显著结果），将打开 TopHits。")
    }
    
    # 3) TopHits（每个 outcome 前 top_n_per_outcome）
    if (exists("top_table") && is.data.frame(top_table) && nrow(top_table) > 0) {
      top_view_cols <- intersect(c(
        "CELL_TYPE", "gene_symbol", "GENE_ID",
        "outcome", "id.outcome",
        "OR_CI", "OR", "OR_lci95", "OR_uci95",
        "b", "se", "pval", "fdr", "nsnp",
        "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate", "steiger_applied"
      ), names(top_table))
      
      top_view <- data.table::as.data.table(top_table)[order(id.outcome, pval)][, ..top_view_cols]
      utils::View(as.data.frame(top_view),
                  title = paste0("MR TopHits (top ", top_n_per_outcome, "/outcome)"))
    }
    
    # 4) 按 CELL_TYPE 的汇总（表格）
    if (exists("mr_all") && is.data.frame(mr_all)) {
      mr_dt <- data.table::as.data.table(mr_all)
      
      # 兜底：如果还没拆 CELL_TYPE
      if (!"CELL_TYPE" %in% names(mr_dt) && "exposure" %in% names(mr_dt)) {
        mr_dt[, CELL_TYPE := data.table::tstrsplit(exposure, "\\|")[[1]]]
      }
      
      ct_sum <- mr_dt[is.finite(fdr), .(
        n_tests = .N,
        n_sig_fdr = sum(fdr <= fdr_cut, na.rm = TRUE),
        best_fdr = min(fdr, na.rm = TRUE),
        best_p = min(pval, na.rm = TRUE)
      ), by = CELL_TYPE][order(-n_sig_fdr, best_fdr)]
      
      utils::View(as.data.frame(ct_sum), title = paste0("Summary by CELL_TYPE (FDR<=", fdr_cut, ")"))
    }
    
    # 5) Windows：打开输出目录和森林图
    if (.Platform$OS.type == "windows") {
      try(shell.exec(getwd()), silent = TRUE)
      if (exists("forest") && !is.null(forest) && !is.null(forest$out_png) && file.exists(forest$out_png)) {
        try(shell.exec(forest$out_png), silent = TRUE)
      }
    }
  }
  
  # 返回：方便你在外部脚本里拿对象继续画图/汇报
  invisible(list(
    base_dir = base_dir,
    keyword = keyword,
    outcomes_meta = top3,
    outcome_ids = top3$id,
    mr_all = mr_all,
    key_table = key_table,
    top_table = top_table,
    out_key = out_key,
    out_key_with_symbol = out_key2,
    out_top = out_top,
    out_top_with_symbol = out_top2,
    forest_by_outcome = forest_by_outcome,
    forest_steiger_by_outcome = forest_steiger_by_outcome,
    forest = forest,
    forest_steiger = forest_steiger,   # ✅ NEW: return steiger-only plot object
    qc = qc,
    results = results
  ))
}
## =========================================================
## 5. 主执行区
## =========================================================
RESULT_DIR <- file.path(PIPELINE_ROOT, slugify(PROJECT_NAME))

if (dir.exists(RESULT_DIR)) {
  stop(
    "结果目录已存在：\n", RESULT_DIR,
    "\n为避免旧结果混入，请执行以下任一操作：",
    "\n1) 改一个新的 PROJECT_NAME",
    "\n2) 手动删除这个旧目录后再跑"
  )
}

res <- run_project(
  keyword = PROJECT_NAME,
  fixed_ids = FIXED_IDS,
  use_search_fill = FALSE,
  n_search = N_SEARCH,
  trait_regex = TRAIT_REGEX,
  id_regex = ID_REGEX,
  population_regex = POPULATION_REGEX,
  pipeline_root = PIPELINE_ROOT,
  rebuild_top3 = REBUILD_TOP3,
  rebuild_mr = REBUILD_MR,
  rebuild_ao = REBUILD_AO,
  do_forest_plot = DO_FOREST_PLOT,
  do_qc_top = DO_QC_TOP,
  verbose = VERBOSE_RUN
)

## 保存完整 res 对象
res_rds <- file.path(res$base_dir, "project_res.rds")
saveRDS(res, res_rds, compress = "xz")

cat("[OK] res 已保存：\n", res_rds, "\n\n")

cat("\n[DONE]\n")
cat("结果目录：\n", res$base_dir, "\n\n")

print(res$outcomes_meta[, .(id, trait, population, ss)])
print(res$outcome_ids)
