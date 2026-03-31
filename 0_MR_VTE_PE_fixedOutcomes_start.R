## =========================================================  # 本行目的：说明步骤或上下文
## 0. 当前项目专用脚本：VTE + PE  # 本行目的：说明步骤或上下文
## 说明：  # 本行目的：说明步骤或上下文
##   - 这是一个单文件脚本  # 本行目的：说明步骤或上下文
##   - 不再 source 其他 core.R  # 本行目的：说明步骤或上下文
##   - 只跑两个 outcome：  # 本行目的：说明步骤或上下文
##       1) finn-b-I9_VTE  # 本行目的：说明步骤或上下文
##       2) finn-b-I9_PULMEMB  # 本行目的：说明步骤或上下文
##   - 先把“项目参数区”放最上面  # 本行目的：说明步骤或上下文
##   - 接下来你把原来 core 的全部函数内容贴到下面指定位置  # 本行目的：说明步骤或上下文
## =========================================================  # 本行目的：说明步骤或上下文


## =========================================================  # 本行目的：说明步骤或上下文
## 1. 基础设置  # 本行目的：说明步骤或上下文
## =========================================================  # 本行目的：说明步骤或上下文
rm(list = ls())  # 本行目的：设置参数或赋值对象
gc()  # 本行目的：执行当前流程语句
options(stringsAsFactors = FALSE)  # 本行目的：设置参数或赋值对象


## =========================================================  # 本行目的：说明步骤或上下文
## 2. 当前项目输入区（以后主要改这里）  # 本行目的：说明步骤或上下文
## =========================================================  # 本行目的：说明步骤或上下文
PROJECT_NAME <- "Venous_thromboembolism_PE_VTE_only"  # 本行目的：设置参数或赋值对象

FIXED_IDS <- c(  # 本行目的：设置参数或赋值对象
  "finn-b-I9_VTE",  # 本行目的：执行当前流程语句
  "finn-b-I9_PULMEMB"  # 本行目的：执行当前流程语句
)  # 本行目的：执行当前流程语句

PIPELINE_ROOT <- "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE"  # 本行目的：设置参数或赋值对象

TRAIT_REGEX <- NULL  # 本行目的：设置参数或赋值对象
ID_REGEX <- NULL  # 本行目的：设置参数或赋值对象
POPULATION_REGEX <- "European|EUR|UKB|UK Biobank|Finnish|FinnGen|Finland"  # 本行目的：设置参数或赋值对象

REBUILD_TOP3 <- TRUE  # 本行目的：设置参数或赋值对象
REBUILD_MR   <- TRUE  # 本行目的：设置参数或赋值对象
REBUILD_AO   <- FALSE  # 本行目的：设置参数或赋值对象

DO_FOREST_PLOT <- TRUE  # 本行目的：设置参数或赋值对象
DO_QC_TOP      <- TRUE  # 本行目的：设置参数或赋值对象
VERBOSE_RUN    <- TRUE  # 本行目的：设置参数或赋值对象

N_SEARCH <- length(FIXED_IDS)  # 本行目的：设置参数或赋值对象


## =========================================================  # 本行目的：说明步骤或上下文
## 3. 下面请粘贴你原来 core 的全部函数内容  # 本行目的：说明步骤或上下文
## 粘贴范围：  # 本行目的：说明步骤或上下文
##   从 slugify <- function(x) {  # 本行目的：说明步骤或上下文
##   一直到 run_project <- function(...) { ... } 结束  # 本行目的：说明步骤或上下文
## =========================================================  # 本行目的：说明步骤或上下文



## ======= 你把原来 core 的内容粘贴在这里下面 =======  # 本行目的：说明步骤或上下文



## =========================================================  # 本行目的：说明步骤或上下文
## 4. 先不要在这里写执行代码  # 本行目的：说明步骤或上下文
## 等你把 core 贴完后，再发给我  # 本行目的：说明步骤或上下文
## 我再给你最后一段“主执行区”  # 本行目的：说明步骤或上下文
## =========================================================  # 本行目的：说明步骤或上下文
slugify <- function(x) {  # 本行目的：定义函数入口
  x <- trimws(x)  # 本行目的：设置参数或赋值对象
  x <- gsub("[^A-Za-z0-9]+", "_", x)  # 本行目的：设置参数或赋值对象
  x <- gsub("^_+|_+$", "", x)  # 本行目的：设置参数或赋值对象
  x  # 本行目的：执行当前流程语句
}  # 本行目的：组织代码块结构

calc_sample_size_num <- function(dt) {  # 本行目的：定义函数入口
  dt <- as.data.table(copy(dt))  # 本行目的：设置参数或赋值对象
  dt[, sample_size_num := suppressWarnings(as.numeric(sample_size))]  # 本行目的：设置参数或赋值对象
  if (all(is.na(dt$sample_size_num)) && all(c("ncase", "ncontrol") %in% names(dt))) {  # 本行目的：组织代码块结构
    dt[, sample_size_num := suppressWarnings(as.numeric(ncase) + as.numeric(ncontrol))]  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  dt  # 本行目的：执行当前流程语句
}  # 本行目的：组织代码块结构

resolve_outcomes_meta <- function(  # 本行目的：定义函数入口
    ao,  # 本行目的：执行当前流程语句
    keyword = NULL,  # 本行目的：设置参数或赋值对象
    fixed_ids = character(),  # 本行目的：设置参数或赋值对象
    mode = c("mixed", "fixed", "search"),  # 本行目的：设置参数或赋值对象
    n = 3,  # 本行目的：设置参数或赋值对象
    population_regex = "European|EUR|UKB|UK Biobank|Finnish|FinnGen|Finland",  # 本行目的：设置参数或赋值对象
    trait_regex = NULL,  # 本行目的：设置参数或赋值对象
    id_regex = NULL,  # 本行目的：设置参数或赋值对象
    min_year = NA_integer_,  # 本行目的：设置参数或赋值对象
    cache_tsv = NULL,  # 本行目的：设置参数或赋值对象
    rebuild = FALSE,  # 本行目的：设置参数或赋值对象
    verbose = TRUE  # 本行目的：设置参数或赋值对象
) {  # 本行目的：组织代码块结构
  mode <- match.arg(mode)  # 本行目的：设置参数或赋值对象
  
  # 目的：有缓存且不重建 -> 直接读缓存  # 本行目的：说明步骤或上下文
  if (!is.null(cache_tsv) && file.exists(cache_tsv) && !rebuild) {  # 本行目的：组织代码块结构
    if (verbose) message("[CACHE] loaded selected outcomes: ", cache_tsv)  # 本行目的：进行条件判断
    out <- fread(cache_tsv)  # 本行目的：读取输入或缓存数据
    return(unique(out, by = "id"))  # 本行目的：返回结果对象
  }  # 本行目的：组织代码块结构
  
  ao_dt <- as.data.table(ao)  # 本行目的：设置参数或赋值对象
  ao_dt <- ao_dt[!is.na(id)]  # 本行目的：设置参数或赋值对象
  ao_dt <- unique(ao_dt, by = "id")  # 本行目的：设置参数或赋值对象
  
  fixed_ids <- unique(fixed_ids[nzchar(fixed_ids)])  # 本行目的：设置参数或赋值对象
  fixed_hit <- ao_dt[id %in% fixed_ids]  # 本行目的：设置参数或赋值对象
  
  # 目的：固定模式仅使用 fixed_ids  # 本行目的：说明步骤或上下文
  if (mode == "fixed") {  # 本行目的：组织代码块结构
    if (length(fixed_ids) > 0 && nrow(fixed_hit) == 0) {  # 本行目的：组织代码块结构
      stop("mode='fixed' 但 fixed_ids 在 available_outcomes 中都未命中，请检查 id。")  # 本行目的：在异常条件下终止并提示
    }  # 本行目的：组织代码块结构
    fixed_hit <- calc_sample_size_num(fixed_hit)  # 本行目的：设置参数或赋值对象
    if (!is.null(cache_tsv)) fwrite(fixed_hit, cache_tsv, sep = "\t")  # 本行目的：进行条件判断
    return(fixed_hit)  # 本行目的：返回结果对象
  }  # 本行目的：组织代码块结构
  
  # 目的：search/mixed 模式需要 keyword；若缺失则 mixed 退化为 fixed  # 本行目的：说明步骤或上下文
  if (is.null(keyword) || !nzchar(keyword)) {  # 本行目的：组织代码块结构
    if (mode == "search") stop("mode='search' 必须提供 keyword。")  # 本行目的：进行条件判断
    fixed_hit <- calc_sample_size_num(fixed_hit)  # 本行目的：设置参数或赋值对象
    if (!is.null(cache_tsv)) fwrite(fixed_hit, cache_tsv, sep = "\t")  # 本行目的：进行条件判断
    return(fixed_hit)  # 本行目的：返回结果对象
  }  # 本行目的：组织代码块结构
  
  # 目的：按 trait keyword 找候选  # 本行目的：说明步骤或上下文
  meta <- ao_dt[!is.na(trait) & grepl(keyword, trait, ignore.case = TRUE)]  # 本行目的：设置参数或赋值对象
  
  # 目的：人群过滤（列存在才用）  # 本行目的：说明步骤或上下文
  if ("population" %in% names(meta) && !is.null(population_regex) && nzchar(population_regex)) {  # 本行目的：组织代码块结构
    meta <- meta[grepl(population_regex, population, ignore.case = TRUE)]  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  # 目的：trait 二次过滤（更精准）  # 本行目的：说明步骤或上下文
  if (!is.null(trait_regex) && nzchar(trait_regex)) {  # 本行目的：组织代码块结构
    meta <- meta[grepl(trait_regex, trait, ignore.case = TRUE)]  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  # 目的：按 id 限定来源（如 ^ukb- / ^finn-）  # 本行目的：说明步骤或上下文
  if (!is.null(id_regex) && nzchar(id_regex)) {  # 本行目的：组织代码块结构
    meta <- meta[grepl(id_regex, id, ignore.case = TRUE)]  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  # 目的：限定年份（列存在才用）  # 本行目的：说明步骤或上下文
  if (!is.na(min_year) && "year" %in% names(meta)) {  # 本行目的：组织代码块结构
    yr <- suppressWarnings(as.integer(meta$year))  # 本行目的：设置参数或赋值对象
    meta <- meta[!is.na(yr) & yr >= min_year]  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  if (nrow(meta) == 0) {  # 本行目的：组织代码块结构
    if (mode == "search") {  # 本行目的：组织代码块结构
      stop(  # 本行目的：在异常条件下终止并提示
        "[STOP] No OpenGWAS outcomes matched keyword = ", keyword,  # 本行目的：设置参数或赋值对象
        "\n建议：换同义词、放宽 population/trait 过滤、或先不限定 EUR。"  # 本行目的：执行当前流程语句
      )  # 本行目的：执行当前流程语句
    }  # 本行目的：组织代码块结构
    fixed_hit <- calc_sample_size_num(fixed_hit)  # 本行目的：设置参数或赋值对象
    if (!is.null(cache_tsv)) fwrite(fixed_hit, cache_tsv, sep = "\t")  # 本行目的：进行条件判断
    return(fixed_hit)  # 本行目的：返回结果对象
  }  # 本行目的：组织代码块结构
  
  # 目的：按样本量排序优先大样本 outcome  # 本行目的：说明步骤或上下文
  meta <- calc_sample_size_num(meta)  # 本行目的：设置参数或赋值对象
  setorder(meta, -sample_size_num, trait)  # 本行目的：执行当前流程语句
  meta <- unique(meta, by = "id")  # 本行目的：设置参数或赋值对象
  
  if (mode == "search") {  # 本行目的：组织代码块结构
    picked <- meta[1:min(n, .N)]  # 本行目的：设置参数或赋值对象
    if (!is.null(cache_tsv)) fwrite(picked, cache_tsv, sep = "\t")  # 本行目的：进行条件判断
    return(picked)  # 本行目的：返回结果对象
  }  # 本行目的：组织代码块结构
  
  # 目的：mixed：固定保留 + 搜索补齐到 n  # 本行目的：说明步骤或上下文
  meta2 <- meta[!id %in% fixed_ids]  # 本行目的：设置参数或赋值对象
  need <- max(0L, n - nrow(fixed_hit))  # 本行目的：设置参数或赋值对象
  picked_search <- if (need > 0L) meta2[1:min(need, .N)] else meta2[0]  # 本行目的：进行条件判断
  
  out <- rbindlist(list(fixed_hit, picked_search), fill = TRUE)  # 本行目的：设置参数或赋值对象
  out <- unique(out, by = "id")  # 本行目的：设置参数或赋值对象
  out <- calc_sample_size_num(out)  # 本行目的：设置参数或赋值对象
  setorder(out, -sample_size_num, trait)  # 本行目的：执行当前流程语句
  
  if (!is.null(cache_tsv)) fwrite(out, cache_tsv, sep = "\t")  # 本行目的：进行条件判断
  out  # 本行目的：执行当前流程语句
}  # 本行目的：组织代码块结构

build_gene_map <- function(vpub_file, vpub_sheet = "S3") {  # 本行目的：定义函数入口
  stopifnot(file.exists(vpub_file))  # 本行目的：执行当前流程语句
  s3 <- as.data.table(readxl::read_xlsx(vpub_file, sheet = vpub_sheet, skip = 1))  # 本行目的：设置参数或赋值对象
  if (nrow(s3) >= 1) s3 <- s3[-.N]  # 本行目的：进行条件判断
  
  gene_map <- s3[, .(  # 本行目的：设置参数或赋值对象
    GENE_ID = trimws(as.character(geneID)),  # 本行目的：设置参数或赋值对象
    gene_symbol = trimws(as.character(geneSymble))  # 本行目的：设置参数或赋值对象
  )]  # 本行目的：执行当前流程语句
  
  gene_map <- gene_map[between(nchar(GENE_ID), 1L, 1e9)]  # 本行目的：设置参数或赋值对象
  gene_map[, GENE_ID := sub("\\.\\d+$", "", GENE_ID)]  # 本行目的：设置参数或赋值对象
  
  gene_map[is.na(gene_symbol), gene_symbol := ""]  # 本行目的：设置参数或赋值对象
  gene_map[, sym_len := nchar(gene_symbol)]  # 本行目的：设置参数或赋值对象
  setorder(gene_map, GENE_ID, -sym_len)  # 本行目的：执行当前流程语句
  gene_map <- unique(gene_map, by = "GENE_ID")  # 本行目的：设置参数或赋值对象
  gene_map[, sym_len := NULL]  # 本行目的：设置参数或赋值对象
  gene_map  # 本行目的：执行当前流程语句
}  # 本行目的：组织代码块结构

make_forest_plot <- function(key_table, top_table, out_dir, theme_pub) {  # 本行目的：定义函数入口
  plot_use <- if (!is.null(key_table) && nrow(key_table) > 0) "significant" else "tophits"  # 本行目的：进行条件判断
  plot_dat <- if (plot_use == "significant") copy(as.data.table(key_table)) else  # 本行目的：进行条件判断
    copy(as.data.table(top_table))  # 本行目的：执行当前流程语句
  
  if (!"gene_symbol" %in% names(plot_dat)) plot_dat[, gene_symbol := NA_character_]  # 本行目的：进行条件判断
  
  need_cols <- c("CELL_TYPE", "GENE_ID", "method", "nsnp", "pval", "fdr", "OR", "OR_lci95", "OR_uci95")  # 本行目的：设置参数或赋值对象
  miss_cols <- setdiff(need_cols, names(plot_dat))  # 本行目的：设置参数或赋值对象
  if (length(miss_cols) > 0) {  # 本行目的：组织代码块结构
    stop("Forest plot missing columns: ", paste(miss_cols, collapse = ", "),  # 本行目的：在异常条件下终止并提示
         "\n请确认主表已生成 OR/CI，且包含这些列。")  # 本行目的：执行当前流程语句
  }  # 本行目的：组织代码块结构
  
  num_cols <- intersect(c("OR", "OR_lci95", "OR_uci95", "pval", "fdr", "nsnp"), names(plot_dat))  # 本行目的：设置参数或赋值对象
  plot_dat[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]  # 本行目的：设置参数或赋值对象
  
  if (!"outcome" %in% names(plot_dat)) plot_dat[, outcome := NA_character_]  # 本行目的：进行条件判断
  if (!"id.outcome" %in% names(plot_dat)) plot_dat[, id.outcome := "outcome"]  # 本行目的：进行条件判断
  plot_dat[, outcome_show := fifelse(!is.na(outcome) & nzchar(outcome), outcome, id.outcome)]  # 本行目的：处理条件分支
  
  plot_dat <- plot_dat[  # 本行目的：设置参数或赋值对象
    is.finite(OR) & is.finite(OR_lci95) & is.finite(OR_uci95) &  # 本行目的：执行当前流程语句
      OR > 0 & OR_lci95 > 0 & OR_uci95 > 0  # 本行目的：执行当前流程语句
  ]  # 本行目的：执行当前流程语句
  if (nrow(plot_dat) == 0) stop("No valid rows left for forest plot after filtering OR/CI.")  # 本行目的：进行条件判断
  
  plot_dat[, label := ifelse(  # 本行目的：处理条件分支
    !is.na(gene_symbol) & nzchar(gene_symbol),  # 本行目的：执行当前流程语句
    paste0(CELL_TYPE, " | ", gene_symbol, " (", GENE_ID, ")"),  # 本行目的：执行当前流程语句
    paste0(CELL_TYPE, " | ", GENE_ID)  # 本行目的：执行当前流程语句
  )]  # 本行目的：执行当前流程语句
  
  plot_dat[, ann := paste0(  # 本行目的：设置参数或赋值对象
    method,  # 本行目的：执行当前流程语句
    " | nsnp=", nsnp,  # 本行目的：设置参数或赋值对象
    " | p=", formatC(pval, format = "e", digits = 2),  # 本行目的：设置参数或赋值对象
    " | FDR=", formatC(fdr, format = "e", digits = 2)  # 本行目的：设置参数或赋值对象
  )]  # 本行目的：执行当前流程语句
  
  plot_dat[, y_fac := paste0(label, "@@", outcome_show)]  # 本行目的：设置参数或赋值对象
  setorder(plot_dat, outcome_show, pval)  # 本行目的：执行当前流程语句
  lev <- plot_dat[, unique(y_fac)]  # 本行目的：设置参数或赋值对象
  plot_dat[, y_fac := factor(y_fac, levels = rev(lev))]  # 本行目的：设置参数或赋值对象
  
  plot_dat[, x_text := pmax(max(OR_uci95, na.rm = TRUE), 1) * 1.5, by = outcome_show]  # 本行目的：设置参数或赋值对象
  
  out_plot_data <- file.path(out_dir, "Fig_MR_forest_data.tsv")  # 本行目的：设置参数或赋值对象
  fwrite(plot_dat, out_plot_data, sep = "\t")  # 本行目的：写出表格文件
  
  p_forest <- ggplot2::ggplot(  # 本行目的：设置参数或赋值对象
    plot_dat,  # 本行目的：执行当前流程语句
    ggplot2::aes(x = OR, y = y_fac, xmin = OR_lci95, xmax = OR_uci95, color = method)  # 本行目的：设置参数或赋值对象
  ) +  # 本行目的：执行当前流程语句
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +  # 本行目的：设置参数或赋值对象
    ggplot2::geom_pointrange(fatten = 1.15) +  # 本行目的：设置参数或赋值对象
    ggplot2::geom_text(ggplot2::aes(x = x_text, label = ann),  # 本行目的：设置参数或赋值对象
                       color = "black", size = 2, hjust = 0) +  # 本行目的：设置参数或赋值对象
    ggplot2::scale_x_log10() +  # 本行目的：执行当前流程语句
    scico::scale_color_scico_d(palette = "batlow", begin = 0.1, end = 0.9) +  # 本行目的：设置参数或赋值对象
    ggplot2::scale_y_discrete(labels = function(x) sub("@@.*$", "", x)) +  # 本行目的：设置参数或赋值对象
    ggplot2::labs(  # 本行目的：执行当前流程语句
      x = "Odds ratio (log10 scale)",  # 本行目的：设置参数或赋值对象
      y = NULL,  # 本行目的：设置参数或赋值对象
      color = "MR method",  # 本行目的：设置参数或赋值对象
      title = "MR forest plot (cell type | gene → outcome)"  # 本行目的：设置参数或赋值对象
    ) +  # 本行目的：执行当前流程语句
    theme_pub(base = 10) +  # 本行目的：设置参数或赋值对象
    ggplot2::theme(  # 本行目的：执行当前流程语句
      panel.grid.minor = ggplot2::element_blank(),  # 本行目的：设置参数或赋值对象
      plot.margin = ggplot2::margin(5.5, 50, 5.5, 5.5),  # 本行目的：设置参数或赋值对象
      axis.text.y = ggplot2::element_text(size = 6)   # <- 加这一行  # 本行目的：设置参数或赋值对象
    ) +  # 本行目的：执行当前流程语句
    ggplot2::coord_cartesian(clip = "off") +  # 本行目的：设置参数或赋值对象
    ggplot2::facet_wrap(~ outcome_show, scales = "free_y", ncol = 1)  # 本行目的：设置参数或赋值对象
  
  out_png <- file.path(out_dir, "Fig_MR_forest.png")  # 本行目的：设置参数或赋值对象
  h <- min(26, max(7, 0.25 * nrow(plot_dat)))  # 本行目的：设置参数或赋值对象
  ggplot2::ggsave(out_png, p_forest, width = 14, height = h, dpi = 300)  # 本行目的：设置参数或赋值对象
  
  list(plot = p_forest, data = plot_dat, plot_use = plot_use, out_png = out_png, out_data = out_plot_data)  # 本行目的：设置参数或赋值对象
}  # 本行目的：组织代码块结构

qc_top_table <- function(top_table, out_dir) {  # 本行目的：定义函数入口
  top_dt <- as.data.table(copy(top_table))  # 本行目的：设置参数或赋值对象
  if (!"gene_symbol" %in% names(top_dt)) top_dt[, gene_symbol := NA_character_]  # 本行目的：进行条件判断
  
  top_dt[, GENE_ID := sub("\\.\\d+$", "", trimws(GENE_ID))]  # 本行目的：设置参数或赋值对象
  top_dt[, CELL_TYPE := trimws(CELL_TYPE)]  # 本行目的：设置参数或赋值对象
  
  if (!"outcome" %in% names(top_dt)) top_dt[, outcome := NA_character_]  # 本行目的：进行条件判断
  if (!"id.outcome" %in% names(top_dt)) top_dt[, id.outcome := NA_character_]  # 本行目的：进行条件判断
  top_dt[, outcome_show := fifelse(!is.na(outcome) & nzchar(outcome), outcome, id.outcome)]  # 本行目的：处理条件分支
  
  top_dt[, gene_show := fifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, GENE_ID)]  # 本行目的：处理条件分支
  
  gene_cell <- unique(top_dt[, .(GENE_ID, gene_show, CELL_TYPE)])  # 本行目的：设置参数或赋值对象
  gene_cell_stat <- gene_cell[, .(  # 本行目的：设置参数或赋值对象
    n_celltypes = uniqueN(CELL_TYPE),  # 本行目的：设置参数或赋值对象
    celltypes = paste(sort(unique(CELL_TYPE)), collapse = "; ")  # 本行目的：设置参数或赋值对象
  ), by = .(GENE_ID, gene_show)]  # 本行目的：设置参数或赋值对象
  setorder(gene_cell_stat, -n_celltypes, gene_show)  # 本行目的：执行当前流程语句
  
  gene_out <- unique(top_dt[, .(GENE_ID, gene_show, outcome_show)])  # 本行目的：设置参数或赋值对象
  gene_out_stat <- gene_out[, .(  # 本行目的：设置参数或赋值对象
    n_outcomes = uniqueN(outcome_show),  # 本行目的：设置参数或赋值对象
    outcomes = paste(sort(unique(outcome_show)), collapse = "; ")  # 本行目的：设置参数或赋值对象
  ), by = .(GENE_ID, gene_show)]  # 本行目的：设置参数或赋值对象
  setorder(gene_out_stat, -n_outcomes, gene_show)  # 本行目的：执行当前流程语句
  
  gene_joint <- merge(gene_cell_stat, gene_out_stat, by = c("GENE_ID", "gene_show"), all = TRUE)  # 本行目的：设置参数或赋值对象
  setorder(gene_joint, -n_outcomes, -n_celltypes, gene_show)  # 本行目的：执行当前流程语句
  
  out_joint <- file.path(out_dir, "QC_topTable_gene_across_celltypes_outcomes.tsv")  # 本行目的：设置参数或赋值对象
  fwrite(gene_joint, out_joint, sep = "\t")  # 本行目的：写出表格文件
  
  ct_stat <- top_dt[, .(n_hits = .N, n_genes = uniqueN(GENE_ID)), by = .(CELL_TYPE)][order(-n_hits)]  # 本行目的：设置参数或赋值对象
  out_ct <- file.path(out_dir, "QC_topTable_hits_by_celltype.tsv")  # 本行目的：设置参数或赋值对象
  fwrite(ct_stat, out_ct, sep = "\t")  # 本行目的：写出表格文件
  
  out_stat <- top_dt[, .(n_hits = .N, n_genes = uniqueN(GENE_ID)), by = .(outcome_show)][order(-n_hits)]  # 本行目的：设置参数或赋值对象
  out_outcome <- file.path(out_dir, "QC_topTable_hits_by_outcome.tsv")  # 本行目的：设置参数或赋值对象
  fwrite(out_stat, out_outcome, sep = "\t")  # 本行目的：写出表格文件
  
  ct_out <- top_dt[, .(n_hits = .N, n_genes = uniqueN(GENE_ID)),  # 本行目的：设置参数或赋值对象
                   by = .(outcome_show, CELL_TYPE)][order(outcome_show, -n_hits)]  # 本行目的：设置参数或赋值对象
  out_ct_out <- file.path(out_dir, "QC_topTable_hits_by_celltype_by_outcome.tsv")  # 本行目的：设置参数或赋值对象
  fwrite(ct_out, out_ct_out, sep = "\t")  # 本行目的：写出表格文件
  
  invisible(list(  # 本行目的：执行当前流程语句
    gene_joint = gene_joint,  # 本行目的：设置参数或赋值对象
    ct_stat = ct_stat,  # 本行目的：设置参数或赋值对象
    out_stat = out_stat,  # 本行目的：设置参数或赋值对象
    ct_out = ct_out,  # 本行目的：设置参数或赋值对象
    out_joint = out_joint,  # 本行目的：设置参数或赋值对象
    out_ct = out_ct,  # 本行目的：设置参数或赋值对象
    out_outcome = out_outcome,  # 本行目的：设置参数或赋值对象
    out_ct_out = out_ct_out  # 本行目的：设置参数或赋值对象
  ))  # 本行目的：组织代码块结构
}  # 本行目的：组织代码块结构

# ============================================================  # 本行目的：说明步骤或上下文
# 主函数：run_project()  # 本行目的：说明步骤或上下文
# 目的：你外部只改 keyword 或 fixed_ids 就能完整跑完并出表/出图  # 本行目的：说明步骤或上下文
# ============================================================  # 本行目的：说明步骤或上下文
run_project <- function(  # 本行目的：定义函数入口
    keyword,  # 本行目的：执行当前流程语句
    fixed_ids = character(),     # 指定 outcome id（不填则走搜索）  # 本行目的：设置参数或赋值对象
    use_search_fill = FALSE,     # TRUE：固定 + 搜索补齐；FALSE：只跑 fixed_ids  # 本行目的：设置参数或赋值对象
    n_search = 3,                # 搜索/补齐时取几个  # 本行目的：设置参数或赋值对象
    trait_regex = NULL,  # 本行目的：设置参数或赋值对象
    id_regex = NULL,  # 本行目的：设置参数或赋值对象
    population_regex = "European|EUR|UKB|UK Biobank|Finnish|FinnGen|Finland",  # 本行目的：设置参数或赋值对象
    min_year = NA_integer_,  # 本行目的：设置参数或赋值对象
    
    # 路径参数（你一般不用改，放这里做默认值）  # 本行目的：说明步骤或上下文
    pipeline_root = "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE",  # 本行目的：设置参数或赋值对象
    exp_rds = "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE/tmp_onek1k/vpub_s3_p5e8_kb10000_r2_0.001/mr_n982_from_s3_F10.rds",  # 本行目的：设置参数或赋值对象
    theme_file = "W:/fangxiaobin-CTRH-omics/000temple/scico_plot_theme.R",  # 本行目的：设置参数或赋值对象
    ao_cache_rds = "D:/fangxiaobin-CTRH-omics/disease_mr_keygene/Reference/MR/cache/available_outcomes.rds",  # 本行目的：设置参数或赋值对象
    vpub_file  = "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE/literature/V_publish.xlsx",  # 本行目的：设置参数或赋值对象
    vpub_sheet = "S3",  # 本行目的：设置参数或赋值对象
    
    # 分析参数  # 本行目的：说明步骤或上下文
    fdr_cut = 0.05,  # 本行目的：设置参数或赋值对象
    top_n_per_outcome = 30,  # 本行目的：设置参数或赋值对象
    
    # 重跑开关  # 本行目的：说明步骤或上下文
    rebuild_top3 = FALSE,  # 本行目的：设置参数或赋值对象
    rebuild_mr   = FALSE,  # 本行目的：设置参数或赋值对象
    rebuild_ao   = FALSE,  # 本行目的：设置参数或赋值对象
    
    # 其它  # 本行目的：说明步骤或上下文
    sleep = 10,  # 本行目的：设置参数或赋值对象
    do_forest_plot = TRUE,  # 本行目的：设置参数或赋值对象
    do_qc_top = TRUE,  # 本行目的：设置参数或赋值对象
    restore_wd = FALSE,  # 本行目的：设置参数或赋值对象
    verbose = TRUE  # 本行目的：设置参数或赋值对象
) {  # 本行目的：组织代码块结构
  suppressPackageStartupMessages({  # 本行目的：组织代码块结构
    library(data.table)  # 本行目的：加载依赖包
    library(TwoSampleMR)  # 本行目的：加载依赖包
    library(ieugwasr)  # 本行目的：加载依赖包
    library(ggplot2)  # 本行目的：加载依赖包
    library(scico)  # 本行目的：加载依赖包
    library(patchwork)  # 本行目的：加载依赖包
    library(readxl)  # 本行目的：加载依赖包
  })  # 本行目的：组织代码块结构
  
  # 目的：可选关闭旧图形设备，避免 ggsave/绘图冲突  # 本行目的：说明步骤或上下文
 # while (length(grDevices::dev.list()) > 0) grDevices::dev.off()  # 本行目的：说明步骤或上下文
  gc()  # 本行目的：执行当前流程语句
  
  fdr_tag <- gsub("\\.", "p", format(fdr_cut, scientific = FALSE, trim = TRUE))  # 本行目的：设置参数或赋值对象
  
  # step1：按 keyword 建目录  # 本行目的：说明步骤或上下文
  key_slug <- slugify(keyword)  # 本行目的：设置参数或赋值对象
  base_dir <- file.path(pipeline_root, key_slug)  # 本行目的：设置参数或赋值对象
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)  # 本行目的：创建输出目录
  
  old_wd <- getwd()  # 本行目的：设置或获取工作目录
  if (restore_wd) on.exit(setwd(old_wd), add = TRUE)  # 本行目的：进行条件判断
  setwd(base_dir)  # 本行目的：设置或获取工作目录
  if (verbose) message("当前项目工作目录：", normalizePath(getwd(), winslash = "/"))  # 本行目的：进行条件判断
  
  # step2：主题  # 本行目的：说明步骤或上下文
  if (!is.null(theme_file) && file.exists(theme_file)) source(theme_file)  # 本行目的：进行条件判断
  if (!exists("theme_pub")) theme_pub <- function(base = 11) ggplot2::theme_minimal(base_size = base)  # 本行目的：定义函数入口
  
  # step3：暴露  # 本行目的：说明步骤或上下文
  stopifnot(file.exists(exp_rds))  # 本行目的：执行当前流程语句
  mr_exp <- readRDS(exp_rds)  # 本行目的：读取输入或缓存数据
  mr_exp <- as.data.table(mr_exp)  # 本行目的：设置参数或赋值对象
  
  exposure_dat <- mr_exp[, .(  # 本行目的：设置参数或赋值对象
    SNP,  # 本行目的：执行当前流程语句
    beta.exposure = beta,  # 本行目的：设置参数或赋值对象
    se.exposure   = se,  # 本行目的：设置参数或赋值对象
    effect_allele.exposure = effect_allele,  # 本行目的：设置参数或赋值对象
    other_allele.exposure  = other_allele,  # 本行目的：设置参数或赋值对象
    eaf.exposure  = eaf,  # 本行目的：设置参数或赋值对象
    pval.exposure = pval,  # 本行目的：设置参数或赋值对象
    exposure      = exposure_id,  # 本行目的：设置参数或赋值对象
    id.exposure   = exposure_id,  # 本行目的：设置参数或赋值对象
    samplesize.exposure = N  # 本行目的：设置参数或赋值对象
  )]  # 本行目的：执行当前流程语句
  # ---- DEDUP exposure: (id.exposure, SNP) 必须唯一 ----  # 本行目的：说明步骤或上下文
  exp_dt <- as.data.table(exposure_dat)  # 本行目的：设置参数或赋值对象
  exp_dt <- unique(exp_dt, by = c("id.exposure", "SNP"))  # 本行目的：设置参数或赋值对象
  exposure_dat <- as.data.frame(exp_dt)  # 本行目的：设置参数或赋值对象
  # -----------------------------------------------  # 本行目的：说明步骤或上下文
  snps <- unique(exposure_dat$SNP)   # <- 这行保留/放在去重后  # 本行目的：设置参数或赋值对象
  
  
  if (verbose) {  # 本行目的：组织代码块结构
    message("Exposure rows: ", nrow(exposure_dat))  # 本行目的：输出运行信息
    message("Unique exposures: ", length(unique(exposure_dat$id.exposure)))  # 本行目的：输出运行信息
    message("Unique SNPs: ", length(snps))  # 本行目的：输出运行信息
  }  # 本行目的：组织代码块结构
  
  # step5：token  # 本行目的：说明步骤或上下文
  token <- ieugwasr::get_opengwas_jwt()  # 本行目的：设置参数或赋值对象
  Sys.setenv(OPENGWAS_JWT = token)  # 本行目的：设置参数或赋值对象
  stopifnot(nzchar(token))  # 本行目的：执行当前流程语句
  
  # step6：available_outcomes 缓存  # 本行目的：说明步骤或上下文
  dir.create(dirname(ao_cache_rds), showWarnings = FALSE, recursive = TRUE)  # 本行目的：创建输出目录
  if (file.exists(ao_cache_rds) && !rebuild_ao) {  # 本行目的：组织代码块结构
    ao <- readRDS(ao_cache_rds)  # 本行目的：读取输入或缓存数据
    if (verbose) message("available_outcomes 从缓存读取：", ao_cache_rds, " 行=", nrow(ao))  # 本行目的：进行条件判断
  } else {  # 本行目的：组织代码块结构
    ao <- TwoSampleMR::available_outcomes()  # 本行目的：设置参数或赋值对象
    saveRDS(ao, ao_cache_rds, compress = TRUE)  # 本行目的：保存RDS结果文件
    if (verbose) message("available_outcomes 已缓存：", ao_cache_rds, " 行=", nrow(ao))  # 本行目的：进行条件判断
  }  # 本行目的：组织代码块结构
  
  # step7：选择 outcomes（只改 fixed_ids/keyword 就能换）  # 本行目的：说明步骤或上下文
  mode <- if (length(fixed_ids) > 0 && !use_search_fill) "fixed" else  # 本行目的：进行条件判断
    if (length(fixed_ids) > 0 && use_search_fill) "mixed" else "search"  # 本行目的：进行条件判断
  
  out_sel <- file.path(getwd(), "opengwas_outcomes_selected.tsv")  # 本行目的：设置或获取工作目录
  
  top3 <- resolve_outcomes_meta(  # 本行目的：设置参数或赋值对象
    ao = ao,  # 本行目的：设置参数或赋值对象
    keyword = keyword,  # 本行目的：设置参数或赋值对象
    fixed_ids = fixed_ids,  # 本行目的：设置参数或赋值对象
    mode = mode,  # 本行目的：设置参数或赋值对象
    n = n_search,  # 本行目的：设置参数或赋值对象
    population_regex = population_regex,  # 本行目的：设置参数或赋值对象
    trait_regex = trait_regex,  # 本行目的：设置参数或赋值对象
    id_regex = id_regex,  # 本行目的：设置参数或赋值对象
    min_year = min_year,  # 本行目的：设置参数或赋值对象
    cache_tsv = out_sel,  # 本行目的：设置参数或赋值对象
    rebuild = rebuild_top3,  # 本行目的：设置参数或赋值对象
    verbose = verbose  # 本行目的：设置参数或赋值对象
  )  # 本行目的：执行当前流程语句
  
  # step8：ss（Steiger 用）  # 本行目的：说明步骤或上下文
  top3 <- as.data.table(top3)  # 本行目的：设置参数或赋值对象
  if (!"sample_size_num" %in% names(top3)) top3 <- calc_sample_size_num(top3)  # 本行目的：进行条件判断
  top3[, ss := fifelse(  # 本行目的：处理条件分支
    !is.na(sample_size_num), as.numeric(sample_size_num),  # 本行目的：执行当前流程语句
    suppressWarnings(as.numeric(ncase) + as.numeric(ncontrol))  # 本行目的：执行当前流程语句
  )]  # 本行目的：执行当前流程语句
  setorder(top3, -ss)  # 本行目的：执行当前流程语句
  
  if (verbose) {  # 本行目的：组织代码块结构
    message("\nSelected outcomes (sorted by ss):")  # 本行目的：输出运行信息
    print(top3[, .(id, trait, population, ss)])  # 本行目的：输出运行信息
  }  # 本行目的：组织代码块结构
  # --- NEW: export selected outcomes meta to CSV ---  # 本行目的：说明步骤或上下文
  out_sel_csv <- file.path(getwd(), "opengwas_outcomes_selected_main.csv")  # 本行目的：设置或获取工作目录
  
  top3_export <- as.data.table(top3)  # 本行目的：设置参数或赋值对象
  keep_cols <- intersect(  # 本行目的：设置参数或赋值对象
    c("id", "trait", "population", "category", "subcategory", "ontology", "year",  # 本行目的：执行当前流程语句
      "sample_size", "ncase", "ncontrol", "sample_size_num", "ss"),  # 本行目的：执行当前流程语句
    names(top3_export)  # 本行目的：执行当前流程语句
  )  # 本行目的：执行当前流程语句
  
  data.table::fwrite(top3_export[, ..keep_cols], out_sel_csv)  # 本行目的：写出表格文件
  if (verbose) message("[OK] wrote selected outcomes (CSV): ", out_sel_csv)  # 本行目的：进行条件判断
  # ============================================================  # 本行目的：说明步骤或上下文
  # step9：跑一个 outcome（闭包：自动拿到 exposure_dat/snps/fdr_cut/...）  # 本行目的：说明步骤或上下文
  # 目的：你的 MR、Steiger、缓存、输出逻辑保持不变  # 本行目的：说明步骤或上下文
  # ============================================================  # 本行目的：说明步骤或上下文
  run_mr_one_outcome <- function(outcome_id, outcome_ss, out_prefix) {  # 本行目的：定义函数入口
    
    # 9.1 输出文件名（MR结果/显著子集/细胞汇总/缓存）  # 本行目的：说明步骤或上下文
    out_all <- paste0(out_prefix, "_MR_all.tsv.gz")  # 本行目的：设置参数或赋值对象
    out_sig <- paste0(out_prefix, "_MR_FDR", fdr_tag, ".tsv.gz")  # 本行目的：设置参数或赋值对象
    out_ct  <- paste0(out_prefix, "_sig_counts_by_celltype.tsv")  # 本行目的：设置参数或赋值对象
    
    cache_mr_rds      <- paste0(out_prefix, "_MR_res.rds")  # 本行目的：设置参数或赋值对象
    cache_dat_rds     <- paste0(out_prefix, "_harmonised_dat2.rds")  # 本行目的：设置参数或赋值对象
    cache_outcome_rds <- paste0(out_prefix, "_outcome_raw.rds")  # 本行目的：设置参数或赋值对象
    
    # 9.2 MR结果缓存：如果已有 MR_all / MR_res，就直接读并返回（不再跑）  # 本行目的：说明步骤或上下文
    if (!rebuild_mr) {  # 本行目的：组织代码块结构
      if (file.exists(cache_mr_rds)) {  # 本行目的：组织代码块结构
        mr_res <- as.data.table(readRDS(cache_mr_rds))  # 本行目的：读取输入或缓存数据
        if (verbose) message("[CACHE] MR loaded: ", cache_mr_rds, " rows=", nrow(mr_res))  # 本行目的：进行条件判断
        return(list(  # 本行目的：返回结果对象
          mr = mr_res,  # 本行目的：设置参数或赋值对象
          dat = if (file.exists(cache_dat_rds)) readRDS(cache_dat_rds) else NULL  # 本行目的：进行条件判断
        ))  # 本行目的：组织代码块结构
      }  # 本行目的：组织代码块结构
      if (file.exists(out_all)) {  # 本行目的：组织代码块结构
        mr_res <- fread(out_all)  # 本行目的：读取输入或缓存数据
        if (verbose) message("[CACHE] MR loaded: ", out_all, " rows=", nrow(mr_res))  # 本行目的：进行条件判断
        return(list(  # 本行目的：返回结果对象
          mr = mr_res,  # 本行目的：设置参数或赋值对象
          dat = if (file.exists(cache_dat_rds)) readRDS(cache_dat_rds) else NULL  # 本行目的：进行条件判断
        ))  # 本行目的：组织代码块结构
      }  # 本行目的：组织代码块结构
    }  # 本行目的：组织代码块结构
    
    # 9.3 outcome 提取（优先读 outcome_raw 缓存；否则联网 extract 并缓存）  # 本行目的：说明步骤或上下文
    if (file.exists(cache_outcome_rds)) {  # 本行目的：组织代码块结构
      outcome_dat <- readRDS(cache_outcome_rds)  # 本行目的：读取输入或缓存数据
    } else {  # 本行目的：组织代码块结构
      args <- list(snps = snps, outcomes = outcome_id, proxies = FALSE)  # 本行目的：设置参数或赋值对象
      if ("splitsize" %in% names(formals(TwoSampleMR::extract_outcome_data))) {  # 本行目的：组织代码块结构
        args$splitsize <- 500L  # 本行目的：设置参数或赋值对象
      }  # 本行目的：组织代码块结构
      
      outcome_dat <- NULL  # 本行目的：设置参数或赋值对象
      for (k in 1:3) {  # 本行目的：组织代码块结构
        outcome_dat <- tryCatch(  # 本行目的：设置参数或赋值对象
          do.call(TwoSampleMR::extract_outcome_data, args),  # 本行目的：执行当前流程语句
          error = function(e) e  # 本行目的：设置参数或赋值对象
        )  # 本行目的：执行当前流程语句
        
        if (!inherits(outcome_dat, "error")) break  # 本行目的：进行条件判断
        
        msg <- conditionMessage(outcome_dat)  # 本行目的：设置参数或赋值对象
        message("[WARN] extract_outcome_data failed (try ", k, "/3): ", msg)  # 本行目的：输出运行信息
        
        # ✅ 如果遇到 401（Unauthorized），自动刷新 token 再重试  # 本行目的：说明步骤或上下文
        if (grepl("\\b401\\b", msg)) {  # 本行目的：组织代码块结构
          message("[AUTH] Detected 401. Refreshing OpenGWAS JWT...")  # 本行目的：输出运行信息
          token <- ieugwasr::get_opengwas_jwt()  # 本行目的：设置参数或赋值对象
          Sys.setenv(OPENGWAS_JWT = token)  # 本行目的：设置参数或赋值对象
        }  # 本行目的：组织代码块结构
        
        Sys.sleep(3 * k)  # 本行目的：执行当前流程语句
      }  # 本行目的：组织代码块结构
      
      if (inherits(outcome_dat, "error")) stop(conditionMessage(outcome_dat))  # 本行目的：进行条件判断
      
      
      saveRDS(outcome_dat, cache_outcome_rds, compress = TRUE)  # 本行目的：保存RDS结果文件
    }  # 本行目的：组织代码块结构
    # ---- DEDUP outcome: (SNP, id.outcome) 必须唯一 ----  # 本行目的：说明步骤或上下文
    out_dt <- as.data.table(outcome_dat)  # 本行目的：设置参数或赋值对象
    if (all(c("SNP", "id.outcome") %in% names(out_dt))) {  # 本行目的：组织代码块结构
      out_dt <- unique(out_dt, by = c("SNP", "id.outcome"))  # 本行目的：设置参数或赋值对象
      outcome_dat <- as.data.frame(out_dt)  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    # --------------------------------------------------  # 本行目的：说明步骤或上下文
    
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {  # 本行目的：组织代码块结构
      warning("No outcome data extracted for outcome_id = ", outcome_id)  # 本行目的：给出警告但不中断
      return(NULL)  # 本行目的：返回结果对象
    }  # 本行目的：组织代码块结构
    
    # 9.4 harmonise：等位基因对齐/回文位点处理  # 本行目的：说明步骤或上下文
    dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2)  # 本行目的：设置参数或赋值对象
    
    # outcome 样本量补齐（Steiger会用）  # 本行目的：说明步骤或上下文
    if (!"samplesize.outcome" %in% names(dat)) dat$samplesize.outcome <- NA  # 本行目的：进行条件判断
    dat$samplesize.outcome[is.na(dat$samplesize.outcome)] <- outcome_ss  # 本行目的：设置参数或赋值对象
    
    # ------------------------------------------------------------  # 本行目的：说明步骤或上下文
    # 关键修复点：dat 是 data.frame，统计时先转 data.table 再用 by=  # 本行目的：说明步骤或上下文
    # ------------------------------------------------------------  # 本行目的：说明步骤或上下文
    dat_dt <- as.data.table(dat)  # 本行目的：设置参数或赋值对象
    # ---- DEDUP#1 harmonise: (SNP, exposure, id.outcome) 必须唯一 ----  # 本行目的：说明步骤或上下文
    if (all(c("SNP","exposure","id.outcome") %in% names(dat_dt))) {  # 本行目的：组织代码块结构
      if ("mr_keep" %in% names(dat_dt)) {  # 本行目的：组织代码块结构
        dat_dt[, mr_keep := as.logical(mr_keep)]  # 本行目的：设置参数或赋值对象
        dat_dt[, mr_keep_int := as.integer(mr_keep)]  # 本行目的：设置参数或赋值对象
        data.table::setorder(dat_dt, exposure, id.outcome, SNP, -mr_keep_int)  # 本行目的：执行当前流程语句
        dat_dt[, mr_keep_int := NULL]  # 本行目的：设置参数或赋值对象
      } else {  # 本行目的：组织代码块结构
        data.table::setorder(dat_dt, exposure, id.outcome, SNP)  # 本行目的：执行当前流程语句
      }  # 本行目的：组织代码块结构
      dat_dt <- unique(dat_dt, by = c("SNP","exposure","id.outcome"))  # 本行目的：设置参数或赋值对象
      dat <- as.data.frame(dat_dt)  # 关键：让 Steiger 用到去重后的 dat  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    # ----------------------------------------------------------------  # 本行目的：说明步骤或上下文
    dat_dt <- as.data.table(dat)   # 确保 harm_sum 用的是去重后的 dat  # 本行目的：设置参数或赋值对象
    
    # 9.5 统计：harmonise 后每个 exposure 实际保留的 SNP 数（mr_keep==TRUE）  # 本行目的：说明步骤或上下文
    harm_sum <- NULL  # 本行目的：设置参数或赋值对象
    if ("mr_keep" %in% names(dat_dt)) {  # 本行目的：组织代码块结构
      harm_sum <- dat_dt[mr_keep %in% TRUE,  # 本行目的：设置参数或赋值对象
                         .(nsnp_harmonised = uniqueN(SNP)),  # 本行目的：设置参数或赋值对象
                         by = "exposure"]  # 本行目的：设置参数或赋值对象
    } else {  # 本行目的：组织代码块结构
      harm_sum <- dat_dt[, .(nsnp_harmonised = uniqueN(SNP)), by = "exposure"]  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    
    # 9.6 Steiger filtering：方向性检验（失败则跳过）  # 本行目的：说明步骤或上下文
    steiger_sum <- NULL  # 本行目的：设置参数或赋值对象
    steiger_applied <- FALSE  # 本行目的：设置参数或赋值对象
    
    dat2_dt <- tryCatch({  # 本行目的：组织代码块结构
      tmp <- TwoSampleMR::steiger_filtering(dat)   # 通常返回 data.frame  # 本行目的：设置参数或赋值对象
      tmp_dt <- as.data.table(tmp)  # 本行目的：设置参数或赋值对象
      
      if ("steiger_dir" %in% names(tmp_dt)) {  # 本行目的：组织代码块结构
        steiger_applied <- TRUE  # 本行目的：设置参数或赋值对象
        
        # 统计：Steiger 通过的 SNP 数（按 exposure）  # 本行目的：说明步骤或上下文
        if ("mr_keep" %in% names(tmp_dt)) {  # 本行目的：组织代码块结构
          steiger_sum <- tmp_dt[mr_keep %in% TRUE,  # 本行目的：设置参数或赋值对象
                                .(nsnp_steiger_keep = uniqueN(SNP[steiger_dir %in% TRUE])),  # 本行目的：设置参数或赋值对象
                                by = "exposure"]  # 本行目的：设置参数或赋值对象
        } else {  # 本行目的：组织代码块结构
          steiger_sum <- tmp_dt[, .(nsnp_steiger_keep = uniqueN(SNP[steiger_dir %in% TRUE])),  # 本行目的：设置参数或赋值对象
                                by = "exposure"]  # 本行目的：设置参数或赋值对象
        }  # 本行目的：组织代码块结构
        
        # 用于 MR 的数据：只保留 steiger_dir==TRUE  # 本行目的：说明步骤或上下文
        tmp_dt <- tmp_dt[steiger_dir %in% TRUE]  # 本行目的：设置参数或赋值对象
      }  # 本行目的：组织代码块结构
      
      tmp_dt  # 本行目的：执行当前流程语句
    }, error = function(e) {  # 本行目的：组织代码块结构
      message("[WARN] Steiger filtering failed, continue without it: ", e$message)  # 本行目的：输出运行信息
      dat_dt  # 本行目的：执行当前流程语句
    })  # 本行目的：组织代码块结构
    # ---- DEDUP#2 steiger: (SNP, exposure, id.outcome) 必须唯一 ----  # 本行目的：说明步骤或上下文
    if (all(c("SNP","exposure","id.outcome") %in% names(dat2_dt))) {  # 本行目的：组织代码块结构
      if ("mr_keep" %in% names(dat2_dt)) {  # 本行目的：组织代码块结构
        dat2_dt[, mr_keep := as.logical(mr_keep)]  # 本行目的：设置参数或赋值对象
        dat2_dt[, mr_keep_int := as.integer(mr_keep)]  # 本行目的：设置参数或赋值对象
        data.table::setorder(dat2_dt, exposure, id.outcome, SNP, -mr_keep_int)  # 本行目的：执行当前流程语句
        dat2_dt[, mr_keep_int := NULL]  # 本行目的：设置参数或赋值对象
      } else {  # 本行目的：组织代码块结构
        data.table::setorder(dat2_dt, exposure, id.outcome, SNP)  # 本行目的：执行当前流程语句
      }  # 本行目的：组织代码块结构
      dat2_dt <- unique(dat2_dt, by = c("SNP","exposure","id.outcome"))  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    # ---------------------------------------------------------------  # 本行目的：说明步骤或上下文
    
    if (nrow(dat2_dt) == 0) {  # 本行目的：组织代码块结构
      warning("All instruments removed after harmonise/Steiger for outcome_id = ", outcome_id)  # 本行目的：给出警告但不中断
      return(NULL)  # 本行目的：返回结果对象
    }  # 本行目的：组织代码块结构
    
    # 9.7 MR：单工具 Wald；多工具 IVW  # 本行目的：说明步骤或上下文
    # 为兼容性，喂给 TwoSampleMR::mr 时用 data.frame  # 本行目的：说明步骤或上下文
    dat2 <- as.data.frame(dat2_dt)  # 本行目的：设置参数或赋值对象
    
    mr_res <- TwoSampleMR::mr(dat2, method_list = c("mr_wald_ratio", "mr_ivw"))  # 本行目的：设置参数或赋值对象
    mr_res <- as.data.table(mr_res)  # 本行目的：设置参数或赋值对象
    
    # 每个 exposure 只保留“该用的那一个方法”  # 本行目的：说明步骤或上下文
    mr_res <- mr_res[  # 本行目的：设置参数或赋值对象
      (nsnp == 1 & method == "Wald ratio") |  # 本行目的：设置参数或赋值对象
        (nsnp > 1 & method == "Inverse variance weighted")  # 本行目的：设置参数或赋值对象
    ]  # 本行目的：执行当前流程语句
    if (nrow(mr_res) == 0) {  # 本行目的：组织代码块结构
      warning("No MR results left after method selection for outcome_id = ", outcome_id)  # 本行目的：给出警告但不中断
      return(NULL)  # 本行目的：返回结果对象
    }  # 本行目的：组织代码块结构
    
    # 9.8 BH-FDR + 拆 CELL_TYPE/GENE_ID  # 本行目的：说明步骤或上下文
    mr_res <- unique(mr_res, by = c("exposure", "id.outcome", "method"))  # 本行目的：设置参数或赋值对象
    mr_res[, fdr := p.adjust(pval, method = "BH")]  # 本行目的：设置参数或赋值对象
    
    mr_res[, c("CELL_TYPE", "GENE_ID") := tstrsplit(exposure, "\\|")]  # 本行目的：设置参数或赋值对象
    
    # 9.9 ★把 harmonise/Steiger 的统计 merge 进 MR 结果  # 本行目的：说明步骤或上下文
    mr_res <- merge(mr_res, harm_sum, by = "exposure", all.x = TRUE)  # 本行目的：设置参数或赋值对象
    if (!is.null(steiger_sum) && nrow(steiger_sum) > 0) {  # 本行目的：组织代码块结构
      mr_res <- merge(mr_res, steiger_sum, by = "exposure", all.x = TRUE)  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    
    # steiger 保留比例（如果 steiger 没跑或没有列，就会是 NA）  # 本行目的：说明步骤或上下文
    if (!"nsnp_steiger_keep" %in% names(mr_res)) mr_res[, nsnp_steiger_keep := NA_real_]  # 本行目的：进行条件判断
    if (!"nsnp_harmonised" %in% names(mr_res)) mr_res[, nsnp_harmonised := NA_real_]  # 本行目的：进行条件判断
    mr_res[, steiger_keep_rate := nsnp_steiger_keep / nsnp_harmonised]  # 本行目的：设置参数或赋值对象
    mr_res[, steiger_applied := steiger_applied]  # 本行目的：设置参数或赋值对象
    
    # 9.10 输出 + 缓存（下次命中缓存就不再跑）  # 本行目的：说明步骤或上下文
    fwrite(mr_res, out_all, sep = "\t", compress = "gzip")  # 本行目的：写出表格文件
    fwrite(mr_res[mr_res$fdr <= fdr_cut], out_sig, sep = "\t", compress = "gzip")  # 本行目的：写出表格文件
    
    sig_ct <- mr_res[fdr <= fdr_cut,  # 本行目的：设置参数或赋值对象
                     .(n_sig = uniqueN(exposure)),  # 本行目的：设置参数或赋值对象
                     by = "CELL_TYPE"][order(-n_sig)]  # 本行目的：设置参数或赋值对象
    fwrite(sig_ct, out_ct, sep = "\t")  # 本行目的：写出表格文件
    
    saveRDS(mr_res, cache_mr_rds, compress = TRUE)  # 本行目的：保存RDS结果文件
    saveRDS(dat2_dt, cache_dat_rds, compress = TRUE)  # 缓存 SNP-level（data.table 也行）  # 本行目的：保存RDS结果文件
    
    if (verbose) {  # 本行目的：组织代码块结构
      message("\n[OK] outcome: ", outcome_id)  # 本行目的：输出运行信息
      message("  harmonised rows: ", nrow(dat2_dt))  # 本行目的：输出运行信息
      message("  MR tests: ", nrow(mr_res))  # 本行目的：输出运行信息
      message("  FDR<= ", fdr_cut, ": ", sum(mr_res$fdr <= fdr_cut))  # 本行目的：输出运行信息
      message("  wrote:\n   ", out_all, "\n   ", out_sig, "\n   ", out_ct)  # 本行目的：输出运行信息
    }  # 本行目的：组织代码块结构
    
    list(mr = mr_res, dat = dat2_dt)  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  # step10：循环跑 outcomes  # 本行目的：说明步骤或上下文
  # step10：循环跑 outcomes（每个 outcome 一个子目录）  # 本行目的：说明步骤或上下文
  results <- list()  # 本行目的：设置参数或赋值对象
  for (i in seq_len(nrow(top3))) {  # 本行目的：组织代码块结构
    outcome_id <- top3$id[i]  # 本行目的：设置参数或赋值对象
    outcome_ss <- top3$ss[i]  # 本行目的：设置参数或赋值对象
    
    outcome_dir <- file.path(getwd(), paste0("outcome_", outcome_id))  # 本行目的：设置或获取工作目录
    dir.create(outcome_dir, showWarnings = FALSE, recursive = TRUE)  # 本行目的：创建输出目录
    
    out_prefix <- file.path(outcome_dir, paste0("MR_", outcome_id))  # 本行目的：设置参数或赋值对象
    
    results[[outcome_id]] <- run_mr_one_outcome(outcome_id, outcome_ss, out_prefix)  # 本行目的：设置参数或赋值对象
    Sys.sleep(sleep)  # 本行目的：执行当前流程语句
  }  # 本行目的：组织代码块结构
  
  # step11：合并 MR_all  # 本行目的：说明步骤或上下文
  # step11：合并 MR_all（递归搜子目录）  # 本行目的：说明步骤或上下文
  mr_files <- list.files(getwd(), pattern = "^MR_.*_MR_all\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)  # 本行目的：设置或获取工作目录
  mr_files <- unique(normalizePath(mr_files, winslash = "/", mustWork = TRUE))  # 本行目的：设置参数或赋值对象
  stopifnot(length(mr_files) > 0)  # 本行目的：执行当前流程语句
  
  mr_all <- rbindlist(lapply(mr_files, fread), fill = TRUE)  # 本行目的：设置参数或赋值对象
  mr_all <- unique(mr_all, by = c("exposure", "id.outcome", "method", "b", "se", "pval"))  # 本行目的：设置参数或赋值对象
  
  if (verbose) message("[INFO] merged MR_all files: ", length(mr_files), " total rows=", nrow(mr_all))  # 本行目的：进行条件判断
  
  # step12：主表 + TopHits  # 本行目的：说明步骤或上下文
  num_cols <- intersect(c("b", "se", "pval", "fdr", "nsnp", "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate"), names(mr_all))  # 本行目的：设置参数或赋值对象
  mr_all[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]  # 本行目的：设置参数或赋值对象
  
  if (!all(c("CELL_TYPE", "GENE_ID") %in% names(mr_all))) {  # 本行目的：组织代码块结构
    src <- if ("exposure" %in% names(mr_all)) "exposure" else "id.exposure"  # 本行目的：进行条件判断
    mr_all[, c("CELL_TYPE", "GENE_ID") := tstrsplit(get(src), "\\|")]  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  mr_all[, `:=`(  # 本行目的：设置参数或赋值对象
    OR       = exp(b),  # 本行目的：设置参数或赋值对象
    OR_lci95 = exp(b - 1.96 * se),  # 本行目的：设置参数或赋值对象
    OR_uci95 = exp(b + 1.96 * se),  # 本行目的：设置参数或赋值对象
    OR_CI    = sprintf("%.3f (%.3f–%.3f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))  # 本行目的：设置参数或赋值对象
  )]  # 本行目的：执行当前流程语句
  
  key_cols <- intersect(c(  # 本行目的：设置参数或赋值对象
    "CELL_TYPE", "GENE_ID",  # 本行目的：执行当前流程语句
    "id.outcome", "outcome",  # 本行目的：执行当前流程语句
    "method", "nsnp",  # 本行目的：执行当前流程语句
    "OR_CI", "OR", "OR_lci95", "OR_uci95",  # 本行目的：执行当前流程语句
    "b", "se", "pval", "fdr",  # 本行目的：执行当前流程语句
    "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate", "steiger_applied"  # 本行目的：执行当前流程语句
  ), names(mr_all))  # 本行目的：执行当前流程语句
  
  key_table <- mr_all[is.finite(fdr) & fdr <= fdr_cut, ..key_cols]  # 本行目的：设置参数或赋值对象
  out_key <- file.path(getwd(), paste0("Table_MR_key_FDR", fdr_tag, ".tsv"))  # 本行目的：设置或获取工作目录
  fwrite(key_table, out_key, sep = "\t")  # 本行目的：写出表格文件
  if (verbose) message("[OK] wrote (MAIN TABLE): ", out_key, " rows=", nrow(key_table), " (FDR<=", fdr_cut, ")")  # 本行目的：进行条件判断
  
  top_table <- mr_all[is.finite(pval)][order(id.outcome, pval)][, head(.SD, top_n_per_outcome), by = id.outcome][, ..key_cols]  # 本行目的：设置参数或赋值对象
  out_top <- file.path(getwd(), "Table_MR_key_TopHits.tsv")  # 本行目的：设置或获取工作目录
  fwrite(top_table, out_top, sep = "\t")  # 本行目的：写出表格文件
  if (verbose) message("[OK] wrote (TOP HITS): ", out_top, " rows=", nrow(top_table))  # 本行目的：进行条件判断
  
  # step12.7：加 gene_symbol  # 本行目的：说明步骤或上下文
  gene_map <- build_gene_map(vpub_file, vpub_sheet)  # 本行目的：设置参数或赋值对象
  
  if (nrow(key_table) > 0) {  # 本行目的：组织代码块结构
    key_table[, GENE_ID := sub("\\.\\d+$", "", GENE_ID)]  # 本行目的：设置参数或赋值对象
    key_table <- merge(key_table, gene_map, by = "GENE_ID", all.x = TRUE)  # 本行目的：设置参数或赋值对象
    front <- intersect(c("CELL_TYPE", "gene_symbol", "GENE_ID"), names(key_table))  # 本行目的：设置参数或赋值对象
    setcolorder(key_table, c(front, setdiff(names(key_table), front)))  # 本行目的：执行当前流程语句
    
    out_key2 <- file.path(getwd(), paste0("Table_MR_key_FDR", fdr_tag, "_with_symbol.tsv"))  # 本行目的：设置或获取工作目录
    fwrite(key_table, out_key2, sep = "\t")  # 本行目的：写出表格文件
    if (verbose) message("[OK] wrote: ", out_key2, " rows=", nrow(key_table))  # 本行目的：进行条件判断
  } else {  # 本行目的：组织代码块结构
    out_key2 <- NA_character_  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  top_table[, GENE_ID := sub("\\.\\d+$", "", GENE_ID)]  # 本行目的：设置参数或赋值对象
  top_table <- merge(top_table, gene_map, by = "GENE_ID", all.x = TRUE)  # 本行目的：设置参数或赋值对象
  front <- intersect(c("CELL_TYPE", "gene_symbol", "GENE_ID"), names(top_table))  # 本行目的：设置参数或赋值对象
  setcolorder(top_table, c(front, setdiff(names(top_table), front)))  # 本行目的：执行当前流程语句
  
  out_top2 <- file.path(getwd(), "Table_MR_key_TopHits_with_symbol.tsv")  # 本行目的：设置或获取工作目录
  fwrite(top_table, out_top2, sep = "\t")  # 本行目的：写出表格文件
  if (verbose) message("[OK] wrote: ", out_top2, " rows=", nrow(top_table))  # 本行目的：进行条件判断
  # ============================================================  # 本行目的：说明步骤或上下文
  # ============================================================  # 本行目的：说明步骤或上下文
  # step12.7.9：Per-outcome forest plots（每个 outcome 一个森林图）  # 本行目的：说明步骤或上下文
  # 位置：gene_symbol merge 完成后、决策表/总森林图之前  # 本行目的：说明步骤或上下文
  # 输出：  # 本行目的：说明步骤或上下文
  #   - <base_dir>/outcome_<id>/Fig_MR_forest.png  # 本行目的：说明步骤或上下文
  #   - <base_dir>/outcome_<id>/Fig_MR_forest_data.tsv  # 本行目的：说明步骤或上下文
  #   - <base_dir>/outcome_<id>/steiger_only/Fig_MR_forest.png （若有）  # 本行目的：说明步骤或上下文
  # 并在 RStudio Plots 里依次显示  # 本行目的：说明步骤或上下文
  # ============================================================  # 本行目的：说明步骤或上下文
  
  forest_by_outcome <- list()  # 本行目的：设置参数或赋值对象
  forest_steiger_by_outcome <- list()  # 本行目的：设置参数或赋值对象
  
  if (isTRUE(do_forest_plot)) {  # 本行目的：组织代码块结构
    
    # 用你选中的 outcomes 顺序画图（更可控）  # 本行目的：说明步骤或上下文
    outcome_ids <- as.character(top3$id)  # 本行目的：设置参数或赋值对象
    
    for (oid in outcome_ids) {  # 本行目的：组织代码块结构
      
      out_dir_i <- file.path(getwd(), paste0("outcome_", oid))  # 本行目的：设置或获取工作目录
      dir.create(out_dir_i, showWarnings = FALSE, recursive = TRUE)  # 本行目的：创建输出目录
      
      # 子集：该 outcome 的 key_table / top_table  # 本行目的：说明步骤或上下文
      kt_i <- data.table::as.data.table(key_table)[id.outcome %in% oid]  # 本行目的：设置参数或赋值对象
      tt_i <- data.table::as.data.table(top_table)[id.outcome %in% oid]  # 本行目的：设置参数或赋值对象
      
      # 如果该 outcome 没有任何行，跳过  # 本行目的：说明步骤或上下文
      if (nrow(kt_i) == 0 && nrow(tt_i) == 0) next  # 本行目的：进行条件判断
      
      # 画该 outcome 的主森林图（优先 kt_i，否则 tt_i）  # 本行目的：说明步骤或上下文
      forest_i <- make_forest_plot(  # 本行目的：设置参数或赋值对象
        key_table = if (nrow(kt_i) > 0) kt_i else NULL,  # 本行目的：进行条件判断
        top_table = tt_i,  # 本行目的：设置参数或赋值对象
        out_dir = out_dir_i,  # 本行目的：设置参数或赋值对象
        theme_pub = theme_pub  # 本行目的：设置参数或赋值对象
      )  # 本行目的：执行当前流程语句
      forest_by_outcome[[oid]] <- forest_i  # 本行目的：设置参数或赋值对象
      
      if (verbose) message("[OK] saved (per-outcome forest): ", forest_i$out_png)  # 本行目的：进行条件判断
      if (interactive() && !is.null(forest_i$plot)) print(forest_i$plot)  # 本行目的：进行条件判断
      
      # ---- Steiger-only（该 outcome）----  # 本行目的：说明步骤或上下文
      # 只在有 steiger_applied 列时做筛选  # 本行目的：说明步骤或上下文
      kt_s <- NULL  # 本行目的：设置参数或赋值对象
      if ("steiger_applied" %in% names(kt_i)) {  # 本行目的：组织代码块结构
        kt_i[, steiger_applied := as.logical(steiger_applied)]  # 本行目的：设置参数或赋值对象
        kt_s <- kt_i[steiger_applied %in% TRUE]  # 本行目的：设置参数或赋值对象
      }  # 本行目的：组织代码块结构
      
      tt_s <- tt_i  # 本行目的：设置参数或赋值对象
      if ("steiger_applied" %in% names(tt_s)) {  # 本行目的：组织代码块结构
        tt_s[, steiger_applied := as.logical(steiger_applied)]  # 本行目的：设置参数或赋值对象
        tt_s <- tt_s[steiger_applied %in% TRUE]  # 本行目的：设置参数或赋值对象
      } else {  # 本行目的：组织代码块结构
        tt_s <- tt_s[0]  # 本行目的：设置参数或赋值对象
      }  # 本行目的：组织代码块结构
      
      if (!((is.null(kt_s) || nrow(kt_s) == 0) && nrow(tt_s) == 0)) {  # 本行目的：组织代码块结构
        out_dir_s <- file.path(out_dir_i, "steiger_only")  # 本行目的：设置参数或赋值对象
        dir.create(out_dir_s, showWarnings = FALSE, recursive = TRUE)  # 本行目的：创建输出目录
        
        forest_s <- make_forest_plot(  # 本行目的：设置参数或赋值对象
          key_table = if (!is.null(kt_s) && nrow(kt_s) > 0) kt_s else NULL,  # 本行目的：进行条件判断
          top_table = tt_s,  # 本行目的：设置参数或赋值对象
          out_dir = out_dir_s,  # 本行目的：设置参数或赋值对象
          theme_pub = theme_pub  # 本行目的：设置参数或赋值对象
        )  # 本行目的：执行当前流程语句
        forest_steiger_by_outcome[[oid]] <- forest_s  # 本行目的：设置参数或赋值对象
        
        if (verbose) message("[OK] saved (per-outcome Steiger-only): ", forest_s$out_png)  # 本行目的：进行条件判断
        if (interactive() && !is.null(forest_s$plot)) print(forest_s$plot)  # 本行目的：进行条件判断
      } else {  # 本行目的：组织代码块结构
        if (verbose) message("[INFO] No Steiger-applied rows for per-outcome Steiger-only: ", oid)  # 本行目的：进行条件判断
      }  # 本行目的：组织代码块结构
    }  # 本行目的：组织代码块结构
  }  # 本行目的：组织代码块结构
  
  # ============================================================  # 本行目的：说明步骤或上下文
  # step12.8：少而稳 - Go/No-go 决策表（基于 key_table with symbol）  # 本行目的：说明步骤或上下文
  # 位置：gene_symbol merge 之后、森林图之前  # 本行目的：说明步骤或上下文
  # ============================================================  # 本行目的：说明步骤或上下文
  decision_dir <- getwd()  # 本行目的：设置或获取工作目录
  
  # 1) 只从“显著表”出发（少而稳核心）  # 本行目的：说明步骤或上下文
  dec_src <- NULL  # 本行目的：设置参数或赋值对象
  if (exists("key_table") && is.data.frame(key_table) && nrow(key_table) > 0) {  # 本行目的：组织代码块结构
    dec_src <- data.table::as.data.table(key_table)  # 本行目的：设置参数或赋值对象
  } else {  # 本行目的：组织代码块结构
    # 少而稳：显著表为空，则 A_go 为空；仍可输出空表做记录  # 本行目的：说明步骤或上下文
    dec_src <- data.table::data.table()  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  # 2) 设定保守阈值（你可后面再调）  # 本行目的：说明步骤或上下文
  min_nsnp <- 3L  # 本行目的：设置参数或赋值对象
  min_steiger_rate <- 0.70  # 本行目的：设置参数或赋值对象
  max_ci_ratio <- 3.0  # 本行目的：设置参数或赋值对象
  min_effect_log_or <- log(1.1)  # 软门槛，用于排序，不一定作为硬过滤  # 本行目的：设置参数或赋值对象
  
  # 3) 生成审计字段与失败原因  # 本行目的：说明步骤或上下文
  decision_audit <- copy(dec_src)  # 本行目的：设置参数或赋值对象
  if ("steiger_applied" %in% names(decision_audit)) decision_audit[, steiger_applied := as.logical(steiger_applied)]  # 本行目的：进行条件判断
  if (nrow(decision_audit) > 0) {  # 本行目的：组织代码块结构
    # 确保数值列是数值  # 本行目的：说明步骤或上下文
    num_cols <- intersect(  # 本行目的：设置参数或赋值对象
      c("fdr", "pval", "OR", "OR_lci95", "OR_uci95", "nsnp",  # 本行目的：执行当前流程语句
        "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate"),  # 本行目的：执行当前流程语句
      names(decision_audit)  # 本行目的：执行当前流程语句
    )  # 本行目的：执行当前流程语句
    decision_audit[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]  # 本行目的：设置参数或赋值对象
    
    # 逐条门槛（硬）  # 本行目的：说明步骤或上下文
    decision_audit[, pass_fdr := is.finite(fdr) & fdr <= fdr_cut]  # 本行目的：设置参数或赋值对象
    decision_audit[, pass_nsnp := is.finite(nsnp) & nsnp >= min_nsnp &  # 本行目的：设置参数或赋值对象
                     is.finite(nsnp_harmonised) & nsnp_harmonised >= min_nsnp]  # 本行目的：设置参数或赋值对象
    
    decision_audit[, pass_steiger := (steiger_applied %in% TRUE) &  # 本行目的：设置参数或赋值对象
                     is.finite(steiger_keep_rate) & steiger_keep_rate >= min_steiger_rate]  # 本行目的：设置参数或赋值对象
    
    decision_audit[, ci_ratio := OR_uci95 / OR_lci95]  # 本行目的：设置参数或赋值对象
    decision_audit[, pass_ci := is.finite(OR_lci95) & is.finite(OR_uci95) &  # 本行目的：设置参数或赋值对象
                     OR_lci95 > 0 & OR_uci95 > 0 & is.finite(ci_ratio) & ci_ratio <= max_ci_ratio]  # 本行目的：设置参数或赋值对象
    
    # 软门槛：效应大小（主要用于排序）  # 本行目的：说明步骤或上下文
    decision_audit[, effect_log_or := abs(log(OR))]  # 本行目的：设置参数或赋值对象
    decision_audit[, pass_effect := is.finite(effect_log_or) & effect_log_or >= min_effect_log_or]  # 本行目的：设置参数或赋值对象
    
    # 汇总失败原因（便于“少而稳”解释）  # 本行目的：说明步骤或上下文
    decision_audit[, fail_reason := ""]  # 本行目的：设置参数或赋值对象
    decision_audit[(pass_fdr == FALSE), fail_reason := paste0(fail_reason, "FDR;")]  # 本行目的：设置参数或赋值对象
    decision_audit[(pass_fdr == TRUE)  & (pass_nsnp == FALSE), fail_reason := paste0(fail_reason, "Low_nsnp;")]  # 本行目的：设置参数或赋值对象
    decision_audit[(pass_fdr == TRUE)  & (pass_nsnp == TRUE) & (pass_steiger == FALSE), fail_reason := paste0(fail_reason, "Steiger;")]  # 本行目的：设置参数或赋值对象
    decision_audit[(pass_fdr == TRUE)  & (pass_nsnp == TRUE) & (pass_steiger == TRUE) & (pass_ci == FALSE), fail_reason := paste0(fail_reason, "Wide_CI;")]  # 本行目的：设置参数或赋值对象
    
    
    # A_go：全部硬门槛都过  # 本行目的：说明步骤或上下文
    decision_audit[, tier := ifelse(pass_fdr & pass_nsnp & pass_steiger & pass_ci, "A_go", "B_hold")]  # 本行目的：处理条件分支
    
    # 给 A_go 做一个排序：先 fdr，再 steiger_keep_rate，再 nsnp_harmonised，再效应  # 本行目的：说明步骤或上下文
    setorder(decision_audit, tier, fdr, -steiger_keep_rate, -nsnp_harmonised, -effect_log_or, pval)  # 本行目的：执行当前流程语句
  }  # 本行目的：组织代码块结构
  
  # 4) 输出（A_go + 审计）  # 本行目的：说明步骤或上下文
  out_audit <- file.path(decision_dir, "Table_MR_decision_audit.tsv")  # 本行目的：设置参数或赋值对象
  data.table::fwrite(decision_audit, out_audit, sep = "\t")  # 本行目的：写出表格文件
  if (verbose) message("[OK] wrote: ", out_audit, " rows=", nrow(decision_audit))  # 本行目的：进行条件判断
  
  decision_A <- decision_audit[0]  # 本行目的：设置参数或赋值对象
  if (nrow(decision_audit) > 0 && "tier" %in% names(decision_audit)) {  # 本行目的：组织代码块结构
    decision_A <- decision_audit[tier == "A_go"]  # 本行目的：设置参数或赋值对象
  }  # 本行目的：组织代码块结构
  
  out_A <- file.path(decision_dir, "Table_MR_decision_A_go.tsv")  # 本行目的：设置参数或赋值对象
  data.table::fwrite(decision_A, out_A, sep = "\t")  # 本行目的：写出表格文件
  if (verbose) message("[OK] wrote: ", out_A, " rows=", nrow(decision_A))  # 本行目的：进行条件判断
  
  # step13：森林图（可选）  # 本行目的：说明步骤或上下文
  forest <- NULL  # 本行目的：设置参数或赋值对象
  if (isTRUE(do_forest_plot)) {  # 本行目的：组织代码块结构
    forest <- make_forest_plot(  # 本行目的：设置参数或赋值对象
      key_table = if (nrow(key_table) > 0) key_table else NULL,  # 本行目的：进行条件判断
      top_table = top_table,  # 本行目的：设置参数或赋值对象
      out_dir = getwd(),  # 本行目的：设置或获取工作目录
      theme_pub = theme_pub  # 本行目的：设置参数或赋值对象
    )  # 本行目的：执行当前流程语句
    if (verbose) message("[OK] saved: ", forest$out_png)  # 本行目的：进行条件判断
    if (interactive() && !is.null(forest)) print(forest$plot)   # ✅ NEW: show in RStudio Plots  # 本行目的：进行条件判断
  }  # 本行目的：组织代码块结构
  
  # TopHits QC（可选）  # 本行目的：说明步骤或上下文
  qc <- NULL  # 本行目的：设置参数或赋值对象
  if (isTRUE(do_qc_top)) {  # 本行目的：组织代码块结构
    qc <- qc_top_table(top_table, out_dir = getwd())  # 本行目的：设置或获取工作目录
  }  # 本行目的：组织代码块结构
  # --- NEW: Forest plot (Steiger-only) ---  # 本行目的：说明步骤或上下文
  forest_steiger <- NULL  # 本行目的：设置参数或赋值对象
  if (isTRUE(do_forest_plot)) {  # 本行目的：组织代码块结构
    
    # 只保留 steiger_applied == TRUE 的结果  # 本行目的：说明步骤或上下文
    key_s <- NULL  # 本行目的：设置参数或赋值对象
    if (exists("key_table") && is.data.frame(key_table) && nrow(key_table) > 0 && "steiger_applied" %in% names(key_table)) {  # 本行目的：组织代码块结构
      key_s <- data.table::as.data.table(key_table)  # 本行目的：设置参数或赋值对象
      key_s[, steiger_applied := as.logical(steiger_applied)]  # 本行目的：设置参数或赋值对象
      key_s <- key_s[steiger_applied %in% TRUE]  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    
    top_s <- data.table::as.data.table(top_table)  # 本行目的：设置参数或赋值对象
    if ("steiger_applied" %in% names(top_s)) {  # 本行目的：组织代码块结构
      top_s[, steiger_applied := as.logical(steiger_applied)]  # 本行目的：设置参数或赋值对象
      top_s <- top_s[steiger_applied %in% TRUE]  # 本行目的：设置参数或赋值对象
    } else {  # 本行目的：组织代码块结构
      top_s <- top_s[0]  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    
    # 若过滤后完全没数据，就不画  # 本行目的：说明步骤或上下文
    if ((is.null(key_s) || nrow(key_s) == 0) && nrow(top_s) == 0) {  # 本行目的：组织代码块结构
      if (verbose) message("[INFO] No Steiger-applied rows available for Steiger-only forest plot.")  # 本行目的：进行条件判断
    } else {  # 本行目的：组织代码块结构
      out_dir_steiger <- file.path(getwd(), "steiger_only")  # 本行目的：设置或获取工作目录
      dir.create(out_dir_steiger, showWarnings = FALSE, recursive = TRUE)  # 本行目的：创建输出目录
      
      forest_steiger <- make_forest_plot(  # 本行目的：设置参数或赋值对象
        key_table = if (!is.null(key_s) && nrow(key_s) > 0) key_s else NULL,  # 本行目的：进行条件判断
        top_table = top_s,  # 本行目的：设置参数或赋值对象
        out_dir = out_dir_steiger,  # 本行目的：设置参数或赋值对象
        theme_pub = theme_pub  # 本行目的：设置参数或赋值对象
      )  # 本行目的：执行当前流程语句
      if (verbose) message("[OK] saved (Steiger-only): ", forest_steiger$out_png)  # 本行目的：进行条件判断
      if (interactive() && !is.null(forest_steiger)) print(forest_steiger$plot)  # ✅ NEW: show in Plots  # 本行目的：进行条件判断
    }  # 本行目的：组织代码块结构
  }  # 本行目的：组织代码块结构
  
  # ============================================================  # 本行目的：说明步骤或上下文
  # stepX：自动 View 最重要的东西（表格优先）  # 本行目的：说明步骤或上下文
  # 目的：跑完立刻在 RStudio 看到结果表；Windows 可自动打开目录/图片  # 本行目的：说明步骤或上下文
  # ============================================================  # 本行目的：说明步骤或上下文
  if (interactive()) {  # 本行目的：组织代码块结构
    
    # 1) 你最终选了哪些 outcome  # 本行目的：说明步骤或上下文
    top3_view <- top3[, .(id, trait, population, ss)]  # 本行目的：设置参数或赋值对象
    utils::View(as.data.frame(top3_view), title = "Selected outcomes")  # 本行目的：设置参数或赋值对象
    
    # 2) 主表（显著结果）  # 本行目的：说明步骤或上下文
    if (exists("key_table") && is.data.frame(key_table) && nrow(key_table) > 0) {  # 本行目的：组织代码块结构
      key_view_cols <- intersect(c(  # 本行目的：设置参数或赋值对象
        "CELL_TYPE", "gene_symbol", "GENE_ID",  # 本行目的：执行当前流程语句
        "outcome", "id.outcome",  # 本行目的：执行当前流程语句
        "OR_CI", "OR", "OR_lci95", "OR_uci95",  # 本行目的：执行当前流程语句
        "b", "se", "pval", "fdr", "nsnp",  # 本行目的：执行当前流程语句
        "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate", "steiger_applied"  # 本行目的：执行当前流程语句
      ), names(key_table))  # 本行目的：执行当前流程语句
      
      key_view <- data.table::as.data.table(key_table)[order(fdr, pval)][, ..key_view_cols]  # 本行目的：设置参数或赋值对象
      utils::View(as.data.frame(key_view), title = paste0("MR KEY (FDR<=", fdr_cut, ")"))  # 本行目的：设置参数或赋值对象
    } else {  # 本行目的：组织代码块结构
      message("[VIEW] key_table 为空（没有显著结果），将打开 TopHits。")  # 本行目的：输出运行信息
    }  # 本行目的：组织代码块结构
    
    # 3) TopHits（每个 outcome 前 top_n_per_outcome）  # 本行目的：说明步骤或上下文
    if (exists("top_table") && is.data.frame(top_table) && nrow(top_table) > 0) {  # 本行目的：组织代码块结构
      top_view_cols <- intersect(c(  # 本行目的：设置参数或赋值对象
        "CELL_TYPE", "gene_symbol", "GENE_ID",  # 本行目的：执行当前流程语句
        "outcome", "id.outcome",  # 本行目的：执行当前流程语句
        "OR_CI", "OR", "OR_lci95", "OR_uci95",  # 本行目的：执行当前流程语句
        "b", "se", "pval", "fdr", "nsnp",  # 本行目的：执行当前流程语句
        "nsnp_harmonised", "nsnp_steiger_keep", "steiger_keep_rate", "steiger_applied"  # 本行目的：执行当前流程语句
      ), names(top_table))  # 本行目的：执行当前流程语句
      
      top_view <- data.table::as.data.table(top_table)[order(id.outcome, pval)][, ..top_view_cols]  # 本行目的：设置参数或赋值对象
      utils::View(as.data.frame(top_view),  # 本行目的：执行当前流程语句
                  title = paste0("MR TopHits (top ", top_n_per_outcome, "/outcome)"))  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    
    # 4) 按 CELL_TYPE 的汇总（表格）  # 本行目的：说明步骤或上下文
    if (exists("mr_all") && is.data.frame(mr_all)) {  # 本行目的：组织代码块结构
      mr_dt <- data.table::as.data.table(mr_all)  # 本行目的：设置参数或赋值对象
      
      # 兜底：如果还没拆 CELL_TYPE  # 本行目的：说明步骤或上下文
      if (!"CELL_TYPE" %in% names(mr_dt) && "exposure" %in% names(mr_dt)) {  # 本行目的：组织代码块结构
        mr_dt[, CELL_TYPE := data.table::tstrsplit(exposure, "\\|")[[1]]]  # 本行目的：设置参数或赋值对象
      }  # 本行目的：组织代码块结构
      
      ct_sum <- mr_dt[is.finite(fdr), .(  # 本行目的：设置参数或赋值对象
        n_tests = .N,  # 本行目的：设置参数或赋值对象
        n_sig_fdr = sum(fdr <= fdr_cut, na.rm = TRUE),  # 本行目的：设置参数或赋值对象
        best_fdr = min(fdr, na.rm = TRUE),  # 本行目的：设置参数或赋值对象
        best_p = min(pval, na.rm = TRUE)  # 本行目的：设置参数或赋值对象
      ), by = CELL_TYPE][order(-n_sig_fdr, best_fdr)]  # 本行目的：设置参数或赋值对象
      
      utils::View(as.data.frame(ct_sum), title = paste0("Summary by CELL_TYPE (FDR<=", fdr_cut, ")"))  # 本行目的：设置参数或赋值对象
    }  # 本行目的：组织代码块结构
    
    # 5) Windows：打开输出目录和森林图  # 本行目的：说明步骤或上下文
    if (.Platform$OS.type == "windows") {  # 本行目的：组织代码块结构
      try(shell.exec(getwd()), silent = TRUE)  # 本行目的：设置或获取工作目录
      if (exists("forest") && !is.null(forest) && !is.null(forest$out_png) && file.exists(forest$out_png)) {  # 本行目的：组织代码块结构
        try(shell.exec(forest$out_png), silent = TRUE)  # 本行目的：设置参数或赋值对象
      }  # 本行目的：组织代码块结构
    }  # 本行目的：组织代码块结构
  }  # 本行目的：组织代码块结构
  
  # 返回：方便你在外部脚本里拿对象继续画图/汇报  # 本行目的：说明步骤或上下文
  invisible(list(  # 本行目的：执行当前流程语句
    base_dir = base_dir,  # 本行目的：设置参数或赋值对象
    keyword = keyword,  # 本行目的：设置参数或赋值对象
    outcomes_meta = top3,  # 本行目的：设置参数或赋值对象
    outcome_ids = top3$id,  # 本行目的：设置参数或赋值对象
    mr_all = mr_all,  # 本行目的：设置参数或赋值对象
    key_table = key_table,  # 本行目的：设置参数或赋值对象
    top_table = top_table,  # 本行目的：设置参数或赋值对象
    out_key = out_key,  # 本行目的：设置参数或赋值对象
    out_key_with_symbol = out_key2,  # 本行目的：设置参数或赋值对象
    out_top = out_top,  # 本行目的：设置参数或赋值对象
    out_top_with_symbol = out_top2,  # 本行目的：设置参数或赋值对象
    forest_by_outcome = forest_by_outcome,  # 本行目的：设置参数或赋值对象
    forest_steiger_by_outcome = forest_steiger_by_outcome,  # 本行目的：设置参数或赋值对象
    forest = forest,  # 本行目的：设置参数或赋值对象
    forest_steiger = forest_steiger,   # ✅ NEW: return steiger-only plot object  # 本行目的：设置参数或赋值对象
    qc = qc,  # 本行目的：设置参数或赋值对象
    results = results  # 本行目的：设置参数或赋值对象
  ))  # 本行目的：组织代码块结构
}  # 本行目的：组织代码块结构
## =========================================================  # 本行目的：说明步骤或上下文
## 5. 主执行区  # 本行目的：说明步骤或上下文
## =========================================================  # 本行目的：说明步骤或上下文
RESULT_DIR <- file.path(PIPELINE_ROOT, slugify(PROJECT_NAME))  # 本行目的：设置参数或赋值对象

if (dir.exists(RESULT_DIR)) {  # 本行目的：组织代码块结构
  stop(  # 本行目的：在异常条件下终止并提示
    "结果目录已存在：\n", RESULT_DIR,  # 本行目的：执行当前流程语句
    "\n为避免旧结果混入，请执行以下任一操作：",  # 本行目的：执行当前流程语句
    "\n1) 改一个新的 PROJECT_NAME",  # 本行目的：执行当前流程语句
    "\n2) 手动删除这个旧目录后再跑"  # 本行目的：执行当前流程语句
  )  # 本行目的：执行当前流程语句
}  # 本行目的：组织代码块结构

res <- run_project(  # 本行目的：调用主流程执行分析
  keyword = PROJECT_NAME,  # 本行目的：设置参数或赋值对象
  fixed_ids = FIXED_IDS,  # 本行目的：设置参数或赋值对象
  use_search_fill = FALSE,  # 本行目的：设置参数或赋值对象
  n_search = N_SEARCH,  # 本行目的：设置参数或赋值对象
  trait_regex = TRAIT_REGEX,  # 本行目的：设置参数或赋值对象
  id_regex = ID_REGEX,  # 本行目的：设置参数或赋值对象
  population_regex = POPULATION_REGEX,  # 本行目的：设置参数或赋值对象
  pipeline_root = PIPELINE_ROOT,  # 本行目的：设置参数或赋值对象
  rebuild_top3 = REBUILD_TOP3,  # 本行目的：设置参数或赋值对象
  rebuild_mr = REBUILD_MR,  # 本行目的：设置参数或赋值对象
  rebuild_ao = REBUILD_AO,  # 本行目的：设置参数或赋值对象
  do_forest_plot = DO_FOREST_PLOT,  # 本行目的：设置参数或赋值对象
  do_qc_top = DO_QC_TOP,  # 本行目的：设置参数或赋值对象
  verbose = VERBOSE_RUN  # 本行目的：设置参数或赋值对象
)  # 本行目的：执行当前流程语句

## 保存完整 res 对象  # 本行目的：说明步骤或上下文
res_rds <- file.path(res$base_dir, "project_res.rds")  # 本行目的：设置参数或赋值对象
saveRDS(res, res_rds, compress = "xz")  # 本行目的：保存RDS结果文件

cat("[OK] res 已保存：\n", res_rds, "\n\n")  # 本行目的：输出运行信息

cat("\n[DONE]\n")  # 本行目的：输出运行信息
cat("结果目录：\n", res$base_dir, "\n\n")  # 本行目的：输出运行信息

print(res$outcomes_meta[, .(id, trait, population, ss)])  # 本行目的：输出运行信息
print(res$outcome_ids)  # 本行目的：输出运行信息
