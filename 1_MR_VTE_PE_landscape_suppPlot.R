## =========================================================
## 0. 清空环境
## =========================================================
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)


## =========================================================
## 1. 配置区
## 以后主要改这里
## =========================================================
PROJECT_DIR <- "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE/Venous_thromboembolism_PE_VTE_only"

RES_RDS <- file.path(PROJECT_DIR, "project_res.rds")

THEME_FILE <- "W:/fangxiaobin-CTRH-omics/000temple/scico_plot_theme.R"

VPUB_FILE <- "W:/fangxiaobin-CTRH-omics/OneK1K_PIPELINE/literature/V_publish.xlsx"
VPUB_SHEET <- "S3"

FDR_CUT <- 0.05
LABEL_TOP_N <- 10L

BLOCK_WIDTH_EQ <- 100
BLOCK_GAP_EQ <- 12

POINT_SIZE_ALL <- 1.0
POINT_SIZE_SIG <- 1.9

COLOR_PALETTE <- "vik"

SAVE_PNG <- TRUE
SAVE_PDF <- TRUE


## =========================================================
## 2. 加载程序包
## =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
  library(scico)
})

if (!file.exists(RES_RDS)) {
  stop("找不到 res 文件：", RES_RDS)
}

res <- readRDS(RES_RDS)

if (file.exists(THEME_FILE)) {
  source(THEME_FILE)
}

if (!exists("theme_pub")) {
  theme_pub <- function(base = 11) ggplot2::theme_minimal(base_size = base)
}


## =========================================================
## 3. 检查 res
## =========================================================
if (is.null(res)) {
  stop("res 为空。")
}

if (is.null(res$mr_all) || nrow(res$mr_all) == 0) {
  stop("res$mr_all 为空。")
}

if (is.null(res$outcome_ids) || length(res$outcome_ids) == 0) {
  stop("res$outcome_ids 为空。")
}

if (is.null(res$base_dir) || !nzchar(res$base_dir)) {
  stop("res$base_dir 为空。")
}


## =========================================================
## 4. 准备主表 manh_dt
## =========================================================
manh_dt <- as.data.table(copy(res$mr_all))

if (!all(c("CELL_TYPE", "GENE_ID") %in% names(manh_dt))) {
  if ("exposure" %in% names(manh_dt)) {
    manh_dt[, c("CELL_TYPE", "GENE_ID") := tstrsplit(exposure, "\\|")]
  }
}

if (!all(c("CELL_TYPE", "GENE_ID") %in% names(manh_dt))) {
  stop("mr_all 中缺少 CELL_TYPE / GENE_ID，且 exposure 也无法拆分。")
}

if (!"id.outcome" %in% names(manh_dt)) {
  stop("mr_all 中缺少 id.outcome。")
}

if (!"pval" %in% names(manh_dt)) {
  stop("mr_all 中缺少 pval。")
}

if (!"fdr" %in% names(manh_dt)) {
  manh_dt[, fdr := NA_real_]
}

manh_dt[, GENE_ID := sub("\\.\\d+$", "", trimws(GENE_ID))]

num_cols <- intersect(c("pval", "fdr"), names(manh_dt))
manh_dt[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]

manh_dt <- manh_dt[
  is.finite(pval) & pval > 0 &
    is.finite(fdr) & fdr > 0
]

if (nrow(manh_dt) == 0) {
  stop("过滤后 manh_dt 为空。")
}


## =========================================================
## 5. 补 gene_symbol
## =========================================================
s3 <- as.data.table(readxl::read_xlsx(VPUB_FILE, sheet = VPUB_SHEET, skip = 1))
if (nrow(s3) >= 1) {
  s3 <- s3[-.N]
}

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

manh_dt <- merge(manh_dt, gene_map, by = "GENE_ID", all.x = TRUE)

if (!"gene_symbol" %in% names(manh_dt)) {
  manh_dt[, gene_symbol := NA_character_]
}

manh_dt[, gene_show := fifelse(
  !is.na(gene_symbol) & nzchar(gene_symbol),
  gene_symbol,
  GENE_ID
)]

if (!"outcome" %in% names(manh_dt)) {
  manh_dt[, outcome := NA_character_]
}

manh_dt[, outcome_show := fifelse(
  !is.na(outcome) & nzchar(outcome),
  outcome,
  id.outcome
)]

manh_dt[, neg_log10_fdr := -log10(pmax(fdr, 1e-300))]
manh_dt[, sig_fdr := is.finite(fdr) & fdr <= FDR_CUT]

setorder(manh_dt, id.outcome, CELL_TYPE, gene_show, fdr, pval)

out_all_tsv <- file.path(res$base_dir, "Fig_MR_FDR_landscape_all_data.tsv")
fwrite(manh_dt, out_all_tsv, sep = "\t")


## =========================================================
## 6. 逐个 outcome 作图
## 一个循环里同时输出：
##   1) 等宽横版
##   2) 等宽翻转版
## =========================================================
for (oid in as.character(res$outcome_ids)) {
  
  dt_i <- copy(manh_dt[id.outcome %in% oid])
  
  if (nrow(dt_i) == 0) {
    next
  }
  
  setorder(dt_i, CELL_TYPE, gene_show, fdr, pval)
  
  ## -------------------------
  ## 每个 CELL_TYPE 固定等宽布局
  ## -------------------------
  ct_levels_i <- unique(dt_i$CELL_TYPE)
  
  layout_i <- data.table(
    CELL_TYPE = ct_levels_i,
    ct_idx = seq_along(ct_levels_i)
  )
  
  layout_i[, block_start := (ct_idx - 1) * (BLOCK_WIDTH_EQ + BLOCK_GAP_EQ) + 1]
  layout_i[, block_end := block_start + BLOCK_WIDTH_EQ - 1]
  layout_i[, x_center := (block_start + block_end) / 2]
  layout_i[, y_center := x_center]
  
  dt_i <- merge(dt_i, layout_i, by = "CELL_TYPE", all.x = TRUE, sort = FALSE)
  
  dt_i[, x := {
    n_i <- .N
    s_i <- block_start[1]
    e_i <- block_end[1]
    if (n_i == 1) (s_i + e_i) / 2 else seq(s_i, e_i, length.out = n_i)
  }, by = CELL_TYPE]
  
  dt_i[, y_pos := {
    n_i <- .N
    s_i <- block_start[1]
    e_i <- block_end[1]
    if (n_i == 1) (s_i + e_i) / 2 else seq(s_i, e_i, length.out = n_i)
  }, by = CELL_TYPE]
  
  sep_pos_i <- layout_i$block_end + BLOCK_GAP_EQ / 2
  if (length(sep_pos_i) > 0) {
    sep_pos_i <- sep_pos_i[-length(sep_pos_i)]
  }
  
  p_line_i <- -log10(FDR_CUT)
  
  label_dt_i <- dt_i[sig_fdr == TRUE][order(fdr, pval)]
  if (nrow(label_dt_i) > LABEL_TOP_N) {
    label_dt_i <- label_dt_i[seq_len(LABEL_TOP_N)]
  }
  if (nrow(label_dt_i) == 0) {
    label_dt_i <- dt_i[order(fdr, pval)][seq_len(min(nrow(dt_i), LABEL_TOP_N))]
  }
  
  title_i <- unique(dt_i$outcome_show)
  title_i <- title_i[!is.na(title_i) & nzchar(title_i)]
  if (length(title_i) == 0) {
    title_i <- oid
  } else {
    title_i <- title_i[1]
  }
  
  out_dir_i <- file.path(res$base_dir, paste0("outcome_", oid))
  dir.create(out_dir_i, showWarnings = FALSE, recursive = TRUE)
  
  
  ## =====================================================
  ## 6A. 等宽横版
  ## =====================================================
  p_i_eq <- ggplot(
    dt_i,
    aes(x = x, y = neg_log10_fdr)
  ) +
    geom_point(
      aes(color = CELL_TYPE),
      size = POINT_SIZE_ALL,
      alpha = 0.80,
      show.legend = FALSE
    ) +
    geom_point(
      data = dt_i[sig_fdr == TRUE],
      color = "red",
      size = POINT_SIZE_SIG,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    scico::scale_color_scico_d(
      palette = COLOR_PALETTE
    ) +
    scale_x_continuous(
      breaks = layout_i$x_center,
      labels = layout_i$CELL_TYPE,
      expand = c(0.01, 0.01)
    ) +
    labs(
      x = NULL,
      y = expression(-log[10](FDR)),
      title = paste0("MR FDR landscape plot (equal width): ", title_i),
      subtitle = paste0(
        "Each CELL_TYPE occupies the same width; dashed line = FDR = ",
        FDR_CUT
      )
    ) +
    theme_pub(base = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 7
      ),
      plot.margin = margin(5.5, 20, 5.5, 5.5)
    ) +
    coord_cartesian(clip = "off") +
    geom_hline(
      yintercept = p_line_i,
      linetype = "dashed",
      color = "red",
      linewidth = 0.4
    )
  
  if (length(sep_pos_i) > 0) {
    p_i_eq <- p_i_eq +
      geom_vline(
        xintercept = sep_pos_i,
        linetype = "solid",
        color = "grey85",
        linewidth = 0.3
      )
  }
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_i_eq <- p_i_eq +
      ggrepel::geom_text_repel(
        data = label_dt_i,
        aes(label = gene_show),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.20,
        min.segment.length = 0,
        show.legend = FALSE
      )
  } else {
    p_i_eq <- p_i_eq +
      geom_text(
        data = label_dt_i,
        aes(label = gene_show),
        size = 2.6,
        vjust = -0.5,
        check_overlap = TRUE,
        show.legend = FALSE
      )
  }
  
  print(p_i_eq)
  
  out_png_i_eq <- file.path(out_dir_i, "Fig_MR_FDR_landscape_equal_width.png")
  out_pdf_i_eq <- file.path(out_dir_i, "Fig_MR_FDR_landscape_equal_width.pdf")
  out_tsv_i_eq <- file.path(out_dir_i, "Fig_MR_FDR_landscape_equal_width_data.tsv")
  
  if (SAVE_PNG) {
    ggsave(
      filename = out_png_i_eq,
      plot = p_i_eq,
      width = 14,
      height = 6.5,
      dpi = 300
    )
  }
  
  if (SAVE_PDF) {
    pdf(out_pdf_i_eq, width = 14, height = 6.5)
    print(p_i_eq)
    dev.off()
  }
  
  fwrite(dt_i, out_tsv_i_eq, sep = "\t")
  
  cat("[OK] 已输出等宽横版：\n")
  if (SAVE_PNG) cat(out_png_i_eq, "\n")
  if (SAVE_PDF) cat(out_pdf_i_eq, "\n")
  cat(out_tsv_i_eq, "\n\n")
  
  
  ## =====================================================
  ## 6B. 等宽翻转版
  ## =====================================================
  p_i_flip <- ggplot(
    dt_i,
    aes(x = neg_log10_fdr, y = y_pos)
  ) +
    geom_point(
      aes(color = CELL_TYPE),
      size = POINT_SIZE_ALL,
      alpha = 0.80,
      show.legend = FALSE
    ) +
    geom_point(
      data = dt_i[sig_fdr == TRUE],
      color = "red",
      size = POINT_SIZE_SIG,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    scico::scale_color_scico_d(
      palette = COLOR_PALETTE
    ) +
    scale_y_continuous(
      breaks = layout_i$y_center,
      labels = layout_i$CELL_TYPE,
      expand = c(0.01, 0.01)
    ) +
    labs(
      x = expression(-log[10](FDR)),
      y = NULL,
      title = paste0("MR FDR landscape plot (equal width, flipped): ", title_i),
      subtitle = paste0(
        "Horizontal axis = -log10(FDR); each CELL_TYPE occupies the same vertical space; line = FDR = ",
        FDR_CUT
      )
    ) +
    theme_pub(base = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 7),
      plot.margin = margin(5.5, 50, 5.5, 5.5)
    ) +
    coord_cartesian(clip = "off") +
    geom_vline(
      xintercept = p_line_i,
      linetype = "dashed",
      color = "red",
      linewidth = 0.4
    )
  
  if (length(sep_pos_i) > 0) {
    p_i_flip <- p_i_flip +
      geom_hline(
        yintercept = sep_pos_i,
        linetype = "solid",
        color = "grey85",
        linewidth = 0.3
      )
  }
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_i_flip <- p_i_flip +
      ggrepel::geom_text_repel(
        data = label_dt_i,
        aes(label = gene_show),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.20,
        min.segment.length = 0,
        show.legend = FALSE
      )
  } else {
    p_i_flip <- p_i_flip +
      geom_text(
        data = label_dt_i,
        aes(label = gene_show),
        size = 2.6,
        hjust = -0.1,
        check_overlap = TRUE,
        show.legend = FALSE
      )
  }
  
  print(p_i_flip)
  
  out_png_i_flip <- file.path(out_dir_i, "Fig_MR_FDR_landscape_equal_width_flipped.png")
  out_pdf_i_flip <- file.path(out_dir_i, "Fig_MR_FDR_landscape_equal_width_flipped.pdf")
  out_tsv_i_flip <- file.path(out_dir_i, "Fig_MR_FDR_landscape_equal_width_flipped_data.tsv")
  
  if (SAVE_PNG) {
    ggsave(
      filename = out_png_i_flip,
      plot = p_i_flip,
      width = 10,
      height = 12,
      dpi = 300
    )
  }
  
  if (SAVE_PDF) {
    pdf(out_pdf_i_flip, width = 10, height = 12)
    print(p_i_flip)
    dev.off()
  }
  
  fwrite(dt_i, out_tsv_i_flip, sep = "\t")
  
  cat("[OK] 已输出等宽翻转版：\n")
  if (SAVE_PNG) cat(out_png_i_flip, "\n")
  if (SAVE_PDF) cat(out_pdf_i_flip, "\n")
  cat(out_tsv_i_flip, "\n\n")
}

cat("[DONE] 全部完成。\n")
cat("总底表：\n", out_all_tsv, "\n")
