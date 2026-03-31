## =========================================================
## 目标：
## 1) 清洗并整合 WHO GDBS + GBD + SDI 数据
## 2) 生成可投稿散点图（PDF + PNG）
## 3) 生成可投稿表格（CSV + PDF）
## 说明：
## - 保持你原始统计口径不变，仅做结构整理与导出增强
## - 直接复制粘贴后，优先只改“配置区”路径即可运行
## =========================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

## =========================================================
## 1. 配置区（只需改这里）
## =========================================================
PATH_ANN6 <- "W:/Current_Work/GBD_h/WHO_GDBS_2021_extracted_clean_csv/reviewed_from_manual/annex6_reviewed_from_manual.csv"
PATH_ANN2 <- "W:/Current_Work/GBD_h/WHO_GDBS_2021_extracted_clean_csv/reviewed_from_manual/annex2_reviewed_from_manual.csv"
PATH_GBD_OUTCOME <- "W:/Current_Work/GBD_h/GBD/IHME_GBD_2023_DATA-2106ed8d-1.csv"
PATH_GBD_POP <- "W:/Current_Work/GBD_h/GBD/IHME_GBD_2023_DATA-e2b105d2-1.csv"
PATH_SDI <- "W:/Current_Work/GBD_h/GBD/IHME_GBD_SDI_2021_SDI_1950_2021_Y2024M05D16.csv"

OUT_DIR <- "./GBD_maternal_hemorrhage_output"

YEAR_MODE <- 2018
YEAR_POOL <- 2014:2018
MIN_COVERAGE <- 80
LABEL_N <- 5
LABEL_N_PER_QUADRANT <- 5

## =========================================================
## 2. 依赖包检查 + 加载
## =========================================================
required_pkgs <- c(
  "readr", "dplyr", "stringr", "countrycode",
  "ggplot2", "ggrepel", "scales", "tibble", "gridExtra"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少必要 R 包：", paste(missing_pkgs, collapse = ", "),
       "\n请先安装后再运行，例如 install.packages(c(...))")
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(countrycode)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(tibble)
  library(gridExtra)
})

## =========================================================
## 3. 通用小工具
## =========================================================
clean_names_simple <- function(x) {
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

as_iso3 <- function(x) {
  x <- as.character(x)

  custom_match <- c(
    "United States of America" = "USA",
    "United States" = "USA",
    "Russian Federation" = "RUS",
    "Russia" = "RUS",
    "Iran (Islamic Republic of)" = "IRN",
    "Iran" = "IRN",
    "Bolivia (Plurinational State of)" = "BOL",
    "Bolivia" = "BOL",
    "Venezuela (Bolivarian Republic of)" = "VEN",
    "Venezuela" = "VEN",
    "Syrian Arab Republic" = "SYR",
    "Syria" = "SYR",
    "United Republic of Tanzania" = "TZA",
    "Tanzania" = "TZA",
    "Lao People's Democratic Republic" = "LAO",
    "Laos" = "LAO",
    "Democratic People's Republic of Korea" = "PRK",
    "North Korea" = "PRK",
    "Republic of Korea" = "KOR",
    "South Korea" = "KOR",
    "Republic of Moldova" = "MDA",
    "Moldova" = "MDA",
    "Viet Nam" = "VNM",
    "Vietnam" = "VNM",
    "Türkiye" = "TUR",
    "Turkey" = "TUR",
    "Brunei Darussalam" = "BRN",
    "Brunei" = "BRN",
    "Cabo Verde" = "CPV",
    "Cape Verde" = "CPV",
    "Micronesia (Federated States of)" = "FSM",
    "The former Yugoslav Republic of Macedonia" = "MKD",
    "North Macedonia" = "MKD",
    "occupied Palestinian territory" = "PSE",
    "occupied Palestinian territory, including east Jerusalem" = "PSE",
    "State of Palestine" = "PSE",
    "Palestine" = "PSE",
    "Côte d'Ivoire" = "CIV",
    "Côte d’Ivoire" = "CIV",
    "Ivory Coast" = "CIV",
    "Eswatini" = "SWZ",
    "Swaziland" = "SWZ",
    "The Gambia" = "GMB",
    "Gambia" = "GMB",
    "The Bahamas" = "BHS",
    "Bahamas" = "BHS",
    "Democratic Republic of the Congo" = "COD",
    "Republic of the Congo" = "COG",
    "Congo" = "COG",
    "Taiwan (Province of China)" = "TWN",
    "Taiwan" = "TWN",
    "Hong Kong Special Administrative Region of China" = "HKG",
    "Hong Kong" = "HKG",
    "Macao Special Administrative Region of China" = "MAC",
    "Macao" = "MAC",
    "Czech Republic" = "CZE",
    "Czechia" = "CZE",
    "United Kingdom of Great Britain and Northern Ireland" = "GBR",
    "United Kingdom" = "GBR"
  )

  countrycode(
    sourcevar = x,
    origin = "country.name",
    destination = "iso3c",
    custom_match = custom_match,
    warn = FALSE
  )
}

read_csv_checked <- function(path) {
  if (!file.exists(path)) {
    stop("文件不存在：", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

save_plot_pdf_png <- function(p, pdf_path, png_path, width = 10, height = 7) {
  ggsave(pdf_path, p, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, p, width = width, height = height, dpi = 320)
}

save_table_pdf <- function(df, pdf_path, title_txt = NULL) {
  g <- gridExtra::tableGrob(df, rows = NULL)
  grDevices::cairo_pdf(pdf_path, width = 14, height = 8)
  if (!is.null(title_txt) && nzchar(title_txt)) {
    grid::grid.newpage()
    grid::grid.text(title_txt, x = 0.01, y = 0.98, just = c("left", "top"), gp = grid::gpar(fontsize = 13, fontface = "bold"))
    grid::grid.draw(g)
  } else {
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  grDevices::dev.off()
}

## =========================================================
## 4. 读取数据
## =========================================================
ann6 <- read_csv_checked(PATH_ANN6)
ann2 <- read_csv_checked(PATH_ANN2)
gbd_outcome <- read_csv_checked(PATH_GBD_OUTCOME)
gbd_population <- read_csv_checked(PATH_GBD_POP)
sdi <- read_csv_checked(PATH_SDI)

names(ann6) <- clean_names_simple(names(ann6))
names(ann2) <- clean_names_simple(names(ann2))

## =========================================================
## 5. 过滤与拼接（统计口径保持原样）
## =========================================================
who_base <- ann6 %>%
  transmute(
    country = as.character(country),
    year = as.integer(data_year),
    units_red_cells = as.numeric(units_red_cells),
    units_note = as.character(units_note),
    units_basis = as.character(units_basis)
  ) %>%
  left_join(
    ann2 %>%
      transmute(
        country = as.character(country),
        year = as.integer(data_year),
        coverage_pct = as.numeric(coverage_pct_value),
        coverage_status = as.character(coverage_status)
      ),
    by = c("country", "year")
  ) %>%
  filter(!is.na(units_red_cells), units_red_cells > 0) %>%
  filter(!is.na(coverage_pct), coverage_pct >= MIN_COVERAGE) %>%
  filter(year == YEAR_MODE) %>%
  mutate(iso3 = as_iso3(country)) %>%
  filter(!is.na(iso3))

pop_df <- gbd_population %>%
  transmute(
    location = as.character(location_name),
    year = as.integer(year),
    sex = as.character(sex_name),
    age = as.character(age_name),
    metric = as.character(metric_name),
    population = as.numeric(val)
  ) %>%
  filter(year %in% YEAR_POOL) %>%
  filter(str_detect(str_to_lower(sex), "both")) %>%
  filter(str_detect(str_to_lower(age), "all")) %>%
  filter(str_detect(str_to_lower(metric), "number")) %>%
  select(location, year, population) %>%
  distinct() %>%
  mutate(iso3 = as_iso3(location)) %>%
  filter(!is.na(iso3)) %>%
  distinct(iso3, year, .keep_all = TRUE)

death_df <- gbd_outcome %>%
  transmute(
    location = as.character(location_name),
    year = as.integer(year),
    sex = as.character(sex_name),
    age = as.character(age_name),
    measure = as.character(measure_name),
    metric = as.character(metric_name),
    cause = as.character(cause_name),
    deaths = as.numeric(val)
  ) %>%
  filter(year %in% YEAR_POOL) %>%
  filter(str_detect(str_to_lower(sex), "female")) %>%
  filter(str_detect(str_to_lower(age), "all")) %>%
  filter(str_detect(str_to_lower(measure), "death")) %>%
  filter(str_detect(str_to_lower(metric), "number")) %>%
  filter(str_detect(str_to_lower(cause), "maternal hemorrhage")) %>%
  select(location, year, deaths) %>%
  distinct() %>%
  mutate(iso3 = as_iso3(location)) %>%
  filter(!is.na(iso3)) %>%
  distinct(iso3, year, .keep_all = TRUE)

inc_df <- gbd_outcome %>%
  transmute(
    location = as.character(location_name),
    year = as.integer(year),
    sex = as.character(sex_name),
    age = as.character(age_name),
    measure = as.character(measure_name),
    metric = as.character(metric_name),
    cause = as.character(cause_name),
    incidence = as.numeric(val)
  ) %>%
  filter(year %in% YEAR_POOL) %>%
  filter(str_detect(str_to_lower(sex), "female")) %>%
  filter(str_detect(str_to_lower(age), "all")) %>%
  filter(str_detect(str_to_lower(measure), "incidence")) %>%
  filter(str_detect(str_to_lower(metric), "number")) %>%
  filter(str_detect(str_to_lower(cause), "maternal hemorrhage")) %>%
  select(location, year, incidence) %>%
  distinct() %>%
  mutate(iso3 = as_iso3(location)) %>%
  filter(!is.na(iso3)) %>%
  distinct(iso3, year, .keep_all = TRUE)

sdi_df <- sdi %>%
  transmute(
    location = as.character(location_name),
    year = as.integer(year_id),
    sdi = as.numeric(mean_value)
  ) %>%
  filter(year %in% YEAR_POOL) %>%
  distinct() %>%
  mutate(iso3 = as_iso3(location)) %>%
  filter(!is.na(iso3)) %>%
  distinct(iso3, year, .keep_all = TRUE)

plot_df0 <- who_base %>%
  left_join(pop_df %>% select(iso3, year, population), by = c("iso3", "year")) %>%
  left_join(death_df %>% select(iso3, year, deaths), by = c("iso3", "year")) %>%
  left_join(inc_df %>% select(iso3, year, incidence), by = c("iso3", "year")) %>%
  left_join(sdi_df %>% select(iso3, year, sdi), by = c("iso3", "year"))

diagnostic_tbl <- tibble::tibble(
  step = c(
    "WHO after filter",
    "joined population",
    "joined deaths",
    "joined incidence",
    "joined SDI",
    "final plot data"
  ),
  n_country = c(
    nrow(who_base),
    sum(!is.na(plot_df0$population)),
    sum(!is.na(plot_df0$deaths)),
    sum(!is.na(plot_df0$incidence)),
    sum(!is.na(plot_df0$sdi)),
    sum(!is.na(plot_df0$population) & !is.na(plot_df0$deaths) & !is.na(plot_df0$incidence) & !is.na(plot_df0$sdi))
  )
)

plot_df <- plot_df0 %>%
  mutate(
    rbc_per_1000 = units_red_cells / population * 1000,
    dir_per_1000 = deaths / incidence * 1000
  ) %>%
  filter(!is.na(rbc_per_1000), rbc_per_1000 > 0) %>%
  filter(!is.na(dir_per_1000), dir_per_1000 > 0) %>%
  filter(!is.na(incidence), incidence > 0) %>%
  filter(!is.na(sdi)) %>%
  mutate(
    mismatch_score =
      as.numeric(scale(log10(dir_per_1000))) -
      as.numeric(scale(log10(rbc_per_1000)))
  ) %>%
  arrange(desc(mismatch_score))

if (nrow(plot_df) == 0) {
  stop("plot_df 为空：请检查年份、覆盖率阈值、或原始文件内容。")
}

## =========================================================
## 6. 图1：全局 mismatch topN 标注
## =========================================================
label_df_1 <- plot_df %>%
  slice_head(n = min(LABEL_N, nrow(plot_df)))

x_med <- median(plot_df$rbc_per_1000, na.rm = TRUE)
y_med <- median(plot_df$dir_per_1000, na.rm = TRUE)
subtitle_txt_1 <- paste0("Year = ", YEAR_MODE, ", coverage ≥ ", MIN_COVERAGE, "%; N = ", nrow(plot_df))

p1 <- ggplot(plot_df, aes(x = rbc_per_1000, y = dir_per_1000)) +
  geom_vline(xintercept = x_med, linetype = 2, linewidth = 0.5, colour = "grey75") +
  geom_hline(yintercept = y_med, linetype = 2, linewidth = 0.5, colour = "grey75") +
  geom_point(aes(size = incidence, colour = sdi), alpha = 0.85) +
  geom_text_repel(
    data = label_df_1,
    aes(label = country),
    size = 3.6,
    seed = 123,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.alpha = 0.5,
    max.overlaps = Inf
  ) +
  scale_x_log10(labels = label_number(accuracy = 0.1)) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  scale_size_continuous(range = c(2.8, 10), labels = label_comma()) +
  scale_colour_viridis_c(option = "C", end = 0.95, labels = label_number(accuracy = 0.01)) +
  labs(
    title = "Maternal hemorrhage fatality vs blood access",
    subtitle = subtitle_txt_1,
    x = "Red cell units per 1,000 population (WHO GDBS)",
    y = "Deaths per 1,000 incident maternal hemorrhage cases (GBD)",
    size = "Incidence\n(number)",
    colour = "SDI",
    caption = "Labels = highest mismatch score (high fatality + low blood access)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey90"),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, colour = "grey30"),
    plot.caption = element_text(size = 10, colour = "grey40"),
    legend.position = "right"
  )

## =========================================================
## 7. 图2：四象限极端国家标注
## =========================================================
plot_df_label <- plot_df %>%
  mutate(
    quadrant = case_when(
      rbc_per_1000 < x_med & dir_per_1000 >= y_med ~ "upper_left",
      rbc_per_1000 >= x_med & dir_per_1000 >= y_med ~ "upper_right",
      rbc_per_1000 < x_med & dir_per_1000 < y_med ~ "lower_left",
      TRUE ~ "lower_right"
    ),
    x_z = as.numeric(scale(log10(rbc_per_1000))),
    y_z = as.numeric(scale(log10(dir_per_1000))),
    quad_score = case_when(
      quadrant == "upper_left"  ~ (-x_z) + y_z,
      quadrant == "upper_right" ~ x_z + y_z,
      quadrant == "lower_left"  ~ (-x_z) + (-y_z),
      quadrant == "lower_right" ~ x_z + (-y_z)
    )
  )

label_df_2 <- plot_df_label %>%
  group_by(quadrant) %>%
  arrange(desc(quad_score), .by_group = TRUE) %>%
  slice_head(n = LABEL_N_PER_QUADRANT) %>%
  ungroup()

subtitle_txt_2 <- paste0("Year = ", YEAR_MODE, ", coverage ≥ ", MIN_COVERAGE, "%; N = ", nrow(plot_df_label))

p2 <- ggplot(plot_df_label, aes(x = rbc_per_1000, y = dir_per_1000)) +
  geom_vline(xintercept = x_med, linetype = 2, linewidth = 0.5, colour = "grey75") +
  geom_hline(yintercept = y_med, linetype = 2, linewidth = 0.5, colour = "grey75") +
  geom_point(aes(size = incidence, colour = sdi), alpha = 0.85) +
  geom_text_repel(
    data = label_df_2,
    aes(label = country),
    size = 3.4,
    seed = 123,
    box.padding = 0.45,
    point.padding = 0.3,
    segment.alpha = 0.5,
    max.overlaps = Inf
  ) +
  scale_x_log10(labels = label_number(accuracy = 0.1)) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  scale_size_continuous(range = c(2.8, 10), labels = label_comma()) +
  scale_colour_viridis_c(option = "C", end = 0.95, labels = label_number(accuracy = 0.01)) +
  labs(
    title = "Maternal hemorrhage fatality vs blood access",
    subtitle = subtitle_txt_2,
    x = "Red cell units per 1,000 population (WHO GDBS)",
    y = "Deaths per 1,000 incident maternal hemorrhage cases (GBD)",
    size = "Incidence\n(number)",
    colour = "SDI",
    caption = "Labels = top 5 most extreme countries in each quadrant"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey90"),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, colour = "grey30"),
    plot.caption = element_text(size = 10, colour = "grey40"),
    legend.position = "right"
  )

## =========================================================
## 8. 生成表格（四象限 Top 表）
## =========================================================
table1_df <- label_df_2 %>%
  mutate(
    quadrant = factor(
      quadrant,
      levels = c("upper_left", "upper_right", "lower_left", "lower_right"),
      labels = c("Upper left", "Upper right", "Lower left", "Lower right")
    )
  ) %>%
  group_by(quadrant) %>%
  arrange(desc(quad_score), .by_group = TRUE) %>%
  mutate(rank_in_quadrant = row_number()) %>%
  ungroup() %>%
  transmute(
    Quadrant = quadrant,
    Rank = rank_in_quadrant,
    Country = country,
    `Red cell units per 1,000 population` = round(rbc_per_1000, 2),
    `Deaths per 1,000 incident cases` = round(dir_per_1000, 2),
    `Incidence (number)` = round(incidence, 1),
    `Deaths (number)` = round(deaths, 1),
    SDI = round(sdi, 3)
  )

## =========================================================
## 9. 导出文件（投稿友好）
## =========================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

write_csv(plot_df0, file.path(OUT_DIR, "data_joined_raw.csv"))
write_csv(plot_df, file.path(OUT_DIR, "data_analysis_ready.csv"))
write_csv(diagnostic_tbl, file.path(OUT_DIR, "table_diagnostic_counts.csv"))
write_csv(table1_df, file.path(OUT_DIR, "table_quadrant_top_countries.csv"))

save_plot_pdf_png(
  p1,
  pdf_path = file.path(OUT_DIR, "figure1_mismatch_topN.pdf"),
  png_path = file.path(OUT_DIR, "figure1_mismatch_topN.png"),
  width = 10,
  height = 7
)

save_plot_pdf_png(
  p2,
  pdf_path = file.path(OUT_DIR, "figure2_quadrant_topN.pdf"),
  png_path = file.path(OUT_DIR, "figure2_quadrant_topN.png"),
  width = 10,
  height = 7
)

save_table_pdf(
  diagnostic_tbl,
  pdf_path = file.path(OUT_DIR, "table_diagnostic_counts.pdf"),
  title_txt = "Diagnostic country counts by merge step"
)

save_table_pdf(
  table1_df,
  pdf_path = file.path(OUT_DIR, "table_quadrant_top_countries.pdf"),
  title_txt = "Top extreme countries in each quadrant"
)

message("完成：结果已导出到 -> ", normalizePath(OUT_DIR, mustWork = FALSE))

