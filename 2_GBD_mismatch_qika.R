## =========================================================
## GBD × WHO 输血可及性失配分析（高质量出图版）
## 目标：
## 1) 输出目录默认与分析脚本同级
## 2) 生成可投稿风格散点图、表格 PDF、图表合并 PDF、Word 附件
## =========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(countrycode)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(gridExtra)
  library(grid)
  library(flextable)
  library(officer)
})

## =========================================================
## 1. 路径与参数
## =========================================================
YEAR_MODE <- 2018
MIN_COVERAGE <- 80
LABEL_N_PER_QUADRANT <- 5

ANN6_FILE <- "W:/Current_Work/GBD_h/WHO_GDBS_2021_extracted_clean_csv/reviewed_from_manual/annex6_reviewed_from_manual.csv"
ANN2_FILE <- "W:/Current_Work/GBD_h/WHO_GDBS_2021_extracted_clean_csv/reviewed_from_manual/annex2_reviewed_from_manual.csv"
GBD_OUTCOME_FILE <- "W:/Current_Work/GBD_h/GBD/IHME_GBD_2023_DATA-2106ed8d-1.csv"
GBD_POP_FILE <- "W:/Current_Work/GBD_h/GBD/IHME_GBD_2023_DATA-e2b105d2-1.csv"
SDI_FILE <- "W:/Current_Work/GBD_h/GBD/IHME_GBD_SDI_2021_SDI_1950_2021_Y2024M05D16.csv"

# 输出目录：尽量和分析脚本在一起
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    normalizePath(dirname(sub("^--file=", "", hit[1])), winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}

SCRIPT_DIR <- get_script_dir()
OUT_DIR <- file.path(SCRIPT_DIR, paste0("outputs_mismatch_", YEAR_MODE))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 2. 工具函数
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
    "United States of America" = "USA", "United States" = "USA",
    "Russian Federation" = "RUS", "Russia" = "RUS",
    "Iran (Islamic Republic of)" = "IRN", "Iran" = "IRN",
    "Bolivia (Plurinational State of)" = "BOL", "Bolivia" = "BOL",
    "Venezuela (Bolivarian Republic of)" = "VEN", "Venezuela" = "VEN",
    "Syrian Arab Republic" = "SYR", "Syria" = "SYR",
    "United Republic of Tanzania" = "TZA", "Tanzania" = "TZA",
    "Lao People's Democratic Republic" = "LAO", "Laos" = "LAO",
    "Democratic People's Republic of Korea" = "PRK", "North Korea" = "PRK",
    "Republic of Korea" = "KOR", "South Korea" = "KOR",
    "Republic of Moldova" = "MDA", "Moldova" = "MDA",
    "Viet Nam" = "VNM", "Vietnam" = "VNM",
    "Türkiye" = "TUR", "Turkey" = "TUR",
    "Brunei Darussalam" = "BRN", "Brunei" = "BRN",
    "Cabo Verde" = "CPV", "Cape Verde" = "CPV",
    "Micronesia (Federated States of)" = "FSM",
    "The former Yugoslav Republic of Macedonia" = "MKD", "North Macedonia" = "MKD",
    "occupied Palestinian territory" = "PSE",
    "occupied Palestinian territory, including east Jerusalem" = "PSE",
    "State of Palestine" = "PSE", "Palestine" = "PSE",
    "Côte d'Ivoire" = "CIV", "Côte d’Ivoire" = "CIV", "Ivory Coast" = "CIV",
    "Eswatini" = "SWZ", "Swaziland" = "SWZ",
    "The Gambia" = "GMB", "Gambia" = "GMB",
    "The Bahamas" = "BHS", "Bahamas" = "BHS",
    "Democratic Republic of the Congo" = "COD",
    "Republic of the Congo" = "COG", "Congo" = "COG",
    "Taiwan (Province of China)" = "TWN", "Taiwan" = "TWN",
    "Hong Kong Special Administrative Region of China" = "HKG", "Hong Kong" = "HKG",
    "Macao Special Administrative Region of China" = "MAC", "Macao" = "MAC",
    "Czech Republic" = "CZE", "Czechia" = "CZE",
    "United Kingdom of Great Britain and Northern Ireland" = "GBR", "United Kingdom" = "GBR"
  )
  countrycode(x, origin = "country.name", destination = "iso3c", custom_match = custom_match, warn = FALSE)
}

## =========================================================
## 3. 读取数据
## =========================================================
for (f in c(ANN6_FILE, ANN2_FILE, GBD_OUTCOME_FILE, GBD_POP_FILE, SDI_FILE)) {
  if (!file.exists(f)) stop("找不到输入文件：", f)
}

ann6 <- read_csv(ANN6_FILE, show_col_types = FALSE)
ann2 <- read_csv(ANN2_FILE, show_col_types = FALSE)
gbd_outcome <- read_csv(GBD_OUTCOME_FILE, show_col_types = FALSE)
gbd_population <- read_csv(GBD_POP_FILE, show_col_types = FALSE)
sdi <- read_csv(SDI_FILE, show_col_types = FALSE)

names(ann6) <- clean_names_simple(names(ann6))
names(ann2) <- clean_names_simple(names(ann2))

## =========================================================
## 4. 构建分析数据
## =========================================================
who_base <- ann6 %>%
  transmute(
    country = as.character(country),
    year = as.integer(data_year),
    units_red_cells = as.numeric(units_red_cells)
  ) %>%
  left_join(
    ann2 %>%
      transmute(
        country = as.character(country),
        year = as.integer(data_year),
        coverage_pct = as.numeric(coverage_pct_value)
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
  filter(year %in% 2014:2018) %>%
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
  filter(year %in% 2014:2018) %>%
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
  filter(year %in% 2014:2018) %>%
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
  filter(year %in% 2014:2018) %>%
  distinct() %>%
  mutate(iso3 = as_iso3(location)) %>%
  filter(!is.na(iso3)) %>%
  distinct(iso3, year, .keep_all = TRUE)

plot_df0 <- who_base %>%
  left_join(pop_df %>% select(iso3, year, population), by = c("iso3", "year")) %>%
  left_join(death_df %>% select(iso3, year, deaths), by = c("iso3", "year")) %>%
  left_join(inc_df %>% select(iso3, year, incidence), by = c("iso3", "year")) %>%
  left_join(sdi_df %>% select(iso3, year, sdi), by = c("iso3", "year"))

plot_df <- plot_df0 %>%
  mutate(
    rbc_per_1000 = units_red_cells / population * 1000,
    dir_per_1000 = deaths / incidence * 1000
  ) %>%
  filter(!is.na(rbc_per_1000), rbc_per_1000 > 0) %>%
  filter(!is.na(dir_per_1000), dir_per_1000 > 0) %>%
  filter(!is.na(incidence), incidence > 0) %>%
  filter(!is.na(sdi))

if (nrow(plot_df) == 0) stop("过滤后 plot_df 为空，请检查输入数据。")

## =========================================================
## 5. 四象限与标签
## =========================================================
x_med <- median(plot_df$rbc_per_1000, na.rm = TRUE)
y_med <- median(plot_df$dir_per_1000, na.rm = TRUE)

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

label_df <- plot_df_label %>%
  group_by(quadrant) %>%
  arrange(desc(quad_score), .by_group = TRUE) %>%
  slice_head(n = LABEL_N_PER_QUADRANT) %>%
  ungroup()

subtitle_txt <- paste0(
  "Year = ", YEAR_MODE,
  ", coverage ≥ ", MIN_COVERAGE, "%; N = ", nrow(plot_df_label)
)

## =========================================================
## 6. 主图（投稿风格）
## =========================================================
p <- ggplot(plot_df_label, aes(x = rbc_per_1000, y = dir_per_1000)) +
  geom_vline(xintercept = x_med, linetype = 2, linewidth = 0.55, colour = "grey70") +
  geom_hline(yintercept = y_med, linetype = 2, linewidth = 0.55, colour = "grey70") +
  geom_point(aes(size = incidence, colour = sdi), alpha = 0.88) +
  geom_text_repel(
    data = label_df,
    aes(label = country),
    size = 3.6,
    seed = 123,
    box.padding = 0.45,
    point.padding = 0.32,
    segment.alpha = 0.55,
    max.overlaps = Inf
  ) +
  scale_x_log10(labels = label_number(accuracy = 0.1)) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  scale_size_continuous(range = c(2.8, 10), labels = label_comma()) +
  scale_colour_viridis_c(option = "C", end = 0.95, labels = label_number(accuracy = 0.01)) +
  labs(
    title = "Maternal hemorrhage fatality vs blood access",
    subtitle = subtitle_txt,
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
    plot.subtitle = element_text(size = 11.2, colour = "grey30"),
    plot.caption = element_text(size = 10, colour = "grey40"),
    legend.position = "right",
    plot.margin = margin(12, 16, 10, 12)
  )

## =========================================================
## 7. 表格（投稿风格）
## =========================================================
phenotype_map <- c(
  upper_left  = "Low RBC / High fatality",
  upper_right = "High RBC / High fatality",
  lower_left  = "Low RBC / Low fatality",
  lower_right = "High RBC / Low fatality"
)

phenotype_levels <- c(
  "Low RBC / High fatality",
  "High RBC / High fatality",
  "Low RBC / Low fatality",
  "High RBC / Low fatality"
)

table1_df <- label_df %>%
  mutate(
    Phenotype = recode(quadrant, !!!phenotype_map),
    Phenotype = factor(Phenotype, levels = phenotype_levels)
  ) %>%
  group_by(Phenotype) %>%
  arrange(desc(quad_score), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup() %>%
  transmute(
    Phenotype,
    Rank,
    Country = country,
    `RBC/1,000` = rbc_per_1000,
    `Deaths/1,000 cases` = dir_per_1000,
    Incidence = incidence,
    Deaths = deaths,
    SDI = sdi
  )

table1_print <- table1_df %>%
  mutate(
    Phenotype = case_when(
      Phenotype == "Low RBC / High fatality"  ~ "Low RBC /\nHigh fatality",
      Phenotype == "High RBC / High fatality" ~ "High RBC /\nHigh fatality",
      Phenotype == "Low RBC / Low fatality"   ~ "Low RBC /\nLow fatality",
      TRUE                                      ~ "High RBC /\nLow fatality"
    ),
    Country = str_wrap(Country, width = 13),
    `RBC/1,000` = sprintf("%.2f", `RBC/1,000`),
    `Deaths/1,000 cases` = sprintf("%.2f", `Deaths/1,000 cases`),
    Incidence = comma(round(Incidence, 0)),
    Deaths = format(round(Deaths, 1), nsmall = 1, trim = TRUE, big.mark = ","),
    SDI = sprintf("%.3f", SDI)
  )

core_fill <- rep(c("white", "#F7F9FC"), length.out = nrow(table1_print))

table_theme <- ttheme_minimal(
  base_size = 7.4,
  core = list(
    fg_params = list(
      col = "#222222",
      fontsize = 7.4,
      hjust = c(0, 0.5, 0, 1, 1, 1, 1, 1),
      x = c(0.02, 0.50, 0.02, 0.98, 0.98, 0.98, 0.98, 0.98)
    ),
    bg_params = list(fill = core_fill, col = NA)
  ),
  colhead = list(
    fg_params = list(col = "white", fontsize = 7.6, fontface = "bold"),
    bg_params = list(fill = "#2F4858", col = NA)
  )
)

tbl_grob <- tableGrob(table1_print, rows = NULL, theme = table_theme)
tbl_grob$widths <- unit(c(0.23, 0.05, 0.16, 0.10, 0.14, 0.10, 0.09, 0.07), "npc")

table_title <- textGrob(
  paste0("Table 1. Top ", LABEL_N_PER_QUADRANT, " countries within each mismatch phenotype in ", YEAR_MODE),
  x = 0, hjust = 0,
  gp = gpar(fontsize = 11, fontface = "bold", col = "#1F1F1F")
)

table_note <- textGrob(
  "Phenotypes were defined by sample medians of RBC units per 1,000 population and maternal hemorrhage deaths per 1,000 incident cases.",
  x = 0, hjust = 0,
  gp = gpar(fontsize = 7.2, col = "#555555")
)

table_page <- arrangeGrob(table_title, tbl_grob, table_note,
  heights = unit.c(unit(0.35, "in"), unit(1, "null"), unit(0.32, "in")))

## =========================================================
## 8. 导出（全部在脚本同级 outputs 目录）
## =========================================================
out_plot_png <- file.path(OUT_DIR, paste0("Figure1_quadrant_scatter_", YEAR_MODE, ".png"))
out_plot_pdf <- file.path(OUT_DIR, paste0("Figure1_quadrant_scatter_", YEAR_MODE, ".pdf"))
out_table_pdf <- file.path(OUT_DIR, paste0("Table1_mismatch_top", LABEL_N_PER_QUADRANT, "_by_phenotype_", YEAR_MODE, ".pdf"))
out_combo_pdf <- file.path(OUT_DIR, paste0("Figure1_and_Table1_mismatch_", YEAR_MODE, ".pdf"))
out_docx <- file.path(OUT_DIR, paste0("Figure1_Table1_Mismatch_", YEAR_MODE, ".docx"))
out_diag <- file.path(OUT_DIR, paste0("diagnostic_table_", YEAR_MODE, ".tsv"))
out_main <- file.path(OUT_DIR, paste0("analysis_plot_df_", YEAR_MODE, ".tsv"))
out_table_tsv <- file.path(OUT_DIR, paste0("table1_top", LABEL_N_PER_QUADRANT, "_", YEAR_MODE, ".tsv"))

diagnostic_tbl <- tibble::tibble(
  step = c("WHO after filter", "joined population", "joined deaths", "joined incidence", "joined SDI", "final plot data"),
  n_country = c(
    nrow(who_base),
    sum(!is.na(plot_df0$population)),
    sum(!is.na(plot_df0$deaths)),
    sum(!is.na(plot_df0$incidence)),
    sum(!is.na(plot_df0$sdi)),
    sum(!is.na(plot_df0$population) & !is.na(plot_df0$deaths) & !is.na(plot_df0$incidence) & !is.na(plot_df0$sdi))
  )
)

write_tsv(diagnostic_tbl, out_diag)
write_tsv(plot_df_label, out_main)
write_tsv(table1_df, out_table_tsv)

ggsave(filename = out_plot_png, plot = p, width = 9, height = 6.4, dpi = 600, bg = "white")
ggsave(filename = out_plot_pdf, plot = p, width = 9, height = 6.4, device = "pdf", bg = "white")

pdf(file = out_table_pdf, width = 12, height = 8.5, onefile = TRUE)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.88, height = 0.84, just = c("center", "center")))
grid.draw(table_page)
popViewport()
dev.off()

pdf(file = out_combo_pdf, width = 12, height = 8.5, onefile = TRUE)
print(p)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.88, height = 0.84, just = c("center", "center")))
grid.draw(table_page)
popViewport()
dev.off()

doc <- read_docx()
doc <- body_add_par(doc, paste0("Figure 1. Maternal hemorrhage fatality versus reported red-cell availability in ", YEAR_MODE), style = "heading 1")
doc <- body_add_fpar(doc, fpar(external_img(out_plot_png, width = 6.6, height = 4.8), fp_p = fp_par(text.align = "center")))
doc <- body_add_par(doc, "Countries were classified into four mismatch phenotypes using sample medians.", style = "Normal")
doc <- body_add_break(doc)
doc <- body_add_par(doc, paste0("Table 1. Top ", LABEL_N_PER_QUADRANT, " countries within each mismatch phenotype in ", YEAR_MODE), style = "heading 1")
doc <- body_add_flextable(doc, value = flextable(table1_df))
print(doc, target = out_docx)

cat("[DONE] 分析完成。输出目录：\n", OUT_DIR, "\n", sep = "")
cat("主要文件：\n", out_plot_pdf, "\n", out_table_pdf, "\n", out_combo_pdf, "\n", out_docx, "\n", sep = "")
