## =========================================================
## 方案一：情景化分层面板预测（2030-2050）
## 仅包含：预测模型代码 + 投稿级图形输出
## =========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(lme4)
  library(splines)
})

## =========================================================
## 0. 输入与输出
## 说明：优先读取 country-year 长面板（推荐你提前整理好）
## 需要至少包含：iso3,country,year,rbc_per_1000,deaths,incidence,sdi
## 可选：dalys,region
## =========================================================
PANEL_FILE <- "W:/Current_Work/GBD_h/panel_country_year_mismatch.tsv"
BASELINE_YEAR <- 2018
FORECAST_YEARS <- 2030:2050
TRAIN_END_YEAR <- 2016

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) dirname(sub("^--file=", "", hit[1])) else getwd()
}
OUT_DIR <- file.path(get_script_dir(), "outputs_scenario_panel_forecast")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PANEL_FILE)) stop("找不到长面板文件：", PANEL_FILE)

panel <- read_tsv(PANEL_FILE, show_col_types = FALSE)

need_cols <- c("iso3", "country", "year", "rbc_per_1000", "deaths", "incidence", "sdi")
miss_cols <- setdiff(need_cols, names(panel))
if (length(miss_cols) > 0) stop("面板缺少必要列：", paste(miss_cols, collapse = ", "))

if (!"region" %in% names(panel)) panel <- panel %>% mutate(region = "Global")
if (!"dalys" %in% names(panel)) panel <- panel %>% mutate(dalys = NA_real_)

panel <- panel %>%
  mutate(
    year = as.integer(year),
    rbc_per_1000 = as.numeric(rbc_per_1000),
    deaths = as.numeric(deaths),
    incidence = as.numeric(incidence),
    sdi = as.numeric(sdi),
    dalys = as.numeric(dalys)
  ) %>%
  filter(!is.na(iso3), !is.na(country), !is.na(year)) %>%
  filter(!is.na(rbc_per_1000), rbc_per_1000 > 0) %>%
  filter(!is.na(deaths), deaths > 0) %>%
  filter(!is.na(incidence), incidence > 0) %>%
  filter(!is.na(sdi)) %>%
  mutate(dir_per_1000 = deaths / incidence * 1000)

## =========================================================
## 1. 构建历史窗口与验证集
## =========================================================
hist_df <- panel %>% arrange(iso3, year)
train_df <- hist_df %>% filter(year <= TRAIN_END_YEAR)
test_df  <- hist_df %>% filter(year > TRAIN_END_YEAR)

if (nrow(train_df) < 50) stop("训练样本太少，请检查 PANEL_FILE。")

# 训练集中心化年份，提升模型稳定性
year_center <- mean(train_df$year, na.rm = TRUE)
train_df <- train_df %>% mutate(year_c = year - year_center)
test_df  <- test_df %>% mutate(year_c = year - year_center)

# 根据训练集时间跨度自动选择年份项（避免 ns 在边界年过少时报错/警告）
n_year_unique <- dplyr::n_distinct(train_df$year)
year_term <- if (n_year_unique >= 6) "ns(year, df = 3)" else "year_c"

# 根据每个国家观测密度自动选择随机效应结构，避免 RE 参数不可识别
obs_by_country <- train_df %>% count(iso3, name = "n_obs")
can_use_random_slope <- nrow(obs_by_country) >= 20 && median(obs_by_country$n_obs, na.rm = TRUE) >= 4
re_term <- if (can_use_random_slope) "(1 + year_c | iso3)" else "(1 | iso3)"

## =========================================================
## 2. 分层面板模型（两个底层量 + incidence）
## =========================================================
m_rbc <- lmer(
  as.formula(paste0("log(rbc_per_1000) ~ ", year_term, " + sdi + ", re_term)),
  data = train_df,
  REML = FALSE,
  control = lmerControl(check.rankX = "ignore", check.conv.singular = "ignore")
)

m_dir <- lmer(
  as.formula(paste0("log(dir_per_1000) ~ ", year_term, " + sdi + log(rbc_per_1000) + ", re_term)),
  data = train_df,
  REML = FALSE,
  control = lmerControl(check.rankX = "ignore", check.conv.singular = "ignore")
)

m_inc <- lmer(
  as.formula(paste0("log(incidence) ~ ", year_term, " + sdi + ", re_term)),
  data = train_df,
  REML = FALSE,
  control = lmerControl(check.rankX = "ignore", check.conv.singular = "ignore")
)

## =========================================================
## 3. 验证（out-of-time）
## =========================================================
pred_test <- test_df %>%
  mutate(
    pred_rbc = exp(predict(m_rbc, newdata = test_df, allow.new.levels = TRUE)),
    pred_dir = exp(predict(m_dir, newdata = test_df, allow.new.levels = TRUE)),
    pred_inc = exp(predict(m_inc, newdata = test_df, allow.new.levels = TRUE))
  )

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
val_tbl <- tibble::tibble(
  metric = c("rbc_per_1000", "dir_per_1000", "incidence"),
  rmse = c(
    rmse(pred_test$rbc_per_1000, pred_test$pred_rbc),
    rmse(pred_test$dir_per_1000, pred_test$pred_dir),
    rmse(pred_test$incidence, pred_test$pred_inc)
  )
)
write_tsv(val_tbl, file.path(OUT_DIR, "validation_rmse_out_of_time.tsv"))

## =========================================================
## 4. 三情景预测：no intervention / reference / improvement
## =========================================================
country_base <- hist_df %>%
  filter(year == BASELINE_YEAR) %>%
  group_by(iso3, country, region) %>%
  summarise(
    sdi_base = median(sdi, na.rm = TRUE),
    rbc_base = median(rbc_per_1000, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(country_base) == 0) stop("BASELINE_YEAR 在数据里无可用记录：", BASELINE_YEAR)

future_grid <- tidyr::expand_grid(
  iso3 = unique(country_base$iso3),
  year = FORECAST_YEARS
) %>%
  left_join(country_base, by = "iso3") %>%
  mutate(year_c = year - year_center) %>%
  mutate(country = as.character(country)) %>%
  select(iso3, country, region, year, year_c, sdi_base, rbc_base)

# 参考情景：按历史趋势前推
ref_pred <- future_grid %>%
  mutate(
    sdi = sdi_base,
    rbc_ref = exp(predict(m_rbc, newdata = mutate(future_grid, sdi = sdi_base, rbc_per_1000 = pmax(rbc_base, 1e-6)), allow.new.levels = TRUE))
  )

# 冻结情景：rbc 固定在 baseline
freeze_pred <- ref_pred %>% mutate(rbc_freeze = rbc_base)

# 改善情景：到2050线性收敛到区域 75 分位
region_target <- country_base %>%
  group_by(region) %>%
  summarise(rbc_q75 = quantile(rbc_base, 0.75, na.rm = TRUE), .groups = "drop")

improve_pred <- ref_pred %>%
  left_join(region_target, by = "region") %>%
  mutate(
    prog = (year - min(FORECAST_YEARS)) / (max(FORECAST_YEARS) - min(FORECAST_YEARS)),
    rbc_improve = pmax(rbc_ref, rbc_ref + (rbc_q75 - rbc_ref) * pmax(0, prog))
  )

# 给定 rbc 场景，预测 dir/incidence/deaths
predict_burden <- function(dt, rbc_col, scenario_name) {
  # 先把模型依赖列显式写入 newdata，避免在同一个 mutate 中预测时找不到列
  dt2 <- dt %>%
    mutate(
      rbc_per_1000 = .data[[rbc_col]],
      sdi = sdi_base
    )
  
  dt2 <- dt2 %>%
    mutate(
      dir_pred = exp(predict(m_dir, newdata = dt2, allow.new.levels = TRUE)),
      inc_pred = exp(predict(m_inc, newdata = dt2, allow.new.levels = TRUE)),
      deaths_pred = dir_pred / 1000 * inc_pred,
      scenario = scenario_name
    )
  if (all(is.na(hist_df$dalys))) {
    dt2 <- dt2 %>% mutate(dalys_pred = NA_real_)
  } else {
    # 简化：用历史 DALYs/Death 比例中位数近似换算
    k <- median(hist_df$dalys / hist_df$deaths, na.rm = TRUE)
    dt2 <- dt2 %>% mutate(dalys_pred = deaths_pred * k)
  }
  dt2
}

pred_ref <- predict_burden(ref_pred, "rbc_ref", "reference")
pred_freeze <- predict_burden(freeze_pred, "rbc_freeze", "no_intervention")
pred_improve <- predict_burden(improve_pred, "rbc_improve", "improvement")

pred_all <- bind_rows(pred_ref, pred_freeze, pred_improve) %>%
  group_by(year) %>%
  mutate(
    z_rbc = (log(rbc_per_1000) - mean(log(rbc_per_1000), na.rm = TRUE)) / sd(log(rbc_per_1000), na.rm = TRUE),
    z_dir = (log(dir_pred) - mean(log(dir_pred), na.rm = TRUE)) / sd(log(dir_pred), na.rm = TRUE),
    mismatch_pred = z_dir - z_rbc
  ) %>%
  ungroup()

write_tsv(pred_all, file.path(OUT_DIR, "scenario_panel_forecast_2030_2050.tsv"))

## =========================================================
## 5. 危害翻译：持续 mismatch 的负担与 excess deaths
## =========================================================
burden_summary <- pred_all %>%
  group_by(year, scenario) %>%
  summarise(
    deaths = sum(deaths_pred, na.rm = TRUE),
    dalys = sum(dalys_pred, na.rm = TRUE),
    n_upper_left = sum(rbc_per_1000 < median(rbc_per_1000, na.rm = TRUE) & dir_pred >= median(dir_pred, na.rm = TRUE), na.rm = TRUE),
    n_upper_right = sum(rbc_per_1000 >= median(rbc_per_1000, na.rm = TRUE) & dir_pred >= median(dir_pred, na.rm = TRUE), na.rm = TRUE),
    .groups = "drop"
  )

excess_tbl <- burden_summary %>%
  select(year, scenario, deaths) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = deaths) %>%
  mutate(
    excess_vs_improve_ref = reference - improvement,
    excess_vs_improve_no_intervention = no_intervention - improvement
  )

write_tsv(burden_summary, file.path(OUT_DIR, "scenario_burden_summary_2030_2050.tsv"))
write_tsv(excess_tbl, file.path(OUT_DIR, "scenario_excess_deaths_vs_improvement.tsv"))

## =========================================================
## 6. 投稿图（美观清晰）
## =========================================================
# 图1：三情景死亡负担轨迹
p1 <- ggplot(burden_summary, aes(x = year, y = deaths, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(
    "no_intervention" = "#B2182B",
    "reference" = "#2166AC",
    "improvement" = "#1B7837"
  )) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Projected maternal hemorrhage deaths under three scenarios",
    subtitle = "Scenario-based hierarchical panel forecast (2030-2050)",
    x = NULL,
    y = "Projected deaths",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 12.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
    plot.title = element_text(face = "bold", size = 14.5),
    plot.subtitle = element_text(size = 10.5, color = "grey35"),
    legend.position = "top",
    plot.margin = margin(12, 12, 10, 12)
  )

# 图2：2050年 no intervention vs improvement 的国家层 excess 风险散点
d2050 <- pred_all %>% filter(year == 2050) %>% select(iso3, country, scenario, deaths_pred, mismatch_pred) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = c(deaths_pred, mismatch_pred)) %>%
  mutate(excess_deaths = deaths_pred_no_intervention - deaths_pred_improvement)

top_label <- d2050 %>% arrange(desc(excess_deaths)) %>% slice_head(n = 15)

p2 <- ggplot(d2050, aes(x = mismatch_pred_no_intervention, y = excess_deaths)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey70", linewidth = 0.4) +
  geom_point(color = "#8C2D04", alpha = 0.80, size = 2.6) +
  geom_text_repel(
    data = top_label,
    aes(label = country),
    size = 3.3,
    seed = 123,
    box.padding = 0.30,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Country-level excess deaths in 2050",
    subtitle = "No intervention vs improvement scenario",
    x = "Predicted mismatch score (no intervention)",
    y = "Excess deaths"
  ) +
  theme_minimal(base_size = 12.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
    plot.title = element_text(face = "bold", size = 14.5),
    plot.subtitle = element_text(size = 10.5, color = "grey35"),
    plot.margin = margin(12, 12, 10, 12)
  )

out_p1_png <- file.path(OUT_DIR, "Fig_scenario_projected_deaths_2030_2050.png")
out_p1_pdf <- file.path(OUT_DIR, "Fig_scenario_projected_deaths_2030_2050.pdf")
out_p2_png <- file.path(OUT_DIR, "Fig_country_excess_deaths_2050.png")
out_p2_pdf <- file.path(OUT_DIR, "Fig_country_excess_deaths_2050.pdf")

ggsave(out_p1_png, p1, width = 9.2, height = 5.8, dpi = 500, bg = "white")
ggsave(out_p1_pdf, p1, width = 9.2, height = 5.8, device = "pdf", bg = "white")
ggsave(out_p2_png, p2, width = 9.2, height = 6.2, dpi = 500, bg = "white")
ggsave(out_p2_pdf, p2, width = 9.2, height = 6.2, device = "pdf", bg = "white")

cat("[DONE] 预测完成。输出目录：\n", OUT_DIR, "\n", sep = "")
cat("关键文件：\n",
    file.path(OUT_DIR, "scenario_panel_forecast_2030_2050.tsv"), "\n",
    file.path(OUT_DIR, "scenario_burden_summary_2030_2050.tsv"), "\n",
    out_p1_pdf, "\n", out_p2_pdf, "\n", sep = "")
