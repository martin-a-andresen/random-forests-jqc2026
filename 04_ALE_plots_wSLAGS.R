# ===============================================================================
# ALE plots — two illustrative cases
#
# Case 1: diversity_25 at the SEGMENT level for TOB
# Case 2: leading CENSUS variable at the CT level for BNER
#
# Both models now include spatial lag predictors, matching the main analysis.
# ===============================================================================

.libPaths(c("/home/andresen/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(ranger)
library(iml)
library(dplyr)

set.seed(42)

impute_medians <- function(d, cols) {
  for (col in cols) {
    if (anyNA(d[[col]])) d[[col]][is.na(d[[col]])] <- median(d[[col]], na.rm = TRUE)
  }
  d
}

# Shared helper — build slag_cols cleanly
clean_slag_cols <- function(df) {
  s <- grep("_slag$", names(df), value = TRUE)
  s <- s[!grepl("_slag\\.[0-9]+$", s)]
  s <- s[!grepl("_slag_slag", s)]
  s
}

# ===============================================================================
# CASE 1: diversity_25, segment level, TOB
# ===============================================================================

cat("=== Case 1: diversity_25, segment level, TOB ===\n")

seg_path <- "/home/andresen/data/segments_final_RF_ready.csv"
df_seg   <- read.csv(seg_path)

network_cols   <- grep("^(segment_length|betweenness|closeness|degree|eigen|tortuosity|dead_end|dist_centre)",
                       names(df_seg), value = TRUE)
network_cols   <- network_cols[!grepl("_slag", network_cols)]
diversity_cols <- grep("^(diversity_|weighted_diversity_)[0-9]{2}$", names(df_seg), value = TRUE)
business_cols  <- grep("^(bar_lounge|cannabis|fast_food|food_retail|liquor_store|restaurant)[0-9]{2}$",
                       names(df_seg), value = TRUE)
da_ct_uid_cols  <- grep("^(DA_uid|CT_uid)_[0-9]{4}$", names(df_seg), value = TRUE)
da_ct_area_cols <- grep("^(DA_area|CT_area|DA_area_pct|CT_area_pct)_[0-9]{4}$", names(df_seg), value = TRUE)
census_cols     <- grep("^(DA_|CT_).*_(2001|2006|2011|2016|2021)$", names(df_seg), value = TRUE)
census_cols     <- setdiff(census_cols, c(da_ct_uid_cols, da_ct_area_cols))
change_cols     <- grep("_(abschange|pctchange|trend|cv)$", names(df_seg), value = TRUE)
slag_cols_seg   <- clean_slag_cols(df_seg)

place_predictor_cols_seg <- c(network_cols, diversity_cols, business_cols,
                              census_cols, change_cols, slag_cols_seg)
place_predictor_cols_seg <- place_predictor_cols_seg[!grepl("_slag_slag", place_predictor_cols_seg)]

all_na <- sapply(df_seg[, place_predictor_cols_seg], function(x) all(is.na(x)))
place_predictor_cols_seg <- place_predictor_cols_seg[!all_na]
cat(sprintf("  Segment predictors: %d (%d slag)\n",
            length(place_predictor_cols_seg),
            sum(grepl("_slag$", place_predictor_cols_seg))))

outcome_col_1 <- "TOB_2025"
rf_df_1 <- df_seg[, c(outcome_col_1, place_predictor_cols_seg)]
rf_df_1 <- rf_df_1[!is.na(rf_df_1[[outcome_col_1]]), ]
rf_df_1 <- impute_medians(rf_df_1, place_predictor_cols_seg)
rf_df_1[, place_predictor_cols_seg] <- lapply(rf_df_1[, place_predictor_cols_seg], as.numeric)

rf1 <- ranger(
  formula     = as.formula(paste(outcome_col_1, "~ .")),
  data        = rf_df_1[, c(outcome_col_1, place_predictor_cols_seg)],
  num.trees   = 500,
  importance  = "permutation",
  num.threads = 1
)

vi1 <- sort(rf1$variable.importance, decreasing = TRUE)
cat(sprintf("  Top predictor: %s (importance = %.5f)\n", names(vi1)[1], vi1[1]))
cat(sprintf("  diversity_25 rank: %d of %d\n",
            which(names(vi1) == "diversity_25"), length(vi1)))
if (names(vi1)[1] != "diversity_25") {
  cat("  NOTE: diversity_25 is not the top predictor in this run.\n")
}

X1 <- rf_df_1[, place_predictor_cols_seg]
y1 <- rf_df_1[[outcome_col_1]]
pred_fun_1 <- function(model, newdata) predict(model, data = newdata)$predictions
predictor1 <- Predictor$new(model = rf1, data = X1, y = y1, predict.fun = pred_fun_1)
ale1 <- FeatureEffect$new(predictor1, feature = "diversity_25", method = "ale")

png("/home/andresen/results/ALE_crimemix2025_segment_TOB.png",
    width = 2400, height = 1800, res = 300)
print(
  plot(ale1) +
    ggplot2::labs(
      title    = "ALE: 2025 crime mix, segment level, TOB (place-only model)",
      subtitle = "Accumulated local effect on predicted 2025 theft of bicycle count",
      x        = "Crime-type mix (2025)",
      y        = "ALE of predicted TOB count"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 14),
      axis.title    = ggplot2::element_text(size = 14),
      axis.text     = ggplot2::element_text(size = 14)
    )
)
dev.off()

write.csv(ale1$results,
          "/home/andresen/results/ALE_diversity25_segment_TOB.csv",
          row.names = FALSE)
cat("  Saved: ALE_crimemix2025_segment_TOB.png / .csv\n")

# ===============================================================================
# CASE 2: leading census variable, CT level, BNER
# ===============================================================================

cat("\n=== Case 2: leading census variable, CT level, BNER ===\n")

ct_path <- "/home/andresen/data/segments_aggregated_CT_2021.csv"
df_ct   <- read.csv(ct_path)

network_cols_ct  <- grep("^(segment_length|betweenness|closeness|degree|eigen|tortuosity|dead_end|dist_centre)",
                         names(df_ct), value = TRUE)
network_cols_ct  <- network_cols_ct[!grepl("_slag", network_cols_ct)]
diversity_cols_ct <- grep("^(diversity_|weighted_diversity_)[0-9]{2}$", names(df_ct), value = TRUE)
business_cols_ct  <- grep("^(bar_lounge|cannabis|fast_food|food_retail|liquor_store|restaurant)[0-9]{2}$",
                          names(df_ct), value = TRUE)
da_ct_uid_cols_ct  <- grep("^(DA_uid|CT_uid)_[0-9]{4}$", names(df_ct), value = TRUE)
da_ct_area_cols_ct <- grep("^(DA_area|CT_area|DA_area_pct|CT_area_pct)_[0-9]{4}$", names(df_ct), value = TRUE)
census_cols_ct     <- grep("^(CT_).*_(2001|2006|2011|2016|2021)$", names(df_ct), value = TRUE)
census_cols_ct     <- setdiff(census_cols_ct, c(da_ct_uid_cols_ct, da_ct_area_cols_ct))
change_cols_ct     <- grep("_(abschange|pctchange|trend|cv)$", names(df_ct), value = TRUE)
slag_cols_ct       <- clean_slag_cols(df_ct)

place_predictor_cols_ct <- c(network_cols_ct, diversity_cols_ct, business_cols_ct,
                             census_cols_ct, change_cols_ct, slag_cols_ct)
place_predictor_cols_ct <- place_predictor_cols_ct[!grepl("_slag_slag", place_predictor_cols_ct)]

all_na_ct <- sapply(df_ct[, place_predictor_cols_ct], function(x) all(is.na(x)))
place_predictor_cols_ct <- place_predictor_cols_ct[!all_na_ct]
cat(sprintf("  CT predictors: %d (%d slag)\n",
            length(place_predictor_cols_ct),
            sum(grepl("_slag$", place_predictor_cols_ct))))

outcome_col_2 <- "BNER_2025"
rf_df_2 <- df_ct[, c(outcome_col_2, place_predictor_cols_ct)]
rf_df_2 <- rf_df_2[!is.na(rf_df_2[[outcome_col_2]]), ]
rf_df_2 <- impute_medians(rf_df_2, place_predictor_cols_ct)
rf_df_2[, place_predictor_cols_ct] <- lapply(rf_df_2[, place_predictor_cols_ct], as.numeric)

rf2 <- ranger(
  formula     = as.formula(paste(outcome_col_2, "~ .")),
  data        = rf_df_2[, c(outcome_col_2, place_predictor_cols_ct)],
  num.trees   = 500,
  importance  = "permutation",
  num.threads = 1
)

vi2 <- sort(rf2$variable.importance, decreasing = TRUE)

# Identify top census variable (non-slag census + change)
census_and_change <- c(census_cols_ct, change_cols_ct)
census_and_change <- census_and_change[!grepl("_slag", census_and_change)]
census_vars_in_model <- intersect(names(vi2), census_and_change)
top_census_var <- census_vars_in_model[which.max(vi2[census_vars_in_model])]

cat(sprintf("  Top overall predictor:        %s (importance = %.5f)\n", names(vi2)[1], vi2[1]))
cat(sprintf("  Top CENSUS-category predictor: %s (importance = %.5f) <- plotting this\n",
            top_census_var, vi2[top_census_var]))

X2 <- rf_df_2[, place_predictor_cols_ct]
y2 <- rf_df_2[[outcome_col_2]]
pred_fun_2 <- function(model, newdata) predict(model, data = newdata)$predictions
predictor2 <- Predictor$new(model = rf2, data = X2, y = y2, predict.fun = pred_fun_2)
ale2 <- FeatureEffect$new(predictor2, feature = top_census_var, method = "ale")

# Clean up variable name for axis label
var_label <- gsub("CT_", "", top_census_var)
var_label <- gsub("_", " ", var_label)

png("/home/andresen/results/ALE_MedIncome_CT_BNER.png",
    width = 2400, height = 1800, res = 300)
print(
  plot(ale2) +
    ggplot2::labs(
      title    = sprintf("ALE: %s, CT level, BNER (place-only model)", top_census_var),
      subtitle = "Accumulated local effect on predicted 2025 residential burglary count",
      x        = var_label,
      y        = "ALE of predicted BNER count"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 14),
      axis.title    = ggplot2::element_text(size = 14),
      axis.text     = ggplot2::element_text(size = 14)
    )
)
dev.off()

write.csv(ale2$results,
          "/home/andresen/results/ALE_topcensus_CT_BNER.csv",
          row.names = FALSE)
cat(sprintf("  Saved: ALE_MedIncome_CT_BNER.png / .csv (variable: %s)\n", top_census_var))

cat("\nDone.\n")
cat(sprintf("Figure 3: diversity_25, segment TOB\n"))
cat(sprintf("Figure 4: %s, CT BNER\n", top_census_var))