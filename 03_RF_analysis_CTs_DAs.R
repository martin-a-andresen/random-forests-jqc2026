# ===============================================================================
# Random Forest Analysis — Crime Count Prediction (Regression)
# DA-level and CT-level robustness check
# DRAC Fir Cluster — JupyterHub RStudio
#
# Mirrors RF_analysis.R (segment-level) exactly in structure/logic, applied to
# the aggregated DA and CT files instead of the segment-level file. Two models
# per crime type per level: (1) place characteristics only, (2) place + lags.
# ===============================================================================
.libPaths(c("/home/andresen/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(ranger)
library(randomForest)  # for OOB convergence curve only
library(dplyr)
library(ggplot2)

# ===============================================================================
# HELPER FUNCTION — median imputation for NA in predictors
# ===============================================================================

impute_medians <- function(df, pred_cols) {
  for (col in pred_cols) {
    if (anyNA(df[[col]])) {
      df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
    }
  }
  df
}

# ===============================================================================
# SETTINGS
# ===============================================================================

CRIME_TYPES  <- c("TOB", "BNEC", "TOV", "BNER", "Theft", "Mischief", "TFV")
OUTCOME_YEAR <- 2025
N_TREES      <- 1000
N_CORES      <- 1

LEVELS <- list(
  DA = list(
    data_path    = "/home/andresen/data/segments_aggregated_DA_2021.csv",
    results_path = "/home/andresen/results/RF_count_DA_SLAGS",
    uid_col      = "DA_uid_2021"
  ),
  CT = list(
    data_path    = "/home/andresen/data/segments_aggregated_CT_2021.csv",
    results_path = "/home/andresen/results/RF_count_CT_SLAGS",
    uid_col      = "CT_uid_2021"
  )
)

# ===============================================================================
# MAIN LOOP — run identical pipeline at both DA and CT levels
# ===============================================================================

all_summary_dfs <- list()

for (LEVEL in names(LEVELS)) {
  
  set.seed(42)  # reset per level for reproducibility, matches segment-level convention
  
  cfg          <- LEVELS[[LEVEL]]
  data_path    <- cfg$data_path
  results_path <- cfg$results_path
  uid_col      <- cfg$uid_col
  dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\n################################################################\n"))
  cat(sprintf("LEVEL: %s\n", LEVEL))
  cat(sprintf("################################################################\n"))
  
  cat("Loading data...\n")
  df <- read.csv(data_path)
  cat(sprintf("Loaded: %d rows x %d columns\n", nrow(df), ncol(df)))
  
  # ===============================================================================
  # DEFINE VARIABLE GROUPS — identical regex logic to segment-level script
  # ===============================================================================
  
  network_cols   <- grep("^(segment_length|betweenness|closeness|degree|eigen|tortuosity|dead_end|dist_centre)",
                         names(df), value = TRUE)
  network_cols   <- network_cols[!grepl("_slag", network_cols)]  # exclude slag and dedup variants
  diversity_cols <- grep("^(diversity_|weighted_diversity_)[0-9]{2}$",
                         names(df), value = TRUE)
  business_cols  <- grep("^(bar_lounge|cannabis|fast_food|food_retail|liquor_store|restaurant)[0-9]{2}$",
                         names(df), value = TRUE)
  da_ct_uid_cols  <- grep("^(DA_uid|CT_uid)_[0-9]{4}$", names(df), value = TRUE)
  da_ct_area_cols <- grep("^(DA_area|CT_area|DA_area_pct|CT_area_pct)_[0-9]{4}$", names(df), value = TRUE)
  census_cols     <- grep("^(DA_|CT_).*_(2001|2006|2011|2016|2021)$", names(df), value = TRUE)
  census_cols     <- setdiff(census_cols, c(da_ct_uid_cols, da_ct_area_cols))
  change_cols     <- grep("_(abschange|pctchange|trend|cv)$", names(df), value = TRUE)
  slag_cols       <- grep("_slag$", names(df), value = TRUE)
  slag_cols       <- slag_cols[!grepl("_slag\\.[0-9]+$", slag_cols)]  # drop R-deduplicated duplicates
  slag_cols       <- slag_cols[!grepl("_slag_slag", slag_cols)]          # drop double-lagged columns
  
  # "Place characteristics" predictors — no lagged crime, no n_segments
  place_predictor_cols <- c(network_cols, diversity_cols, business_cols,
                            census_cols, change_cols, slag_cols)
  place_predictor_cols <- place_predictor_cols[!grepl("_slag_slag", place_predictor_cols)]  # belt-and-suspenders
  
  # Drop all-NA columns as a final safety check
  all_na_cols <- sapply(df[, place_predictor_cols], function(x) all(is.na(x)))
  if (any(all_na_cols)) {
    cat(sprintf("Dropping %d all-NA predictor columns\n", sum(all_na_cols)))
    place_predictor_cols <- place_predictor_cols[!all_na_cols]
  }
  
  cat(sprintf("Place-characteristic predictors: %d\n", length(place_predictor_cols)))
  cat(sprintf("  network: %d | diversity: %d | business: %d | census: %d | change: %d | slag: %d\n",
              length(network_cols), length(diversity_cols), length(business_cols),
              length(census_cols), length(change_cols),
              sum(grepl("_slag$", place_predictor_cols))))
  cat(sprintf("n_segments range: %d - %d (NOT included as predictor)\n",
              min(df$n_segments), max(df$n_segments)))
  
  # ===============================================================================
  # OOB CONVERGENCE DIAGNOSTIC — TFV only, using randomForest package
  # ===============================================================================
  
  cat(sprintf("\n--- OOB Convergence Diagnostic (%s level, TFV, full model with lags) ---\n", LEVEL))
  
  ct_diag <- "TFV"
  outcome_col_diag <- paste0(ct_diag, "_", OUTCOME_YEAR)
  
  lag_years <- 2003:(OUTCOME_YEAR - 1)
  lag_cols_diag <- paste0(ct_diag, "_", lag_years)
  lag_cols_diag <- lag_cols_diag[lag_cols_diag %in% names(df)]
  
  full_predictor_cols_diag <- c(place_predictor_cols, lag_cols_diag)
  
  diag_df <- df[, c(outcome_col_diag, full_predictor_cols_diag)]
  diag_df <- diag_df[!is.na(diag_df[[outcome_col_diag]]), ]
  diag_df <- impute_medians(diag_df, full_predictor_cols_diag)
  
  cat(sprintf("Outcome: %s (n = %d, mean = %.2f, sd = %.2f, max = %.2f)\n",
              outcome_col_diag, nrow(diag_df),
              mean(diag_df[[outcome_col_diag]]), sd(diag_df[[outcome_col_diag]]),
              max(diag_df[[outcome_col_diag]])))
  
  cat(sprintf("  Fitting randomForest for convergence curve (%d trees)...\n", N_TREES))
  rf_diag <- randomForest(
    x          = diag_df[, full_predictor_cols_diag],
    y          = diag_df[[outcome_col_diag]],
    ntree      = N_TREES,
    importance = FALSE,
    do.trace   = 50
  )
  
  oob_df <- data.frame(n_trees = 1:N_TREES, oob_mse = rf_diag$mse)
  
  p_conv <- ggplot(oob_df, aes(x = n_trees, y = oob_mse)) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    labs(title    = sprintf("OOB MSE Convergence — %s (%s level)", outcome_col_diag, LEVEL),
         subtitle = sprintf("randomForest package, %d trees, full model with lags", N_TREES),
         x        = "Number of Trees", y = "OOB Mean Squared Error") +
    theme_bw()
  
  ggsave(file.path(results_path, sprintf("oob_convergence_%s_%s.png", LEVEL, outcome_col_diag)),
         p_conv, width = 8, height = 5, dpi = 150)
  cat(sprintf("  Saved: oob_convergence_%s_%s.png\n", LEVEL, outcome_col_diag))
  cat("  OOB MSE, last 10 tree counts:\n")
  print(tail(oob_df, 10))
  
  # ===============================================================================
  # MAIN RF ANALYSIS — all crime types, ranger regression
  # ===============================================================================
  
  cat(sprintf("\n--- Main Random Forest Analysis (%s level, ranger, regression) ---\n", LEVEL))
  
  results_list <- list()
  
  for (ct in CRIME_TYPES) {
    
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Level: %s | Crime type: %s\n", LEVEL, ct))
    cat(sprintf("========================================\n"))
    
    outcome_col <- paste0(ct, "_", OUTCOME_YEAR)
    
    if (!outcome_col %in% names(df)) {
      cat(sprintf("  SKIPPED — outcome column %s not found\n", outcome_col))
      next
    }
    
    lag_cols <- paste0(ct, "_", lag_years)
    lag_cols <- lag_cols[lag_cols %in% names(df)]
    
    full_predictor_cols <- c(place_predictor_cols, lag_cols)
    
    rf_df <- df[, c(outcome_col, full_predictor_cols)]
    rf_df <- rf_df[!is.na(rf_df[[outcome_col]]), ]
    rf_df <- impute_medians(rf_df, full_predictor_cols)
    
    cat(sprintf("Outcome: %s\n", outcome_col))
    cat(sprintf("  n = %d, mean = %.2f, sd = %.2f, min = %.2f, max = %.2f\n",
                nrow(rf_df), mean(rf_df[[outcome_col]]), sd(rf_df[[outcome_col]]),
                min(rf_df[[outcome_col]]), max(rf_df[[outcome_col]])))
    cat(sprintf("Predictors: %d place characteristics + %d lag years = %d total\n",
                length(place_predictor_cols), length(lag_cols), length(full_predictor_cols)))
    
    # ── MODEL 1: Place characteristics only ──────────────────────────────────
    cat(sprintf("\n[Model 1] Place characteristics only (%d predictors)...\n",
                length(place_predictor_cols)))
    
    t_start <- proc.time()
    rf_place <- ranger(
      formula     = as.formula(paste(outcome_col, "~ .")),
      data        = rf_df[, c(outcome_col, place_predictor_cols)],
      num.trees   = N_TREES,
      importance  = "permutation",
      num.threads = N_CORES,
      verbose     = TRUE
    )
    t_elapsed <- proc.time() - t_start
    
    rsq_place <- rf_place$r.squared
    cat(sprintf("Done in %.1f sec | OOB R-squared: %.4f | OOB MSE: %.4f\n",
                t_elapsed["elapsed"], rsq_place, rf_place$prediction.error))
    
    # ── MODEL 2: Place characteristics + lagged crime ────────────────────────
    cat(sprintf("\n[Model 2] Place characteristics + %d lag years...\n", length(lag_cols)))
    
    t_start <- proc.time()
    rf_full <- ranger(
      formula     = as.formula(paste(outcome_col, "~ .")),
      data        = rf_df[, c(outcome_col, full_predictor_cols)],
      num.trees   = N_TREES,
      importance  = "permutation",
      num.threads = N_CORES,
      verbose     = TRUE
    )
    t_elapsed <- proc.time() - t_start
    
    rsq_full <- rf_full$r.squared
    cat(sprintf("Done in %.1f sec | OOB R-squared: %.4f | OOB MSE: %.4f\n",
                t_elapsed["elapsed"], rsq_full, rf_full$prediction.error))
    
    cat(sprintf("\nR-squared improvement from adding lags: %.4f -> %.4f (+%.4f)\n",
                rsq_place, rsq_full, rsq_full - rsq_place))
    
    # ── Variable importance ───────────────────────────────────────────────────
    
    vi_place <- data.frame(
      variable   = names(rf_place$variable.importance),
      importance = rf_place$variable.importance
    ) |> dplyr::arrange(desc(importance))
    
    vi_full <- data.frame(
      variable   = names(rf_full$variable.importance),
      importance = rf_full$variable.importance
    ) |> dplyr::arrange(desc(importance))
    
    cat(sprintf("\nTop 15 variables — Model 1 (place only):\n"))
    print(head(vi_place, 15))
    cat(sprintf("\nTop 15 variables — Model 2 (place + lags):\n"))
    print(head(vi_full, 15))
    
    # ── Save outputs ───────────────────────────────────────────────────────────
    
    write.csv(vi_place,
              file.path(results_path, sprintf("varimp_%s_placeOnly_%s.csv", LEVEL, ct)),
              row.names = FALSE)
    write.csv(vi_full,
              file.path(results_path, sprintf("varimp_%s_withLags_%s.csv", LEVEL, ct)),
              row.names = FALSE)
    
    p_vi_place <- ggplot(head(vi_place, 30), aes(x = reorder(variable, importance), y = importance)) +
      geom_col(fill = "steelblue") + coord_flip() +
      labs(title = sprintf("Variable Importance (Place Only) — %s (%s level)", outcome_col, LEVEL),
           subtitle = sprintf("OOB R-squared: %.3f", rsq_place), x = NULL, y = "Permutation Importance") +
      theme_bw(base_size = 10)
    ggsave(file.path(results_path, sprintf("varimp_%s_placeOnly_%s.png", LEVEL, ct)),
           p_vi_place, width = 8, height = 8, dpi = 150)
    
    p_vi_full <- ggplot(head(vi_full, 30), aes(x = reorder(variable, importance), y = importance)) +
      geom_col(fill = "darkorange") + coord_flip() +
      labs(title = sprintf("Variable Importance (With Lags) — %s (%s level)", outcome_col, LEVEL),
           subtitle = sprintf("OOB R-squared: %.3f", rsq_full), x = NULL, y = "Permutation Importance") +
      theme_bw(base_size = 10)
    ggsave(file.path(results_path, sprintf("varimp_%s_withLags_%s.png", LEVEL, ct)),
           p_vi_full, width = 8, height = 8, dpi = 150)
    
    results_list[[ct]] <- list(
      crime_type = ct, outcome_col = outcome_col, n_obs = nrow(rf_df),
      n_lag_years = length(lag_cols),
      rsq_place = rsq_place, rsq_full = rsq_full, rsq_improvement = rsq_full - rsq_place,
      mse_place = rf_place$prediction.error, mse_full = rf_full$prediction.error,
      top10_place = head(vi_place$variable, 10), top10_full = head(vi_full$variable, 10)
    )
    
    cat(sprintf("Saved outputs for %s level, %s\n", LEVEL, ct))
  }
  
  # ===============================================================================
  # SUMMARY TABLE — this level
  # ===============================================================================
  
  cat(sprintf("\n--- Summary across all crime types (%s level) ---\n", LEVEL))
  
  summary_df <- data.frame(
    level           = LEVEL,
    crime_type      = names(results_list),
    outcome         = sapply(results_list, function(x) x$outcome_col),
    n_obs           = sapply(results_list, function(x) x$n_obs),
    n_lag_years     = sapply(results_list, function(x) x$n_lag_years),
    rsq_place_only  = round(sapply(results_list, function(x) x$rsq_place), 4),
    rsq_with_lags   = round(sapply(results_list, function(x) x$rsq_full), 4),
    rsq_improvement = round(sapply(results_list, function(x) x$rsq_improvement), 4),
    mse_place_only  = round(sapply(results_list, function(x) x$mse_place), 4),
    mse_with_lags   = round(sapply(results_list, function(x) x$mse_full), 4)
  )
  
  print(summary_df)
  write.csv(summary_df, file.path(results_path, sprintf("RF_summary_%s.csv", LEVEL)), row.names = FALSE)
  cat(sprintf("\n%s level done! Results saved to: %s\n", LEVEL, results_path))
  
  all_summary_dfs[[LEVEL]] <- summary_df
}

# ===============================================================================
# COMBINED SUMMARY
# ===============================================================================

combined_summary <- bind_rows(all_summary_dfs)
cat("\n--- Combined DA + CT summary ---\n")
print(combined_summary)
write.csv(combined_summary,
          "/home/andresen/results/RF_summary_DA_CT_combined_SLAGS.csv",
          row.names = FALSE)

cat("\nAll levels done.\n")
