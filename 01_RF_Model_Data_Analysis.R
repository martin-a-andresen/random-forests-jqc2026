# ===============================================================================
# RELOAD segments_final_RF_ready.csv — READY FOR RF MODELING
# ===============================================================================

library(dplyr)

data_dir <- "mydatadirectory"
df <- read.csv(file.path(data_dir, "segments_final_RF_ready.csv"))

cat(sprintf("Loaded: %d rows x %d columns\n", nrow(df), ncol(df)))

# ===============================================================================
# REBUILD VARIABLE GROUPS
# ===============================================================================

id_cols      <- c("segment_id", "hblock", "streetuse", "street_source")

network_cols <- grep("^(segment_length|betweenness|closeness|degree|eigen|tortuosity|dead_end|dist_centre)",
                     names(df), value = TRUE)

crime_cols   <- grep("^(TFV|Theft|Mischief|BNER|BNEC|TOV|TOB)_[0-9]{4}$",
                     names(df), value = TRUE)

diversity_cols <- grep("^(diversity_|weighted_diversity_)[0-9]{2}$",
                       names(df), value = TRUE)

business_cols  <- grep("^(bar_lounge|cannabis|fast_food|food_retail|liquor_store|restaurant)[0-9]{2}$",
                       names(df), value = TRUE)

da_ct_uid_cols  <- grep("^(DA_uid|CT_uid)_[0-9]{4}$", names(df), value = TRUE)
da_ct_area_cols <- grep("^(DA_area|CT_area|DA_area_pct|CT_area_pct)_[0-9]{4}$", names(df), value = TRUE)

census_cols_clean <- grep("^(DA_|CT_).*_(2001|2006|2011|2016|2021)$", names(df), value = TRUE)
census_cols_clean <- setdiff(census_cols_clean, c(da_ct_uid_cols, da_ct_area_cols))

change_cols <- grep("_(abschange|pctchange|trend|cv)$", names(df), value = TRUE)

# Final predictor set (excludes IDs, crime outcomes, DA/CT UIDs and area variables)
predictor_cols <- c(network_cols, diversity_cols, business_cols, census_cols_clean, change_cols)

cat(sprintf("\n=== Variable group summary ===\n"))
cat(sprintf("ID columns:      %d\n", length(id_cols)))
cat(sprintf("Network:         %d\n", length(network_cols)))
cat(sprintf("Crime (outcome): %d\n", length(crime_cols)))
cat(sprintf("Diversity:       %d\n", length(diversity_cols)))
cat(sprintf("Business:        %d\n", length(business_cols)))
cat(sprintf("DA/CT UIDs:      %d\n", length(da_ct_uid_cols)))
cat(sprintf("DA/CT Area:      %d\n", length(da_ct_area_cols)))
cat(sprintf("Census levels:   %d\n", length(census_cols_clean)))
cat(sprintf("Census change:   %d\n", length(change_cols)))
cat(sprintf("TOTAL PREDICTORS:%d\n", length(predictor_cols)))
cat(sprintf("Total columns:   %d\n", ncol(df)))

# ===============================================================================
# MEDIAN IMPUTATION FOR REMAINING PARTIAL NAs
# ===============================================================================

impute_medians <- function(df, pred_cols) {
  for (col in pred_cols) {
    if (anyNA(df[[col]])) {
      df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
    }
  }
  df
}

df <- impute_medians(df, predictor_cols)

cat(sprintf("\nAny remaining NAs in predictors after imputation: %d\n",
            sum(sapply(df[, predictor_cols], anyNA))))

cat("\nReady for RF modeling.\n")
cat(sprintf("Dataset: %d rows, %d predictors, %d crime outcome columns available\n",
            nrow(df), length(predictor_cols), length(crime_cols)))

####################################################################################

