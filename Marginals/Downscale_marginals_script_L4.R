# =====================================================================
# MODULAR DOWNSCALING SCRIPT (ALL-IN-ONE)
# =====================================================================

# ---------------------------------------------------------------------
# 1. LIBRARIES & CONFIGURATION
# ---------------------------------------------------------------------
library(dplyr)
library(lubridate)
library(lmomco)
library(DEoptim)
library(tools)
library(ggplot2)
library(tidyr)

# Auto-detect Google Drive path
gdrive_path <- paste(LETTERS[file.exists(paste0(LETTERS, ":/My Drive"))], ":/My Drive", sep = "")

APP_CONFIG <- list(
 paths = list(
  data_dir = file.path(gdrive_path, "General_Data/GSDR", "QC_d data - Germany"),
  export_dir = file.path(gdrive_path, "Academic_git/DownScaling/Marginals", "QC_d data - Germany", "Station_Results")
 ),
 parameters = list(
  k_star = 24,
  all_k = c(1:12, 24, 24*2, 24*3, 24*4, 24*5, 24*9, 24*10),
  max_H0 = 1,
  global_expected_models = c("Weibull_1p", "PowerLaw_1p", "Weibull_2p", "PowerLaw_2p"),
  global_expected_data = c("Calib. data", "Valid. data")
 ),
 deoptim_ctrl = DEoptim.control(
  trace = FALSE, itermax = 1000, NP = 50, reltol = 1e-11, steptol = 500
 )
)

# ---------------------------------------------------------------------
# 2. MATHEMATICAL MODELS & LOSS FUNCTIONS
# ---------------------------------------------------------------------
H_W_2p <- function(k, H0, k_star, q_k_star, a) { H0 * ((q_k_star / H0) ^ ((k / k_star)^a)) }
H_L_2p <- function(k, H0, k_star, q_k_star, b) { H0 * ((q_k_star / H0) ^ (log(1 + b * k) / log(1 + b * k_star))) }
H_W_1p_pdry <- function(k, k_star, q_k_star, a) { 0.99 * ((q_k_star / 0.99) ^ ((k / k_star)^a)) }
H_L_1p_pdry <- function(k, k_star, q_k_star, b) { 0.99 * ((q_k_star / 0.99) ^ (log(1 + b * k) / log(1 + b * k_star))) }

mse_W_2p <- function(par, k, y, k_star, q_k_star) mean((y - H_W_2p(k, par[1], k_star, q_k_star, par[2]))^2, na.rm=TRUE)
mse_L_2p <- function(par, k, y, k_star, q_k_star) mean((y - H_L_2p(k, par[1], k_star, q_k_star, par[2]))^2, na.rm=TRUE)
mse_W_1p <- function(par, k, y, k_star, q_k_star) mean((y - H_W_1p_pdry(k, k_star, q_k_star, par[1]))^2, na.rm=TRUE)
mse_L_1p <- function(par, k, y, k_star, q_k_star) mean((y - H_L_1p_pdry(k, k_star, q_k_star, par[1]))^2, na.rm=TRUE)

calculate_mse <- function(actual, predicted) { mean((actual - predicted)^2, na.rm = TRUE) }

# ---------------------------------------------------------------------
# 3. DATA PREP & AGGREGATION
# ---------------------------------------------------------------------
load_station_data <- function(filepath) {
 df <- read.table(filepath, skip = 21, col.names = "precip")
 df$precip[df$precip == -999] <- NA
 
 header_lines <- readLines(filepath, n = 21)
 start_date_str <- sub("Start datetime: ", "", grep("^Start datetime:", header_lines, value = TRUE))
 end_date_str   <- sub("End datetime: ", "", grep("^End datetime:", header_lines, value = TRUE))
 
 start_ts <- ymd_h(start_date_str)
 end_ts   <- ymd_h(end_date_str)
 df$date  <- seq(start_ts, end_ts, by = "hour")
 
 theoretical_records <- as.numeric(difftime(end_ts, start_ts, units = "hours")) + 1
 metadata <- list(
  Start_Datetime = start_ts, End_Datetime = end_ts,
  Theoretical_Records = theoretical_records,
  Explicit_NAs = sum(is.na(df$precip)),
  Missing_Rows = max(0, theoretical_records - nrow(df))
 )
 
 return(list(data = df, metadata = metadata))
}

aggregate_time_scales <- function(df, all_k) {
 df$precip_rate <- df$precip 
 df$hour_index <- 0:(nrow(df) - 1)
 clean_bins_list <- list()
 
 for (agg_length in all_k) {
  scale_name <- paste0("scale_", agg_length, "Hours")
  
  if (agg_length < 6) max_na_allowed <- 0 
  else if (agg_length >= 6 & agg_length <= 12) max_na_allowed <- 1
  else max_na_allowed <- min(ceiling(agg_length / 5), agg_length - 1) 
  
  valid_bins <- df %>%
   mutate(bin_id = hour_index %/% agg_length) %>%
   group_by(bin_id) %>%
   mutate(
    hours_in_bin = n(),
    na_count = sum(is.na(precip_rate)), 
    bin_mean_precip = ifelse(na_count <= max_na_allowed, mean(precip_rate, na.rm = TRUE), NA_real_),
    is_valid_bin = (!is.na(bin_mean_precip) & hours_in_bin == agg_length)
   ) %>%
   filter(is_valid_bin) %>%
   summarise(agg_date = min(date), bin_mean_precip = first(bin_mean_precip), .groups = "drop")
  
  clean_bins_list[[scale_name]] <- valid_bins
 }
 return(clean_bins_list)
}

compute_scale_statistics <- function(clean_bins_list) {
 summary_stats_list <- list()
 
 for (scale_name in names(clean_bins_list)) {
  x_all <- clean_bins_list[[scale_name]]$bin_mean_precip
  x_pos <- x_all[x_all > 0]
  
  p_zero <- sum(x_all == 0) / length(x_all)
  
  get_lmoms <- function(x) {
   if (length(x) >= 4) {
    lm <- lmoms(x)
    return(c(lm$lambdas[1:3], lm$ratios[2:4]))
   } else {
    return(rep(NA, 6))
   }
  }
  
  lm_all <- get_lmoms(x_all)
  lm_pos <- get_lmoms(x_pos)
  
  summary_stats_list[[scale_name]] <- data.frame(
   scale = scale_name, p_zero = p_zero,
   mean_all = mean(x_all, na.rm = TRUE), var_all = var(x_all, na.rm = TRUE),
   l1_all = lm_all[1], l2_all = lm_all[2], l3_all = lm_all[3], 
   t2_all = lm_all[4], t3_all = lm_all[5], t4_all = lm_all[6],
   mean_pos = ifelse(length(x_pos) > 0, mean(x_pos, na.rm = TRUE), NA), 
   var_pos = ifelse(length(x_pos) > 1, var(x_pos, na.rm = TRUE), NA),
   l1_pos = lm_pos[1], l2_pos = lm_pos[2], l3_pos = lm_pos[3], 
   t2_pos = lm_pos[4], t3_pos = lm_pos[5], t4_pos = lm_pos[6]
  )
 }
 return(bind_rows(summary_stats_list) %>% mutate(k_hours = as.numeric(gsub("[^0-9.]", "", scale))))
}

# ---------------------------------------------------------------------
# 4. MODEL FITTING (DEoptim)
# ---------------------------------------------------------------------
run_optimization <- function(model_name, loss_func, lower_bounds, upper_bounds, k_obs, y_obs, k_star, q_k_star, de_ctrl, stat_name) {
 fit_result <- DEoptim(fn = loss_func, lower = lower_bounds, upper = upper_bounds, k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
 
 best_params <- fit_result$optim$bestmem
 is_1p <- grepl("1p", model_name)
 is_weibull <- grepl("Weibull", model_name)
 
 return(data.frame(
  Statistic = stat_name, Model = model_name,
  H0 = if (is_1p) 1 else best_params[1],
  Param_a = if (is_weibull) (if (is_1p) best_params[1] else best_params[2]) else NA,
  Param_b = if (!is_weibull) (if (is_1p) best_params[1] else best_params[2]) else NA,
  MSE = fit_result$optim$bestval
 ))
}

fit_scaling_models <- function(final_summary_df, APP_CONFIG) {
 k_star <- APP_CONFIG$parameters$k_star
 de_ctrl <- APP_CONFIG$deoptim_ctrl
 max_H0 <- APP_CONFIG$parameters$max_H0
 
 fit_data <- final_summary_df %>% filter(k_hours >= k_star) %>% arrange(k_hours)
 
 # EXCLUDE t_pos from DEoptim fitting
 cols_to_exclude <- c("scale", "k_hours", "mean_all", "var_all", "mean_pos", "var_pos")
 stats_to_fit <- setdiff(names(fit_data), cols_to_exclude)
 
 all_fits_list <- list()
 
 for (stat in stats_to_fit) {
  y_obs <- fit_data[[stat]]
  if(all(is.na(y_obs))) next
  q_k_star <- fit_data[[stat]][fit_data$k_hours == k_star][1]
  if(is.na(q_k_star)) next
  
  stat_results <- list(
   run_optimization("Weibull_2p", mse_W_2p, c(1e-6, 1e-6), c(max_H0, 1000), fit_data$k_hours, y_obs, k_star, q_k_star, de_ctrl, stat),
   run_optimization("PowerLaw_2p", mse_L_2p, c(1e-6, 1e-6), c(max_H0, 1000), fit_data$k_hours, y_obs, k_star, q_k_star, de_ctrl, stat)
  )
  
  if (stat == "p_zero") {
   stat_results <- append(stat_results, list(
    run_optimization("Weibull_1p", mse_W_1p, 1e-6, 10, fit_data$k_hours, y_obs, k_star, q_k_star, de_ctrl, stat),
    run_optimization("PowerLaw_1p", mse_L_1p, 1e-6, 10, fit_data$k_hours, y_obs, k_star, q_k_star, de_ctrl, stat)
   ))
  }
  all_fits_list[[stat]] <- bind_rows(stat_results)
 }
 return(bind_rows(all_fits_list))
}

# ---------------------------------------------------------------------
# 5. PREDICTIONS & VALIDATION
# ---------------------------------------------------------------------
predict_scaling_values <- function(model_name, k_target, k_star, q_k_star, H0, param_a, param_b) {
 switch(model_name,
        "Weibull_2p"  = H_W_2p(k_target, H0, k_star, q_k_star, param_a),
        "PowerLaw_2p" = H_L_2p(k_target, H0, k_star, q_k_star, param_b),
        "Weibull_1p"  = H_W_1p_pdry(k_target, k_star, q_k_star, param_a),
        "PowerLaw_1p" = H_L_1p_pdry(k_target, k_star, q_k_star, param_b),
        rep(NA, length(k_target))
 )
}

generate_base_predictions <- function(final_summary_df, optimized_parameters_df, APP_CONFIG, custom_k = NULL) {
 k_star <- APP_CONFIG$parameters$k_star
 k_predict <- if (is.null(custom_k)) c(0.25, sort(unique(final_summary_df$k_hours))) else custom_k 
 
 predicted_list <- list()
 for (stat in unique(optimized_parameters_df$Statistic)) {
  q_k_star_val <- final_summary_df[[stat]][final_summary_df$k_hours == k_star]
  df_stat <- data.frame(Scale_k = k_predict, Statistic = stat, Actual = final_summary_df[[stat]][match(k_predict, final_summary_df$k_hours)])
  params <- optimized_parameters_df %>% filter(Statistic == stat)
  
  for (i in seq_len(nrow(params))) {
   df_stat[[params$Model[i]]] <- predict_scaling_values(params$Model[i], k_predict, k_star, q_k_star_val, params$H0[i], params$Param_a[i], params$Param_b[i])
  }
  predicted_list[[stat]] <- df_stat
 }
 return(bind_rows(predicted_list))
}

evaluate_models <- function(predictions_df, optimized_parameters_df, APP_CONFIG) {
 val_preds <- predictions_df %>% filter(Scale_k < APP_CONFIG$parameters$k_star)
 eval_results <- list()
 available_models <- setdiff(names(predictions_df), c("Scale_k", "Statistic", "Actual", "Data_Type"))
 
 for (stat in unique(val_preds$Statistic)) {
  stat_data <- val_preds %>% filter(Statistic == stat)
  if (all(is.na(stat_data$Actual))) next
  
  for (mod in available_models) {
   if (mod %in% names(stat_data) && !all(is.na(stat_data[[mod]]))) {
    mse_val <- calculate_mse(stat_data$Actual, stat_data[[mod]])
    orig_param <- optimized_parameters_df %>% filter(Statistic == stat, Model == mod)
    
    eval_results[[length(eval_results) + 1]] <- data.frame(
     Statistic = stat, Model = mod,
     H0 = if(nrow(orig_param) > 0) orig_param$H0[1] else NA,
     Param_a = if(nrow(orig_param) > 0) orig_param$Param_a[1] else NA,
     Param_b = if(nrow(orig_param) > 0) orig_param$Param_b[1] else NA,
     MSE_Calibration = if(nrow(orig_param) > 0) orig_param$MSE[1] else NA,
     mse_validation = mse_val
    )
   }
  }
 }
 evaluated_df <- bind_rows(eval_results)
 best_df <- evaluated_df %>% group_by(Statistic) %>% slice_min(order_by = mse_validation, n = 1, with_ties = FALSE) %>% ungroup()
 return(list(evaluated = evaluated_df, best = best_df))
}

# =====================================================================
# PART A: HIGH-SPEED PROCESSING LOOP (DATA ONLY)
# =====================================================================
# Run this block to churn through all your text files incredibly fast.

station_files <- list.files(path = APP_CONFIG$paths$data_dir, pattern = "\\.txt$", full.names = TRUE)

for (file in station_files) {
 
 station_name <- tools::file_path_sans_ext(basename(file))
 station_folder <- file.path(APP_CONFIG$paths$export_dir, station_name)
 
 tryCatch({
  raw_data <- load_station_data(file)
  clean_bins <- aggregate_time_scales(raw_data$data, APP_CONFIG$parameters$all_k)
  summary_df <- compute_scale_statistics(clean_bins)
  
  params_df <- fit_scaling_models(summary_df, APP_CONFIG)
  
  preds_df <- generate_base_predictions(summary_df, params_df, APP_CONFIG)
  base_eval  <- evaluate_models(preds_df, params_df, APP_CONFIG)
  
  eval_res <- list(
   evaluated = base_eval$evaluated,
   best      = base_eval$best
  )
  
  data_to_save <- list(
   final_summary = summary_df, optimized_params = params_df, 
   eval_results = eval_res, predictions = preds_df, metadata = raw_data$metadata
  )
  
  # Assuming 'save_station_data' is defined elsewhere in your environment
  save_station_data(station_folder, data_to_save)
  cat("Successfully processed data for:", station_name, "\n")
  
 }, error = function(e) {
  message("FAILED: Station ", station_name, " encountered an error. Detail: ", conditionMessage(e))
 })
}