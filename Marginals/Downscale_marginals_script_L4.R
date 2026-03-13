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
# 3. ALGEBRAIC DERIVATIONS (Derive t_pos from t_all)
# ---------------------------------------------------------------------
calc_t2_pos <- function(p0, t2_all) {
  p1 <- 1 - p0
  ifelse(p1 == 0, NA, (t2_all - p0) / p1)
}

calc_t3_pos <- function(p0, t2_all, t3_all) {
  p1 <- 1 - p0
  num <- t3_all * t2_all - 3 * p0 * t2_all + (p0^2) + p0
  den <- p1 * (t2_all - p0)
  ifelse(den == 0 | is.na(den), NA, num / den)
}

calc_t4_pos <- function(p0, t2_all, t3_all, t4_all) {
  p1 <- 1 - p0
  num <- t4_all * t2_all - 5 * p0 * t3_all * t2_all + (6 * (p0^2) + 3 * p0) * t2_all - ((p0^3) + 3 * (p0^2) + p0)
  den <- (p1^2) * (t2_all - p0)
  ifelse(den == 0 | is.na(den), NA, num / den)
}

# ---------------------------------------------------------------------
# 4. DATA PREP & AGGREGATION
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
# 5. MODEL FITTING (DEoptim)
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
  cols_to_exclude <- c("scale", "k_hours", "mean_all", "var_all", "l1_all", "l2_all", "l3_all", "mean_pos", "var_pos", "l1_pos", "l2_pos", "l3_pos", "t2_pos", "t3_pos", "t4_pos")
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
# 6. PREDICTIONS & VALIDATION (Base + Derived Best t_pos)
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

derive_tpos_predictions <- function(final_summary_df, base_preds_df, best_models_df, custom_k = NULL) {
  k_predict <- if (is.null(custom_k)) c(0.25, sort(unique(final_summary_df$k_hours))) else custom_k 
  
  get_best_curve <- function(stat_name) {
    best_mod <- best_models_df$Model[best_models_df$Statistic == stat_name][1]
    if (is.na(best_mod)) return(rep(NA, length(k_predict)))
    stat_df <- base_preds_df %>% filter(Statistic == stat_name)
    if (nrow(stat_df) > 0 && best_mod %in% names(stat_df)) return(stat_df[[best_mod]])
    return(rep(NA, length(k_predict)))
  }
  
  p0_curve <- get_best_curve("p_zero")
  t2_all_curve <- get_best_curve("t2_all")
  t3_all_curve <- get_best_curve("t3_all")
  t4_all_curve <- get_best_curve("t4_all")
  
  res_list <- list()
  
  if (!all(is.na(p0_curve)) && !all(is.na(t2_all_curve))) {
    df_t2 <- data.frame(Scale_k = k_predict, Statistic = "t2_pos", Actual = final_summary_df$t2_pos[match(k_predict, final_summary_df$k_hours)])
    df_t2$Calculated_Best <- calc_t2_pos(p0_curve, t2_all_curve)
    res_list[["t2_pos"]] <- df_t2
  }
  if (!all(is.na(p0_curve)) && !all(is.na(t2_all_curve)) && !all(is.na(t3_all_curve))) {
    df_t3 <- data.frame(Scale_k = k_predict, Statistic = "t3_pos", Actual = final_summary_df$t3_pos[match(k_predict, final_summary_df$k_hours)])
    df_t3$Calculated_Best <- calc_t3_pos(p0_curve, t2_all_curve, t3_all_curve)
    res_list[["t3_pos"]] <- df_t3
  }
  if (!all(is.na(p0_curve)) && !all(is.na(t2_all_curve)) && !all(is.na(t3_all_curve)) && !all(is.na(t4_all_curve))) {
    df_t4 <- data.frame(Scale_k = k_predict, Statistic = "t4_pos", Actual = final_summary_df$t4_pos[match(k_predict, final_summary_df$k_hours)])
    df_t4$Calculated_Best <- calc_t4_pos(p0_curve, t2_all_curve, t3_all_curve, t4_all_curve)
    res_list[["t4_pos"]] <- df_t4
  }
  return(bind_rows(res_list))
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

# ---------------------------------------------------------------------
# 7. PLOTTING FUNCTIONS
# ---------------------------------------------------------------------
plot_statistic_dynamic <- function(stat_name, final_summary_df, smooth_preds_df, APP_CONFIG, y_label = stat_name, show_legend = FALSE, fixed_y = TRUE) {
  emp_data <- final_summary_df %>% select(k_hours, value = !!sym(stat_name)) %>% filter(!is.na(value)) %>% mutate(Data_Type = ifelse(k_hours >= APP_CONFIG$parameters$k_star, "Calib. data", "Valid. data"))
  if(nrow(emp_data) == 0) return(NULL)
  
  df_stat_smooth <- smooth_preds_df %>% filter(Statistic == stat_name)
  if(nrow(df_stat_smooth) == 0) return(NULL)
  
  existing_models <- intersect(c(APP_CONFIG$parameters$global_expected_models, "Calculated_Best"), names(df_stat_smooth))
  smooth_long <- df_stat_smooth %>% select(Scale_k, all_of(existing_models)) %>% pivot_longer(cols = all_of(existing_models), names_to = "Model", values_to = "value") %>% filter(!is.na(value))
  
  global_colors <- c("Weibull_1p" = "blue", "PowerLaw_1p" = "red", "Weibull_2p" = "blue", "PowerLaw_2p" = "red", "Calculated_Best" = "darkgreen")
  global_lines  <- c("Weibull_1p" = "solid", "PowerLaw_1p" = "solid", "Weibull_2p" = "dashed", "PowerLaw_2p" = "dashed", "Calculated_Best" = "solid")
  custom_breaks <- c(1, 6, 12, 24, 120, 240); custom_breaks <- custom_breaks[custom_breaks <= max(emp_data$k_hours, na.rm=TRUE)] 
  
  p <- ggplot() +
    geom_line(data = smooth_long, aes(x = Scale_k, y = value, color = Model, linetype = Model), linewidth = 1) +
    geom_point(data = emp_data, aes(x = k_hours, y = value, fill = Data_Type), shape = 21, color = "white", size = 3, stroke = 0.5) +
    geom_vline(xintercept = APP_CONFIG$parameters$k_star, linetype = "dashed", color = "darkgray", linewidth = 0.8) +
    scale_x_log10(breaks = custom_breaks, labels = custom_breaks) +
    scale_fill_manual(name = "", values = c("Calib. data" = "black", "Valid. data" = "orange")) +
    scale_color_manual(name = "", values = global_colors) +
    scale_linetype_manual(name = "", values = global_lines) +
    labs(x = "Temporal scale, k [h]", y = y_label) + theme_bw() + theme(legend.position = if(show_legend) "bottom" else "none")
  
  if (fixed_y) p <- p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
  return(p)
}

plot_discrete_points <- function(stat_name, preds_df, APP_CONFIG, y_label = stat_name, show_legend = TRUE, fixed_y = TRUE) {
  df_stat <- preds_df %>% filter(Statistic == stat_name)
  if(nrow(df_stat) == 0) return(NULL)
  
  existing_models <- intersect(c(APP_CONFIG$parameters$global_expected_models, "Calculated_Best"), names(df_stat))
  cols_to_gather <- c("Actual", existing_models)
  df_long <- df_stat %>% select(Scale_k, all_of(cols_to_gather)) %>% pivot_longer(cols = all_of(cols_to_gather), names_to = "Series", values_to = "value") %>% filter(!is.na(value))
  
  my_colors <- c("Actual" = "black", "Weibull_1p" = "blue", "PowerLaw_1p" = "red", "Weibull_2p" = "blue", "PowerLaw_2p" = "red", "Calculated_Best" = "darkgreen")
  my_shapes <- c("Actual" = 16, "Weibull_1p" = 17, "PowerLaw_1p" = 15, "Weibull_2p" = 17, "PowerLaw_2p" = 15, "Calculated_Best" = 18)
  
  custom_breaks <- c(1, 6, 12, 24, 120, 240); custom_breaks <- custom_breaks[custom_breaks <= max(df_stat$Scale_k, na.rm = TRUE)]
  
  p <- ggplot(df_long, aes(x = Scale_k, y = value, color = Series, shape = Series)) +
    geom_point(size = 3.5, alpha = 0.8) +
    geom_vline(xintercept = APP_CONFIG$parameters$k_star, linetype = "dashed", color = "darkgray", linewidth = 0.8) +
    scale_x_log10(breaks = custom_breaks, labels = custom_breaks) +
    scale_color_manual(values = my_colors) + scale_shape_manual(values = my_shapes) +
    labs(x = "Temporal scale, k [h]", y = y_label) + theme_bw() + theme(legend.position = if(show_legend) "bottom" else "none", legend.title = element_blank())
  
  if (fixed_y) p <- p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
  return(p)
}

plot_qq_dynamic <- function(stat_name, all_predictions_df, y_label = stat_name) {
  df_stat <- all_predictions_df %>% filter(Statistic == stat_name, !is.na(Actual)) %>% mutate(Data_Type = ifelse(Scale_k >= 24, "Calib. data", "Valid. data"))
  if(nrow(df_stat) == 0) return(NULL)
  
  model_cols <- setdiff(names(df_stat), c("Scale_k", "Statistic", "Actual", "Data_Type"))
  df_long <- df_stat %>% pivot_longer(cols = all_of(model_cols), names_to = "Model", values_to = "Predicted") %>% filter(!is.na(Predicted))
  if(nrow(df_long) == 0) return(NULL)
  
  min_val <- min(c(df_long$Actual, df_long$Predicted), na.rm = TRUE)
  max_val <- max(c(df_long$Actual, df_long$Predicted), na.rm = TRUE)
  buffer <- max(0.05, (max_val - min_val) * 0.05)
  
  p <- ggplot(df_long, aes(x = Actual, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgray", linewidth = 0.8) +
    geom_point(aes(fill = Data_Type), shape = 21, color = "white", size = 2.5, stroke = 0.5) +
    facet_wrap(~ Model) + coord_fixed(ratio = 1) + 
    scale_fill_manual(name = "", values = c("Calib. data" = "black", "Valid. data" = "orange")) +
    scale_x_continuous(limits = c(min_val - buffer, max_val + buffer)) +
    scale_y_continuous(limits = c(min_val - buffer, max_val + buffer)) +
    labs(title = paste("Empirical vs Predicted:", y_label), x = paste("Empirical (Actual)", y_label), y = paste("Predicted", y_label)) +
    theme_bw() + theme(panel.grid.major = element_line(linetype = "dotted", color = "lightgray"), panel.grid.minor = element_blank(), strip.background = element_rect(fill = "whitesmoke"), strip.text = element_text(face = "bold"), legend.position = "bottom")
  return(p)
}

# ---------------------------------------------------------------------
# 8. I/O MANAGER (Saving Data & Plots)
# ---------------------------------------------------------------------
save_station_data <- function(station_folder, df_list) {
  dir.create(station_folder, recursive = TRUE, showWarnings = FALSE)
  saveRDS(df_list$final_summary, file.path(station_folder, "final_summary.rds"))
  saveRDS(df_list$optimized_params, file.path(station_folder, "optimized_parameters.rds"))
  saveRDS(df_list$eval_results$evaluated, file.path(station_folder, "final_evaluated_models.rds"))
  saveRDS(df_list$eval_results$best, file.path(station_folder, "best_models.rds"))
  saveRDS(df_list$predictions, file.path(station_folder, "all_predictions.rds"))
  saveRDS(df_list$metadata, file.path(station_folder, "extra_metadata.rds"))
}

save_station_plots <- function(station_folder, summary_df, smooth_preds_df, preds_df, APP_CONFIG) {
  dir.create(station_folder, recursive = TRUE, showWarnings = FALSE)
  stats_all <- c("p_zero", "t2_all", "t3_all", "t4_all")
  stats_pos <- c("t2_pos", "t3_pos", "t4_pos")
  
  for (stat in stats_all) {
    p_fixed <- plot_statistic_dynamic(stat, summary_df, smooth_preds_df, APP_CONFIG, y_label = stat, show_legend = TRUE, fixed_y = TRUE)
    if (!is.null(p_fixed)) ggsave(file.path(station_folder, paste0(stat, "_scaling.png")), p_fixed + ggtitle(paste("Scaling Curve (Overall Process):", stat)), width = 8, height = 6, bg = "white")
  }
  
  for (stat in stats_pos) {
    p_fixed <- plot_discrete_points(stat, preds_df, APP_CONFIG, y_label = stat, show_legend = TRUE, fixed_y = TRUE)
    if (!is.null(p_fixed)) ggsave(file.path(station_folder, paste0(stat, "_scaling.png")), p_fixed + ggtitle(paste("Scaling Scatter (Positive Process):", stat)), width = 8, height = 6, bg = "white")
  }
  
  for (stat in stats_pos) {
    p_qq <- plot_qq_dynamic(stat, preds_df, y_label = stat)
    if (!is.null(p_qq)) ggsave(file.path(station_folder, paste0(stat, "_qq_validation.png")), p_qq + labs(title = paste("Positive Process (Calculated vs Empirical):", stat), x = paste("Empirical Actual", stat), y = paste("Calculated_Best", stat)), width = 10, height = 8, bg = "white")
  }
  
  for (stat in stats_all) {
    p_qq_all <- plot_qq_dynamic(stat, preds_df, y_label = stat)
    if (!is.null(p_qq_all)) ggsave(file.path(station_folder, paste0(stat, "_qq_validation.png")), p_qq_all + labs(title = paste("Overall Process (Fitted vs Empirical):", stat), x = paste("Empirical Actual", stat), y = paste("DEoptim Fitted", stat)), width = 10, height = 8, bg = "white")
  }
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
    
    base_preds <- generate_base_predictions(summary_df, params_df, APP_CONFIG)
    base_eval  <- evaluate_models(base_preds, params_df, APP_CONFIG)
    
    tpos_preds <- derive_tpos_predictions(summary_df, base_preds, base_eval$best)
    tpos_eval  <- evaluate_models(tpos_preds, params_df, APP_CONFIG)
    
    preds_df <- bind_rows(base_preds, tpos_preds)
    eval_res <- list(
      evaluated = bind_rows(base_eval$evaluated, tpos_eval$evaluated),
      best      = bind_rows(base_eval$best, tpos_eval$best)
    )
    
    data_to_save <- list(
      final_summary = summary_df, optimized_params = params_df, 
      eval_results = eval_res, predictions = preds_df, metadata = raw_data$metadata
    )
    
    save_station_data(station_folder, data_to_save)
    cat("Successfully processed data for:", station_name, "\n")
    
  }, error = function(e) {
    message("FAILED: Station ", station_name, " encountered an error. Detail: ", conditionMessage(e))
  })
}


# =====================================================================
# PART B: ON-DEMAND PLOTTER
# =====================================================================
# Run this block separately whenever you want to generate images for a station.

target_station <- "Station_001"   # <--- TYPE YOUR STATION NAME HERE

station_folder <- file.path(APP_CONFIG$paths$export_dir, target_station)

if (dir.exists(station_folder)) {
  summary_df <- readRDS(file.path(station_folder, "final_summary.rds"))
  params_df  <- readRDS(file.path(station_folder, "optimized_parameters.rds"))
  preds_df   <- readRDS(file.path(station_folder, "all_predictions.rds"))
  best_df    <- readRDS(file.path(station_folder, "best_models.rds"))
  
  k_smooth <- exp(seq(log(1), log(max(summary_df$k_hours, na.rm=TRUE)), length.out = 500))
  
  base_smooth <- generate_base_predictions(summary_df, params_df, APP_CONFIG, custom_k = k_smooth)
  tpos_smooth <- derive_tpos_predictions(summary_df, base_smooth, best_df, custom_k = k_smooth)
  smooth_preds_df <- bind_rows(base_smooth, tpos_smooth)
  
  save_station_plots(station_folder, summary_df, smooth_preds_df, preds_df, APP_CONFIG)
  cat("Done! Plots for", target_station, "are saved in", station_folder, "\n")
} else {
  message("Skipping plotting. Folder not found for station: ", target_station)
}