# =====================================================================
# 1. MATHEMATICAL & OBJECTIVE FUNCTIONS
# =====================================================================

H_W_2p <- function(k, H0, k_star, q_k_star, a) {
 power_exponent <- (k / k_star)^a
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

H_L_2p <- function(k, H0, k_star, q_k_star, b) {
 ln_numerator <- log(1 + b * k)
 ln_denominator <- log(1 + b * k_star)
 power_exponent <- ln_numerator / ln_denominator
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

H_W_1p_pdry <- function(k, k_star, q_k_star, a) {
 H0 <- 0.99
 power_exponent <- (k / k_star)^a
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

H_L_1p_pdry <- function(k, k_star, q_k_star, b) {
 H0 <- 0.99
 ln_numerator <- log(1 + b * k)
 ln_denominator <- log(1 + b * k_star)
 power_exponent <- ln_numerator / ln_denominator
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

mse_W_2p <- function(par, k, y, k_star, q_k_star) mean((y - H_W_2p(k, par[1], k_star, q_k_star, par[2]))^2, na.rm=TRUE)
mse_L_2p <- function(par, k, y, k_star, q_k_star) mean((y - H_L_2p(k, par[1], k_star, q_k_star, par[2]))^2, na.rm=TRUE)
mse_W_1p <- function(par, k, y, k_star, q_k_star) mean((y - H_W_1p_pdry(k, k_star, q_k_star, par[1]))^2, na.rm=TRUE)
mse_L_1p <- function(par, k, y, k_star, q_k_star) mean((y - H_L_1p_pdry(k, k_star, q_k_star, par[1]))^2, na.rm=TRUE)

calculate_mse <- function(actual, predicted) {
 mean((actual - predicted)^2, na.rm = TRUE)
}

# --- NEW: CONSTRAINED MSE FOR t2_all (BLIND EXTRAPOLATION BOUNDS) ---
mse_W_2p_t2_constrained <- function(par, k, y, k_star, q_k_star, p0_k1_pred) {
 pred_t2_all <- H_W_2p(k, par[1], k_star, q_k_star, par[2])
 base_mse <- mean((y - pred_t2_all)^2, na.rm=TRUE)
 
 # Blind mathematical extrapolation to k=1
 pred_t2_all_k1 <- H_W_2p(1, par[1], k_star, q_k_star, par[2])
 
 # Guard against division by zero if p0 hits or exceeds 1 during extrapolation
 safe_p0 <- min(p0_k1_pred, 0.9999)
 
 derived_t2_pos_k1 <- (pred_t2_all_k1 - safe_p0) / (1 - safe_p0)
 
 # Theoretical boundaries for L-CV (t2) are [0, 1]
 if (is.na(derived_t2_pos_k1) || derived_t2_pos_k1 > 1.0 || derived_t2_pos_k1 < 0) {
  return(base_mse + 1e6) 
 }
 return(base_mse)
}

mse_L_2p_t2_constrained <- function(par, k, y, k_star, q_k_star, p0_k1_pred) {
 pred_t2_all <- H_L_2p(k, par[1], k_star, q_k_star, par[2])
 base_mse <- mean((y - pred_t2_all)^2, na.rm=TRUE)
 
 # Blind mathematical extrapolation to k=1
 pred_t2_all_k1 <- H_L_2p(1, par[1], k_star, q_k_star, par[2])
 
 safe_p0 <- min(p0_k1_pred, 0.9999)
 derived_t2_pos_k1 <- (pred_t2_all_k1 - safe_p0) / (1 - safe_p0)
 
 if (is.na(derived_t2_pos_k1) || derived_t2_pos_k1 > 1.0 || derived_t2_pos_k1 < 0) {
  return(base_mse + 1e6) 
 }
 return(base_mse)
}

# --- L-MOMENT CONVERSION FUNCTIONS ---
calc_t2_all <- function(p0, t2_pos) {
 p1 <- 1 - p0
 return(p0 + p1 * t2_pos)
}

calc_t3_all <- function(p0, t2_pos, t3_pos) {
 p1 <- 1 - p0
 numerator <- t2_pos * (p1^2) * t3_pos + 3 * p0 * p1 * t2_pos + 2 * (p0^2) - p0
 denominator <- p0 + p1 * t2_pos
 return(numerator / denominator)
}

calc_t4_all <- function(p0, t2_pos, t3_pos, t4_pos) {
 p1 <- 1 - p0
 numerator <- p0 * (5 * p0^2 - 5 * p0 + 1) + 
  3 * p0 * p1 * (3 * p0 - 1) * t2_pos + 
  5 * p0 * (p1^2) * t2_pos * t3_pos + 
  (p1^3) * t2_pos * t4_pos
 denominator <- p0 + p1 * t2_pos
 return(numerator / denominator)
}

calc_t2_pos <- function(p0, t2_all) {
 p1 <- 1 - p0
 return((t2_all - p0) / p1)
}

calc_t3_pos <- function(p0, t2_all, t3_all) {
 p1 <- 1 - p0
 numerator <- t3_all * t2_all - 3 * p0 * t2_all + (p0^2) + p0
 denominator <- p1 * (t2_all - p0)
 return(numerator / denominator)
}

calc_t4_pos <- function(p0, t2_all, t3_all, t4_all) {
 p1 <- 1 - p0
 numerator <- t4_all * t2_all - 
  5 * p0 * t3_all * t2_all + 
  3 * p0 * (2 * p0 + 1) * t2_all - 
  p0 * (p0^2 + 3 * p0 + 1)
 denominator <- (p1^2) * (t2_all - p0)
 return(numerator / denominator)
}


# =====================================================================
# 1. MATHEMATICAL & OBJECTIVE FUNCTIONS
# =====================================================================

H_W_2p <- function(k, H0, k_star, q_k_star, a) {
 power_exponent <- (k / k_star)^a
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

H_L_2p <- function(k, H0, k_star, q_k_star, b) {
 ln_numerator <- log(1 + b * k)
 ln_denominator <- log(1 + b * k_star)
 power_exponent <- ln_numerator / ln_denominator
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

H_W_1p_pdry <- function(k, k_star, q_k_star, a) {
 H0 <- 0.99
 power_exponent <- (k / k_star)^a
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

H_L_1p_pdry <- function(k, k_star, q_k_star, b) {
 H0 <- 0.99
 ln_numerator <- log(1 + b * k)
 ln_denominator <- log(1 + b * k_star)
 power_exponent <- ln_numerator / ln_denominator
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

# --- NEW: HYPERBOLIC FUNCTION FOR p_zero (STRICTLY CONVEX) ---
H_Hyperbolic_1p <- function(k, par_H0, k_star, q_k_star) {
 # Calculate 'c' using the 24h anchor point to ensure the curve hits q_k_star
 c_param <- (par_H0 - q_k_star) / (k_star * q_k_star)
 
 # Ensure c_param isn't negative to maintain physical boundaries
 c_param <- max(c_param, 1e-6) 
 
 # The pure hyperbolic curve
 return(par_H0 / (1 + c_param * k))
}

mse_W_2p <- function(par, k, y, k_star, q_k_star) mean((y - H_W_2p(k, par[1], k_star, q_k_star, par[2]))^2, na.rm=TRUE)
mse_L_2p <- function(par, k, y, k_star, q_k_star) mean((y - H_L_2p(k, par[1], k_star, q_k_star, par[2]))^2, na.rm=TRUE)
mse_W_1p <- function(par, k, y, k_star, q_k_star) mean((y - H_W_1p_pdry(k, k_star, q_k_star, par[1]))^2, na.rm=TRUE)
mse_L_1p <- function(par, k, y, k_star, q_k_star) mean((y - H_L_1p_pdry(k, k_star, q_k_star, par[1]))^2, na.rm=TRUE)
mse_Hyperbolic_1p <- function(par, k, y, k_star, q_k_star) mean((y - H_Hyperbolic_1p(k, par[1], k_star, q_k_star))^2, na.rm=TRUE)

calculate_mse <- function(actual, predicted) {
 mean((actual - predicted)^2, na.rm = TRUE)
}

# --- L-MOMENT CONVERSION FUNCTIONS ---
calc_t2_all <- function(p0, t2_pos) {
 p1 <- 1 - p0
 return(p0 + p1 * t2_pos)
}

calc_t3_all <- function(p0, t2_pos, t3_pos) {
 p1 <- 1 - p0
 numerator <- t2_pos * (p1^2) * t3_pos + 3 * p0 * p1 * t2_pos + 2 * (p0^2) - p0
 denominator <- p0 + p1 * t2_pos
 return(numerator / denominator)
}

calc_t4_all <- function(p0, t2_pos, t3_pos, t4_pos) {
 p1 <- 1 - p0
 numerator <- p0 * (5 * p0^2 - 5 * p0 + 1) + 
  3 * p0 * p1 * (3 * p0 - 1) * t2_pos + 
  5 * p0 * (p1^2) * t2_pos * t3_pos + 
  (p1^3) * t2_pos * t4_pos
 denominator <- p0 + p1 * t2_pos
 return(numerator / denominator)
}

calc_t2_pos <- function(p0, t2_all) {
 p1 <- 1 - p0
 return((t2_all - p0) / p1)
}

calc_t3_pos <- function(p0, t2_all, t3_all) {
 p1 <- 1 - p0
 numerator <- t3_all * t2_all - 3 * p0 * t2_all + (p0^2) + p0
 denominator <- p1 * (t2_all - p0)
 return(numerator / denominator)
}

calc_t4_pos <- function(p0, t2_all, t3_all, t4_all) {
 p1 <- 1 - p0
 numerator <- t4_all * t2_all - 
  5 * p0 * t3_all * t2_all + 
  3 * p0 * (2 * p0 + 1) * t2_all - 
  p0 * (p0^2 + 3 * p0 + 1)
 denominator <- (p1^2) * (t2_all - p0)
 return(numerator / denominator)
}


# =====================================================================
# 2. CORE STATION PROCESSING FUNCTION
# =====================================================================

process_station <- function(current_station_file, individual_dir) {
 
 station_name <- tools::file_path_sans_ext(basename(current_station_file))
 
 tryCatch({
  # --- 2.1 File Reading (Optimized with fread) ---
  df <- data.table::fread(current_station_file, skip = 21, col.names = "precip", data.table = FALSE)
  df$precip[df$precip == -999] <- NA
  
  header_lines <- readLines(current_station_file, n = 21)
  start_date_str <- sub("Start datetime: ", "", grep("^Start datetime:", header_lines, value = TRUE))
  end_date_str <- sub("End datetime: ", "", grep("^End datetime:", header_lines, value = TRUE))
  
  start_ts <- ymd_h(start_date_str)
  end_ts <- ymd_h(end_date_str)
  df$date <- seq(start_ts, end_ts, by = "hour")
  
  # --- 2.2 Metadata Calculation ---
  metadata_records <- as.numeric(sub("Number of records: ", "", grep("^Number of records:", header_lines, value = TRUE)))
  metadata_missing_pct <- as.numeric(sub("Percent missing data: ", "", grep("^Percent missing data:", header_lines, value = TRUE)))
  
  theoretical_records <- as.numeric(difftime(end_ts, start_ts, units = "hours")) + 1
  explicit_nas <- sum(is.na(df$precip)) 
  missing_rows <- max(0, theoretical_records - nrow(df))
  total_missing_data_points <- explicit_nas + missing_rows
  calculated_missing_pct <- (total_missing_data_points / theoretical_records) * 100
  
  # --- 2.3 & 2.4 ULTRA-FAST AGGREGATION & L-MOMENTS ---
  setDT(df)
  df[, hour_index := 0:(.N - 1)]
  
  all_k <- c(1:12, 24, 24*2, 24*3, 24*4, 24*5, 24*9, 24*10)
  
  summary_stats_list <- vector("list", length(all_k))
  scale_na_summary_list <- vector("list", length(all_k))
  
  for (i in seq_along(all_k)) {
   agg_length <- all_k[i]
   scale_name <- paste0("scale_", agg_length, "Hours")
   
   if (agg_length < 6) { max_na <- 0 } 
   else if (agg_length >= 6 & agg_length <= 12) { max_na <- 1 } 
   else { max_na <- min(ceiling(agg_length / 5), agg_length - 1) }
   
   agg_dt <- df[, .(
    n_total = .N,
    n_na = sum(is.na(precip)),
    sum_precip = sum(precip, na.rm = TRUE)
   ), by = .(bin_id = hour_index %/% agg_length)]
   
   valid_bins <- agg_dt[n_total == agg_length & n_na <= max_na, sum_precip / (n_total - n_na)]
   
   theor_bins <- theoretical_records %/% agg_length
   valid_bins_count <- length(valid_bins)
   
   scale_na_summary_list[[i]] <- data.frame(
    Scale_k_hours = agg_length, Max_NA_Allowed = max_na, Theoretical_Bins = theor_bins,
    Valid_Bins = valid_bins_count, Dropped_Bins = theor_bins - valid_bins_count,
    Dropped_Pct = ((theor_bins - valid_bins_count) / theor_bins) * 100
   )
   
   x_all <- valid_bins
   x_pos <- x_all[x_all > 0]
   
   get_lmoms <- function(x) {
    if (length(x) >= 4) {
     lm <- tryCatch(lmomco::lmoms(x), error = function(e) list(lambdas = rep(NA, 4), ratios = rep(NA, 6)))
     return(c(lm$lambdas[1:4], lm$ratios[2:4]))
    } else {
     return(rep(NA, 7))
    }
   }
   
   lm_all <- get_lmoms(x_all)
   lm_pos <- get_lmoms(x_pos)
   
   summary_stats_list[[i]] <- data.frame(
    scale = scale_name, k_hours = agg_length,
    p_zero = if(length(x_all) > 0) sum(x_all == 0) / length(x_all) else NA,
    mean_all = mean(x_all, na.rm = TRUE), var_all = var(x_all, na.rm = TRUE),
    l1_all = lm_all[1], l2_all = lm_all[2], l3_all = lm_all[3], l4_all = lm_all[4],
    t2_all = lm_all[5], t3_all = lm_all[6], t4_all = lm_all[7],
    mean_pos = if(length(x_pos) > 0) mean(x_pos, na.rm = TRUE) else NA,
    var_pos = if(length(x_pos) > 1) var(x_pos, na.rm = TRUE) else NA,
    l1_pos = lm_pos[1], l2_pos = lm_pos[2], l3_pos = lm_pos[3], l4_pos = lm_pos[4],
    t2_pos = lm_pos[5], t3_pos = lm_pos[6], t4_pos = lm_pos[7]
   )
  }
  
  scale_na_summary_df <- data.table::rbindlist(scale_na_summary_list)
  final_summary_df <- data.table::rbindlist(summary_stats_list)
  
  station_metadata_obj <- list(
   General_Info = list(Station = station_name, Start_Datetime = start_ts, End_Datetime = end_ts, Theoretical_Records = theoretical_records, Metadata_Records = metadata_records, Actual_Rows = nrow(df)),
   Missing_Data_Summary = list(Explicit_NAs = explicit_nas, Missing_Rows = missing_rows, Total_Missing_Points = total_missing_data_points, Calculated_Missing_Pct = calculated_missing_pct, Metadata_Missing_Pct = metadata_missing_pct),
   Scale_Data_Loss = as.data.frame(scale_na_summary_df)
  )
  
  # --- 2.5 DEoptim Model Fitting ---
  fit_data <- final_summary_df %>% filter(k_hours >= 24) %>% arrange(k_hours)
  k_obs <- fit_data$k_hours
  k_star <- 24
  
  stats_to_fit <- setdiff(names(fit_data), c("scale", "k_hours", "mean_all", "var_all", "mean_pos", "var_pos"))
  de_ctrl <- DEoptim.control(trace = FALSE, itermax = 1000, NP = 50, reltol = 1e-11, steptol = 500)
  all_fits_list <- list()
  
  for (stat in stats_to_fit) {
   y_obs <- fit_data[[stat]]
   if (all(is.na(y_obs))) next
   
   q_k_star <- fit_data[[stat]][fit_data$k_hours == k_star][1]
   if (is.na(q_k_star)) next
   max_H0 <- 1 
   
   fit_W_2p <- DEoptim(mse_W_2p, lower = c(1e-6, 1e-6), upper = c(max_H0, 1000), k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
   fit_L_2p <- DEoptim(mse_L_2p, lower = c(1e-6, 1e-6), upper = c(max_H0, 1000), k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
   
   res_row <- data.frame(Statistic = stat, Model = c("Weibull_2p", "PowerLaw_2p"), H0 = c(fit_W_2p$optim$bestmem[1], fit_L_2p$optim$bestmem[1]), Param_a = c(fit_W_2p$optim$bestmem[2], NA), Param_b = c(NA, fit_L_2p$optim$bestmem[2]), MSE = c(fit_W_2p$optim$bestval, fit_L_2p$optim$bestval))
   
   if (stat == "p_zero") {
    fit_W_1p <- DEoptim(mse_W_1p, lower = 1e-6, upper = 10, k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
    fit_L_1p <- DEoptim(mse_L_1p, lower = 1e-6, upper = 10, k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
    
    # --- NEW: Fit the Hyperbolic 1p model ---
    fit_Hyp_1p <- DEoptim(mse_Hyperbolic_1p, lower = 1e-6, upper = 1, k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
    
    res_1p <- data.frame(Statistic = stat, Model = c("Weibull_1p", "PowerLaw_1p", "Hyperbolic_1p"), 
                         H0 = c(1, 1, fit_Hyp_1p$optim$bestmem[1]), 
                         Param_a = c(fit_W_1p$optim$bestmem[1], NA, NA), 
                         Param_b = c(NA, fit_L_1p$optim$bestmem[1], NA), 
                         MSE = c(fit_W_1p$optim$bestval, fit_L_1p$optim$bestval, fit_Hyp_1p$optim$bestval))
    res_row <- bind_rows(res_row, res_1p)
   }
   all_fits_list[[stat]] <- res_row
  }
  
  optimized_parameters_df <- bind_rows(all_fits_list)
  
  # --- 2.6 Evaluate Validation Data ---
  validation_data <- final_summary_df %>% filter(k_hours < 24) %>% arrange(k_hours)
  empirical_24h <- final_summary_df %>% filter(k_hours == k_star)
  
  validation_results_list <- lapply(seq_len(nrow(optimized_parameters_df)), function(i) {
   current_fit <- optimized_parameters_df[i, ]
   actual_vals <- validation_data[[current_fit$Statistic]]
   q_k_star_val <- empirical_24h[[current_fit$Statistic]]
   
   if (all(is.na(actual_vals))) { current_fit$mse_validation <- NA; return(current_fit) }
   
   predicted_vals <- switch(current_fit$Model,
                            "Weibull_2p"    = H_W_2p(validation_data$k_hours, current_fit$H0, k_star, q_k_star_val, current_fit$Param_a),
                            "PowerLaw_2p"   = H_L_2p(validation_data$k_hours, current_fit$H0, k_star, q_k_star_val, current_fit$Param_b),
                            "Weibull_1p"    = H_W_1p_pdry(validation_data$k_hours, k_star, q_k_star_val, current_fit$Param_a),
                            "PowerLaw_1p"   = H_L_1p_pdry(validation_data$k_hours, k_star, q_k_star_val, current_fit$Param_b),
                            "Hyperbolic_1p" = H_Hyperbolic_1p(validation_data$k_hours, current_fit$H0, k_star, q_k_star_val),
                            rep(NA, length(validation_data$k_hours))
   )
   current_fit$mse_validation <- calculate_mse(actual_vals, predicted_vals)
   return(current_fit)
  })
  
  final_evaluated_models <- bind_rows(validation_results_list)
  best_models <- final_evaluated_models %>% group_by(Statistic) %>% slice_min(order_by = mse_validation, n = 1) %>% ungroup()
  
  # --- 2.7 Generate All Predictions & Derived L-Moments ---
  k_predict <- c(0.25, sort(unique(final_summary_df$k_hours))) 
  
  predicted_list <- lapply(unique(optimized_parameters_df$Statistic), function(stat) {
   q_k_star_val <- final_summary_df[[stat]][final_summary_df$k_hours == k_star]
   actual_vals <- final_summary_df[[stat]][match(k_predict, final_summary_df$k_hours)]
   params <- optimized_parameters_df %>% filter(Statistic == stat)
   df_stat <- data.frame(Scale_k = k_predict, Statistic = stat, Actual = actual_vals)
   
   for (i in seq_len(nrow(params))) {
    m_name <- params$Model[i]
    predicted_vals <- switch(m_name,
                             "Weibull_2p"    = H_W_2p(k_predict, params$H0[i], k_star, q_k_star_val, params$Param_a[i]),
                             "PowerLaw_2p"   = H_L_2p(k_predict, params$H0[i], k_star, q_k_star_val, params$Param_b[i]),
                             "Weibull_1p"    = H_W_1p_pdry(k_predict, k_star, q_k_star_val, params$Param_a[i]),
                             "PowerLaw_1p"   = H_L_1p_pdry(k_predict, k_star, q_k_star_val, params$Param_b[i]),
                             "Hyperbolic_1p" = H_Hyperbolic_1p(k_predict, params$H0[i], k_star, q_k_star_val),
                             rep(NA, length(k_predict)))
    df_stat[[m_name]] <- predicted_vals
   }
   return(df_stat)
  })
  
  all_predictions_df <- bind_rows(predicted_list)
  
  # Calculate derived L-moments based on best predicted lines AND Actual values
  derived_l_moments_df <- all_predictions_df %>%
   inner_join(best_models %>% select(Statistic, Best_Model = Model), by = "Statistic") %>%
   mutate(Pred = case_when(
    Best_Model == "Weibull_2p"    ~ Weibull_2p,
    Best_Model == "PowerLaw_2p"   ~ PowerLaw_2p,
    Best_Model == "Weibull_1p"    ~ Weibull_1p,
    Best_Model == "PowerLaw_1p"   ~ PowerLaw_1p,
    Best_Model == "Hyperbolic_1p" ~ Hyperbolic_1p,
    TRUE ~ NA_real_
   )) %>%
   select(Scale_k, Statistic, Pred, Actual) %>%
   pivot_wider(names_from = Statistic, values_from = c(Pred, Actual)) %>%
   mutate(
    derived_Pred_t2_all = if (all(c("Pred_p_zero", "Pred_t2_pos") %in% names(.))) calc_t2_all(Pred_p_zero, Pred_t2_pos) else NA,
    derived_Pred_t3_all = if (all(c("Pred_p_zero", "Pred_t2_pos", "Pred_t3_pos") %in% names(.))) calc_t3_all(Pred_p_zero, Pred_t2_pos, Pred_t3_pos) else NA,
    derived_Pred_t4_all = if (all(c("Pred_p_zero", "Pred_t2_pos", "Pred_t3_pos", "Pred_t4_pos") %in% names(.))) calc_t4_all(Pred_p_zero, Pred_t2_pos, Pred_t3_pos, Pred_t4_pos) else NA,
    
    derived_Pred_t2_pos = if (all(c("Pred_p_zero", "Pred_t2_all") %in% names(.))) calc_t2_pos(Pred_p_zero, Pred_t2_all) else NA,
    derived_Pred_t3_pos = if (all(c("Pred_p_zero", "Pred_t2_all", "Pred_t3_all") %in% names(.))) calc_t3_pos(Pred_p_zero, Pred_t2_all, Pred_t3_all) else NA,
    derived_Pred_t4_pos = if (all(c("Pred_p_zero", "Pred_t2_all", "Pred_t3_all", "Pred_t4_all") %in% names(.))) calc_t4_pos(Pred_p_zero, Pred_t2_all, Pred_t3_all, Pred_t4_all) else NA,
    
    derived_Actual_t2_all = if (all(c("Actual_p_zero", "Actual_t2_pos") %in% names(.))) calc_t2_all(Actual_p_zero, Actual_t2_pos) else NA,
    derived_Actual_t3_all = if (all(c("Actual_p_zero", "Actual_t2_pos", "Actual_t3_pos") %in% names(.))) calc_t3_all(Actual_p_zero, Actual_t2_pos, Actual_t3_pos) else NA,
    derived_Actual_t4_all = if (all(c("Actual_p_zero", "Actual_t2_pos", "Actual_t3_pos", "Actual_t4_pos") %in% names(.))) calc_t4_all(Actual_p_zero, Actual_t2_pos, Actual_t3_pos, Actual_t4_pos) else NA,
    
    derived_Actual_t2_pos = if (all(c("Actual_p_zero", "Actual_t2_all") %in% names(.))) calc_t2_pos(Actual_p_zero, Actual_t2_all) else NA,
    derived_Actual_t3_pos = if (all(c("Actual_p_zero", "Actual_t2_all", "Actual_t3_all") %in% names(.))) calc_t3_pos(Actual_p_zero, Actual_t2_all, Actual_t3_all) else NA,
    derived_Actual_t4_pos = if (all(c("Actual_p_zero", "Actual_t2_all", "Actual_t3_all", "Actual_t4_all") %in% names(.))) calc_t4_pos(Actual_p_zero, Actual_t2_all, Actual_t3_all, Actual_t4_all) else NA
   )
  
  # --- 2.8 Export Results ---
  station_folder <- file.path(individual_dir, station_name)
  dir.create(station_folder, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(final_summary_df, file.path(station_folder, "final_summary.rds"))
  saveRDS(optimized_parameters_df, file.path(station_folder, "optimized_parameters.rds"))
  saveRDS(final_evaluated_models, file.path(station_folder, "final_evaluated_models.rds"))
  saveRDS(best_models, file.path(station_folder, "best_models.rds"))
  saveRDS(all_predictions_df, file.path(station_folder, "all_predictions.rds"))
  saveRDS(derived_l_moments_df, file.path(station_folder, "derived_l_moments.rds")) 
  saveRDS(station_metadata_obj, file.path(station_folder, "extra_metadata.rds"))
  
  return(paste("Success:", station_name))
  
 }, error = function(e) {
  return(paste("FAILED:", station_name, "-", conditionMessage(e)))
 }) 
}


# =====================================================================
# 3. DIRECTORY SETUP & SEQUENTIAL EXECUTION
# =====================================================================

gdrive_path <- paste(LETTERS[file.exists(paste0(LETTERS, ":/My Drive"))], ":/My Drive", sep = "")
project_path <- file.path(gdrive_path, "Academic_git/DownScaling")

data_dir_name <- "QC_d data - Germany" 
data_dir <- file.path(gdrive_path, "General_Data/GSDR", data_dir_name)
individual_dir <- file.path(project_path, "Marginals", data_dir_name, "Testing") 

station_files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

cat("Starting sequential processing of", length(station_files), "files...\n\n")
station_files<-station_files[8]
results <- list()
for (file in station_files) {
 cat("Processing file:", basename(file), "...")
 res <- process_station(file, individual_dir)
 results[[basename(file)]] <- res
 cat(ifelse(grepl("^Success", res), " Done\n", paste("\n -> Error:", res, "\n")))
}

cat("\nProcessing Complete! Summary:\n")
print(table(grepl("^Success", unlist(results))))

# =====================================================================
# 4. PLOTTING
# =====================================================================
library(ggplot2)

plot_moments <- function(data, metric, x_var = "Scale_k") {
 
 col_pred           <- paste0("Pred_", metric)
 col_actual         <- paste0("Actual_", metric)
 col_derived_pred   <- paste0("derived_Pred_", metric)
 col_derived_actual <- paste0("derived_Actual_", metric)
 
 req_cols <- c(x_var, col_pred, col_actual, col_derived_pred, col_derived_actual)
 missing_cols <- setdiff(req_cols, names(data))
 if (length(missing_cols) > 0) {
  stop("Missing columns in your dataset: ", paste(missing_cols, collapse = ", "))
 }
 
 lbl_pred           <- paste("Predicted", metric)
 lbl_actual         <- paste("Actual", metric)
 lbl_derived_pred   <- paste("Derived Pred", metric)
 lbl_derived_actual <- paste("Derived Actual", metric)
 
 legend_colors <- setNames(
  c("blue", "black", "red", "orange"),
  c(lbl_pred, lbl_actual, lbl_derived_pred, lbl_derived_actual)
 )
 
 p <- ggplot(data = data) +
  
  geom_line(aes(x = .data[[x_var]], y = .data[[col_pred]], color = lbl_pred),
            linewidth = 1, alpha = 0.5, na.rm = TRUE) +
  
  geom_point(aes(x = .data[[x_var]], y = .data[[col_derived_pred]], color = lbl_derived_pred),
             shape = 16, size = 3, alpha = 0.5, na.rm = TRUE) +
  
  geom_point(aes(x = .data[[x_var]], y = .data[[col_actual]], color = lbl_actual),
             shape = 3, size = 3, stroke = 1, alpha = 0.8, na.rm = TRUE) +
  
  geom_point(aes(x = .data[[x_var]], y = .data[[col_derived_actual]], color = lbl_derived_actual),
             shape = 4, size = 3, stroke = 1, alpha = 0.8, na.rm = TRUE) +
  
  scale_color_manual(
   name = "Legend",
   values = legend_colors
  ) +
  
  labs(
   x = paste(x_var, "(hours)"),
   y = paste(metric, "values"),
   title = paste("Comparison of Actual vs Predicted and Derived:", metric)
  ) +
  
  theme_bw() +
  scale_x_log10()
 
 return(p)
}

library(ggplot2)
library(dplyr)

# Function to plot specifically the p_zero curves to check for the S-shape
plot_pzero_comparison <- function(predictions_df) {
 
 # Isolate only the p_zero data
 p_zero_data <- predictions_df %>% filter(Statistic == "p_zero")
 
 p <- ggplot(data = p_zero_data, aes(x = Scale_k)) +
  
  # 1. Actual Data Points (The truth)
  geom_point(aes(y = Actual, color = "Actual p_zero"), 
             shape = 3, size = 3, stroke = 1.2, na.rm = TRUE) +
  
  # 2. Weibull 2-parameter Fit (Prone to flattening out / S-Shape)
  geom_line(aes(y = Weibull_2p, color = "Weibull_2p (S-Shape)"), 
            linewidth = 1.2, linetype = "dashed", na.rm = TRUE) +
  
  # 3. Hyperbolic 1-parameter Fit (Strictly Convex)
  geom_line(aes(y = Hyperbolic_1p, color = "Hyperbolic_1p (Convex)"), 
            linewidth = 1.2, na.rm = TRUE) +
  
  # 4. Calibration Boundary (Show where extrapolation begins)
  geom_vline(xintercept = 24, linetype = "dotted", color = "darkgray", linewidth = 1) +
  annotate("text", x = 24, y = max(p_zero_data$Actual, na.rm=TRUE) * 0.9, 
           label = "Calibration Threshold (24h)", angle = 90, vjust = -0.5, color = "darkgray") +
  
  scale_color_manual(
   name = "Legend",
   values = c("Actual p_zero" = "black", 
              "Weibull_2p (S-Shape)" = "red", 
              "Hyperbolic_1p (Convex)" = "blue")
  ) +
  
  labs(
   x = "Scale k (hours)",
   y = "Probability of Zero Rainfall (p_zero)",
   title = "Comparison of p_zero Extrapolation: Weibull vs Hyperbolic"
  ) +
  
  theme_bw() +
  scale_x_log10()
 
 return(p)
}

# If you have the all_predictions_df in your environment, just run:
my_pzero_plot <- plot_pzero_comparison(all_predictions)
print(my_pzero_plot)

# Or, if you need to load a saved station's data first:
# my_data <- readRDS("path_to_your_testing_folder/YourStationName/all_predictions.rds")
# print(plot_pzero_comparison(my_data))