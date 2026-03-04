
library(dplyr)
library(lubridate)
library(lmomco)
library(DEoptim)
library(tools)


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
 H0 <- 1 
 power_exponent <- (k / k_star)^a
 base_fraction <- q_k_star / H0
 H0 * (base_fraction ^ power_exponent)
}

H_L_1p_pdry <- function(k, k_star, q_k_star, b) {
 H0 <- 1 
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

# =====================================================================
# =====================================================================


gdrive_path <- paste(LETTERS[file.exists(paste0(LETTERS, ":/My Drive"))], ":/My Drive", sep = "")
project_path <- file.path(gdrive_path, "Academic_git/DownScaling")
generalData_path <- file.path(gdrive_path, "General_Data")
individual_dir <- file.path(project_path, "Station_Outputs") # Export directory


data_dir <- file.path(generalData_path, "GSDR/QC_d data - Germany")
station_files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

# =====================================================================
# =====================================================================
station_files<-station_files[1:3]

for (current_station_file in station_files) {
 
 # Extract current station name for export later
 station_name <- tools::file_path_sans_ext(basename(current_station_file))

 df <- read.table(current_station_file, skip = 21, col.names = "precip")
 df$precip[df$precip == -999] <- NA
 
 header_lines <- readLines(current_station_file, n = 21)
 start_date_str <- sub("Start datetime: ", "", grep("^Start datetime:", header_lines, value = TRUE))
 end_date_str <- sub("End datetime: ", "", grep("^End datetime:", header_lines, value = TRUE))
 
 start_ts <- ymd_h(start_date_str)
 end_ts <- ymd_h(end_date_str)
 df$date <- seq(start_ts, end_ts, by = "hour")
 

 Base_time_chunks_per_hour <- 4
 df$precip_rate <- df$precip 
 df$hour_index <- 0:(nrow(df) - 1)
 
 all_k <- c(empirical_k_validation, empirical_k_training)
 
 inspection_list <- list()

for (i in seq_along(all_k)) {
  
  agg_length <- all_k[i]
  
  # 1. Define the NA tolerance based on the scale
  if (agg_length <= 12) {
    max_na_allowed <- 0
  } else {
    # Linear interpolation: starts at 1 NA for 24h, ends at 3 NAs for 240h
    # Formula: y = y1 + (x - x1) * ((y2 - y1) / (x2 - x1))
    max_na_allowed <- round(1 + (agg_length - 24) * ((3 - 1) / (240 - 24)))
  }
  
  # 2. Apply the aggregation with the new logic
  df_detailed <- df %>%
    mutate(bin_id = hour_index %/% agg_length) %>%
    group_by(bin_id) %>%
    mutate(
      hours_in_bin = n(),
      # Count exactly how many NAs are in this specific bin
      na_count = sum(is.na(precip_rate)), 
      
      # Calculate mean conditionally
      bin_mean_precip = ifelse(
        na_count <= max_na_allowed, 
        mean(precip_rate, na.rm = TRUE), # If within tolerance, ignore NAs to get the mean
        NA_real_                         # If tolerance exceeded, the bin is invalid (NA)
      ),
      
      is_valid_bin = (!is.na(bin_mean_precip) & hours_in_bin == agg_length)
    ) %>%
    ungroup() %>%
    # Clean up the temporary na_count column so it doesn't clutter your dataframe
    select(-na_count) 
  
  # Save the dataframe to our list
  list_name <- paste0("scale_", agg_length, "Hours")
  inspection_list[[list_name]] <- df_detailed
}
 
 clean_bins_list <- list()
 for (scale_name in names(inspection_list)) {
  df_detailed <- inspection_list[[scale_name]]
  valid_bins <- df_detailed %>%
   filter(is_valid_bin) %>%
   group_by(bin_id) %>%
   summarise(
    agg_date = min(date),
    bin_mean_precip = first(bin_mean_precip),
    .groups = "drop"
   )
  clean_bins_list[[scale_name]] <- valid_bins
 }
 
 summary_stats_list <- list()
 for (scale_name in names(clean_bins_list)) {
  df_clean <- clean_bins_list[[scale_name]]
  x_all <- df_clean$bin_mean_precip
  x_pos <- x_all[x_all > 0]
  
  p_zero <- sum(x_all == 0) / length(x_all)
  mean_all <- mean(x_all, na.rm = TRUE)
  var_all  <- var(x_all, na.rm = TRUE)
  mean_pos <- ifelse(length(x_pos) > 0, mean(x_pos, na.rm = TRUE), NA)
  var_pos  <- ifelse(length(x_pos) > 1, var(x_pos, na.rm = TRUE), NA)
  
  if (length(x_all) >= 4) {
   lm_all <- lmoms(x_all)
   l1_all <- lm_all$lambdas[1]
   l2_all <- lm_all$lambdas[2]
   l3_all <- lm_all$lambdas[3]
   t2_all <- lm_all$ratios[2]
   t3_all <- lm_all$ratios[3]
   t4_all <- lm_all$ratios[4]
  } else {
   l1_all <- l2_all <- l3_all <- t1_all <- t2_all <- t3_all <- t4_all <- NA
  }
  
  if (length(x_pos) >= 4) {
   lm_pos <- lmoms(x_pos)
   l1_pos <- lm_pos$lambdas[1]
   l2_pos <- lm_pos$lambdas[2]
   l3_pos <- lm_pos$lambdas[3]
   t2_pos <- lm_pos$ratios[2]
   t3_pos <- lm_pos$ratios[3]
   t4_pos <- lm_pos$ratios[4]
  } else {
   l1_pos <- l2_pos <- l3_pos <- t1_pos <- t2_pos <- t3_pos <- t4_pos <- NA
  }
  
  scale_summary <- data.frame(
   scale = scale_name, p_zero = p_zero,
   mean_all = mean_all, var_all = var_all,
   l1_all = l1_all, l2_all = l2_all, l3_all = l3_all, t2_all = t2_all, t3_all = t3_all, t4_all = t4_all,
   mean_pos = mean_pos, var_pos = var_pos,
   l1_pos = l1_pos, l2_pos = l2_pos, l3_pos = l3_pos, t2_pos = t2_pos, t3_pos = t3_pos, t4_pos = t4_pos
  )
  summary_stats_list[[scale_name]] <- scale_summary
 }
 
 final_summary_df <- bind_rows(summary_stats_list)
 
 # --- Filter Data & Optimize ---
 final_summary_df <- final_summary_df %>%
  mutate(k_hours = as.numeric(gsub("[^0-9.]", "", scale)))
 
 fit_data <- final_summary_df %>% 
  filter(k_hours >= 24) %>% 
  arrange(k_hours)
 
 k_obs <- fit_data$k_hours
 k_star <- 24
 
 stats_to_fit <- setdiff(names(fit_data), c(
  "scale", "k_hours", "mean_all", "var_all", "l1_all", "l2_all", "l3_all",
  "mean_pos", "var_pos", "l1_pos", "l2_pos", "l3_pos"
 ))
 
 all_fits_list <- list()
 de_ctrl <- DEoptim.control(
  trace = FALSE, 
  itermax = 1000,   
  NP = 50,         
  reltol = 1e-11,    
  steptol = 500     
 )
 
 for (stat in stats_to_fit) {
  y_obs <- fit_data[[stat]]
  if(all(is.na(y_obs))) next
  
  q_k_star <- fit_data[[stat]][fit_data$k_hours == k_star][1]
  if(is.na(q_k_star)) next
  
  max_H0 <- ifelse(stat == "p_zero", 1, 1)
  
  fit_W_2p <- DEoptim(mse_W_2p, lower = c(1e-6, 1e-6), upper = c(max_H0, 1000), 
                      k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
  
  fit_L_2p <- DEoptim(mse_L_2p, lower = c(1e-6, 1e-6), upper = c(max_H0, 1000), 
                      k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
  
  res_row <- data.frame(
   Statistic = stat,
   Model = c("Weibull_2p", "PowerLaw_2p"),
   H0 = c(fit_W_2p$optim$bestmem[1], fit_L_2p$optim$bestmem[1]),
   Param_a = c(fit_W_2p$optim$bestmem[2], NA),
   Param_b = c(NA, fit_L_2p$optim$bestmem[2]),
   MSE = c(fit_W_2p$optim$bestval, fit_L_2p$optim$bestval)
  )
  
  if (stat == "p_zero") {
   fit_W_1p <- DEoptim(mse_W_1p, lower = 1e-6, upper = 10, 
                       k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
   
   fit_L_1p <- DEoptim(mse_L_1p, lower = 1e-6, upper = 10, 
                       k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
   
   res_1p <- data.frame(
    Statistic = stat,
    Model = c("Weibull_1p", "PowerLaw_1p"),
    H0 = c(1, 1), 
    Param_a = c(fit_W_1p$optim$bestmem[1], NA),
    Param_b = c(NA, fit_L_1p$optim$bestmem[1]),
    MSE = c(fit_W_1p$optim$bestval, fit_L_1p$optim$bestval)
   )
   res_row <- bind_rows(res_row, res_1p)
  }
  all_fits_list[[stat]] <- res_row
 }
 
 optimized_parameters_df <- bind_rows(all_fits_list)
 
 # --- Evaluate Validation Data ---
 validation_data <- final_summary_df %>% filter(k_hours < 24) %>% arrange(k_hours)
 k_val <- validation_data$k_hours
 empirical_24h <- final_summary_df %>% filter(k_hours == k_star)
 
 validation_results_list <- list()
 for (i in seq_len(nrow(optimized_parameters_df))) {
  current_fit <- optimized_parameters_df[i, ]
  stat <- current_fit$Statistic
  model_type <- current_fit$Model
  
  H0_val <- current_fit$H0
  a_val <- current_fit$Param_a
  b_val <- current_fit$Param_b
  
  actual_vals <- validation_data[[stat]]
  q_k_star_val <- empirical_24h[[stat]]
  
  if (all(is.na(actual_vals))) {
   current_fit$mse_validation <- NA
   validation_results_list[[i]] <- current_fit
   next
  }
  
  predicted_vals <- switch(model_type,
                           "Weibull_2p"  = H_W_2p(k_val, H0_val, k_star, q_k_star_val, a_val),
                           "PowerLaw_2p" = H_L_2p(k_val, H0_val, k_star, q_k_star_val, b_val),
                           "Weibull_1p"  = H_W_1p_pdry(k_val, k_star, q_k_star_val, a_val),
                           "PowerLaw_1p" = H_L_1p_pdry(k_val, k_star, q_k_star_val, b_val),
                           rep(NA, length(k_val)))
  
  current_fit$mse_validation <- calculate_mse(actual_vals, predicted_vals)
  validation_results_list[[i]] <- current_fit
 }
 
 final_evaluated_models <- bind_rows(validation_results_list)
 
 best_models <- final_evaluated_models %>%
  group_by(Statistic) %>%
  slice_min(order_by = mse_validation, n = 1) %>%
  ungroup()
 
 # --- Generate All Predictions ---
 k_empirical <- sort(unique(final_summary_df$k_hours))
 k_predict <- c(0.25, k_empirical) 
 predicted_list <- list()
 
 for (stat in unique(optimized_parameters_df$Statistic)) {
  q_k_star_val <- final_summary_df[[stat]][final_summary_df$k_hours == k_star]
  actual_vals <- final_summary_df[[stat]][match(k_predict, final_summary_df$k_hours)]
  params <- optimized_parameters_df %>% filter(Statistic == stat)
  
  df_stat <- data.frame(Scale_k = k_predict, Statistic = stat, Actual = actual_vals)
  
  for (i in seq_len(nrow(params))) {
   m_name <- params$Model[i]
   H0_val <- params$H0[i]
   a_val <- params$Param_a[i]
   b_val <- params$Param_b[i]
   
   predicted_vals <- switch(m_name,
                            "Weibull_2p"  = H_W_2p(k_predict, H0_val, k_star, q_k_star_val, a_val),
                            "PowerLaw_2p" = H_L_2p(k_predict, H0_val, k_star, q_k_star_val, b_val),
                            "Weibull_1p"  = H_W_1p_pdry(k_predict, k_star, q_k_star_val, a_val),
                            "PowerLaw_1p" = H_L_1p_pdry(k_predict, k_star, q_k_star_val, b_val),
                            rep(NA, length(k_predict)))
   df_stat[[m_name]] <- predicted_vals
  }
  predicted_list[[stat]] <- df_stat
 }
 
 all_predictions_df <- bind_rows(predicted_list)
 
 # =====================================================================
 # EXPORT STATION RESULTS
 # =====================================================================
 station_folder <- file.path(individual_dir, station_name)
 dir.create(station_folder, recursive = TRUE, showWarnings = FALSE)
 
 saveRDS(final_summary_df, file.path(station_folder, "final_summary.rds"))
 saveRDS(optimized_parameters_df, file.path(station_folder, "optimized_parameters.rds"))
 saveRDS(final_evaluated_models, file.path(station_folder, "final_evaluated_models.rds"))
 saveRDS(best_models, file.path(station_folder, "best_models.rds"))
 saveRDS(all_predictions_df, file.path(station_folder, "all_predictions.rds"))
 
 cat("Successfully exported results for:", station_name, "\n")
 
} 


subset_data <- all_predictions[all_predictions$Statistic == "p_zero" & 
                                all_predictions$Scale_k >= 1 & 
                                all_predictions$Scale_k <= 23, ]


mse_weibull_2p <- mean((subset_data$Actual - subset_data$Weibull_2p)^2, na.rm = TRUE)


mse_weibull_1p <- mean((subset_data$Actual - subset_data$Weibull_1p)^2, na.rm = TRUE)










