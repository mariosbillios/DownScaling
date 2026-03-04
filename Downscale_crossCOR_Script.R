library(dplyr)
library(lubridate)
library(lmomco)
library(DEoptim)
library(tools)

# =====================================================================
# 1. INITIAL SETUP & PATHS
# =====================================================================
gdrive_path <- paste(LETTERS[file.exists(paste0(LETTERS, ":/My Drive"))], ":/My Drive", sep = "")
project_path <- file.path(gdrive_path, "Academic_git/DownScaling")
generalData_path <- file.path(gdrive_path, "General_Data")

data_dir <- file.path(generalData_path, "GSDR/QC_d data - Germany")

folder_name <- basename(data_dir)
folder_name <- sub(".*- ", "", folder_name)
export_folder_name <- paste("Cross_Corr_Downs", folder_name)
individual_dir <- file.path(project_path, export_folder_name)

if (!dir.exists(individual_dir)) {
 dir.create(individual_dir, recursive = TRUE)
 message("Created new directory: ", individual_dir)
}

station_files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

 station_files <- station_files[1:3] 

 
 
 
 
 # =====================================================================
 # CROSS-CORRELATION FUNCTIONS (Monotonically Increasing)
 # =====================================================================
 
 # -- 2-Parameter Functions --
 # Equation (10): Weibull-H>
 H_W_2p_cross <- function(k, H_inf, k_star, q_k_star, a) {
  exponent <- (k / k_star)^a
  base_fraction <- 1 - (q_k_star / H_inf)
  base_fraction <- ifelse(base_fraction <= 0, 1e-6, base_fraction) 
  H_inf * (1 - (base_fraction ^ exponent))
 }
 
 # Equation (11): Lomax-H>
 H_L_2p_cross <- function(k, H_inf, k_star, q_k_star, b) {
  ln_numerator <- log(1 + b * k)
  ln_denominator <- log(1 + b * k_star)
  exponent <- ln_numerator / ln_denominator
  base_fraction <- 1 - (q_k_star / H_inf)
  base_fraction <- ifelse(base_fraction <= 0, 1e-6, base_fraction)
  H_inf * (1 - (base_fraction ^ exponent))
 }
 
 # -- 1-Parameter Functions (Assuming H_inf = 1) --
 H_W_1p_cross <- function(k, k_star, q_k_star, a) {
  H_inf <- 1 
  exponent <- (k / k_star)^a
  base_fraction <- 1 - (q_k_star / H_inf)
  base_fraction <- ifelse(base_fraction <= 0, 1e-6, base_fraction)
  H_inf * (1 - (base_fraction ^ exponent))
 }
 
 H_L_1p_cross <- function(k, k_star, q_k_star, b) {
  H_inf <- 1 
  ln_numerator <- log(1 + b * k)
  ln_denominator <- log(1 + b * k_star)
  exponent <- ln_numerator / ln_denominator
  base_fraction <- 1 - (q_k_star / H_inf)
  base_fraction <- ifelse(base_fraction <= 0, 1e-6, base_fraction)
  H_inf * (1 - (base_fraction ^ exponent))
 }
 
 # -- MSE Functions for Optimization --
 mse_W_2p_cross <- function(par, k, y, k_star, q_k_star) mean((y - H_W_2p_cross(k, par[1], k_star, q_k_star, par[2]))^2, na.rm=TRUE)
 mse_L_2p_cross <- function(par, k, y, k_star, q_k_star) mean((y - H_L_2p_cross(k, par[1], k_star, q_k_star, par[2]))^2, na.rm=TRUE)
 
 mse_W_1p_cross <- function(par, k, y, k_star, q_k_star) mean((y - H_W_1p_cross(k, k_star, q_k_star, par[1]))^2, na.rm=TRUE)
 mse_L_1p_cross <- function(par, k, y, k_star, q_k_star) mean((y - H_L_1p_cross(k, k_star, q_k_star, par[1]))^2, na.rm=TRUE)
 
 
 
 
 
 
 
 
 
# =====================================================================
# 2. READ, FILTER & AGGREGATE DATA (Builds 'all_clean_bins_lists')
# =====================================================================
raw_data_list <- list()
all_clean_bins_lists <- list()

cat("\n--- Phase 1: Reading Raw Data ---\n")
for (file_path in station_files) {
 station_name <- tools::file_path_sans_ext(basename(file_path))
 cat("Reading:", station_name, "\n")
 
 df <- read.table(file_path, skip = 21, col.names = "precip")
 df$precip[df$precip == -999] <- NA
 
 header_lines <- readLines(file_path, n = 21)
 start_date_str <- sub("Start datetime: ", "", grep("^Start datetime:", header_lines, value = TRUE))
 end_date_str <- sub("End datetime: ", "", grep("^End datetime:", header_lines, value = TRUE))
 
 start_ts <- ymd_h(start_date_str)
 end_ts <- ymd_h(end_date_str)
 df$date <- seq(start_ts, end_ts, by = "hour")
 
 raw_data_list[[station_name]] <- df
}

# Find Common Date Window
common_start <- max(do.call(c, lapply(raw_data_list, function(df) min(df$date))))
common_end <- min(do.call(c, lapply(raw_data_list, function(df) max(df$date))))

Base_time_chunks_per_hour <- 4
all_k <- c(c(1,2,4,6,8,10,12) * Base_time_chunks_per_hour, seq(1:10) * 24 * Base_time_chunks_per_hour)

cat("\n--- Phase 2: Aggregating to Time Scales ---\n")
for (station_name in names(raw_data_list)) {
 df <- raw_data_list[[station_name]] %>% filter(date >= common_start & date <= common_end)
 df$precip_rate <- df$precip / Base_time_chunks_per_hour
 df$hour_index <- 0:(nrow(df) - 1)
 
 station_clean_bins_list <- list()
 
 for (i in seq_along(all_k)) {
  agg_length <- all_k[i] / Base_time_chunks_per_hour
  list_name <- paste0("scale_", agg_length, "Hours")
  
  valid_bins <- df %>%
   mutate(bin_id = hour_index %/% agg_length) %>%
   group_by(bin_id) %>%
   summarise(
    hours_in_bin = n(),
    agg_date = min(date),
    bin_mean_precip = mean(precip_rate, na.rm = FALSE),
    .groups = "drop"
   ) %>%
   filter(!is.na(bin_mean_precip) & hours_in_bin == agg_length)
  
  station_clean_bins_list[[list_name]] <- valid_bins
 }
 all_clean_bins_lists[[station_name]] <- station_clean_bins_list
}

# =====================================================================
# 3. WRAPPER FUNCTION: PAIRWISE DOWNSCALING & EXPORT
# =====================================================================
process_station_pair <- function(stat1, stat2, clean_bins_list, export_base_dir) {
 
 cat("Processing Pair: ", stat1, " & ", stat2, "\n")
 scale_names <- names(clean_bins_list[[stat1]])
 
 correlation_results <- data.frame(Scale = character(), Numeric_Hours = numeric(), Valid_Bin_Pairs = integer(), Correlation_Lag_0 = numeric(), stringsAsFactors = FALSE)
 
 for (scale_name in scale_names) {
  df1 <- clean_bins_list[[stat1]][[scale_name]]
  df2 <- clean_bins_list[[stat2]][[scale_name]]
  
  joined_df <- inner_join(df1, df2, by = "bin_id", suffix = c("_1", "_2"))
  n_pairs <- nrow(joined_df)
  corr_val <- if(n_pairs > 1) cor(joined_df$bin_mean_precip_1, joined_df$bin_mean_precip_2, method = "pearson") else NA
  
  num_hours <- as.numeric(gsub("[^0-9]", "", scale_name))
  correlation_results <- rbind(correlation_results, data.frame(Scale = scale_name, Numeric_Hours = num_hours, Valid_Bin_Pairs = n_pairs, Correlation_Lag_0 = corr_val))
 }
 
 correlation_results <- correlation_results %>% arrange(Numeric_Hours)
 
 fit_data <- data.frame(k_hours = correlation_results$Numeric_Hours, cross_corr = correlation_results$Correlation_Lag_0) %>% arrange(k_hours)
 
 k_star <- 24
 training_data <- fit_data %>% filter(k_hours >= k_star)
 validation_data <- fit_data %>% filter(k_hours < k_star & k_hours >= 1)
 
 k_obs <- training_data$k_hours
 y_obs <- training_data$cross_corr
 q_k_star <- training_data$cross_corr[training_data$k_hours == k_star][1]
 
 if(is.na(q_k_star)) {
  cat("  -> Skipped: Missing 24h correlation.\n")
  return(NULL)
 }
 
 de_ctrl <- DEoptim.control(trace = FALSE, itermax = 5000, NP = 500, reltol = 1e-11, steptol = 500)
 
 fit_W_2p <- DEoptim(mse_W_2p_cross, lower = c(1e-5, 1e-6), upper = c(1.0, 1000), k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
 fit_L_2p <- DEoptim(mse_L_2p_cross, lower = c(1e-5, 1e-6), upper = c(1.0, 1000), k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
 fit_W_1p <- DEoptim(mse_W_1p_cross, lower = 1e-6, upper = 1000, k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
 fit_L_1p <- DEoptim(mse_L_1p_cross, lower = 1e-6, upper = 1000, k = k_obs, y = y_obs, k_star = k_star, q_k_star = q_k_star, control = de_ctrl)
 
 optimized_params <- data.frame(
  Statistic = "cross_corr_lag0",
  Model = c("Weibull_2p", "PowerLaw_2p", "Weibull_1p", "PowerLaw_1p"),
  H_inf = c(fit_W_2p$optim$bestmem[1], fit_L_2p$optim$bestmem[1], 1, 1),
  Param_shape = c(fit_W_2p$optim$bestmem[2], fit_L_2p$optim$bestmem[2], fit_W_1p$optim$bestmem[1], fit_L_1p$optim$bestmem[1]),
  Training_MSE = c(fit_W_2p$optim$bestval, fit_L_2p$optim$bestval, fit_W_1p$optim$bestval, fit_L_1p$optim$bestval)
 )
 
 k_val <- validation_data$k_hours
 y_val_actual <- validation_data$cross_corr
 
 final_evaluated_models <- optimized_params
 final_evaluated_models$Validation_MSE <- c(
  mean((y_val_actual - H_W_2p_cross(k_val, optimized_params$H_inf[1], k_star, q_k_star, optimized_params$Param_shape[1]))^2, na.rm = TRUE),
  mean((y_val_actual - H_L_2p_cross(k_val, optimized_params$H_inf[2], k_star, q_k_star, optimized_params$Param_shape[2]))^2, na.rm = TRUE),
  mean((y_val_actual - H_W_1p_cross(k_val, k_star, q_k_star, optimized_params$Param_shape[3]))^2, na.rm = TRUE),
  mean((y_val_actual - H_L_1p_cross(k_val, k_star, q_k_star, optimized_params$Param_shape[4]))^2, na.rm = TRUE)
 )
 
 best_model_info <- final_evaluated_models %>% slice_min(order_by = Validation_MSE, n = 1)
 
 k_all <- c(0.25, validation_data$k_hours, training_data$k_hours)
 all_predictions <- data.frame(k_hours = k_all, Actual_Corr = c(NA, validation_data$cross_corr, training_data$cross_corr))
 
 all_predictions$Weibull_2p <- H_W_2p_cross(k_all, optimized_params$H_inf[1], k_star, q_k_star, optimized_params$Param_shape[1])
 all_predictions$PowerLaw_2p <- H_L_2p_cross(k_all, optimized_params$H_inf[2], k_star, q_k_star, optimized_params$Param_shape[2])
 all_predictions$Weibull_1p <- H_W_1p_cross(k_all, k_star, q_k_star, optimized_params$Param_shape[3])
 all_predictions$PowerLaw_1p <- H_L_1p_cross(k_all, k_star, q_k_star, optimized_params$Param_shape[4])
 
 # Export handling
 pair_folder_name <- paste0(stat1, "_vs_", stat2)
 pair_export_dir <- file.path(export_base_dir, pair_folder_name)
 dir.create(pair_export_dir, recursive = TRUE, showWarnings = FALSE)
 
 saveRDS(correlation_results, file.path(pair_export_dir, "final_summary.rds"))
 saveRDS(optimized_params, file.path(pair_export_dir, "optimized_parameters.rds"))
 saveRDS(final_evaluated_models, file.path(pair_export_dir, "final_evaluated_models.rds"))
 saveRDS(best_model_info, file.path(pair_export_dir, "best_models.rds"))
 saveRDS(all_predictions, file.path(pair_export_dir, "all_predictions.rds"))
 
 return(TRUE)
}

# =====================================================================
# 4. EXECUTION LOOP: GENERATE PAIRS & RUN DOWNSCALING
# =====================================================================
cat("\n--- Phase 3: Executing Downscaling for Station Pairs ---\n")
station_names <- names(all_clean_bins_lists)
station_pairs <- combn(station_names, 2, simplify = FALSE)

cat("Found", length(station_pairs), "unique station pairs to process.\n\n")

for (pair in station_pairs) {
 process_station_pair(
  stat1 = pair[1], 
  stat2 = pair[2], 
  clean_bins_list = all_clean_bins_lists, 
  export_base_dir = individual_dir 
 )
}

cat("\nDone! All pairs processed and exported.\n")