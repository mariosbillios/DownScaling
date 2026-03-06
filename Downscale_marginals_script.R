
library(dplyr)
library(lubridate)
library(lmomco)
library(DEoptim)
library(tools)
library(patchwork)
library(ggplot2)
library(tidyr)



# Define expected models and data globally
global_expected_models <- c("Weibull_1p", "PowerLaw_1p", "Weibull_2p", "PowerLaw_2p", "Weibull_3p", "PowerLaw_3p")
global_expected_data   <- c("Calib. data", "Valid. data")

# Define the plotting function
plot_statistic_dynamic <- function(stat_name, y_label = stat_name, show_legend = FALSE, fixed_y = TRUE) {
  
  # Prepare Empirical Data points (Using the _df variables from your loop)
  emp_data <- final_summary_df %>%
    select(k_hours, value = !!sym(stat_name)) %>%
    filter(!is.na(value)) %>%
    mutate(Data_Type = ifelse(k_hours >= 24, "Calib. data", "Valid. data"))
  
  max_k <- max(emp_data$k_hours, na.rm = TRUE)
  q_k_star_val <- emp_data$value[emp_data$k_hours == 24]
  
  if(length(q_k_star_val) == 0) {
    warning(paste("No 24h anchor found for", stat_name, "- skipping plot."))
    return(NULL) 
  }
  
  k_smooth <- exp(seq(log(1), log(max_k), length.out = 500))
  params <- optimized_parameters_df %>% filter(Statistic == stat_name)
  
  smooth_lines_list <- list()
  for (i in 1:nrow(params)) {
    m_name <- params$Model[i]
    H0_val <- params$H0[i]
    a_val  <- params$Param_a[i]
    b_val  <- params$Param_b[i]
    
    pred_smooth <- switch(m_name,
                          "Weibull_3p"  = H_W_3p(k_smooth, H0_val, a_val, b_val),
                          "PowerLaw_3p" = H_L_3p(k_smooth, H0_val, a_val, b_val),
                          "Weibull_2p"  = H_W_2p(k_smooth, H0_val, 24, q_k_star_val, a_val),
                          "PowerLaw_2p" = H_L_2p(k_smooth, H0_val, 24, q_k_star_val, b_val),
                          "Weibull_1p"  = H_W_1p_pdry(k_smooth, 24, q_k_star_val, a_val),
                          "PowerLaw_1p" = H_L_1p_pdry(k_smooth, 24, q_k_star_val, b_val),
                          rep(NA, length(k_smooth))
    )
    smooth_lines_list[[i]] <- data.frame(Scale_k = k_smooth, Model = m_name, value = pred_smooth)
  }
  
  smooth_data <- bind_rows(smooth_lines_list) %>% filter(!is.na(value))
  
  line_colors <- c("Weibull_1p" = "blue", "PowerLaw_1p" = "red",
                   "Weibull_2p" = "blue", "PowerLaw_2p" = "red",
                   "Weibull_3p" = "darkblue", "PowerLaw_3p" = "darkred")
  
  line_types <- c("Weibull_1p" = "solid", "PowerLaw_1p" = "solid",
                  "Weibull_2p" = "dashed", "PowerLaw_2p" = "dashed",
                  "Weibull_3p" = "dotted", "PowerLaw_3p" = "dotted")
  
  custom_breaks <- c(1, 6, 12, 24, 120, 240)
  custom_breaks <- custom_breaks[custom_breaks <= max_k] 
  
  fill_guide  <- if(show_legend) guide_legend(order = 1, ncol = 1) else "none"
  line_guide  <- if(show_legend) guide_legend(order = 2, nrow = 2) else "none"
  
  p <- ggplot() +
    geom_line(data = smooth_data, aes(x = Scale_k, y = value, color = Model, linetype = Model), linewidth = 1) +
    geom_point(data = emp_data, aes(x = k_hours, y = value, fill = Data_Type), 
               shape = 21, color = "white", size = 3, stroke = 0.5) +
    geom_vline(xintercept = 24, linetype = "dashed", color = "darkgray", linewidth = 0.8) +
    scale_x_log10(breaks = custom_breaks, labels = custom_breaks) +
    scale_fill_manual(name = "", values = c("Calib. data" = "black", "Valid. data" = "orange"),
                      limits = global_expected_data, drop = FALSE, guide = fill_guide) +
    scale_color_manual(name = "", values = line_colors, limits = global_expected_models, 
                       drop = FALSE, guide = line_guide) +
    scale_linetype_manual(name = "", values = line_types, limits = global_expected_models, 
                          drop = FALSE, guide = line_guide) +
    labs(x = "Temporal scale, k [h]", y = y_label) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(linetype = "dashed", color = "lightgray"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
  
  if (fixed_y) {
    p <- p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
  }
  
  return(p)
}

plot_qq_dynamic <- function(stat_name, y_label = stat_name) {
  
  # 1. Filter predictions for the specific statistic
  df_stat <- all_predictions_df %>% 
    filter(Statistic == stat_name) %>%
    filter(!is.na(Actual)) %>% # Remove extrapolated k=0.25 if it has no empirical data
    mutate(Data_Type = ifelse(Scale_k >= 24, "Calib. data", "Valid. data"))
  
  if(nrow(df_stat) == 0) return(NULL)
  
  # 2. Pivot data to long format so we can facet by Model
  # We gather all columns EXCEPT Scale_k, Statistic, Actual, and Data_Type
  model_cols <- setdiff(names(df_stat), c("Scale_k", "Statistic", "Actual", "Data_Type"))
  
  df_long <- df_stat %>%
    pivot_longer(cols = all_of(model_cols), names_to = "Model", values_to = "Predicted") %>%
    filter(!is.na(Predicted))
  
  if(nrow(df_long) == 0) return(NULL)
  
  # 3. Find global min/max so the X and Y axes perfectly match (true 1:1 plot)
  min_val <- min(c(df_long$Actual, df_long$Predicted), na.rm = TRUE)
  max_val <- max(c(df_long$Actual, df_long$Predicted), na.rm = TRUE)
  buffer <- (max_val - min_val) * 0.05
  if(buffer == 0) buffer <- 0.05
  
  # 4. Generate the Plot
  p <- ggplot(df_long, aes(x = Actual, y = Predicted)) +
    # The 1:1 perfect agreement line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgray", linewidth = 0.8) +
    
    # The actual data points
    geom_point(aes(fill = Data_Type), shape = 21, color = "white", size = 2.5, stroke = 0.5) +
    
    # Create separate panels for each model (e.g., Weibull_1p, Weibull_2p, etc.)
    facet_wrap(~ Model) + 
    
    # Formatting
    scale_fill_manual(name = "", values = c("Calib. data" = "black", "Valid. data" = "orange")) +
    scale_x_continuous(limits = c(min_val - buffer, max_val + buffer)) +
    scale_y_continuous(limits = c(min_val - buffer, max_val + buffer)) +
    coord_fixed(ratio = 1) + # Forces the plot to be perfectly square
    
    labs(
      title = bquote("Empirical vs Predicted:" ~ .(y_label[[1]])),
      x = bquote("Empirical (Actual)" ~ .(y_label[[1]])),
      y = bquote("Predicted" ~ .(y_label[[1]]))
    )  +
    theme_bw() +
    theme(
      panel.grid.major = element_line(linetype = "dotted", color = "lightgray"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "whitesmoke"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}












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
station_files<-station_files[64:65]

for (current_station_file in station_files) {
 
 
 # Extract current station name for export and logging
 station_name <- tools::file_path_sans_ext(basename(current_station_file))
 
 # Wrap the entire processing logic in tryCatch
 tryCatch({

 df <- read.table(current_station_file, skip = 21, col.names = "precip")
 df$precip[df$precip == -999] <- NA
 
 header_lines <- readLines(current_station_file, n = 21)
 start_date_str <- sub("Start datetime: ", "", grep("^Start datetime:", header_lines, value = TRUE))
 end_date_str <- sub("End datetime: ", "", grep("^End datetime:", header_lines, value = TRUE))
 
 start_ts <- ymd_h(start_date_str)
 end_ts <- ymd_h(end_date_str)
 df$date <- seq(start_ts, end_ts, by = "hour")
 
 
 
 # =====================================================================
 # 1. INITIAL METADATA CHECK (Place this right after defining start_ts / end_ts)
 # =====================================================================
 
 meta_records_str <- grep("^Number of records:", header_lines, value = TRUE)
 meta_missing_str <- grep("^Percent missing data:", header_lines, value = TRUE)
 
 metadata_records <- as.numeric(sub("Number of records: ", "", meta_records_str))
 metadata_missing_pct <- as.numeric(sub("Percent missing data: ", "", meta_missing_str))
 
 # Calculate Theoretical Data (Total hours between start and end)
 theoretical_records <- as.numeric(difftime(end_ts, start_ts, units = "hours")) + 1
 
 explicit_nas <- sum(is.na(df$precip)) 
 missing_rows <- max(0, theoretical_records - nrow(df))
 total_missing_data_points <- explicit_nas + missing_rows
 calculated_missing_pct <- (total_missing_data_points / theoretical_records) * 100
 
 
 
 
 Base_time_chunks_per_hour <- 4
 df$precip_rate <- df$precip 
 df$hour_index <- 0:(nrow(df) - 1)
 
 
 
 all_k<-c(1:12,24,24*2,24*3,24*4,24*5,24*9,24*10)
 
 inspection_list <- list()
 na_thresholds_list <- list() # To store thresholds for the .rds metadata
 
 for (i in seq_along(all_k)) {
   
   agg_length <- all_k[i] 
   scale_name <- paste0("scale_", agg_length, "Hours")
   
   # 1. Apply the conditional threshold rule
   if (agg_length < 6) {
     max_na_allowed <- 0  # Strict zero-tolerance strictly under 6h
     
   } else if (agg_length >= 6 & agg_length <= 12) {
     max_na_allowed <- 1  # Rule for 6h to 12h (easy to change to 1 later if needed)
     
   } else {
     # Fast rule for larger scales: max 1/3 of the values missing, rounded up
     max_na_allowed <- ceiling(agg_length / 5)
     
     # Safety check: never allow a completely empty bin to pass
     max_na_allowed <- min(max_na_allowed, agg_length - 1) 
   }
   
   # 2. Store the calculated threshold for your metadata export
   na_thresholds_list[[scale_name]] <- max_na_allowed
   
   # 3. Apply the aggregation
   df_detailed <- df %>%
     mutate(bin_id = hour_index %/% agg_length) %>%
     group_by(bin_id) %>%
     mutate(
       hours_in_bin = n(),
       na_count = sum(is.na(precip_rate)), 
       
       # Calculate mean conditionally
       bin_mean_precip = ifelse(
         na_count <= max_na_allowed, 
         mean(precip_rate, na.rm = TRUE),
         NA_real_                        
       ),
       
       is_valid_bin = (!is.na(bin_mean_precip) & hours_in_bin == agg_length)
     ) %>%
     ungroup() %>%
     select(-na_count) 
   
   # Save the dataframe to our list
   inspection_list[[scale_name]] <- df_detailed
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
 
 
 
 
 # Initialize an empty list to store the scale data frames
 scale_na_summary_list <- list()
 
 for (scale_name in names(clean_bins_list)) {
   # Extract the numerical scale from the name (e.g., "scale_24Hours" -> 24)
   agg_length <- as.numeric(gsub("[^0-9.]", "", scale_name))
   
   # Calculate how many bins there WOULD be if the gauge never broke down
   theor_bins <- theoretical_records %/% agg_length
   
   # Count how many bins actually survived the NA tolerance check
   valid_bins_count <- nrow(clean_bins_list[[scale_name]])
   
   # Calculate losses
   dropped_bins <- theor_bins - valid_bins_count
   dropped_pct <- (dropped_bins / theor_bins) * 100
   
   # Store the results as a data frame row rather than text
   scale_na_summary_list[[scale_name]] <- data.frame(
     Scale_k_hours = agg_length,
     Theoretical_Bins = theor_bins,
     Valid_Bins = valid_bins_count,
     Dropped_Bins = dropped_bins,
     Dropped_Pct = dropped_pct
   )
 }
 
 
 
 scale_na_summary_list <- list()
 
 for (scale_name in names(clean_bins_list)) {
   agg_length <- as.numeric(gsub("[^0-9.]", "", scale_name))
   theor_bins <- theoretical_records %/% agg_length
   valid_bins_count <- nrow(clean_bins_list[[scale_name]])
   dropped_bins <- theor_bins - valid_bins_count
   dropped_pct <- (dropped_bins / theor_bins) * 100
   
   # Retrieve the threshold we stored earlier
   current_max_na <- na_thresholds_list[[scale_name]]
   
   # Store results, now including the dynamic threshold
   scale_na_summary_list[[scale_name]] <- data.frame(
     Scale_k_hours = agg_length,
     Max_NA_Allowed = current_max_na,    # NEW COLUMN ADDED HERE
     Theoretical_Bins = theor_bins,
     Valid_Bins = valid_bins_count,
     Dropped_Bins = dropped_bins,
     Dropped_Pct = dropped_pct
   )
 }
 
 # Combine all scales into your final data frame
 scale_na_summary_df <- bind_rows(scale_na_summary_list)
 
 
 
 
 
 
 # Create a comprehensive metadata list object
 station_metadata_obj <- list(
   General_Info = list(
     Station = station_name,
     Start_Datetime = start_ts,
     End_Datetime = end_ts,
     Theoretical_Records = theoretical_records,
     Metadata_Records = metadata_records,
     Actual_Rows = nrow(df)
   ),
   Missing_Data_Summary = list(
     Explicit_NAs = explicit_nas,
     Missing_Rows = missing_rows,
     Total_Missing_Points = total_missing_data_points,
     Calculated_Missing_Pct = calculated_missing_pct,
     Metadata_Missing_Pct = metadata_missing_pct
   ),
   Scale_Data_Loss = scale_na_summary_df
 )
 
 
 
 
 
 
 
 
 
 
 
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
 saveRDS(station_metadata_obj, file.path(station_folder,"extra_metadata.rds"))
 
 
 
 # 1. Fixed Y-Axis Plots (0 to 1)
 plot_p_zero_fixed <- plot_statistic_dynamic("p_zero", y_label = expression(p[0]^{(k)}), show_legend = TRUE, fixed_y = TRUE)
 plot_t2_fixed     <- plot_statistic_dynamic("t2_pos", y_label = expression(t[2]^{(k)}), fixed_y = TRUE)
 plot_t3_fixed     <- plot_statistic_dynamic("t3_pos", y_label = expression(t[3]^{(k)}), fixed_y = TRUE)
 plot_t4_fixed     <- plot_statistic_dynamic("t4_pos", y_label = expression(t[4]^{(k)}), fixed_y = TRUE)
 
 
 
 # Check if plots generated successfully before combining
 if (!is.null(plot_p_zero_fixed)) {
   combined_plot_fixed <- (plot_p_zero_fixed | plot_t2_fixed) / (plot_t3_fixed | plot_t4_fixed) + 
     plot_layout(guides = "collect") & 
     theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "left", legend.spacing.x = unit(0.5, "cm"))
   
   # Export Fixed Plot
   ggsave(file.path(station_folder, "scaling_plots_fixed.png"), combined_plot_fixed, width = 10, height = 8, bg = "white")
 }
 
 # 2. Dynamic Y-Axis Plots (Unconstrained)
 plot_p_zero_dyn <- plot_statistic_dynamic("p_zero", y_label = expression(p[0]^{(k)}), show_legend = TRUE, fixed_y = FALSE)
 plot_t2_dyn     <- plot_statistic_dynamic("t2_pos", y_label = expression(t[2]^{(k)}), fixed_y = FALSE)
 plot_t3_dyn     <- plot_statistic_dynamic("t3_pos", y_label = expression(t[3]^{(k)}), fixed_y = FALSE)
 plot_t4_dyn     <- plot_statistic_dynamic("t4_pos", y_label = expression(t[4]^{(k)}), fixed_y = FALSE)
 
 if (!is.null(plot_p_zero_dyn)) {
   combined_plot_dynamic <- (plot_p_zero_dyn | plot_t2_dyn) / (plot_t3_dyn | plot_t4_dyn) + 
     plot_layout(guides = "collect") & 
     theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "left", legend.spacing.x = unit(0.5, "cm"))
   
   # Export Dynamic Plot
   ggsave(file.path(station_folder, "scaling_plots_dynamic.png"), combined_plot_dynamic, width = 10, height = 8, bg = "white")
 }
 
 
 
 # =====================================================================
 # 4. GENERATE EMPIRICAL VS PREDICTED (Q-Q) SCATTER PLOTS
 # =====================================================================
 
 qq_p_zero <- plot_qq_dynamic("p_zero", y_label = expression(p[0]^{(k)}))
 qq_t2     <- plot_qq_dynamic("t2_pos", y_label = expression(t[2]^{(k)}))
 qq_t3     <- plot_qq_dynamic("t3_pos", y_label = expression(t[3]^{(k)}))
 qq_t4     <- plot_qq_dynamic("t4_pos", y_label = expression(t[4]^{(k)}))
 
 # Combine them into one big grid using patchwork
 # We use wrap_plots to handle any NULLs safely if a stat failed to compute
 qq_list <- list(qq_p_zero, qq_t2, qq_t3, qq_t4)
 qq_list <- qq_list[!sapply(qq_list, is.null)] # Remove empty plots
 
 if(length(qq_list) > 0) {
   combined_qq_plots <- wrap_plots(qq_list, ncol = 2) + 
     plot_layout(guides = "collect") & 
     theme(legend.position = "bottom")
   
   # Save the huge validation grid
   ggsave(file.path(station_folder, "qq_validation_plots.png"), 
          combined_qq_plots, width = 12, height = 10, bg = "white")
 }
 
 
 
 
 
 
 
 cat("Successfully exported results and plots for:", station_name, "\n")
 
 }, error = function(e) {
  
  # This block ONLY executes if an error occurs in the code above
  # It prints a warning to the console and gracefully continues to the next file
  message("--------------------------------------------------")
  message("FAILED: Station ", station_name, " encountered an error.")
  message("Error detail: ", conditionMessage(e))
  message("Skipping to the next station...")
  message("--------------------------------------------------")
  
 }) 
 
} # End of the main station loop









