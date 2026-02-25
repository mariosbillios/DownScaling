Lm_Pd_downscaleFun <- function(file_path, intervals_per_day = 96) {
 # --- 1. Load Data ---
 header_lines <- readLines(file_path, n = 16)
 start_date_str <- sub("Start datetime: ", "", grep("^Start datetime:", header_lines, value = TRUE))
 end_date_str <- sub("End datetime: ", "", grep("^End datetime:", header_lines, value = TRUE))
 
 df <- read.table(file_path, skip = 16, col.names = "precip")
 df$precip[df$precip == -999] <- NA
 df$date <- seq(ymd(start_date_str), ymd(end_date_str), by = "day")
 
 # --- 2. Aggregation Logic ---
 df$precip_rate <- df$precip / intervals_per_day
 df$day_index <- as.numeric(df$date - min(df$date))
 
 aggregation_days <- 1:20
 inspection_list <- list()
 
 for (agg_length in aggregation_days) {
  valid_bins <- df %>%
   mutate(bin_id = day_index %/% agg_length) %>%
   group_by(bin_id) %>%
   summarise(days_in_bin = n(), bin_mean_precip = mean(precip_rate, na.rm = FALSE), .groups = 'drop') %>%
   filter(!is.na(bin_mean_precip) & days_in_bin == agg_length)
  
  get_stats <- function(data) {
   if(length(data) > 3) {
    m <- lmoms(data)$ratios
    return(c(m[2], m[3]))
   } 
   return(c(NA, NA))
  }
  
  l_all <- get_stats(valid_bins$bin_mean_precip)
  l_pos <- get_stats(valid_bins$bin_mean_precip[valid_bins$bin_mean_precip > 0])
  
  inspection_list[[agg_length]] <- data.frame(
   station = basename(file_path), scale_days = agg_length, scale_k = agg_length * intervals_per_day,
   prob_dry = mean(valid_bins$bin_mean_precip == 0),
   L_var_all = l_all[1], L_skew_all = l_all[2], L_var_pos = l_pos[1], L_skew_pos = l_pos[2]
  )
 }
 inspection_df <- do.call(rbind, inspection_list)
 
 # --- 3. Optimization Logic ---
 bounded_model <- function(k, p) p[1] ^ ((1 + (p[2]^(-1/p[3]) - 1) * (k - 1))^p[3])
 
 run_full_opt <- function(empirical_stat, k_vals, stat_label) {
  valid <- !is.na(empirical_stat)
  if(sum(valid) < 3) return(rep(NA, 3))
  
  obj <- function(p) sum((1 - (bounded_model(k_vals[valid], p) / empirical_stat[valid]))^2, na.rm = TRUE)
  fit <- DEoptim(fn = obj, lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1), 
                 control = DEoptim.control(itermax = 500, trace = FALSE))
  
  res <- fit$optim$bestmem
  names(res) <- paste0(stat_label, c("_m_basic", "_xi", "_eta"))
  return(res)
 }
 
 # Run optimization for all 5 stats
 p1 <- run_full_opt(inspection_df$prob_dry, inspection_df$scale_k, "p_dry")
 p2 <- run_full_opt(inspection_df$L_var_all, inspection_df$scale_k, "L_var_all")
 p3 <- run_full_opt(inspection_df$L_skew_all, inspection_df$scale_k, "L_skew_all")
 p4 <- run_full_opt(inspection_df$L_var_pos, inspection_df$scale_k, "L_var_pos")
 p5 <- run_full_opt(inspection_df$L_skew_pos, inspection_df$scale_k, "L_skew_pos")
 
 # Store Parameters
 param_df <- data.frame(Station = basename(file_path), t(c(p1, p2, p3, p4, p5)))
 
 # Store 15-min Estimates (m_basic values)
 estimates_df <- data.frame(
  Station = basename(file_path),
  p_dry_15min = p1[1], L_var_all_15min = p2[1], L_skew_all_15min = p3[1],
  L_var_pos_15min = p4[1], L_skew_pos_15min = p5[1]
 )
 
 return(list(estimates = estimates_df, inspection = inspection_df, params = param_df))
}