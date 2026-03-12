# Load required libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(lubridate)
library(DEoptim)
library(lmomco)
library(tidyr)

# =====================================================================
# 1. DEFINE PATHS & DIRECTORIES (Your Setup)
# =====================================================================
gdrive_path   <- "G:/My Drive"
project_path  <- file.path(gdrive_path, "Academic_git/DownScaling")
data_dir_name <- "QC_d data - Germany" 

metadata_dir <- file.path(project_path, "Marginals", data_dir_name, "Metadata")
plots_dir    <- file.path(project_path, "Marginals", data_dir_name, "Plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 2. LOAD REGIONAL MASTER DATA
# =====================================================================
master_rds_path <- file.path(metadata_dir, paste0(data_dir_name, "_Master_List_Filtered.rds"))
cat("Loading Master Regional Data...\n")
master_regional_data <- readRDS(master_rds_path)

# =========================================================================
# 3. SETUP TRAIN & TEST STATIONS
# =========================================================================
set.seed(42) # For reproducibility
all_stations <- unique(master_regional_data$Final_Summary$Station_ID)

# Pick 20 stations for training
train_stations <- sample(all_stations, 20)

# Pick 20 DIFFERENT stations for testing
remaining_stations <- setdiff(all_stations, train_stations)
test_stations <- sample(remaining_stations, min(20, length(remaining_stations)))

cat(sprintf("Selected %d Training Stations.\n", length(train_stations)))
cat(sprintf("Selected %d Testing Stations.\n\n", length(test_stations)))

# =========================================================================
# 4. TRAINING PHASE: FIND THE REGIONAL RATIO (A / alpha) & PLOT
# =========================================================================
cat("Training Phase: Calculating Regional Ratio and Plotting...\n")

surge_4param <- function(k, A, alpha, beta, gamma) {
 return(A * (k^alpha) * exp(-(beta * k)^gamma))
}

ratios <- c()

for (stn in train_stations) {
 stn_data <- master_regional_data$Final_Summary %>% filter(Station_ID == stn)
 
 obj_train <- function(p) {
  est <- surge_4param(stn_data$k_hours, p[1], p[2], p[3], p[4])
  return(mean(abs((stn_data$t2_pos - est) / stn_data$t2_pos)))
 }
 
 # Fit full 4-parameter model
 fit <- suppressWarnings(DEoptim(obj_train, lower = c(0.1, 0.01, 0.001, 0.1), 
                                 upper = c(10, 5, 1, 3), 
                                 control = list(itermax = 400, trace = FALSE)))
 
 opt_A <- fit$optim$bestmem[1]
 opt_alpha <- fit$optim$bestmem[2]
 opt_beta <- fit$optim$bestmem[3]
 opt_gamma <- fit$optim$bestmem[4]
 
 ratios <- c(ratios, opt_A / opt_alpha)
 
 # PLOTTING THE TRAINING FIT
 k_seq <- exp(seq(log(min(stn_data$k_hours)), log(max(stn_data$k_hours)), length.out = 500))
 t2_pred <- surge_4param(k_seq, opt_A, opt_alpha, opt_beta, opt_gamma)
 
 p_train <- ggplot() +
  geom_point(data = stn_data, aes(x = k_hours, y = t2_pos), size = 3, shape = 21, fill = "gray70") +
  geom_line(data = data.frame(k = k_seq, t2 = t2_pred), aes(x = k, y = t2), color = "steelblue", linewidth = 1.2) +
  scale_x_log10() +
  labs(title = sprintf("Training Fit: Station %s", stn),
       subtitle = sprintf("4-Parameter Model | Ratio A/alpha = %.2f", opt_A / opt_alpha),
       x = "Temporal scale, k [h] (Log Scale)", 
       y = expression(t[2])) +
  theme_minimal()
 
 # Save plot
 ggsave(file.path(plots_dir, sprintf("Train_Station_%s.png", stn)), plot = p_train, width = 7, height = 5, bg = "white")
}

R_mean <- mean(ratios)
cat(sprintf("=== Regional Parameter Link ===\nMean Ratio (A / alpha) = %.4f\n\n", R_mean))

# =========================================================================
# 5. TESTING PHASE: EVALUATE ON 20 UNSEEN STATIONS & PLOT
# =========================================================================
cat("Testing Phase: Fitting Tails, Extrapolating Heads, and Plotting...\n")

# 3-Parameter Linked Model
surge_linked <- function(k, alpha, beta, gamma) {
 linked_A <- R_mean * alpha
 return(linked_A * (k^alpha) * exp(-(beta * k)^gamma))
}

# Pre-allocate list for results
test_results_list <- list()

for (i in seq_along(test_stations)) {
 stn <- test_stations[i]
 
 if (i %% 5 == 0) cat(sprintf("Processing test station %d of %d...\n", i, length(test_stations)))
 
 test_data <- master_regional_data$Final_Summary %>% filter(Station_ID == stn)
 tail_data <- test_data %>% filter(k_hours >= 24)
 
 obj_test_tail <- function(p) {
  est <- surge_linked(tail_data$k_hours, p[1], p[2], p[3])
  return(mean(abs((tail_data$t2_pos - est) / tail_data$t2_pos)))
 }
 
 # Fit ONLY to the tail of the current test station
 fit_test <- suppressWarnings(DEoptim(obj_test_tail, lower = c(0.01, 0.001, 0.1), upper = c(5, 1, 3), 
                                      control = list(itermax = 400, trace = FALSE)))
 
 opt_p <- fit_test$optim$bestmem
 
 # Predict values for ALL scales
 predicted_t2 <- surge_linked(test_data$k_hours, opt_p[1], opt_p[2], opt_p[3])
 
 # Calculate Absolute Percentage Error (APE) for each scale
 ape <- abs(test_data$t2_pos - predicted_t2) / test_data$t2_pos * 100
 
 # Store in dataframe
 test_results_list[[i]] <- data.frame(
  Station_ID = stn,
  k_hours = test_data$k_hours,
  APE = ape,
  Scale_Type = ifelse(test_data$k_hours >= 24, "Tail (>= 24h) [Fitted]", "Head (< 24h) [Extrapolated]")
 )
 
 # PLOTTING THE VALIDATION FIT
 k_smooth <- exp(seq(log(min(test_data$k_hours)), log(max(test_data$k_hours)), length.out = 500))
 t2_pred_smooth <- surge_linked(k_smooth, opt_p[1], opt_p[2], opt_p[3])
 
 df_pred <- data.frame(k = k_smooth, t2 = t2_pred_smooth)
 df_orig <- test_data %>%
  mutate(Scale = ifelse(k_hours >= 24, "Tail (>= 24h) [Fitted]", "Head (< 24h) [Hidden/Extrapolated]"))
 
 p_test <- ggplot() +
  geom_line(data = df_pred, aes(x = k, y = t2), color = "#d35400", linewidth = 1.2) +
  geom_point(data = df_orig, aes(x = k_hours, y = t2_pos, fill = Scale), size = 3, shape = 21, color="black") +
  scale_x_log10() +
  scale_fill_manual(values = c("Head (< 24h) [Hidden/Extrapolated]" = "orange", "Tail (>= 24h) [Fitted]" = "midnightblue")) +
  labs(title = sprintf("Validation Fit: Station %s", stn),
       subtitle = sprintf("Trained entirely on tail data | Linked Ratio A/alpha = %.2f", R_mean),
       x = "Temporal scale, k [h] (Log Scale)", 
       y = expression(t[2])) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
 
 # Save plot
 ggsave(file.path(plots_dir, sprintf("Test_Station_%s.png", stn)), plot = p_test, width = 7, height = 5, bg = "white")
}

# Combine all results
all_test_results <- bind_rows(test_results_list)

cat("\nCross-Validation Complete! Plots exported to directory.\n")

# =========================================================================
# 6. VISUALIZE MAPE DISTRIBUTION ACROSS SCALES
# =========================================================================
# Convert k_hours to a factor so ggplot treats them as distinct categorical boxes
p_boxplot <- ggplot(all_test_results, aes(x = factor(k_hours), y = APE, fill = Scale_Type)) +
 geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.alpha = 0.5, alpha = 0.8) +
 scale_fill_manual(values = c("Head (< 24h) [Extrapolated]" = "#f39c12", 
                              "Tail (>= 24h) [Fitted]" = "#2c3e50")) +
 labs(title = "Cross-Validation Error by Temporal Scale",
      subtitle = sprintf("Evaluated on %d unseen stations | Linked Ratio (A/alpha) = %.2f", length(test_stations), R_mean),
      x = "Temporal scale, k [h]", 
      y = "Absolute Percentage Error (%)") +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
       legend.position = "bottom",
       legend.title = element_blank())

# Save the final summary boxplot
ggsave(file.path(plots_dir, "Summary_CrossValidation_Error.png"), plot = p_boxplot, width = 9, height = 6, bg = "white")