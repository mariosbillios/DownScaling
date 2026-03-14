library(ggplot2)
library(dplyr)
library(tidyr)

# ==========================================
# HELPER: FORMAT STATISTICS AS FORMULAS
# ==========================================
get_stat_expr <- function(stat_name) {
 # Using quote() instead of expression() allows seamless injection into bquote() later
 switch(stat_name,
        "p_zero"   = quote(p[0]),
        "mean_all" = quote(mu[all]),
        "var_all"  = quote(sigma[all]^2),
        "mean_pos" = quote(mu["+"]),
        "var_pos"  = quote(sigma["+"]^2),
        "l1_pos"   = quote(lambda[1]^("+")),
        "l2_all"   = quote(lambda[2]^{(all)}),
        "l2_pos"   = quote(lambda[2]^("+")),
        "l3_all"   = quote(lambda[3]^{(all)}),
        "l3_pos"   = quote(lambda[3]^("+")),
        "l4_all"   = quote(lambda[4]^{(all)}),
        "l4_pos"   = quote(lambda[4]^("+")),
        "t2_all"   = quote(tau[2]^{(all)}),
        "t2_pos"   = quote(tau[2]^("+")),
        "t3_all"   = quote(tau[3]^{(all)}),
        "t3_pos"   = quote(tau[3]^("+")),
        "t4_all"   = quote(tau[4]^{(all)}),
        "t4_pos"   = quote(tau[4]^("+")),
        parse(text = stat_name)[[1]] # Fallback if not listed
 )
}

# ==========================================
# CUSTOM MOMENTS PLOTTING FUNCTION
# ==========================================
plot_moments <- function(data, metric, x_var = "Scale_k") {
 
 col_pred         <- paste0("Pred_", metric)
 col_actual       <- paste0("Actual_", metric)
 col_derived_pred   <- paste0("derived_Pred_", metric)
 col_derived_actual <- paste0("derived_Actual_", metric)
 
 req_cols <- c(x_var, col_pred, col_actual, col_derived_pred, col_derived_actual)
 missing_cols <- setdiff(req_cols, names(data))
 if (length(missing_cols) > 0) {
  stop("Missing columns in your dataset: ", paste(missing_cols, collapse = ", "))
 }
 
 metric_expr <- get_stat_expr(metric)
 
 lbl_pred         <- "Predicted"
 lbl_actual       <- "Actual"
 lbl_derived_pred   <- "Derived Pred"
 lbl_derived_actual <- "Derived Actual"
 
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
  scale_x_log10(breaks = c(1, 2, 4, 7, 12, 24, 48, 96, 240)) + 
  scale_color_manual(name = "Legend", values = legend_colors) +
  labs(
   x = "Aggregation Scale k (hours)",
   y = metric_expr,
   title = bquote("Comparison of Actual vs Predicted and Derived:" ~ .(metric_expr))
  ) +
  theme_bw() +
  theme(
   panel.grid.major = element_line(linetype = "dotted", color = "gray70"),
   panel.grid.minor = element_blank(),
   legend.position = "right",
   legend.background = element_blank(),
   legend.key = element_blank()
  )
 
 return(p)
}

# ==========================================
# MASTER GENERATOR FUNCTION
# ==========================================
generate_station_plots <- function(station_name, individual_dir) {
 
 station_folder <- file.path(individual_dir, station_name)
 plots_folder <- file.path(station_folder, "plots")
 dir.create(plots_folder, showWarnings = FALSE) 
 
 tryCatch({
  preds_df <- readRDS(file.path(station_folder, "all_predictions.rds"))
  models_df <- readRDS(file.path(station_folder, "final_evaluated_models.rds"))
  derived_df <- readRDS(file.path(station_folder, "derived_l_moments.rds"))
  meta_obj <- readRDS(file.path(station_folder, "extra_metadata.rds"))
  best_models <- readRDS(file.path(station_folder, "best_models.rds"))
  
  # ---------------------------------------------------------
  # RESTRICT DOMAIN & REMOVE l1_all
  # ---------------------------------------------------------
  preds_df <- preds_df %>% filter(Scale_k >= 1 & Scale_k <= 240, Statistic != "l1_all")
  derived_df <- derived_df %>% filter(Scale_k >= 1 & Scale_k <= 240)
  models_df <- models_df %>% filter(Statistic != "l1_all")
  best_models <- best_models %>% filter(Statistic != "l1_all")
  
  all_stats <- unique(preds_df$Statistic)
  
  # ==========================================
  # 1. SCALING PLOTS
  # ==========================================
  for (stat in all_stats) {
   df_stat <- preds_df %>% filter(Statistic == stat)
   metric_expr <- get_stat_expr(stat)
   
   # Extract Best Model
   best_mod <- best_models %>% filter(Statistic == stat) %>% pull(Model)
   best_mod_label <- gsub("Weibull_", "Wei-", gsub("PowerLaw_", "Lom-", best_mod))
   
   # Define active models based on statistic
   if (stat == "p_zero") {
    active_models <- c("Weibull_1p", "Weibull_2p", "PowerLaw_1p", "PowerLaw_2p")
   } else {
    active_models <- c("Weibull_2p", "PowerLaw_2p")
   }
   
   # Pivot models into long format for line aesthetics
   df_lines <- df_stat %>%
    select(Scale_k, all_of(active_models)) %>%
    pivot_longer(cols = -Scale_k, names_to = "Model", values_to = "Value") %>%
    mutate(Model = factor(Model, 
                          levels = c("PowerLaw_1p", "PowerLaw_2p", "Weibull_1p", "Weibull_2p"),
                          labels = c("Lom-1p", "Lom-2p", "Wei-1p", "Wei-2p"))) %>%
    filter(!is.na(Value))
   
   # Phase points
   df_points <- df_stat %>%
    mutate(Phase = factor(ifelse(Scale_k >= 24, "Calib. data", "Valid. data"),
                          levels = c("Calib. data", "Valid. data")))
   
   p_scale <- ggplot() +
    geom_vline(xintercept = 24, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_line(data = df_lines, aes(x = Scale_k, y = Value, color = Model, linetype = Model), linewidth = 1) +
    geom_point(data = df_points, aes(x = Scale_k, y = Actual, fill = Phase), 
               shape = 21, color = "transparent", size = 2.5) +
    scale_x_log10(breaks = c(1, 2, 4, 7, 12, 24, 48, 96, 240)) +
    scale_color_manual(name = "", values = c("Lom-1p" = "red", "Lom-2p" = "red", "Wei-1p" = "blue", "Wei-2p" = "blue")) +
    scale_linetype_manual(name = "", values = c("Lom-1p" = "solid", "Lom-2p" = "dashed", "Wei-1p" = "solid", "Wei-2p" = "dashed")) +
    scale_fill_manual(name = "", values = c("Calib. data" = "black", "Valid. data" = "orange")) +
    labs(
     title = bquote("Scaling Behavior of" ~ .(metric_expr)), 
     subtitle = paste("Best Model:", best_mod_label), 
     x = "Aggregation Scale k (hours)", 
     y = metric_expr
    ) +
    theme_bw() + 
    theme(
     panel.grid.major = element_line(linetype = "dotted", color = "gray70"),
     panel.grid.minor = element_blank(),
     legend.position = "right",
     legend.background = element_blank(), 
     legend.key = element_blank(),        
     legend.spacing.y = unit(0, "pt"),
     legend.margin = margin(t = 2, r = 5, b = 2, l = 5)
    ) +
    guides(
     fill = guide_legend(order = 1, override.aes = list(size = 3)),
     color = guide_legend(order = 2),
     linetype = guide_legend(order = 2)
    )
   
   ggsave(file.path(plots_folder, paste0("1_Scaling_", stat, ".png")), plot = p_scale, width = 7, height = 6)
  }
  
  # ==========================================
  # 2. QQ PLOTS: FACETED MODEL VS EMPIRICAL
  # ==========================================
  for (stat in all_stats) {
   metric_expr <- get_stat_expr(stat)
   
   # Define active models based on statistic
   if (stat == "p_zero") {
    active_models <- c("Weibull_1p", "Weibull_2p", "PowerLaw_1p", "PowerLaw_2p")
    plot_width <- 12 # Wider to fit 4 facets
   } else {
    active_models <- c("Weibull_2p", "PowerLaw_2p")
    plot_width <- 7  # Standard width for 2 facets
   }
   
   # Pivot to evaluate dynamic models simultaneously via facets
   df_qq <- preds_df %>% 
    filter(Statistic == stat) %>%
    select(Scale_k, Actual, all_of(active_models)) %>%
    pivot_longer(cols = all_of(active_models),
                 names_to = "Model", values_to = "Predicted") %>%
    mutate(
     Model = factor(Model, 
                    levels = c("PowerLaw_1p", "PowerLaw_2p", "Weibull_1p", "Weibull_2p"),
                    labels = c("Lom-1p", "Lom-2p", "Wei-1p", "Wei-2p")),
     Phase = factor(ifelse(Scale_k >= 24, "Calib", "Valid"),
                    levels = c("Calib", "Valid"))
    ) %>%
    filter(!is.na(Predicted) & !is.na(Actual))
   
   p_qq <- ggplot(df_qq, aes(x = Actual, y = Predicted, color = Phase)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
    geom_point(size = 2.5, alpha = 0.8) +
    facet_wrap(~ Model, nrow = 1) +
    scale_color_manual(name = "", values = c("Calib" = "black", "Valid" = "orange")) +
    labs(
     title = bquote("Model vs empirical values for statistic:" ~ .(metric_expr)),
     x = "Empirical values",
     y = "Model values"
    ) +
    theme_bw() +
    theme(
     panel.grid.major = element_line(linetype = "dotted", color = "gray70"),
     panel.grid.minor = element_blank(),
     strip.background = element_rect(fill = "gray90"),
     legend.position = "right",
     legend.background = element_blank(), 
     legend.key = element_blank()         
    )
   
   ggsave(file.path(plots_folder, paste0("2_QQ_", stat, ".png")), plot = p_qq, width = plot_width, height = 4)
  }
  
  # ==========================================
  # 3. MODEL PERFORMANCE DUEL (LOLLIPOP / SMAPE)
  # ==========================================
  duel_data <- models_df %>%
   filter(!is.na(smape_validation), Model %in% c("Weibull_2p", "PowerLaw_2p")) %>%
   select(Statistic, Model, smape_validation) %>%
   pivot_wider(names_from = Model, values_from = smape_validation) %>%
   mutate(Statistic = factor(Statistic))
  
  p_duel <- ggplot(duel_data) +
   geom_segment(aes(x = Weibull_2p, xend = PowerLaw_2p, y = Statistic, yend = Statistic), 
                color = "gray60", linetype = "dashed", linewidth = 0.8) +
   geom_point(aes(x = Weibull_2p, y = Statistic, color = "Wei-2p"), size = 3.5) +
   geom_point(aes(x = PowerLaw_2p, y = Statistic, color = "Lom-2p"), size = 3.5) +
   scale_color_manual(values = c("Wei-2p" = "blue", "Lom-2p" = "red")) +
   scale_y_discrete(labels = function(breaks) {
    do.call(expression, lapply(breaks, get_stat_expr))
   }) +
   labs(title = "Validation SMAPE Comparison by Statistic", subtitle = "Lower is better",
        x = "SMAPE Validation Error (%)", y = "Statistic", color = "Model") +
   theme_bw() + 
   theme(
    panel.grid.major = element_line(linetype = "dotted", color = "gray70"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.background = element_blank(), 
    legend.key = element_blank(),        
    axis.text.y = element_text(size = 12) 
   )
  
  ggsave(file.path(plots_folder, "3_SMAPE_Duel.png"), plot = p_duel, width = 8, height = 6)
  
  # ==========================================
  # 4. ACTUAL VS PREDICTED & DERIVED PLOTS
  # ==========================================
  moment_metrics <- c("t2_all", "t3_all", "t4_all", "t2_pos", "t3_pos", "t4_pos")
  
  for (m_stat in moment_metrics) {
   req_cols <- c("Scale_k", paste0("Pred_", m_stat), paste0("Actual_", m_stat), 
                 paste0("derived_Pred_", m_stat), paste0("derived_Actual_", m_stat))
   
   if (all(req_cols %in% names(derived_df))) {
    p_moments <- plot_moments(data = derived_df, metric = m_stat, x_var = "Scale_k")
    ggsave(file.path(plots_folder, paste0("4_Moments_", m_stat, ".png")), plot = p_moments, width = 8, height = 6)
   }
  }
  
  # ==========================================
  # 5. DATA LOSS WATERFALL
  # ==========================================
  data_loss <- meta_obj$Scale_Data_Loss %>% filter(Scale_k_hours >= 1 & Scale_k_hours <= 240)
  p_loss <- ggplot(data_loss, aes(x = as.factor(Scale_k_hours), y = Dropped_Pct)) +
   geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
   geom_text(aes(label = round(Dropped_Pct, 1)), vjust = -0.5, size = 3) +
   labs(title = "Data Loss per Aggregation Scale", subtitle = "Percentage of theoretical bins dropped due to missing data",
        x = "Aggregation Scale k (Hours)", y = "Dropped Bins (%)") +
   theme_minimal()
  
  ggsave(file.path(plots_folder, "5_Data_Loss.png"), plot = p_loss, width = 8, height = 5)
  
  return(paste("SUCCESS: Plots generated and saved for", station_name))
  
 }, error = function(e) {
  return(paste("FAILED to generate plots for", station_name, "-", conditionMessage(e)))
 })
}

# ==========================================
# HOW TO USE IT:
# ==========================================

gdrive_path <- paste(LETTERS[file.exists(paste0(LETTERS, ":/My Drive"))], ":/My Drive", sep = "")
project_path <- file.path(gdrive_path, "Academic_git/DownScaling")

data_dir_name <- "QC_d data - Germany" 
data_dir <- file.path(gdrive_path, "General_Data/GSDR", data_dir_name)
individual_dir <- file.path(project_path, "Marginals", data_dir_name, "Testing") 
station_code <- "DE_00003"

generate_station_plots(station_name = "DE_00003", individual_dir = individual_dir)