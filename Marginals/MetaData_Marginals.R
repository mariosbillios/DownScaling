# Load necessary libraries
library(dplyr)
library(lubridate)
library(tools)
library(writexl)


# =====================================================================
# 1. DEFINE PATHS & DIRECTORIES
# =====================================================================
gdrive_path <- "G:/My Drive"
project_path <- file.path(gdrive_path, "Academic_git/DownScaling")

# Define Data directory
data_dir_name <- "QC_d data - Germany" 
data_dir <- file.path(gdrive_path, "General_Data/GSDR", data_dir_name)

# UPDATED: Station Results Directory
individual_dir <- file.path(project_path, "Marginals", data_dir_name, "Station_Results")

# Define and create the Metadata export directory
metadata_export_dir <- file.path(project_path, "Marginals", data_dir_name, "Metadata" )
dir.create(metadata_export_dir, recursive = TRUE, showWarnings = FALSE)



# =====================================================================
# 2. HELPER FUNCTION TO EXTRACT FIELDS FROM TXT HEADER
# =====================================================================

extract_field <- function(lines, key) {
 matched_line <- grep(paste0("^", key, ":"), lines, value = TRUE)
 if (length(matched_line) > 0) {
  val <- sub(paste0("^", key, ":\\s*"), "", matched_line[1])
  return(trimws(val))
 } else {
  return(NA)
 }
}

# =====================================================================
# 3. RUN EXTRACTION LOOP (TXT + RDS)
# =====================================================================
station_files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)
metadata_list <- list()

cat("Starting master metadata extraction for", length(station_files), "files...\n")




for (current_file in station_files) {
 
 
 station_name <- tools::file_path_sans_ext(basename(current_file))
 
 # --- A. READ TXT HEADER ---
 header_lines <- readLines(current_file, n = 21, warn = FALSE)
 
 # --- B. INITIALIZE DEFAULT RESULTS ---
 ds_successful      <- FALSE
 explicit_nas       <- NA
 missing_rows       <- NA
 total_missing      <- NA
 calc_missing_pct   <- NA
 
 # --- C. READ DOWNSCALING RESULTS (If they exist) ---
 station_folder <- file.path(individual_dir, station_name)
 extra_meta_file <- file.path(station_folder, "extra_metadata.rds")
 
 if (dir.exists(station_folder) && file.exists(extra_meta_file)) {
  tryCatch({
   meta_obj <- readRDS(extra_meta_file)
   ds_successful <- TRUE
   
   mds <- meta_obj$Missing_Data_Summary
   if (!is.null(mds)) {
    explicit_nas       <- mds$Explicit_NAs
    missing_rows       <- mds$Missing_Rows
    total_missing      <- mds$Total_Missing_Points
    calc_missing_pct   <- mds$Calculated_Missing_Pct
   }
  }, error = function(e) {
   warning(paste("Failed to read RDS for", station_name))
  })
 }
 
 # --- D. BUILD RAW ROW ---
 metadata_list[[station_name]] <- data.frame(
  Station_ID              = extract_field(header_lines, "Station ID"),
  Country                 = extract_field(header_lines, "Country"),
  Original_Station_Number = extract_field(header_lines, "Original Station Number"),
  Original_Station_Name   = extract_field(header_lines, "Original Station Name"),
  Original_Data_Path      = extract_field(header_lines, "Path to original data"),
  
  Latitude                = as.numeric(extract_field(header_lines, "Latitude")),
  Longitude               = as.numeric(extract_field(header_lines, "Longitude")),
  Elevation_m             = as.numeric(gsub("m", "", extract_field(header_lines, "Elevation"))), 
  
  Start_Datetime          = ymd_h(extract_field(header_lines, "Start datetime")),
  End_Datetime            = ymd_h(extract_field(header_lines, "End datetime")),
  Time_Zone               = extract_field(header_lines, "Time Zone"),
  Daylight_Saving_Info    = extract_field(header_lines, "Daylight Saving info"),
  
  Original_Timestep       = extract_field(header_lines, "Original Timestep"),
  New_Timestep            = extract_field(header_lines, "New Timestep"),
  Original_Units          = extract_field(header_lines, "Original Units"),
  New_Units               = extract_field(header_lines, "New Units"),
  
  Header_Records          = as.numeric(extract_field(header_lines, "Number of records")),
  Header_Missing_Pct      = as.numeric(extract_field(header_lines, "Percent missing data")),
  No_Data_Value           = as.numeric(extract_field(header_lines, "No data value")),
  Resolution              = as.numeric(extract_field(header_lines, "Resolution")),
  Other_Info              = extract_field(header_lines, "Other"),
  
  Downscaling_Successful  = ds_successful,
  Explicit_NAs            = explicit_nas,
  Missing_Rows            = missing_rows,
  Total_Missing_Points    = total_missing,
  Calculated_Missing_Pct  = calc_missing_pct,
  
  stringsAsFactors = FALSE
 )
}

# 4. Combine all stations
global_metadata_df <- bind_rows(metadata_list)

# =====================================================================
# 5. POST-PROCESSING: Calculate Years & Check Mismatches
# =====================================================================
global_metadata_df <- global_metadata_df %>%
 mutate(
  # 1. Calculate the Years
  Start_Year     = Start_Datetime,
  End_Year       = End_Datetime,
  Duration_Years = round(time_length(interval(Start_Datetime, End_Datetime), "years"), 2),
  
  # 2. Round calculated missing to 1 decimal
  Calculated_Missing_Pct = round(Calculated_Missing_Pct, 2),
  
  # 3. Create the mismatch flag (Comparing Header to Calculated)
  # Using case_when to safely handle NA values if downscaling failed
  Missing_Pct_Mismatch = case_when(
   is.na(Calculated_Missing_Pct) ~ NA, 
   Header_Missing_Pct != Calculated_Missing_Pct ~ TRUE,
   TRUE ~ FALSE
  )
 )



# =====================================================================
# 6. EXPORT TO RDS AND XLSX
# =====================================================================
export_base_name <- paste0(data_dir_name, "_Master_Metadata")

rds_path  <- file.path(metadata_export_dir, paste0(export_base_name, ".rds"))
xlsx_path <- file.path(metadata_export_dir, paste0(export_base_name, ".xlsx"))

# Save Master Files
saveRDS(global_metadata_df, rds_path)
write_xlsx(global_metadata_df, xlsx_path)


#--------------------------------------------------------------------------------

# Load necessary libraries
library(dplyr)
library(readxl)

# =====================================================================
# 1. DEFINE PATHS
# =====================================================================
gdrive_path <- "G:/My Drive"
project_path <- file.path(gdrive_path, "Academic_git/DownScaling")

data_dir_name <- "QC_d data - Germany" 
metadata_dir <- file.path(project_path, "Marginals", data_dir_name, "Metadata")
results_dir <- file.path(project_path, "Marginals", data_dir_name, "Station_Results")

# =====================================================================
# 2. LOAD AND FILTER METADATA
# =====================================================================
# We can read the .rds version of your master metadata since it retains strict data types
meta_path <- file.path(metadata_dir, paste0(data_dir_name, "_Master_Metadata.rds"))
meta_df <- readRDS(meta_path)

# Apply the strict quality filters
quality_stations <- meta_df %>%
 filter(
  Downscaling_Successful == TRUE,
  Duration_Years >= 10,
  Calculated_Missing_Pct <= 2
 )

cat("Filtering Complete:\n")
cat(" - Original Stations:", nrow(meta_df), "\n")
cat(" - High-Quality Stations (Passed criteria):", nrow(quality_stations), "\n\n")

# =====================================================================
# 3. INITIALIZE LISTS FOR COMBINING
# =====================================================================
# We will collect the dataframes from all stations into these lists
list_final_summary   <- list()
list_opt_params      <- list()
list_eval_models     <- list()
list_best_models     <- list()
list_all_predictions <- list()

cat("Merging RDS results for valid stations...\n")

# =====================================================================
# 4. LOOP THROUGH QUALITY STATIONS & COMBINE RESULTS
# =====================================================================
for (i in seq_len(nrow(quality_stations))) {
 
 # Assuming your folder names match the Station_ID (e.g., "DE_00003")
 st_id <- quality_stations$Station_ID[i]
 st_folder <- file.path(results_dir, st_id)
 
 # Define paths to the station's RDS files
 file_summary     <- file.path(st_folder, "final_summary.rds")
 file_params      <- file.path(st_folder, "optimized_parameters.rds")
 file_eval_models <- file.path(st_folder, "final_evaluated_models.rds")
 file_best_models <- file.path(st_folder, "best_models.rds")
 file_predictions <- file.path(st_folder, "all_predictions.rds")
 
 # Safely read and append each file if it exists
 # We use mutate(Station_ID = st_id) so you can track the data after merging!
 
 if (file.exists(file_summary)) {
  list_final_summary[[st_id]] <- readRDS(file_summary) %>% mutate(Station_ID = st_id)
 }
 
 if (file.exists(file_params)) {
  list_opt_params[[st_id]] <- readRDS(file_params) %>% mutate(Station_ID = st_id)
 }
 
 if (file.exists(file_eval_models)) {
  list_eval_models[[st_id]] <- readRDS(file_eval_models) %>% mutate(Station_ID = st_id)
 }
 
 if (file.exists(file_best_models)) {
  list_best_models[[st_id]] <- readRDS(file_best_models) %>% mutate(Station_ID = st_id)
 }
 
 if (file.exists(file_predictions)) {
  list_all_predictions[[st_id]] <- readRDS(file_predictions) %>% mutate(Station_ID = st_id)
 }
}

# =====================================================================
# 5. BIND ROWS AND CREATE MASTER OBJECT
# =====================================================================
# bind_rows combines the hundreds of tiny dataframes into massive regional dataframes
master_regional_data <- list(
 Quality_Metadata         = quality_stations, # Keeps the metadata of only the stations used
 Final_Summary            = bind_rows(list_final_summary),
 Optimized_Parameters     = bind_rows(list_opt_params),
 Final_Evaluated_Models   = bind_rows(list_eval_models),
 Best_Models              = bind_rows(list_best_models),
 All_Predictions          = bind_rows(list_all_predictions)
)

# =====================================================================
# 6. EXPORT MASTER RDS
# =====================================================================
# Name it exactly as the data_dir_name (e.g., "QC_d data - Germany.rds")
export_master_path <- file.path(metadata_dir, paste0(data_dir_name, "_Master_List_Filtered" ,".rds"))

saveRDS(master_regional_data, export_master_path)

cat("\n--- MERGE COMPLETE ---\n")
cat("Successfully created Regional Master Data Object.\n")
cat("Saved to:", export_master_path, "\n\n")

# Print a quick preview of the resulting sizes
cat("Master Object Contents:\n")
cat(" - Final_Summary Rows:          ", nrow(master_regional_data$Final_Summary), "\n")
cat(" - Optimized_Parameters Rows:   ", nrow(master_regional_data$Optimized_Parameters), "\n")
cat(" - Best_Models Rows:            ", nrow(master_regional_data$Best_Models), "\n")