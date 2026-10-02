clean_headers <- function(file_path, remove_second_lactose_peak = FALSE) {
  # Read the first two rows to get headers
  header1 <- read_excel(file_path, sheet = "Rawdata", range = "A1:ZZ1", col_names = FALSE)
  header2 <- read_excel(file_path, sheet = "Rawdata", range = "A2:ZZ2", col_names = FALSE)
  batch_number = gsub("Batch", "", basename(file_path))
  batch_number = gsub(".xlsx", "", batch_number, fixed = TRUE)
  
  # Read the actual data 
  data <- read_excel(file_path, sheet = "Rawdata")
  
  # Convert to character vectors
  h1 <- as.character(header1[1, ])
  h2 <- as.character(header2[1, ])
  
  # Step 1: Fill in "Sample" for columns 2-7 in row 1
  for(i in 2:7) {
    if(h1[i] == "NA" || h1[i] == "") {
      h1[i] <- "Sample"
    }
  }
  
  # Step 2: Fill in compound names for their paired RT/Resp columns
  # Find all compound names that end with "Results" and fill their next column
  i <- 8  # Start from column 8 as specified
  while(i <= length(h1)) {
    if(!is.na(h1[i]) && grepl("Results$", h1[i])) {
      # Fill the next column with the same compound name
      if(i + 1 <= length(h1)) {
        h1[i + 1] <- h1[i]
      }
      i <- i + 2  # Skip to next potential compound (they come in pairs)
    } else {
      i <- i + 1
    }
  }
  
  # cat("After filling row 1:\n")
  # cat("First 20 values of h1:", h1[1:min(20, length(h1))], "\n")
  
  # Ensure headers match data length
  if(length(h1) > ncol(data)) {
    h1 <- h1[1:ncol(data)]
    h2 <- h2[1:ncol(data)]
  }
  
  # Step 4: Create final headers
  final_headers <- character(length(h2))
  
  for(i in 1:length(h2)) {
    # If row 1 has "Sample", leave row 2 as is
    if(h1[i] != "NA" && h1[i] == "Sample") {
      final_headers[i] <- h2[i]
    }
    # For Results columns, combine compound name with RT/Resp
    else if(h1[i] != "NA" && grepl("Results$", h1[i]) && 
            h2[i] != "NA" && h2[i] %in% c("RT", "Resp.")) {
      # Remove "Results" from compound name and combine
      compound_name <- gsub(" Results$", "", h1[i])
      final_headers[i] <- paste(compound_name, h2[i], sep = " ")
    }
    # For other cases, use row 2 if available
    else if(h2[i] != "NA" && h2[i] != "") {
      final_headers[i] <- h2[i]
    }
  }
  
  # Remove empty headers and corresponding columns
  keep_cols <- which(final_headers != "" & final_headers != "NA")
  data <- data[, keep_cols]
  final_headers <- final_headers[keep_cols]
  
  # Apply the processed headers
  colnames(data) <- final_headers
  
  # Remove h2 and empty columns
  data <- data[-1, ]
  # data <- data[, -(1:2)]
  if(remove_second_lactose_peak) {
    lactose_peak_1_to_remove = grep("Lactose Pk1 361", colnames(data))
    lactose_peak_2_to_remove = grep("Lactose Pk2 361", colnames(data))
    to_rm = c(lactose_peak_1_to_remove, lactose_peak_2_to_remove)
    if(length(to_rm) > 0) {
      data = data[, -to_rm]
    }
  }
  
  data$batch = batch_number
  
  return(data)
}

clean_headers_new <- function(file_path, remove_second_lactose_peak = FALSE,
                              measures = c("RT", "Resp.")) {
  batch_number <- gsub("Batch|\\.xlsx$", "", basename(file_path))
  
  # Read the whole sheet once, so the headers and data always line up
  raw  <- read_excel(file_path, sheet = "Rawdata", col_names = FALSE, col_types = "text")
  h1   <- unlist(raw[1, ], use.names = FALSE)   # compound names ("... Results")
  h2   <- unlist(raw[2, ], use.names = FALSE)   # Data File, Type, ..., Resp., RT
  data <- raw[-(1:2), ]
  
  is_blank <- function(x) is.na(x) | x %in% c("", "NA")
  
  # Step 1: Compound name for each column
  cmp <- ifelse(!is_blank(h1) & grepl("Results$", h1), sub(" Results$", "", h1), NA)
  
  # Step 2: Carry the compound name into the empty cells that follow it,
  # but only where row 2 holds a measure (e.g. the RT column after Resp.)
  for (i in seq_along(cmp)[-1]) {
    if (is.na(cmp[i]) && is_blank(h1[i]) && !is.na(cmp[i - 1]) && h2[i] %in% measures) {
      cmp[i] <- cmp[i - 1]
    }
  }
  
  # Step 3: Final headers: "Compound Measure" for compound columns, row 2 for everything else
  final_headers <- ifelse(!is.na(cmp) & h2 %in% measures, paste(cmp, h2),
                          ifelse(!is_blank(h2), h2, ""))
  
  # Remove columns without a header
  keep <- final_headers != ""
  data <- data[, keep]
  final_headers <- final_headers[keep]
  
  if (anyDuplicated(final_headers)) {
    warning("Duplicated column names in ", basename(file_path), ": ",
            paste(unique(final_headers[duplicated(final_headers)]), collapse = ", "))
  }
  colnames(data) <- final_headers
  
  if (remove_second_lactose_peak) {
    to_rm <- grep("Lactose Pk1 361|Lactose Pk2 361", colnames(data))
    if (length(to_rm) > 0) data <- data[, -to_rm]
  }
  
  data$batch <- batch_number
  data
}


introduce_known_concentrations = function(df, custom_standard_map = NULL) {

    galactose_known_con = c(
        0,
        0.6372,
        0.8496,
        1.062,
        1.6992,
        2.3364
    )

    lactose_known_con = c(
        0,
        0.856,
        1.284,
        1.926,
        2.889,
        3.852
    )

    galactitol_known_con = c(
        0,
        0.11585,
        0.173775,
        0.2317,
        0.34755,
        0.4634
    )

    galactonate_known_con = c(
        0,
        1.0385,
        1.55775,
        2.077,
        3.1155,
        4.154
    )
    
    if (!is.null(custom_standard_map)) {
      
      known_concentration_idx = (nrow(df)-11):nrow(df)
      df[known_concentration_idx, "galactose_concentration_ug_per_100ul"] = rep(galactose_known_con, 2)
      df[known_concentration_idx, "lactose_concentration_ug_per_100ul"] = rep(lactose_known_con, 2)
      df[known_concentration_idx, "galactitol_concentration_ug_per_100ul"] = rep(galactitol_known_con, 2)
      df[known_concentration_idx, "galactonate_concentration_ug_per_100ul"] = rep(galactonate_known_con, 2)
      
    } else {
      known_concentration_idx = (nrow(df)-5):nrow(df)
      df[known_concentration_idx, "galactose_concentration_ug_per_100ul"] = galactose_known_con
      df[known_concentration_idx, "lactose_concentration_ug_per_100ul"] = lactose_known_con
      df[known_concentration_idx, "galactitol_concentration_ug_per_100ul"] = galactitol_known_con
      df[known_concentration_idx, "galactonate_concentration_ug_per_100ul"] = galactonate_known_con
    }

    return(df)

}

introduce_known_concentrations_special_batch = function(df, custom_standard_map = NULL) {
  
  galactose_known_con = c(
    0,
    0.6372,
    0.8496,
    1.062,
    1.6992,
    2.3364,
    4.6728,
    9.3456,
    21.24
  )
  
  lactose_known_con = c(
    0,
    0.856,
    1.284,
    1.926,
    2.889,
    3.852,
    7.704,
    19.26,
    42.8
  )
  
  galactitol_known_con = c(
    0,
    0.11585,
    0.173775,
    0.2317,
    0.34755,
    0.4634,
    0.9268,
    2.317,
    4.634
  )
  
  galactonate_known_con = c(
    0,
    1.0385,
    1.55775,
    2.077,
    3.1155,
    4.154,
    8.308,
    20.77,
    41.54
  )
  
  if (!is.null(custom_standard_map)) {
    
    known_concentration_idx = (nrow(df)-17):nrow(df)
    df[known_concentration_idx, "galactose_concentration_ug_per_100ul"] = rep(galactose_known_con, 2)
    df[known_concentration_idx, "lactose_concentration_ug_per_100ul"] = rep(lactose_known_con, 2)
    df[known_concentration_idx, "galactitol_concentration_ug_per_100ul"] = rep(galactitol_known_con, 2)
    df[known_concentration_idx, "galactonate_concentration_ug_per_100ul"] = rep(galactonate_known_con, 2)
    
  } else {
    known_concentration_idx = (nrow(df)-8):nrow(df)
    df[known_concentration_idx, "galactose_concentration_ug_per_100ul"] = galactose_known_con
    df[known_concentration_idx, "lactose_concentration_ug_per_100ul"] = lactose_known_con
    df[known_concentration_idx, "galactitol_concentration_ug_per_100ul"] = galactitol_known_con
    df[known_concentration_idx, "galactonate_concentration_ug_per_100ul"] = galactonate_known_con
  }
  
  return(df)
  
}

subract_matrix_blanks = function(df, custom_standard_map = NULL) {

    # Lookup column index of the first target compound (from left to right) with "final" values:
    # which(colnames(df) == "galactose_ion_area") 
    # fst_col = 25 #Define fst_col as column index of first target compound
    fst_col = which(colnames(df) == "galactose_ion_area") 

    if (!is.null(custom_standard_map)) {
      
      # Subtract matrix blank signal from standard addition signals
      galactose_blank_areas = df %>% filter(grepl("Std_0", Name)) %>% pull(galactose_ion_area) %>% rep(each = 5)
      # print(galactose_blank_area)
      df[(nrow(df) - 11):nrow(df), fst_col] = df[(nrow(df) - 11):nrow(df), fst_col] - galactose_blank_areas
      
      lactose_blank_areas = df %>% filter(grepl("Std_0", Name)) %>% pull(lactose_ion_area) %>% rep(each = 5)
      # print(lactose_blank_area)
      df[(nrow(df) - 11):nrow(df), fst_col + 1] = df[(nrow(df) - 11):nrow(df), fst_col + 1] - lactose_blank_areas
      
      galactitol_blank_areas = df %>% filter(grepl("Std_0", Name)) %>% pull(galactitol_ion_area) %>% rep(each = 5)
      # print(galactitol_blank_area)
      df[(nrow(df) - 11):nrow(df), fst_col + 2] = df[(nrow(df) - 11):nrow(df), fst_col + 2] - galactitol_blank_areas
      
      galactonate_blank_areas = df %>% filter(grepl("Std_0", Name)) %>% pull(galactonate_ion_area) %>% rep(each = 5)
      # print(galactonate_blank_area)
      df[(nrow(df) - 11):nrow(df), fst_col + 3] = df[(nrow(df) - 11):nrow(df), fst_col + 3] - galactonate_blank_areas
      
    } else {
      
      # Subtract matrix blank signal from standard addition signals
      galactose_blank_area = df[nrow(df) - 5, fst_col]
      # print(galactose_blank_area)
      for (i in (nrow(df)-5) : nrow(df)) {
        df[i, "galactose_ion_area"] = df[i, "galactose_ion_area"] - galactose_blank_area
      }
      
      lactose_blank_area = df[nrow(df) - 5, fst_col+1]
      # print(lactose_blank_area)
      for (i in (nrow(df)-5) : nrow(df)) {
        df[i, "lactose_ion_area"] = df[i, "lactose_ion_area"] - lactose_blank_area
      }
      
      galactitol_blank_area = df[nrow(df) - 5, fst_col+2]
      # print(galactitol_blank_area)
      for (i in (nrow(df)-5) : nrow(df)) {
        df[i, "galactitol_ion_area"] = df[i, "galactitol_ion_area"] - galactitol_blank_area
      }
      
      galactonate_blank_area = df[nrow(df) - 5, fst_col+3]
      # print(galactonate_blank_area)
      for (i in (nrow(df)-5) : nrow(df)) {
        df[i, "galactonate_ion_area"] = df[i, "galactonate_ion_area"] - galactonate_blank_area
      }
      
    }
    

    return(df)

}

subract_matrix_blanks_special_batch = function(df, custom_standard_map = NULL) {
  
  # Lookup column index of the first target compound (from left to right) with "final" values:
  # which(colnames(df) == "galactose_ion_area") 
  # fst_col = 25 #Define fst_col as column index of first target compound
  fst_col = which(colnames(df) == "galactose_ion_area") 
  
  if (!is.null(custom_standard_map)) {
    
    # Subtract matrix blank signal from standard addition signals
    galactose_blank_areas = df %>% filter(grepl("Std_0", Name)) %>% pull(galactose_ion_area) %>% rep(each = 9)
    # print(galactose_blank_area)
    df[(nrow(df) - 17):nrow(df), fst_col] = df[(nrow(df) - 17):nrow(df), fst_col] - galactose_blank_areas
    
    lactose_blank_areas = df %>% filter(grepl("Std_0", Name)) %>% pull(lactose_ion_area) %>% rep(each = 9)
    # print(lactose_blank_area)
    df[(nrow(df) - 17):nrow(df), fst_col + 1] = df[(nrow(df) - 17):nrow(df), fst_col + 1] - lactose_blank_areas
    
    galactitol_blank_areas = df %>% filter(grepl("Std_0", Name)) %>% pull(galactitol_ion_area) %>% rep(each = 9)
    # print(galactitol_blank_area)
    df[(nrow(df) - 17):nrow(df), fst_col + 2] = df[(nrow(df) - 17):nrow(df), fst_col + 2] - galactitol_blank_areas
    
    galactonate_blank_areas = df %>% filter(grepl("Std_0", Name)) %>% pull(galactonate_ion_area) %>% rep(each = 9)
    # print(galactonate_blank_area)
    df[(nrow(df) - 17):nrow(df), fst_col + 3] = df[(nrow(df) - 17):nrow(df), fst_col + 3] - galactonate_blank_areas
    
  } else {
    
    # Subtract matrix blank signal from standard addition signals
    galactose_blank_area = df[nrow(df) - 8, fst_col]
    # print(galactose_blank_area)
    for (i in (nrow(df)-8) : nrow(df)) {
      df[i, "galactose_ion_area"] = df[i, "galactose_ion_area"] - galactose_blank_area
    }
    
    lactose_blank_area = df[nrow(df) - 8, fst_col+1]
    # print(lactose_blank_area)
    for (i in (nrow(df)-8) : nrow(df)) {
      df[i, "lactose_ion_area"] = df[i, "lactose_ion_area"] - lactose_blank_area
    }
    
    galactitol_blank_area = df[nrow(df) - 8, fst_col+2]
    # print(galactitol_blank_area)
    for (i in (nrow(df)-8) : nrow(df)) {
      df[i, "galactitol_ion_area"] = df[i, "galactitol_ion_area"] - galactitol_blank_area
    }
    
    galactonate_blank_area = df[nrow(df) - 8, fst_col+3]
    # print(galactonate_blank_area)
    for (i in (nrow(df)-8) : nrow(df)) {
      df[i, "galactonate_ion_area"] = df[i, "galactonate_ion_area"] - galactonate_blank_area
    }
    
  }
  
  
  return(df)
  
}

compute_standard_specific_model = function(df, target_compound, r2_output_file) {
  
  # Define standard concentrations of target compound
  standard_concentrations = df %>%
    filter(grepl("Std", Name)) %>%
    pull(paste0(target_compound, "_concentration_ug_per_100ul")) %>%
    unique() %>%
    na.omit()
  # print(standard_concentrations)
  # Define standard response values of target compound
  # Model for samples using standard 1
  standard_responses = df %>% 
    filter(grepl("Std_[012345678]_1", Name)) %>%
    pull(paste0(target_compound, "_ion_area"))
  
  # Check if we need to remove any point
  standard_responses_lag = lag(standard_responses)
  standard_responses_lag[1] = 0
  diffs = standard_responses - standard_responses_lag
  idx_to_rm = which(diffs < 0)
  if (length(idx_to_rm) > 0) {
    standard_concentrations = standard_concentrations[-idx_to_rm]
    standard_responses = standard_responses[-idx_to_rm]
  }
  
  sample_responses = df %>%
    filter(!grepl("Std", Name)) %>%
    filter(closest_standard == 1) %>%
    pull(paste0(target_compound, "_ion_area"))
  
  quantified_values = generate_standard_curve(standard_concentrations, standard_responses, sample_responses, r2_output_file)
  # Negative values are converted to zero
  quantified_values[quantified_values < 0] = 0
  
  col_idx_to_update = which(colnames(df) == paste0(target_compound, "_concentration_ug_per_100ul"))
  df[!grepl("Std", df$Name) & df$closest_standard == 1, col_idx_to_update] = quantified_values
  
  # Model for samples using standard 2
  standard_concentrations = df %>%
    filter(grepl("Std", Name)) %>%
    pull(paste0(target_compound, "_concentration_ug_per_100ul")) %>%
    unique() %>%
    na.omit()
  
  standard_responses = df %>% 
    filter(grepl("Std_[012345678]_2", Name)) %>%
    pull(paste0(target_compound, "_ion_area"))
  
  # Check if we need to remove any point
  standard_responses_lag = lag(standard_responses)
  standard_responses_lag[1] = 0
  diffs = standard_responses - standard_responses_lag
  idx_to_rm = which(diffs < 0)
  if (length(idx_to_rm) > 0) {
    standard_concentrations = standard_concentrations[-idx_to_rm]
    standard_responses = standard_responses[-idx_to_rm]
  }
  
  sample_responses = df %>%
    filter(!grepl("Std", Name)) %>%
    filter(closest_standard == 2) %>%
    pull(paste0(target_compound, "_ion_area"))
  
  quantified_values = generate_standard_curve(standard_concentrations, standard_responses, sample_responses, r2_output_file)
  quantified_values[quantified_values < 0] = 0
  
  df[!grepl("Std", df$Name) & df$closest_standard == 2, col_idx_to_update] = quantified_values
  
  return(df)
  
}

get_closest_standard = function(df) {
  
  standard_info = df %>% 
    filter(grepl("Std", Name)) %>% 
    arrange(injection_order, decreasing = FALSE) %>% 
    select(injection_order) %>%
    mutate(
      standard_order = rep(1:2, each = 6),
      midpoint = c(
        rep(mean(injection_order[1:6]), 6),
        rep(mean(injection_order[7:12]), 6)
      )
    )
  df = df %>% 
    left_join(standard_info, by = "injection_order")
  
  df$closest_standard = sapply(1:nrow(df), function(x) {
    which.min(abs(df$injection_order[x] - unique(standard_info$midpoint)))
  })
  
  out = df %>%
    select(batch, Name, closest_standard, injection_order)
  
  return(out)
  
}

calculate_sample_concentrations = function(df, r2_output_file, custom_standard_map = NULL) {

    if (!is.null(custom_standard_map)) {
      
        df = df %>%
          inner_join(custom_standard_map, by = "Name")
      
        # Figure out closest standard curve using injection orders
        # standard_info = df %>% 
        #     filter(grepl("Std", Name)) %>% 
        #     arrange(injection_order, decreasing = FALSE) %>% 
        #     select(injection_order) %>%
        #     mutate(
        #         standard_order = rep(1:2, each = 6),
        #         midpoint = c(
        #             rep(mean(injection_order[1:6]), 6),
        #             rep(mean(injection_order[7:12]), 6)
        #         )
        #     )
        # df = df %>% 
        #     left_join(standard_info, by = "injection_order")
        # 
        # df$closest_standard = sapply(1:nrow(df), function(x) {
        #     which.min(abs(df$injection_order[x] - unique(standard_info$midpoint)))
        # })
        
        df = compute_standard_specific_model(df, "galactose", r2_output_file = r2_output_file)
        df = compute_standard_specific_model(df, "lactose", r2_output_file = r2_output_file)
        df = compute_standard_specific_model(df, "galactitol", r2_output_file = r2_output_file)
        df = compute_standard_specific_model(df, "galactonate", r2_output_file = r2_output_file)
      
    } else {
      
        # Galactose
        # Define standard concentrations of target compound
        standard_concentrations = df$galactose_concentration_ug_per_100ul[(nrow(df)-5):nrow(df)]
        # print(standard_concentrations)
        # Define standard response values of target compound
        standard_responses = df$galactose_ion_area[(nrow(df)-5):nrow(df)]
        # print(standard_responses)
        # Define sample response values to be quantified
        sample_responses = df$galactose_ion_area[1:(nrow(df)-6)]
        # print(sample_responses)
        # Run function
        quantified_values = generate_standard_curve(standard_concentrations, standard_responses, sample_responses)
        # print(quantified_values) 
        #Update data
        df$galactose_concentration_ug_per_100ul[1:(nrow(df)-6)] = quantified_values
        
        # Lactose
        # Define standard concentrations of target compound
        standard_concentrations = df$lactose_concentration_ug_per_100ul[(nrow(df)-5):nrow(df)]
        # print(standard_concentrations)
        # Define standard response values of target compound
        standard_responses = df$lactose_ion_area[(nrow(df)-5):nrow(df)]
        # print(standard_responses)
        # Define sample response values to be quantified
        sample_responses = df$lactose_ion_area[1:(nrow(df)-6)]
        # print(sample_responses)
        # Run function
        quantified_values = generate_standard_curve(standard_concentrations, standard_responses, sample_responses)
        # print(quantified_values) 
        #Update data
        df$lactose_concentration_ug_per_100ul[1:(nrow(df)-6)] = quantified_values
        
        # Galactitol
        # Define standard concentrations of target compound
        standard_concentrations = df$galactitol_concentration_ug_per_100ul[(nrow(df)-5):nrow(df)]
        # print(standard_concentrations)
        # Define standard response values of target compound
        standard_responses = df$galactitol_ion_area[(nrow(df)-5):nrow(df)]
        # print(standard_responses)
        # Define sample response values to be quantified
        sample_responses = df$galactitol_ion_area[1:(nrow(df)-6)]
        # print(sample_responses)
        # Run function
        quantified_values = generate_standard_curve(standard_concentrations, standard_responses, sample_responses)
        # print(quantified_values) 
        #Update data
        df$galactitol_concentration_ug_per_100ul[1:(nrow(df)-6)] = quantified_values
        
        # Galactonate
        # Define standard concentrations of target compound
        standard_concentrations = df$galactonate_concentration_ug_per_100ul[(nrow(df)-5):nrow(df)]
        # print(standard_concentrations)
        # Define standard response values of target compound
        standard_responses = df$galactonate_ion_area[(nrow(df)-5):nrow(df)]
        # print(standard_responses)
        # Define sample response values to be quantified
        sample_responses = df$galactonate_ion_area[1:(nrow(df)-6)]
        # print(sample_responses)
        # Run function
        quantified_values = generate_standard_curve(standard_concentrations, standard_responses, sample_responses)
        # print(quantified_values) 
        #Update data
        df$galactonate_concentration_ug_per_100ul[1:(nrow(df)-6)] = quantified_values
      
    }

    return(df)

}

# Create function to generate and plot standard addition curves, and calculate unknown concentrations 
generate_standard_curve = function(standard_concentrations, standard_responses, sample_responses, r2_output_file) {
  # Combine standard concentrations and responses into a data frame
  df = data.frame(standard_concentrations,standard_responses)
  # Fit linear model
  standard_curve = lm(standard_responses ~ standard_concentrations, data = df)
  # print(summary(standard_curve))
  # Calculate concentrations of unknowns from sample responses
  calculated_concentrations = (sample_responses - coef(standard_curve)[1]) / coef(standard_curve)[2]
  r_squared = summary(standard_curve)$adj.r.squared
  if (!file.exists(r2_output_file)) {
    r_square_df = data.frame(r2 = r_squared)
  } else {
    r_square_df = read.delim(r2_output_file)
    tmp_df = data.frame(r2 = r_squared)
    r_square_df = rbind(r_square_df, tmp_df)
  }
  write.table(
    r_square_df,
    r2_output_file,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  # TODO add calculated points in plot
  # Generate plot
  plot_range = range(standard_concentrations) + c(-1, 1)
  p = ggplot(df, aes(x = standard_concentrations, y = standard_responses)) +
    geom_point(color = "blue") +
    geom_point(data = data.frame(calculated_concentrations = calculated_concentrations, sample_responses = sample_responses), aes(x = calculated_concentrations, y = sample_responses)) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    labs(title = paste0("Standard Addition Curve", r_squared),
         x = "Added Concentration",
         y = "Response") +
    xlim(plot_range)
  # print(p)
  # Return calculated concentrations
  return(calculated_concentrations)
}

normalize_by_urine_volume = function(df) {

    # Create "timepoint traces" for timepoint-specific correction by whole urine volume
    # b1 = grepl("AH", df$Name)
    b2 = grepl("BL", df$Name)
    p1 = grepl("0-3", df$Name)
    p2 = grepl("3-6", df$Name)

    #Galactose
    df$galactose_ug = NA_real_ #Create new column for absolute quantities
    # df$galactose_ug[b1] = df$galactose_concentration_ug_per_100ul[b1] * df$urine_b2_ml[b1] / 10 # shouldn't this be * 10 instead of / 10?
    df$galactose_ug[b2] = df$galactose_concentration_ug_per_100ul[b2] * df$urine_b2_ml[b2] / 10
    df$galactose_ug[p1] = df$galactose_concentration_ug_per_100ul[p1] * df$urine_p1_ml[p1] / 10
    df$galactose_ug[p2] = df$galactose_concentration_ug_per_100ul[p2] * df$urine_p2_ml[p2] / 10

    #Lactose
    df$lactose_ug = NA_real_
    df$lactose_ug[b2] = df$lactose_concentration_ug_per_100ul[b2] * df$urine_b2_ml[b2] / 10
    df$lactose_ug[p1] = df$lactose_concentration_ug_per_100ul[p1] * df$urine_p1_ml[p1] / 10
    df$lactose_ug[p2] = df$lactose_concentration_ug_per_100ul[p2] * df$urine_p2_ml[p2] / 10

    #Galactitol
    df$galactitol_ug = NA_real_
    df$galactitol_ug[b2] = df$galactitol_concentration_ug_per_100ul[b2] * df$urine_b2_ml[b2] / 10
    df$galactitol_ug[p1] = df$galactitol_concentration_ug_per_100ul[p1] * df$urine_p1_ml[p1] / 10
    df$galactitol_ug[p2] = df$galactitol_concentration_ug_per_100ul[p2] * df$urine_p2_ml[p2] / 10

    #Galactonate
    df$galactonate_ug = NA_real_
    df$galactonate_ug[b2] = df$galactonate_concentration_ug_per_100ul[b2] * df$urine_b2_ml[b2] / 10
    df$galactonate_ug[p1] = df$galactonate_concentration_ug_per_100ul[p1] * df$urine_p1_ml[p1] / 10
    df$galactonate_ug[p2] = df$galactonate_concentration_ug_per_100ul[p2] * df$urine_p2_ml[p2] / 10

    return(df)

}


normalize_by_signal_per_volume <- function(df, per_n_ml = 1, signal_suffix = "_ion_area") {

  b2 <- grepl("BL", df$Name)
  p1 <- grepl("0-3", df$Name)
  p2 <- grepl("3-6", df$Name)
  
  volume_ml <- rep(NA_real_, nrow(df))
  volume_ml[b2] <- df$urine_b2_ml[b2]
  volume_ml[p1] <- df$urine_p1_ml[p1]
  volume_ml[p2] <- df$urine_p2_ml[p2]
  
  # Find every raw signal column to normalize (e.g. galactose_ion_area, lactose_ion_area, ...)
  signal_cols <- names(df)[endsWith(names(df), signal_suffix)]
  
  for (col in signal_cols) {
    compound <- sub(signal_suffix, "", col)
    new_col <- paste0(compound, "_per_", per_n_ml, "ml")
    df[[new_col]] <- df[[col]] / volume_ml * per_n_ml
  }
  df
}

prepare_and_normalize <- function(df, volume_data, signal_suffix, per_n_ml = 1, id_pattern = "^LB[0-9]+") {
  
  # Derive study_id_lb from sample Name
  df <- df %>%
    mutate(study_id_lb = if_else(
      Type == "Sample",
      str_extract(Name, id_pattern),
      NA_character_
    ))
  
  # Keep only actual samples -- QC/Std/Blank rows have no urine volume
  # and shouldn't go through volume normalization
  df <- df %>%
    filter(Type == "Sample")
  
  # Join volume/metadata onto compound data
  df <- df %>%
    left_join(volume_data, by = "study_id_lb")
  
  # Calculate whole urine volumes (assumed urine density ~1.02 g/ml)
  df <- df %>%
    mutate(
      urine_b2_ml = (int_wt_urine_b2 - int_wt_empty_urine_b2) / 1.02,
      urine_p1_ml = (int_urine_3_00_wt_pool1 - int_urine_wt_empty_pool1) / 1.02,
      urine_p2_ml = (int_urine_6_00_wt_pool2 - int_urine_wt_empty_pool2) / 1.02
    )
  
  # Calculate fluid consumption (water density ~1.00 g/ml)
  df <- df %>%
    mutate(
      water_consumed_b2_ml = 650,
      water_consumed_p1_ml = int_water_wt_full - int_water_wt_3_00,
      water_consumed_p2_ml = int_water_wt_full - int_water_wt_6_00
    )
  
  # Normalize signal by urine volume
  df <- normalize_by_signal_per_volume(df, per_n_ml = per_n_ml, signal_suffix = signal_suffix)
  
  df
}

plot_by_compound = function(comp) {

  df = samples_cr_clean %>% filter(compound == comp)
  y_label = gsub("_ug$", " / µg", comp)
  y_label = paste0(toupper(substr(y_label, 1, 1)), substr(y_label, 2, nchar(y_label)))
  
  ggplot(
    df,
    aes(
      x = timepoint,
      y = value,
      group = interaction(study_id, treatment),
      linetype = treatment,
      color = screen_group,
      shape = screen_group
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    xlab("Timepoints") +
    ylab(y_label) +
    scale_color_manual(
      name = "Group",
      values = c(
        "Control" = "#984EA3FF",
        "LM without symptoms" = "#a7c193",
        "LM with symptoms"   = "#ED685E"
      )
    ) +
    scale_shape_manual(
      name = "Group",
      values = c(
        "Control" = 15,
        "LM without symptoms" = 17,
        "LM with symptoms"   = 16
      )
    ) +
    theme_consistent
}

plot_box_by_compound = function(comp) {

  df = samples_cr_clean %>% filter(compound == comp)
  y_label = gsub("_ug$", " / µg", comp)
  y_label = paste0(toupper(substr(y_label, 1, 1)), substr(y_label, 2, nchar(y_label)))
  
  ggplot(df, aes(x = timepoint, y = value, fill = screen_group)) +
    geom_boxplot(outlier.shape = 21, color = "black", width = 0.7) +
    xlab("Timepoints") +
    ylab(y_label) +
    scale_fill_manual(
      name = "Group",
      values = c(
        "Control" = "#984EA3FF",
        "LM without symptoms" = "#a7c193",
        "LM with symptoms"   = "#ED685E"
      )
    ) +
    facet_wrap(treatment ~ .) +
    theme_consistent
}

build_posthoc_plot = function(data, posthoc_results, test, sig_thr = 0.05, plot_col = "screen_group") {

  if (test == "timewise") {
    
    sig_res = posthoc_results %>%
      filter(p_adj <= sig_thr) %>%
      mutate(plot_col = paste0(compound, ":", timepoint, ":", treatment))
    
    i = 1
    plots_ls = list()
    comparison_ls = list()
    for(group in unique(sig_res$plot_col)) {
      
      if (plot_col == "screen_group") {
        
        group_compound = sapply(strsplit(group, ":"), "[[", 1)
        group_timepoint = sapply(strsplit(group, ":"), "[[", 2)
        group_treatment = sapply(strsplit(group, ":"), "[[", 3)
        
        data_sub = data %>%
          filter(compound == group_compound) %>%
          filter(timepoint == group_timepoint) %>%
          filter(treatment == group_treatment) %>%
          mutate(screen_group = factor(screen_group))
        
        max_y = max(data_sub$value)
        padding = 0.1 * max_y
        
        comparison_ls[[i]] = sig_res %>%
          filter(plot_col == group) %>%
          separate(comparison, c("g1", "g2"), sep = " - ") %>%
          mutate(
            nrow = 1:length(compound),
            g1 = factor(g1, levels = c("Control", "LM-A", "LM-S")),
            g2 = factor(g2, levels = c("Control", "LM-A", "LM-S")),
            g1 = as.numeric(g1),
            g1 = case_when(
              g1 == 1 ~ -0.3,
              g1 == 2 ~ 0,
              g1 == 3 ~ 0.3
            ),
            g2 = as.numeric(g2),
            g2 = case_when(
              g2 == 1 ~ -0.3,
              g2 == 2 ~ 0,
              g2 == 3 ~ 0.3
            ),
            xlabel = as.numeric(timepoint) + ((g1 + g2) / 2),
            xstart = as.numeric(timepoint) + g1,
            xend = as.numeric(timepoint) + g2,
            yval = max_y + padding * nrow,
            ylabel = yval + 0.01 * yval
          )
        
      } else {
        group_compound = sapply(strsplit(group, ":"), "[[", 1)
        group_timepoint = sapply(strsplit(group, ":"), "[[", 2)
        group_treatment = sapply(strsplit(group, ":"), "[[", 3)
        
        data_sub = data %>%
          filter(compound == group_compound) %>%
          filter(timepoint == group_timepoint) %>%
          filter(treatment == group_treatment) %>%
          mutate(gen_test_rs4988235 = factor(gen_test_rs4988235))
        
        max_y = max(data_sub$value)
        padding = 0.1 * max_y
        
        comparison_ls[[i]] = sig_res %>%
          filter(plot_col == group) %>%
          separate(comparison, c("g1", "g2"), sep = " - ") %>%
          mutate(
            nrow = 1:length(compound),
            g1 = factor(g1, levels = c("TT", "CT", "CC")),
            g2 = factor(g2, levels = c("TT", "CT", "CC")),
            g1 = as.numeric(g1),
            g1 = case_when(
              g1 == 1 ~ -0.3,
              g1 == 2 ~ 0,
              g1 == 3 ~ 0.3
            ),
            g2 = as.numeric(g2),
            g2 = case_when(
              g2 == 1 ~ -0.3,
              g2 == 2 ~ 0,
              g2 == 3 ~ 0.3
            ),
            xlabel = as.numeric(timepoint) + ((g1 + g2) / 2),
            xstart = as.numeric(timepoint) + g1,
            xend = as.numeric(timepoint) + g2,
            yval = max_y + padding * nrow,
            ylabel = yval + 0.01 * yval
          )
      }
      
      i = i + 1
      
    }
    
  } else {
    
    sig_res = posthoc_results %>%
      filter(p_adj <= sig_thr) %>%
      mutate(plot_col = paste0(compound, ":", !!sym(plot_col), ":", treatment))
    
    i = 1
    plots_ls = list()
    comparison_ls = list()
    for(group in unique(sig_res$plot_col)) {
      
      if (plot_col == "screen_group") {
        
        group_compound = sapply(strsplit(group, ":"), "[[", 1)
        group_screen_group = sapply(strsplit(group, ":"), "[[", 2)
        group_treatment = sapply(strsplit(group, ":"), "[[", 3)
        
        data_sub = data %>%
          filter(compound == group_compound) %>%
          filter(screen_group == group_screen_group) %>%
          filter(treatment == group_treatment) %>%
          mutate(
            timepoint = factor(timepoint)
          )
        
        max_y = max(data_sub$value)
        padding = 0.1 * max_y
        
        comparison_ls[[i]] = sig_res %>%
          filter(plot_col == group) %>%
          separate(comparison, c("t1", "t2"), sep = " - ") %>%
          mutate(
            nrow = 1:length(compound),
            t1 = factor(t1, levels = c("Baseline", "0-3h", "3-6h")),
            t2 = factor(t2, levels = c("Baseline", "0-3h", "3-6h")),
            g1 = as.numeric(t1),
            g1 = case_when(
              g1 == 1 ~ -0.3,
              g1 == 2 ~ 0,
              g1 == 3 ~ 0.3
            ),
            g2 = as.numeric(t2),
            g2 = case_when(
              g2 == 1 ~ -0.3,
              g2 == 2 ~ 0,
              g2 == 3 ~ 0.3
            ),
            xlabel = as.numeric(screen_group) + ((g1 + g2) / 2),
            xstart = as.numeric(screen_group) + g1,
            xend = as.numeric(screen_group) + g2,
            yval = max_y + padding * nrow,
            ylabel = yval + 0.01 * yval
          )
        
      } else {
        
        group_compound = sapply(strsplit(group, ":"), "[[", 1)
        group_genotype = sapply(strsplit(group, ":"), "[[", 2)
        group_treatment = sapply(strsplit(group, ":"), "[[", 3)
        
        data_sub = data %>%
          filter(compound == group_compound) %>%
          filter(gen_test_rs4988235 == group_genotype) %>%
          filter(treatment == group_treatment) %>%
          mutate(
            timepoint = factor(timepoint)
          )
        
        max_y = max(data_sub$value)
        padding = 0.1 * max_y
        
        comparison_ls[[i]] = sig_res %>%
          filter(plot_col == group) %>%
          separate(comparison, c("t1", "t2"), sep = " - ") %>%
          mutate(
            nrow = 1:length(compound),
            t1 = factor(t1, levels = c("Baseline", "0-3h", "3-6h")),
            t2 = factor(t2, levels = c("Baseline", "0-3h", "3-6h")),
            g1 = as.numeric(t1),
            g1 = case_when(
              g1 == 1 ~ -0.3,
              g1 == 2 ~ 0,
              g1 == 3 ~ 0.3
            ),
            g2 = as.numeric(t2),
            g2 = case_when(
              g2 == 1 ~ -0.3,
              g2 == 2 ~ 0,
              g2 == 3 ~ 0.3
            ),
            xlabel = as.numeric(gen_test_rs4988235) + ((g1 + g2) / 2),
            xstart = as.numeric(gen_test_rs4988235) + g1,
            xend = as.numeric(gen_test_rs4988235) + g2,
            yval = max_y + padding * nrow,
            ylabel = yval + 0.01 * yval
          )
        
      }
      
      i = i + 1
      
    }
    
  }
  
  

  # return(plots_ls)
  return(comparison_ls)

}
