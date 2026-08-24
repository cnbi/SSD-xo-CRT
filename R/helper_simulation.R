#################### HELPERS IN SIMULATION #########################

run_sim <- function(row, design_matrix,
                    ndatasets,
                    Max,
                    batch_size,
                    results_folder) {
    
    fixed <- design_matrix[row, "fixed"]
    fixed <- toupper(fixed)
    if (fixed == "N1") {
        finding <- "N2"
    } else if (fixed == "N2") {
        finding <- "N1"
    }
    
    # Start time
    start_time <- Sys.time()
    
    if (fixed == "N1") {
        ssd_results <- ssd_crt_xo(eff_size = design_matrix[row, "eff_size"],
                                  cac = design_matrix[row, "cac"],
                                  wp_icc = design_matrix[row, "wp_icc"],
                                  periods = design_matrix[row, "periods"],
                                  n2 = design_matrix[row, "n2"],
                                  n1 = design_matrix[row, "n1"],
                                  treatment_n = 2, pattern = "chess",
                                  seed = design_matrix[row, "seed"],
                                  BF_thresh = design_matrix[row, "BF_thresh"],
                                  eta = design_matrix[row, "eta"],
                                  batch_size = batch_size, max_sample = Max, ndatasets = ndatasets,
                                  fixed = as.character(design_matrix[row, "fixed"]))
    } else if (fixed == "N2") {
        ssd_results <- ssd_crt_xo(eff_size = design_matrix[row, "eff_size"],
                                  cac = design_matrix[row, "cac"],
                                  wp_icc = design_matrix[row, "wp_icc"],
                                  periods = design_matrix[row, "periods"],
                                  n2 = design_matrix[row, "n2"],
                                  n1 = design_matrix[row, "n1"],
                                  treatment_n = 2, pattern = "chess",
                                  seed = design_matrix[row, "seed"],
                                  BF_thresh = design_matrix[row, "BF_thresh"],
                                  eta = design_matrix[row, "eta"],
                                  batch_size = batch_size, max_sample = Max, ndatasets = ndatasets,
                                  fixed = as.character(design_matrix[row, "fixed"]))
    }
    
    # End time and save results
    end_time <- Sys.time()
    file_name <- file.path(results_folder,
                           paste0("Results", finding, "Row", row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    
    # Save running time
    running_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("time", finding, "Row", row, ".RDS"))
    saveRDS(running_time, file = time_name)
    
    # Clean up memory
    rm(ssd_results)
    NULL
    gc()
}

# Check simulation--------------------
missing_rows <- function(folder_path,
                         name_pattern = NULL,
                         check_numbers,
                         underscore = TRUE) {
    files_names <- list.files(folder_path)
    
    # Filter by name pattern
    if (!is.null(name_pattern)) {
        filtered_names <- files_names[grep(name_pattern, files_names)]
    }
    
    # Extract the number of row
    if (underscore) {
        row_numbers <- sapply(filtered_names, function(names){
            parts <- unlist(strsplit(names, "_"))
            number <- parts[length(parts)]
            
            # Remove extension
            number_only <- sub("\\.[^.]+$", "", number)
            return(number_only)
        })   
    } else if (underscore == FALSE) {
        
        row_numbers <- sapply(filtered_names, function(names){
            # Remove extension
            name_no_ext <- tools::file_path_sans_ext(names)
            
            number <- regmatches(name_no_ext, gregexpr("\\d+$", name_no_ext))[[1]]
            
            return(number)
        })
    }
    
    
    difference <- setdiff(check_numbers, row_numbers)
    print(difference)
}

# Collect results -------------------------------------------
collect_results <- function(design_matrix,
                            results_folder,
                            pair,
                            results_name,
                            b = 1, rows, 
                            save = TRUE,
                            file_name) {
    if (missing(rows)) {
        rows <-  seq(nrow(design_matrix))
    }
    
    # Create results' name
    fixed <- design_matrix[, "fixed"]
    fixed <- toupper(fixed)
    finding <- ifelse(fixed == "N1", "N2", "N1")
    results_name_full <- paste0(results_name, finding, "Row", rows, ".RDS")
    
    if (pair == 2) {
        # When evaluating set 2 of hypotheses (H1vsH2)
        new_matrix <- matrix(NA, ncol = 7, nrow = nrow(design_matrix))
        for (row_design in rows) {
            # Open file and extract the information
            
            stored_result <- readRDS(file.path(results_folder, results_name_full[row_design]))
            median.BF1c <- median(stored_result[[4]][, "BF.1c"], na.rm = TRUE)
            median.BFc1 <- median(stored_result[[4]][, "BF.c1"], na.rm = TRUE)
            mean.PMP1 <- mean(stored_result[[4]][, "PMP.1"], na.rm = TRUE)
            mean.PMPc <- mean(stored_result[[4]][, "PMP.c"], na.rm = TRUE)
            eta.BF12 <- stored_result$Proportion.BF1c
            n2 <- stored_result$n2
            n1 <- stored_result$n1
            
            # Save in matrix
            new_matrix[row_design, ] <- c(mean.PMP1,
                                          median.BF1c,
                                          median.BFc1,
                                          mean.PMPc,
                                          eta.BF12,
                                          n2,
                                          n1)
        }
        new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
        colnames(new_matrix) <- c(
            names(design_matrix),
            "mean.PMP1",
            "median.BF1c",
            "median.BFc1",
            "mean.PMPc",
            "eta.BF12",
            "n2.final",
            "n1.final"
        )
        
        saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, "_set2.RDS")))
        
    } else if (pair == 1) {
        #TODO:
        # Hypothesis set 1 (H0vsH1)
        new_matrix <- matrix(NA, ncol = 11, nrow = nrow(design_matrix)*3)
        row_res <- 1
        for (row_design in rows) {
            for (b_frac in seq(b)) {
                stored_result <- readRDS(file.path(results_folder, paste0(results_name,  results_name_full[row_design])))
                median.BF01 <- median(stored_result[[b_frac]][[6]][, "BF.01"], na.rm = TRUE)       # 6: data_H0
                median.BF10 <- median(stored_result[[b_frac]][[7]][, "BF.10"], na.rm = TRUE)       # 7: data_H1
                mean.PMP0.H0 <- mean(stored_result[[b_frac]][[6]][, "PMP.0"], na.rm = TRUE)        # 6: data_H0
                mean.PMP1.H0 <- mean(stored_result[[b_frac]][[6]][, "PMP.1"], na.rm = TRUE)        # 6: data_H0
                mean.PMP0.H1 <- mean(stored_result[[b_frac]][[7]][, "PMP.0"], na.rm = TRUE)        # 7: data_H1
                mean.PMP1.H1 <- mean(stored_result[[b_frac]][[7]][, "PMP.1"], na.rm = TRUE)        # 7: data_H1
                n2 <- stored_result[[b_frac]]$n2
                eta.BF01 <- stored_result[[b_frac]]$Proportion.BF01
                eta.BF10 <- stored_result[[b_frac]]$Proportion.BF10
                n1 <- stored_result[[b_frac]]$n1
                new_matrix[row_res, ] <- c(b_frac,
                                           median.BF01,
                                           mean.PMP0.H0,
                                           mean.PMP1.H0,
                                           median.BF10,
                                           mean.PMP0.H1,
                                           mean.PMP1.H1,
                                           eta.BF01,
                                           eta.BF10,
                                           n1,
                                           n2
                )
                row_res <- row_res + 1}
        }
        design_matrix <- design_matrix %>% slice(rep(1:n(), each = b))
        new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
        colnames(new_matrix) <- c(
            names(design_matrix),
            "b",
            "median.BF01",
            "mean.PMP0.H0",
            "mean.PMP1.H0",
            "median.BF10",
            "mean.PMP0.H1",
            "mean.PMP1.H1",
            "eta.BF01",
            "eta.BF10",
            "n1.final",
            "n2.final"
        )
        if (save == TRUE) {
            saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, "_set1.RDS")))
        }
    }
    return(new_matrix)
}

# Collect times in a matrix ----
collect_times <- function(design_matrix,
                          rows,
                          pair,
                          times_name,
                          results_folder,
                          file_name) {
    if (missing(rows)) {
        rows <-  seq(nrow(design_matrix))
    }
    
    # Create results' name
    fixed <- design_matrix[, "fixed"]
    fixed <- toupper(fixed)
    finding <- ifelse(fixed == "N1", "N2", "N1")
    times_name_full <- paste0(times_name, finding, "Row", rows, ".RDS")
    
    new_matrix <- matrix(NA, nrow = nrow(design_matrix), ncol = 1)
    
    # Pair of hypotheses to compare
    if (pair == 1) {
        type <- "_set1.RDS"
    } else if (pair == 2) {
        type <- "_set2.RDS"
    }
    
    for (row_result in rows) {
        stored_result <- readRDS(file.path(results_folder, times_name_full[row_result]))
        new_matrix[row_result, ] <- stored_result
    }
    new_matrix <- as.data.frame(cbind(design_matrix, new_matrix))
    colnames(new_matrix) <- c(names(design_matrix), "total.time.min")
    saveRDS(new_matrix, file = file.path(results_folder, paste0(file_name, type)))
    
    return(new_matrix)
}
