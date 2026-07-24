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
                   fixed = design_matrix[row, "fixed"])
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
                   fixed = design_matrix[row, "fixed"])
    }
    
    # End time and save results
    end_time <- Sys.time()
    file_name <- file.path(results_folder,
                           paste0("Results", finding, "Row", Row, ".RDS"))
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