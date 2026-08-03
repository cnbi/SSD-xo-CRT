################# SAMPLE SIZE DETERMINATION FOR CROSSOVER DESIGN ##############




ssd_crt_xo <- function(eff_size, cac, wp_icc, bp_icc, periods, n2, n1, treatment_n, 
                       pattern, seed, BF_thresh, eta, batch_size, max_sample, 
                       ndatasets, fixed, gc_every = 10) {
    
    # Libraries
    library(lme4)
    library(bain)
    # library(BFpack)
    
    # Warnings
    if (is.numeric(c(eff_size, n1, n2, ndatasets, BF_thresh, max_sample, batch_size)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    if (eff_size < 0) stop("The effect size must be a positive value ")
    if (eta > 1) stop("The probability of exceeding Bayes Factor threshold cannot be larger than 1")
    if (is.character(fixed) == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    
    # Functions
    source("data_generation_xo.R")
    source("helper_functions.R")
    source("get_bf.R")
    source("print_results.R")
    
    # Initial values
    n2_seq <- n2/2 #Number of clusters per sequence
    previous_high <- 0
    previous_eta <- 0
    current_eta <- 0
    singular_warn <- 0
    ultimate_sample_size <- FALSE
    
    if (fixed == "n1") {
        min_sample <- 6                     # Minimum number of clusters
        low <- min_sample                   #lower bound
    } else if (fixed == "n2") {
        min_sample <- 4                     # Minimum cluster size
        low <- min_sample                   #lower bound
    }
    high <- max_sample                    #higher bound
    
    results_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)
    final_SSD <- vector(mode = "list")
    
    while (ultimate_sample_size == FALSE) {
        # Generate data ------------------------------------------------------------
        if (missing(bp_icc)) {
            data_H1 <- do.call(data_generation, list(eff_size = eff_size, 
                                                     cac = cac, wp_icc = wp_icc,
                                                     periods = periods, n2 = n2, 
                                                     n1 = n1,
                                                     treatment_n = treatment_n, 
                                                     seed = seed, 
                                                     ndatasets = ndatasets, 
                                                     batch_size = batch_size,
                                                     gc_every = gc_every))
        } else {
            data_H1 <- do.call(data_generation, list(eff_size = eff_size,
                                                     bp_icc = bp_icc,
                                                     wp_icc = wp_icc,
                                                     periods = periods, n2 = n2, 
                                                     n1 = n1,
                                                     treatment_n = treatment_n, 
                                                     seed = seed,
                                                     ndatasets = ndatasets, 
                                                     batch_size = batch_size,
                                                     gc_every = gc_every))
        }
        
        
        # If H0 is true
        # data_H0 <- do.call(data_generation, list(eff_size = 0, cac,wp_icc, periods, n2, 
        #                                          n1, treatment_n, 
        #                                          seed, ndatasets, 
        #                                          batch_size = batch_size))
        
        # Compute Bayes factor --------------------------------------------------
        n_eff <- ((n1 * n2) / (1 + (n1 - 1) * data_H1$emp_wp_icc)) / 2 # Effecttive sample size
        output_bf <- Map(get_BF, n_eff, data_H1$treatment_eff,  data_H1$var_treat, 2)
        
        # Evaluate power criterion ----------------------------------------------
        results_H1[, 1] <- unlist(lapply(output_bf, `[`, "bf.12")) # Bayes factor H1vsHc
        results_H1[, 2] <- unlist(lapply(output_bf, `[`, "pmp1")) #posterior model probabilities of H1
        results_H1[, 3] <- unlist(lapply(output_bf, `[`, "bf.21")) # Bayes factor HcvsH1
        results_H1[, 4] <- unlist(lapply(output_bf, `[`, "pmp2")) #posterior model probabilities of Hc
        
        # results_H1[, 1] <- unlist(lapply(output_bf, extract_results, 1)) # Bayes factor H1vsHc
        # results_H1[, 2] <- unlist(lapply(output_bf, extract_results, 4)) #posterior model probabilities of H1
        # results_H1[, 3] <- unlist(lapply(output_bf, extract_results, 2)) # Bayes factor HcvsH1
        # results_H1[, 4] <- unlist(lapply(output_bf, extract_results, 3)) #posterior model probabilities of Hc
        
        colnames(results_H1) <- c("BF.1c", "PMP.1", "BF.c1", "PMP.c")
        
        # Proportion
        previous_eta <- current_eta
        current_eta <- length(which(results_H1[, "BF.1c"] > BF_thresh)) / ndatasets
        
        # Evaluation
        ifelse(current_eta > eta, condition_met <- TRUE, condition_met <- FALSE)
        
        # Update sample size ---------------------------------------------------
        binary_search <- update_sample(results = results_H1, eta = eta, 
                                       current_eta = current_eta, 
                                       n2 = n2, n1 = n1, 
                                       condition_met = condition_met,
                                       fixed = fixed, hypotheses_set = 2,
                                       results_H1 = results_H1, data_H1 = data_H1,
                                       final_SSD = final_SSD, 
                                       min_sample = min_sample, 
                                       max_sample = max_sample, 
                                       singular_warn = singular_warn,
                                       low = low, high = high, 
                                       previous_high = previous_high, 
                                       previous_eta = previous_eta, 
                                       ultimate_sample_size = ultimate_sample_size)
        low <- binary_search$low
        high <- binary_search$high
        n1 <- binary_search$n1
        n2 <- binary_search$n2
        ultimate_sample_size <- binary_search$ultimate_sample_size
        previous_high <- binary_search$previous_high
        
        rm(data_H1)
        gc(full = TRUE)
    }
    
    rm(output_bf)
    if (n2 < 30) warning("The number of groups is less than 30.
                                             This may cause problems in convergence and singularity.")
    
    # Display results
    print_results(list_results = binary_search$final_SSD, hypotheses_set = 2, 
                  BF_thresh = BF_thresh)
    if (any(singular_warn > 0)) warning("At least one of the fitted models is singular. For more information about singularity see help('isSingular').
                               The number of models that are singular can be found in the output object.")
    invisible(binary_search$final_SSD)
}

# #Test
# a <- ssd_crt_xo(eff_size = 0.4, cac = 0.8, wp_icc = 0.01, periods = 2, n2 = 60,
#                 n1 = 10, treatment_n = 2,
#                 pattern = "chess", seed = 23, BF_thresh = 2, eta = 0.8,
#                 batch_size = 10, max_sample = 80,
#                 ndatasets = 10, fixed = "n1")
