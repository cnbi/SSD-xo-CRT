################################### GET POWER ##############################

get_power <- function(eff_size, cac, wp_icc, bp_icc, periods, n2, n1, treatment_n, 
                      pattern, seed, BF_thresh, batch_size, 
                      ndatasets, gc_every = 10){
    
    # Libraries
    library(lme4)
    library(bain)
    # library(BFpack)
    
    # Warnings
    if (is.numeric(c(eff_size, n1, n2, ndatasets, BF_thresh, batch_size)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    if (eff_size < 0) stop("The effect size must be a positive value ")
    
    # Functions
    source("data_generation_xo.R")
    source("helper_functions.R")
    source("get_bf.R")
    source("print_results.R")
    
    # Initial values
    n2_seq <- n2/2 #Number of clusters per sequence
    singular_warn <- 0

    
    results_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)
    final_SSD <- vector(mode = "list")
    
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
    current_eta <- length(which(results_H1[, "BF.1c"] > BF_thresh)) / ndatasets
    
    result <- list("n1" = n1,
                   "n2" = n2,
                   "Proportion.BF1c" = current_eta)
    
    # Display results
    print_results(list_results = result, hypotheses_set = 2, 
                  BF_thresh = BF_thresh)
}

# get_power(eff_size = 0.4, cac = 0.8, wp_icc = 0.01, periods = 2, n2 = 60,
#                           n1 = 10, treatment_n = 2,
#                           pattern = "chess", seed = 23, BF_thresh = 2,
#                           batch_size = 10,
#                           ndatasets = 10)
