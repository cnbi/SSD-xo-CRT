############################ HELPERS XO CRT ###################################

# Fit multilevel model ---------------------------------------------------------
fit_lmer <- function(x) {
    
    fitted_model <- tryCatch({
        suppressMessages(suppressWarnings({
            lmer(resp ~ treat_cond + period + (1|cluster) + (1|cluster:period), 
                 data = x)
        }))},  error = function(e) return(NULL))
    
    return(fitted_model)
}

# Variance-covariance fixed effects---------------------------------------------
vcov_ <- function(output_lmer, name) {
    if (is.null(output_lmer)) return(matrix(NA_real_, 1, 1))
    var_cov_fix <- tryCatch(
        suppressWarnings(vcov.merMod(output_lmer)[name, name]),
        error = function(e) {
            message("vcov() failed (non-positive definite matrix): returning NA")
            NA_real_
        }
    )
    return(var_cov_fix)
}

# Random effects variances -----------------------------------------------------
get_var <- function(output_lmer) {
    variances <- as.data.frame(VarCorr(output_lmer))
    var_g <- variances[1, 4]
    var_c <- variances[2, 4]
    var_e <- variances[3, 4]
    
    return(c("var_g" = var_g, "var_c" = var_c, "var_e" = var_e))
}

# Extract results Bayes factor -------------------------------------------------

extract_results <- function(output_bf, number) {
    results <- output_bf[[number]]
    return(results)
}

# Rounding half away from zero -------------------------------------------------
round2 <- function(number, decimals = 0) {
    sign_number <- sign(number)
    number <- abs(number) * 10^decimals
    number <- number + 0.5 + sqrt(.Machine$double.eps)
    number <- trunc(number)
    number <- number / 10 ^ decimals
    number * sign_number
}
# Source: https://stackoverflow.com/questions/66600344/commercial-rounding-in-r-i-e-always-round-up-from-5/66600470#66600470

# Binary search algorithm ------------------------------------------------------
update_sample <- function(results, eta, current_eta, n2, n1, condition_met,
                          fixed, hypotheses_set, b, prop_BF01, prop_BF10, results_H0,
                          results_H1, data_H0, data_H1, final_SSD, min_sample, 
                          max_sample, singular_warn, low, high, previous_high, 
                          previous_eta, ultimate_sample_size) {
    if (hypotheses_set == 1) {
        if (condition_met == FALSE) {
            # print(c("Using cluster size:", n1, "and number of clusters:", n2,
            #         "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10, "b:", b,
            #         "low:", low, "high:", high))
            message("Increasing sample")
            if (fixed == "n1") {
                if ((n2 == max_sample) | (n2 > max_sample))    { # If the sample size reaches the maximum
                    final_SSD[[b]] <- list("n1" = n1,
                                           "n2" = n2,
                                           "Proportion.BF01" = prop_BF01,
                                           "Proportion.BF10" = prop_BF10,
                                           "b.frac" = b,
                                           "data_H0" = results_H0,
                                           "data_H1" = results_H1,
                                           "singularity" = cbind(H0 = data_H0$singularity,
                                                                 H1 = data_H1$singularity))
                    b <- b + 1
                    low <- min_sample
                    previous_eta <- 0
                    previous_high <- 0
                    high <- max_sample
                    next
                } else {
                    # Increase the number of clusters since eta is too small
                    low <- n2                         #lower bound
                    high <- high                      #higher bound
                    n2 <- round2((low + high) / 2)     #point in the middle
                    if (n2 %% 2 != 0) n2 <- n2 + 1 # To ensure number of clusters is even
                    
                    # Adjust higher bound when there is a ceiling effect
                    if (low + n2 == high * 2) {
                        low <- n2                         #lower bound
                        if (previous_high > 0) {
                            high <- previous_high
                        } else {
                            high <- max_sample                   #higher bound
                        }
                        n2 <- round2((low + high) / 2)     #point in the middle
                    }
                }
            } else if (fixed == "n2") {
                if ((n1 == max_sample) | (n1 > max_sample))    {# If the sample size reaches the maximum
                    final_SSD[[b]] <- list("n1" = n1,
                                           "n2" = n2,
                                           "Proportion.BF01" = prop_BF01,
                                           "Proportion.BF10" = prop_BF10,
                                           "b.frac" = b,
                                           "data_H0" = results_H0,
                                           "data_H1" = results_H1,
                                           "singularity" = cbind(H0 = data_H0$singularity,
                                                                 H1 = data_H1$singularity))
                    b <- b + 1
                    low <- min_sample
                    previous_eta <- 0
                    previous_high <- 0
                    high <- max_sample
                    next
                } else {
                    # Increase the cluster sizes since eta is too small
                    low <- n1                        #lower bound
                    high <- high                     #higher bound
                    n1 <- round2((low + high) / 2)    #point in the middle
                    
                    # Adjust higher bound when there is a ceiling effect
                    if ((low + n1 == high * 2) | (current_eta == previous_eta)) {
                        low <- n1                        #lower bound
                        # Set the higher bound based on the previous high or the maximum
                        if (previous_high > 0 ) {
                            high <- previous_high
                        } else {
                            high <- max_sample
                        }
                        n1 <- round2((low + high) / 2)    #point in the middle
                    }
                }
            }
            break
        } else if (condition_met == TRUE) {
            # print(c("Using cluster size:", n1,
            #         "and number of clusters:", n2,
            #         "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10,
            #         "low: ", low, "high: ", high, "b:", b))
            previous_high <- high
            SSD_object <- list("n1" = n1,
                               "n2" = n2,
                               "Proportion.BF01" = prop_BF01,
                               "Proportion.BF10" = prop_BF10,
                               "b.frac" = b,
                               "data_H0" = results_H0,
                               "data_H1" = results_H1,
                               "singularity" = cbind(H0 = data_H0$singularity,
                                                     H1 = data_H1$singularity))
            # print("Lowerign sample")
            # print(c("previous:", previous_eta))
            previous_eta <- current_eta
            
            if (fixed == "n1") {
                # Eta is close enough to the desired eta
                if (current_eta - eta < 0.1 && n2 - low == 2) {
                    final_SSD[[b]] <- SSD_object
                    singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                    b <- b + 1
                    low <- min_sample
                    previous_eta <- 0
                    previous_high <- 0
                    high <- max_sample
                    next
                    
                } else if (previous_eta == current_eta && n2 - low == 2) {
                    # If there is no change in eta and the lower bound is close to the middle point
                    final_SSD[[b]] <- SSD_object
                    singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                    b <- b + 1
                    low <- min_sample
                    previous_eta <- 0
                    previous_high <- 0
                    high <- max_sample
                    next
                    
                } else {
                    # Decreasing to find the ultimate number of clusters
                    low <- low                         #lower bound
                    high <- n2                         #higher bound
                    n2 <- round2((low + high) / 2)      #point in the middle
                    if (n2 %% 2 != 0) n2 <- n2 + 1
                    # print("Lowering") # Eliminate later
                    break
                    
                }
            } else if (fixed == "n2") {
                # Eta is close enough to the desired eta
                if (current_eta - eta < 0.1 && n1 - low == 1) {
                    final_SSD[[b]] <- SSD_object
                    singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                    b <- b + 1
                    low <- min_sample
                    previous_eta <- 0
                    previous_high <- 0
                    high <- max_sample
                    next
                    
                } else if (current_eta == previous_eta && n1 - low == 1) {
                    # If there is no change in eta and the lower bound is close to the middle point
                    final_SSD[[b]] <- SSD_object
                    singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                    b <- b + 1
                    low <- min_sample
                    previous_eta <- 0
                    previous_high <- 0
                    high <- max_sample
                    next
                    
                } else if (current_eta == previous_eta && low + n1 == high * 2) {
                    # Reached the minimum number that meets the Bayesian power condition
                    final_SSD[[b]] <- SSD_object
                    b <- b + 1
                    singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                    low <- min_sample
                    previous_eta <- 0
                    previous_high <- 0
                    high <- max_sample
                    next
                    
                } else {
                    # Decreasing the cluster size to find the ultimate sample size
                    low <- low                         #lower bound
                    high <- n1                         #higher bound
                    n1 <- round2((low + high) / 2)      #point in the middle
                    # print("Lowering") # Eliminate later
                    break
                }
            }
        } # Finish condition met
        
    } else if (hypotheses_set == 2) {
        if (condition_met == FALSE) {
            print(c("Using cluster size:", n1, "/", "and number of clusters:", n2, "/",
                    "prop_BF1c: ", current_eta,  "/",
                    "low:", low, " / high:", high,  "/",
                    " ultimate sample: ", ultimate_sample_size))
            if (fixed == "n1") {
                if ((n2 == max_sample) | (n2 > max_sample))    { # If the sample size reaches the maximum
                    final_SSD <- list("n1" = n1,
                                      "n2" = n2,
                                      "Proportion.BF1c" = current_eta,
                                      "data_H1" = results_H1,
                                      "singularity" = data_H1$singularity)
                    ultimate_sample_size <- TRUE
                } else {
                    message("Increasing sample")
                    # Increase the number of clusters since eta is too small
                    low <- n2                         #lower bound
                    high <- high                      #higher bound
                    n2 <- round2((low + high) / 2)     #point in the middle
                    if (n2 %% 2 != 0) n2 <- n2 + 1 # To ensure number of clusters is even
                    
                    # Adjust higher bound when there is a ceiling effect
                    if (low + n2 == high * 2) {
                        low <- n2                         #lower bound
                        if (previous_high > 0) {
                            high <- previous_high
                        } else {
                            high <- max_sample                   #higher bound
                        }
                        n2 <- round2((low + high) / 2)     #point in the middle
                    }
                }
            } else if (fixed == "n2") {
                if ((n1 == max_sample) | (n1 > max_sample))    {# If the sample size reaches the maximum
                    final_SSD <- list("n1" = n1,
                                      "n2" = n2,
                                      "Proportion.BF1c" = current_eta,
                                      "data_H1" = results_H1,
                                      "singularity" = data_H1$singularity)
                    ultimate_sample_size <- TRUE
                } else {
                    # Increase the cluster sizes since eta is too small
                    message("Increasing sample")
                    low <- n1                        #lower bound
                    high <- high                     #higher bound
                    n1 <- round2((low + high) / 2)    #point in the middle
                    
                    # Adjust higher bound when there is a ceiling effect
                    if ((low + n1 == high * 2) | (current_eta == previous_eta)) {
                        low <- n1                        #lower bound
                        # Set the higher bound based on the previous high or the maximum
                        if (previous_high > 0 ) {
                            high <- previous_high
                        } else {
                            high <- max_sample
                        }
                        n1 <- round2((low + high) / 2)    #point in the middle
                    }
                }
            }
            
        } else if (condition_met == TRUE) {
            print(c("Using cluster size:", n1,
                    "and number of clusters:", n2,
                    "prop_BF1c: ", current_eta,
                    "low: ", low, "high: ", high))
            previous_high <- high
            SSD_object <- list("n1" = n1,
                               "n2" = n2,
                               "Proportion.BF1c" = current_eta,
                               "data_H1" = results_H1,
                               "singularity" = data_H1$singularity)
            # print("Lowerign sample")
            # print(c("previous:", previous_eta))
            previous_eta <- current_eta
            
            if (fixed == "n1") {
                # Eta is close enough to the desired eta
                if (current_eta - eta < 0.1 && n2 - low == 2) {
                    final_SSD <- SSD_object
                    ultimate_sample_size <- TRUE
                    
                } else if (previous_eta == current_eta && n2 - low == 2) {
                    # If there is no change in eta and the lower bound is close to the middle point
                    final_SSD <- SSD_object
                    ultimate_sample_size <- TRUE
                    
                } else {
                    # Decreasing to find the ultimate number of clusters
                    message("Decreasing sample size")
                    low <- low                         #lower bound
                    high <- n2                         #higher bound
                    n2 <- round2((low + high) / 2)      #point in the middle
                    if (n2 %% 2 != 0) n2 <- n2 + 1

                }
            } else if (fixed == "n2") {
                # Eta is close enough to the desired eta
                if (current_eta - eta < 0.1 && n1 - low == 1) {
                    final_SSD <- SSD_object
                    ultimate_sample_size <- TRUE
                    
                } else if (current_eta == previous_eta && n1 - low == 1) {
                    # If there is no change in eta and the lower bound is close to the middle point
                    final_SSD <- SSD_object
                    ultimate_sample_size <- TRUE
                    
                } else if (current_eta == previous_eta && low + n1 == high * 2) {
                    # Reached the minimum number that meets the Bayesian power condition
                    final_SSD <- SSD_object
                    ultimate_sample_size <- TRUE
                    
                } else {
                    # Decreasing the cluster size to find the ultimate sample size
                    message("Decreasing sample size")
                    low <- low                         #lower bound
                    high <- n1                         #higher bound
                    n1 <- round2((low + high) / 2)      #point in the middle
                }
            }
        } # Finish condition met
    }

    return(list("final_SSD" = final_SSD,
                "n1" = n1, "n2" = n2,
                "ultimate_sample_size" = ultimate_sample_size,
                "previous_high" = previous_high, "low" = low, "high" = high))
}
