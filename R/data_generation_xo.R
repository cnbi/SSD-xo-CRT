############################### DATA GENERATION ################################

# From: https://github.com/klgrantham/bayesian-SW/blob/master/sim.R
# From: https://github.com/KMTanvir/RedundantCrossovers/tree/master

data_generation <- function(eff_size, cac, wp_icc, bp_icc, periods, n2, n1, treatment_n, 
                            seed, ndatasets, batch_size, gc_every = 10){
    
    treatment_eff <- eff_size
    
    # Create objects for results
    data_list <- vector(mode = "list")
    output_lmer <- vector(mode = "list")
    
    # Design matrix with sequences
    pattern <- "chess"
    sequences_matrix <- matrix(0, nrow = n2, ncol = periods)
    i <- 1
    j <- 0
    small_matrix <- matrix(0, 2, ncol = periods)
    
    if (pattern == "chess") {
        large_pattern <- rep(c(1, 0), periods)
        while (i < (periods + 1)) {
            small_matrix[i, ] <- large_pattern[i:(periods + j)]
            i <- i + 1
            j <- j + 1
        }
        sequences_matrix <- do.call(rbind, replicate(n2/2, small_matrix, simplify = FALSE))
    }
    
    # Period effects
    period_effects <- 0.1 * matrix(data = seq(1:periods), 
                                   nrow = n2, ncol = periods, byrow = TRUE)
    period_eff_vec <- as.data.frame(as.table(period_effects)) # In long format
    period_eff_vec <- period_eff_vec[order(period_eff_vec[, 1]), ][ , 3]
    period_eff_vec <- rep(period_eff_vec, n1)
    
    # Treatment effects vector
    treat_cond_vec <- as.data.frame(as.table(sequences_matrix))[, 3] # In long format
    treat_cond_vec <- rep(treat_cond_vec, n1)
    
    
    # Variances
    total_var <- 1
    sigma_e <- 1 - wp_icc
    if (missing(cac)) {
        sigma_c <- bp_icc
    } else {
        sigma_c <- wp_icc * cac
    }
    
    sigma_g <- wp_icc - sigma_c
    
    # Seeds
    if (missing(seed)) {
        seeds <- sample(2^32 / 2, ndatasets)
    } else {
        set.seed(seed)
        seeds <- sample(2^32 / 2, ndatasets)
    }
    
    # Create objects for results
    treatments   <- vector("numeric", ndatasets)
    var_fix_eff  <- vector("numeric", ndatasets)
    emp_var_g    <- vector("numeric", ndatasets)
    emp_var_c    <- vector("numeric", ndatasets)
    emp_var_e    <- vector("numeric", ndatasets)
    
    cluster <- rep(seq(1, n2), each = periods * n1)
    period <- rep(seq(1, periods), n1 * n2)
    id <- rep(seq(1, (n1*n2)), each = periods) # Variable to organise g
    
    # Create batches
    nbatches <- ceiling(ndatasets / batch_size)
    
    # Creating and fitting model by batches
    for (batch in seq(nbatches)) {
        #Indexes
        start_index <- (batch_size * (batch - 1)) + 1
        end_index <- min(batch * batch_size, ndatasets)
        #Seeds
        batch_seeds <- seeds[start_index:end_index]
        
        batch_data <- lapply(batch_seeds, function(s) {
            set.seed(s)
            e <- rnorm(n1 * periods * n2, mean = 0, sd = sqrt(sigma_e)) # Residual
            c <- rnorm(n2, mean = 0, sd = sqrt(sigma_c))                # Cluster-level random intercept
            c <- rep(c, each = n1 * periods)
            g <- rnorm(n2 * periods, mean = 0, sd = sqrt(sigma_g))      # Cluster-period random effect
            g <- rep(g, each = n1)
            
            
            order_g <- cbind(id, cluster, period)
            order_g <- order_g[order(order_g[, "cluster"], order_g[, "period"]), ]
            order_g <- cbind(order_g, g)
            order_g <- order_g[order(order_g[ , "id"]), ]
            g <- order_g[, "g"]
            
            # Response
            resp <- period_eff_vec + treatment_eff * treat_cond_vec + c + g + e           # In long format
            
            # data.frame(
            #     cluster = as.factor(cluster),
            #     period = as.factor(period),
            #     period_eff = period_eff_vec, 
            #     treat_cond = treat_cond_vec, 
            #     c = c,
            #     g = g,
            #     e = e, 
            #     resp = resp
            # )
            # 
            df <- data.frame(
                cluster = as.factor(cluster),
                period = as.factor(period),
                treat_cond = treat_cond_vec,
                resp = resp
            )
            
            fitted_model <- fit_lmer(df)
            
            # Extract estimations
            if (is.null(fitted_model)) {
                return(list(
                    treat_cond = NA_real_,
                    var_treat  = NA_real_,
                    var_g = NA_real_, var_c = NA_real_, var_e = NA_real_,
                    singular = NA
                ))
            }
            
            treatment_eff_ <- fixef(fitted_model)["treat_cond"]
            var_treatment <- vcov_(fitted_model, "treat_cond")
            var_rand_eff <- get_var(fitted_model)
            singular_marker <- tryCatch(lme4::isSingular(fitted_model), error = function(e) NA)
            
            output_est <- list(
                treat_cond = unname(treatment_eff_),
                var_treat  = unname(var_treatment),
                var_g = unname(var_rand_eff["var_g"]),
                var_c = unname(var_rand_eff["var_c"]),
                var_e = unname(var_rand_eff["var_e"]),
                singular = singular_marker
            )
            
            rm(fitted_model, df, e, c, g, order_g)
            return(output_est)
        })
        
        # Save results 
        i_results <- start_index:end_index
        for (k in seq(i_results)) {
            treatments[i_results[k]] <- batch_data[[k]]$treat_cond
            var_fix_eff[i_results[k]]  <- batch_data[[k]]$var_treat
            emp_var_g[i_results[k]]    <- batch_data[[k]]$var_g
            emp_var_c[i_results[k]]    <- batch_data[[k]]$var_c
            emp_var_e[i_results[k]]    <- batch_data[[k]]$var_e
        }
        
        rm(batch_data)
        if (batch %% gc_every == 0 || batch == nbatches) gc()
    }
    
    # Calculate ICC and total variance
    emp_total_var <- emp_var_c + emp_var_e + emp_var_g
    emp_wp_icc <- (emp_var_c + emp_var_g) / emp_total_var
    emp_cac <- emp_var_c/(emp_var_c + emp_var_g)
    
    
    # Return
    return(list("treatment_eff" = treatments,
                "var_treat" = var_fix_eff,
                "emp_total_var" = emp_total_var,
                "emp_wp_icc" = emp_wp_icc,
                "emp_cac" = emp_cac))
}


# # 1) empirical variance of simulated cluster-period effects (one per cluster-period)
# cpvals <- aggregate(g ~ cluster + period, data = data_, FUN = unique)
# cat(sprintf("Empirical var(g) across cluster-periods: %.6f (target sigma_g^2 = %.6f)\n",
#             var(cpvals$g), sigma_g))
# 
# # 2) empirical variance of simulated cluster intercepts (one per cluster)
# cvals <- aggregate(c ~ cluster, data = data_, FUN = unique)
# cat(sprintf("Empirical var(c) across clusters: %.6f (target sigma_c^2 = %.6f)\n",
#             var(cvals$c), sigma_c))

# Test -----------
# data_generation(eff_size = 0.4, cac = 0.8, wp_icc = 0.1, periods = 2, n2 = 30, 
#                 n1 = 20, treatment_n = 2, seed = 2207, ndatasets = 10, batch_size = 10)