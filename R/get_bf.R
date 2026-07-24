########################## GET THE BAYES  FACTOR ##############################
# Approximated Adjusted Fractional Bayes Factor



get_BF <- function(n_eff, treatment, var_treatm, hypotheses, b){
    # Unequality constrained hypotheses---------------------------------
    if (hypotheses == 2) {
        # Complexities
        comp1 <- .5
        comp2 <- .5
        
        # Fit
        fit2 <- pnorm(0, mean = treatment, sd = sqrt(var_treatm))
        fit1 <- 1 - fit2
        
        # Calculation BFs
        AAFBF1u <- (fit1 / comp1)
        AAFBF2u <- (fit2 / comp2)
        AAFBF12 <- AAFBF1u / AAFBF2u
        AAFBF21 <- 1 / AAFBF12
        
        #Calculation of PMPs
        pmp1 <- AAFBF1u / (AAFBF1u + AAFBF2u)
        pmp2 <- 1 - pmp1
        output <- list(bf.12 = AAFBF12, bf.21 = AAFBF21, pmp1 = pmp1, pmp2 = pmp2)
        
    } else if (hypotheses == 1) {
        # Equality constrained hypothesis-------------------------------------------
        b_calc <- b * 1 / n_eff                   # Calculate b
        
        # Complexities
        comp0 <- dnorm(0, mean = 0, sd = sqrt(var_treatm / b_calc))     # overlap of parameter under H0 and unconstrained prior -> density of the prior under Hu at the focal point 0
        comp1 <- 1 - pnorm(0, mean = 0, sd = sqrt(var_treatm / b_calc))
        
        # Fit
        fit0 <- dnorm(0, mean = treatment, sd = sqrt(var_treatm)) # overlap of parameter under H0 and posterior -> density of the posterior at focal point 0
        fit1 <- 1 - pnorm(0, mean = treatment, sd = sqrt(var_treatm)) # the fit is equal to 1 - the fit of the complement
        
        # Calculation of BFs
        AAFBF0u <- fit0 / comp0                    # AAFBF of H0 vs Hu
        AAFBF1u <- fit1 / comp1                    # AAFBF of H1 vs Hu
        AAFBF01 <- AAFBF0u / AAFBF1u
        AAFBF10 <- 1 / AAFBF01
        
        # Calculation of PMPs
        pmp0 <- AAFBF0u / (AAFBF0u + AAFBF1u)
        pmp1 <- 1 - pmp0
        output <- list(bf.10 = AAFBF10, bf.01 = AAFBF01, pmp0 = pmp0, pmp1 = pmp1)
    }
    
    # Return results
    return(output)
    
}