############################ FUNCTIONS FOR PLOTS ##############################

# Density plot of Bayes factors
density_plot <- function(eta, data_plot, b, BF_threshold){
    
    library(ggplot2)
    library(ggpubr)
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#661100")
    
    # Draw plots
    # When working with equality constrains
    if (input$hypotheses == "Equality") {
        plot_left <- ggplot(data_plot, aes(x = log(BF_01), color = as.factor(b), fill = as.factor(b))) + 
            geom_histogram(binwidth = 1, aes(y = after_stat(density)), alpha = 0.5, position = "identity") +
            geom_density(alpha = .2) +
            geom_vline(aes(xintercept = log(BF_threshold), colour = cbbPalette[9]), linetype = "dashed") +
            scale_fill_manual(values = c("1" = cbbPalette[2], "2" = cbbPalette[3], "3" = cbbPalette[8]), name = "Fraction b") +
            scale_color_manual(values = c("1" = cbbPalette[2], "2" = cbbPalette[3], "3" = cbbPalette[8]), name = "Fraction b") +
            ylab("Density") + xlab(bquote("log Bayes factor"["01"])) +
            theme(legend.position = "bottom", axis.title = element_text(size = 16),
                  axis.text = element_text(size = 15),
                  plot.caption = element_text(hjust = 0, size = 12),
                  plot.title = element_text(size = 14, face = "bold")) + 
            labs(title = "Null hypothesis is true")
        
        plot_right <- ggplot(data_plot, aes(x = log(BF_10), color = as.factor(b), fill = as.factor(b))) + 
            geom_histogram(binwidth = 1, aes(y = after_stat(density)), alpha = 0.5, position = "identity") +
            geom_density(alpha = .2) +
            geom_vline(aes(xintercept = log(BF_threshold), colour = cbbPalette[9]), linetype = "dashed") +
            scale_fill_manual(values = c("1" = cbbPalette[2], "2" = cbbPalette[3], "3" = cbbPalette[8]), name = "Fraction b") +
            scale_color_manual(values = c("1" = cbbPalette[2], "2" = cbbPalette[3], "3" = cbbPalette[8]), name = "Fraction b") +
            ylab("Density") + xlab(bquote("log Bayes factor"["10"])) + labs(title = "Alternative hypothesis is true") +
            theme(legend.position = "bottom", axis.title = element_text(size = 16),
                  axis.text = element_text(size = 15),
                  plot.title = element_text(size = 14, face = "bold"))
        
    } # When working with inequality constraints
    else if (input$hypotheses == "Informative") {
        plot_BF <- ggplot(data_plot, aes(x = log(BF_12))) + 
            geom_histogram(binwidth = 1, colour = "#69b3a2", fill = "#69b3a2", 
                           aes(y = after_stat(density)), alpha = 0.8, position = "identity") +
            geom_density(alpha = .2, colour = "#69b3a2") + 
            geom_vline(aes(xintercept = log(BF_threshold), colour = cbbPalette[9]), linetype = "dashed") +
            ylab("Density") + xlab(bquote("log Bayes factor"["12"])) + 
            labs(caption = paste0("Eta (\u03B7) = ", eta)) +
            theme(legend.position = "none", axis.title = element_text(size = 16),
                  axis.text = element_text(size = 15),
                  plot.caption = element_text(hjust = 0, size = 12))
    }
    # Arrange plots
    if (input$hypotheses == "Equality") {
        grid_plots <- ggarrange(plot_left, plot_right, ncol = 2, common.legend = TRUE,
                                legend = "bottom")
        annotate_figure(grid_plots, bottom = text_grob(as.character(paste0("Eta (\u03B7) = ", eta))))
    } else if (input$hypotheses == "Informative") {
        plot_BF
    }
}


# Plot increasing sample
increasing_sample_plot <- function(eff_size, cac, wp_icc, bp_icc, periods, n2, n1, treatment_n, 
                                   pattern, seed, BF_thresh, batch_size, max_sample, 
                                   ndatasets, fixed, gc_every = 10, set, b){
    
    # Libraries
    library(lme4)
    library(bain)
    library(BFpack)
    library(ggplot2)
    
    # Warnings
    if (is.numeric(c(eff_size, n1, n2, ndatasets, BF_thresh, max_sample, batch_size)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    if (eff_size < 0) stop("The effect size must be a positive value ")
    if (is.character(fixed) == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    
    # Functions
    source("data_generation_xo.R")
    source("helper_functions.R")
    source("get_bf.R")
    source("print_results.R")
    
    # Initial values
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#661100")
    n2_seq <- n2/2 #Number of clusters per sequence
    current_eta <- 0
    i <- 1
    
    results_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)
    if (fixed == "n1") {
        rows <- ((max_sample - n2)/2) + 1
        final_SSD <- matrix(NA, nrow = rows, ncol = 5)
    } else if (fixed == "n2")
        rows <- (max_sample - n1) + 1
        final_SSD <- matrix(NA, nrow = rows, ncol = 5)
    
    # Generate data ------------------------------------------------------------
    while (i < (rows + 1)) {
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
        
        
        
        
        # Compute Bayes factor --------------------------------------------------
        n_eff <- ((n1 * n2) / (1 + (n1 - 1) * data_H1$emp_wp_icc)) / 2 # Effecttive sample size
        output_bf <- Map(get_BF, n_eff, data_H1$treatment_eff,  data_H1$var_treat, 2)
        
        # Proportion
        
        if (set == 1) { #TODO
            # If H0 is true
            if (missing(bp_icc)) {
            data_H0 <- do.call(data_generation, list(eff_size = eff_size, 
                                                     cac = cac, wp_icc = wp_icc,
                                                     periods = periods, n2 = n2, 
                                                     n1 = n1,
                                                     treatment_n = treatment_n,
                                                     seed = seed,
                                                     ndatasets = ndatasets,
                                                     batch_size = batch_size,
                                                     gc_every = gc_every, 
                                                     b = b))
            } else {
                data_H0 <- do.call(data_generation, list(eff_size = eff_size,
                                                         bp_icc = bp_icc,
                                                         wp_icc = wp_icc,
                                                         periods = periods, n2 = n2, 
                                                         n1 = n1,
                                                         treatment_n = treatment_n,
                                                         seed = seed,
                                                         ndatasets = ndatasets,
                                                         batch_size = batch_size,
                                                         gc_every = gc_every,
                                                         b = b))
            }
            
        } else {
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
            current_eta <- length(which(results_H1[, "BF.1c"] > BF_thresh)) / ndatasets
            median_Bf12 <- median(results_H1[, 1])
            mean_pmp1 <- mean(results_H1[, 2])
        }

        final_SSD[i, ] <- c(n1, n2, current_eta, median_Bf12, mean_pmp1)
        
        i <- i + 1
        # Increasing sample
        if (fixed == "n1") {
            n2 <- n2 + 2
        } else if (fixed == "n2") {
            n1 <- n1 + 1
        }
    }
    
    final_SSD <- as.data.frame(final_SSD)
    colnames(final_SSD) <- c("n1", "n2", "eta", "median_Bf12", "mean_pmp1")
    # Plots
    if (set == 2) {
        if (fixed == "n1") {
            plot_output <- ggplot(final_SSD, aes(x = n2, y = current_eta)) +
                geom_line(colour = cbbPalette[6], linewidth = 1.2) +
                geom_point(colour = cbbPalette[3], size = 2) +
                labs(x = "Number of clusters", y = "Probability of exceeding threshold") +
                theme_minimal()
        } else {
            plot_output <- ggplot(final_SSD, aes(x = n1, y = current_eta)) +
                geom_line(colour = cbbPalette[6], linewidth = 1.2) +
                geom_point(colour = cbbPalette[3], size = 2) +
                labs(x = "Cluster size", y = "Probability of exceeding threshold") +
                theme_minimal()
        }
    } else if (set == 1) { #TODO
        if (fixed == "n1") {
            plot_output <- ggplot(final_SSD, aes(x = n2, y = current_eta)) +
                geom_line(colour = cbbPalette[6], linewidth = 1.2) +
                geom_point(colour = cbbPalette[3], size = 2) +
                labs(x = "Number of clusters", y = "Probability of exceeding threshold") +
                theme_minimal()
        } else {
            plot_output <- ggplot(final_SSD, aes(x = n1, y = current_eta)) +
                geom_line(colour = cbbPalette[6], linewidth = 1.2) +
                geom_point(colour = cbbPalette[3], size = 2) +
                labs(x = "Cluster size", y = "Probability of exceeding threshold") +
                theme_minimal()
        }
    }
    
    return(plot_output)
}



# # Tests
# increasing_sample_plot(eff_size = 0.4, cac = 0.2, wp_icc = 0.1, periods = 2, 
#                        n2 = 10, n1 = 10, treatment_n = 2, 
#                        pattern = "chess", seed = 408, BF_thresh = 3,
#                        batch_size = 100, max_sample = 30, 
#                        ndatasets = 100, fixed = "n2", gc_every = 10, set = 2)
