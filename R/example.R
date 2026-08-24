############################################# EXAMPLE IN MANUSCRIPT #################################

source("ssd_crt_xo.R")
source("helper_simulation.R")

# Initial values
effect_size <- 0.4
eta <- 0.7
BF_thres <- 5
ndatasets <- 1000
wp_icc <- 0.2
bp_icc <- 0.08

# Find cluster size
start_time <- Sys.time()
example_n1 <- ssd_crt_xo(
    eff_size = effect_size,
    wp_icc = wp_icc,
    bp_icc = bp_icc,
    periods = 2,
    n2 = 30,
    n1 = 10,
    treatment_n = 2,
    pattern = "chess",
    BF_thresh = BF_thres,
    seed = 853,
    eta = eta,
    batch_size = 100,
    max_sample = 30,
    ndatasets = ndatasets,
    fixed = "n2",
    gc_every = 5
)
end_time <- Sys.time()
# Save results
example$data_H1


# Plot
library(ggplot2)
library(ggpubr) # arrange plots
library(ggsci) #
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#661100")
right <- ggplot(example_n1$data_H1, aes(x = log(BF.1c))) + 
    geom_histogram(binwidth = 1, colour = "#69b3a2", fill = "#69b3a2", 
                   aes(y = after_stat(density)), alpha = 0.8, position = "identity") +
    geom_density(alpha = 0.3, colour = "#69b3a2", fill = "#69b3a2") + 
    geom_vline(aes(xintercept = log(BF_thres), colour = cbbPalette[9]), linetype = "dashed",
               size = 0.85) +
    ylab("Density") + xlab(bquote("log Bayes factor"["1c"])) + 
    theme(legend.position = "none", axis.title = element_text(size = 16),
          axis.text = element_text(size = 3),
          plot.caption = element_text(hjust = 0, size = 12)) +
    annotate("text", x = c(log(BF_thres), 3.2), y = c(0.23, 0.05), parse = TRUE, 
             hjust = -0.1,
             label = c("BF[t][h][r][e][s] == 5", as.character(expression(~eta == 0.7))),
             colour = "grey30", size = 4) +
    ylim(c(0, 0.25)) + theme_minimal() + guides(color = "none")



# Find number of clusters
example_n2 <- ssd_crt_xo(
    eff_size = effect_size,
    wp_icc = wp_icc,
    bp_icc = bp_icc,
    periods = 2,
    n2 = 30,
    n1 = 4,
    treatment_n = 2,
    pattern = "chess",
    BF_thresh = BF_thres,
    seed = 853,
    eta = eta,
    batch_size = 100,
    max_sample = 30,
    ndatasets = ndatasets,
    fixed = "n1",
    gc_every = 5
)

# Save results
example_n2$data_H1

# Plot
left <- ggplot(example_n2$data_H1, aes(x = log(BF.1c))) + 
    geom_histogram(binwidth = 1, colour = "#69b3a2", fill = "#69b3a2",  
                   aes(y = after_stat(density)), alpha = 0.8, position = "identity") +
    geom_density(alpha = 0.4, colour = "#69b3a2", fill = "#69b3a2") + 
    geom_vline(aes(xintercept = log(BF_thres), colour = cbbPalette[9]), linetype = "dashed",
               size = 0.85) +
    ylab("Density") + xlab(bquote("log Bayes factor"["1c"])) + 
    theme(legend.position = "none", axis.title = element_text(size = 16),
          axis.text = element_text(size = 3),
          plot.caption = element_text(hjust = 0, size = 12)) + 
    annotate("text", x = c(log(BF_thres), 1.9), y = c(0.23, 0.05), hjust = -0.1,
             parse = TRUE, label = c("BF[t][h][r][e][s] == 5", as.character(expression(~eta == 0.7))),
             colour = "grey30", size = 4) + ylim(c(0, 0.25)) + theme_minimal() +
    guides(color = "none")

example_comb <- ggarrange(left, right, ncol = 2, legend = "none")
ggsave("sampling_dis_example.pdf", example_comb,
       device = cairo_pdf,
       width = 180,
       height = 130,
       units = "mm",
       dpi    = 300,
       path = "C:/Users/barra006/OneDrive - Universiteit Utrecht/Documents/GitHub/SSD-xo-CRT/output/plots")
# ---------------------------------------------
