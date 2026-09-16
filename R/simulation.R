

# Project set up
dir.create("R") # Folder for results
dir.create("data")
dir.create("results")


########################### SIMULATION ########################################
library(tidyr)
library(nanoparquet) # Save design matrix as parquet file
library(dplyr)    # Data frame
library(ggplot2)
library(parallel)
library(tidyverse)   # To check simulation
library(tools)

source("R/ssd_crt_xo.R")
source("R/helper_simulation.R")

ndatasets <- 5000
max_sample <- 400
batch_size <- 250

#============= Figure 1: n1, n2, treatment effects =============
treat_eff <- c(0.2, 0.4, 0.6, 0.8)
wp_icc <- 0.05
cac <- 0.5
BF_thresh <- 10 
eta <- 0.8
periods <- c(2, 3, 4)

## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_f1.1 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_f1.2 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

# Design matrix
sim_design_f1 <- rbind(sim_design_f1.1, sim_design_f1.2)
sim_design_f1 <- mutate(sim_design_f1, seed = as.integer(sample(2^32 / 2, n())))
nrow_design <- nrow(sim_design_f1)

## Figure 1 folder
path <- "~/"
results_folder <- "data/figure1_data_3"
if (!dir.exists(results_folder)) {
    dir.create(results_folder)
}
write_parquet(sim_design_f1, paste0(results_folder, "/design_matrix_f1"))
# Read design matrix
sim_design_f1 <- read_parquet("data/figure1_data/design_matrix_f1")

# Run simulation
run_sim_wrapper <- function(Row) {
    run_sim(
        row = Row,
        design_matrix = sim_design_f1,
        ndatasets = ndatasets,
        Max = max_sample,
        batch_size = batch_size,
        results_folder = results_folder
    )
}

clusters <- makeForkCluster(detectCores()*0.5)
output <- parallel::parLapply(cl = clusters,
                              X = 1:nrow_design,
                              fun = run_sim_wrapper)

stopCluster(clusters)

# =================== Figure 2: n1, n2, wp_icc, cac =========================

treat_eff <- 0.2
BF_thresh <- 10
eta <- 0.8
periods <- 3
wp_icc <- c(0.01, 0.05, 0.1)
cac <- c(0.5, 0.7, 0.9)

## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_f2.1 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_f2.2 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

# Design matrix
sim_design_f2 <- rbind(sim_design_f2.1, sim_design_f2.2)
sim_design_f2 <- mutate(sim_design_f2, seed = as.integer(sample(2^32 / 2, n())))
nrow_design <- nrow(sim_design_f2)

# Figure 2 folder
path <- "~/"
results_folder <- "data/figure2_data_3"
if (!dir.exists(results_folder)) {
    dir.create(results_folder)
}
write_parquet(sim_design_f1, paste0(results_folder, "/design_matrix_f2"))
# Read design matrix
sim_design_f2 <- read_parquet("data/figure2_data/design_matrix_f2")

# Run simulation
run_sim_wrapper <- function(Row) {
    run_sim(
        row = Row,
        design_matrix = sim_design_f2,
        ndatasets = ndatasets,
        Max = max_sample,
        batch_size = batch_size,
        results_folder = results_folder
    )
}

clusters <- makeForkCluster(detectCores()*0.5)
output <- parallel::parLapply(cl = clusters,
                              X = 1:nrow_design,
                              fun = run_sim_wrapper)

stopCluster(clusters)

# ================== Figure 3: n1, n2, Bayes factor, eta =========================

treat_eff <- 0.2
wp_icc <- 0.05
cac <- 0.5
periods <- 3
BF_thresh <- c(3, 5, 10, 20)
eta <- c(0.7, 0.8, 0.9, 0.95)

## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_f3.1 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_f3.2 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

# Design matrix
sim_design_f3 <- rbind(sim_design_f3.1, sim_design_f3.2)
sim_design_f3 <- mutate(sim_design_f3, seed = as.integer(sample(2^32 / 2, n())))
nrow_design <- nrow(sim_design_f3)

# Figure 3 folder
path <- "~/"
results_folder <- "data/figure3_data_3"
if (!dir.exists(results_folder)) {
    dir.create(results_folder)
}
write_parquet(sim_design_f1, paste0(results_folder, "/design_matrix_f3"))
# Read design matrix
sim_design_f3 <- read_parquet("data/figure3_data/design_matrix_f3")

# Run simulation
run_sim_wrapper <- function(Row) {
    run_sim(
        row = Row,
        design_matrix = sim_design_f3,
        ndatasets = ndatasets,
        Max = max_sample,
        batch_size = batch_size,
        results_folder = results_folder
    )
}

clusters <- makeForkCluster(detectCores()*0.5)
output <- parallel::parLapply(cl = clusters,
                              X = 1:nrow_design,
                              fun = run_sim_wrapper)

stopCluster(clusters)

############################ FULL SIMULATION ###################################

library(tidyr)
library(nanoparquet) # Save design matrix as parquet file
library(dplyr)    # Data frame
library(ggplot2)
library(parallel)
library(tidyverse)
library(tools)   # To check simulation

source("ssd_crt_xo.R")
source("helper_simulation.R")

ndatasets <- 5000
max_sample <- 400
batch_size <- 625

treat_eff <- c(0.2, 0.4, 0.6, 0.8)
wp_icc <- c(0.01, 0.05, 0.1)
cac <- c(0.5, 0.7, 0.9)
BF_thresh <- c(3, 5, 10, 20)
eta <- c(0.7, 0.8, 0.9, 0.95)
periods <- 2

## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_set2.1 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_set2.2 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

# Design matrix
sim_design_set2 <- rbind(sim_design_set2.1, sim_design_set2.2)
sim_design_set2 <- mutate(sim_design_set2, seed = as.integer(sample(2^32 / 2, n())))
nrow_design <- nrow(sim_design_set2)

## Hypothesis set 2 folder
path <- "~/"
results_folder <- "set2_data"
if (!dir.exists(results_folder)) {
    dir.create(results_folder)
}
write_parquet(sim_design_set2,
              paste0(results_folder, "/design_matrix_set2"))
# Read design matrix
sim_design_set2 <- read_parquet("data/set2_data/design_matrix_set2")

# Run simulation
run_sim_wrapper <- function(Row) {
    run_sim(
        row = Row,
        design_matrix = sim_design_set2,
        ndatasets = ndatasets,
        Max = max_sample,
        batch_size = batch_size,
        results_folder = results_folder
    )
}

clusters <- makeForkCluster(detectCores() * 0.5)
output <- parallel::parLapply(cl = clusters, X = missing_until_now, fun = run_sim_wrapper)

a <- missing_rows(
    results_folder,
    name_pattern = "ResultsN1Row",
    check_numbers = 1:nrow_design,
    underscore = F
)
b <- missing_rows(
    results_folder,
    name_pattern = "ResultsN2Row",
    check_numbers = 1:nrow_design,
    underscore = F
)
missing_until_now <- intersect(a, b)

stopCluster(clusters)

####################### SIMULATION AFTER CHANGING CAC ###############################

library(tidyr)
library(nanoparquet) # Save design matrix as parquet file
library(dplyr)    # Data frame
library(ggplot2)
library(parallel)
library(tidyverse)
library(tools)   # To check simulation

source("ssd_crt_xo.R")
source("helper_simulation.R")

ndatasets <- 5000
max_sample <- 400
batch_size <- 625

treat_eff <- c(0.2, 0.4, 0.6, 0.8)
wp_icc <- c(0.01, 0.05, 0.1)
bp_icc <- c(0.005, 0.025, 0.04)
BF_thresh <- c(3, 5, 10, 20)
eta <- c(0.7, 0.8, 0.9, 0.95)
periods <- 2

## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_set2.1 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_set2.2 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

# Design matrix
sim_design_set2 <- rbind(sim_design_set2.1, sim_design_set2.2)
sim_design_set2 <- mutate(sim_design_set2, seed = as.integer(sample(2^32 / 2, n())))
nrow_design <- nrow(sim_design_set2)

## Hypothesis set 2 folder
path <- "~/"
results_folder <- "set2_data"
if (!dir.exists(results_folder)) {
    dir.create(results_folder)
}
write_parquet(sim_design_set2,
              paste0(results_folder, "/design_matrix_set2"))

# Run simulation
run_sim_wrapper <- function(Row) {
    run_sim(
        row = Row,
        design_matrix = sim_design_set2,
        ndatasets = ndatasets,
        Max = max_sample,
        batch_size = batch_size,
        results_folder = results_folder
    )
}

clusters <- makeForkCluster(detectCores() * 0.5)
output <- parallel::parLapply(cl = clusters, X = missing_until_now, fun = run_sim_wrapper)

a <- missing_rows(
    results_folder,
    name_pattern = "ResultsN1Row",
    check_numbers = 1:nrow_design,
    underscore = F
)
b <- missing_rows(
    results_folder,
    name_pattern = "ResultsN2Row",
    check_numbers = 1:nrow_design,
    underscore = F
)
missing_until_now <- intersect(a, b)

stopCluster(clusters)

####################### COLLECT RESULTS ################################
# Collect results: Figure 1
## Results
figure1_data <- collect_results(
    sim_design_f1,
    results_folder = "data/figure1_data",
    pair = 2,
    results_name = "Results",
    save = T,
    file_name = "figure1_data"
)
## Times
figure1_times <- collect_times(
    design_matrix = sim_design_f1,
    pair = 2,
    times_name = "time",
    results_folder = "data/figure1_data_2",
    file_name = "figure1_times"
)

# Collect results: Figure 2
## Results
figure2_data <- collect_results(
    sim_design_f2,
    results_folder = "data/figure2_data_2",
    pair = 2,
    results_name = "Results",
    save = T,
    file_name = "figure2_data"
)
## Times
figure2_times <- collect_times(
    design_matrix = sim_design_f2,
    pair = 2,
    times_name = "time",
    results_folder = "data/figure2_data",
    file_name = "figure2_times"
)

# Collect results: Figure 3
## Results
figure3_data <- collect_results(
    sim_design_f3,
    results_folder = "data/figure3_data_2",
    pair = 2,
    results_name = "Results",
    save = T,
    file_name = "figure3_data2"
)
## Times
figure3_times <- collect_times(
    design_matrix = sim_design_f3,
    pair = 2,
    times_name = "time",
    results_folder = "data/figure3_data",
    file_name = "figure3_times"
)

# Full simulation
## Results
full_simulation_data <- collect_results(
    design_matrix = sim_design_set2,
    results_folder = "data/set2_data",
    pair = 2,
    results_name = "Results",
    file_name = "full_sim_data",
    save = T
)

## Times
full_simulation_times <- collect_times(
    design_matrix = sim_design_set2,
    pair = 2,
    times_name = "time",
    results_folder = "data/set2_data",
    file_name = "full_sim_times"
)

######################### PLOTS #############################
## Figure 1--------
library(ggplot2)
# Plot median BF
ggplot(figure1_data,
       aes(
           y = log(median.BF1c),
           x = n1.final,
           color = as.factor(n2.final),
           shape = as.factor(n2.final)
       )) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(cols = vars(eff_size)) +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "Number of \nclusters",
         shape = "Number of \nclusters") +
    xlab("Cluster sizes") + ylab("Bayes Factor") + theme(legend.position = "bottom")
# Plot sample sizes
ggplot(figure1_data,
       aes(
           y = n1.final,
           x = n2.final,
           color = as.factor(eff_size),
           shape = as.factor(eff_size)
       )) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "Effect sizes",
         shape = "Effect sizes") +
    xlab("Cluster sizes") + ylab("Bayes Factor") + theme(legend.position = "bottom")

ggplot(figure1_data,
       aes(
           y = n2.final,
           x = n1.final,
           color = as.factor(eff_size),
           shape = as.factor(eff_size)
       )) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Number of clusters",
         color = "Effect sizes",
         shape = "Effect sizes") +
    xlab("Cluster sizes") + ylab("Bayes Factor") + theme(legend.position = "bottom")

ggplot(figure1_data,
       aes(
           y = n2.final,
           x = n1.final,
           color = as.factor(eff_size),
           shape = as.factor(eff_size)
       )) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Number of clusters",
         color = "Effect sizes",
         shape = "Effect sizes") +
    xlab("Cluster sizes") + ylab("Bayes Factor") + theme(legend.position = "bottom")

ggplot(figure1_data[figure1_data$fixed == "n2",],
       aes(
           y = eta.BF12,
           x = n2.final,
           color = as.factor(eff_size),
           shape = as.factor(eff_size)
       )) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "Effect size",
         shape = "Effect size") +
    xlab("Number of clusters") + ylab("Probability") + theme(legend.position = "bottom")

ggplot(figure1_data[figure1_data$fixed == "n1",],
       aes(
           y = eta.BF12,
           x = n1.final,
           color = as.factor(eff_size),
           shape = as.factor(eff_size)
       )) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "Effect size",
         shape = "Effect size") +
    xlab("Cluster size") + ylab("Probability") + theme(legend.position = "bottom")

## Figure 2-------------
ggplot(figure2_data[figure2_data$fixed == "n1",],
       aes(
           y = log(median.BF1c),
           x = n1.final,
           color = as.factor(cac),
           shape = as.factor(cac)
       )) +
    geom_point() + geom_line() +
    facet_grid(cols = vars(wp_icc)) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "CAC",
         shape = "CAC") +
    xlab("Cluster sizes") + ylab("Bayes Factor") + theme(legend.position = "bottom") 


ggplot(figure2_data[figure2_data$fixed == "n2",],
       aes(
           y = log(median.BF1c),
           x = n2.final,
           color = as.factor(cac),
           shape = as.factor(cac)
       )) +
    geom_point() + geom_line() +
    facet_grid(cols = vars(wp_icc)) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "CAC",
         shape = "CAC") +
    xlab("Cluster sizes") + ylab("Bayes Factor") + theme(legend.position = "bottom")

ggplot(figure2_data[figure2_data$fixed == "n2",],
       aes(
           y = eta.BF12,
           x = n2.final,
           color = as.factor(cac),
           shape = as.factor(cac)
       )) +
    geom_point() + geom_line() +
    facet_grid(cols = vars(wp_icc)) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "CAC",
         shape = "CAC") +
    xlab("Number of clusters") + ylab("Probability") + theme(legend.position = "bottom")

ggplot(figure2_data[figure2_data$fixed == "n1",],
       aes(
           y = eta.BF12,
           x = n1.final,
           color = as.factor(cac),
           shape = as.factor(cac)
       )) +
    geom_point() + geom_line() +
    facet_grid(cols = vars(wp_icc)) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "CAC",
         shape = "CAC") +
    xlab("Cluster size") + ylab("Probability") + theme(legend.position = "bottom")



## Figure 3---------------
ggplot(figure3_data[figure3_data$fixed == "n1",],
       aes(
           y = log(median.BF1c),
           x = n1.final,
           color = as.factor(BF_thresh),
           shape = as.factor(BF_thresh)
       )) +
    geom_point() + geom_line() +
    facet_grid(cols = vars(eta)) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "BF threshold",
         shape = "BF threshold") +
    xlab("Cluster size") + ylab("Bayes Factor") + theme(legend.position = "bottom")


ggplot(figure3_data[figure3_data$fixed == "n2",],
       aes(
           y = log(median.BF1c),
           x = n2.final,
           color = as.factor(BF_thresh),
           shape = as.factor(BF_thresh)
       )) +
    geom_point() + geom_line() +
    facet_grid(cols = vars(eta)) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "BF threshold",
         shape = "BF threshold") +
    xlab("Number of clusters") + ylab("Bayes Factor") + theme(legend.position = "bottom")

ggplot(figure3_data[figure3_data$fixed == "n2",],
       aes(
           y = eta.BF12,
           x = n2.final,
           color = as.factor(BF_thresh),
           shape = as.factor(BF_thresh)
       )) +
    geom_point() + geom_line() +
    facet_grid(cols = vars(eta)) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "BF threshold",
         shape = "BF threshold") +
    xlab("Number of clusters") + ylab("Probability") + theme(legend.position = "bottom")

ggplot(figure3_data[figure3_data$fixed == "n1",],
       aes(
           y = eta.BF12,
           x = n1.final,
           color = as.factor(BF_thresh),
           shape = as.factor(BF_thresh)
       )) +
    geom_point() + geom_line() +
    facet_grid(cols = vars(eta)) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs Hc",
         color = "BF threshold",
         shape = "BF threshold") +
    xlab("Cluster size") + ylab("Probability") + theme(legend.position = "bottom")


# Explore full simulation

library(ggplot2)

# Focus on a specific power target (e.g., eta = 0.8)
plot_data <- subset(full_sim_data_set2, fixed == "n2" & eta == 0.8)

ggplot(plot_data, aes(
    x = n2.final, 
    y = eta.BF12, 
    color = factor(wp_icc), 
    linetype = factor(cac),
    group = interaction(wp_icc, cac)
)) +
    geom_hline(yintercept = 0.80, linetype = "dotted", color = "grey40") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_grid(
        rows = vars(BF_thresh), 
        cols = vars(eff_size), 
        labeller = label_both,
        scales = "free_x"
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(
        title = "Probability of Reaching BF Threshold (Target Power η = 0.8)",
        x = "Number of Clusters (n2)",
        y = "P(BF > Threshold)",
        color = "Within-Period ICC",
        linetype = "CAC"
    ) +
    theme_bw(base_size = 11) +
    theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95")
    )


library(ggplot2)

plot_simulation_grid <- function(data, 
                                 fixed_param = "n2", 
                                 target_eta = 0.80, 
                                 palette = "Dark2") {
    
    # Select the correct x-axis variable based on fixed parameter design
    x_var <- if (fixed_param == "n2") "n2.final" else "n1.final"
    x_label <- if (fixed_param == "n2") "Number of Clusters (n2)" else "Cluster Size (n1)"
    
    plot_df <- subset(data, fixed == fixed_param & eta == target_eta)
    
    ggplot(plot_df, aes(
        x = .data[[x_var]], 
        y = eta.BF12, 
        color = factor(wp_icc), 
        linetype = factor(cac),
        shape = factor(cac),
        group = interaction(wp_icc, cac)
    )) +
        geom_hline(yintercept = target_eta, linetype = "dotted", color = "grey40", linewidth = 0.7) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2) +
        facet_grid(
            rows = vars(BF_thresh), 
            cols = vars(eff_size), 
            labeller = label_both,
            scales = "free_x"
        ) +
        scale_color_brewer(palette = palette) +
        scale_y_continuous(limits = c(0.7, 1), breaks = seq(0, 1, 0.2)) +
        labs(
            title = paste0("Probability of Exceeding BF Threshold (Fixed: ", fixed_param, ", Target η = ", target_eta, ")"),
            x = x_label,
            y = "P(BF > Threshold)",
            color = "Within-Period ICC",
            linetype = "CAC",
            shape = "CAC"
        ) +
        theme_bw(base_size = 11) +
        theme(
            legend.position = "bottom",
            legend.box = "horizontal",
            strip.background = element_rect(fill = "grey92"),
            panel.grid.minor = element_blank()
        )
}

# Example usage:

plot_simulation_grid(full_sim_data_set2, fixed_param = "n2", target_eta = 0.70)
plot_simulation_grid(full_sim_data_set2, fixed_param = "n1", target_eta = 0.70)
plot_simulation_grid(full_sim_data_set2, fixed_param = "n2", target_eta = 0.80)
plot_simulation_grid(full_sim_data_set2, fixed_param = "n1", target_eta = 0.80)
plot_simulation_grid(full_sim_data_set2, fixed_param = "n1", target_eta = 0.90)
plot_simulation_grid(full_sim_data_set2, fixed_param = "n2", target_eta = 0.90)
plot_simulation_grid(full_sim_data_set2, fixed_param = "n2", target_eta = 0.95)
plot_simulation_grid(full_sim_data_set2, fixed_param = "n1", target_eta = 0.95)



#############################################################################
####### EXTRA SIMULATION WITH MORE PERIODS AND LARGER THRESHOLDS ############


treat_eff <- c(0.2, 0.4, 0.6, 0.8)
wp_icc <- c(0.01, 0.05, 0.1)
bp_icc <- c(0.005, 0.025, 0.04)
BF_thresh <- c(3, 5, 10, 20, 50, 100)
eta <- c(0.7, 0.8, 0.9, 0.95)
periods <- c(2, 3, 4)


## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_set2.1 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_set2.2 <- expand.grid(
    "eff_size" = treat_eff,
    "fixed" = fixed,
    "n2" = n2,
    "n1" = n1,
    "wp_icc" = wp_icc,
    "cac" = cac,
    "BF_thresh" = BF_thresh,
    "eta" = eta,
    "periods" = periods
)

# Design matrix
sim_design_set2 <- rbind(sim_design_set2.1, sim_design_set2.2)
sim_design_set2 <- mutate(sim_design_set2, seed = as.integer(sample(2^32 / 2, n())))
nrow_design <- nrow(sim_design_set2)

# Filter rows to include only new conditions
big_simulation <- read_parquet("data/set2_data/design_matrix_set2")
columns_to_compare <- names(big_simulation)[1:9]
key1 <- do.call(paste, c(big_simulation[columns_to_compare], sep = "___"))
key2 <- do.call(paste, c(sim_design_set2[columns_to_compare], sep = "___"))
common_index <- length(which(key2 %in% key1))

additional_combinations <- anti_join(sim_design_set2[, 1:9], big_simulation[, 1:9])

## Hypothesis set 2 
path <- "~/"
results_folder <- "set2_data_extra"
if (!dir.exists(results_folder)) {
    dir.create(results_folder)
}
write_parquet(sim_design_set2,
              paste0(results_folder, "/design_matrix_set2"))

# Run simulation
run_sim_wrapper <- function(Row) {
    run_sim(
        row = Row,
        design_matrix = sim_design_set2,
        ndatasets = ndatasets,
        Max = max_sample,
        batch_size = batch_size,
        results_folder = results_folder
    )
}

clusters <- makeForkCluster(detectCores() * 0.5)
output <- parallel::parLapply(cl = clusters, X = missing_until_now, fun = run_sim_wrapper)
