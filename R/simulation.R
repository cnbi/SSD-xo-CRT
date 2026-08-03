

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

source("ssd_crt_xo.R")
source("helper_simulation.R")

ndatasets <- 5000
max_sample <- 400
batch_size <- 1000

#============= First figure: n1, n2, treatment effects =============
treat_eff <- c(0.2, 0.4, 0.6, 0.8)
wp_icc <- 0.05
cac <- 0.8
BF_thresh <- 5
eta <- 0.8
periods <- 2

## Find n1
fixed <- "n2"
n2 <- c(20, 40, 60)
n1 <- 10
sim_design_f1.1 <- expand.grid("eff_size" = treat_eff, "fixed" = fixed, "n2" = n2,
                               "n1" = n1, "wp_icc" = wp_icc, "cac" = cac,
                               "BF_thresh" = BF_thresh, "eta" = eta, "periods" = periods)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_f1.2 <- expand.grid("eff_size" = treat_eff, "fixed" = fixed, "n2" = n2,
                               "n1" = n1, "wp_icc" = wp_icc, "cac" = cac,
                               "BF_thresh" = BF_thresh, "eta" = eta, "periods" = periods)

# Design matrix
sim_design_f1 <- rbind(sim_design_f1.1, sim_design_f1.2)
sim_design_f1 <- mutate(sim_design_f1, seed = as.integer(sample(2 ^ 32 / 2, n())))
nrow_design <- nrow(sim_design_f1)

## Figure 1 folder
path <- "~/"
results_folder <- "figure1_data"
if (!dir.exists(results_folder)) {dir.create(results_folder)}
write_parquet(sim_design_f1, paste0(results_folder,"/design_matrix_f1"))

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

clusters <- makeForkCluster(nrow_design)
output <- parallel::parLapply(cl = clusters, 
                              X = 1:nrow_design, 
                              fun = run_sim_wrapper)

stopCluster(clusters)

# =================== Figure 2: n1, n2, wp_icc, cac =========================

treat_eff <- 0.4
BF_thresh <- 5
eta <- 0.8
periods <- 2
wp_icc <- c(0.01, 0.05, 0.1)
cac <- c(0.5, 0.8, 0.95)

## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_f2.1 <- expand.grid("eff_size" = treat_eff, "fixed" = fixed, "n2" = n2,
                               "n1" = n1, "wp_icc" = wp_icc, "cac" = cac,
                               "BF_thresh" = BF_thresh, "eta" = eta, "periods" = periods)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_f2.2 <- expand.grid("eff_size" = treat_eff, "fixed" = fixed, "n2" = n2,
                               "n1" = n1, "wp_icc" = wp_icc, "cac" = cac,
                               "BF_thresh" = BF_thresh, "eta" = eta, "periods" = periods)

# Design matrix
sim_design_f2 <- rbind(sim_design_f2.1, sim_design_f2.2)
sim_design_f2 <- mutate(sim_design_f2, seed = as.integer(sample(2 ^ 32 / 2, n())))
nrow_design <- nrow(sim_design_f2)

# Figure 2 folder
path <- "~/"
results_folder <- "figure2_data"
if (!dir.exists(results_folder)) {dir.create(results_folder)}
write_parquet(sim_design_f1, paste0(results_folder,"/design_matrix_f2"))

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

clusters <- makeForkCluster(nrow_design)
output <- parallel::parLapply(cl = clusters, 
                              X = 1:nrow_design, 
                              fun = run_sim_wrapper)

stopCluster(clusters)

# ================== Figure 3: n1, n2, Bayes factor, eta =========================

treat_eff <- 0.4
wp_icc <- 0.05
cac <- 0.8
periods <- 2
BF_thresh <- c(3, 5, 10, 20)
eta <- c(0.7, 0.8, 0.9, 0.95)

## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_f3.1 <- expand.grid("eff_size" = treat_eff, "fixed" = fixed, "n2" = n2,
                               "n1" = n1, "wp_icc" = wp_icc, "cac" = cac,
                               "BF_thresh" = BF_thresh, "eta" = eta, "periods" = periods)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_f3.2 <- expand.grid("eff_size" = treat_eff, "fixed" = fixed, "n2" = n2,
                               "n1" = n1, "wp_icc" = wp_icc, "cac" = cac,
                               "BF_thresh" = BF_thresh, "eta" = eta, "periods" = periods)

# Design matrix
sim_design_f3 <- rbind(sim_design_f3.1, sim_design_f3.2)
sim_design_f3 <- mutate(sim_design_f3, seed = as.integer(sample(2 ^ 32 / 2, n())))
nrow_design <- nrow(sim_design_f3)

# Figure 3 folder
path <- "~/"
results_folder <- "figure3_data"
if (!dir.exists(results_folder)) {dir.create(results_folder)}
write_parquet(sim_design_f1, paste0(results_folder,"/design_matrix_f3"))

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

clusters <- makeForkCluster(nrow_design)
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
cac <- c(0.5, 0.8, 0.95)
BF_thresh <- c(3, 5, 10, 20)
eta <- c(0.7, 0.8, 0.9, 0.95)
periods <- 2

## Find n1
fixed <- "n2"
n2 <- c(10, 20, 40, 60)
n1 <- 10
sim_design_set2.1 <- expand.grid("eff_size" = treat_eff, "fixed" = fixed, "n2" = n2,
                                 "n1" = n1, "wp_icc" = wp_icc, "cac" = cac,
                                 "BF_thresh" = BF_thresh, "eta" = eta, "periods" = periods)

## Find n2
fixed <- "n1"
n1 <- c(5, 10, 20, 30)
n2 <- 40
sim_design_set2.2 <- expand.grid("eff_size" = treat_eff, "fixed" = fixed, "n2" = n2,
                                 "n1" = n1, "wp_icc" = wp_icc, "cac" = cac,
                                 "BF_thresh" = BF_thresh, "eta" = eta, "periods" = periods)

# Design matrix
sim_design_set2 <- rbind(sim_design_set2.1, sim_design_set2.2)
sim_design_set2 <- mutate(sim_design_set2, seed = as.integer(sample(2 ^ 32 / 2, n())))
nrow_design <- nrow(sim_design_set2)

## Hypothesis set 2 folder
path <- "~/"
results_folder <- "set2_data"
if (!dir.exists(results_folder)) {dir.create(results_folder)}
write_parquet(sim_design_set2, paste0(results_folder,"/design_matrix_set2"))

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
output <- parallel::parLapply(cl = clusters, 
                              X = missing_until_now, 
                              fun = run_sim_wrapper)

a <- missing_rows(results_folder, name_pattern = "ResultsN1Row", check_numbers = 1:nrow_design, underscore = F)
b <- missing_rows(results_folder, name_pattern = "ResultsN2Row", check_numbers = 1:nrow_design, underscore = F)
missing_until_now <- intersect(a, b)

stopCluster(clusters)
