#########################
# scripts and libraries #
#########################
remove(list = ls())
options(warn = -1)
pacman::p_load(data.table,magrittr,ggplot2,latex2exp,xtable,patchwork)
source("./source/vectorial_methods.R")
source("./source/auxiliar_methods.R") 
source("./source/simulations.R")

##############
# Parameters #
##############
m <- 11
r_values <- c(10,9,8)
i_values <- matrix(c(11,0,10,0,9,1,8,2),nrow=4,ncol=2,byrow = TRUE)
rownames(i_values) <- c('Case 1','Case 2','Case 3','Case 4')
Tt <- 100 # series length
S <- 300 # number of simulation
persistence <- "low" ; dist <- "t" # persistence and innovation process distribution
dependence <- TRUE
seeds <- c(1,1) # seeds for reproducibility
spca_sparse <- "varnum"  # type of sparsity: "penalty" or "varnum"
spca_engine <- "elasticnet" # type of sparsity and engine for SPCA
spca_eta <- 0.6
spls_eta <- spca_eta  # e.g. reuse SPCA penalty; or set manually, e.g. 0.6

###############################################
# Produce table for low-dimensional scenarios #
###############################################
df_low_dimension <- run_simulation(seeds, m, r_values, i_values, Tt, S,
                    dist = dist, persistence = persistence, dependence = dependence,
                    spca_sparse = spca_sparse, spca_engine = spca_engine,
                    spca_eta = spca_eta, spls_eta = spls_eta,
                    parallel = TRUE,  # Enable parallelization
                    n_cores = NULL)   # Use all available cores (minus 1)

                
dt_low_dimension <- as.data.table(df_low_dimension)
dt_low_dimension[, mean_n_coint := sprintf("%.3f (%.3f)", mean_n_coint, sd_n_coint)]
dt_low_dimension[, mean_n_norms := sprintf("%.3f (%.3f)", mean_n_norms, sd_n_norms)]

# Remove SD columns as they are now embedded
dt_low_dimension[, c("sd_n_coint", "sd_n_norms") := NULL]
dt_low_dimension[, c("i1", "i2") := NULL]
setnames(dt_low_dimension, old = c("mean_n_coint","mean_n_norms"),
                     new = c("Dimension", "Subspace"))


# Prepare data
dt_aux <- copy(dt_low_dimension)
# dt_aux[, c("m", "Dimension", "X") := NULL]
dt_aux[, c("m", "X") := NULL]

# Order data
setorder(dt_aux, -r, Method)

# Format into a single long table
# dt_table <- dt_aux[, .(Method, r, Case, Subspace)]
dt_table <- dt_aux[, .(Method, r, Case, Subspace, Dimension)]

# Optionally reshape wide to combine cases into a single row per Method & r
dt_wide_subspace <- dcast(dt_table, Method + r ~ Case, value.var = "Subspace")
dt_wide_dimension <- dcast(dt_table, Method + r ~ Case, value.var = "Dimension")

setorder(dt_wide_subspace,-r)
setorder(dt_wide_dimension,-r)

# Print as LaTeX table
print(xtable(dt_wide_subspace, align = paste0("lrl", paste(rep("r", ncol(dt_wide_subspace) - 2), collapse = ""))),
      include.rownames = FALSE)

print(xtable(dt_wide_dimension, align = paste0("lrl", paste(rep("r", ncol(dt_wide_dimension) - 2), collapse = ""))),
      include.rownames = FALSE)
