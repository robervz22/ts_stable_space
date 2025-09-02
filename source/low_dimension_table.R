#########################
# scripts and libraries #
#########################
remove(list = ls())
options(warn = -1)
library(data.table)
library(magrittr)
library(ggplot2)
library(latex2exp)
library(xtable)
library(patchwork)

source("./source/simulations.R")
source("./source/vectorial_methods.R")
source("./source/auxiliar_methods.R")

##############
# Parameters #
##############
m <- 11
r_values <- c(10,9,8)
i_values <- matrix(c(11,0,10,0,9,1,8,2),nrow=4,ncol=2,byrow = TRUE)
rownames(i_values) <- c('Case 1','Case 2','Case 3','Case 4')
Tt <- 100 # series length
S <- 500 # number of simulation
persistence <- "low" ; dist <- "t" # persistence and innovation process distribution
dependence <- TRUE
seeds <- c(1,1) # seeds for reproducibility
spca_sparse <- "penalty"  # type of sparsity: "penalty" or "varnum"
spca_para <- 0.25 # sparsity parameter for SPCA
spca_engine <- "elasticnet" # type of sparsity and engine for SPCA

###############################################
# Produce table for low-dimensional scenarios #
###############################################
df_low_dimension <- run_simulation(seeds,m,r_values,i_values,Tt,S,
                    dist = dist, persistence = persistence, dependence = dependence,
                    spca_sparse = spca_sparse, spca_para = spca_para, spca_engine = spca_engine)

                
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
dt_aux[, c("m", "Dimension", "X") := NULL]

# Order data
setorder(dt_aux, -r, Method)

# Format into a single long table
dt_table <- dt_aux[, .(Method, r, Case, Subspace)]

# Optionally reshape wide to combine cases into a single row per Method & r
dt_wide <- dcast(dt_table, Method + r ~ Case, value.var = "Subspace")
setorder(dt_wide,-r)

# Print as LaTeX table
print(xtable(dt_wide, align = paste0("lrl", paste(rep("r", ncol(dt_wide) - 2), collapse = ""))),
      include.rownames = FALSE)




