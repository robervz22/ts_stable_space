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
# parameters #
##############
Tt <- 100 # series length
values_pairs <- list(c(6,4), c(11,9)) # series' dimension (list makes indexing easier)
case_values <- c(1,2) # scenario to consider 1 cointegration, 2 mixed
S <- 500 # number of simulations
crit <- "SC(n)" # criterio Johansen estimation
ec_det <- "none" # non trend and seasonality for Johansen procedure
pct <- "5pct" # Johansen's significance level
test <- "kpss" # test to use in PLS o PCA options: "adf" and "kpss"
persistence <- "low" ; dist <- "t" # persistence and innovation process distribution
trend <- FALSE # consider a simple linear trend
dependence <- TRUE # if we consider a covariance structure in the innovation process
seeds <- c(1,1) # seeds for reproducibility
sigma <- 1 # equal variance for the non-correlated case
spca_sparse <- "penalty"  # type of sparsity: "penalty" or "varnum"
spca_para <- 0.25 # sparsity parameter for SPCA
spca_engine <- "elasticnet" # type of sparsity and engine for SPCA
methods_colors <- c("Johansen" = "#E69F00","PLS" = "#56B4E9", "PCA" = "#009E73", "SPCA" = "#CC79A7") 

##################
# Collect plots  #
##################
bar_plots <- list()
violin_plots <- list()
plot_idx <- 1

for (pair in values_pairs) {
  ###############################
  # Theoretical VAR(1) structure#
  ###############################
  m <- pair[1]; r <- pair[2]

  # Shared Sigma_eps
  set.seed(seeds[1])  # for reproducibility
  if (dependence) {
    Sigma_eps <- clusterGeneration::genPositiveDefMat(m, covMethod = "eigen")$Sigma
  } else {
    Sigma_eps <- sigma * diag(1, m, m)
  }

  # Shared stable basis
  set.seed(seeds[2] + r)  # shared seed for beta, beta_orto, alpha, gamma
  I_m <- qr.Q(qr(matrix(rnorm(m^2), m, m)))
  beta <- I_m[, 1:r]
  beta_orto <- if (r < m) I_m[, (r + 1):m] else matrix(0, m, m)

  # Alpha and gamma
  alpha <- if (persistence == "low") runif(r, 0.1, 0.3) else runif(r, 0.3, 0.7)
  gamma <- runif(ncol(beta_orto), -0.7, 0.7)

  # ---- IMPORTANT: iterate over 'case' INSIDE the pair loop ----
  for (case in case_values) {

    # cointegration / mixed orders by case
    orders <- switch(
      as.character(case),
      "1" = c(m, 0),
      "2" = c(m - 3, 2),
      stop("Unsupported case")
    )
    i1 <- orders[1]; i2 <- orders[2]

    ###################
    # Simulation loop #
    ###################
    n_coint <- matrix(NA, S, 4)
    n_norms <- matrix(NA, S, 4)
    colnames(n_coint) <- c("Johansen", "PLS", "PCA", "SPCA")
    colnames(n_norms) <- c("Johansen", "PLS", "PCA", "SPCA")

    for (s in 1:S) {
      # Epsilon innovations
      set.seed(s)
      epsilon <- if (dist == "normal") {
        if (trend) {
          portes::varima.sim(n = Tt + 100, k = m, sigma = Sigma_eps, trend = rep(0.05, m), demean = rep(0, m))
        } else {
          portes::varima.sim(n = Tt + 100, k = m, sigma = Sigma_eps, demean = rep(0, m))
        }
      } else if (dist == "t") {
        if (trend) {
          portes::varima.sim(n = Tt + 100, k = m, sigma = Sigma_eps, trend = rep(0.05, m), innov.dist = "t", dft = 3, demean = rep(0, m))
        } else {
          portes::varima.sim(n = Tt + 100, k = m, sigma = Sigma_eps, innov.dist = "t", dft = 3, demean = rep(0, m))
        }
      }

      # Simulate X using the improved core
      X_sim <- X_simulation(seed = s,
                            m = m, r = r, Tt = Tt + 200,
                            beta = beta, beta_orto = beta_orto,
                            alpha = alpha, gamma = gamma, Sigma_eps = Sigma_eps,
                            mix = TRUE, dist = dist, trend = trend)(i1, i2)

      X <- X_sim$X[201:nrow(X_sim$X), ]  # remove burn-in

      # PLS, PCA, SPCA
      basis_PLS <- basis_stable(X, method = "pls", test = test)
      basis_PCA <- basis_stable(scale(X), method = "pca", test = test)
      basis_SPCA <- basis_stable(scale(X), method = "spca", test = test, 
                                    spca_sparse = spca_sparse, 
                                    spca_engine = spca_engine, 
                                    spca_para = spca_para)

      if (!is.null(ncol(basis_PLS$basis_S))){ 
        n_coint[s, "PLS"] <- ncol(basis_PLS$basis_S) - r
      }
      else {
        n_coint[s, "PLS"] <- -r
      }

      if (!is.null(ncol(basis_PCA$basis_S))) {
        n_coint[s, "PCA"] <- ncol(basis_PCA$basis_S) - r
      }
      else{
        n_coint[s, "PCA"] <- -r
      }

      if (!is.null(ncol(basis_SPCA$basis_S))) {
        n_coint[s, "SPCA"] <- ncol(basis_SPCA$basis_S) - r
      }
      else{
        n_coint[s, "SPCA"] <- -r
      }
      n_norms[s, "PLS"] <- grassmann_distance(beta, basis_PLS$basis_S)
      n_norms[s, "PCA"] <- grassmann_distance(beta, basis_PCA$basis_S)
      n_norms[s, "SPCA"] <- grassmann_distance(beta, basis_SPCA$basis_S)


      # Johansen
      if (m <= 11) {
        basis_johansen <- basis_stable(X, method = "johansen", ec_det = ec_det)
        if (!is.null(ncol(basis_johansen$basis_S))) {
          n_coint[s, "Johansen"] <- ncol(basis_johansen$basis_S) - r
        }
        else {
          n_coint[s, "Johansen"] <- -r
        }
        n_norms[s, "Johansen"] <- grassmann_distance(beta, basis_johansen$basis_S)
      }
    }

    ###########
    # Barplot #
    ###########
    data_n_coint <- as.data.table(n_coint)
    data_n_coint <- data_n_coint[, which(colSums(is.na(data_n_coint)) == 0), with = FALSE]

    data_n_coint_melt <- melt(
      data_n_coint,
      variable.name = "Method",
      value.name   = "AnnError"
    )
    data_n_coint_melt[, AnnError := as.integer(AnnError)]

    gg_barplot <- ggplot(data_n_coint_melt, aes(x = AnnError, fill = Method)) +
      geom_bar(alpha = 0.7) +
      facet_grid(~ Method) +
      ggtitle(paste0("Scenario ", case, " (m = ", m, ", r = ", r, ")")) +
      xlab(TeX(r"($\hat{r}-r$)")) + mytheme +
      scale_x_continuous(breaks = seq(min(data_n_coint_melt$AnnError, na.rm = TRUE),
                                      max(data_n_coint_melt$AnnError, na.rm = TRUE), by = 1)) +
      ylab("Counts") +
      scale_fill_manual(values = methods_colors) +
      theme(legend.position = "none")

    ###############
    # Violin plot #
    ###############
    data_norms <- as.data.table(n_norms)
    data_norms_melt <- melt(data_norms, variable.name = "Method", value.name = "Norm")

    gg_violin <- ggplot(data_norms_melt, aes(x = Method, y = Norm, fill = Method)) +
      geom_violin(alpha = 0.7) +
      labs(
        title = paste0("Scenario ", case, " (m = ", m, ", r = ", r, ")"),
        x = "Method",
        y = TeX(r"($\delta(\hat{\beta},\beta)$)")
      ) + mytheme +
      scale_fill_manual(values = methods_colors) +
      theme(legend.position = "none")

    # store plots
    bar_plots[[plot_idx]]    <- gg_barplot
    violin_plots[[plot_idx]] <- gg_violin

    # save individual PDFs (square, e.g., 5in × 5in)
    # ggsave(
    #   filename = file.path(out_dir, sprintf("barplot_case%d_m%d_r%d.pdf", case, m, r)),
    #   plot     = gg_barplot,
    #   width    = 5, height = 5, units = "in", device = "pdf", dpi = 300, useDingbats = FALSE
    # )
    # ggsave(
    #   filename = file.path(out_dir, sprintf("violin_case%d_m%d_r%d.pdf", case, m, r)),
    #   plot     = gg_violin,
    #   width    = 5, height = 5, units = "in", device = "pdf", dpi = 300, useDingbats = FALSE
    # )

    plot_idx <- plot_idx + 1
  }
}

#####################################
# 2×2 grids (patchwork) + save PDF  #
#####################################
# if you ever vary lengths, pad/trim accordingly; here we expect 4 plots each
bar_grid    <- wrap_plots(bar_plots,    ncol = 2)
violin_grid <- wrap_plots(violin_plots, ncol = 2)

# Save the grids as square PDFs (e.g., 8in × 8in)
ggsave(
  filename = "./images/Figure_1.pdf",
  plot     = bar_grid,
  width    = 32, height = 16, units = "in", device = "pdf", dpi = 300, useDingbats = FALSE
)
ggsave(
  filename = "./images/Figure_2.pdf",
  plot     = violin_grid,
  width    = 32, height = 16, units = "in", device = "pdf", dpi = 300, useDingbats = FALSE
)




