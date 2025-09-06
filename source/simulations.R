###############
# simulations #
###############

# Grassmann distance for subspaces (allows different dimensions and NULL)
grassmann_distance <- function(S1, S2) {
  
  # Treat NULL as zero-dimensional subspaces
  if (is.null(S1) && is.null(S2)) return(0)

  if (is.null(S1)) {
    if (!is.matrix(S2)) stop("S2 must be a matrix or NULL.")
    return(sqrt(ncol(S2)) * (pi / 2))
  }
  if (is.null(S2)) {
    if (!is.matrix(S1)) stop("S1 must be a matrix or NULL.")
    return(sqrt(ncol(S1)) * (pi / 2))
  }

  # Ensure the inputs are matrices
  if (!is.matrix(S1) || !is.matrix(S2)) {
    stop("Inputs must be matrices (orthonormal bases) or NULL.")
  }

  # Ambient dimension check
  if (nrow(S1) != nrow(S2)) {
    stop("Both bases must be in the same ambient space (same number of rows).")
  }

  n1 <- ncol(S1)
  n2 <- ncol(S2)
  d  <- min(n1, n2)
  k  <- abs(n1 - n2)

  # If one (or both) subspaces are zero-dimensional, skip SVD
  if (d == 0) {
    return(sqrt(k) * (pi / 2))
  }

  # Compute principal angles via SVD of cross-basis matrix
  C <- t(S1) %*% S2
  svd_result <- svd(C, nu = 0, nv = 0)
  cos_theta <- svd_result$d[1:d]
  cos_theta <- pmin(pmax(cos_theta, -1), 1)  # numerical stability
  theta <- acos(cos_theta)

  # Pad with π/2 for the dimension mismatch
  if (k > 0) {
    theta <- c(theta, rep(pi / 2, k))
  }

  return(sqrt(sum(theta^2)))
}



# Optimised simulation function
X_simulation <- function(seed, m, r, Tt,
                         alpha, beta, beta_orto, gamma,
                         Sigma_eps, dist = "normal",
                         trend = FALSE, mix = FALSE) {
  Z0 <- rep(0, m)
  D_alpha <- diag(alpha, ncol(beta), ncol(beta))
  D_gamma <- diag(gamma, ncol(beta_orto), ncol(beta_orto))

  # Generate innovations
  set.seed(seed)
  epsilon <- if (dist == "normal") {
    if (trend) {
      portes::varima.sim(n = Tt, k = m, sigma = Sigma_eps, trend = rep(0.05, m), demean = rep(0, m))
    } else {
      portes::varima.sim(n = Tt, k = m, sigma = Sigma_eps, demean = rep(0, m))
    }
  } else if (dist == "t") {
    if (trend) {
      portes::varima.sim(n = Tt, k = m, sigma = Sigma_eps, trend = rep(0.05, m), innov.dist = "t", dft = 3, demean = rep(0, m))
    } else {
      p###############
# simulations #
###############

# Grassmann distance for subspaces (allows different dimensions and NULL)
grassmann_distance <- function(S1, S2) {
  ortes::varima.sim(n = Tt, k = m, sigma = Sigma_eps, innov.dist = "t", dft = 3, demean = rep(0, m))
    }
  }

  # Preallocate
  ZS_sim <- matrix(0, Tt, m)
  ZN_sim <- matrix(0, Tt, m)
  Zaux <- Z0

  # Efficient loop
  for (t in 1:Tt) {
    ZS <- beta %*% D_alpha %*% t(beta) %*% Zaux
    ZN <- beta_orto %*% D_gamma %*% t(beta_orto) %*% Zaux
    Zaux <- ZS + ZN + epsilon[t, ]
    ZS_sim[t, ] <- ZS
    ZN_sim[t, ] <- ZN
  }

  # Integration order mix
  integration_order <- function(i1, i2) {
    X_aux <- ZN_sim
    if (i1 > 0) {
      X_aux[, 1:i1] <- apply(as.matrix(ZN_sim[, 1:i1]), 2, cumsum)
    }
    if (i2 > 0) {
      X_aux[, (i1 + 1):(i1 + i2)] <- apply(apply(as.matrix(ZN_sim[, (i1 + 1):(i1 + i2)]), 2, cumsum), 2, cumsum)
    }
    X_final <- ts(X_aux + ZS_sim)
    list(X = X_final, beta = beta, beta_orto = beta_orto, Sigma = Sigma_eps)
  }

  if (mix) {
    return(integration_order)
  } else {
    X_final <- ts(apply(ZN_sim, 2, cumsum) + ZS_sim)
    return(list(X = X_final, beta = beta, beta_orto = beta_orto, Sigma = Sigma_eps))
  }
}

# Optimised simulation study
run_simulation <- function(seeds, m, r_values, i_values, Tt, S,
                           test = "kpss", dist = "normal",
                           persistence = "low", dependence = FALSE,
                           trend = FALSE,burnin = 200, methods = c("Johansen", "PLS", "PCA", "SPCA"),                          
                           # --- new: spca controls with safe defaults ---
                            spca_engine = c("elasticnet", "PMA"),
                            spca_sparse = c("varnum", "penalty"), # only for elasticnet
                            spca_k = NULL, # number of sparse components (defaults to ncol(X))
                            spca_para = NULL, # sparsity control (engine-specific; see spca_alg docs below)
                            spca_center = FALSE, spca_scale = FALSE) {

  
  results_list <- list()
  i1_i2_names <- rownames(i_values)


  set.seed(seeds[1])  # ensure reproducibility
  Sigma_eps <- if (dependence) {
    clusterGeneration::genPositiveDefMat(m, covMethod = "eigen")$Sigma
  } else {
    diag(1, m)
  }

  for (r in r_values) {
    set.seed(seeds[2] + r)  # reproducible beta/beta_orto
    I_m <- qr.Q(qr(matrix(rnorm(m^2), m, m)))
    beta <- as.matrix(I_m[, 1:r])
    beta_orto <- if (r < m) as.matrix(I_m[, (r + 1):m]) else matrix(0, m, m)

    alpha <- switch(persistence,
                    "low" = runif(r, 0.1, 0.3),
                    "high" = runif(r, 0.3, 0.7))
    gamma <- runif(ncol(beta_orto), -0.7, 0.7)


    for (idx in 1:nrow(i_values)) {
      pair <- i_values[idx, ]
      i1 <- pair[1]; i2 <- pair[2]
      case_name <- i1_i2_names[idx]

      for (s in 1:S) {
        # Simulate with precomputed components
        start_simulation <- Sys.time()
        sim_fun <- X_simulation(
          seed = s, m = m, r = r, Tt = Tt + burnin,
          alpha = alpha, beta = beta, beta_orto = beta_orto, gamma = gamma,
          Sigma_eps = Sigma_eps, dist = dist, trend = trend, mix = TRUE
        )
        sim_data <- sim_fun(i1, i2)

        X <- sim_data$X[(burnin+1):nrow(sim_data$X), ]  # remove burn-in
        end_simulation <- Sys.time()
        elapsed_simulation <- difftime(end_simulation, start_simulation, units = "secs")

        cat("\nSimulation time: ", round(elapsed_simulation), "(s)")

        iter_results <- data.frame(Method = methods, m = m, r = r,
                                   i1 = i1, i2 = i2, Case = case_name,
                                   n_coint = NA, n_norms = NA)

        for (method in intersect(c("PLS", "PCA", "SPCA"),methods)) {
          start_method <- Sys.time()
          if(method == "PLS"){
            basis <- basis_stable(X, method = "pls", test = test)
          } else if(method == "PCA"){
            basis <- basis_stable(scale(X), method = "pca", test = test)
          } else if(method == "SPCA"){
            basis <- basis_stable(scale(X), method = "spca", test = test, 
                                      spca_sparse = spca_sparse, 
                                      spca_engine = spca_engine, 
                                      spca_para = spca_para)
          }
          end_method <- Sys.time()
          elapsed_method <- difftime(end_method, start_method, units = "secs")
          cat("\n",method, " time: ",round(elapsed_method,2),"(s)")

          if (!is.null(ncol(basis$basis_S))) {
            iter_results[iter_results$Method == method, "n_coint"] <- ncol(basis$basis_S) - r
          } else {
            iter_results[iter_results$Method == method, "n_coint"] <- -r
          }
          iter_results[iter_results$Method == method, "n_norms"] <- grassmann_distance(beta, basis$basis_S)
        }

        if (m <= 11 && "Johansen" %in% methods) {
          basis_johansen <- basis_stable(X, method = "johansen")
          if (!is.null(ncol(basis_johansen$basis_S))) {
            iter_results[iter_results$Method == method, "n_coint"] <- ncol(basis_johansen$basis_S) - r
          } else {
            iter_results[iter_results$Method == method, "n_coint"] <- -r
          }
          iter_results[iter_results$Method == method, "n_norms"] <- grassmann_distance(beta, basis_johansen$basis_S)

        }
        results_list[[length(results_list) + 1]] <- iter_results
      }
    }
  }

  # Aggregate simulation results
  results_df <- do.call(rbind, results_list)

  summary_table <- results_df %>%
    dplyr::group_by(Method, m, r, i1, i2, Case) %>%
    dplyr::summarise(
      mean_n_coint = mean(n_coint, na.rm = TRUE),
      sd_n_coint = sd(n_coint, na.rm = TRUE),
      mean_n_norms = mean(n_norms, na.rm = TRUE),
      sd_n_norms = sd(n_norms, na.rm = TRUE),
      .groups = "drop"
    )

  return(summary_table)
}




