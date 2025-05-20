###############
# simulations #
###############

# Grassman distance for subspaces with different dimensions
grassmann_distance <- function(S1, S2) {
  # Ensure the inputs are matrices
  if (!is.matrix(S1) || !is.matrix(S2)) {
    stop("Inputs must be matrices (orthonormal bases).")
  }

  n1 <- ncol(S1)
  n2 <- ncol(S2)
  d <- min(n1, n2)
  k <- abs(n1 - n2)

  # Check dimensional compatibility
  if (nrow(S1) != nrow(S2)) {
    stop("Both bases must be in the same ambient space (same number of rows).")
  }

  # Compute the principal angles using SVD
  svd_result <- svd(t(S1) %*% S2)
  cos_theta <- svd_result$d[1:d]
  cos_theta <- pmin(pmax(cos_theta, -1), 1)  # numerical stability
  theta <- acos(cos_theta)

  # Pad with π/2 if dimensions differ
  if (k > 0) {
    theta <- c(theta, rep(pi / 2, k))
  }

  # Grassmann distance
  d_grassmann <- sqrt(sum(theta^2))
  return(d_grassmann)
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
      portes::varima.sim(n = Tt, k = m, sigma = Sigma_eps, innov.dist = "t", dft = 3, demean = rep(0, m))
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
                           trend = FALSE) {
  results_list <- list()
  methods <- c("Johansen", "PLS", "PCA")
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
        sim_fun <- X_simulation_new(
          seed = s, m = m, r = r, Tt = Tt + 100,
          alpha = alpha, beta = beta, beta_orto = beta_orto, gamma = gamma,
          Sigma_eps = Sigma_eps, dist = dist, trend = trend, mix = TRUE
        )
        sim_data <- sim_fun(i1, i2)

        X <- ts(scale(sim_data$X[101:(Tt + 100), ]))
        XX <- X[1:(nrow(X) - 1), ]; Y <- X[2:nrow(X), ]

        beta_teo <- beta

        iter_results <- data.frame(Method = methods, m = m, r = r,
                                   i1 = i1, i2 = i2, Case = case_name,
                                   n_coint = NA, n_norms = NA)

        # Estimate cointegration basis with all 3 methods
        basis_PLS <- basis_stable(X, method = "pls", test = test)
        basis_PCA <- basis_stable(XX, method = "pca", test = test)

        for (method in c("PLS", "PCA")) {
          basis <- get(paste0("basis_", method))
          if (!is.null(ncol(basis$basis_S))) {
            iter_results[iter_results$Method == method, "n_coint"] <- ncol(basis$basis_S) - r
            iter_results[iter_results$Method == method, "n_norms"] <- grassmann_distance(beta_teo, basis$basis_S)
          }
        }

        if (m <= 11) {
          method <- "Johansen"
          basis_johansen <- basis_stable(XX, method = "johansen")
          if (!is.null(ncol(basis_johansen$basis_S))) {
            iter_results[iter_results$Method == method, "n_coint"] <- ncol(basis_johansen$basis_S) - r
            iter_results[iter_results$Method == method, "n_norms"] <- grassmann_distance(beta_teo, basis_johansen$basis_S)
          }
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




