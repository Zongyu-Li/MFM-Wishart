# The code in this section is adapted from the code for Cappozzo & Casa. (2025)
# The code is adapted to allow shared degrees of freedom for all clusters
# See https://github.com/AndreaCappozzo/sparsemixwishart
# Andrea Cappozzo, Alessandro Casa,
# Model-based clustering for covariance matrices via penalized Wishart mixture models,
# Computational Statistics & Data Analysis, Volume 212, 2025, 108232

EM_controls <- function(tol = c(1e-03),
                        max_iter = 1e02,
                        type_start = c("hc", "random"),
                        linkage_hc_start= "ward.D",
                        # n_subset_start = NULL,
                        n_random_start = 50,
                        nu_max=NULL){
  list(tol = tol,
       max_iter = max_iter,
       type_start = match.arg( type_start, choices = eval(formals(EM_controls)$type_start) ),
       linkage_hc_start= linkage_hc_start,
       # n_subset_start = n_subset_start,
       n_random_start = n_random_start,
       nu_max=nu_max)
}


nu_k_root_finder_equation <- function(data, p, z_k, n_k, nu_k, Sigma_inv_k) {
  # Equation (5) of \cite{Hidot2010} Pattern Recognition letter
  sum(z_k * apply(data, 3, function(s)
    determinant(x = (s %*% Sigma_inv_k) / 2, logarithm = TRUE)$modulus)) - n_k * sum(digamma(x = (nu_k -
                                                                                                    (1:p) + 1) / 2))
}

weighted_sample_S_calculator <- function(data, z, p, n_K) {
  K <- ncol(z)
  weighted_data <- vector(mode = "list", length = K)
  S_weighted <- array(dim = c(p, p, K))
  for (k in 1:K) {
    # weighted_data[[k]] <- sweep(x = data,MARGIN = 3,STATS = z[,k],FUN = "*")
    weighted_data <- sweep(
      x = data,
      MARGIN = 3,
      STATS = z[, k],
      FUN = "*"
    )
    S_weighted[, , k] <- apply(weighted_data, MARGIN = c(1, 2), FUN = sum) /
      n_K[k]
  }
  S_weighted
  # out <- list(weighted_data=weighted_data, S_k_weighted=S_k_weighted)
  # out
}

log_dens_pro_calculator <- function(data, nu, Sigma, pro, N, K) {
  
  log_dens <- matrix(NA_real_, N, K)
  
  for (k in 1:K) {
    log_dens[, k] <- apply(
      data,
      MARGIN = 3,
      FUN = function(Omega) {
        LaplacesDemon::dwishart(
          Omega = Omega,
          nu = nu,
          S = Sigma[, , k],
          log = TRUE
        )
      }
    )
  }
  
  sweep(log_dens, 2, log(pro), "+")
}

Q_function_calculator <- function(data, nu, Sigma, pro, z, N, K,LAMBDA) {
  # Calculate log density
  log_dens_pro <- log_dens_pro_calculator(
    data = data,
    nu = nu,
    Sigma = Sigma,
    pro = pro,
    N = N,
    K = K
  )
  
  penalty_term <- sum(sapply(1:K, function(k)  sum(abs(LAMBDA[,,1]*Sigma[,,k]))))
  
  # Calculate Q function
  Q_function <- sum(log_dens_pro * z) - penalty_term
  
  return(Q_function)
}


nu_root_finder_equation_shared <- function(nu, data, p, z, Sigma_inv) {
  K <- dim(Sigma_inv)[3]
  lhs <- 0
  
  for (k in 1:K) {
    log_det_term_k <- apply(
      data,
      3,
      function(s) {
        as.numeric(
          determinant(
            x = (s %*% Sigma_inv[, , k]) / 2,
            logarithm = TRUE
          )$modulus
        )
      }
    )
    lhs <- lhs + sum(z[, k] * log_det_term_k)
  }
  
  rhs <- sum(z) * sum(digamma((nu - (1:p) + 1) / 2))
  lhs - rhs
}


M_step_sigma_and_nu <- function(data, z, n_K, p, K, nu, LAMBDA, tol, max_iter, N, pro, nu_max) {
  
  crit_Q_M <- TRUE
  iter_Q_M <- 0
  Q_M <- Q_M_prev <- -.Machine$double.xmax
  
  Sigma <- array(0, dim = c(p, p, K))
  Sigma_inv <- array(0, dim = c(p, p, K))
  
  while (crit_Q_M) {
    
    weighted_S <- weighted_sample_S_calculator(
      data = data,
      z = z,
      n_K = n_K,
      p = p
    )
    
    S_tilde <- weighted_S / nu
    
    for (k in 1:K) {
      covglasso_tmp <- covglasso::covglasso(
        S = S_tilde[, , k],
        n = n_K[k],
        lambda = 2 * LAMBDA / (n_K[k] * nu),
        penalize.diag = TRUE
      )
      
      Sigma[, , k] <- covglasso_tmp$sigma
      Sigma_inv[, , k] <- covglasso_tmp$omega
    }
    
    nu <- tryCatch(
      stats::uniroot(
        f = nu_root_finder_equation_shared,
        interval = c((p - 1) + sqrt(.Machine$double.eps), nu_max),
        data = data,
        p = p,
        z = z,
        Sigma_inv = Sigma_inv
      )$root,
      error = function(e) {
        stop("Error in numerically finding shared nu; try increasing nu_max", call. = FALSE)
      }
    )
    
    nu <- max(nu, p)
    
    Q_M <- Q_function_calculator(data, nu, Sigma, pro, z, N, K, LAMBDA)
    
    err_Q_M <- abs(Q_M - Q_M_prev) / (1 + abs(Q_M))
    Q_M_prev <- Q_M
    
    crit_Q_M <- (err_Q_M > tol & iter_Q_M < max_iter)
    iter_Q_M <- iter_Q_M + 1
  }
  
  list(Sigma = Sigma, nu = nu)
}


E_step_z_log_lik <- function(data, nu, Sigma, pro, N, K) {
  
  # Calculate log density with proportionality constant
  log_dens_pro <- log_dens_pro_calculator(
    data = data,
    nu = nu,
    Sigma = Sigma,
    pro = pro,
    N = N,
    K = K
  )
  
  # Calculate the maximum log density for each row
  zMax <- apply(log_dens_pro, 1, max)
  
  # Compute the log-likelihood for each obs
  loglik_i <- zMax + log(rowSums(exp(log_dens_pro - zMax)))
  
  # Update the responsibility matrix
  z <- exp(log_dens_pro - loglik_i)
  
  return(list(z = z, loglik_i = loglik_i))
}


