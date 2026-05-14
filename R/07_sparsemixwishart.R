# The code in this section is adapted from the code for Cappozzo & Casa. (2025)
# The code is adapted to allow shared degrees of freedom for all clusters
# See https://github.com/AndreaCappozzo/sparsemixwishart
# Andrea Cappozzo, Alessandro Casa,
# Model-based clustering for covariance matrices via penalized Wishart mixture models,
# Computational Statistics & Data Analysis, Volume 212, 2025, 108232

em_sparse_mix_wishart <- function(data,
                                  K,
                                  penalty,
                                  P,
                                  control,
                                  data_dim,
                                  hc_init) {
  p <- data_dim[1]
  N <- data_dim[3]
  
  P <- if(is.null(P)){
    matrix(data = 1,nrow = p,ncol = p)
  } else{
    P
  }
  
  LAMBDA <- array(P*penalty,dim = c(p,p,1)) # an array of dims p x p x 1, to be used with the covglasso package
  
  # store EM parameters
  tol <- control$tol
  max_iter <- control$max_iter
  n_random_start <-  control$n_random_start
  
  # Define upper bound for nu uniroot search
  nu_max <- ifelse(is.null(control$nu_max), max(N,3*p), control$nu_max)
  
  # initialization of z -----------------------------------------------------------------
  class_init <- if(is.null(hc_init)) {
    sample(1:K, size = N, replace = TRUE)
  } else {
    stats::cutree(hc_init, k = K)
  }
  z <- mclust::unmap(class_init)
  nu <- p
  n_K <- colSums(z)
  weighted_S <- weighted_sample_S_calculator(data = data,z = z,n_K = n_K,p = p)
  S_tilde <- sweep(weighted_S, MARGIN = 3,STATS = nu,FUN = "/")
  Sigma <- Sigma_inv <- array(dim = c(p,p,K))
  
  for (k in 1:K) {
    Sigma[, , k] <- covglasso::covglasso(
      S = S_tilde[, , k],
      n = n_K[k], # not used
      lambda = 2 * LAMBDA / (n_K[k] * nu),
      penalize.diag = TRUE
    )$sigma
  }
  
  
  iter <- 0
  loglik_pen <- loglik_pen_prev <- -.Machine$integer.max / 2
  loglik_pen_vec <- NULL
  crit <- TRUE
  
  while (crit) {
    
    # M step ------------------------------------------------------------------
    
    #tau
    n_K <- colSums(z)
    pro <- n_K / N
    
    # Sigma and nu
    m_step_output_Sigma_nu <- M_step_sigma_and_nu(data, z, n_K, p, K, nu, LAMBDA, tol, max_iter, N, pro, nu_max)
    
    Sigma <- m_step_output_Sigma_nu$Sigma
    nu <- m_step_output_Sigma_nu$nu
    
    # E step ------------------------------------------------------------------
    e_step_output <- E_step_z_log_lik(data, nu, Sigma, pro, N, K)
    z <- e_step_output$z
    
    # Log lik and objective function calculation ------------------------------
    loglik <- sum(e_step_output$loglik_i)
    
    penalty_term <- sum(sapply(1:K, function(k)  sum(abs(LAMBDA[,,1]*Sigma[,,k]))))
    
    loglik_pen <- loglik - penalty_term
    
    # Check convergence -------------------------------------------------------
    
    err <- abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
    loglik_pen_prev <- loglik_pen
    loglik_pen_vec <- c(loglik_pen_vec, loglik_pen)
    iter <- iter + 1
    crit <- (err > tol & iter < max_iter)
    
  }
  
  # Results collection
  
  parameters <- list(pro=pro, Sigma=Sigma, nu=nu)
  n_par_pro <- K - 1
  n_par_nu <- 1
  n_par_Sigma <- p * K + sum(apply(Sigma,
                                   3, function(A) A[upper.tri(A)] != 0))
  bic_final <- 2 * loglik - (n_par_pro + n_par_nu + n_par_Sigma) * log(N)
  
  OUT <-
    list(
      loglik = loglik,
      loglik_pen = loglik_pen,
      parameters = parameters,
      z = z,
      K=K,
      classification = mclust::map(z),
      bic = bic_final,
      n_params = c(pro = n_par_pro, nu = n_par_nu, Sigma = n_par_Sigma),
      penalty_term = penalty_term,
      obj_func_trace = loglik_pen_vec,
      iter = iter
    )
  
  return( OUT )
}

## ----- Wrapper function -----
sparsemixwishart <- function(data,
                             K = 2,
                             penalty = 10,
                             P=NULL, # arbitrary p x p matrix with nonnegative elements (see Sparse estimation of a covariance matrix Bien 2011)
                             control = EM_controls(),
                             verbose = interactive()) {
  
  # The best model is the one that maximizes the BIC
  
  call <- match.call()
  
  start_time <- Sys.time()
  
  data_dim <- dim(data)
  
  p <- data_dim[1]
  N <- data_dim[3]
  
  all_hyperparameters <-
    expand.grid(K = K,
                penalty = penalty)
  
  n_different_models <- nrow(all_hyperparameters)
  
  models_container <-
    vector(mode = "list", length = n_different_models)
  
  if (verbose) {
    cat("Fitting:\n")
    utils::flush.console()
    pbar <- utils::txtProgressBar(min = 0,
                                  max = n_different_models,
                                  style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  
  # Initialization performed via hierarchical clustering calcolating distances between covariance matrices, following
  # Non-Euclidean statistics for covariance matrices, with applications to diffusion tensor imaging (Dryden2009)
  hc_init <- if ( control$type_start == "hc" ) {
    data_vec <- matrix(data, N, p*p, byrow = TRUE)
    Dist_cov_init <- stats::dist(data_vec, method = "euclidean") # Euclidean metric is considered
    stats::hclust(d = Dist_cov_init, method = control$linkage_hc_start)
  } else NULL
  
  for (model in 1:n_different_models) {
    models_container[[model]] <-
      tryCatch(em_sparse_mix_wishart(
        data = data,
        K = all_hyperparameters[model, "K"],
        penalty = all_hyperparameters[model, "penalty"],
        P = P,
        control = control,
        data_dim=data_dim,
        hc_init=hc_init
      ),error = function(e) {
        list(bic = NA)
      })
    
    if (verbose) {
      ipbar <- ipbar + 1
      utils::setTxtProgressBar(pbar, ipbar)
    }
  }
  
  models_BICS <- sapply(models_container, "[[", "bic")
  max_bic_model <- which.max(models_BICS)
  
  if(length(max_bic_model)==0){ # Exiting the function if no model converged
    return(stop("None of the considered models converged",call. = FALSE))
  }
  
  selected_model <- models_container[[max_bic_model]]
  selected_model$BIC <- cbind(all_hyperparameters, bic=models_BICS)
  
  end_time <- Sys.time()
  selected_model$elapsed_time <- end_time - start_time
  
  if (verbose) {
    cat(
      "\nModel with K=",
      selected_model$K,
      ", penalty=",
      all_hyperparameters[max_bic_model, "penalty"],
      " returns the highest BIC.",
      sep = ""
    )
  }
  selected_model
}