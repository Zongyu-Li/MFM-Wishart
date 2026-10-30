fit_moewishart_bic <- function(
    data,
    K = 1:10,
    niter = 3000,
    estimate_nu = TRUE,
    n_restarts = 3,
    restart_iters = 20,
    tol = 1e-6,
    verbose = FALSE,
    z_true = NULL
) {
  
  S_list <- asplit(data, 3L)
  
  fits <- lapply(K, function(k) {
    
    tryCatch(
      moewishart::mixturewishart(
        S_list = S_list,
        K = k,
        niter = niter,
        method = "em",
        init_pi = NULL,
        init_nu = NULL,
        init_Sigma = NULL,
        estimate_nu = estimate_nu,
        n_restarts = n_restarts,
        restart_iters = restart_iters,
        tol = tol,
        verbose = verbose
      ),
      
      error = function(e) {
        stop(
          sprintf(
            "moewishart failed for K=%s: %s",
            k, conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
  })
  
  names(fits) <- paste0("K_", K)
  
  BIC_values <- vapply(
    fits,
    function(fit) {
      as.numeric(moewishart::computeIC(fit)$BIC)
    },
    numeric(1)
  )
  
  if (any(!is.finite(BIC_values))) {
    stop(
      "Non-finite BIC for K: ",
      paste(K[!is.finite(BIC_values)], collapse = ", "),
      call. = FALSE
    )
  }
  
  classifications <- lapply(seq_along(K), function(j) {
    
    tau <- fits[[j]]$tau
    
    if (!is.matrix(tau) || !is.numeric(tau) ||
        nrow(tau) != n || ncol(tau) != K[j] ||
        any(!is.finite(tau))) {
      stop(
        sprintf("Invalid tau matrix for K=%s.", K[j]),
        call. = FALSE
      )
    }
    
    max.col(tau)
  })
  
  ARI_values <- rep(NA_real_, length(K))
  
  if (!is.null(z_true)) {
    ARI_values <- vapply(
      classifications,
      function(z) {
        mclust::adjustedRandIndex(z_true, z)
      },
      numeric(1)
    )
  }
  
  best_idx <- which.min(BIC_values)
  
  comparison <- data.frame(
    K = K,
    BIC = unname(BIC_values),
    ARI = unname(ARI_values),
    selected = seq_along(K) == best_idx
  )
  
  list(
    K = K[best_idx],
    classification = classifications[[best_idx]],
    BIC = unname(BIC_values[best_idx]),
    ARI = unname(ARI_values[best_idx]),
    summary = comparison,
    best_fit = fits[[best_idx]],
    fits = fits
  )
}