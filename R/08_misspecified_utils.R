# Helper functions for simulations in mis-specified cases

## finite-sample covariance-level effective sample size
neff_cov_var1 <- function(T, phi){
  if(T < 2) stop("T must be >= 2.")
  r <- phi^2
  infl <- 1 + 2 * sum((1 - (1:(T-1))/T) * r^(1:(T-1)))
  T / infl
}

## choose T so that neff_cov_var1(T, phi) is close to target_df
choose_T_var1 <- function(phi, target_df, T_max = 500){
  Ts <- 2:T_max
  neff_vals <- sapply(Ts, neff_cov_var1, phi = phi)
  Ts[which.min(abs(neff_vals - target_df))]
}

## generate n_draws SPD matrices from stationary VAR(1),
## then rescale to match Wishart(target_df, Sigma) mean level
rVAR1_scatter_match <- function(n_draws, Sigma, phi, T_len, target_df){
  p <- nrow(Sigma)
  R <- chol(Sigma)  # upper triangular, t(R)%*%R = Sigma
  out <- array(0, dim = c(p, p, n_draws))
  
  innov_scale <- sqrt(1 - phi^2)
  
  for(m in 1:n_draws){
    X <- matrix(0, nrow = T_len, ncol = p)
    
    ## stationary start: X1 ~ N(0, Sigma)
    X[1, ] <- drop(rnorm(p) %*% R)
    
    ## VAR(1): Xt = phi Xt-1 + eps_t, eps_t ~ N(0, (1-phi^2) Sigma)
    if(T_len > 1){
      for(tt in 2:T_len){
        eps_t <- innov_scale * drop(rnorm(p) %*% R)
        X[tt, ] <- phi * X[tt - 1, ] + eps_t
      }
    }
    
    ## raw scatter: sum Xt Xt^T
    S <- crossprod(X)
    
    ## rescale so E[W] = target_df * Sigma, matching Wishart(target_df, Sigma)
    out[,,m] <- (target_df / T_len) * S
  }
  
  out
}