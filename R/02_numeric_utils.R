logdet_spd <- function(A){ 
  R <- chol(A)
  2 * sum(log(diag(R)))
}

log_multigamma <- function(a, p){ 
  # log Γ_p(a) = p(p-1)/4 log(pi) + sum_{j=1}^p log Γ(a-(j-1)/2)
  out <- (p*(p-1)/4) * log(pi)
  for(j in 1:p){
    out <- out + lgamma(a - (j-1)/2)
  }
  out
}

symm <- function(A) (A + t(A))/2