log_pred_existing <- function(Wi, Si, ni, nu, Psi0, kappa0, logdetWi){
  ## ni = n_{c,-i}, Si = S_{c,-i}
  p <- nrow(Wi)
  logdet_Psi_S   <- logdet_spd(Psi0 + Si)
  logdet_Psi_S_W <- logdet_spd(Psi0 + Si + Wi)
  
  out <- ( log_multigamma((kappa0 + (ni+1)*nu)/2, p)
           - log_multigamma((kappa0 + ni*nu)/2, p)
           - log_multigamma(nu/2, p)
           + ((nu - p - 1)/2) * logdetWi
           + ((kappa0 + ni*nu)/2) * logdet_Psi_S
           - ((kappa0 + (ni+1)*nu)/2) * logdet_Psi_S_W )
  out
}

log_pred_new <- function(Wi, nu, Psi0, kappa0, logdetWi){
  ## m(Wi|nu) = prior predictive (new cluster)
  p <- nrow(Wi)
  out <- ( log_multigamma((nu + kappa0)/2, p)
           - log_multigamma(nu/2, p)
           - log_multigamma(kappa0/2, p)
           + ((nu - p - 1)/2) * logdetWi
           + (kappa0/2) * logdet_spd(Psi0)
           - ((nu + kappa0)/2) * logdet_spd(Psi0 + Wi) )
  out
}



log_target_nu <- function(nu, Sarr, A, nvec, Psi0, kappa0){ 
  #fully-collapsed nu: log target
  p <- dim(Sarr)[1]
  K <- length(nvec)
  out <- 0
  for(c in 1:K){
    nc <- nvec[c]
    Sc <- Sarr[,,c]
    out <- out +
      log_multigamma((kappa0 + nc*nu)/2, p) -
      nc * log_multigamma(nu/2, p) +
      ((nu - p - 1)/2) * A[c] -
      ((kappa0 + nc*nu)/2) * logdet_spd(Psi0 + Sc)
  }
  out
}
