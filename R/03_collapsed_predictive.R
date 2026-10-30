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
  ## m(Wi|nu) = prior predictive (new cluster), with Sigma integrated out
  p <- nrow(Wi)
  out <- ( log_multigamma((nu + kappa0)/2, p)
           - log_multigamma(nu/2, p)
           - log_multigamma(kappa0/2, p)
           + ((nu - p - 1)/2) * logdetWi
           + (kappa0/2) * logdet_spd(Psi0)
           - ((nu + kappa0)/2) * logdet_spd(Psi0 + Wi) )
  out
}



## The collapsed log target for one occupied cluster's nu_c.
log_target_nu_cluster <- function(nu, Sc, Ac, nc, Psi0, kappa0){
  p <- nrow(Sc)
  out <- ( log_multigamma((kappa0 + nc*nu)/2, p)
           - nc * log_multigamma(nu/2, p)
           + ((nu - p - 1)/2) * Ac
           - ((kappa0 + nc*nu)/2) * logdet_spd(Psi0 + Sc) )
  out
}