user_lib <- "/myRlibs"
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

repo_root <- normalizePath(".")

source(file.path(repo_root, "R", "01_postprocess.R"))
source(file.path(repo_root, "R", "02_numeric_utils.R"))
source(file.path(repo_root, "R", "03_collapsed_predictive.R"))
source(file.path(repo_root, "R", "04_collapsed_samplers.R"))
source(file.path(repo_root, "R", "05_distance_clustering.R"))
source(file.path(repo_root, "R", "06_moewishart.R"))

library(sna)
library(MCMCpack)
library(mclust)
library(cluster)
library(shapes)
library(sparsemixwishart)
library(moewishart)
library(future)
library(future.apply)

run_one <- function(idx, repo_root){
  set.seed(idx)
  
  ## ====================================
  ## 1) Create cluster-specific Sigma's
  ## ====================================
  p <- 3  # For small-matrix simulation, we set p = 3
  #p <- 6 # For medium-matrix simulation, we set p = 6
  
  make_ar1 <- function(rho, p){
    outer(1:p, 1:p, function(i,j) rho^abs(i-j))
  }
  
  ## Cluster 1: AR(1), rho=0.8
  Sigma1 <- make_ar1(0.8, p)
  
  ## Cluster 2: compound symmetry, rho=0.3
  rho2 <- 0.3
  Sigma2 <- (1 - rho2) * diag(p) + rho2 * matrix(1, p, p)
  
  ## Cluster 3: AR(1) with negative rho=-0.9
  Sigma3 <- make_ar1(-0.9, p)
  
  ## Cluster 4: identity
  Sigma4 <- make_ar1(0, p)
  
  ## Cluster 3: AR(1) with negative rho=-0.5
  Sigma5 <- make_ar1(-0.5, p)
  
  ## sanity check
  eig1 <- eigen(Sigma1, symmetric=TRUE, only.values=TRUE)$values
  eig2 <- eigen(Sigma2, symmetric=TRUE, only.values=TRUE)$values
  eig3 <- eigen(Sigma3, symmetric=TRUE, only.values=TRUE)$values
  eig4 <- eigen(Sigma4, symmetric=TRUE, only.values=TRUE)$values
  eig5 <- eigen(Sigma5, symmetric=TRUE, only.values=TRUE)$values
  print(round(eig1, 6)); print(round(eig2, 6)); print(round(eig3, 6)); print(round(eig4, 6)); print(round(eig5, 6))
  
  stopifnot(min(eig1) > 1e-10, min(eig2) > 1e-10, min(eig3) > 1e-10, min(eig4) > 1e-10, min(eig5) > 1e-10)
  
  ## =================================================
  ## 2) Generating data：heterogeneous Wishart samples
  ## =================================================
  n_vec <- c(10, 10, 10, 10, 10)     # for n=50,  balanced setting
  # n_vec <- c(20, 20, 20, 20, 20)   # for n=100, balanced setting
  # n_vec <- c(40, 40, 40, 40, 40)   # for n=200, balanced setting
  # n_vec <- c(5, 5, 10, 15, 15)     # for n=50,  unbalanced setting
  # n_vec <- c(10, 10, 20, 30, 30)   # for n=100, unbalanced setting
  # n_vec <- c(20, 20, 40, 60, 60)   # for n=200, unbalanced setting
  
  n <- sum(n_vec) # sample size
  K_true <- 5 # true number of clusters
  
  nu1_true <- 10
  nu2_true <- 15
  nu3_true <- 20
  nu4_true <- 25
  nu5_true <- 30
  
  W1 <- rWishart(n_vec[1], df = nu1_true, Sigma = Sigma1)  
  W2 <- rWishart(n_vec[2], df = nu2_true, Sigma = Sigma2)
  W3 <- rWishart(n_vec[3], df = nu3_true, Sigma = Sigma3)
  W4 <- rWishart(n_vec[4], df = nu4_true, Sigma = Sigma4)
  W5 <- rWishart(n_vec[5], df = nu5_true, Sigma = Sigma5)
  
  data <- array(NA_real_, dim=c(p,p,n))
  data[,, 1:n_vec[1]] <- W1
  data[,, (n_vec[1]+1):(n_vec[1]+n_vec[2])] <- W2
  data[,, (n_vec[1]+n_vec[2]+1):(n_vec[1]+n_vec[2]+n_vec[3])] <- W3
  data[,, (n_vec[1]+n_vec[2]+n_vec[3]+1):(n_vec[1]+n_vec[2]+n_vec[3]+n_vec[4])] <- W4
  data[,, (n_vec[1]+n_vec[2]+n_vec[3]+n_vec[4]+1):n] <- W5
  
  z_true <- c(rep(1, n_vec[1]), rep(2, n_vec[2]), rep(3, n_vec[3]), rep(4, n_vec[4]), rep(5, n_vec[5]))
  perm <- sample.int(n)
  data <- data[,,perm]
  z_true <- z_true[perm]
  
  ## =============================================
  ## 3) Running MFM-Wishart and baseline methods
  ## =============================================
  GAMMA <- 1 # MFM concentration parameter
  
  ## IW prior: Sigma ~ IW(Psi0, kappa0)
  ## Using kappa0=p+2 and Psi0=I => prior mean = I (since: mean = Psi0/(kappa0-p-1))
  kappa0 <- p + 2
  Psi0 <- diag(p)
  
  ## nu prior bounds
  nu_L <- p + 1 
  nu_U <- 50
  nu_init <- 20
  
  niterations = 10000
  
  ## ----- MFM-Wishart -----
  T_MFM <- system.time({
    MFM_fit <- MFM_Wishart_collapsed(
      data = data,
      niterations = niterations,
      GAMMA = GAMMA,
      Psi0 = Psi0,
      kappa0 = kappa0,
      nu_init = nu_init,
      nu_L = nu_L,
      nu_U = nu_U,
      initNClusters = n,
      proposal_sd = 1.0,
      lambda = 1,
      k_max = 500,
      m_aux = 10
    )
  })
  
  ## ----- DPM-Wishart -----
  T_DPM <- system.time({
    DPM_fit <- DPM_Wishart_collapsed(
      data = data,
      niterations = niterations,
      GAMMA = GAMMA,
      Psi0 = Psi0,
      kappa0 = kappa0,
      nu_init = nu_init,
      nu_L = nu_L,
      nu_U = nu_U,
      initNClusters = n,
      proposal_sd = 1,
      m_aux = 10
    )
  })
  
  ## ----- Finite Mixture of Wishart: moewishart + BIC selection -----
  SW_fit <- fit_moewishart_bic(
    data = data,
    K = 1:10,
    niter = 3000,
    estimate_nu = TRUE,
    n_restarts = 3,
    restart_iters = 20,
    tol = 1e-6,
    verbose = FALSE,
    z_true = z_true
  )
  
  ## ----- Penalized finite Mixture of Wishart (library(sparsemixwishart)) -----
  SW_pen_fit <- sparsemixwishart(data, K=c(1,2,3,4,5,6,7,8,9,10), penalty=c(0,2,4,6,8,10))
  
  ## ===========================================
  ## 4.1) Dahl's postprocessing for MFM-Wishart
  ## ===========================================
  burn <- 4000
  dahl <- getDahl(MFM_fit, burn = burn)
  
  z_hat  <- dahl[[1]]
  nu_hat <- dahl[[2]]
  
  K_MFM <- length(unique(z_hat))
  ARI_MFM <- adjustedRandIndex(z_true, z_hat)
  K_trace  <- sapply(MFM_fit$Iterates, function(x) length(unique(x$zout)))
  P_MFM <- mean(K_trace[(burn+1):niterations]==K_true)
  
  
  ## ----- Distance clustering (using the true number of clusters) -----
  ## ---- Riemannian distance matrix
  Dmat <- riem_dist_mat(data)
  
  #### ---- Hierarchical clustering
  z_hc <- hc_riemann_ward(Dmat, K = K_true)
  ARI_HC <- adjustedRandIndex(z_true, z_hc)
  
  #### ---- PAM (k-medoids)
  z_pam <- kmedoids_pam_riemann(Dmat, K = K_true)
  ARI_PAM <- adjustedRandIndex(z_true, z_pam)
  
  
  ## ===========================================
  ## 4.2) Dahl's postprocessing for DPM-Wishart
  ## ===========================================
  burn <- 4000
  dahl <- getDahl(DPM_fit, burn = burn)
  
  z_hat  <- dahl[[1]]
  nu_hat <- dahl[[2]]
  
  K_DPM <- length(unique(z_hat))
  ARI_DPM <- adjustedRandIndex(z_true, z_hat)
  K_trace  <- sapply(DPM_fit$Iterates, function(x) length(unique(x$zout)))
  P_DPM <- mean(K_trace[(burn+1):niterations]==K_true)
  
  ## ==========================================
  ## 4.3) Evaluate FMM (original + penalized)
  ## ==========================================
  K_SW <- SW_fit$K
  ARI_SW <- adjustedRandIndex(z_true, SW_fit$classification)
  
  K_SW_pen <- SW_pen_fit$K
  ARI_SW_pen <- adjustedRandIndex(z_true, SW_pen_fit$classification)
  
  res_list <- list(K_MFM=K_MFM,
                   K_DPM=K_DPM,
                   K_SW=K_SW,
                   K_SW_pen=K_SW_pen,
                   P_MFM=P_MFM,
                   P_DPM=P_DPM,
                   T_MFM=T_MFM,
                   T_DPM=T_DPM,
                   ARI_MFM=ARI_MFM,
                   ARI_DPM=ARI_DPM,
                   ARI_HC=ARI_HC,
                   ARI_SW=ARI_SW,
                   ARI_SW_pen=ARI_SW_pen,
                   ARI_PAM=ARI_PAM)
  
  out_dir <- file.path(repo_root, "results", "sim_small_5_50") # change the filename according to the setting
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_file <- file.path(out_dir, sprintf("small_5_%d.rds", idx)) # change the filename according to the setting
  saveRDS(res_list, out_file)
}


## ============================


plan(
  multisession,
  workers = 15
)

files <- future_lapply(
  1:100,
  run_one,
  repo_root = repo_root,
  future.seed = TRUE,
  future.packages = c(
    "sna",
    "MCMCpack",
    "mclust",
    "cluster",
    "shapes",
    "sparsemixwishart",
    "moewishart"
  ),
  future.stdout = FALSE
)

plan(sequential)