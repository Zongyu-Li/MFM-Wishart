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
source(file.path(repo_root, "R", "06_sparsemixwishart_helpers.R"))
source(file.path(repo_root, "R", "07_sparsemixwishart.R"))

library(sna)
library(MCMCpack)
library(mclust)
library(cluster)
library(shapes)

idx = 1

set.seed(idx)

## ====================================
## 1) Create cluster-specific Sigma's
## ====================================
p <- 12

make_ar1 <- function(rho, p) outer(1:p, 1:p, function(i,j) rho^abs(i-j))

make_block_cs <- function(block_sizes, rho_in=0.5, rho_out=0.1){ # Generating Sigma with sparse structure; fixed acorss replicates
  p <- sum(block_sizes)
  S <- matrix(rho_out, p, p); diag(S) <- 1
  idx_start <- cumsum(c(1, head(block_sizes, -1)))
  for(b in seq_along(block_sizes)){
    idx <- idx_start[b]:(idx_start[b] + block_sizes[b] - 1)
    S[idx, idx] <- rho_in
    diag(S[idx, idx]) <- 1
  }
  S
}

make_random_corr <- function(p, df=p+10, seed=NULL){ # random correlation matrix; changing across replicates
  if(!is.null(seed)) set.seed(seed)
  A <- matrix(rnorm(df*p), df, p)
  cov2cor(crossprod(A))
}

Sigma0 <- make_ar1(0.2, p) 

Sigma1 <- make_block_cs(c(2,3,7), rho_in=0.5, rho_out=0)
Sigma2 <- make_block_cs(c(3,4,5), rho_in=0.5, rho_out=0)
Sigma3 <- make_random_corr(p, df=p+12, seed=123+idx)

## Sanity check
all.equal(diag(Sigma1), rep(1, p))
all.equal(diag(Sigma2), rep(1, p))
all.equal(diag(Sigma3), rep(1, p))
eig1 <- eigen(Sigma1, symmetric=TRUE, only.values=TRUE)$values
eig2 <- eigen(Sigma2, symmetric=TRUE, only.values=TRUE)$values
eig3 <- eigen(Sigma3, symmetric=TRUE, only.values=TRUE)$values
print(round(eig1, 6)); print(round(eig2, 6)); print(round(eig3, 6))

stopifnot(min(eig1) > 1e-10, min(eig2) > 1e-10, min(eig3) > 1e-10)

## =================================================
## 2) Generating data：heterogeneous Wishart samples
## =================================================
n_vec <- c(16, 17, 17)      # n = 50
# n_vec <- c(33, 33, 34)    # n = 100
# n_vec <- c(66, 67, 67)    # n = 200

n <- sum(n_vec)  #sample size
K_true <- 3 # true number of clusters

nu_true <- 15 # shared degrees of freedom

W1 <- rWishart(n_vec[1], df = nu_true, Sigma = Sigma1)  
W2 <- rWishart(n_vec[2], df = nu_true, Sigma = Sigma2)
W3 <- rWishart(n_vec[3], df = nu_true, Sigma = Sigma3)

data <- array(NA_real_, dim=c(p,p,n))
data[,, 1:n_vec[1]] <- W1
data[,, (n_vec[1]+1):(n_vec[1]+n_vec[2])] <- W2
data[,, (n_vec[1]+n_vec[2]+1):(n_vec[1]+n_vec[2]+n_vec[3])] <- W3

z_true <- c(rep(1, n_vec[1]), rep(2, n_vec[2]), rep(3, n_vec[3]))
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
    k_max = 500
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
    proposal_sd = 1.0
  )
})

## ----- Penalized finite Mixture of Wishart -----
SW_pen_fit <- sparsemixwishart(data, K=c(1,2,3,4,5), penalty=c(0,2,4,6,8,10))

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
P_MFM <- mean(K_trace[(burn+1):niterations]==3)


## ----- Distance clustering (using the estimated number of clusters from MFM-Wishart) -----
## ---- Riemannian distance matrix
Dmat <- riem_dist_mat(data)

#### ---- Hierarchical clustering
z_hc <- hc_riemann_ward(Dmat, K = K_MFM)
ARI_HC <- adjustedRandIndex(z_true, z_hc)

#### ---- PAM (k-medoids)
z_pam <- kmedoids_pam_riemann(Dmat, K = K_MFM)
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
P_DPM <- mean(K_trace[(burn+1):niterations]==3)


## ==========================================
## 4.3) Evaluate FMM (original + penalized)
## ==========================================
K_SW_pen <- SW_pen_fit$K
ARI_SW_pen <- adjustedRandIndex(z_true, SW_pen_fit$classification)


res_list <- list(K_MFM=K_MFM,
                 K_DPM=K_DPM,
                 K_SW_pen=K_SW_pen,
                 P_MFM=P_MFM,
                 P_DPM=P_DPM,
                 T_MFM=T_MFM,
                 T_DPM=T_DPM,
                 ARI_MFM=ARI_MFM,
                 ARI_DPM=ARI_DPM,
                 ARI_SW_pen=ARI_SW_pen,
                 ARI_HC=ARI_HC,
                 ARI_PAM=ARI_PAM)

out_dir <- file.path(repo_root, "results", "sim_large_3_50") # change the filename according to the setting
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, sprintf("large_3_%d.rds", idx)) # change the filename according to the setting
saveRDS(res_list, out_file)