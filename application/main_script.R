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

## ====================================
## 1) Import the data
## ====================================
dir_csv <- "./data/processed_data"   # put the preprocessed data here
files <- list.files(dir_csv, pattern = "\\.csv$", full.names = TRUE)

get_id <- function(f) as.integer(sub("\\.csv$", "", basename(f)))
files <- files[order(sapply(files, get_id))]

length(files)  # should be 99

## ====================================
## 2) Load channels 40:46 
## ====================================
idx <- 40:46
p_sub <- length(idx)

data <- array(NA_real_, dim = c(p_sub, p_sub, length(files)))

for (i in seq_along(files)) {
  R <- as.matrix(read.csv(files[i], header = FALSE))
  
  if (!all(dim(R) >= max(idx))) {
    stop(sprintf("File %s has dim %s, cannot slice 40:46",
                 basename(files[i]), paste(dim(R), collapse="x")))
  }
  
  data[,,i] <- R[idx, idx]
}

dim(data)  # 7 7 99


## ====================================
## 3) Running MCMC for MFM--Wishart
## ====================================
library(sna)      
library(MCMCpack) 
library(mclust) 
library(cluster)


p <- 7
nu0 <- 55

kappa0 <- p + 5     # 12
Psi0   <- (kappa0 - p - 1) * (1/nu0) * diag(p)   # (4/53)*I
nu_init <- nu0
nu_L <- 10
nu_U <- 100

GAMMA <- 1

niterations = 10000

set.seed(123)

## MFM time
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
    initNClusters = 99,
    proposal_sd = 4.0,
    lambda = 1,
    k_max = 500
  )
})

## ===========================================
## 4) Dahl's postprocessing for MFM-Wishart
## ===========================================
burn <- 4000
dahl <- getDahl(MFM_fit, burn = burn)

z_hat  <- dahl[[1]]
nu_hat <- dahl[[2]]

K_MFM <- length(unique(z_hat))
K_MFM
z_hat
#P_MFM <- mean(K_trace[(burn+1):niterations]==3)

## ===========================================
## 5) Summary and Plots
## ===========================================
stopifnot(length(z_hat) == dim(data)[3])

clusters <- sort(unique(z_hat))
p <- dim(data)[1]

cluster_mean <- lapply(clusters, function(k){
  idx <- which(z_hat == k)
  Mk <- apply(data[,,idx, drop = FALSE], c(1,2), mean, na.rm = TRUE)
  symm(Mk)
})
names(cluster_mean) <- paste0("Cluster_", clusters)


table(z_hat)

library(fields)

chan_labels <- 40:46
stopifnot(length(chan_labels) == p)

col_full <- fields::tim.colors(64)
col_pos  <- col_full[(floor(length(col_full) / 2) + 1):length(col_full)]

plot_corrmat2 <- function(M, main = "", zlim = c(0, 1),
                          show_numbers = TRUE, labels = NULL){
  M <- symm(M)
  p <- nrow(M)
  if(is.null(labels)) labels <- 1:p
  
  Z <- t(M[p:1, ])
  
  image(
    1:p, 1:p, Z,
    zlim = zlim,
    col = col_pos,
    xaxt = "n", yaxt = "n",
    xlab = "", ylab = "",
    main = main
  )
  
  axis(1, at = 1:p, labels = labels, tick = FALSE)
  axis(2, at = 1:p, labels = rev(labels), tick = FALSE)
  box()
  
  if(show_numbers){
    for(i in 1:p){
      for(j in 1:p){
        text(i, p-j+1, sprintf("%.2f", M[i,j]), cex = 0.8)
      }
    }
  }
}


op <- par(no.readonly = TRUE)

layout(
  matrix(c(1, 0, 2, 3), nrow = 1),
  widths = c(1, 0.05, 1, 0.18)
)

for(k in clusters){
  Mk <- cluster_mean[[paste0("Cluster_", k)]]
  nk <- sum(z_hat == k)
  
  par(mar = c(4, 4, 3, 1))
  
  plot_corrmat2(
    Mk,
    main = sprintf("Cluster %d (n=%d)", k, nk),
    zlim = c(0, 1),
    labels = chan_labels
  )
}

par(mar = c(4, 0, 3, 4))
plot.new()

fields::image.plot(
  legend.only = TRUE,
  zlim = c(0, 1),
  legend.width = 5, 
  col = col_pos,
  axis.args = list(
    at = seq(0, 1, by = 0.2),
    labels = sprintf("%.1f", seq(0, 1, by = 0.2))
  )
)

layout(1)
par(op)



## ===========================================
## 6) Post hoc analysis: sex
## ===========================================

# https://osf.io/7fzkm/files/n6dsr
# sex: 1 = Female, 2 = Male
gender <- c(1, 1, 1, 2, 1, 1, 2, 2, 2,
            1, 2, 1, 1, 2, 2, 2, 1, 1,
            1, 1, 2, 1, 1, 1, 2, 2, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 1, 1, 1, 2, 1, 2, 1, 2,
            2, 1, 2, 2, 1, 1, 1, 2, 1,
            1, 2, 1, 2, 1, 2, 2, 2, 2,
            2, 2, 2, 2, 1, 1, 2, 1, 2,
            2, 1, 2, 1, 1, 1, 1, 2, 1,
            1, 2, 2, 1, 2, 1, 2, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 2)
adjustedRandIndex(z_hat, gender)

gender_factor <- factor(gender, levels = c(1, 2),
                        labels = c("Female", "Male"))

cluster_factor <- factor(z_hat)

tab_from_data <- table(cluster_factor, gender_factor)
tab_from_data

fisher.test(tab_from_data, alternative = "two.sided")$p.value
