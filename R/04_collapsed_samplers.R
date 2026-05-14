# The code adapted from the code for Yin et al. (2023)
# See MFM_MN.R in https://github.com/fyin-stats/MFM-MxN
# Fan Yin, Guanyu Hu & Weining Shen (2023) Analysis of Professional
# Basketball Field Goal Attempts via a Bayesian Matrix Clustering Approach, Journal of
# Computational and Graphical Statistics, 32:1, 49-60, DOI: 10.1080/10618600.2022.2085727

MFM_Wishart_collapsed <- function(
    data, niterations,
    GAMMA,            # Dirichlet gamma
    Psi0, kappa0,     # IW prior hyperparams
    nu_init, nu_L, nu_U,  # nu prior bounds + init
    initNClusters,    # Number of initial clusters
    proposal_sd = 1,  # RW-MH sd for nu
    lambda = 1,       # Poisson^+(lambda) for K
    k_max = 500       # truncation for VN
){
  
  ## data: p x p x n array, each slice SPD
  n <- dim(data)[3]
  p <- nrow(data[,,1])
  
  ## Precompute logdet(Wi)(for updating z and nu)
  logdetW <- rep(NA_real_, n)
  for(i in 1:n){
    logdetW[i] <- logdet_spd(data[,,i])
  }
  
  ## Precompute VN in the MFM prior (k truncated at a large integer)
  N <- n
  VN <- 0
  tmax <- min(n + 10, k_max)  
  VN <- rep(NA_real_, tmax)
  
  for (t in 1:tmax){
    r <- log(0)
    for (k in t:k_max){
      b <- sum(log((k-t+1):k)) -
        sum(log((k*GAMMA):(k*GAMMA+N-1))) +
        dpois(k-1, lambda, log = TRUE)
      m <- max(b, r)
      r <- log(exp(r-m) + exp(b-m)) + m
    }
    VN[t] <- r
  }
  
  ## Initialize clusterAssign
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  ## for collapsed sampling
  ## For fast computation, we keep the sufficient stats in each cluster
  ##                S_c = sum W_i,  A_c = sum log|W_i|
  clusterSizes <- table(as.factor(clusterAssign))
  nClusters <- length(clusterSizes)
  
  Sarr <- array(0, dim=c(p,p,nClusters))
  A    <- rep(0, nClusters)
  for(i in 1:n){
    c <- clusterAssign[i]
    Sarr[,,c] <- Sarr[,,c] + data[,,i]
    A[c] <- A[c] + logdetW[i]
  }
  
  nu <- nu_init
  
  History <- vector("list", niterations)
  
  ## start MCMC
  for(iter in 1:niterations){
    
    #if(iter %% 1000 == 0){
    #  print(iter)
    #}
    
    ## update z, the labels for each sample
    clusterSizes <- table(as.factor(clusterAssign))
    nClusters <- length(clusterSizes)
    
    for(i in 1:n){
      
      cur.cluster.i <- clusterAssign[i]
      Wi <- data[,,i]
      logdetWi <- logdetW[i]
      
      if(clusterSizes[cur.cluster.i] > 1){
        
        ## i is a non-singleton: |C_-i| = nClusters
        c.counts.noi <- clusterSizes
        c.counts.noi[cur.cluster.i] <- c.counts.noi[cur.cluster.i] - 1
        
        ### existing clusters
        log_clusterProbs <- sapply(1:nClusters, function(x){
          
          ni <- as.numeric(c.counts.noi[x])
          if(ni <= 0) return(-Inf)
          
          ## Si = S_{x,-i}
          Si <- Sarr[,,x]
          if(x == cur.cluster.i){
            Si <- Si - Wi
          }
          
          log_prior <- log(GAMMA + ni)
          log_like  <- log_pred_existing(Wi, Si, ni, nu, Psi0, kappa0, logdetWi)
          log_prior + log_like
        })
        
        ### new cluster
        log_prior_new <- log(GAMMA) + (VN[nClusters+1] - VN[nClusters])
        log_like_new  <- log_pred_new(Wi, nu, Psi0, kappa0, logdetWi)
        log_clusterProbs[nClusters+1] <- log_prior_new + log_like_new
        
        ### sample
        cluster.i <- sample(1:(nClusters+1), size=1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs)))
        
        ### Update sufficient stats 
        if(cluster.i == cur.cluster.i){
          #### stay: do nothing
        } else if(cluster.i <= nClusters){
          #### move to existing
          Sarr[,,cur.cluster.i] <- Sarr[,,cur.cluster.i] - Wi
          A[cur.cluster.i] <- A[cur.cluster.i] - logdetWi
          Sarr[,,cluster.i] <- Sarr[,,cluster.i] + Wi
          A[cluster.i] <- A[cluster.i] + logdetWi
          clusterAssign[i] <- cluster.i
        } else {
          #### open new cluster
          Sarr[,,cur.cluster.i] <- Sarr[,,cur.cluster.i] - Wi
          A[cur.cluster.i] <- A[cur.cluster.i] - logdetWi
          
          Snew <- array(0, dim=c(p,p,nClusters+1))
          Snew[,,1:nClusters] <- Sarr
          Snew[,,nClusters+1] <- Wi
          Sarr <- Snew
          
          A <- c(A, logdetWi)
          
          clusterAssign[i] <- nClusters + 1
        }
        
        clusterSizes <- table(as.factor(clusterAssign))
        nClusters <- length(clusterSizes)
        
      } else {
        
        ## i is a  singleton: |C_-i| = nClusters - 1
        ### existing clusters
        log_clusterProbs <- sapply(1:nClusters, function(x){
          
          if(x == cur.cluster.i) return(-Inf)
          
          ni <- as.numeric(clusterSizes[x])  # n_{x,-i}=n_x
          Si <- Sarr[,,x]
          
          log_prior <- log(GAMMA + ni)
          log_like  <- log_pred_existing(Wi, Si, ni, nu, Psi0, kappa0, logdetWi)
          log_prior + log_like
        })
        
        ### new cluster
        log_prior_new <- log(GAMMA) + (VN[nClusters] - VN[nClusters-1])
        
        log_like_new  <- log_pred_new(Wi, nu, Psi0, kappa0, logdetWi)
        log_clusterProbs[nClusters+1] <- log_prior_new + log_like_new
        
        cluster.i <- sample(1:(nClusters+1), size=1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs)))
        
        if(cluster.i > nClusters){
          #### stay alone: keep label cur.cluster.i, no stats change
          clusterAssign[i] <- cur.cluster.i
        } else {
          #### move i into existing cluster, then delete empty cur cluster and relabel
          #### update target stats first (indices are still old)
          Sarr[,,cluster.i] <- Sarr[,,cluster.i] + Wi
          A[cluster.i] <- A[cluster.i] + logdetWi
          
          ## drop cur.cluster.i slice from Sarr/A
          Sarr <- Sarr[,,-cur.cluster.i, drop=FALSE]
          A    <- A[-cur.cluster.i]
          
          ## relabel clusterAssign
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign)
        }
        
        clusterSizes <- table(as.factor(clusterAssign))
        nClusters <- length(clusterSizes)
      }
    } ## end i loop
    
    ## ---- update nu by RW-MH (fully collapsed)
    ## recompute nvec from current z (consistent with Sarr/A)
    nvec <- as.numeric(table(factor(clusterAssign, levels=1:max(clusterAssign))))
    nvec <- nvec[nvec > 0]  # occupied clusters
    
    nu_prop <- nu + rnorm(1, mean=0, sd=proposal_sd)
    
    if(nu_prop >= nu_L && nu_prop <= nu_U){
      
      log_old <- log_target_nu(nu,      Sarr, A, nvec, Psi0, kappa0)
      log_new <- log_target_nu(nu_prop, Sarr, A, nvec, Psi0, kappa0)
      
      log_alpha <- log_new - log_old
      if(log(runif(1)) < log_alpha){
        nu <- nu_prop
      }
    }## else: reject automatically
    
    ## store
    History[[iter]] <- list(zout = clusterAssign, nuout = nu)
    
  } # end iter loop
  
  list(Iterates = History, VN = VN)
}



#########################################
DPM_Wishart_collapsed <- function(
    data, niterations,
    GAMMA,            # DP concentration (i.e., alpha)
    Psi0, kappa0,     # IW prior hyperparams
    nu_init, nu_L, nu_U,  # nu prior bounds + init
    initNClusters,   # number of initial clusters
    proposal_sd = 1  # RW-MH sd for nu
){
  
  ## data: p x p x n array, each slice SPD
  n <- dim(data)[3]
  p <- nrow(data[,,1])
  
  ## Precompute logdet(Wi)(for updating z and nu)
  logdetW <- rep(NA_real_, n)
  for(i in 1:n){
    logdetW[i] <- logdet_spd(data[,,i])
  }
  
  ## initialize clusterAssign
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  ## for collapsed sampling
  ## For fast computation, we keep the sufficient stats in each cluster
  ##                S_c = sum W_i,  A_c = sum log|W_i|
  clusterSizes <- table(as.factor(clusterAssign))
  nClusters <- length(clusterSizes)
  
  Sarr <- array(0, dim=c(p,p,nClusters))
  A    <- rep(0, nClusters)
  for(i in 1:n){
    c <- clusterAssign[i]
    Sarr[,,c] <- Sarr[,,c] + data[,,i]
    A[c] <- A[c] + logdetW[i]
  }
  
  nu <- nu_init
  
  History <- vector("list", niterations)
  
  ## start MCMC
  for(iter in 1:niterations){
    
    #if(iter %% 1000 == 0){
    #  print(iter)
    #}
    
    ## update z
    clusterSizes <- table(as.factor(clusterAssign))
    nClusters <- length(clusterSizes)
    
    for(i in 1:n){
      
      cur.cluster.i <- clusterAssign[i]
      Wi <- data[,,i]
      logdetWi <- logdetW[i]
      
      if(clusterSizes[cur.cluster.i] > 1){
        
        ## i is a non-singleton: |C_-i| = nClusters
        c.counts.noi <- clusterSizes
        c.counts.noi[cur.cluster.i] <- c.counts.noi[cur.cluster.i] - 1
        
        ### existing clusters
        log_clusterProbs <- sapply(1:nClusters, function(x){
          
          ni <- as.numeric(c.counts.noi[x])
          if(ni <= 0) return(-Inf)
          
          ### Si = S_{x,-i}
          Si <- Sarr[,,x]
          if(x == cur.cluster.i){
            Si <- Si - Wi
          }
          
          log_prior <- log(ni)
          log_like  <- log_pred_existing(Wi, Si, ni, nu, Psi0, kappa0, logdetWi)
          log_prior + log_like
        })
        
        ### new cluster
        log_prior_new <- log(GAMMA)
        log_like_new  <- log_pred_new(Wi, nu, Psi0, kappa0, logdetWi)
        log_clusterProbs[nClusters+1] <- log_prior_new + log_like_new
        
        ## sample
        cluster.i <- sample(1:(nClusters+1), size=1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs)))
        
        ## update sufficient stats
        if(cluster.i == cur.cluster.i){
          ### stay: do nothing
        } else if(cluster.i <= nClusters){
          ### move to existing
          Sarr[,,cur.cluster.i] <- Sarr[,,cur.cluster.i] - Wi
          A[cur.cluster.i] <- A[cur.cluster.i] - logdetWi
          Sarr[,,cluster.i] <- Sarr[,,cluster.i] + Wi
          A[cluster.i] <- A[cluster.i] + logdetWi
          clusterAssign[i] <- cluster.i
        } else {
          ### open new cluster
          Sarr[,,cur.cluster.i] <- Sarr[,,cur.cluster.i] - Wi
          A[cur.cluster.i] <- A[cur.cluster.i] - logdetWi
          
          Snew <- array(0, dim=c(p,p,nClusters+1))
          Snew[,,1:nClusters] <- Sarr
          Snew[,,nClusters+1] <- Wi
          Sarr <- Snew
          
          A <- c(A, logdetWi)
          
          clusterAssign[i] <- nClusters + 1
        }
        
        clusterSizes <- table(as.factor(clusterAssign))
        nClusters <- length(clusterSizes)
        
      } else {
        
        ## i is a  singleton: |C_-i| = nClusters - 1
        ### existing clusters
        log_clusterProbs <- sapply(1:nClusters, function(x){
          
          if(x == cur.cluster.i) return(-Inf)
          
          ni <- as.numeric(clusterSizes[x])  # n_{x,-i}=n_x
          Si <- Sarr[,,x]
          
          log_prior <- log(ni)
          log_like  <- log_pred_existing(Wi, Si, ni, nu, Psi0, kappa0, logdetWi)
          log_prior + log_like
        })
        
        
        ### new cluster
        log_prior_new <- log(GAMMA)
        
        log_like_new  <- log_pred_new(Wi, nu, Psi0, kappa0, logdetWi)
        log_clusterProbs[nClusters+1] <- log_prior_new + log_like_new
        
        cluster.i <- sample(1:(nClusters+1), size=1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs)))
        
        if(cluster.i > nClusters){
          ## stay alone: keep label cur.cluster.i, no stats change
          clusterAssign[i] <- cur.cluster.i
        } else {
          ## move i into existing cluster, then delete empty cur cluster and relabel
          ## update target stats first (indices are still old)
          Sarr[,,cluster.i] <- Sarr[,,cluster.i] + Wi
          A[cluster.i] <- A[cluster.i] + logdetWi
          
          ## drop cur.cluster.i slice from Sarr/A
          Sarr <- Sarr[,,-cur.cluster.i, drop=FALSE]
          A    <- A[-cur.cluster.i]
          
          ## relabel clusterAssign 
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign)
        }
        
        clusterSizes <- table(as.factor(clusterAssign))
        nClusters <- length(clusterSizes)
      }
    } # end i loop
    
    ## ---- update nu by RW-MH (fully collapsed)
    ## recompute nvec from current z (consistent with Sarr/A)
    nvec <- as.numeric(table(factor(clusterAssign, levels=1:max(clusterAssign))))
    nvec <- nvec[nvec > 0]  # occupied clusters
    
    nu_prop <- nu + rnorm(1, mean=0, sd=proposal_sd)
    
    if(nu_prop >= nu_L && nu_prop <= nu_U){
      
      log_old <- log_target_nu(nu,      Sarr, A, nvec, Psi0, kappa0)
      log_new <- log_target_nu(nu_prop, Sarr, A, nvec, Psi0, kappa0)
      
      log_alpha <- log_new - log_old
      if(log(runif(1)) < log_alpha){
        nu <- nu_prop
      }
    }## else: reject automatically
    
    ## store
    History[[iter]] <- list(zout = clusterAssign, nuout = nu)
    
  } # end iter loop
  
  list(Iterates = History)
}