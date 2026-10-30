dist_cov_creator <- function(data, method = "Euclidean") {
  # Compute Non-Euclidean distances for covariance matrices, according to Dryden AOAS 2009
  # Useful for subsequently performing hierarchical clustering with different distance types
  N <- dim(data)[3]
  dist_matrix <- matrix(NA, nrow = N, ncol = N)
  
  for (row in 1:N) {
    for (col in row:N) {
      dist_row_col <- shapes::distcov(
        S1 = data[, , row],
        S2 = data[, , col],
        method = method
      )
      dist_matrix[row, col] <- dist_row_col
      dist_matrix[col, row] <- dist_row_col
    }
  }
  stats::as.dist(dist_matrix)
}


riem_dist_mat <- function(data){ # construct n x n distance matrix
  n <- dim(data)[3]
  D <- matrix(0, n, n)
  for(i in 1:(n-1)){
    Ai <- data[,,i]
    for(j in (i+1):n){
      d <- distcov(Ai, data[,,j]) # need library(shapes)
      D[i,j] <- d
      D[j,i] <- d
    }
  }
  D
}

hc_riemann_ward <- function(Dmat, K){
  ## Dmat: n x n distance matrix
  D <- as.dist(Dmat)
  hc <- hclust(D, method = "ward.D")  
  cutree(hc, k = K)
}


kmedoids_pam_riemann <- function(Dmat, K){
  ## need library(cluster)
  D <- as.dist(Dmat)
  pam_fit <- cluster::pam(D, k = K, diss = TRUE)
  pam_fit$clustering
}


