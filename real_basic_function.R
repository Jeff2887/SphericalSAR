haversine <- function(lon1, lat1, lon2, lat2, R = 6371) {
  # Convert to radians
  to_rad <- pi / 180
  lon1 <- lon1 * to_rad; lat1 <- lat1 * to_rad
  lon2 <- lon2 * to_rad; lat2 <- lat2 * to_rad
  
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  R * c  # distance in km (Earth radius ≈ 6371 km)
}



SAR_direct_sim_function <- function(seed_set,sample_size,neighbor_set, rho0=0.5){
  set.seed(seed_set)
  # Parameters
  n <- sample_size      # number of observations
  
  
  # Generate row-standardized weight matrix W:: sparse
  k <- neighbor_set         # number of non-zero neighbors per row
  
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  for (i in 1:n) {
    # Randomly select k columns (excluding the diagonal)
    neighbors <- sample(setdiff(1:n, i), k)
    
    # Assign random weights to these neighbors
    W[i, neighbors] <- runif(k)
    
    # Row-standardize
    W[i, ] <- W[i, ] / sum(W[i, ])
  }
  
  
  # Noise terms
  # Generate Gaussian noise
  epsilon <- rnorm(n, mean=0, sd=1)
  
  I_n <- diag(n)
  ##y <- solve(I_n - rho0 * W) %*% epsilon
  y <- rep(0,n)
  # Compute residuals
  Q_resid <- solve(I_n - rho0 * W)   # operator on residuals
  # Apply to each element in E
  for (i in 1:n) {
    yi <- 0
    for (j in 1:n) {
      yi <- yi + Q_resid[i,j] * epsilon[j]
    }
    y[i] <- yi
  }
  
  ## cost function
  objective <- function(rho, W, y) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- t(y) %*% t(S) %*% W %*% S %*% y
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W, y=y)$minimum
  
  
  ##hypothesis test
  sigma_est <- n * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  ##sigma_est <- n * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho0 * W))))^2
  wald_test <- n * rho_hat^2 / sigma_est
  
  
  return(c(rho_hat,wald_test))
}

## PSSAR
# Function to generate random skew-symmetric matrix
rSkewSym <- function(d) {
  M <- matrix(rnorm(d^2), d, d)
  M <- M - t(M)
  M <- M / sum(M^2)
  return(M) # skew-symmetric
}

# Function to generate a random zero-mean skew-symmetric matrix of size d x d
random_skew_symmetric_zero_mean <- function(d) {
  A <- matrix(rnorm(d^2,0,1), nrow=d)
  A <- A - t(A)  # Make skew-symmetric
  # Subtract the mean of all entries to make the "overall mean" zero
  A <- A - mean(A)
  A <- A / sum(A^2)
  return(A)
}

# Function to generate random skew-symmetric matrix
rSkewSym <- function(d) {
  m1 <- abs(rnorm(d,0,1))
  m2 <- rep(1,d)
  opt_trans <- log_rotation(m1,m2)
  M <- opt_trans$L
  return(M) 
}

# Function to generate a random skew-symmetric matrix of size d x d
random_skew_symmetric_zero_mean <- function(d) {
  m1 <- abs(rnorm(d,0,1))
  m2 <- rep(1,d)
  opt_trans <- log_rotation(m1,m2)
  M <- opt_trans$L
  return(M) 
}




PSSAR_direct_sim_function <- function(seed_set,sample_size,dim_set,neighbor_set, rho0=0.5,dim_red='PCA'){
  set.seed(seed_set)
  # Parameters
  n <- sample_size      # number of observations
  d <- dim_set          # dimension of skew-symmetric matrices
  
  # Generate row-standardized weight matrix W:: sparse
  k <- neighbor_set         # number of non-zero neighbors per row
  
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  
  for (i in 1:n) {
    # choose k neighbors for node i
    neighbors <- sample(setdiff(1:n, i), k)
    
    for (j in neighbors) {
      # assign a random weight
      w <- runif(1)
      W[i, j] <- w
      W[j, i] <- w   # enforce symmetry
    }
  }
  
  # Row-standardize
  row_sums <- rowSums(W)
  W <- W / row_sums
  
  # Mean matrix (skew-symmetric)
  q_bar <- random_skew_symmetric_zero_mean(d)
  
  # Noise terms
  epsilon <- lapply(1:n, function(i) random_skew_symmetric_zero_mean(d))
  epsilon_mean <- Reduce("+", epsilon) / n
  
  # Solve SAR-like system for q_i
  # Write in stacked form: Q - 1⊗q_bar = rho0 W (Q - 1⊗q_bar) + E
  # Equivalent: (I - rho0 W) (Q - 1⊗q_bar) = E
  I_n <- diag(n)
  
  # Stack epsilon as vector of skew-symmetric matrices
  E <- epsilon
  
  # Compute residuals
  Q_resid <- solve(I_n - rho0 * W)  # operator on residuals
  # Apply to each element in E
  Q <- vector("list", n)
  for (i in 1:n) {
    Qi <- matrix(0, d, d)
    for (j in 1:n) {
      Qi <- Qi + Q_resid[i,j] * (E[[j]] - epsilon_mean) ##E[[j]] ##
    }
    Q[[i]] <- q_bar + Qi
  }
  
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q) / n
  G_hat <- matrix(0, n, n)
  for (j in 1:n) {
    for (k in 1:n) {
      A <- Q[[j]] - q_hat ##q_bar
      B <- Q[[k]] - q_hat ##q_bar
      G_hat[j,k] <- sum(diag(t(A) %*% B))   # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
  
  
  ##hypothesis test
  
  if (dim_red=='PCA'){
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
    }
    # Step 1. Vectorize each epsilon_i
    U <- t(sapply(E_est_list, function(E) as.vector(E)))  # n x d^2
    
    # Step 2. Center columns
    Uc <- scale(U, center = TRUE, scale = FALSE)
    
    # Step 3. Truncated SVD (to avoid forming huge covariance)
    if (dim_set<=10){
      Sigma_hat <- cov(Uc)  # d^2 x d^2 matrix
      eig <- eigen(Sigma_hat, symmetric = TRUE)
      eigs <- eig$values
    } else {
      max_m = 50
      m <- min(max_m, n - 1)
      sv <- irlba(Uc, nv = m, nu = 0)
      eigs <- (sv$d^2) / n   # eigenvalues of sample covariance
    }
    
    
    # Step 4. Choose k by variance explained
    cumvar <- cumsum(eigs) / sum(eigs)
    k <- which(cumvar >= 0.9)[1]
    
    # Step 5. Compute varpi’s using top-k eigenvalues
    lambda <- eigs[1:k]
    
    # Step 5. Compute varpi_1, varpi_3, varpi_4
    varpi_1 <- (sum(lambda))^2
    varpi_2 <- sum(lambda^2)
  } else{
    varpi_1 <- 0
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
      cov_est <- sum(E_est*E_est) ##sum(diag(t(E_est) %*% E_est)) ##\|\epsilon_i\|^2
      varpi_1 <- varpi_1 + cov_est
    }
    varpi_1 <- (varpi_1/n)^2
    
    
    varpi_2 <- 0
    
    for (i in 1:n) {
      for (j in setdiff(1:n, i)) {       # j goes over 1:n excluding i
        ip <- sum(E_est_list[[i]] * E_est_list[[j]])  # Frobenius inner product
        varpi_2 <- varpi_2 + ip^2
      }
    }
    
    # divide by n*(n-1) for unbiased estimator
    varpi_2 <- varpi_2 / (n * (n - 1))
  }
  
  
  
  sigma_cst <- varpi_2 / varpi_1
  
  sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  ##sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W)))
  wald_test <- n * rho_hat^2 / sigma_est
  
  
  return(c(rho_hat,wald_test))
}


PSSAR_fit_function <- function(y,W,dim_red='PCA'){
  # Parameters
  n <- nrow(y)      # number of observations
  d <- ncol(y)          # dimension of skew-symmetric matrices
  dim_set <- ncol(y)  
  
  mu_intrinsic <- intrinsic_mean_sphere(y)
  Q <- list()
  for (i in 1:nrow(y)) {
    opt_tran <- log_rotation(mu_intrinsic,y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q) / n
  G_hat <- matrix(0, n, n)
  for (j in 1:n) {
    for (k in 1:n) {
      A <- Q[[j]] - q_hat ##q_bar
      B <- Q[[k]] - q_hat ##q_bar
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
  
  ##hypothesis test
  
  if (dim_red=='PCA'){
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
    }
    # Step 1. Vectorize each epsilon_i
    U <- t(sapply(E_est_list, function(E) as.vector(E)))  # n x d^2
    
    # Step 2. Center columns
    Uc <- scale(U, center = TRUE, scale = FALSE)
    
    # Step 3. Truncated SVD (to avoid forming huge covariance)
    if (dim_set<=10){
      Sigma_hat <- cov(Uc)  # d^2 x d^2 matrix
      eig <- eigen(Sigma_hat, symmetric = TRUE)
      eigs <- eig$values
    } else {
      max_m = 50
      m <- min(max_m, n - 1)
      sv <- irlba(Uc, nv = m, nu = 0)
      eigs <- (sv$d^2) / n   # eigenvalues of sample covariance
    }
    
    
    # Step 4. Choose k by variance explained
    cumvar <- cumsum(eigs) / sum(eigs)
    k <- which(cumvar >= 0.9)[1]
    
    # Step 5. Compute varpi’s using top-k eigenvalues
    lambda <- eigs[1:k]
    
    # Step 5. Compute varpi_1, varpi_3, varpi_4
    varpi_1 <- (sum(lambda))^2
    varpi_2 <- sum(lambda^2)
  } else{
    varpi_1 <- 0
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
      cov_est <- sum(E_est*E_est) ##sum(diag(t(E_est) %*% E_est)) ##\|\epsilon_i\|^2
      varpi_1 <- varpi_1 + cov_est
    }
    varpi_1 <- (varpi_1/n)^2
    
    
    varpi_2 <- 0
    
    for (i in 1:n) {
      for (j in setdiff(1:n, i)) {       # j goes over 1:n excluding i
        ip <- sum(E_est_list[[i]] * E_est_list[[j]])  # Frobenius inner product
        varpi_2 <- varpi_2 + ip^2
      }
    }
    
    # divide by n*(n-1) for unbiased estimator
    varpi_2 <- varpi_2 / (n * (n - 1))
  }
  
  
  
  sigma_cst <- varpi_2 / varpi_1
  
  sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(diag(n) - rho_hat * W))))^2
  ##sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W)))
  wald_test <- n * rho_hat^2 / sigma_est
  
  
  return(c(rho_hat,wald_test))
}

PSSAR_fit_function_fast <- function(y,W,dim_red='PCA'){
  # Parameters
  n <- nrow(y)      # number of observations
  d <- ncol(y)          # dimension of skew-symmetric matrices
  dim_set <- ncol(y)  
  I_n <- diag(n)
  
  mu_intrinsic <- intrinsic_mean_sphere(y)
  
  # 1. Compute all log rotations and store in 3D array
  Q <- array(0, dim = c(d, d, n))
  Q_LIST <- list()
  for (i in 1:n) {
    qsave <- log_rotation(mu_intrinsic, y[i,])$L
    Q[,,i] <- qsave
    Q_LIST[[i]] <- qsave
  }
  # 2. q_hat
  q_hat <- apply(Q, c(1,2), mean)
  
  # 3. Center matrices
  Qc <- sweep(Q, c(1,2), q_hat)
  
  # 4. Vectorize
  d2 <- d*d
  Qc_mat <- matrix(0, d2, n)
  for (i in 1:n) {
    Qc_mat[,i] <- as.vector(Qc[,,i])
  }
  
  # 5. Gram matrix (fast)
  G_hat <- crossprod(Qc_mat)
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
  
  
  ##hypothesis test
  
  if (dim_red=='PCA'){
    
    # Convert back to list of d × d matrices
    ##Q <- lapply(1:n, function(i) matrix(Q_mat[i, ], d, d))
    Q <- Q_LIST
    q_hat <- Reduce("+", Q) / n
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
    }
    
    
    # Step 1. Vectorize each epsilon_i
    U <- t(sapply(E_est_list, function(E) as.vector(E)))  # n x d^2
    
    # Step 2. Center columns
    Uc <- scale(U, center = TRUE, scale = FALSE)
    
    # Step 3. Truncated SVD (to avoid forming huge covariance)
    if (dim_set<=10){
      Sigma_hat <- cov(Uc)  # d^2 x d^2 matrix
      eig <- eigen(Sigma_hat, symmetric = TRUE)
      eigs <- eig$values
    } else {
      ##max_m = dim_set
      max_m = d
      m <- min(max_m, n - 1)
      sv <- irlba(Uc, nv = m, nu = 0)
      eigs <- (sv$d^2) / n   # eigenvalues of sample covariance
    }
    
    
    
    # Step 4. Choose k by variance explained
    cumvar <- cumsum(eigs) / sum(eigs)
    k <- which(cumvar >= 0.9)[1]
    
    
    # Step 5. Compute varpi’s using top-k eigenvalues
    lambda <- eigs[1:k]
    
    # Step 5. Compute varpi_1, varpi_3, varpi_4
    varpi_1 <- (sum(lambda))^2
    varpi_2 <- sum(lambda^2)
    sigma_cst <- varpi_2 / varpi_1
    sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  } else if (dim_red=='bootstrap'){
    Q_resid <- I_n
    
    ## faster
    ## ---- Precompute E_est_list ----
    # Q_centered_mat = Q_mat - q_hat_vec
    ##Q_centered_mat <- sweep(Q_mat, 2, q_hat_vec)
    Q_centered_mat <- t(Qc_mat)
    q_hat_vec <- as.vector(q_hat)
    # Compute E_est_mat under the null
    E_est_mat <- Q_centered_mat ##- rho_hat * (W %*% Q_centered_mat)
    
    # Convert back to list of d × d matrices
    E_est_list <- lapply(1:n, function(i) matrix(E_est_mat[i, ], d, d))
    
    
    
    rho_est_boots <- rep(0,500)
    for (boots_ind in 1:500) {
      print(boots_ind)
      
      ## faster 
      # 2. Bootstrap residual sample
      idx <- sample.int(n, n, TRUE)
      eps_boots_mat  <- E_est_mat[idx, , drop=FALSE]
      eps_boots_mean <- colMeans(eps_boots_mat)
      
      # 3. Compute Q_i
      eps_center_mat <- sweep(eps_boots_mat, 2, eps_boots_mean)
      Q_mat <- sweep(Q_resid %*% eps_center_mat, 2, q_hat_vec, "+")
      
      # 4. Compute G_hat
      q_hat_boots_vec <- colMeans(Q_mat)
      Q_centered_mat <- sweep(Q_mat, 2, q_hat_boots_vec)
      G_hat <- Q_centered_mat %*% t(Q_centered_mat)
      
      ## cost function
      objective <- function(rho, W, G_hat) {
        n <- nrow(W)
        S <- diag(n) - rho * W
        val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
        return(val^2)
      }
      
      rho_est_boots[boots_ind] <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
    }
    
    ##sigma_est <- var(rho_est_boots) * n
    ##plot(density(rho_est_boots))
    ci <- quantile(rho_est_boots, c(0.025, 0.975))
  } else if(dim_red=='No test'){
    wald_test <- NULL
  }else{
    Q <- Q_LIST
    q_hat <- Reduce("+", Q) / n
    
    varpi_1 <- 0
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
      cov_est <- sum(E_est*E_est) ##sum(diag(t(E_est) %*% E_est)) ##\|\epsilon_i\|^2
      varpi_1 <- varpi_1 + cov_est
    }
    varpi_1 <- (varpi_1/n)^2
    
    
    varpi_2 <- 0
    
    for (i in 1:n) {
      for (j in setdiff(1:n, i)) {       # j goes over 1:n excluding i
        ip <- sum(E_est_list[[i]] * E_est_list[[j]])  # Frobenius inner product
        varpi_2 <- varpi_2 + ip^2
      }
    }
    
    # divide by n*(n-1) for unbiased estimator
    varpi_2 <- varpi_2 / (n * (n - 1))
    sigma_cst <- varpi_2 / varpi_1
    sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  }
  

  if (dim_red=='bootstrap'){
    rej_result <- (rho_hat < ci[1] || rho_hat > ci[2])
  } else if(dim_red=='No test'){
    rej_result <- NULL
  }else{
    threshold_set <- qchisq(0.05, df=1, lower.tail=FALSE)
    wald_test <- n * rho_hat^2 / sigma_est
    rej_result <- (wald_test>=threshold_set)
  }
  
  return(c(rho_hat,rej_result))
}



## optimal transport
# -----------------------------
# Utilities for S^{d-1}
# -----------------------------
unit_norm <- function(x) x / sqrt(sum(x^2))

geo_dist <- function(x, y) {
  # Geodesic distance on the unit sphere (radians)
  z <- sum(x * y)
  z <- max(min(z, 1), -1)
  acos(z)
}

rand_on_sphere <- function(d) {
  z <- rnorm(d)
  unit_norm(z)
}

# Orthonormal basis of span{x,y}
span_basis_xy <- function(x, y) {
  u1 <- unit_norm(x)
  proj <- sum(y * u1) * u1
  v <- y - proj
  nv <- sqrt(sum(v^2))
  if (nv < 1e-12) {
    # y == +/- x: pick any orthonormal u2 orthogonal to u1
    tmp <- rnorm(length(x))
    tmp <- tmp - sum(tmp * u1) * u1
    u2 <- unit_norm(tmp)
  } else {
    u2 <- v / nv
  }
  list(u1 = u1, u2 = u2)
}



log_rotation <- function(x, y, tol = 1e-12) {
  x <- as.numeric(x); y <- as.numeric(y)
  stopifnot(length(x) == length(y))
  
  # ensure unit length
  x <- x / sqrt(sum(x^2))
  y <- y / sqrt(sum(y^2))
  
  cxy <- sum(x * y)
  cxy <- max(min(cxy, 1), -1)
  
  # identical or antipodal cases
  if (1 - abs(cxy) < tol) {
    if (cxy > 0) {
      return(list(L = matrix(0, length(x), length(x)), theta = 0, u1 = x, u2 = NA))
    } else {
      # antipodal: choose any unit u2 perpendicular to x
      j <- which.min(abs(x))
      e <- rep(0, length(x)); e[j] <- 1
      v <- e - sum(e * x) * x
      vn <- sqrt(sum(v^2))
      if (vn < tol) {
        if (length(x) < 2) stop("Need dimension >= 2 for antipodal case")
        v <- c(-x[2], x[1], if (length(x) > 2) rep(0, length(x) - 2) else NULL)
        vn <- sqrt(sum(v^2))
      }
      u2 <- v / vn
      theta <- pi
      ##Q <- tcrossprod(x, u2) - tcrossprod(u2, x)
      Q <- tcrossprod(u2, x) - tcrossprod(x, u2)
      return(list(L = theta * Q, theta = theta, u1 = x, u2 = u2))
    }
  }
  
  # general case
  u1 <- x
  v <- y - cxy * u1
  u2 <- v / sqrt(sum(v^2))
  theta <- acos(cxy)
  ##Q <- tcrossprod(u1, u2) - tcrossprod(u2, u1)  # skew-symmetric generator
  Q <- tcrossprod(u2, u1) - tcrossprod(u1, u2)
  
  list(L = theta * Q, theta = theta, u1 = u1, u2 = u2)
}


exp_from_log_rotation_predict <- function(q,theta) {
  #th <- attr(L, "theta")
  th <- theta
  L <- q
  if (is.null(th)) {
    # generic skew-symmetric: use series-safe expm via base eigen (skew-sym -> pure imaginary)
    eig <- eigen(L, symmetric = FALSE)
    V <- eig$vectors; D <- diag(exp(eig$values))
    Re(V %*% D %*% solve(V))
  } else {
    I <- diag(nrow(L))
    K <- L / max(th, 1e-16)
    I + sin(th) * K + (1 - cos(th)) * (K %*% K)
  }
}


exp_from_log_rotation <- function(L,a) {
  #th <- attr(L, "theta")
  th <- L$theta
  L <- L$L
  if (is.null(th)) {
    # generic skew-symmetric: use series-safe expm via base eigen (skew-sym -> pure imaginary)
    eig <- eigen(L, symmetric = FALSE)
    V <- eig$vectors; D <- diag(exp(eig$values))
    Re(V %*% D %*% solve(V))
  } else {
    I <- diag(nrow(L))
    K <- a * L / max(th, 1e-16)
    I + sin(th) * K + (1 - cos(th)) * (K %*% K)
  }
}


#### intrinsic mean
intrinsic_mean_sphere <- function(y, tol = 1e-8, max_iter = 1000) {
  n <- nrow(y)
  d <- ncol(y)
  
  # Initialize with extrinsic mean
  mu <- colSums(y)
  mu <- mu / sqrt(sum(mu^2))
  
  for (iter in 1:max_iter) {
    # Compute tangent vectors
    tangent_sum <- rep(0, d)
    for (i in 1:n) {
      theta <- acos(sum(mu * y[i,]))
      if (theta > 1e-12) {
        tangent <- (theta / sin(theta)) * (y[i,] - cos(theta) * mu)
        tangent_sum <- tangent_sum + tangent
      }
    }
    tangent_mean <- tangent_sum / n
    norm_t <- sqrt(sum(tangent_mean^2))
    
    # Check convergence
    if (norm_t < tol) break
    
    # Update mu using exponential map
    mu <- cos(norm_t) * mu + sin(norm_t) * (tangent_mean / norm_t)
  }
  
  return(mu)
}

##### function to generate z and normalize
PSSAR_trans_sim_function <- function(seed_set,sample_size, rho0=0.5){
  set.seed(seed_set)
  n <- sample_size
  d <- 10
  k <- 10 
  alpha <- rho0    # neighbor correlation strength (0 < alpha < 1)
  mu <- rep(1,d)   # shared mean vector
  mu <- mu / sqrt(sum(mu^2))  # normalize
  
  # ----------------------
  # Define neighbor adjacency matrix W
  # ----------------------
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  for (i in 1:n) {
    # Randomly select k columns (excluding the diagonal)
    neighbors <- sample(setdiff(1:n, i), k)
    
    # Assign random weights to these neighbors
    W[i, neighbors] <- runif(k)
    
    # Row-standardize
    W[i, ] <- W[i, ] / sum(W[i, ])
  }
  
  # ----------------------
  # Generate latent Gaussian z with neighbor correlation
  # ----------------------
  z <- matrix(0, n, d)
  for (k in 1:d) {
    eps <- rnorm(n)                     # independent noise
    z[,k] <- mu[k] + solve(diag(n) - alpha*W) %*% eps
  }
  
  # ----------------------
  # Normalize to unit sphere
  # ----------------------
  y <- t(apply(z, 1, function(v) v / sqrt(sum(v^2))))
  
  mu_intrinsic <- intrinsic_mean_sphere(y)
  
  Q <- list()
  for (i in 1:nrow(y)) {
    opt_tran <- log_rotation(mu_intrinsic,y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q) / n
  G_hat <- matrix(0, n, n)
  for (j in 1:n) {
    for (k in 1:n) {
      A <- Q[[j]] - q_hat ##q_bar
      B <- Q[[k]] - q_hat ##q_bar
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
  
  ##hypothesis test
  sigma_est <- n * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(diag(n) - rho_hat * W))))^2
  ##sigma_est <- n * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho0 * W))))^2
  wald_test <- n * rho_hat^2 / sigma_est
  
  
  return(c(rho_hat,wald_test))
  
}


##### function to generate z and normalize
PSSAR_trans_sim_prediction_function <- function(seed_set,sample_size, rho0=0.5){
  set.seed(seed_set)
  n <- sample_size
  d <- 10           # dimension of spherical vectors, previous setting 3
  k <- 10
  alpha <- rho0    # neighbor correlation strength (0 < alpha < 1)
  mu <- rep(1,d)   # shared mean vector
  mu <- mu / sqrt(sum(mu^2))  # normalize
  
  # ----------------------
  # Define neighbor adjacency matrix W
  # ----------------------
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  for (i in 1:n) {
    # Randomly select k columns (excluding the diagonal)
    neighbors <- sample(setdiff(1:n, i), k)
    
    # Assign random weights to these neighbors
    W[i, neighbors] <- runif(k)
    
    # Row-standardize
    W[i, ] <- W[i, ] / sum(W[i, ])
  }
  
  # ----------------------
  # Generate latent Gaussian z with neighbor correlation
  # ----------------------
  z <- matrix(0, n, d)
  for (k in 1:d) {
    eps <- rnorm(n)                     # independent noise
    z[,k] <- mu[k] + solve(diag(n) - alpha*W) %*% eps
  }
  
  ## consider use 1,...,n-1 to predict n
  n_train <- n-1
  z_train <- z[1:n_train,]
  W_train <- W[1:n_train,1:n_train]
  # ----------------------
  # Normalize to unit sphere
  # ----------------------
  y <- t(apply(z_train, 1, function(v) v / sqrt(sum(v^2))))
  
  mu_intrinsic <- intrinsic_mean_sphere(y)
  
  Q <- list()
  for (i in 1:nrow(y)) {
    opt_tran <- log_rotation(mu_intrinsic,y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q) / n_train
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- Q[[j]] - q_hat ##q_bar
      B <- Q[[k]] - q_hat ##q_bar
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(1:n_train, function(j) {
    W[n, j] * (Q[[j]] - q_hat)
  }))
  
  W_row <- W[n, ]        # n-th row
  w_nn <- W[n, n]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## project back
  angle_n <- sqrt(sum((q_n_pred %*% mu_intrinsic)^2))
  y_pred <- exp_from_log_rotation_predict(q_n_pred,angle_n) %*% mu_intrinsic
  
  z_n <- z[n,]
  y_true <- z_n / sqrt(sum(z_n^2))
  pred_dist <- as.numeric(acos(t(y_pred) %*% y_true))
  
  ##ignore spherical, use extrinsic way
  y_bar <- colMeans(y)
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- y[j,] - y_bar 
      B <- y[k,] - y_bar 
      G_hat[j,k] <- as.numeric(t(A) %*% B)  # Frobenius inner product
    }
  }
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(1:n_train, function(j) {
    W[n, j] * (y[j] - y_bar)
  }))
  
  W_row <- W[n, ]        # n-th row
  w_nn <- W[n, n]
  y_n_pred <- y_bar + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  y_n_pred <-  y_n_pred / sqrt(sum(y_n_pred^2))
  pred_dist_ignore <- as.numeric(acos(t(y_n_pred) %*% y_true))
  
  return(c(pred_dist,pred_dist_ignore))
}


predict_function <- function(y,W,ind,measure_set){
  ## consider use 1,...,n-1 to predict n
  n = nrow(y)
  n_train <- n-1
  ind_set <- setdiff(1:n,ind)
  y_train <- y[ind_set,]
  W_train <- W[ind_set,ind_set]
  
  ## SSAR
  mu_intrinsic <- intrinsic_mean_sphere(y_train)
  
  Q <- list()
  for (i in 1:n) {
    opt_tran <- log_rotation(mu_intrinsic,y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  Q_train <- Q[ind_set]
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q_train) / n_train
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- Q_train[[j]] - q_hat ##q_bar
      B <- Q_train[[k]] - q_hat ##q_bar
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(ind_set, function(j) {
    W[ind, j] * (Q[[j]] - q_hat)
  }))
  
  w_nn <- W[ind, ind]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## project back
  #angle_n <- sqrt(sum((q_n_pred %*% mu_intrinsic)^2))
  #y_pred1 <- exp_from_log_rotation_predict(q_n_pred,angle_n) %*% mu_intrinsic
  y_pred1 <- expm::expm(q_n_pred) %*% mu_intrinsic
  print(sum(y_pred1^2))
  
  y_true <- y[ind,]
  pred_dist <- as.numeric(acos(t(y_pred1) %*% y_true))
  
  y_true_density <- y_true^2
  y_pred1_density <- y_pred1^2
  
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD1 <- JSD(y_pred1_density,y_true_density)
  } else{
    pred_distJSD1 <- sum((y_pred1_density-y_true_density)^2) 
  }

  
  ##ignore spherical, use extrinsic way
  y_bar <- colMeans(y_train)
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- y_train[j,] - y_bar 
      B <- y_train[k,] - y_bar 
      G_hat[j,k] <- as.numeric(t(A) %*% B)  
    }
  }
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(ind_set, function(j) {
    W[ind, j] * (y[j,] - y_bar)
  }))
  
  w_nn <- W[ind, ind]
  y_n_pred2 <- y_bar + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  y_n_pred2 <-  y_n_pred2 / sqrt(sum(y_n_pred2^2))
  pred_dist_ignore <- as.numeric(acos(t(y_n_pred2) %*% y_true))
  
  y_pred2_density <- y_n_pred2^2
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD2 <- JSD(y_pred2_density,y_true_density)
  } else{
    pred_distJSD2 <- sum((y_pred2_density-y_true_density)^2) 
  }
  

  return(list(pred_dist,pred_dist_ignore,pred_distJSD1,pred_distJSD2,y_true_density,y_pred1_density,y_pred2_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}



KL <- function(p, q) {
  # Only indices where p > 0
  idx <- which(p > 0)
  sum(p[idx] * log(p[idx] / q[idx]))
}

JSD <- function(p, q) {
  m <- 0.5 * (p + q)
  0.5 * KL(p, m) + 0.5 * KL(q, m)
}

clr_trans <- function(mat) {
  if (is.null(dim(mat))) stop("mat must be a matrix with rows = compositions")
  gmean <- exp(rowMeans(log(mat)))
  
  # apply CLR
  clr <- log(mat / gmean)
  return(clr)
}

inv_clr <- function(clr_row) {
  x <- exp(clr_row)
  x / sum(x)
}

clr_predict_function <- function(y,W,ind,measure_set){
  y <- y^2
  ## deal with zero
  eps <- 1e-6
  y_pos <- y
  y_pos[y_pos == 0] <- eps
  y_pos <- y_pos / rowSums(y_pos)
  ## consider use 1,...,n-1 to predict n
  n = nrow(y_pos)
  n_train <- n-1
  ind_set <- setdiff(1:n,ind)
  y_train <- y_pos[ind_set,]
  W_train <- W[ind_set,ind_set]
  
  ## use CLR transfromation, ignore spherical, use extrinsic way
  clr_mat <- CLR(y_train)$LR ##clr_trans(y_train)
  clr_mat_whole <- CLR(y_pos)$LR ##clr_trans(y_train)
  
  y_bar <- colMeans(clr_mat)
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- clr_mat[j,] - y_bar 
      B <- clr_mat[k,] - y_bar 
      G_hat[j,k] <- as.numeric(t(A) %*% B)  
    }
  }
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(ind_set, function(j) {
    W[ind, j] * (clr_mat_whole[j,] - y_bar)
  }))
  
  w_nn <- W[ind, ind]
  y_n_pred2 <- y_bar + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  y_n_pred2 <- sqrt(inv_clr(y_n_pred2)) ## first back clr, then square root to sphere
  

  y_true_density <- y[ind,]
  y_true <- sqrt(y_true_density)
  
  pred_dist_ignore <- as.numeric(acos(t(y_n_pred2) %*% y_true)) ## angle 
  
  y_pred2_density <- y_n_pred2^2
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD2 <- JSD(y_pred2_density,y_true_density)
  } else{
    pred_distJSD2 <- sum((y_pred2_density-y_true_density)^2) 
  }
  
  
  return(list(pred_dist_ignore,pred_distJSD2,y_pred2_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}

########### kriging ########### 
krig_predict_function <- function(y,W,ind,measure_set){
  n = nrow(y)
  ind_set <- setdiff(1:n,ind)
  y_train <- y[ind_set,]

  
  
  neighbors_i <- which(W[ind, ] != 0)
  if (length(neighbors_i)==0){
    y_n_pred2 <- intrinsic_mean_sphere(y_train) ## if no neighour, use sample mean
  } else if (length(neighbors_i)==1){
    y_n_pred2 <- y[neighbors_i,]
  } else{
    y_neighbour <- y[neighbors_i,]
    y_n_pred2 <- intrinsic_mean_sphere(y_neighbour)
  }
  
  
  y_true <- y[ind,]
  y_true_density <- y_true^2
  
  pred_dist_ignore <- as.numeric(acos(t(y_n_pred2) %*% y_true)) ## angle 
  
  y_pred2_density <- y_n_pred2^2
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD2 <- JSD(y_pred2_density,y_true_density)
  } else{
    pred_distJSD2 <- sum((y_pred2_density-y_true_density)^2) 
  }
  
  
  return(list(pred_dist_ignore,pred_distJSD2,y_pred2_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}


########### SSAREXR ########### 

conditional_intrinsic_mean_sphere <- function(Y, w = NULL, max_iter = 100, tol = 1e-8) {
  # Y: n x d matrix, rows are y_i on the unit sphere
  # w: weights (length n), default = equal
  n <- nrow(Y)
  d <- ncol(Y)
  if (is.null(w)) w <- rep(1/n, n)
  w <- w / sum(w)
  
  # init: normalized extrinsic mean
  nu <- colSums(w * Y)
  nu <- nu / sqrt(sum(nu^2))
  
  for (iter in 1:max_iter) {
    # log map: tangent vectors
    theta <- acos(pmin(1, pmax(-1, Y %*% nu)))  # angles
    log_vecs <- matrix(0, n, d)
    
    for (i in 1:n) {
      if (theta[i] > 1e-12) {
        log_vecs[i, ] <- (theta[i] / sin(theta[i])) * 
          (Y[i, ] - cos(theta[i]) * nu)
      }
    }
    
    # weighted average in tangent space
    v <- colSums(w * log_vecs)
    norm_v <- sqrt(sum(v^2))
    if (norm_v < tol) break
    
    # exp map back to sphere
    nu <- cos(norm_v) * nu + sin(norm_v) * (v / norm_v)
    nu <- nu / sqrt(sum(nu^2))  # re-normalize
  }
  
  return(nu)
}


SSAREXR_predict_function <- function(y,covariate,W,ind,measure_set){
  ## consider use 1,...,n-1 to predict n
  n = nrow(y)
  n_train <- n-1
  ind_set <- setdiff(1:n,ind)
  y_train <- y[ind_set,]
  covariate_train <- covariate[ind_set,]
  W_train <- W[ind_set,ind_set]
  
  ## SSAREXR
  ## cost function
  objective <- function(rho, W, covariate, y) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    S_inv <- solve(S)
    
    cov_weight_matrix <- S_inv %*% covariate %*% solve(t(covariate) %*% covariate) %*% t(covariate) %*% S
    
    
    mean_matrix <- matrix(NA, nrow = n, ncol = ncol(y))
    for (i in 1:n) {
      weight_vec <- cov_weight_matrix[i,]
      current_mean <- (t(weight_vec) %*% y) / sum(weight_vec)
      mean_matrix[i,] <- current_mean / sqrt(sum(current_mean^2))
    }
    
    Q <- list()
    for (i in 1:nrow(y)) {
      opt_tran <- log_rotation(mean_matrix[i,],y[i,])
      Q[[i]] <- opt_tran$L
    }
    
    ## compute G_n
    # Compute Gram matrix G_hat
    q_hat <- Reduce("+", Q) / n
    G_hat <- matrix(0, n, n)
    for (j in 1:n) {
      for (k in 1:n) {
        A <- Q[[j]] - q_hat ##q_bar
        B <- Q[[k]] - q_hat ##q_bar
        G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
      }
    }
    
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, covariate=covariate_train, y=y_train)$minimum
  
  
  ##prediction
  S_train <- diag(n_train) - rho_hat * W_train
  S_inv_train <- solve(S_train)
  
  cov_weight_matrix <- S_inv_train %*% covariate_train %*% solve(t(covariate_train) %*% covariate_train) %*% t(covariate_train) %*% S_train
  
  
  mean_matrix <- matrix(NA, nrow = n_train, ncol = ncol(y))
  for (i in 1:n_train) {
    weight_vec <- cov_weight_matrix[i,]
    current_mean <- (t(weight_vec) %*% y_train) / sum(weight_vec)
    mean_matrix[i,] <- current_mean / sqrt(sum(current_mean^2))
  }
  
  Q <- list()
  for (i in 1:nrow(y_train)) {
    opt_tran <- log_rotation(mean_matrix[i,],y_train[i,])
    Q[[i]] <- opt_tran$L
  }
  
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q) / n
  # Weighted sum of deviations
  W_vec <- W[ind, -ind]
  
  sum_dev <- Reduce("+", lapply(1:nrow(y_train), function(j) {
    W_vec[j] * (Q[[j]] - q_hat)
  }))
  
  
  w_nn <- W[ind, ind]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## mean prediction
  S <- diag(n) - rho_hat * W
  S_inv <- solve(S)
  
  
  cov_weight_matrix <- S_inv %*% covariate %*% solve(t(covariate_train) %*% covariate_train) %*% t(covariate_train) %*% S_train
  weight_vec <- cov_weight_matrix[ind,]
  current_mean <- (t(weight_vec) %*% y_train) / sum(weight_vec)
  mu_pred <- as.numeric(current_mean / sqrt(sum(current_mean^2)))
  
  ## project back
  angle_n <- sqrt(sum((q_n_pred %*% mu_pred)^2))
  y_pred1 <- exp_from_log_rotation_predict(q_n_pred,angle_n) %*% mu_pred
  
  
  y_true <- y[ind,]
  pred_dist <- as.numeric(acos(t(y_pred1) %*% y_true))
  
  y_true_density <- y_true^2
  y_pred1_density <- y_pred1^2
  
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD1 <- JSD(y_pred1_density,y_true_density)
  } else{
    pred_distJSD1 <- sum((y_pred1-y_true)^2) 
  }
  
  
  
  return(list(pred_dist,pred_distJSD1,y_true_density,y_pred1_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}

intr_SSAREXR_predict_function <- function(y,covariate,W,ind,measure_set){
  ## consider use 1,...,n-1 to predict n
  n = nrow(y)
  n_train <- n-1
  ind_set <- setdiff(1:n,ind)
  y_train <- y[ind_set,]
  covariate_train <- covariate[ind_set,]
  W_train <- W[ind_set,ind_set]
  
  ## SSAREXR
  ## cost function
  objective <- function(rho, W, covariate, y) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    S_inv <- solve(S)
    
    cov_weight_matrix <- S_inv %*% covariate %*% solve(t(covariate) %*% covariate) %*% t(covariate) %*% S
    
    
    mean_matrix <- matrix(NA, nrow = n, ncol = ncol(y))
    for (i in 1:n) {
      weight_vec <- cov_weight_matrix[i,]
      current_mean <- conditional_intrinsic_mean_sphere(y,weight_vec)
      mean_matrix[i,] <- current_mean
    }
    
    Q <- list()
    for (i in 1:nrow(y)) {
      opt_tran <- log_rotation(mean_matrix[i,],y[i,])
      Q[[i]] <- opt_tran$L
    }
    
    ## compute G_n
    # Compute Gram matrix G_hat
    q_hat <- Reduce("+", Q) / n
    G_hat <- matrix(0, n, n)
    for (j in 1:n) {
      for (k in 1:n) {
        A <- Q[[j]] - q_hat ##q_bar
        B <- Q[[k]] - q_hat ##q_bar
        G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
      }
    }
    
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, covariate=covariate_train, y=y_train)$minimum
  
  
  ##prediction
  S_train <- diag(n_train) - rho_hat * W_train
  S_inv_train <- solve(S_train)
  
  cov_weight_matrix <- S_inv_train %*% covariate_train %*% solve(t(covariate_train) %*% covariate_train) %*% t(covariate_train) %*% S_train
  
  
  mean_matrix <- matrix(NA, nrow = n_train, ncol = ncol(y))
  for (i in 1:n_train) {
    weight_vec <- cov_weight_matrix[i,]
    current_mean <- conditional_intrinsic_mean_sphere(y_train,weight_vec)
    mean_matrix[i,] <- current_mean
  }
  
  Q <- list()
  for (i in 1:nrow(y_train)) {
    opt_tran <- log_rotation(mean_matrix[i,],y_train[i,])
    Q[[i]] <- opt_tran$L
  }
  
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q) / n
  # Weighted sum of deviations
  W_vec <- W[ind, -ind]
  
  sum_dev <- Reduce("+", lapply(1:nrow(y_train), function(j) {
    W_vec[j] * (Q[[j]] - q_hat)
  }))
  
  
  w_nn <- W[ind, ind]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## mean prediction
  S <- diag(n) - rho_hat * W
  S_inv <- solve(S)
  
  
  cov_weight_matrix <- S_inv %*% covariate %*% solve(t(covariate_train) %*% covariate_train) %*% t(covariate_train) %*% S_train
  weight_vec <- cov_weight_matrix[ind,]
  current_mean <- (t(weight_vec) %*% y_train) / sum(weight_vec)
  mu_pred <- as.numeric(current_mean / sqrt(sum(current_mean^2)))
  
  ## project back
  angle_n <- sqrt(sum((q_n_pred %*% mu_pred)^2))
  y_pred1 <- exp_from_log_rotation_predict(q_n_pred,angle_n) %*% mu_pred
  
  
  y_true <- y[ind,]
  pred_dist <- as.numeric(acos(t(y_pred1) %*% y_true))
  
  y_true_density <- y_true^2
  y_pred1_density <- y_pred1^2
  
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD1 <- JSD(y_pred1_density,y_true_density)
  } else{
    pred_distJSD1 <- sum((y_pred1-y_true)^2) 
  }
  
  
  
  return(list(pred_dist,pred_distJSD1,y_true_density,y_pred1_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}



#### Frechet regression ####
l2norm <- function(x){
  #sqrt(sum(x^2))
  as.numeric(sqrt(crossprod(x)))
}

SpheGeoDist <- function(y1,y2) {
  if (abs(length(y1) - length(y2)) > 0) {
    stop("y1 and y2 should be of the same length.")
  }
  if ( !isTRUE( all.equal(l2norm(y1),1) ) ) {
    stop("y1 is not a unit vector.")
  }
  if ( !isTRUE( all.equal(l2norm(y2),1) ) ) {
    stop("y2 is not a unit vector.")
  }
  y1 = y1 / l2norm(y1)
  y2 = y2 / l2norm(y2)
  if (sum(y1 * y2) > 1){
    return(0)
  } else if (sum(y1*y2) < -1){
    return(pi)
  } else return(acos(sum(y1 * y2)))
}


SpheGeoGrad <- function(x,y) { 
  tmp <- 1 - sum(x * y) ^ 2
  return(- (tmp) ^ (-0.5) * x)
  # if (tmp < tol) {
  #   return(- Inf * x)
  # } else {
  #   return(- (tmp) ^ (-0.5) * x)
  # }
}

SpheGeoHess <- function(x,y) { #,tol = 1e-10){
  return(- sum(x * y) * (1 - sum(x * y) ^ 2) ^ (-1.5) * x %*% t(x))
}

GloSpheReg <- function(xin=NULL, yin=NULL, xout=NULL){
  
  if (is.null(xout))
    xout <- xin
  
  
  yout <- GloSpheGeoReg(xin = xin, yin = yin, xout = xout)
  
  res <- list(xout = xout, yout = yout, xin = xin, yin = yin)
  class(res) <- "spheReg"
  return(res)
}




GloSpheGeoReg <- function(xin, yin, xout) {
  
  k = nrow(xout)
  n = nrow(xin)
  m = ncol(yin)
  
  xbar <- colMeans(xin)
  Sigma <- cov(xin) * (n-1) / n
  invSigma <- solve(Sigma)
  
  yout = sapply(1:k, function(j){
    s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xout[j,] - xbar)
    s <- as.vector(s)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if ( any( sapply( 1:n, function(i) isTRUE( all.equal( sum(yin[i,]*y0), 1 ) ) ) ) ){
      # if (sum(sapply(1:n, function(i) sum(yin[i,]*y0)) > 1-1e-8)){
      #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) + 
      #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
      y0[1] = y0[1] + 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument)
  })
  return(t(yout))
}


GloSpheRegNew <- function(xin=NULL, yin=NULL, xout=NULL, add_matrix){
  
  if (is.null(xout))
    xout <- xin
  
  yout <- GloSpheGeoRegNew(xin = xin, yin = yin, xout = xout, add_matrix)
  res <- list(xout = xout, yout = yout, xin = xin, yin = yin)
  class(res) <- "spheReg"
  return(res)
}


GloSpheGeoRegNew <- function(xin, yin, xout, add_matrix) {
  
  k = nrow(xout)
  n = nrow(xin)
  m = ncol(yin)
  
  xin <- cbind(1,xin)
  xout <- cbind(1,xout)
  
  yout = sapply(1:k, function(j){
    s <- t(xout[j,]) %*% solve(t(xin) %*% t(add_matrix) %*% add_matrix %*% xin) %*% t(xin) %*% t(add_matrix) %*% add_matrix

    s <- as.vector(s)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if ( any( sapply( 1:n, function(i) isTRUE( all.equal( sum(yin[i,]*y0), 1 ) ) ) ) ){
      # if (sum(sapply(1:n, function(i) sum(yin[i,]*y0)) > 1-1e-8)){
      #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) + 
      #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
      y0[1] = y0[1] + 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument)
  })
  return(t(yout))
}

######## new SSAREXR ######
new_SSAREXR_predict_function <- function(y,covariate,W,ind,measure_set){
  ## consider use 1,...,n-1 to predict n
  n = nrow(y)
  n_train <- n-1
  ind_set <- setdiff(1:n,ind)
  y_train <- y[ind_set,]
  
  if (dim(covariate)[2]==1){
    covariate_train <- matrix(covariate[ind_set,],ncol = 1)
  } else{
    covariate_train <- covariate[ind_set,]
  }
  W_train <- W[ind_set,ind_set]
  
  mean_est <- GloSpheReg(covariate_train,y_train,covariate)$yout
  
  ## SSAR
  
  Q <- list()
  for (i in 1:n) {
    opt_tran <- log_rotation(mean_est[i,],y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  Q_train <- Q[ind_set]
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q_train) / n_train
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- Q_train[[j]] - q_hat ##q_bar
      B <- Q_train[[k]] - q_hat ##q_bar
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  ##update beta using GLSE
  ##S_n_new <- diag(n_train) - rho_hat * W_train
  ##mean_est <- GloSpheRegNew(covariate_train,y_train,covariate,add_matrix=S_n_new)$yout
  ##Q <- list()
  ##for (i in 1:n) {
  ##  opt_tran <- log_rotation(mean_est[i,],y[i,])
  ##  Q[[i]] <- opt_tran$L
  ##}
  
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(ind_set, function(j) {
    W[ind, j] * (Q[[j]] - q_hat)
  }))
  
  w_nn <- W[ind, ind]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## project back
  #angle_n <- sqrt(sum((q_n_pred %*% mean_est[ind,])^2))
  #y_pred1 <- exp_from_log_rotation_predict(q_n_pred,angle_n) %*% mean_est[ind,]
  
  y_pred1 <- expm::expm(q_n_pred) %*% mean_est[ind,]
  
  print(sum(y_pred1^2))
  
  y_true <- y[ind,]
  pred_dist <- as.numeric(acos(t(y_pred1) %*% y_true))
  
  y_true_density <- y_true^2
  y_pred1_density <- y_pred1^2
  
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD1 <- JSD(y_pred1_density,y_true_density)
  } else{
    pred_distJSD1 <- sum((y_pred1-y_true)^2) 
  }
  
  
  
  return(list(pred_dist,pred_distJSD1,y_true_density,y_pred1_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}




###### block prediction #####
PSSAR_predict_block <- function(y,W,test_idx,measure_set){
  n <- nrow(y)
  train_idx <- setdiff(1:n, test_idx)
  
  # Optional: remove training nodes that are direct neighbors of test nodes
  ##adjacency <- (W != 0)
  ##neighbors <- apply(adjacency[train_idx, test_idx, drop=FALSE], 1, any)
  ##train_idx <- train_idx[!neighbors]
  n_train <- length(train_idx)
  n_test <- length(test_idx)
  
  # Submatrices
  W11 <- W[train_idx, train_idx]
  W21 <- W[test_idx, train_idx]
  W22 <- W[test_idx, test_idx]
  y_train <- y[train_idx,]
  y_test <- y[test_idx,]
  
  
  mu_intrinsic <- intrinsic_mean_sphere(y_train)
  
  Q <- list()
  for (i in 1:n) {
    opt_tran <- log_rotation(mu_intrinsic,y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  Q_train <- Q[train_idx]
  ## compute G_n
  # Compute Gram matrix G_hat
  q_hat <- Reduce("+", Q_train) / n_train
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- Q_train[[j]] - q_hat ##q_bar
      B <- Q_train[[k]] - q_hat ##q_bar
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W11, G_hat=G_hat)$minimum
  
  ##prediction
  weight_matrix <- rho_hat * solve(diag(n_test)-rho_hat * W22) %*% W21
  y_pred <- y_test
  if (n_test==1){
    Q_test <- Reduce("+", lapply(1:n_train, function(j) {
      weight_matrix[1, j] * (Q_train[[j]] - q_hat)
    }))
    Q_test <- Q_test + q_hat
    ## project back
    angle_n <- sqrt(sum((Q_test %*% mu_intrinsic)^2))
    y_pred1 <- exp_from_log_rotation_predict(Q_test,angle_n) %*% mu_intrinsic
    y_pred <- y_pred1
    
    angle_vec <- as.numeric(acos(t(y_pred) %*% y_test))
    y_true_density <- y_test^2
    y_pred1_density <- y_pred^2
    
    if (measure_set=='JSD'){
      ##JSD divergence
      mse_vec <- JSD(y_pred1_density,y_true_density)
    } else{
      mse_vec <- sum((y_pred1_density-y_true_density)^2) 
    }
  } else{
    for (i in 1:n_test) {
      Q_test <- Reduce("+", lapply(1:n_train, function(j) {
        weight_matrix[i, j] * (Q_train[[j]] - q_hat)
      }))
      Q_test <- Q_test + q_hat
      ## project back
      angle_n <- sqrt(sum((Q_test %*% mu_intrinsic)^2))
      y_pred1 <- exp_from_log_rotation_predict(Q_test,angle_n) %*% mu_intrinsic
      y_pred[i,] <- y_pred1
    }
    angle_vec <- rep(0,n_test)
    mse_vec <- rep(0,n_test)
    for (i in 1:n_test) {
      angle_vec[i] <- as.numeric(acos(t(y_pred[i,]) %*% y_test[i,]))
      y_true_density <- y_test[i,]^2
      y_pred1_density <- y_pred[i,]^2
      
      if (measure_set=='JSD'){
        ##JSD divergence
        mse_vec[i] <- JSD(y_pred1_density,y_true_density)
      } else{
        mse_vec[i] <- sum((y_pred1_density-y_true_density)^2) 
      }
    }
  }
  
  
  
  return(list(angle_vec,mse_vec))
}


Extrinsic_predict_block <- function(y,W,test_idx,measure_set){
  n <- nrow(y)
  train_idx <- setdiff(1:n, test_idx)
  
  # Optional: remove training nodes that are direct neighbors of test nodes
  ##adjacency <- (W != 0)
  ##neighbors <- apply(adjacency[train_idx, test_idx, drop=FALSE], 1, any)
  ##train_idx <- train_idx[!neighbors]
  n_train <- length(train_idx)
  n_test <- length(test_idx)
  
  # Submatrices
  W11 <- W[train_idx, train_idx]
  W21 <- W[test_idx, train_idx]
  W22 <- W[test_idx, test_idx]
  y_train <- y[train_idx,]
  y_test <- y[test_idx,]
  
  
  y_bar <- colMeans(y_train)
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- y_train[j,] - y_bar 
      B <- y_train[k,] - y_bar 
      G_hat[j,k] <- as.numeric(t(A) %*% B)  
    }
  }
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W11, G_hat=G_hat)$minimum
  
  
  ##prediction
  weight_matrix <- rho_hat * solve(diag(n_test)-rho_hat * W22) %*% W21
  y_pred <- y_test
  
  if (n_test==1){
    y_pred_i <- Reduce("+", lapply(1:n_train, function(j) {
      weight_matrix[1, j] * (y_train[j,] - y_bar)
    }))
    y_pred_i <- y_pred_i + y_bar
    ## project back
    
    y_pred <- y_pred_i / sqrt(sum(y_pred_i^2))
    
    angle_vec <- as.numeric(acos(t(y_pred) %*% y_test))
    y_true_density <- y_test^2
    y_pred1_density <- y_pred^2
    
    if (measure_set=='JSD'){
      ##JSD divergence
      mse_vec <- JSD(y_pred1_density,y_true_density)
    } else{
      mse_vec <- sum((y_pred1_density-y_true_density)^2) 
    }
  } else{
    for (i in 1:n_test) {
      y_pred_i <- Reduce("+", lapply(1:n_train, function(j) {
        weight_matrix[i, j] * (y_train[j,] - y_bar)
      }))
      y_pred_i <- y_pred_i + y_bar
      ## project back
      
      y_pred[i,] <- y_pred_i / sqrt(sum(y_pred_i^2))
    }
    
    
    angle_vec <- rep(0,n_test)
    mse_vec <- rep(0,n_test)
    for (i in 1:n_test) {
      angle_vec[i] <- as.numeric(acos(t(y_pred[i,]) %*% y_test[i,]))
      y_true_density <- y_test[i,]^2
      y_pred1_density <- y_pred[i,]^2
      
      if (measure_set=='JSD'){
        ##JSD divergence
        mse_vec[i] <- JSD(y_pred1_density,y_true_density)
      } else{
        mse_vec[i] <- sum((y_pred1_density-y_true_density)^2) 
      }
    }
  }
  
  return(list(angle_vec,mse_vec))
}


CLR_predict_block <- function(y,W,test_idx,measure_set){
  
  y <- y^2
  ## deal with zero
  eps <- 1e-6
  y_pos <- y
  y_pos[y_pos == 0] <- eps
  y_pos <- y_pos / rowSums(y_pos)
  
  n <- nrow(y)
  train_idx <- setdiff(1:n, test_idx)
  
  # Optional: remove training nodes that are direct neighbors of test nodes
  ##adjacency <- (W != 0)
  ##neighbors <- apply(adjacency[train_idx, test_idx, drop=FALSE], 1, any)
  ##train_idx <- train_idx[!neighbors]
  n_train <- length(train_idx)
  n_test <- length(test_idx)
  
  # Submatrices
  W11 <- W[train_idx, train_idx]
  W21 <- W[test_idx, train_idx]
  W22 <- W[test_idx, test_idx]
  y_train <- y_pos[train_idx,]
  y_test <- sqrt(y[test_idx,])
  
  ## use CLR transfromation, ignore spherical, use extrinsic way
  clr_mat <- CLR(y_train)$LR ##clr_trans(y_train)
  
  y_bar <- colMeans(clr_mat)
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- clr_mat[j,] - y_bar 
      B <- clr_mat[k,] - y_bar 
      G_hat[j,k] <- as.numeric(t(A) %*% B)  
    }
  }
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W11, G_hat=G_hat)$minimum
  
  
  ##prediction
  weight_matrix <- rho_hat * solve(diag(n_test)-rho_hat * W22) %*% W21
  y_pred <- y_test
  
  if (n_test==1){
    y_pred_i <- Reduce("+", lapply(1:n_train, function(j) {
      weight_matrix[1, j] * (y_train[j,] - y_bar)
    }))
    y_pred_i <- y_pred_i + y_bar
    ## project back
    
    y_pred <- sqrt(inv_clr(y_pred_i))
    
    angle_vec <- as.numeric(acos(t(y_pred) %*% y_test))
    y_true_density <- y_test^2
    y_pred1_density <- y_pred^2
    
    if (measure_set=='JSD'){
      ##JSD divergence
      mse_vec <- JSD(y_pred1_density,y_true_density)
    } else{
      mse_vec <- sum((y_pred1_density-y_true_density)^2) 
    }
  } else{
    for (i in 1:n_test) {
      y_pred_i <- Reduce("+", lapply(1:n_train, function(j) {
        weight_matrix[i, j] * (clr_mat[j,] - y_bar)
      }))
      y_pred_i <- y_pred_i + y_bar
      ## project back
      y_pred[i,] <- sqrt(inv_clr(y_pred_i)) ## first back clr, then square root to sphere
    }
    
    
    angle_vec <- rep(0,n_test)
    mse_vec <- rep(0,n_test)
    for (i in 1:n_test) {
      angle_vec[i] <- as.numeric(acos(t(y_pred[i,]) %*% y_test[i,]))
      y_true_density <- y_test[i,]^2
      y_pred1_density <- y_pred[i,]^2
      
      if (measure_set=='JSD'){
        ##JSD divergence
        mse_vec[i] <- JSD(y_pred1_density,y_true_density)
      } else{
        mse_vec[i] <- sum((y_pred1_density-y_true_density)^2) 
      }
    }
  }
  
  
  return(list(angle_vec,mse_vec))
}

##### Multivariate SAR model #####
clr_MSAR_predict_function <- function(y,W,ind,measure_set){
  y_square <- y^2
  ## deal with zero
  eps <- 1e-6
  y_pos <- y_square
  y_pos[y_pos == 0] <- eps
  y_pos <- y_pos / rowSums(y_pos)
  ## consider use 1,...,n-1 to predict n
  n = nrow(y_pos)
  n_train <- n-1
  ind_set <- setdiff(1:n,ind)
  y_train <- y_pos[ind_set,]
  W_train <- W[ind_set,ind_set]
  
  ## use CLR transfromation, ignore spherical, use extrinsic way
  clr_mat <- CLR(y_train)$LR ##clr_trans(y_train)
  clr_mat_whole <- CLR(y_pos)$LR ##clr_trans(y_train)
  
  
  y_n_pred2 <- rep(0,ncol(clr_mat))
  for (j in 1:ncol(clr_mat)) {
    y_uni <- as.vector(clr_mat[,j])
    y_uni_bar <- mean(y_uni)
    
    ## cost function
    objective <- function(rho, W, y, y_bar) {
      n <- nrow(W)
      S <- diag(n) - rho * W
      val <- t(y-y_bar) %*% t(S) %*% W %*% S %*% (y-y_bar)
      return(val^2)
    }
    
    rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, y=y_uni, y_bar=y_uni_bar)$minimum
    
    ##prediction
    # Weighted sum of deviations
    sum_dev <- Reduce("+", lapply(ind_set, function(k) {
      W[ind, k] * (clr_mat_whole[k,j] - y_uni_bar)
    }))
    
    w_nn <- W[ind, ind]
    y_n_pred2[j] <- y_uni_bar + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  }
  
  y_n_pred2 <- sqrt(inv_clr(y_n_pred2)) ## first back clr, then square root to sphere
  
  
  y_true_density <- y_square[ind,]
  y_true <- sqrt(y_true_density)
  
  pred_dist_ignore <- as.numeric(acos(t(y_n_pred2) %*% y_true)) ## angle 
  
  y_pred2_density <- y_n_pred2^2
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD2 <- JSD(y_pred2_density,y_true_density)
  } else{
    pred_distJSD2 <- sum((y_pred2_density-y_true_density)^2) 
  }
  
  
  return(list(pred_dist_ignore,pred_distJSD2,y_pred2_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}

clr_MSAR_fit_function <- function(y,W,ind){
  y_square <- y^2
  ## deal with zero
  eps <- 1e-6
  y_pos <- y_square
  y_pos[y_pos == 0] <- eps
  y_pos <- y_pos / rowSums(y_pos)
  
  ## use CLR transfromation, ignore spherical, use extrinsic way
  clr_mat <- CLR(y_pos)$LR ##clr_trans(y_train)
  
  
  p <- ncol(clr_mat)
  q <- 1 ##ncol(covariate)
  ### Initial value setup
  D0 = Matrix(0, p, p)
  B0 = Matrix(0, q, p)
  Sige0 = Diagonal(n = p, x = rep(0.4,p))
  ##Sige0[1,2] = 0.15
  ##Sige0[2,1] = 0.15
  
  Ymat <- clr_mat
  X <- rep(1,nrow(clr_mat))##covariate
  
  lse_res = MSAR.Lse(Ymat, X, W, D = D0, B = B0, Sige = Sige0) # LSE estimation
  
}


clr_MSAR_predict_function2 <- function(y,W,ind,measure_set){
  y_square <- y^2
  ## deal with zero
  eps <- 1e-6
  y_pos <- y_square
  y_pos[y_pos == 0] <- eps
  y_pos <- y_pos / rowSums(y_pos)
  ## consider use 1,...,n-1 to predict n
  n = nrow(y_pos)
  n_train <- n-1
  ind_set <- setdiff(1:n,ind)
  y_train <- y_pos[ind_set,]
  W_train <- W[ind_set,ind_set]
  
  ## use CLR transfromation, ignore spherical, use extrinsic way
  clr_mat <- CLR(y_train)$LR ##clr_trans(y_train)

  
  p <- ncol(clr_mat)
  q <- 1 ##ncol(covariate)
  ### Initial value setup
  D0 = Matrix(0, p, p)
  B0 = Matrix(0, q, p)
  Sige0 = Diagonal(n = p, x = rep(0.4,p))
  ##Sige0[1,2] = 0.15
  ##Sige0[2,1] = 0.15
  
  Ymat <- clr_mat
  X <- rep(1,nrow(clr_mat))##covariate
  
  lse_res = MSAR.Lse(Ymat, X, W_train, D = D0, B = B0, Sige = Sige0) # LSE estimation
  
  w_n_vec <- W[ind,ind_set]
  y_pred_MSAR <- as.vector(w_n_vec %*% clr_mat %*% lse_res$D + lse_res$B)
  
  y_n_pred2 <- sqrt(inv_clr(y_pred_MSAR)) ## first back clr, then square root to sphere
  
  
  y_true_density <- y_square[ind,]
  y_true <- sqrt(y_true_density)
  
  pred_dist_ignore <- as.numeric(acos(t(y_n_pred2) %*% y_true)) ## angle 
  
  y_pred2_density <- y_n_pred2^2
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD2 <- JSD(y_pred2_density,y_true_density)
  } else{
    pred_distJSD2 <- sum((y_pred2_density-y_true_density)^2) 
  }
  
  
  return(list(pred_dist_ignore,pred_distJSD2,y_pred2_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}


CLR_MSAR_predict_block <- function(y,W,test_idx,measure_set){
  
  y_square <- y^2
  ## deal with zero
  eps <- 1e-6
  y_pos <- y_square
  y_pos[y_pos == 0] <- eps
  y_pos <- y_pos / rowSums(y_pos)
  
  n <- nrow(y_square)
  train_idx <- setdiff(1:n, test_idx)
  
  # Optional: remove training nodes that are direct neighbors of test nodes
  ##adjacency <- (W != 0)
  ##neighbors <- apply(adjacency[train_idx, test_idx, drop=FALSE], 1, any)
  ##train_idx <- train_idx[!neighbors]
  n_train <- length(train_idx)
  n_test <- length(test_idx)
  
  # Submatrices
  W11 <- W[train_idx, train_idx]
  W21 <- W[test_idx, train_idx]
  W22 <- W[test_idx, test_idx]
  y_train <- y_pos[train_idx,]
  y_test <- sqrt(y_square[test_idx,])
  
  ## use CLR transfromation, ignore spherical, use extrinsic way
  clr_mat <- CLR(y_train)$LR ##clr_trans(y_train)
  
  y_pred_i <- NULL
  for (j in 1:ncol(clr_mat)) {
    y_uni <- as.vector(clr_mat[,j])
    y_uni_bar <- mean(y_uni)
    
    ## cost function
    objective <- function(rho, W, y, y_bar) {
      n <- nrow(W)
      S <- diag(n) - rho * W
      val <- t(y-y_bar) %*% t(S) %*% W %*% S %*% (y-y_bar)
      return(val^2)
    }
    
    rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W11, y=y_uni, y_bar=y_uni_bar)$minimum
    
    ##prediction
    # Weighted sum of deviations
    y_pred <- y_uni_bar + rho_hat * solve(diag(n_test)-rho_hat * W22) %*% W21 %*% (y_uni-y_uni_bar)
    
    y_pred_i <- cbind(y_pred_i,y_pred)
  }
  
  
  
  if (n_test==1){
    
    
    y_pred <- sqrt(inv_clr(y_pred_i))
    
    angle_vec <- as.numeric(acos(t(y_pred) %*% y_test))
    y_true_density <- y_test^2
    y_pred1_density <- y_pred^2
    
    if (measure_set=='JSD'){
      ##JSD divergence
      mse_vec <- JSD(y_pred1_density,y_true_density)
    } else{
      mse_vec <- sum((y_pred1_density-y_true_density)^2) 
    }
  } else{
    y_pred <- y_pred_i
    for (i in 1:n_test) {
      ## project back
      y_pred[i,] <- sqrt(inv_clr(y_pred_i[i,])) ## first back clr, then square root to sphere
    }
    
    
    angle_vec <- rep(0,n_test)
    mse_vec <- rep(0,n_test)
    for (i in 1:n_test) {
      angle_vec[i] <- as.numeric(acos(t(y_pred[i,]) %*% y_test[i,]))
      y_true_density <- y_test[i,]^2
      y_pred1_density <- y_pred[i,]^2
      
      if (measure_set=='JSD'){
        ##JSD divergence
        mse_vec[i] <- JSD(y_pred1_density,y_true_density)
      } else{
        mse_vec[i] <- sum((y_pred1_density-y_true_density)^2) 
      }
    }
  }
  
  
  return(list(angle_vec,mse_vec))
}


clr_MSAR_covariate_predict_function <- function(y,W,covariate,ind,measure_set){
  y_square <- y^2
  ## deal with zero
  eps <- 1e-6
  y_pos <- y_square
  y_pos[y_pos == 0] <- eps
  y_pos <- y_pos / rowSums(y_pos)
  ## consider use 1,...,n-1 to predict n
  n = nrow(y_pos)
  n_train <- n-1
  ind_set <- setdiff(1:n,ind)
  y_train <- y_pos[ind_set,]
  W_train <- W[ind_set,ind_set]
  covariate_train <- covariate[ind_set,]
  covariate_test <- covariate[ind,]
  ## use CLR transfromation, ignore spherical, use extrinsic way
  clr_mat <- CLR(y_train)$LR ##clr_trans(y_train)
  clr_mat_whole <- CLR(y_pos)$LR ##clr_trans(y_train)
  
  
  y_n_pred2 <- rep(0,ncol(clr_mat))
  for (j in 1:ncol(clr_mat)) {
    y_uni <- as.vector(clr_mat[,j])
    beta_est <- solve(t(covariate_train)%*%covariate_train)%*%t(covariate_train)%*%y_uni
    y_uni_bar <- mean(covariate_train %*% beta_est)
    
    ## cost function
    objective <- function(rho, W, y, y_bar) {
      n <- nrow(W)
      S <- diag(n) - rho * W
      val <- t(y-y_bar) %*% t(S) %*% W %*% S %*% (y-y_bar)
      return(val^2)
    }
    
    rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, y=y_uni, y_bar=y_uni_bar)$minimum
    
    ##prediction
    # Weighted sum of deviations
    sum_dev <- Reduce("+", lapply(ind_set, function(k) {
      W[ind, k] * (clr_mat_whole[k,j] - y_uni_bar)
    }))
    
    y_bar_pred <- as.numeric(covariate_test %*% beta_est)
    w_nn <- W[ind, ind]
    y_n_pred2[j] <- y_bar_pred + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  }
  
  y_n_pred2 <- sqrt(inv_clr(y_n_pred2)) ## first back clr, then square root to sphere
  
  
  y_true_density <- y_square[ind,]
  y_true <- sqrt(y_true_density)
  
  pred_dist_ignore <- as.numeric(acos(t(y_n_pred2) %*% y_true)) ## angle 
  
  y_pred2_density <- y_n_pred2^2
  if (measure_set=='JSD'){
    ##JSD divergence
    pred_distJSD2 <- JSD(y_pred2_density,y_true_density)
  } else{
    pred_distJSD2 <- sum((y_pred2_density-y_true_density)^2) 
  }
  
  
  return(list(pred_dist_ignore,pred_distJSD2,y_pred2_density))
  ##return(list(pred_distJSD1,pred_distJSD2,y_true,y_pred1,y_n_pred2))
}

#### SSARR


SRMSAR_fit_function<- function(y,W,X,dim_red='PCA'){
  # Parameters
  n <- nrow(y)      # number of observations
  d <- ncol(y)          # dimension of skew-symmetric matrices
  dim_set <- ncol(y)  
  I_n <- diag(n)
  
  mean_est <- GloSpheReg(X,y,X)$yout
  
  # 1. Compute all log rotations and store in 3D array
  Q <- array(0, dim = c(d, d, n))
  Q_LIST <- list()
  for (i in 1:n) {
    qsave <- log_rotation(mean_est[i,], y[i,])$L
    Q[,,i] <- qsave
    Q_LIST[[i]] <- qsave
  }
  # 2. q_hat
  q_hat <- apply(Q, c(1,2), mean)
  
  # 3. Center matrices
  Qc <- sweep(Q, c(1,2), q_hat)
  
  # 4. Vectorize
  d2 <- d*d
  Qc_mat <- matrix(0, d2, n)
  for (i in 1:n) {
    Qc_mat[,i] <- as.vector(Qc[,,i])
  }
  
  # 5. Gram matrix (fast)
  G_hat <- crossprod(Qc_mat)
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
  
  
  ##hypothesis test
  
  if (dim_red=='PCA'){
    
    # Convert back to list of d × d matrices
    ##Q <- lapply(1:n, function(i) matrix(Q_mat[i, ], d, d))
    Q <- Q_LIST
    q_hat <- Reduce("+", Q) / n
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
    }
    
    
    # Step 1. Vectorize each epsilon_i
    U <- t(sapply(E_est_list, function(E) as.vector(E)))  # n x d^2
    
    # Step 2. Center columns
    Uc <- scale(U, center = TRUE, scale = FALSE)
    
    # Step 3. Truncated SVD (to avoid forming huge covariance)
    if (dim_set<=10){
      Sigma_hat <- cov(Uc)  # d^2 x d^2 matrix
      eig <- eigen(Sigma_hat, symmetric = TRUE)
      eigs <- eig$values
    } else {
      ##max_m = dim_set
      max_m = d
      m <- min(max_m, n - 1)
      sv <- irlba(Uc, nv = m, nu = 0)
      eigs <- (sv$d^2) / n   # eigenvalues of sample covariance
    }
    
    
    
    # Step 4. Choose k by variance explained
    cumvar <- cumsum(eigs) / sum(eigs)
    k <- which(cumvar >= 0.9)[1]
    
    
    # Step 5. Compute varpi’s using top-k eigenvalues
    lambda <- eigs[1:k]
    
    # Step 5. Compute varpi_1, varpi_3, varpi_4
    varpi_1 <- (sum(lambda))^2
    varpi_2 <- sum(lambda^2)
    sigma_cst <- varpi_2 / varpi_1
    sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  } else if (dim_red=='bootstrap'){
    #E_est_list <- list()
    
    #for (i in 1:n) {
    #  E_est <- Q[[i]]-q_hat
    #  for (j in 1:n) {
    #    E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
    #  }
    #  E_est_list[[i]] <- E_est
    #}
    
    # Compute residuals
    ##Q_resid <- solve(I_n - rho_hat * W)  # operator on residuals
    Q_resid <- I_n
    
    ## faster
    ## ---- Precompute E_est_list ----
    # Q_centered_mat = Q_mat - q_hat_vec
    ##Q_centered_mat <- sweep(Q_mat, 2, q_hat_vec)
    Q_centered_mat <- t(Qc_mat)
    q_hat_vec <- as.vector(q_hat)
    # Compute E_est_mat under the null
    E_est_mat <- Q_centered_mat ##- rho_hat * (W %*% Q_centered_mat)
    
    # Convert back to list of d × d matrices
    E_est_list <- lapply(1:n, function(i) matrix(E_est_mat[i, ], d, d))
    
    
    
    rho_est_boots <- rep(0,500)
    for (boots_ind in 1:500) {
      print(boots_ind)
      
      ## faster 
      # 2. Bootstrap residual sample
      idx <- sample.int(n, n, TRUE)
      eps_boots_mat  <- E_est_mat[idx, , drop=FALSE]
      eps_boots_mean <- colMeans(eps_boots_mat)
      
      # 3. Compute Q_i
      eps_center_mat <- sweep(eps_boots_mat, 2, eps_boots_mean)
      Q_mat <- sweep(Q_resid %*% eps_center_mat, 2, q_hat_vec, "+")
      
      # 4. Compute G_hat
      q_hat_boots_vec <- colMeans(Q_mat)
      Q_centered_mat <- sweep(Q_mat, 2, q_hat_boots_vec)
      G_hat <- Q_centered_mat %*% t(Q_centered_mat)
      
      ## cost function
      objective <- function(rho, W, G_hat) {
        n <- nrow(W)
        S <- diag(n) - rho * W
        val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
        return(val^2)
      }
      
      rho_est_boots[boots_ind] <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
    }
    
    ##sigma_est <- var(rho_est_boots) * n
    ##plot(density(rho_est_boots))
    ci <- quantile(rho_est_boots, c(0.025, 0.975))
  } else if(dim_red=='No test'){
    wald_test <- NULL
  }else{
    Q <- Q_LIST
    q_hat <- Reduce("+", Q) / n
    
    varpi_1 <- 0
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
      cov_est <- sum(E_est*E_est) ##sum(diag(t(E_est) %*% E_est)) ##\|\epsilon_i\|^2
      varpi_1 <- varpi_1 + cov_est
    }
    varpi_1 <- (varpi_1/n)^2
    
    
    varpi_2 <- 0
    
    for (i in 1:n) {
      for (j in setdiff(1:n, i)) {       # j goes over 1:n excluding i
        ip <- sum(E_est_list[[i]] * E_est_list[[j]])  # Frobenius inner product
        varpi_2 <- varpi_2 + ip^2
      }
    }
    
    # divide by n*(n-1) for unbiased estimator
    varpi_2 <- varpi_2 / (n * (n - 1))
    sigma_cst <- varpi_2 / varpi_1
    sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  }
  
  
  if (dim_red=='bootstrap'){
    rej_result <- (rho_hat < ci[1] || rho_hat > ci[2])
  } else if(dim_red=='No test'){
    rej_result <- NULL
  }else{
    threshold_set <- qchisq(0.05, df=1, lower.tail=FALSE)
    wald_test <- n * rho_hat^2 / sigma_est
    rej_result <- (wald_test>=threshold_set)
  }
  
  return(c(rho_hat,rej_result))
}


conformal_PI_function <- function(y,W,ind,dist_matrix,alpha = 0.1){
  ## consider use 1,...,n-1 to predict n
  n = nrow(y)
  d = ncol(y)
  ind_set <- setdiff(1:n,ind)
  
  idx <- sample(ind_set, (n-1))
  train <- idx[1:round((n-1)/2)]
  calib <- idx[(round((n-1)/2) + 1):(n-1)]
  
  
  dist_numeric <- units::drop_units(dist_matrix)
  calib_distances <- dist_numeric[n, calib]
  eta <- median(dist_numeric[dist_numeric > 0])
  w_raw_initial <- exp(-(calib_distances^2)/(2 * eta^2))
  w_raw <- w_raw_initial
  sum_w <- sum(w_raw) + 1
  w_tilde <- w_raw / sum_w
  w_tilde_n <- 1 / sum_w
  
  # 3. Fit SAR on training
  W_train <- W[train, train]
  y_train <- y[train,]
  mu_intrinsic <- intrinsic_mean_sphere(y_train)
  
  Q <- list()
  for (k in 1:n) {
    opt_tran <- log_rotation(mu_intrinsic,y[k,])
    Q[[k]] <- opt_tran$L
  }
  
  Q_train <- Q[train]
  train_len <- length(train)
  
  q_hat <- Reduce("+", Q_train) / train_len
  G_hat <- matrix(0, train_len, train_len)
  for (j in 1:train_len) {
    for (k in 1:train_len) {
      A <- Q_train[[j]] - q_hat 
      B <- Q_train[[k]] - q_hat 
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  ## cost function
  objective <- function(rho, W, G_hat) {
    S <- diag(nrow(W)) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  # 4. Compute calibration residuals using full W 
  R_vec <- rep(0,(n-1))
  res_weight <- diag(n-1) - rho_hat * W[ind_set,ind_set]
  row_ind <- 1
  for (i in ind_set) {
    Ri <- matrix(0, d, d)
    col_ind <- 1
    for (j in ind_set) {
      Ri <- Ri + res_weight[row_ind,col_ind] * (Q[[j]]-q_hat) 
      col_ind <- col_ind + 1
    }
    R_vec[row_ind] <- sqrt(sum(Ri^2))
    row_ind <- row_ind + 1
  }
  
  R_vec_full <- rep(0,n)
  R_vec_full[ind_set] <- R_vec
  scores <- R_vec_full[calib]
  
  ord <- order(scores)
  sort_scores <- scores[ord]
  sort_w <- w_tilde[ord]
  
  vals <- c(sort_scores, Inf)
  cum_probs <- cumsum(c(sort_w, w_tilde_n))
  
  q_alpha <- vals[which(cum_probs >= (1 - alpha))[1]]
  
  
  # For a new residual R_new, check coverage
  y_true <- y[ind,]
  R_new <- log_rotation(mu_intrinsic,y_true)
  
  sum_dev <- Reduce("+", lapply(train, function(j) {
    W[n, j] * (Q[[j]] - q_hat)
  }))
  
  q_n_pred <- q_hat + rho_hat * sum_dev 
  
  predict_score <- sqrt(sum((R_new$L-q_n_pred)^2))
  
  inside <- (predict_score <= q_alpha)
  
  
  return(inside)
}
