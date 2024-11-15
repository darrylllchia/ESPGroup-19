# Jiayi Wu - s2664441, Ruobing Dai - s2655029, Darryl Chia - s2740198
# Darryl Chia writed the LMMsetup function, Jiayi Wu completed LMMprof function and Ruobing Dai write lmm function.
# Then we checked the LMMprof function carafully and modify it for efficiency. 
# All of us explained our parts to each other to understand and improve our work. 
# Finally, we test our function and agreed on this final piece of work. 

# Aim of this task: Write a function to estimate simple linear mixed model by optimizing log likelihood.

# Approach:
# 1. Use cholesky decomposition and QR decomposition to calculate beta_hat efficiently.
# 2. Write loglikelihood function as objective function of optimization.
# 3. Use Nelder-Mead method to optimize objective function with respect to theta and estimate beta_hat and theta.


LMMsetup <- function(form, dat, ref){
  #' The function is to set up our fixed effect and random effect matrices X and Z, along with the y vector. 
  #' The vector f contains the number of columns of each block in each term in ref. 
  
  # Convert the input data into a data frame format
  df <- data.frame(dat) 
  
  # Set up X (Fixed effects)
  X = model.matrix(form, df)
  
  # Set up y (Response variable)
  y <- df[, names(df) == as.character(form[[2]])]
  
  # Set up Z (Random effects)
  Z <- matrix(nrow = dim(df)[1], ncol = 0)
  
  # Initilise vector f to store the number of columns of each block
  f <- c()
  
  # Update random effect matrix Z if we have random effect
  if (length(ref) > 0){
    for (i in 1:length(ref)){
      
      # Convert random effect into formula format
      fomula <- paste('~', paste(paste(ref[[i]], collapse = ':'), '-1', sep = ''), sep = '')
      # Set up model matrix
      model <- model.matrix(as.formula(fomula), df)
      
      # Combine model matrix of each block together
      Z <- cbind(Z, model)
      # Number of columns of each block
      f <- c(f, ncol(model))
    }
  }
  
  return(list(X = X, Z = Z, y = y, f = f))
  
}



LMMprof <- function(theta,X,Z,y,f){
  #'
  #' Evaluate the (negative) log likelihood for a given theta
  #' Compute the corresponding beta_hat
  #'
  #' @param theta: vector of log(std_devs). First element of theta is the log(std_dev) of residuals.
  #'               Subsequent elements correspond to the log(std_devs) of each block in Z.
  #' @return negative log likelihood function for optim to minimise and the calculated beta_hat vector given a theta vector.
  
  # Number of columns in Z matrix
  p <- sum(f)
  # Number of observations
  n <- dim(X)[1]
  # Number of fixed effects
  p1 <- dim(X)[2]
  
  # Set up sigma
  sigma <- exp(theta[1])
  
  # If there is no random effect
  if (length(f) == 0) {
    # Compute beta_hat and negative log likelihood function
    qrx <- qr(X)
    beta_hat <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p1])
    neg_loglike <- t(y - X %*% beta_hat) %*% (y - X %*% beta_hat) /
      (2*sigma^2) + n*log(sigma)
  }
  
  # If there is random effect
  else {
    # theta = c(log(sigma), ...) each block have the same standard deviation
    # Extend theta to diagonal of psi
    ex_theta <- rep(exp(theta[-1])^2, f)
    
    # QR decompose Z
    qrz <- qr(Z)
    R <- qr.R(qrz)
    
    # Calculate first term of W (R(psi)R^{T} + sigma^2)
    first_term <- R %*% (ex_theta * t(R)) + diag(sigma^2, p)
    
    # Decompose first_term using Cholesky decomposition
    S <- chol(first_term)
    qty <- qr.qty(qrz, y)
    qtx <- qr.qty(qrz, X)
    
    M <- forwardsolve(t(S), qty[1:p]) # first p row of S^{-T}Q^{T}y
    N <- forwardsolve(t(S), qtx[1:p, ]) # first p row of S^{-T}Q^{T}X
    
    # Calculate XTWy and XTWX
    XTWy <- t(N) %*% M + t(qtx[(p+1):n, ]) %*% qty[(p+1):n] / sigma^2
    XTWX <- t(N) %*% N + t(qtx[(p+1):n, ]) %*% qtx[(p+1):n, ] / sigma^2
    
    # Check if XTWX is positive definite before Cholesky decomposition
    is_positive_definite_chol <- function(matrix) {
      # If successfully decompose, matrix is positive definite, directly use Cholesky decomposition
      tryCatch({chol(matrix)}, 
      # If fail, matrix isn't positive definite, add a small regularization term before Cholesky decomposition
      error = function(e) {
        chol(matrix + 1e-6 * diag(ncol(matrix))) 
      })
    }
    
    # Decompose (adjusted)XTWX with Cholesky decomposition
    l <- is_positive_definite_chol(XTWX)
    
    beta_hat <- backsolve(l, forwardsolve(t(l), XTWy))
    
    # Calculate log likelihood
    qtyx <- qr.qty(qrz, y - X %*% beta_hat)
    P <- forwardsolve(t(S), qtyx[1:p]) # first p row of S^{-T}Q^{T}(y-X*beta)
    
    neg_loglike <- (sum(P^2) + t(qtyx[(p+1):n]) %*% qtyx[(p+1):n] / sigma^2)/2 +
      sum(log(diag(S))) + (n-p)*log(sigma)
  }
  
  # Create beta_hat as attribute of negative log likelihood
  attr(neg_loglike, 'beta_hat') <- beta_hat
  
  return(neg_loglike)
}



lmm <- function(form, dat, ref=list()) {
  #'
  #' Optimize loglikelihood function with respect to theta
  #' Estimate the model parameters beta and theta
  #'
  #' @param form: A model formula
  #' @param dat: The data frame containing all the variables needed in the model
  #' @param ref: A list of vectors of variable names specifying random effects for Zb part of the model
  #' @return A list containing two elements:
  #'         - theta: the MLEs in theta
  #'         - beta: the MLEs in beta
  #' @approach "BFGS" method to optimize linear model
  #'           "Nelder-Mead" method to optimize linear mixed model
  
  # Set start point of theta
  theta_start <- rep(1, length(ref)+1)
  
  # Set up Z, X, y and f with LMMsetup function
  setup <- LMMsetup(form, dat, ref)
  X <- setup$X
  Z <- setup$Z
  n <- dim(X)[1]
  
  # Calculate the rank of fixed effect matrix
  rank_X <- rankMatrix(X)[1]
  
  # If fixed-effect model matrix is not full rank, raise warning
  if (rank_X < ncol(X)) {
    print("Fixed-effect model matrix is rank deficient")
  }
  
  # Find the max df in linear and linear mixed model
  p <- max(sum(setup$f), dim(X)[2])
  
  # If p > n, stop, can not use this method
  if (p > n) {
    stop("Unable to run, number of columns more than number of observations")
  }

  # If n > p, use the method
  # If there is no random effect, use BFGS method
  if (length(theta_start) == 1) {method = "BFGS"}
  
  # If there are random effects, check the linearly correlation, use Nelder-Mead method
  else {
    # Calculate the rank of random effect matrix
    rank_Z <- rankMatrix(Z)[1]
    # If random-effect model matrix is not full rank, raise warning
    if (rank_Z < ncol(Z)) {
      print("With linearly correlated random effects")
    }
    method = "Nelder-Mead"
  }
  
  # Optimization
  fit <- optim(par = theta_start, fn = LMMprof, X=setup$X, Z=setup$Z, y=setup$y, f=setup$f, 
                method = method)
  
  # Get optimized theta and beta_hat
  theta <- fit$par
  beta_hat <- attr(LMMprof(X=setup$X, Z=setup$Z, y=setup$y, f=setup$f , theta), "beta_hat")
    
  return(list(theta = theta, beta_hat = beta_hat))
  
}