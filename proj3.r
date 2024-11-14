



LMMsetup <- function(form, dat, ref){
  # The main purpose of this function is to initialise and set up our fixed effect and random effect matrices X and Z,
  # along with the y vector. The vector f contains the number of columns of each block in each term in ref. So the
  # first term in ref will correspond to the first f[1] columns in Z and so on. Each block will have the same standard deviation
  
  
  # Convert the input data into a data frame format
  df <- data.frame(dat) 
  
  # Set up X (Fixed effects)
  X = model.matrix(form, df)
  
  # Set up y (Response variable)
  y <- df[, names(df) == as.character(form[[2]])]
  
  # Set up Z (Random Effects)
  Z <- matrix(nrow = dim(df)[1], ncol = 0)
  f <- c()

  if(length(ref) > 0){
    for(i in 1:length(ref)){
      # Paste each element in ref to get a string in formula format
      fomula <- paste('~', paste(paste(ref[[i]], collapse = ':'), '-1', sep = ''), sep = '')
      model <- model.matrix(as.formula(fomula), df)
      # Combine columns of model matrix of each block together
      Z <- cbind(Z, model)
      # Number of columns of each block
      f <- c(f, ncol(model))
    }
  }
  
  # Set up sigma 
  return(list(X = X, Z = Z, y = y, f = f))
  
}


# theta <- c(-0.03932876,1.42168890,1.22317100)
LMMprof <- function(form, dat, ref, theta){
  # Evaluate the (negative) log likelihood for a given theta
  # Compute the corresponding beta_hat
  # @param theta: vector of log(std_devs). First element of theta is the log(std_dev) of residuals and subsequent elements
  # correspond to the log(std_devs) of each block in Z
  # returns: negative log likelihood equation for optim to minimise and the calculated beta_hat vector given a theta vector
  
  # Set up X, Z and y with LMMsetup function
  result <- LMMsetup(form, dat, ref)
  Z <- result$Z
  X <- result$X
  y <- result$y
  f <- result$f
  # Number of columns in the Z matrix
  p <- sum(f)
  # Number of data samples
  n <- dim(X)[1]
  # Number of fixed effects
  p1 <- dim(X)[2]
  
  sigma <- exp(theta[1])
  
  # If there are no random effects 
  if(length(f) == 0){
    # Compute beta_hat using ordinary least squares estimates
    qrx <- qr(X)
    beta_hat <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p1])
    neg_loglike <- t(y-X%*%beta_hat) %*% (y-X%*%beta_hat) / (2*sigma^2) + n*log(sigma)
  }
  else{
    # QR decompose Z
    qrz <- qr(Z)
    R <- qr.R(qrz)
    
    # extend theta
    # theta <- c(log(sigma), log()) it doesn't matter what theta is
    
    # Calculate first term of RHS of W
    first_term <- R %*% t(R) * rep(exp(theta[-1])^2, f) + diag(sigma^2, p)
    
    # Decompose first_term using Cholesky decomposition
    S <- chol(first_term)
    qty <- qr.qty(qrz, y)
    qtx <- qr.qty(qrz, X)
    M <- forwardsolve(t(S), qty[1:p])
    N <- forwardsolve(t(S), qtx[1:p, ])
    
    # Calculate XTWy and XTWX
    XTWy <- t(N) %*% M + t(qtx[(p+1):n, ]) %*% qty[(p+1):n] / sigma^2
    XTWX <- t(N) %*% N + t(qtx[(p+1):n, ]) %*% qtx[(p+1):n, ] / sigma^2
    # XTWX <- XTWX + 1e-6 * diag(ncol(XTWX))
    
    # Decompose XTWX with Cholesky decomposition
    l <- chol(XTWX)
    
    beta_hat <- backsolve(l, forwardsolve(t(l), XTWy))
    
    # Calcualte log likelihood
    qtyx <- qr.qty(qrz, y - X%*%beta_hat)
    P <- forwardsolve(t(S), qtyx[1:p])
    neg_loglike <- (sum(P^2) + t(qtyx[(p+1):n]) %*% qtyx[(p+1):n] / sigma^2)/2 +
      sum(log(diag(S))) + (n-p)*log(sigma)
    
    #neg_loglike <- (t(P) %*% P + t(qtyx[(p+1):n]) %*% qtyx[(p+1):n] * sigma^2)/2 +
    #sum(log(diag(S))) + (n-p)*log(sigma) +  n*log(2*pi)/2
    
    #neg_loglike <- (sum(P^2)*2 + t(qtyx[(p+1):n]) %*% qtyx[(p+1):n] * sigma^2)/2 +
    #sum(log(diag(S))) + (n-p)*log(sigma) +  n*log(2*pi)/2
    #print(theta)
  }
  # Append beta_hat to the attributes of the negative log likelihood
  attr(neg_loglike, 'beta_hat') <- beta_hat
  return(neg_loglike)
}



lmm <- function(form,dat,ref=list()){
  # estimate the model parameters beta and theta
  # return a list containing the MLEs in beta and theta.
  # input:
  # form: a model formula
  # dat: the data frame containing all the variables needed in the model
  # ref: a list of vectors of variable names specifying random effects for Zb part of the model
  
  # Set the initial theta vector (number of blocks + sigma)
  block_num <- length(ref)
  theta_start <- rep(0, length(ref)+1)

  
  setup <- LMMsetup(form, dat, ref)
    X <- LMMsetup(form, dat, ref)$X
    n <- dim(X)[1]
    p <- max(sum(setup$f), dim(X)[2])
    if(p>n) {
      return("Unable to run, number of columns more than number of observations")
    }
    else {
      # If there is no random effect, use the BFGS method
      if (length(theta_start) == 1) {
        method = "BFGS"} 
      # If there are random effects, use the Nelder-Mead method
      else {method = "Nelder-Mead"}
      
      # Use Nelder-Mead or another method to do optimization
      fit <- optim(theta_start, LMMprof, form = form, dat = dat, ref = ref, 
                   method = method)
      # Get optimized theta and beta_hat
      theta <- fit$par
      beta_hat <- attr(LMMprof(form, dat, ref, theta), "beta_hat")
      
      return(list(fit = fit, theta = theta, beta_hat = beta_hat))
    }

}

library(nlme);library(lme4)
lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine),data=Machines,
     REML=FALSE)
result <- lmm(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
result
exp(result$theta)
lmm(score ~ Machine,Machines,list())
reg <- lm(score ~ Machine,Machines)
summary(reg)
summary(reg)$sigma




data(package = 'lme4')
head(Milk)
Milk$Diet[1:100]
Milk$Cow
lmer(protein ~ Diet + (1|Cow:Diet),data=Milk,REML=FALSE)
result <- lmm(protein ~ Diet,Milk,list(c("Cow","Diet")))
exp(result$theta)
