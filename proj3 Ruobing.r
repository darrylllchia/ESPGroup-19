library(nlme);library(lme4)

LMMsetup <- function(form, dat, ref) {
  # Create the model matrix X and the response variable y
  X <- model.matrix(form, dat)
  y <- model.response(model.frame(form, dat))
  
  # Create the model matrix Z
  Z <- NULL
  ref_len <- rep(0, length(ref))
  
  for (i in 1:length(ref)) {
    vec <- ref[[i]]
    if (length(vec) == 1) {
      z <- model.matrix(as.formula(paste('~',vec,'-1')),dat)
    } else {
      z <- model.matrix(as.formula(paste('~',paste(vec,collapse = ':'),'-1')),dat)
    }
    Z <- cbind(Z,z)
    ref_len[i] <- ncol(z) 
  }
  # Return 
  return(list(X = X, y = y, Z = Z, ref_len = ref_len))
}

LMMprof <- function(theta, form, dat, ref) {
  setup <- LMMsetup(form, dat, ref)
  X <- setup$X
  Z <- setup$Z
  y <- setup$y
  ref_len <- setup$ref_len
  p <- sum(ref_len)
  n <- dim(X)[1]
  
  qrz <- qr(Z)
  R <- qr.R(qrz)
  qty <- qr.qty(qrz, y)
  qtx <- qr.qty(qrz, X)

  sigma <- exp(theta[1])
  rsirt <- R %*% t(R) * rep(exp(theta[-1])^2, ref_len) + diag(sigma^2, p)
  
  S <- chol(rsirt)
  U <- forwardsolve(t(S), qtx[1:p, ])
  L <- forwardsolve(t(S), qty[1:p])
  XTWX <- t(U) %*% U + t(qtx[(p+1):n, ]) %*% qtx[(p+1):n, ] * sigma^2
  XTWy <- t(U) %*% L + t(qtx[(p+1):n, ]) %*% qty[(p+1):n] * sigma^2
  
  B <- chol(XTWX)
  beta_hat <- backsolve(B, forwardsolve(t(B), XTWy))
  
  qty_xb <- qr.qty(qrz, y - X %*% beta_hat)
  D <- forwardsolve(t(S), qty_xb[1:p])
  
  logll <- (t(D) %*% D + t(qty_xb[(p+1):n]) %*% qty_xb[(p+1):n] * sigma^2) / 
    (2 * sigma^2) + sum(log(diag(S)))  + (n / 2) * log(2 * pi)
  
  attr(logll, 'beta_hat') <- beta_hat
  
  return(logll)
}


lmm <- function(form, dat, ref=list()) {
  
  theta_init <- c(log(sd(y)), rep(log(1), length(ref)))
  
  #result <- optim(par = theta_init,fn = LMMprof,method = "BFGS")
  result <- optim(par = theta_init, LMMprof,form = form, dat = dat, ref = ref,method = "BFGS")
  
  theta <- result$par
  beta_hat <- attr(LMMprof(theta,form, dat, ref), "beta_hat")
  
  return(list(beta_hat = beta_hat, theta = theta))
}

form <- score ~ Machine
dat <- Machines
ref <- list("Worker", c("Worker", "Machine"))
result <- lmm(form, dat, ref)

print(result$theta)
print(result$beta_hat)

lmm(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine),data=Machines,
     REML=FALSE)
