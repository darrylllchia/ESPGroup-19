library(nlme);library(lme4)
lmm(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
exp(lmm(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))$theta)
lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine),data=Machines,
REML=FALSE)

form <- score ~ Machine
dat <- Machines
ref <- list("Worker",c("Worker","Machine"))


LMMsetup <- function(form, dat, ref){
    # set up Z, X and other things that need to be setup only o
    df <- data.frame(dat) 
    # colnames(df) <- c("Column1", "Column2")

    # set up X
    X = model.matrix(form, df)

    # set up y
    y <- df[, names(df) == as.character(form[[2]])]

    # set up Z
    Z <- matrix(nrow = dim(df)[1], ncol = 0)
    f <- c()
    for(i in 1:length(ref)){
        # paste each element in ref to get a string in formula format
        fomula <- paste('~', paste(paste(ref[[i]], collapse = ':'), '-1', sep = ''), sep = '')
        
        model <- model.matrix(as.formula(fomula), df)
        # combine columns of model matrix of each block together
        Z <- cbind(Z, model)
        
        f <- c(f, ncol(model))
    }

    # set up sigma 
    return(list(X = X, Z = Z, y = y, f = f))

}



LMMprof <- function(form, dat, ref, theta){
    # evaluate the (negative) log likelihood for a given theta
    # compute the corresponding beta_hat

    # set up X, Z and y with LMMsetup function
    result <- LMMsetup(form, dat, ref)
    Z <- result$Z
    X <- result$X
    y <- result$y
    f <- result$f
    f_num <- sum(f)
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    
    # QR decompose Z
    qrz <- qr(Z)
    R <- qr.R(qrz)
    
    # extend theta
    # theta <- c(log(sigma), log()) it doesn't matter what theta is
    sigma <- exp(theta[1])
    # calculate first term of RHS of w
    first_term <- R %*% t(R) * rep(exp(theta[2:(length(ref)+1)])^2, f) + diag(sigma^2, f_num)

    # decompose first_term using Cholesky decomposition
    S <- chol(first_term)
    qty <- qr.qty(qrz, y)
    qtx <- qr.qty(qrz, X)
    M <- forwardsolve(t(S), qty[1:f_num])
    N <- forwardsolve(t(S), qtx[1:f_num, ])

    # calculate XTWy and XTWX
    XTWy <- t(N) %*% M + t(qtx[(f_num+1):n, ]) %*% qty[(f_num+1):n] * sigma^2
    XTWX <- t(N) %*% N + t(qtx[(f_num+1):n, ]) %*% qtx[(f_num+1):n, ] * sigma^2

    # decompose XTWX with Cholesky decomposition
    l <- chol(XTWX)
    
    beta_hat <- backsolve(l, forwardsolve(t(l), XTWy))

    # calcualte log likelihood
    qtyx <- qr.qty(qrz, y - X%*%beta_hat)
    P <- forwardsolve(t(S), qtyx[1:f_num])

    loglike <- (t(P) %*% P + t(qtyx[(f_num+1):n]) %*% qtyx[(f_num+1):n] * sigma^2)/
    (sigma^2) / 2 + sum(log(diag(S))) + (n-p)*log(sigma)
    
    attr(loglike, 'beta_hat') <- beta_hat
    
    return(loglike)
}



lmm <- function(form,dat,ref=list()){
    # estimate the model parameters beta and theta
    # return a list containing the MLEs in beta and theta.
    # input:
    # form: a model formula
    # dat: the data frame containing all the variables needed in the model
    # ref: a list of vectors of variable names specifying random effects for Zb part of the model
    block_num <- length(ref)
    theta_start <- rep(0, length(ref)+1)
    
    fit <- optim(theta_start, LMMprof, form = form, dat = dat, ref = ref)
    
    theta <- fit$par
    beta_hat <- attr(LMMprof(form, dat, ref, theta), "beta_hat")
    
    return(list(theta = theta, beta_hat = beta_hat))
}


