library(nlme);library(lme4)
lmm(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
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
    for(i in length(ref)){
        # paste each element in ref to get a string in formula format
        fomula <- paste('~', paste(paste(ref[[i]], collapse = ':'), '-1', sep = ''), sep = '')

        # combine columns of model matrix of each block together
        Z <- cbind(Z, model.matrix(as.formula(fomula), df))
    }

    # set up sigma 
    # reg <- lm(form, df)
    # sigma <- sqrt(sum(reg$residuals^2) / (reg$df.residual+1))
    return(list(X = X, Z = Z, y = y, sigma = sigma))

}



LMMprof <- function(form, dat, ref, theta){
    # evaluate the (negative) log likelihood for a given theta
    # compute the corresponding beta_hat

    # set up X, Z and y with LMMsetup function
    result <- LMMsetup(form, dat, ref)
    Z <- result$Z
    X <- result$X
    y <- result$y
    n <- dim(X)[1]
    p <- dim(X)[2]

    # QR decompose Z
    qrz <- qr(Z)
    R <- qr.R(qrz)

    # calculate first term of RHS of w
    first_term <- R %*% t(R) * rep(exp(theta)^2, p) + diag(sigma^2, p)

    # decompose first_term using Cholesky decomposition
    S <- chol(first_term)
    qty <- qr.qty(qrz, y)
    qtx <- qr.qty(qrz, X)
    M <- forwardsolve(t(S), qty[1:p,])
    N <- forwardsolve(t(S), qtx[1:p,])

    # calculate XTWy and XTWX
    XTWy <- rbind(t(N) %*% M, t(qtx[(p+1):n,]) %*% qty[(p+1):n,] * sigma^2)
    XTWX <- rbind(t(N) %*% N, t(qtx[(p+1):n,]) %*% qtx[(p+1):n,] * sigma^2)

    # decompose XTWX with Cholesky decomposition
    l <- chol(XTWX)

    beta_hat <- function(form, dat, ref, theta){
        beta_hat <- backwardsolve(l, forwardsolve(t(l), XTWy))
        beta_l <- c(beta_l, list(bata_hat))
        return(bata_hat)
    }

    # calcualte log likelihood
    qtyx <- qr.qty(qrz, y - X%*%bata_hat)
    P <- forwardsolve(t(S), qtyx[1:p,])

    neg_loglike <- function(form, dat, ref, theta){
        loglike <- rbind(t(P) %*% P, t(qtyx[(p+1):n,]) %*% qtyx[(p+1):n,] * sigma^2)/
        (sigma^2) / 2 + sum(log(diag(S))) + (n-p)*log(sigma)
        return(loglike)3
    }

    return(list(neg_loglike = neg_loglike, beta_hat = beta_hat))
}



lmm <- function(form,dat,ref=list()){
    # estimate the model parameters beta and theta
    # return a list containing the MLEs in beta and theta.
    # input:
    # form: a model formula
    # dat: the data frame containing all the variables needed in the model
    # ref: a list of vectors of variable names specifying random effects for Zb part of the model
    theta_l <- list()
    beta_l <- list()
    tar_fun <- function(form, dat, ref, theta){
        neg_loglike <- LMMprof(form, dat, ref, theta)$neg_loglike
        attr(neg_loglike, 'beta_hat') <- LMMprof(form, dat, ref, theta)$beta_hat
        theta_l <- c(params, list(par))
        beta_l <- c(beta_l, list(beta_hat))
    }
    sigma <- LMMsetup(form,dat,ref)$sigma
    fit <- optim(log(sigma), LMMprof, form = form, dat = dat, ref = ref)
    return(list(beta = beta_l, theta = theta_l))
}


