# Aim of this task: Use simulation based method to infer fatal incidence rates from Covid deaths in English hospitals.

# 1. Assign each victim a guessed time of infection
# 2. add a random draw from the infection-to-death distribution to each time of infection to get the implied times of death
# Initially the resulting distribution of times of death will poorly match the real distribution.
# 3. randomly propose to move each time of infection a few days, only accepting the proposed move if the change improves the simulation’s fit to the real death time distribution. 
# 4. This process is simply iterated until we match the death time distribution well.

##read table to get  the first 150 rows of data
setwd("/Users/macbook/Statistical-Programming/assignment")
data <- read.table('engcov.txt')[1:150, ] ## read data, store it as data.frame

t <- data$julian
deaths <- data$nhs

deconv <- function(t,deaths,n.rep=100,bs=FALSE,t0=NULL){
    n <- 29442
    P_vec <- c()
    t0_mat <- matrix(nrow = n, ncol = n.rep)
    inft <- matrix(nrow = 310, ncol = n.rep)
    intr <- matrix(nrow = 310, ncol = length(deaths))

    # 1
    ##find the probability vector of each disease duration based on given log normal density.
    dis_dur <- seq(1, 80) # create disease duration vector
    # calculate the probability vector, meanlog and sdlog are mean and standard deviation of the distribution on the log scale
    prob_vec <- dlnorm(dis_dur, meanlog = 3.152, sdlog = 0.451) 
    norm_prob_vec <- prob_vec / sum(prob_vec) # normalise the probability vector so that the probabilities sum to 1

    if(is.null(t0)){
        # 2
        #create a vector of death day for each individual fatalities
        death_d <- rep(t, deaths)
        # randomly generate n infection-to-death durations with replacement based on probability vector
        duration <- sample(c(1:80), n, replace = TRUE, prob = norm_prob_vec)
        # death day substract infection-to-death duration to get initial guesses for the days of infection
        t0 <- death_d - duration 
    }

    d <- tabulate(death_d, nbins = 310) # create a vector of number of deaths on each day
    plot(c(1:310), d, ylim = c(0, 1200), type = 'l', col = "black", xlab = 'days of the year')
    legend("topright", legend = c('real deaths', ' simulated deaths', 'estimated incidence(10e-4)'), col = c("black", "blue", "red"), lwd = 3)

    for(k in 1:n.rep){
        # 4
        # generate n new draws from the infection-to-death distribution as simulation duaration
        simu_duat <- sample(c(1:80), n, replace = TRUE, prob = norm_prob_vec)
        death_d_s <- t0 + simu_duat # add simulation duaration to t0 to get simulated death days
        d_s <- tabulate(death_d_s, nbins = 310) # create a vector of simulated deaths on each day
        P <- sum((d - d_s)^2) / max(1, d_s) # calculate modified Pearson statistic P
        # Problem: ignore death_d_s<0?
        
        # 5
        t0_dura <- cbind(t0, simu_duat) # combine t0 to duration time to keep infection to death duration fixed
        t0_dura_ran <- t0_dura[sample(length(t0)), ] # randomly order t0
        step <- c()
        if (k >= 1 && k <= 50) {
        step <- c(-8, -4, -2, -1, 1, 2, 4, 8)
        } else if (k >= 51 && k <= 75) {
        step <- c(-4, -2, -1, 1, 2, 4)
        } else {
        step <- c(-2, -1, 1, 2)
        }
        moving <- sample(step, length(t0), replace = TRUE) # create moving vector by randomly choose from step
        death_d_m <- t0_dura_ran[ ,1] + t0_dura_ran[ ,2] + moving # create death day vector(contains death day of each individual fatalities) after moving
        for (i in 1:length(death_d_m)) {
        d_s[death_d_m[i]] <- d_s[death_d_m[i]] + 1 # increase deaths on day i by 1
        P_m <- sum((d - d_s)^2) / max(1, d_s) # calculate P after moving
        if(P_m < P) { #if P after moving decrease, update P and accept move
            P <- P_m
            t0_dura_ran[i,1] <- t0_dura_ran[i,1] + moving[i]
        } 
        else { #if P after moving do not decrease, leave t0 and P unchanged
            d_s[death_d_m[i]] <- d_s[death_d_m[i]] - 1
        }
        }
        t0 <- t0_dura_ran[ ,1]
        lines(c(1:310), d_s/n*10000, col = 'red')
        lines(c(1:310), d_s, col = 'blue')
        inft[ ,k] <- d_s
        P_vec <- c(P_vec,P)
        t0_mat[ ,k] <- t0
    }


    if(bs){
        for(k in 1:length(deaths)){
        deaths <- rpois(length(deaths), deaths)
        if(is.null(t0)){
            # 2
            #create a vector of death day for each individual fatalities
            death_d <- rep(t, deaths)
            # randomly generate n infection-to-death durations with replacement based on probability vector
            duration <- sample(c(1:80), n, replace = TRUE, prob = norm_prob_vec)
            # death day substract infection-to-death duration to get initial guesses for the days of infection
            t0 <- death_d - duration 
        }
        
        # 4
        # generate n new draws from the infection-to-death distribution as simulation duaration
        simu_duat <- sample(c(1:80), n, replace = TRUE, prob = norm_prob_vec)
        death_d_s <- t0 + simu_duat # add simulation duaration to t0 to get simulated death days
        d <- tabulate(death_d, nbins = 310) # create a vector of number of deaths on each day
        d_s <- tabulate(death_d_s, nbins = 310) # create a vector of simulated deaths on each day
        P <- sum((d - d_s)^2) / max(1, d_s) # calculate modified Pearson statistic P
        # Problem: ignore death_d_s<0?
        
        # 5
        t0_dura <- cbind(t0, simu_duat) # combine t0 to duration time to keep infection to death duration fixed
        t0_dura_ran <- t0_dura[sample(length(t0)), ] # randomly order t0
        step <- c()
        if(k>=1 && k<=50){step <- c(-8, -4, -2, -1, 1, 2, 4, 8)}
        else if (k>=51 && k<=75) {step <- c(-4, -2, -1, 1, 2, 4)}
        else {step <- c(-2, -1, 1, 2)}
        moving <- sample(step, length(t0), replace = TRUE) # create moving vector by randomly choose from step
        death_d_m <- rowSums(t0_dura_ran) + moving # create death day vector(contains death day of each individual fatalities) after moving
        for (i in 1:length(death_d_m)) {
            d_s[death_d_m[i]] <- d_s[death_d_m[i]] + 1 # increase deaths on day i by 1
            P_m <- sum((d - d_s)^2) / max(1, d_s) # calculate P after moving
            if(P_m < P) { #if P after moving decrease, update P and accept move
            P <- P_m
            t0_dura_ran[i,1] <- t0_dura_ran[i,1] + moving[i]
            } 
            else { #if P after moving do not decrease, leave t0 and P unchanged
            d_s[death_d_m[i]] <- d_s[death_d_m[i]] - 1
            }
        }
        t0 <- t0_dura_ran[ ,1]
        # plot(c(1:310), deaths, type = 'l', col = "black", xlab = 'days of the year')
        # lines(c(1:310), deaths, col ='black')
        # lines(c(1:310), d_s/n, col = 'red')
        # lines(c(1:310), d_s, col = 'blue')
        intr[ ,k] <- d_s/n
        # inft2 <- cbind(inft, d_s)
        # P_vec2 <- c(P_vec,P)
        # t0_mat2 <- cbind(t0_mat2, t0)
        }
        # calculate the confidence interval of P, inft, t0 after boostrapping to quantify uncertainty
        intr_mean <- apply(intr, 1, mean)
        intr_se <- apply(intr, 1, function(x){sd(x) / sqrt(length(x))})
        lower <- intr_mean - intr_se*1.96
        upper <- intr_mean + intr_se*1.96
    }
    plot(1:310, inft[ ,ncol(inft)]/n, type = 'l', col = 'blue', xlab = "days of the year", ylab = "estimated incidence")
    lines(1:310, lower, col = 'orange')
    lines(1:310, upper, col = 'orange')
    polygon(c(1:310, rev(1:310)), c(lower, rev(upper)), col = rgb(0.5, 0.5, 1, 0.5), border = NA)
    P <- P_vec
    t0 <- t0_mat[, ncol(t0_mat)]
    return(list(P = P, inft = inft, t0 = t0))
}

output <- deconv(t, deaths, n.rep = 100, bs = TRUE, t0 = NULL)
output$P
output$inft
output$t0