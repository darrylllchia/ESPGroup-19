# Jiayi Wu - s2664441, Ruobing Dai - s2655029, Darryl Chia - s2740198
# Darryl Chia and Ruobing Dai completed the process of iteratively fitting the death day, Jiayi Wu was handled the approximate bootstrapping approach part.Then we check the parts completed by other members.
# Finally, all of us explained our parts to each other to understand and improve our work. We discussed on ways to optimise our code and agreed on this final piece of work. 

# Aim of this task: Use simulation based method to infer fatal incidence rates from Covid deaths in English hospitals.

# 1. Assign each victim a guessed time of infection
# 2. Add a random draw from the infection-to-death distribution to each time of infection to get the implied times of death
# 3. Randomly propose to move each time of infection a few days, only accepting the proposed move if the change improves the simulation’s fit to the real death time distribution. 
# 4. Iterate this process until we have a good match for the death time distribution.
# 5. Plot a graph showing the actual death time, simulated death time, and infection rate, adding each new death time line as it iterates.
# 6. Use an approximate bootstrapping approach applied once the method has converged to quantify its uncertainty.
# (Use real death data as expected values for Poisson variables, then simulate Poisson data based on these values to replace the real data.)
# 7. Plot a graph to illustrate the fit of the simulation model to the real death data.

setwd("/Users/macbook/Statistical-Programming/assignment")
data <- read.table('engcov.txt')[1:150, ] 

t <- data$julian
deaths <- data$nhs

deconv <- function(t,deaths,n.rep=100,bs=FALSE,t0=NULL,n.times = 100){
  
# The purpose of this function is to display the estimated death and incidence rate of covid 19 compared to the actual deaths in England. 
# Data on daily death numbers is compared against simulated data sampled from a distribution with set parameters.
# The loss function (P) involved is proportional to the squared difference between the simulated deaths and the actual deaths each day. 
# The process repeats for a set number of times (n.rep). Over time, P is expected to decrease.
# Bootstrapping also helps to see any uncertainty with our simulation and can be activated with the logical input to the function 'bs'.
# It compares the converged vector of days of death(t0) with data generated from a poisson distribution, using 
# the actual number of deaths of each day as the lambda parameter for each day.
# The function takes in the number of deaths on each day, the days to be analysed and hyperparameters bs, n.rep and t0.
# The function outputs a list containing a matrix(inft) of the simulated number of deaths after each iteration,
# a vector of Ps after each iteration and the converged t0.

  death_d <- rep(t, deaths) # create a vector of death day for each individual fatalities
  n <- length(death_d) # number of fatalities
  max_day <- 310 # set the max observation day
  P_vec <- c() # initialise vector to store each P
  t0_mat <- matrix(nrow = n, ncol = n.rep) # initialise matrix to store each t0 after each full update
  inft <- matrix(nrow = max_day, ncol = n.rep) # initialise matrix to store the number of new infections per day 
  intr <- matrix(nrow = max_day, ncol = n.times) # initialise matrix to store the number of new infections per day after bootstrapping
  
  # 1 find the probability vector of each disease duration based on given log normal density
  dis_dur <- seq(1, 80) # create disease duration vector
  # calculate the probability vector, meanlog and sdlog are mean and standard deviation of the distribution on the log scale
  prob_vec <- dlnorm(dis_dur, meanlog = 3.152, sdlog = 0.451) 
  norm_prob_vec <- prob_vec / sum(prob_vec) # normalise the probability vector so that the probabilities sum to 1
  
  if(is.null(t0)){
    # 2
    #create a vector of day of death for each individual fatality
    death_d <- rep(t, deaths)
    # randomly generate n infection-to-death durations with replacement based on probability vector
    duration <- sample(c(1:80), n, replace = TRUE, prob = norm_prob_vec)
    # infection-to-death duration subtracted from death day to get initial estimates for the days of infection
    t0 <- death_d - duration
    t0[t0 <= 0] <- 1 # set the non-positive infection day as 1
  }
  
  d <- tabulate(death_d, nbins = max_day) # create a vector of number of deaths on each day
  # draw a plot of relationship between real deaths and days of the year
  plot(c(1:max_day), d, ylim = c(0, 1700), type = 'l', lwd = 3, col = "black", xlab = 'Days of the year', ylab = "Number of deaths/incidence", main = "Deaths/incidence each day until converge") 
  legend("topright", legend = c('Real deaths', 'Simulated deaths', 'Estimated incidence'), col = c("black", "blue", "red"), lwd = 3) 
  
  for(k in 1:n.rep){
    # 4 
    # generate n new draws from the infection-to-death distribution as simulation duaration
    simu_duat <- sample(c(1:80), n, replace = TRUE, prob = norm_prob_vec)
    death_d_s <- t0 + simu_duat # add simulation duration to t0 to get simulated death days
    d_s <- tabulate(death_d_s, nbins = max_day) # create a vector of simulated deaths on each day
    P <- sum((d - d_s)^2 / pmax(1, d_s)) # calculate modified Pearson statistic P
    
    # 5. randomly order each time of infection and move it a few days
    t0_dura <- cbind(t0, simu_duat) # combine t0 to duration time to keep infection to death duration fixed
    t0_dura_ran <- t0_dura[sample(length(t0)), ] # randomly order t0
    
    # change the moving step to improve convergence rate
    step <- c()
    if (k >= 1 && k <= 50) {
      step <- c(-8, -4, -2, -1, 1, 2, 4, 8)
    } else if (k >= 51 && k <= 75) {
      step <- c(-4, -2, -1, 1, 2, 4)
    } else {
      step <- c(-2, -1, 1, 2)
    }
    
    moving <- sample(step, length(t0), replace = TRUE) # create moving vector by randomly choosing from step
    death_d_bm <- t0_dura_ran[ ,1] + t0_dura_ran[ ,2] # create death day of each individual before moving
    death_d_m <- t0_dura_ran[ ,1] + t0_dura_ran[ ,2] + moving # create death day of each individual after moving
    
    # for each elements in t0, accept the proposed move if the change improves the simulation’s fit to the real death time distribution, ie. reduce P
    for (i in 1:length(death_d_m)) {
      dead_bef_move <- death_d_bm[i] # the death day before move
      dead_aft_move <- death_d_m[i] # the death day after move
      # calculate change in P: only need to calculate the difference in deaths with the actual days for the 2 days involved in the moving
      # P_pri is before moving, P_for is after moving
      P_pri <- (d_s[dead_bef_move] - d[dead_bef_move])^2 / max(1, d_s[dead_bef_move]) + (d_s[dead_aft_move] - d[dead_aft_move])^2 / max(1, d_s[dead_aft_move])
      d_s[dead_aft_move] <- d_s[dead_aft_move] + 1 # increase deaths on day after moving by 1
      d_s[dead_bef_move] <- d_s[dead_bef_move] - 1 # decrease deaths on day before moving by 1
      P_for <- (d_s[dead_bef_move] - d[dead_bef_move])^2 / max(1, d_s[dead_bef_move]) + (d_s[dead_aft_move] - d[dead_aft_move])^2 / max(1, d_s[dead_aft_move])
      P_m <- P + P_for - P_pri # calculate P after moving
      if(P_m < P) { # if P decreases after moving, update P and accept move
        P <- P_m
        t0_dura_ran[i,1] <- t0_dura_ran[i,1] + moving[i]
      } 
      else { #if P does not decrease after moving, leave t0 and P unchanged
        d_s[dead_aft_move] <- d_s[dead_aft_move] - 1
        d_s[dead_bef_move] <- d_s[dead_bef_move] + 1
      }
    }
    
    t0 <- t0_dura_ran[ ,1] 
    inci <- tabulate(t0, nbins = max_day)
    P_vec <- c(P_vec,P)
    inft[ ,k] <- inci
    t0_mat[ ,k] <- t0
    lines(c(1:max_day), inci, col = 'red')
    lines(c(1:max_day), d_s, col = 'blue')
    
  }
  
  #
  if(bs){
    plot(1:max_day, inft[ ,ncol(inft)], ylim = c(0, 1700), type = 'l', col = 'blue', lwd = 4, xlab = "Days of the year", ylab = "Number of deaths/incidence", main = "Deaths/incidence each day in bootstrapping")
    legend("topright", legend = c('Real incidence', 'Simulated incidence', 'Real deaths', 'The first day of UK lockdown'), col = c("blue", "red", "grey", 'black'), lwd = 3)
    lines(1:max_day, d, col = 'grey', lwd = 3)
    abline(v = 84, col = "black", lwd = 2, lty = 2)
    text(x = 84, y = -50, labels = "84", pos = 3, col = "black")
    death_simulation <- sapply(deaths, function(x) rpois(n.times, x))
    
    for(k in 1:n.times){
      
      deaths_new <- death_simulation[k, ]
      
      if(is.null(t0)){
        # 2
        #create a vector of death day for each individual fatality
        death_d <- rep(t, deaths_new)
        # randomly generate n infection-to-death durations with replacement based on probability vector
        duration <- sample(c(1:80), n, replace = TRUE, prob = norm_prob_vec)
        # infection-to-death duration subtracted from death day to get initial estimates for the days of infection
        t0 <- death_d - duration 
      }
      
      death_d <- rep(t, deaths_new)
      # 4 
      # generate n new draws from the infection-to-death distribution as simulation duaration
      simu_duat <- sample(c(1:80), n, replace = TRUE, prob = norm_prob_vec)
      death_d_s <- t0 + simu_duat # add simulation duration to t0 to get simulated death days
      d_s <- tabulate(death_d_s, nbins = max_day) # create a vector of simulated deaths on each day
      P <- sum((d - d_s)^2 / pmax(1, d_s)) # calculate modified Pearson statistic P
      
      # 5. randomly order each time of infection and move it a few days
      t0_dura <- cbind(t0, simu_duat) # combine t0 to duration time to keep infection to death duration fixed
      t0_dura_ran <- t0_dura[sample(length(t0)), ] # randomly order t0
      
      # change the moving step to improve convergence rate
      step <- c()
      if (k >= 1 && k <= 50) {
        step <- c(-8, -4, -2, -1, 1, 2, 4, 8)
      } else if (k >= 51 && k <= 75) {
        step <- c(-4, -2, -1, 1, 2, 4)
      } else {
        step <- c(-2, -1, 1, 2)
      }
      
      moving <- sample(step, length(t0), replace = TRUE) # create moving vector by randomly choosing from step
      death_d_bm <- t0_dura_ran[ ,1] + t0_dura_ran[ ,2] # create death day of each individual before moving
      death_d_m <- t0_dura_ran[ ,1] + t0_dura_ran[ ,2] + moving # create death day of each individual after moving
      
      # for each elements in t0, accept the proposed move if the change improves the simulation’s fit to the real death time distribution, ie. reduce P
      for (i in 1:length(death_d_m)) {
        dead_bef_move <- death_d_bm[i] # the death day before move
        dead_aft_move <- death_d_m[i] # the death day after move
        # calculate change in P: only need to calculate the difference in deaths with the actual days for the 2 days involved in the moving
        # P_pri is before moving, P_for is after moving
        P_pri <- (d_s[dead_bef_move] - d[dead_bef_move])^2 / max(1, d_s[dead_bef_move]) + (d_s[dead_aft_move] - d[dead_aft_move])^2 / max(1, d_s[dead_aft_move])
        d_s[dead_aft_move] <- d_s[dead_aft_move] + 1 # increase deaths on day after moving by 1
        d_s[dead_bef_move] <- d_s[dead_bef_move] - 1 # decrease deaths on day before moving by 1
        P_for <- (d_s[dead_bef_move] - d[dead_bef_move])^2 / max(1, d_s[dead_bef_move]) + (d_s[dead_aft_move] - d[dead_aft_move])^2 / max(1, d_s[dead_aft_move])
        P_m <- P + P_for - P_pri # calculate P after moving
        if(P_m < P) { # if P decreases after moving, update P and accept move
          P <- P_m
          t0_dura_ran[i,1] <- t0_dura_ran[i,1] + moving[i]
        } 
        else { #if P does not decrease after moving, leave t0 and P unchanged
          d_s[dead_aft_move] <- d_s[dead_aft_move] - 1
          d_s[dead_bef_move] <- d_s[dead_bef_move] + 1
        }
      }
      
      inci <- tabulate(t0_dura_ran[ ,1], nbins = max_day)
      lines(c(1:max_day), inci, col = 'red')
      intr[ ,k] <- inci
      
    }
  }
  P <- P_vec
  t0 <- t0_mat[, ncol(t0_mat)]
  return(list(P = P, inft = inft, t0 = t0))
}

output <- deconv(t, deaths, n.rep = 100, bs = TRUE, t0 = NULL, n.times = 100)
output$P
output$inft
output$t0
