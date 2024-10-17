setwd("C:/Users/HUAWEI/Desktop/sds/extend statisting programming/作业2")
getwd()
data <- read.table('engcov.txt')[1:150,] # Read the first 150 rows of data

# 1
infect_to_death <- dlnorm(1:80, meanlog = 3.152, sdlog = 0.451) # Probability distribution of time from infection to death
itd_prob <- infect_to_death / sum(infect_to_death) # Normalize the probabilities to sum to 1

# 2
day <- data$julian # Extract the number of days to death (julian column)
num <- data$nhs # Extract the number of deaths (nhs column)
death_days <- rep(day, num) # Each death corresponds to a date of death
n <- length(death_days) # Total number of deaths
infect_to_death_days <- sample(1:80, size = n, replace = TRUE, prob = itd_prob) # Randomly sample infection-to-death days
t0 <- death_days - infect_to_death_days # Initial guess of infection time

# 4
pearson <- function(actual_death, simulate_death) {
  p <- sum((actual_death - simulate_death)^2 / pmax(1, simulate_death))
  return(p)
} # Define the Pearson Function

# 5
t0_new <- function(t0, p, num, itd_prob, death_days, move_sample) {
  t0_n <- t0 # Store the updated t0 outside the loop
  for (i in seq_along(t0)) {
    t0_n[i] <- pmax(1, pmin(t0[i] + move_sample[i], 310)) # Adjust infection time
  }
  # Generate simulated death days for the entire new t0
  simulate_death_days_new <- t0_n + sample(1:80, size = length(t0), replace = TRUE, prob = itd_prob)
  simulate_death_num_new <- tabulate(simulate_death_days_new, nbins = max(death_days)) # New simulated number of deaths per day
  p_new <- pearson(num, simulate_death_num_new[death_days]) # New Pearson statistic
  
  return(list(t0 = t0_n, p = p_new))
}

deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL) {
  death_d <- rep(t, deaths)
  if (is.null(t0)) {
    infect_to_death_days <- sample(1:80, size = length(death_d), replace = TRUE, prob = itd_prob)
    t0 <- death_d - infect_to_death_days
  }

  itd_days <- sample(1:80, size = length(death_d), replace = TRUE, prob = itd_prob) # sampled infection-to-death days

  P <- numeric(n.rep)
  inft <- matrix(NA, nrow = 310, ncol = n.rep)
  
  for (j in 1:n.rep) {
    if (bs) {
      deaths <- rpois(length(deaths), lambda = deaths)
    }
    
    move_sample <- sample(c(-4, -2, -1, 1, 2, 4), size = length(t0), replace = TRUE) # move samples
    
    # Update t0 and p
    tp_list <- t0_new(t0, p, deaths, itd_prob, death_d, move_sample)
    t0 <- tp_list$t0
    P[j] <- tp_list$p
    inft[, j] <- tabulate(t0, nbins = 310) # Store the result
    #cat("Iteration:", j, " Pearson statistic:", P[j], "\n") # Output progress
  }
  
  return(list(p = P, inft = inft, t0 = t0))
}

tp <- deconv(day, num, n.rep = 100)

mean_inft <- rowMeans(tp$inft) 
plot(1:310, mean_inft, type = "l", ylim = range(c(mean_inft, ci)), 
     ylab = "Estimated Infections", xlab = "Days", col = "blue", lwd = 2)
points(day, num, col = "red", pch = 20)
abline(v = 84, col = "black", lty = 2)
legend("topright", legend = c("Estimated Infections", "95% CI", "Actual Deaths", "UK Lockdown"), 
       col = c("blue", rgb(0.1, 0.1, 0.9, 0.2), "red", "black"), 
       lty = c(1, NA, NA, 2), pch = c(NA, 15, 20, NA), lwd = 2, pt.cex = 1.5)
