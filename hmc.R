library(numDeriv)
library(mvtnorm)

# We will sample from bivariate MVN(c(0,0), Sigma), with a very high correlation.

rho <- 0.95
Sigma <- matrix(c(1, rho, rho, 1), nrow=2)
Sigma_inv <- solve(Sigma)

# As in Neal's tutorial
       
# U is minus the log likelihood
U <- function(q) {
  return(as.numeric(matrix(q, nrow=1) %*% Sigma_inv %*% matrix(q, nrow=2))/2 + log(2*pi) + 0.5*log(det(Sigma)))
}

# Gradient of U
grad_U <- function(q) {
  return(as.numeric(Sigma_inv %*% matrix(q, nrow=2)))
}

# HMC Sampler
HMC <- function(U, grad_U, epsilon, L, current_q) {
  
  # Current "position"
  q <- current_q
  
  # Sample a momentum and initialize
  p <- rnorm(length(q), 0, 1)
  current_p <- p
  
  # Make the first half step
  p <- p - epsilon*grad_U(q) / 2
  
  # Do the leapfrog part
  for (i in 1:L) {
    q <- q + epsilon*p
    if (i < L) p <- p - epsilon*grad_U(q)
  }
  
  # Make the last half step
  p <- p - epsilon*grad_U(q) / 2
  p <- -p
  
  # Accept or reject the proposal (q)
  current_U <- U(current_q)
  current_K <- sum(current_p^2) / 2
  proposed_U <- U(q)
  proposed_K <- sum(p^2) / 2
  if (runif(1) < exp(proposed_U+proposed_K-current_U-current_K))
    return(list(new_val=q, accepted=1))
  else
    return(list(new_val = current_q, accepted=0))
}

# Random walk Metropolis Sampler
Metropolis_sampler <- function(current_q, proposal_sd) {
  
  q <- current_q + rnorm(2, 0, proposal_sd)

  # Accept or reject the proposal (q)
  num <- -U(q)
  den <- -U(current_q)
  
  if (runif(1) < exp(num-den))
    return(list(new_val=q, accepted=1))
  else
    return(list(new_val = current_q, accepted=0))
}


# Step size
epsilon <- 1/30

# Number ofleapfrog jumps
L <- 30

init <- c(0,0)

# Number of posterior sample
n_mcmc <- 10000

# Initialize chain
sample_hmc <- matrix(0, nrow=n_mcmc, ncol=2)
sample_hmc[1,] <- init
accepted_hmc <- 0

# Sample from target distribution
for(i in 2:n_mcmc) {
  temp <- HMC(U=U, grad_U=grad_U, epsilon=epsilon, L=L, current_q = sample_hmc[i-1,])
  sample_hmc[i,] <- temp$new_val
  accepted_hmc <- accepted_hmc + temp$accepted
}

# Verify that the values drawn are actually from the target distribution
plot(sample_hmc, xlim=c(-3, 3), ylim=c(-3, 3))
cor(sample_hmc)

# Acceptance rate (should be very high)
accepted_hmc / n_mcmc

## Do the same thing for the random walk Metropolis sampler
# Initialize chain
sample <- matrix(0, nrow=n_mcmc, ncol=2)
sample[1,] <- init
accepted <- 0
proposal_sd <- 0.1

# Sample from target distribution
for(i in 2:n_mcmc) {
  temp <- Metropolis_sampler(sample[i-1,], proposal_sd)
  sample[i,] <- temp$new_val
  accepted <- accepted + temp$accepted
}

# Verify that the sample correlation is quite different from the true rho (0.95)
plot(sample, xlim=c(-3, 3), ylim=c(-3, 3))
cor(sample)

# Acceptance rate (lower than in HMC)
accepted / n_mcmc