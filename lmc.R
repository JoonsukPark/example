set.seed(100)

## Langevin Monte Carlo vs Random-walk Metropolis-Hastings MCMC
# We will sample from bivariate MVN(c(0,0), Sigma), with a very high correlation.

rho <- 0.95
Sigma <- matrix(c(1, rho, rho, 1), nrow=2)
Sigma_inv <- solve(Sigma)

# U = loglik (compare this with the case of HMC where U = -loglik)
U <- function(q) {
  return(matrix(q, nrow=1) %*% Sigma_inv %*% matrix(q, nrow=2))/2 + log(2*pi) + 0.5*log(det(Sigma))
}

# Gradient of U
grad_U <- function(q) {
  return(-Sigma_inv %*% matrix(q, nrow=2))
}

# Langevin Monte Carlo sampler
LMC <- function(current_q, epsilon, proposal_sd) {
  
  q <- current_q + epsilon*grad_U(current_q) + sqrt(2*epsilon)*rnorm(length(current_q), 0, proposal_sd)
  
  # Accept or reject the proposal (q)
  num <- -U(q)
  den <- -U(current_q)
  
  if (runif(1) < exp(num-den))
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

# Number of posterior sample
n_mcmc <- 100

## Random walk MCMC
# Initialize chain
init <- c(-10, 10)
sample <- matrix(0, nrow=n_mcmc, ncol=2)
sample[1,] <- init
accepted <- 0

# Sample from the target distribution
for(i in 2:n_mcmc) {
  temp <- Metropolis_sampler(sample[i-1,], proposal_sd)
  sample[i,] <- temp$new_val
  accepted <- accepted + temp$accepted
}

# Verify that the sample correlation is quite different from the true rho (0.95)
sample <- sample[(n_mcmc/2+1):n_mcmc,]
plot(sample, xlim=c(-3, 3), ylim=c(-3, 3))
cor(sample)

# Acceptance rate of Random-walk MCMC
accepted / n_mcmc

## LMC
# Initialize chain
sample_LMC <- matrix(0, nrow=n_mcmc, ncol=2)
sample_LMC[1,] <- init
accepted_LMC <- 0
epsilon <- 0.05
proposal_sd <- 0.5

# Sample from target distribution
for(i in 2:n_mcmc) {
  temp <- LMC(current_q = sample_LMC[i-1,], epsilon=epsilon, proposal_sd=proposal_sd)
  sample_LMC[i,] <- temp$new_val
  accepted_LMC <- accepted_LMC + temp$accepted
}

# Not as good as HMC, but better than random walk MCMC
sample_LMC <- sample_LMC[(n_mcmc/2+1):n_mcmc,]
plot(sample_LMC, xlim=c(-3, 3), ylim=c(-3, 3))
cor(sample_LMC)

# Acceptance rate of LMC
accepted_LMC / n_mcmc
