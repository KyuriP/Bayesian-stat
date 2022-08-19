##============================================================
## 1. Preparation
##============================================================
## Required packages to run the following code
library(scales) # used in section 4: assess convergence (to tone down the color in the plots)
library(MASS) # used in section 8: Bayes factor
library(brms) # used in section 9: Bayesian causal inference

## Load "Birthweight" data
dat <- read.csv("Birthweight.csv")
## Subset the data: only using 3 variables for the main analysis
## (i.e., birth weight of babies, number of cigarettes smoked by mothers, and mother's age)
BW <- dat[,c(3,8,7)]
## Get the descriptive statistics
descriptives <- psych::describe(BW)[,c(2,3,4,5,8,9)]
## Center the predictors (it helps with reducing dependencies between them)
BW[,2:3] <- apply(BW[,-1], 2, function(x) x - mean(x))
## Rename the DV and IVs for convenience (DV = birthweight, x1 = number of cigarettes, x2 = mother's age)
colnames(BW) <- c("y", "x1", "x2")

##============================================================
## 2. Gibbs Sampler (with Random walk MH step) function
##============================================================
## First, I built mini-functions that sample each of the parameters,
## which will be nested under a big sampling function that loops over these mini-functions.

## Gibbs sampling mini-function for tau (precision)
sample_tau <- function(data, b0, b1 = 0, b2 = 0, alpha0, beta0){
  rgamma(1, # draw a random sample from the gamma distribution
         shape = alpha0 + nrow(data)*0.5, # specify shape parameter 
         rate = beta0 + ((1/2) * sum((data$y - (b0 + as.matrix(data$x1) %*% b1 + as.matrix(data$x2) %*% b2))^2))) # specify rate parameter
}

## Gibbs sampling mini-function for b0
sample_b0 <- function(data, b1 = 0, b2 = 0, tau, mu00, tau00) {
  prec <- tau00 + tau * nrow(data) # specify the precision (tau)
  mean <- ((tau00 * mu00) + (tau * sum(data$y - as.matrix(data$x1) %*% b1 - as.matrix(data$x2) %*% b2))) / prec # specify the mean
  rnorm(1, mean = mean, sd = 1 / sqrt(prec)) # draw a random sample from the normal dist.
}

## Random Walk MH-sampling mini-function for b1 : proposal density = normal, prior = t-distribution
MH_sample_b1 <- function(data, b0, b1, b2, tau, nu10, mu10, tau10){
  proposed.b1 <- rnorm(1, b1, 0.01) # draw a candidate value from the (normal) proposal density with mean of (current) b1 and sd = 0.01 (based on trial and error)
  prec <- tau10 + tau * sum(data$x1 * data$x1) # specify the precision
  mean <- ((tau10*mu10) + tau * sum((data$y - b0 - as.matrix(data$x2) %*% b2) * data$x1)) / prec # specify the mean
  prop.like <- dnorm(proposed.b1, mean = mean, sd = 1 / sqrt(prec)) # likelihood for the proposed b1
  current.like <- dnorm(b1, mean = mean, sd = 1 / sqrt(prec)) # likelihood for the current b1
  prop.prior <- (1 + ((proposed.b1 - mu10)^2 * tau10 /nu10))^(-(nu10 + 1) / 2) # prior for proposed b1
  current.prior <- (1 + ((b1 - mu10)^2 * tau10 /nu10))^(-(nu10 + 1) / 2) # prior for current b1
  A <- (prop.like * prop.prior)/(current.like * current.prior) # as proposal density is symmetric, the acceptance ratio (A) is simplified to: P(proposal)/P(current), where P = posterior probability density
  ifelse( A > runif(1), proposed.b1, b1) # accept proposed b1, if acceptance ratio > randomly drawn number; reject otherwise.
}

## Gibbs sampling mini-function for b1
sample_b1 <- function(data, b0, b2 = 0, tau, mu10, tau10) {
  prec <- tau10 + tau * sum(data$x1 * data$x1) # specify the precision
  mean <- ((tau10*mu10) + tau * sum((data$y - b0 - as.matrix(data$x2) %*% b2) * data$x1)) / prec # specify the mean
  rnorm(1, mean = mean, sd = 1 / sqrt(prec)) # draw a random sample from the normal dist.
}

## Gibbs sampling mini-function for b2
sample_b2 <- function(data, b0, b1 = 0, tau, mu20, tau20) {
  prec <- tau20 + tau * sum(data$x2 * data$x2) # specify the precision
  mean <- ((tau20*mu20) + tau * sum((data$y - b0 - b1*data$x1) * data$x2)) / prec # specify the mean
  rnorm(1, mean = mean, sd = 1 / sqrt(prec)) # draw a random sample from the normal dist.
}

## Big sampling function: loop through the pre-specified mini-functions for 'n.iter' times (default = 10000)
gibbs_sampler <- function(data, n.iter = 10000, init.tau, init.b0, init.b1, init.b2, alpha_tau, beta_tau, mu_b0, tau_b0, mu_b1, tau_b1, nu_b1, mu_b2, tau_b2) { # function takes initial values for each parameter and hyperparameters of each prior distribution
  tau <- numeric(n.iter); b0 <- numeric(n.iter); b1 <- numeric(n.iter); b2 <- numeric(n.iter) # storage for sampled parameters
  tau[1] <- init.tau; b0[1] <- init.b0; b1[1] <- init.b1; b2[1] <- init.b2 # plug-in the initial values
  # start the loop: start from 2nd iteration (as 1st value is already set)
  for (i in 2:n.iter) {
    tau[i] <- sample_tau(data, b0[i-1], b1[i-1], b2[i-1], alpha_tau, beta_tau)
    b0[i] <- sample_b0(data, b1[i-1], b2[i-1], tau[i], mu_b0, tau_b0)
    b1[i] <- MH_sample_b1(data, b0[i], b1[i-1], b2[i-1], tau[i], nu_b1, mu_b1, tau_b1)
    b2[i] <- sample_b2(data, b0[i], b1[i], tau[i], mu_b2, tau_b2)
  }
  # collect the resulting samples and store it in a dataframe
  sampled.pars <- data.frame(tau, b0, b1, b2)
  return(sampled.pars)
}

## Final function to implement multiple chains: use "gibbs_sampler" function "nchains" times
multiple_chains <- function(data, nchains, niter, initial_values, prior_values){
  chains <- list() # storage for each chain
  for (i in 1:nchains){ # run "gibbs_sampler" for "nchains" times
    chains[[i]] <- gibbs_sampler(data = data, n.iter = niter,
      # initial_values need to be given in a list format
      init.tau = initial_values$tau[i], init.b0 = initial_values$b0[i], init.b1 = initial_values$b1[i], init.b2 = initial_values$b2[i],
      # prior_values need to be given in a list format
      alpha_tau = prior_values$alpha_tau, # shape parameter for gamma dist.
      beta_tau = prior_values$beta_tau, # rate parameter for gamma dist.
      mu_b0 = prior_values$mu_b0, # mu for normal dist. (b0)
      tau_b0 = prior_values$tau_b0, # tau for normal dist. (b0)
      mu_b1 = prior_values$mu_b1, # mu for t-dist. (b1)
      tau_b1 = prior_values$tau_b1, # tau for t-dist. (b1)
      nu_b1 = prior_values$nu_b1, # nu (degrees of freedom) for t-dist. (b1)
      mu_b2 = prior_values$mu_b2, # mu for normal dist. (b2)
      tau_b2 = prior_values$tau_b2 # tau for normal dist. (b2)
    )
  }
  return(chains)
}

##============================================================
## 3. Perform MCMC sampling (run Gibbs sampler)
##============================================================
## Specify initial values (for two chains)
init <- list(tau = c(0.5, 0.3), b0 = c(1, 5),  b1 = c(0, -0.5), b2 = c(0, 0.3))

## Specify parameter values for prior distributions (uninformative)
prior_values <- list(alpha_tau = 0.001, beta_tau = 0.001, mu_b0 = 0, tau_b0 = 0.01,  mu_b1 = 0, 
                     tau_b1 = 0.01, nu_b1 = 3, mu_b2 = 0, tau_b2 = 0.01)

## Call the final function ("multiple_chains") and store the posterior samples
set.seed(123) # set the seed 
gibbs_results <- multiple_chains(data = BW, nchains = 2, niter = 1e5, initial_values = init, prior_values = prior_values) # run the function (two chains and number of iteration is 100000)

## Specify the burn-in = 1000
chain1 <- gibbs_results[[1]][-c(1:1000),]; chain2 <- gibbs_results[[2]][-c(1:1000),]

## Combine the samples from both chains
combined.chain <-rbind(chain1, chain2)

##============================================================
## 4. Assess convergence
##============================================================
## 1) Trace plot, Density plot, and Autocorrelation plot
for (i in 1:(ncol(combined.chain))){ # for all parameters
  layout(matrix(c(2,2,2,3,1), 1, 5, byrow=T)) # set the plotting layout
  # Density plot
  plot(density(combined.chain[,i]), col = "navy", main ="", xlab="Density", ylab="", cex.axis=2, cex.lab=2)
  # Trace plot for each chain
  plot(1:nrow(chain1), chain1[,i], type="l", col="salmon", main = "", xlab="Iteration", ylab="", cex.axis=2, cex.lab=2)
  lines(chain2[,i], col=alpha("skyblue", 0.8))
  # Autocorrelation plot
  auto.corr <- c() # storage for autocorrelation
  y <- combined.chain[,i] # extract the parameter value that the autocorr. is computed for 
  for (j in 1:50){ # compute the autocorrelation for 50 lags
    auto.corr[j] <- sum((y - mean(y)) * (dplyr::lag(x=y, j-1) - mean(y)),na.rm = T) / sum((y - mean(y))^2)
  }
  plot(auto.corr, type="h",main = "", ylab="", xlab= "Lag", ylim=c(0,1), col="navy", cex.axis=2,cex.lab=2) # plot the autocorrelation
  mtext(paste0("Diagnostics for ", colnames(combined.chain)[i]), side = 3, line = -2.5, outer = TRUE, cex = 1.7, font = 2) # give a common title for the diagnostic plots
}
## 2) Compute Gelman-Rubin statistic
nObs = nrow(chain1) # number of iteration for each chain (same for all chains)
nCols = 2  # number of chains (two in this case)
n1 = floor(nObs * 0.5) # divide the chain by two parts (first half)
n2 = nObs - n1 # second half

Gelman <- c() # storage for the resulting statistic values
for (i in 1:ncol(combined.chain)){ 
  halfchain1 = chain1[,i][-(1:n1)]  # subset the second half from each chain
  halfchain2 = chain2[,i][-(1:n1)] 
  half <- cbind(halfchain1, halfchain2) 
  vars = apply(half, 2, var) # compute the variance of each half
  means = apply(half, 2, mean) # compute the mean of each half
  mBar = mean(means) # compute the grand mean
  B = n2 * sum((means - mBar)^2)/(nCols - 1) # compute Between-chain variance
  W = sum(vars)/nCols # compute Within-chain variance
  vHat = ((n2-1) * W + B)/(n2)  + B/(n2 * nCols) # weighted sum of W and B
  Gelman[i] <- sqrt(vHat/W * (n2/(n2 - 2))) # store Gelman-rubin statistic
}
## 3) Compute MC error (Naive SE) = SD / sqrt(number of iterations)
combined.chain$resid.sd <- 1/sqrt(combined.chain$tau) # create another column of residual SD
MCerror <- apply(combined.chain, 2, function(x) sd(x)/sqrt(nrow(combined.chain)))

##==============================================================
## 5. Posterior Predictive Check (Homoscedasticity & Normality)
##==============================================================
set.seed(123) # set the seed to reproduce the result

## Step1: sample 10,000 times from the posterior distributions 
sampled.b0 <- sample(combined.chain$b0, 10000)
sampled.b1 <-sample(combined.chain$b1, 10000)
sampled.b2 <-sample(combined.chain$b2, 10000)
sampled.tau <-sample(combined.chain$tau, 10000)
sampled.par <- cbind(sampled.b0, sampled.b1, sampled.b2, sampled.tau) # all sampled parameter values

## Step2: generate Y using sampled parameters
sim.y <- matrix(NA, nrow(sampled.par), nrow(BW)) # storage for the simulated Y
for (i in 1:nrow(sampled.par)){ # for each set of parameters
  for(j in 1:nrow(BW)){ # for each row (person) in the data 
    error = rnorm(nrow(BW), 0, 1/sqrt(sampled.tau[i])) # generate a set of error 
    sim.y[i,j] <- sampled.b0[i] + sampled.b1[i]*BW$x1[j] + sampled.b2[i]*BW$x2[j] + error[j] # generate Y value
  }
}

## Step3: compute the simulated and observed residuals
sim.residuals <- matrix(NA, nrow(sampled.par), nrow(BW)) # storage for simulated residuals
obs.residuals <- matrix(NA, nrow(sampled.par), nrow(BW)) # storage for observed residuals
for (i in 1:nrow(sampled.par)){ # for each set of sampled parameters
  for(j in 1:nrow(BW)){ # for each row (person) in the data 
    # simulated residual = simulated Y - predicted value
    sim.residuals[i,j] <- sim.y[i,j] - (sampled.b0[i] + sampled.b1[i]*BW$x1[j] + sampled.b2[i]*BW$x2[j])
    # observed residual = observed Y - predicted value
    obs.residuals[i,j] <- BW$y[j] - (sampled.b0[i] + sampled.b1[i]*BW$x1[j] + sampled.b2[i]*BW$x2[j]) 
  }
}

## Step4-1: compute discrepancy measure for homoscedasticity = absolute value of difference between the variance of the first half of ordered residual and the second half of the ordered residuals
sim.diff <- numeric(nrow(sampled.par)) # storage for discrepancy measure value of simulated data
obs.diff <- numeric(nrow(sampled.par)) # storage for discrepancy measure value of observed data
p <- 0 # p-value shall be updated
for(i in 1:nrow(sampled.par)){ # for each set of simulated/observed residuals
  # order the simulated residuals by the simulated Y value
  sim.residuals[i,] <- sim.residuals[i,][order(sim.y[i,])]
  # order the observed residuals by observed Y value
  obs.residuals[i,] <- obs.residuals[i,][order(BW$y)] 
  # absolute difference between the first half and the second half of simulated residual variance
  sim.diff[i] <- abs(var(sim.residuals[i, 1:21]) - var(sim.residuals[i, 22:42]))
  # absolute difference between the first half and the second half of observed residual variance
  obs.diff[i] <- abs(var(obs.residuals[i, 1:21]) - var(obs.residuals[i, 22:42]))
  # if absolute difference from simulated residuals is larger than observed residuals, update p-value 
  if(sim.diff[i] > obs.diff[i]) p <- p + 1
}
## p-value (for homoscedasticity):
p.bayesian = p/nrow(sampled.par); p.bayesian

## Step4-2: compute discrepancy measure for normality = absolute value of difference between 0.95 and the area under the probability density curve within 1.96 SDs from its mean
sim.diff.95 <- numeric(nrow(sampled.par)) # storage for discrepancy measure value of simulated data
obs.diff.95 <- numeric(nrow(sampled.par)) # storage for discrepancy measure value of observed data
p <- 0 # p-value shall be updated
for(i in 1:nrow(sampled.par)){
  # specify empirical cumulative distribution function (CDF) for each of the simulated and observed residuals using ecdf() function
  sim.cum.density <- ecdf(sim.residuals[i,])
  obs.cum.density <- ecdf(obs.residuals[i,])
  # compute the boundary values for the simulated residuals
  sim.sd <- sd(sim.residuals[i,])
  sim.mean <- mean(sim.residuals[i,])
  sim.boundary <- c(sim.mean - 1.96*sim.sd, sim.mean + 1.96*sim.sd)
  # compute the boundary values for the observed residuals
  obs.sd <- sd(obs.residuals[i,])
  obs.mean <- mean(obs.residuals[i,])
  obs.boundary <- c(obs.mean - 1.96*obs.sd, obs.mean + 1.96*obs.sd)
  # compute the absolute value of difference between 0.95 and the area under the curve(CDF[upper.boundary] - CDF[lower boundary])
  sim.diff.95[i] <- abs(0.95 - (sim.cum.density(sim.boundary[2]) - sim.cum.density(sim.boundary[1])))
  obs.diff.95[i] <- abs(0.95 - (obs.cum.density(obs.boundary[2]) - obs.cum.density(obs.boundary[1])))
  # if absolute difference from simulated residuals is larger than observed residuals, update p-value 
  if(sim.diff.95[i] > obs.diff.95[i]) p <- p +1
}
## p-value (for normality):
p.bayesian <- p/nrow(sampled.par); p.bayesian 

##==============================================================
## 6. Parameter estimates
##==============================================================
## Compute the means
apply(combined.chain, 2, mean)
## Compute the medians
apply(combined.chain, 2, median)
## Compute the SDs
apply(combined.chain, 2, sd)
## Compute the 95% credible intervals
apply(combined.chain,2, function(x) quantile(x, probs=c(0.025, 0.975)))
## Compute the 95% Highest Density intervals (HDI)
computeHDI <- function(x, credlevel = 0.95){ # define a function that computes HDI: "computeHDI"
  sorted.x <- sort(x) # order x in an ascending order
  n <- length(sorted.x) # get the length of x
  # compute the number of values to exclude 
  # (e.g. if there are 1000 values and computing 95% hdi, then we need to exclude 1000 - 950 = 50 values)
  exclude <- n - floor(n * credlevel) 
  # compute the possible lower limits and corresponding upper limits
  # using brute-force method: we take all possible lower bound values and corresponding upper bound values
  possible.lower <- sorted.x[1:exclude]           
  possible.upper <- sorted.x[(n - exclude + 1):n]  
  # find the combination that gives the narrowest interval
  narrowest <- which.min(possible.upper - possible.lower)
  HDI <- c(possible.lower[narrowest], possible.upper[narrowest]) 
  names(HDI) <- c("lower", "upper")
  return(HDI)
}
# apply the "computeHDI" function and compute 95% HDI
apply(combined.chain, 2, function(x) computeHDI(x, 0.95))

##===================================================================
## 7. Model comparison: WAIC (Wantanabe-Akaike Information Criterion)
##===================================================================
## Models that are compared:
# M1 <- y = b0
# M2 <- y = b0 + b1*X1
# M3 <- y = b0 + b2*X2
# M4 <- y = b0 + b1*X1 + b2*X2 (the original model that I've been using)

## First, generate samples from the posterior of each model using Gibbs sampler (M4 is the original model which we sampled already above.)
## Gibbs_sampler function for M1
gibbs_sampler_M1 <- function(data, n.iter = 10000, init.tau, init.b0, alpha_tau, beta_tau, mu_b0, tau_b0) { 
  tau <- numeric(n.iter); b0 <- numeric(n.iter) # storage for sampled parameters
  tau[1] <- init.tau; b0[1] <- init.b0 # plug-in the initial values
  # start the loop: start from 2nd iteration (as 1st value is already set)
  for (i in 2:n.iter) {
    tau[i] <- sample_tau(data = data, b0 = b0[i-1],  alpha0 = alpha_tau, beta0 =beta_tau)
    b0[i] <- sample_b0(data=data, tau=tau[i], mu00=mu_b0, tau00=tau_b0)
  }
  # collect the resulting samples and store it in a dataframe
  sampled.pars <- data.frame(tau, b0)
  return(sampled.pars)
}
## Function to implement multiple chains for M1
multiple_chains_M1 <- function(data, nchains, niter, initial_values, prior_values){
  chains <- list() # storage for each chain
  for (i in 1:nchains){ # run "gibbs_sampler_M1" for the number of chains
    chains[[i]] <- gibbs_sampler_M1(data = data, n.iter = niter,
      # initial_values need to be given in a list format
      init.tau = initial_values[[1]][i], init.b0 = initial_values[[2]][i],
      # prior_values need to be given in a list format
      alpha_tau = prior_values[[1]], beta_tau = prior_values[[2]], mu_b0 = prior_values[[3]], tau_b0 = prior_values[[4]])
  }
  return(chains)
}
## Generate samples from the posterior of M1
set.seed(123) # set the seed 
gibbs_results_M1 <- multiple_chains_M1(data = BW, nchains = 2, niter = 1e5, initial_values = init, prior_values = prior_values) # use the same initial values and prior parameters as specified above.
## Combine the samples from both chains (burn-in = 1000)
combined.chain_M1 <-rbind(gibbs_results_M1[[1]][-c(1:1000),], gibbs_results_M1[[2]][-c(1:1000),])

## Gibbs_sampler function for M2
gibbs_sampler_M2 <- function(data, n.iter = 10000, init.tau, init.b0, init.b1, alpha_tau, beta_tau, mu_b0, tau_b0, mu_b1, tau_b1) { 
  tau <- numeric(n.iter); b0 <- numeric(n.iter); b1 <- numeric(n.iter) # storage for sampled parameters
  tau[1] <- init.tau; b0[1] <- init.b0; b1[1] <- init.b1 # plug-in the initial values
  # start the loop: start from 2nd iteration (as 1st value is already set)
  for (i in 2:n.iter) {
    tau[i] <- sample_tau(data =data, b0=b0[i-1], b1=b1[i-1], alpha0=alpha_tau, beta0=beta_tau)
    b0[i] <- sample_b0(data=data, b1=b1[i-1], tau=tau[i], mu00=mu_b0, tau00=tau_b0)
    b1[i] <- sample_b1(data=data, b0=b0[i],  tau=tau[i], mu10=mu_b1, tau10=tau_b1)
  }
  # collect the resulting samples and store it in a dataframe
  sampled.pars <- data.frame(tau, b0, b1)
  return(sampled.pars)
}
## Function to implement multiple chains for M2
multiple_chains_M2 <- function(data, nchains, niter, initial_values, prior_values){
  chains <- list() # storage for each chain
  for (i in 1:nchains){ # run "gibbs_sampler_M2" for the number of chains
    chains[[i]] <- gibbs_sampler_M2(data = data, n.iter = niter,
      # initial_values need to be given in a list format
      init.tau = initial_values$tau[i], init.b0 = initial_values$b0[i], init.b1 = initial_values$b1[i],
      # prior_values need to be given in a list format
      alpha_tau = prior_values$alpha_tau, beta_tau = prior_values$beta_tau, mu_b0 = prior_values$mu_b0, tau_b0 = prior_values$tau_b0, mu_b1 = prior_values$mu_b1, tau_b1 = prior_values$tau_b1)
  }
  return(chains)
}
## Generate samples from the posterior of M2
set.seed(123) # set the seed 
gibbs_results_M2 <- multiple_chains_M2(data = BW, nchains = 2, niter = 1e5, initial_values = init, prior_values = prior_values)  ## use the same initial values and prior parameters as specified above.
## Combine the samples from both chains (burn-in = 1000)
combined.chain_M2 <-rbind(gibbs_results_M2[[1]][-c(1:1000),], gibbs_results_M2[[2]][-c(1:1000),])

## Gibbs_sampler function for M3
gibbs_sampler_M3 <- function(data, n.iter = 10000, init.tau, init.b0, init.b2, alpha_tau, beta_tau, mu_b0, tau_b0, mu_b2, tau_b2) {
  tau <- numeric(n.iter); b0 <- numeric(n.iter); b2 <- numeric(n.iter) # storage for parameters
  tau[1] <- init.tau; b0[1] <- init.b0; b2[1] <- init.b2 # plug-in the initial values
  # start the loop: start from 2nd iteration (as 1st value is already set)
  for (i in 2:n.iter) {
    tau[i] <- sample_tau(data=data, b0=b0[i-1], b2=b2[i-1], alpha0=alpha_tau, beta0=beta_tau)
    b0[i] <- sample_b0(data=data, b2=b2[i-1], tau=tau[i], mu00=mu_b0, tau00=tau_b0)
    b2[i] <- sample_b2(data=data, b0=b0[i], tau=tau[i], mu20=mu_b2, tau20=tau_b2)
  }
  # collect the resulting samples and store it in a dataframe
  sampled.pars <- data.frame(tau, b0, b2)
  return(sampled.pars)
}
## Function to implement multiple chains for M3
multiple_chains_M3 <- function(data, nchains, niter, initial_values, prior_values){
  chains <- list() # storage for each chain
  for (i in 1:nchains){ # run "gibbs_sampler_M3" for the number of chains
    chains[[i]] <- gibbs_sampler_M3(data = data, n.iter = niter,
      # initial_values need to be given in a list format
      init.tau = initial_values$tau[i], init.b0 = initial_values$b0[i], init.b2 = initial_values$b2[i],
      # prior_values need to be given in a list format
      alpha_tau = prior_values$alpha_tau, beta_tau = prior_values$beta_tau, mu_b0 = prior_values$mu_b0, tau_b0 = prior_values$tau_b0, mu_b2 = prior_values$mu_b2, tau_b2 = prior_values$tau_b2)
  }
  return(chains)
}
## Generate samples from the posterior of M3
set.seed(123) # set the seed 
gibbs_results_M3 <- multiple_chains_M3(data = BW, nchains = 2, niter = 1e5, initial_values = init, prior_values = prior_values) 
## Combine the samples from both chains (burn-in = 1000)
combined.chain_M3 <-rbind(gibbs_results_M3[[1]][-c(1:1000),], gibbs_results_M3[[2]][-c(1:1000),])

## Define a function that computes WAIC: "getWAIC"
getWAIC <- function(posterior){ # take the posterior samples as input argument
  N <- nrow(BW) # number of samples in data
  M <- nrow(posterior) # number of posterior samples
  loglike <- matrix(NA, M, N) # storage for loglikelihood of each observation
  mu <- "posterior$b0" # define the mu (every model has an intercept)
  if(any(grepl("b1", colnames(posterior)))) mu <- paste(mu, "+posterior$b1 * BW$x1[i]") # if a model has "b1" parameter, then add this to the function of mu
  if(any(grepl("b2", colnames(posterior)))) mu <- paste(mu,"+posterior$b2 * BW$x2[i]") # if a model has "b2" parameter, then add this to the function of mu
  for (i in 1:N){ 
    # marginal loglikelihood for each of N observation
    loglike[,i] <- dnorm(BW$y[i], mean = eval(parse(text=mu)), sd = sqrt(1/posterior$tau), log=T)
  }
  # pWAIC (number of effective parameter) = variance of loglikelihood across M posterior samples for each observation
  pWAIC <- apply(loglike, 2, var)
  # ppd = average likelihoods across N obs
  ppd <- apply(exp(loglike), 2, mean)
  # WAIC = -2(sum of averaged log likelihoods - number of effective parameter)
  WAIC <- -2 * sum(log(ppd) - pWAIC)
  return(WAIC)
}

## Compute WAIC for each model
Model1<- getWAIC(combined.chain_M1)
Model2 <- getWAIC(combined.chain_M2)
Model3 <- getWAIC(combined.chain_M3)
Model4 <- getWAIC(combined.chain) # the model that we run above (combined.chain: previously obtained posterior samples)
model.comparison <- data.frame(Model1, Model2, Model3, Model4) # store the results in a dataframe

##===================================================================
## 8. Hypothesis evaluation: Bayes Factor
##===================================================================
## Data needs to be standardized first, so that we can compare b1 and b2
scaled_data <- as.data.frame(lapply(BW, scale))

## Generate the posterior samples given the scaled data
set.seed(123) # set the seed to reproduce the results
# Use the same initial values and parameters for the prior distributions as previously specified
gibbs_results_scaled <- multiple_chains(data = scaled_data, nchains = 2, niter = 1e5, initial_values = init, prior_values = prior_values) 
# Specify the burn-in = 1000
chain1_scaled <- gibbs_results_scaled[[1]][-c(1:1000),] 
chain2_scaled <- gibbs_results_scaled[[2]][-c(1:1000),] 
# Combine the samples from both chains
combined_chain_scaled <-rbind(chain1_scaled, chain2_scaled)

## Considered hypotheses:
# H1: b1 < 0 & b2 <0
# H2: abs(b1) > abs(b2)
# Hu: b1, b2 (unconstrained)

# Extract the parameters of interest (b1, b2)
par_interest <- combined_chain_scaled[,3:4] 
# Variance-Covariance between the parameter estimates
varcov <- var(par_interest)     
# Means of the parameter estimates
means <- colMeans(par_interest) 
# Number of samples that will be used
n_samples <- 1e6

## Compute Bayes Factor for H1: "b1 < 0 and b2 < 0"
set.seed(123) # set the seed to reproduce the results
n_constraint <- 2 # number of independent constraints in H1
# Specify the minimum training sample size = # of independent constraints / sample size
min_train_sample <- n_constraint/nrow(scaled_data) 
# Draw random samples from the posterior distribution
fit_samples <- data.frame(mvrnorm(n_samples, means, varcov)) 
# Draw random samples from the prior distribution
complexity_samples <- data.frame(mvrnorm(n_samples, rep(0, length(means)), varcov / min_train_sample))
fit <- matrix(NA, n_samples,1) # storage for fit values
complexity <- matrix(NA, n_samples,1) # storage for complexity values
for (i in 1:n_samples){ # for each sampled value,
  # Check if posterior sample is in agreement with H1
  fit[i,] <- fit_samples$b1[i] <0 & fit_samples$b2[i] <0 
  # Check if prior sample is in agreement with H1
  complexity[i,] <- complexity_samples$b1[i] <0 & complexity_samples$b2[i] <0 
}
fit_H1 <- mean(fit) # fit = proportion of posterior samples that are in accordance with H1
comp_H1 <- mean(complexity) # complexity = proportion of prior samples that are accordance with H1
BF.u_H1 <- fit_H1 / comp_H1 # Bayes factor of H1 versus the unconstrained hypothesis Hu
Bf.c_H1 <- (fit_H1 / comp_H1) / ((1-fit_H1)/(1-comp_H1)) # Bayes factor of H1 versus its complement Hc

## Compute Bayes Factor for H2: "abs(b1) > abs(b2)"
set.seed(123) # set the seed to reproduce the results
n_constraint <- 1 # number of independent constraints in H2
# Specify the minimum training sample size = # of independent constraints / sample size
min_train_sample <- n_constraint/nrow(scaled_data)
# Draw random samples from the posterior distribution
fit_samples <- data.frame(mvrnorm(n_samples, means, varcov)) 
# Draw random samples from the prior distribution
complexity_samples <- data.frame(mvrnorm(n_samples, rep(0, length(means)), varcov / min_train_sample))
fit <- matrix(NA, n_samples,1) # storage for fit values
complexity <- matrix(NA, n_samples,1) # storage for complexity values
for (i in 1:n_samples){ # for each sampled values,
  # Check if posterior sample is in agreement with H2
  fit[i,] <- abs(fit_samples$b1[i]) > abs(fit_samples$b2[i]) 
  # Check if prior sample is in agreement with H2
  complexity[i,] <- abs(complexity_samples$b1[i]) > abs(complexity_samples$b2[i])
}
fit_H2 <- mean(fit) # fit = proportion of posterior samples that are in accordance with H2
comp_H2 <- mean(complexity) # complexity = proportion of prior samples that are in accordance with H2
BF.u_H2 <- fit_H2 / comp_H2 # Bayes factor of H2 versus the unconstrained hypothesis Hu
Bf.c_H2 <- (fit_H2 / comp_H2) / ((1-fit_H2)/(1-comp_H2)) # Bayes factor of H2 versus its complement Hc

## Compute posterior model probabilities (PMP)
PMP_H1 <- (fit_H1/comp_H1) / ((fit_H1/comp_H1)+(fit_H2/comp_H2) + 1) # 1 is BF of Hu
PMP_H2 <- (fit_H2/comp_H2) / ((fit_H1/comp_H1)+(fit_H2/comp_H2) + 1) # 1 is BF of Hu
PMP_Hu <- 1 - (PMP_H1 + PMP_H2)

##===================================================================
## 9. Bayesian Causal Inference (Inverse Propensity Weighting)
##===================================================================
## Question of interest: what is the average causal effect (ACE) of maternal smoking on baby's birthweight?

## Step1: Predict mother's smoking (cause) based on covariates available in the original dataset
## (e.g., mother's age, father's smoking, father's education years) using logistic regression.
## Here I use "brm" function from "brms" package to fit Bayesian logistic regression model 
## with non-informative default priors.
model_smoking <- brm(
  bf(msmoker ~ mage + fnocig + fedyrs), # "msmoker" = binary variable (0: smoker, 1: non-smoker)
  data = dat, family = bernoulli(), iter = 1000, chains = 2,
  seed = 123) # set the seed to reproduce the results

## Step2: Compute posterior predicted propensity scores 
## using posterior_epred() function from brms package.
## This will give 1000 propensity scores for each of 47 (data sample size) individuals.
set.seed(123)
pred_probs_chains <- posterior_epred(model_smoking)
dim(pred_probs_chains) # rows are each iteration; columns are each individual

## Step3: Compute IPW using each of the propensity scores, and
## run the outcome model 1000 times to obtain ACE estimate along with its standard deviation
ACE <- matrix(NA, nrow(pred_probs_chains),2) # storage for ACE estimates
for (i in 1:nrow(pred_probs_chains)){ # for each of 1000 sets of propensity scores,
  propensity_scores <- pred_probs_chains[i,]
  Y <- dat$Birthweight # outcome variable
  X <- dat$msmoker # cause variable
  # compute the inverse propensity weight using the propensity score
  ipw <- (X/propensity_scores) +((1-X)/(1-propensity_scores)) 
  # run the outcome model using those weights
  outcome_model <- lm(Y ~ X, weights = ipw) 
  ACE[i,1] <- outcome_model$coefficient[2] # obtain the ACE estimates
  ACE[i,2] <- sqrt(diag(vcov(outcome_model)))[2] # obtain the SD of ACE estimates
}

## Step4: Combine the resulting estimates
# Density of ACE (we get a whole distribution, not a single point estimate) 
plot(density(ACE[,1]), main = "Distribution of Average Causal Effect")
# Compute mean and SD of ACE
meanACE <- colMeans(ACE)[1]; seACE <- colMeans(ACE)[2]
# Compute credible interval of ACE
quantile(ACE[,1], probs = c(0.025, 0.975))

####################################################################
## Total lines of code =                                          ##
## 531 lines - 197 lines (comments + empty space) = 334 lines :)  ##
####################################################################