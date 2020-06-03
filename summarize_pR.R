
# Load packages and data --------------------------------------------------

library(tidyverse)
library(MCMCvis)
library(boot)
library(doSNOW)
library(foreach)
library(abind)

results_cjspop <- readRDS("results_cjspop.rds")

mcmc_chains <- results_cjspop$mcmc_chains
chdata <- results_cjspop$chdata
spcode <- results_cjspop$spcode

# Functions ---------------------------------------------------------------

prep.wdata <- function(chdata) {
  
  wdata <- map(chdata$weather_pca, function(x) array(x, dim = c(chdata$npop, chdata$nyear, 1)))
  wdata <- do.call(abind, list(wdata, along = 3))
  
  return(wdata)
  
}

estimate.pR <- function(chains, wdata, time, pop, iter) {
  
  alpha1 <- logit(chains[,"survival_juv"])
  alpha2 <- logit(chains[,"survival_ad"])
  beta <- chains[,"beta"] 
  theta <- log(chains[,"fecundity"]) 
  zeta <- chains[,"zeta"]
  sigma_s <- chains[,"sigma_s"]
  sigma_f <- chains[,"sigma_f"]
  
  eff_s <- array(rnorm(time*iter, 0, sigma_s), dim = c(time, iter))
  eff_f <- array(rnorm(time*iter, 0, sigma_f), dim = c(time, iter))
  
  index1 <- which(grepl("beta_l", colnames(chains)))
  index2 <- which(grepl("zeta_l", colnames(chains)))
  swchains <- chains[,index1]
  fwchains <- chains[,index2]
  
  R <- array(dim = c(pop, time, nrow(chains), iter))
  
  a <- 0
  for (k in 1:pop) {
    
    for (t in 1:time) {
      
      for (i in 1:iter) {
        
        s1 <- inv.logit(alpha1 + beta*(-1) + swchains %*% wdata[k,t,] + eff_s[t,i])
        s2 <- inv.logit(alpha2 + beta*(-1) + swchains %*% wdata[k,t,] + eff_s[t,i])
        fec <-  exp(theta + zeta*(-1) + fwchains %*% wdata[k,t,] + eff_f[t,i])
        
        R[k,t,,i] <- s2 + s1*fec 
        
      }
      
    }#t
    
  }#k
  
  f1 <- function(R) {
    
    avgR <- apply(R, c(2,3), function(x) prod(x)^(1/length(x)))
    pR <- length(which(avgR>1))/length(avgR)
    
    return(pR)
    
  }
  
  pR <- apply(R, 1, f1)
  
  
  return(pR)
  
}


# Estimate P(R>1) ---------------------------------------------------------

wdata <- map(chdata, prep.wdata)

cl <- makeCluster(4, type = "SOCK")
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(spcode), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

pR <- foreach (i=1:length(mcmc_chains), .options.snow = opts, 
               .packages = c("abind", "boot")) %dopar% {
  
  estimate.pR(chains = mcmc_chains[[i]],
              wdata = wdata[[i]],
              time = chdata[[i]]$nyear,
              pop = chdata[[i]]$npop,
              iter = 10)
  
}

names(pR) <- spcode

stopCluster(cl)

saveRDS(pR, "results_pR.rds")


