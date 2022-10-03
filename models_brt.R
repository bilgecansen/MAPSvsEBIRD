# Boosted Regression Trees on selected species

library(dismo)
library(gbm)
library(doSNOW)
library(foreach)
library(tictoc)
library(tidyverse)

data_sdm <- readRDS("data_sdm_pca.rds")
chdata <- readRDS("results_cjspop.rds")$chdata


# Boosted Regression Trees ------------------------------------------------

cl <- makeCluster(5, types = "SOCK")
registerDoSNOW(cl)

pb <- txtProgressBar(max = 17, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic()

results <- foreach(i=1:17, .packages = "dismo", .options.snow = opts) %dopar% {
  
  rate <- 0.001
  
  res <- gbm.step(data = data_sdm[[i]], 
                  gbm.x = 1:(ncol(data_sdm[[i]])-1), 
                  gbm.y = ncol(data_sdm[[i]]), 
                  tree.complexity = 3, 
                  family = "bernoulli", 
                  learning.rate = rate, 
                  bag.fraction = 0.75, 
                  step.size = 250, 
                  n.trees = 1000,
                  max.trees = 15000)
  
  return(res)
          
}

toc()

stopCluster(cl)

names(results) <- names(data_sdm)

saveRDS(results, file = "results_brt.rds")

  