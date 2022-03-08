
data {
  int<lower=0> N;
  //int<lower=0> L;
  vector[N] X;
  //vector[N] X_pred;
  vector[N] y;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  mu = X*beta + alpha;
}

model {
  y ~ normal(mu, sigma);
}

generated quantities {
  real Rsq;
  vector[N] res;
  vector[N] error;
  
  for (i in 1:N) {
    res[i] = (mu[i]- mean(y))^2;
    error[i] = (y[i] - mu[i])^2;
  }

  Rsq = (sum(res)/(N-1))/((sum(res)/(N-1)) + (sum(error)/(N-1)));
}