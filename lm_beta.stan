
data {
  int<lower=0> N;
  //int<lower=0> L;
  vector[N] X;
  //vector[N] X_pred;
  vector[N] y;
}

parameters {
  real alpha;
  real<lower=0> phi;
  real beta;
}

transformed parameters {
  vector[N] mu;
  vector[N] a;
  vector[N] b;
  mu = inv_logit(X*beta + alpha);
  
  for (i in 1:N) {
    a[i] = mu[i]*phi;
    b[i] = (1-mu[i])*phi;
  }
}

model {
  for (i in 1:N) {
    y[i] ~ beta(a[i], b[i]);
  }
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