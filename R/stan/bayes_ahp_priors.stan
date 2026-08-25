// R/stan/bayes_ahp_priors.stan
// ===========================================================================
//  Bayesian AHP -- SAME log-linear model as bayes_ahp.stan, but with the prior
//  standard deviations passed in as DATA so a prior-sensitivity sweep can vary
//  them without editing/recompiling the model.
//
//  Model:  log(a_ij) = theta_i - theta_j + e_ij ,  e_ij ~ Normal(0, tau^2)
//          theta_1 = 0 (reference); w = softmax(theta).
//  Priors: theta_k ~ Normal(0, theta_prior_sd)   (k = 2..n)
//          tau     ~ Normal(0, tau_prior_sd) T[0,]   (half-normal)
//
//  Setting theta_prior_sd = 3 and tau_prior_sd = 1 reproduces bayes_ahp.stan
//  EXACTLY, so the baseline run is recoverable as one row of the sweep.
// ===========================================================================
data {
  int<lower=2> n;
  int<lower=1> J;
  int<lower=1, upper=n> ii[J];
  int<lower=1, upper=n> jj[J];
  vector[J] y;
  real<lower=0> tau_prior_sd;     // swept (e.g. 0.5, 1, 2)
  real<lower=0> theta_prior_sd;   // usually fixed at 3; can also be swept
}
parameters {
  vector[n-1] theta_free;
  real<lower=0> tau;
}
transformed parameters {
  vector[n] theta;
  theta[1] = 0;
  theta[2:n] = theta_free;
}
model {
  theta_free ~ normal(0, theta_prior_sd);
  tau ~ normal(0, tau_prior_sd);   // half-normal via <lower=0>
  for (k in 1:J)
    y[k] ~ normal(theta[ii[k]] - theta[jj[k]], tau);
}
generated quantities {
  simplex[n] w = softmax(theta);
}
