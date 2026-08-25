// R/stan/bayes_ahp.stan
// ===========================================================================
//  Bayesian AHP -- log-linear (Crawford-Williams) model for RQ3, Part 1.
//  Same model as the hand-written Metropolis-Hastings fallback in
//  11_bayesian_ahp.R, expressed for Stan's NUTS sampler.
//
//  Model:  log(a_ij) = theta_i - theta_j + e_ij ,  e_ij ~ Normal(0, tau^2)
//          theta = log-weights ; theta_1 fixed at 0 for identifiability.
//          Final weights w = softmax(theta).
//  Priors: theta_k ~ Normal(0, 3)   (k = 2..n, weak)
//          tau     ~ Normal(0, 1) T[0,]  (half-normal; data estimate the noise)
//
//  Lineage: Gargallo, Moreno-Jimenez & Salvador (2007). AHP judgements are
//  Saaty 1-9 RATIOS, so this log-normal model is the correct estimator -- not
//  the Bradley-Terry / bpcs family (which is for win/loss data).
// ===========================================================================
data {
  int<lower=2> n;              // number of factors (5)
  int<lower=1> J;              // number of upper-triangle judgements (10)
  int<lower=1, upper=n> ii[J]; // row index i (i < j)
  int<lower=1, upper=n> jj[J]; // col index j
  vector[J] y;                 // y[k] = log(a_ij)
}
parameters {
  vector[n-1] theta_free;      // theta_2..theta_n ; theta_1 = 0
  real<lower=0> tau;           // judgement-noise sd (estimated, not fixed)
}
transformed parameters {
  vector[n] theta;
  theta[1] = 0;                // reference factor fixes identifiability
  theta[2:n] = theta_free;
}
model {
  theta_free ~ normal(0, 3);   // weak prior on log-weights
  tau ~ normal(0, 1);          // half-normal via <lower=0>
  for (k in 1:J)
    y[k] ~ normal(theta[ii[k]] - theta[jj[k]], tau);
}
generated quantities {
  simplex[n] w = softmax(theta);   // posterior weight draw, sums to 1
}
