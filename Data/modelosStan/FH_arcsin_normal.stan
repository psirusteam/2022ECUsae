data {
  int<lower=0> N1;   // number of data items
  int<lower=0> N2;   // number of data items for prediction
  int<lower=0> p;   // number of predictors
  matrix[N1, p] X;   // predictor matrix
  matrix[N2, p] Xs;   // predictor matrix
  vector[N1] y;      // predictor matrix 
  vector[N1] sigma_e; // known variances
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[p] beta;       // coefficients for predictors
  real<lower=0> sigma2_v;
  vector[N1] v;
}

transformed parameters{
  vector[N1] theta;
  vector[N1] lp;
  real<lower=0> sigma_v;
  lp = X * beta + v;
  sigma_v = sqrt(sigma2_v);
  for(k in 1:N1){
    theta[k] = pow(sin(lp[k]), 2);
  }
}

model {
  // likelihood
  y ~ normal(lp, sigma_e); 
  // priors
  beta ~ normal(0, 100);
  v ~ normal(0, sigma_v);
  sigma2_v ~ inv_gamma(0.0001, 0.0001);
}

generated quantities{
  vector[N2] theta_pred;
  vector[N2] lppred;
  for(j in 1:N2) {
    lppred[j] = normal_rng(Xs[j] * beta, sigma_v);
    theta_pred[j] = pow(sin(lppred[j]), 2);
  }
}
