data {
  int<lower=0> n;   // Número de observaciones
  int<lower=0> K;   // Número de predictores
  matrix[n, K] x;   // Matrix de predictores
  int<lower=0,upper=1> y[n];      // Vector respuesta
}
parameters {
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
transformed parameters {
    vector[n] inv_eta;
   inv_eta = inv_logit(x * beta);
}

model {
  to_vector(beta) ~ normal(0, 10000);
  y ~ bernoulli(inv_eta);  // likelihood
}
generated quantities {
    real ypred[n];                    // vector de longitud n
    ypred = bernoulli_rng(inv_eta);
}

