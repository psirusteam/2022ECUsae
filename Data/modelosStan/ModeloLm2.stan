data {
  int<lower=0> n;   // Número de observaciones
  int<lower=0> K;   // Número de predictores
  matrix[n, K] x;   // Matrix de predictores
  vector[n] y;      // Vector respuesta
}
parameters {
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
model {
  to_vector(beta) ~ normal(250, 10000);
  sigma ~ cauchy(1, 100); 
  y ~ normal(x * beta, sigma);  // likelihood
}
generated quantities {
    real ypred[n];                    // vector de longitud n
    ypred = normal_rng(x * beta, sigma);
}

