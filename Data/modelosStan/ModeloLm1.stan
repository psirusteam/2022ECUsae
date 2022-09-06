data {
  int<lower=0> n;   // Número de observaciones
  vector[n] x;      // Variable predictora
  vector[n] y;      // Variable respuesta
}
parameters {
  real b0;            // Intercepto
  real b1;            // Pendiente
  real<lower=0> sigma;  // error scale
}
model {
  b0 ~ normal(250, 100);
  b1 ~ normal(-1000, 100);
  sigma ~ cauchy(60, 15); 
  y ~ normal(b0 + b1*x, sigma);  // likelihood
}
generated quantities {
    real ypred[n];                    // vector de longitud n
    ypred = normal_rng(b0 + b1*x, sigma);

}
