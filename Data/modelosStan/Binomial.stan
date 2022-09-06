data {
  int<lower=0> n;              // Número de ensayos 
  int<lower=0> s;              // Número de éxitos
}
parameters {
  real<lower=0, upper=1> theta;     // theta|s1,s2..sD
}
model {

  s ~ binomial(n, theta);

  theta ~ beta(1600, 5000);

}

generated quantities {
    real spred;                    // vector de longitud D
    spred = binomial_rng(n, theta);

}
