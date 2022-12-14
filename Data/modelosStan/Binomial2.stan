data {
  int<lower=0> K;                 // Número de provincia  
  int<lower=0> n[K];              // Número de ensayos 
  int<lower=0> s[K];              // Número de éxitos
}
parameters {
  real<lower=0, upper=1> theta; // theta|s1,s2,...sK
}
model {
  for(kk in 1:K) {
  s[kk] ~ binomial(n[kk], theta);
  }
  theta ~ beta(0.5, 0.5);
}

generated quantities {
    real spred[K];                    // vector de longitud D
    for(kk in 1:K){
    spred[kk] = binomial_rng(n[kk],theta);
}
}
