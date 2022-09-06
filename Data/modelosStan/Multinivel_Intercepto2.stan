data {
  int <lower = 0>N;                       // Número de observaciones
  int <lower = 0>k;                       // Número de covariables 
  vector[N] y;                          // Variables respuesta
  matrix [N,k] X;                        // Variables regresoras
  // efecto aleatorio 
  int <lower = 0>kz; 
  matrix [N,kz] Z;                        // Variables regresoras
}

parameters {
  vector[k] beta;         // coeficientes del modelo
  vector[kz] u;         // coeficientes del modelo
  real <lower = 0> sigma;
}

transformed parameters{
 vector[N] mu;
 mu =  X*beta + Z * u;
}

model {
  beta ~ normal(0,100);
  u ~ normal(0,1000);
  sigma ~ cauchy(0,100);
  y ~ normal(mu , sigma);
}
