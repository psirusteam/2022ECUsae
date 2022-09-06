data {
  int <lower = 0>N;                       // Número de observaciones
  int <lower = 0>k;                       // Número de covariables 
  int  Nv1;                               // Número de grupos
  int  Grupo[N];                         // Ids de grupos
  vector[N] y;                          // Variables respuesta
  matrix [N,k] X;                        // Variables regresoras
}

parameters {
  vector[Nv1] beta0;       // coeficientes del modelo
  vector[k] beta1;         // coeficientes del modelo
  real <lower = 0> sigma;
}

transformed parameters{
 vector[N] mu;
 mu =  beta0[Grupo] + X*beta1;
}

model {
  beta0 ~ normal(0,100);
  beta1 ~ normal(0,100);
  sigma ~ cauchy(0,100);
  y ~ normal(mu , sigma);
  
}
