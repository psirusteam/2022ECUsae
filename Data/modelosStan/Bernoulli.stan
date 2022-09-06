data {                         // entrada el modelo 
  int<lower=0> n;              // numero de observaciones  
  int y[n];                    // vector de longitud n
}
parameters {                   // definir parámetro
  real<lower=0, upper=1> theta;
}
model {                        // definir modelo
  y ~ bernoulli(theta);
  theta ~ beta(350, 650);      // distribución previa 
}
generated quantities {
    real ypred[n];                    // vector de longitud n
    for (ii in 1:n){
    ypred[ii] = bernoulli_rng(theta);
    }
}
