data {
  int<lower=1> N1;                    // sample size
  int<lower=1> N2;                  // sample size
  int<lower=1> p;                     // p predictors
  vector<lower=0,upper=1>[N1] y;      // response 
  matrix[N1,p] X;
  matrix[N2,p] Xs;
  vector<lower=0>[N1] phi;                   // dispersion parameter
}

parameters {
  vector[p] beta;
  real<lower=0> sigma2_v;                      // K predictors
  vector[N1] v;
// reg coefficients
}

transformed parameters{
  vector[N1] LP;
  real<lower=0> sigma_v;
  vector[N1] theta;                        // linear predictor
  LP = X * beta + v;
  sigma_v = sqrt(sigma2_v); 
  for (i in 1:N1) { 
    theta[i] = inv_logit(LP[i]); 
  }
}

model {
  // model calculations
  vector[N1] a;                         // parameter for beta distn
  vector[N1] b;                         // parameter for beta distn

  for (i in 1:N1) { 
    a[i] = theta[i] * phi[i];
    b[i] = (1 - theta[i]) * phi[i];
  }

  // priors
  beta ~ normal(0, 100);
  v ~ normal(0, sigma_v);
  sigma2_v ~ inv_gamma(0.0001, 0.0001);
  //phi ~ cauchy(0, 5);                  // different options for phi  
  //phi ~ inv_gamma(.001, .001);
  //phi ~ uniform(0, 500);             // put upper on phi if using this

  // likelihood
  y ~ beta(a, b);
}

generated quantities {
  vector[N2] y_pred;
  vector[N2] thetapred;

  for (i in 1:N2) {
    y_pred[i] = normal_rng(Xs[i] * beta, sigma_v);
    thetapred[i] = inv_logit(y_pred[i]);
  }
}
  




// data {
//   int<lower=0> N1;   // number of data items
//   int<lower=0> N2;   // number of data items for prediction
//   int<lower=0> p;   // number of predictors
//   matrix[N1, p] X;   // predictor matrix
//   matrix[N2, p] Xs;   // predictor matrix
//   vector[N1] y;      // predictor matrix
//   vector[N1] phi;    // known variances
// }
// 
// // The parameters accepted by the model. Our model
// // accepts two parameters 'mu' and 'sigma'.
// parameters {
//   vector[p] beta;       // coefficients for predictors
//   real<lower=0> sigma2_v;
//   vector[N1] v;
// }
// 
// transformed parameters{
//   vector<lower=0>[N1] a;
//   vector<lower=0>[N1] b;
//   real<lower=0> sigma_v;
//   vector[N1] theta;
//   vector[N1] f;
// 
//   for(i in 1:N1) {
//     a[i] = theta[i] * phi[i];
//     b[i] = (1 - theta[i]) * phi[i];
//   }
//   sigma_v = sqrt(sigma2_v);
//   theta = inv_logit(f);
// }
// 
// model {
//   f ~ normal(X * beta, sigma_v);
//   // likelihood
//   for(i in 1:N1){
//     y[i] ~ beta(a[i], b[i]);
//   }
//   
//   // priors
//   beta ~ normal(0, 100);
//   sigma2_v ~ inv_gamma(0.0001, 0.0001);
// }
// 
// generated quantities{
//   vector[N2] theta_pred;
//   vector[N2] mupred;
//   for(j in 1:N2) {
//     theta_pred[j] = normal_rng(Xs[j] * beta, sigma_v);
//     mupred[j] = inv_logit(theta_pred[j]);
//   }
// }
