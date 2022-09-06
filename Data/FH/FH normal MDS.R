rm(list = ls())

library(cmdstanr)
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(emdi)
library(posterior)
library(patchwork)

rstan_options(auto_write = TRUE)
options(width = 90)
set.seed(2018)

# data and simulated model ----------------------------------------------------------

data <- readRDS("Data/FH/base_sae_2009.rds")
data$vardir[data$vardir == 0] <- 0.00001
data_dir <- data %>% filter(!is.na(nd))
data_syn <- data %>% filter(is.na(nd))

thetahat <- data_dir$pobreza
vhat.dir <- data_dir$vardir

X <- data_dir %>%
  dplyr::select(log_rem_mintrab,
         porcentaje_depoblacinrural,
         porcentajedeasistenciaescolarcom,
         analfabeto,
         pad2009,
         pnna2009,
         region1,
         region2,
         region3,
         region4,
         region5,
         region6,
         region7,
         region8,
         region9,
         region10,
         region11,
         region12,
         region13,
         region14,
         region15) 

Xs <- data_syn %>%
  dplyr::select(log_rem_mintrab,
         porcentaje_depoblacinrural,
         porcentajedeasistenciaescolarcom,
         analfabeto,
         pad2009,
         pnna2009,
         region1,
         region2,
         region3,
         region4,
         region5,
         region6,
         region7,
         region8,
         region9,
         region10,
         region11,
         region12,
         region13,
         region14,
         region15) 

N1 <- nrow(data_dir)
N2 <- nrow(data_syn)

normal_data <- list(N1 = N1,
                    N2 = N2,
                    p = ncol(X), 
                    X = X, 
                    y = thetahat, 
                    sigma_e = sqrt(vhat.dir),
                    Xs = Xs)

# STAN fit ----------------------------------------------------------------

#' # Draw from posterior distribution
#+ results='hide'
fit_FH_Nornal <- cmdstan_model("Data/FH/FH_normal.stan")

model_FH_Nornal <-
  fit_FH_Nornal$sample(
    data = normal_data,
    chains = 4,
    parallel_chains = 4,
    seed = 1234,
    refresh = 1000
    
  )

model_FH_Nornal$summary(
  variables = c("sigma2_v", "beta"))

(mcmc_dens_chains(model_FH_Nornal$draws("sigma2_v")) +
    mcmc_areas(model_FH_Nornal$draws("sigma2_v")))/ 
  mcmc_trace(model_FH_Nornal$draws("sigma2_v")) 

mcmc_areas(model_FH_Nornal$draws(c("beta[1]","beta[2]","beta[3]","beta[4]")))+ 
mcmc_areas(model_FH_Nornal$draws(c("beta[5]","beta[6]","beta[7]","beta[8]"))) 

mcmc_trace(model_FH_Nornal$draws(c("beta[1]","beta[2]","beta[3]","beta[4]")))+ 
mcmc_trace(model_FH_Nornal$draws(c("beta[5]","beta[6]","beta[7]","beta[8]"))) 


N1 <- nrow(data_dir)
y_pred_B <- model_FH_Nornal$draws(variables = "theta", format = "matrix")
rowsrandom <- sample(nrow(y_pred_B), 100)
y_pred2 <- y_pred_B[rowsrandom, 1:N1]
ppc_dens_overlay(y = as.numeric(thetahat), y_pred2)


theta_FH <- model_FH_Nornal$summary(variables =  "theta")
plot(theta_FH$mean, thetahat)
abline(b=1,a=0, col = "red")


# EMDI estimations --------------------------------------------------------

fh_normal <- fh(fixed = pobreza ~ log_rem_mintrab + prom_pob + porcentaje_depoblacinrural + 
                  porcentajedeasistenciaescolarcom + analfabeto + pad2009 + 
                  pnna2009 + region1 + region2 + region3 + region4 + region5 + 
                  region6 + region7 + region8 + region9 + region10 + region11 + 
                  region12 + region13 + region14 + region15 - 1, 
                vardir = "vardir",
                combined_data = data,
                domains = "com_id", 
                method = "reml",
                transformation = "no",
                MSE = TRUE, mse_type = "analytical")

# Summarising estimations (Bayes + Emdi) ----------------------------------

#--- sigma_u ------
model_FH_Nornal$summary(
  variables = c("sigma2_v"))
summary(fh_normal) #0.001123079
#--- theta ------

estimacionesFH <- estimators(fh_normal,
                             indicator = "All",
                             MSE = TRUE, CV = TRUE) %>%
  as.data.frame()

fitsummary <- model_FH_Nornal$summary(variables =  "theta")
#fitsummary[, 1]

estimaciones <- data.frame(comuna = data_dir$com_str,
                           directo = data_dir$pobreza,
                           FH = estimacionesFH$FH[!is.na(estimacionesFH$Direct)],
                           Bayes = fitsummary[, 2],
                           directo_rmse = data_dir$vardir,
                           FH_rmse = sqrt(estimacionesFH$FH_MSE[!is.na(estimacionesFH$Direct)]),
                           bayes_rmse = fitsummary[, 4])


fitsummary <- summary(fit, pars = "y_pred")$summary

predicciones <- data.frame(comuna = data_syn$com_str,
                           directo = data_syn$pobreza,
                           FH = estimacionesFH$FH[is.na(estimacionesFH$Direct)],
                           Bayes = fitsummary[, 1],
                           directo_rmse = data_syn$vardir,
                           FH_rmse = sqrt(estimacionesFH$FH_MSE[is.na(estimacionesFH$Direct)]),
                           bayes_rmse = fitsummary[, 3])


EstimacionesBayes <- rbind.data.frame(estimaciones, 
                                      predicciones)


plot(EstimacionesBayes$FH, EstimacionesBayes$Bayes)

# Plotting the MCMC -------------------------------------------------------
posterior <- as.array(fit)
dim(posterior)
dimnames(posterior)

idbeta <- str_detect(dimnames(posterior)$parameters, "beta")
betapars <- dimnames(posterior)$parameters[idbeta]
idtheta <- str_detect(dimnames(posterior)$parameters, "theta")
thetapars <- dimnames(posterior)$parameters[idtheta]
idpred <- str_detect(dimnames(posterior)$parameters, "y_pred")
ypredpars <- dimnames(posterior)$parameters[idpred]

color_scheme_set("red")
plot(fit, pars = betapars)
plot(fit, pars = thetapars)
plot(fit, pars = ypredpars)

mcmc_intervals(posterior, pars = betapars)
mcmc_intervals(posterior, pars = thetapars)
mcmc_intervals(posterior, pars = ypredpars)

mcmc_areas(
  posterior, 
  pars = betapars,
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

mcmc_areas(
  posterior, 
  pars = ypredpars,
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

color_scheme_set("green")
mcmc_hist(posterior, pars = betapars)
mcmc_hist(posterior, pars = ypredpars)

color_scheme_set("brightblue")
mcmc_hist_by_chain(posterior, pars = betapars)

color_scheme_set("purple")
mcmc_dens(posterior, pars = betapars)
mcmc_dens_overlay(posterior, pars = ypredpars)

color_scheme_set("teal")
mcmc_violin(posterior, pars = thetapars, probs = c(0.1, 0.5, 0.9))

color_scheme_set("gray")
mcmc_scatter(posterior, pars = c("beta0", "beta1"), 
             size = 1.5, alpha = 0.5)

# requires hexbin package
if (requireNamespace("hexbin", quietly = TRUE)) {
  mcmc_hex(posterior, pars = c("beta0", "beta1"))
}


color_scheme_set("pink")
mcmc_pairs(posterior, pars = thetapars,
           off_diag_args = list(size = 1.5))

color_scheme_set("blue")
mcmc_trace(posterior, pars = thetapars)

color_scheme_set("mix-blue-red")
mcmc_trace(posterior, pars = thetapars, 
           facet_args = list(ncol = 1, strip.position = "left"))

mcmc_trace_highlight(posterior, pars = thetapars, highlight = 3)


# Visual MCMC diagnostics -------------------------------------------------

color_scheme_set("darkgray")
np <- nuts_params(fit)
mcmc_pairs(fit, pars = thetapars)

color_scheme_set("red")
mcmc_nuts_energy(np)

rhats <- rhat(fit)
print(rhats)

color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

mcmc_acf(posterior, pars = thetapars, lags = 10)


# Graphical posterior predictive checks -----------------------------------

sims <- as.data.frame(fit)
yrep <- as.matrix(sims[, c(4: 33)])
rowsrandom <- sample(nrow(yrep), 20)
yrep2 <- yrep[rowsrandom, ]

color_scheme_set("brightblue")

ppc_dens_overlay(y, yrep2)
ppc_hist(y, yrep2)
ppc_ecdf_overlay(y, yrep2)

prop_gzero <- function(x) mean(x >= 0)
prop_gzero(y) # check proportion of values greater tha zero in y
ppc_stat(y, yrep, stat = "prop_gzero")
ppc_stat(y, yrep, stat = "max")
ppc_stat(y, yrep, stat = "min")


# Andres' own and inneficient graphs --------------------------------------

#' ## Plot the conditional mean function
sims <- as.data.frame(fit)
plot(x, y, las = 1, pch = 20)

for (s in 1:nrow(fit)) {
  with(sims, curve(a[s] + b[s] * x, lwd = 0.5, add = TRUE, col = "red"))
}

summary_fit <- get_posterior_mean(fit)
a.hat <- summary_fit[1,1]
b.hat <- summary_fit[2,1]
curve(a.hat + b.hat * x, add = TRUE, col = "green")
curve(a + b * x, add = TRUE, col = "blue")

sims1 <- sims %>% 
  select(-c(a, b, sigma, lp__)) 

colnames(sims1) <- as.character(1:N)

sims2 <- sims1  %>%
  pivot_longer(everything(),
               names_to = "iter",
               values_to = "y_pred") %>%
  mutate(k = sort(rep(1:nrow(sims), N)))

y.orig <- data.frame(
  y_pred = y,
  k = 0,
  i = 1:N
)

g1 <- sims2 %>% ggplot() +
  geom_density(
    aes(x = y_pred, 
        group = as.factor(k)),
    col = "grey"
  ) + 
  theme_classic() +
  theme(legend.position = "none")

g1 + geom_density(
  data = y.orig,
  aes(x = y_pred)
)

g2 <- sims2 %>% ggplot() +
  geom_boxplot(
    outlier.shape = NA,
    aes(y = y_pred, 
        x = iter)) + 
  theme_classic() +
  theme(legend.position = "none")

g2 + geom_point(
  data = y.orig,
  aes(x = reorder(i, i),
      y = y_pred),
  col = 2,
  size = 4
)

launch_shinystan(fit)



