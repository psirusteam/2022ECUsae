#########################################################
# Proyecto IPM                   #
# Mapas de pobreza CEPAL                                #
# Lectura y preparación de las bases de datos           #
# Autor: Stalyn Guerrero & Andrés Gutiérrez             #
#########################################################
rm(list =ls())

# Loading required libraries ----------------------------------------------

library(patchwork)
library(trafo)
library(normtest)
library(nortest)
library(lme4)
library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
source("Frecuentista_depto/0Funciones/funciones_mrp.R", encoding = "UTF-8")

# Loading data ------------------------------------------------------------
memory.limit(10000000000000)
fit_educacion <- readRDS( file = "Frecuentista_depto_b/COL/Data/fit_freq_educacion.rds")
fit_empleo <- readRDS( file = "Frecuentista_depto_b/COL/Data/fit_freq_empleo.rds")

encuesta_ipm <- readRDS("Frecuentista_depto_b/COL/Data/encuesta_ipm.rds")
censo_ipm <- readRDS("Frecuentista_depto_b/COL/Data/censo_ipm2.rds") 
tasa_desocupados <- readRDS("Frecuentista_depto_b/COL/Data/tasa_desocupacion.rds")

# Agregando encuesta ------------------------------------------------------
statelevel_predictors_df <- tasa_desocupados
byAgrega <- c("depto", "area", "sexo", "edad",
              "ipm_Material",
              "ipm_Hacinamiento",
              "ipm_Agua",
              "ipm_Saneamiento",
              "ipm_Energia",
              "ipm_Internet" )

poststrat_df <- inner_join(censo_ipm %>% mutate_at(all_of(byAgrega), as.character),
                           statelevel_predictors_df)


poststrat_df$prob_educacion <- predict(
  fit_educacion,
  newdata = poststrat_df,
  type = "response",
  allow.new.levels = TRUE
)
poststrat_df$prob_empleo <- predict(
  fit_empleo,
  newdata = poststrat_df,
  type = "response",
  allow.new.levels = TRUE
)

poststrat_MC <- poststrat_df %>% data.frame() %>% 
  dplyr::select(depto, n,matches("prob|ipm_") ) %>% 
  tibble()

poststrat_MC %<>% mutate_at(vars(matches("ipm_")), as.numeric) 

fipm <- function(){
iter_ipm2 <-  pmap(as.list(poststrat_MC),
       function(
       depto,n,
       ipm_Material,
       ipm_Saneamiento,
       ipm_Energia,
       ipm_Internet,
       ipm_Agua,
       ipm_Hacinamiento,
       prob_educacion,
       prob_empleo
        ) {
         y_empleo =  rbinom(n, 1, prob = prob_empleo)
         y_educacion =  rbinom(n, 1, prob = prob_educacion)
         ipm <- 0.1 * (
           ipm_Material+
           ipm_Saneamiento+
           ipm_Energia+
           ipm_Internet+
           ipm_Agua+
           ipm_Hacinamiento) +
  0.2 * (y_educacion +   y_empleo) 

ipm_dummy <- ifelse(ipm < 0.4, 0, 1) 
mean(ipm_dummy)
                          })

poststrat_MC$ipm <- unlist(iter_ipm2) 

poststrat_MC %>% group_by(depto) %>%
summarise(ipm = sum((n*ipm))/sum(n))
}

ipm_MC <- replicate(100,fipm())

ipm_estimado <- map_df(1:ncol(ipm_MC), function(x) data.frame(t(ipm_MC)[x,])) %>%
group_by(depto) %>% summarise(ipm_estimado_MC = mean(ipm))

saveRDS(ipm_estimado, 
    file = "Frecuentista_depto_b/COL/Data/ipm_MC.rds")


