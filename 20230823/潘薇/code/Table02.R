library(tidyverse)
library(lavaan)
library(lavaanExtra)
library(semTools)
rawdf <- readr::read_rds("./data/rawdf.rds")

model <- '
  # latent variables
  DigPC =~ DigPC1 + DigPC2 + DigPC3 + DigPC4 + DigPC5 + DigPC6 + DigPC7
'

fit_cfa1 <- cfa(model = model, data = rawdf)

model <- '
   # latent variables
   DigPC =~ DigPC1 + DigPC2 + DigPC3 + DigPC4 + DigPC5 + DigPC6 + DigPC7
   ImpC  =~ ImpC1 + ImpC2 + ImpC3 + ImpC4 + ImpC5 + ImpC6 + ImpC7 + ImpC8 + ImpC9 + ImpC10 + ImpC11 + ImpC12
'

fit_cfa2 <- cfa(model = model, data = rawdf)

model <- '
   # latent variables
   DigPC =~ DigPC1 + DigPC2 + DigPC3 + DigPC4 + DigPC5 + DigPC6 + DigPC7
   ImpC  =~ ImpC1 + ImpC2 + ImpC3 + ImpC4 + ImpC5 + ImpC6 + ImpC7 + ImpC8 + ImpC9 + ImpC10 + ImpC11 + ImpC12
   OrgR  =~ OrgR1 + OrgR2 + OrgR3 + OrgR4
'

fit_cfa3 <- cfa(model = model, data = rawdf)

model <- '
   # latent variables
   DigPC =~ DigPC1 + DigPC2 + DigPC3 + DigPC4 + DigPC5 + DigPC6 + DigPC7 
   ImpC  =~ ImpC1 + ImpC2 + ImpC3 + ImpC4 + ImpC5 + ImpC6 + ImpC7 + ImpC8 + ImpC9 + ImpC10 + ImpC11 + ImpC12
   OrgR  =~ OrgR1 + OrgR2 + OrgR3 + OrgR4
   InnP  =~ InnP1 + InnP2 + InnP3 + InnP4 + InnP5
'

fit_cfa4 <- cfa(model = model, data = rawdf)

table02 <- list(fit_cfa4, fit_cfa3, fit_cfa2, fit_cfa1) %>% 
  nice_fit(nice_table = TRUE) |> 
  flextable::fontsize(size = 9, part = "all") 

