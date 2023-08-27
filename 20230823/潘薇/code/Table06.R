library(tidyverse)
library(lavaan)
library(semTools)
rawdf <- readr::read_rds("./data/rawdf.rds")

model <- '
  DigPC =~ DigPC1 + DigPC2 + DigPC3 + DigPC4 + DigPC5 + DigPC6 + DigPC7
  ImpC  =~ ImpC1 + ImpC2 + ImpC3 + ImpC4 + ImpC5 + ImpC6 + ImpC7 + ImpC8 + ImpC9 + ImpC10 + ImpC11 + ImpC12
  OrgR  =~ OrgR1 + OrgR2 + OrgR3 + OrgR4
  InnP  =~ InnP1 + InnP2 + InnP3 + InnP4 + InnP5
  
  OrgR ~ a1 * DigPC + a2 * ImpC
  InnP ~ b * OrgR + cprime1 * DigPC + cprime2 * ImpC
  
  c1         := cprime1 + a1 * b
  c2         := cprime2 + a2 * b
'

fit_med <- sem(model, data = rawdf)

table06 <- fit_med |> 
  parameterEstimates(standardized = TRUE) |> 
  filter(op == ":=") |> 
  filter(label == "c2") |> 
  select(!lhs:rhs) |> 
  mutate(label = "IC → OR → IP") |> 
  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "lllr", 
    caption = 'direct effect B.',
    linesep = "",
    digits = 3
  )
