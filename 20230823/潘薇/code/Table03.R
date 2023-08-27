library(tidyverse)
library(lavaan)
library(semTools)
rawdf <- readr::read_rds("./data/rawdf.rds")

model <- '
  DigPC =~ DigPC1 + DigPC2 + DigPC3 + DigPC4 + DigPC5 + DigPC6 + DigPC7
  ImpC  =~ ImpC1 + ImpC2 + ImpC3 + ImpC4 + ImpC5 + ImpC6 + ImpC7 + ImpC8 + ImpC9 + ImpC10 + ImpC11 + ImpC12
  OrgR  =~ OrgR1 + OrgR2 + OrgR3 + OrgR4
  InnP  =~ InnP1 + InnP2 + InnP3 + InnP4 + InnP5
'

fit_cfa <- cfa(model, data = rawdf)

Alpha <- reliability(fit_cfa)[1, ]
CR    <- compRelSEM(fit_cfa)
AVE   <- AVE(fit_cfa)

df_ACA <- tibble(
  Details = names(Alpha),
  Alpha   = Alpha,
  CR      = CR,
  AVE     = AVE
)

df_correlation <- fit_cfa |> 
  lavInspect(what = "cor.lv") |> 
  as.data.frame() |> 
  rownames_to_column("Variable")

df <- rawdf |> 
  rowwise() |> 
  transmute(
    DigPC = mean(c(DigPC1, DigPC2, DigPC3, DigPC4, DigPC5, DigPC6, DigPC7)),
    ImpC  = mean(c(ImpC1, ImpC2, ImpC3, ImpC4, ImpC5, ImpC6, ImpC7, ImpC8, ImpC9, ImpC10, ImpC11, ImpC12)),
    OrgR  = mean(c(OrgR1, OrgR2, OrgR3, OrgR4)),
    InnP  = mean(c(InnP1, InnP2, InnP3, InnP4, InnP5))
  ) |> 
  ungroup()

table03 <- df |> 
  summarise(
    across(everything(), list(Mean = mean, SD = sd))
  ) |> 
  pivot_longer(
    everything(),
    names_to      = c("Variable", ".value"),
    names_pattern = "(.*)_(Mean|SD)"
  ) |>
  left_join(
    df_ACA |> select(Details, Alpha),
    by = join_by("Variable" == Details)
  ) |> 
  left_join(
    df_correlation,
    by = join_by("Variable")
  ) |> 
  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "lllr", 
    caption = 'mean sd and correlation coefficient.',
    linesep = "",
    digits = 3
  )
