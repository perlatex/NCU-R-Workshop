library(tidyverse)
library(lavaan)
library(semTools)
options(knitr.kable.NA = '')

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

order <- c("DigPC", "DigPC1", "DigPC2", "DigPC3", "DigPC4", "DigPC5", "DigPC6", "DigPC7", "ImpC", "ImpC1", "ImpC2", "ImpC3", "ImpC4", "ImpC5", "ImpC6", "ImpC7", "ImpC8", "ImpC9", "ImpC10", "ImpC11", "ImpC12", "OrgR", "OrgR1", "OrgR2", "OrgR3", "OrgR4", "InnP", "InnP1", "InnP2", "InnP3", "InnP4", "InnP5")

table01 <- fit_cfa |> 
  parameterEstimates(standardized = TRUE) |> 
  filter(op == "=~") |> 
  select("Details" = rhs, "F.L" = std.all, "t.value" = z) |> 
  add_row(tibble(Details = c("DigPC", "ImpC", "OrgR", "InnP")), .before = 1) |> 
  left_join(
    df_ACA, by = join_by("Details")
  ) |> 
  arrange(factor(Details, levels = order)) |> 
  
  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "lrrrrr",  
    caption = 'reliability and validity.',
    linesep = "",
    digits = 3
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  )
