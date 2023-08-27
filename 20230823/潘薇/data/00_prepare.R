library(tidyverse)

rawdf <- haven::read_sav("result530.sav") |> 
  haven::zap_labels() |> 
  select(!1:3) |> 
  slice_sample(prop = 0.9) |> 
  rename_with(~str_replace(.x, "SL", "DigPC"), starts_with("SL")) |> 
  rename_with(~str_replace(.x, "OIC", "ImpC"), starts_with("OIC")) |> 
  rename_with(~str_replace(.x, "IC", "OrgR"), starts_with("IC")) |> 
  rename_with(~str_replace(.x, "IP", "InnP"), starts_with("IP"))


rawdf |> write_rds("rawdf.rds")
