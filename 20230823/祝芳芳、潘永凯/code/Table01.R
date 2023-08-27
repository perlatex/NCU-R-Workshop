library(tidyverse)

rawdf1 <- haven::read_dta("./data/Dataset1.dta") 


d1 <- rawdf1 %>% 
  select(DPSS_P, totalSDS, DPSS_S, Negative_Body_Image) 


table01 <- d1 %>% 
  summarise(
    across(everything(), list(M = mean, SD = sd))
  ) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(M|SD)"
  ) %>% 
  left_join(
    d1 %>% corrr::correlate() %>% corrr::shave(),
    by = join_by(variable == term)
  ) %>% 
  select(-last_col()) %>% 
  flextable::flextable() %>% 
  flextable::autofit() %>% 
  flextable::colformat_double(digits = 3)



