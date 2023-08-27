library(tidyverse)

rawdf2 <- haven::read_dta("./data/Dataset2.dta") 

d2 <- rawdf2 %>% 
  select(DPSS_P, TDSSmor, TDSSsex, TDSSpath, totalSDES, DPSS_S, Negative_Body_Image)


table04 <- d2 %>% 
  summarise(
    across(everything(), list(M = mean, SD = sd))
  ) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(M|SD)"
  ) %>% 
  left_join(
    d2 %>% corrr::correlate() %>% corrr::shave(),
    by = join_by(variable == term)
  ) %>% 
  select(-last_col()) %>% 
  flextable::flextable() %>% 
  flextable::colformat_double(digits = 3) %>%  
  flextable::fontsize(size = 9, part = "all") %>% 
  flextable::autofit() 
