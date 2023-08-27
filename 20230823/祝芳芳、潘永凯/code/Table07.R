library(tidyverse)

rawdf2 <- haven::read_dta("./data/Dataset2.dta") 


d2 <- rawdf2 %>% 
  select(DPSS_P, TDSSmor, TDSSsex, TDSSpath, totalSDES, DPSS_S, Negative_Body_Image)


d2_std <- d2 %>% 
  mutate(across(everything(), ~ (.x - mean(.x)) /sd(.x)))


model1 <- lm(Negative_Body_Image ~ TDSSpath + TDSSsex + TDSSmor, data = d2_std)

model2 <- lm(Negative_Body_Image ~ TDSSpath + DPSS_P, data = d2_std)

model3 <- lm(Negative_Body_Image ~ TDSSpath + DPSS_P + DPSS_S, data = d2_std)



models <- list("Model 1" = model1, "Model 2" = model2, "Model 3" = model3)

table07 <- models %>% 
  purrr::imap_dfr(
    ~ broom::tidy(.x, conf.int = TRUE),
    .id = "model"
  ) %>% 
  filter(term != "(Intercept)") %>% 
  rename(Predictor = term) %>% 
  
  flextable::as_grouped_data(groups = c("model")) %>% 
  flextable::flextable() %>% 
  flextable::colformat_double(digits = 3)
