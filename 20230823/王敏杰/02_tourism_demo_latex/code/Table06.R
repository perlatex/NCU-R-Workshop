library(tidyverse)
library(gtsummary)

rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>% 
  drop_na()


d6 <- rawdf %>% 
  select(ex1:ex18) 


library(lavaan)

model_cfa <- '

  romance_and_relax =~ ex1 + ex11 + ex14 + ex2 + ex5 + ex15 + ex4
  chance_encounter  =~ ex10 + ex9 + ex8 + ex6 + ex7
  sense_of_loss     =~ ex17 + ex18 + ex16 
  aberration        =~ ex3 + ex13 + ex12
  
'

fit_cfa6 <- cfa(model = model_cfa, data = d6) 



table06 <- fit_cfa6 %>% 
  semTools::discriminantValidity() %>%  
  select(lhs, rhs, est, ci.lower, ci.upper) %>%
  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "l", 
    caption = 'Discriminant validity test of sub-dimensions of liminal experience (confidence interval test).',
    linesep = "",
    digits = 3
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  ) 




