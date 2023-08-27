library(tidyverse)
library(lavaan)


rawdf2 <- haven::read_dta("./data/Dataset2.dta") 

d2 <- rawdf2 %>% 
  select(DPSS_P, TDSSmor, TDSSsex, TDSSpath, totalSDES, DPSS_S, Negative_Body_Image)



d2_c <- d2 %>% 
  mutate(
    across(c(totalSDES, DPSS_S), ~.x - mean(.x))
  )


model <- '

   totalSDES ~ a * DPSS_P
   Negative_Body_Image ~ b1 * totalSDES + cprime * DPSS_P + b2 * DPSS_S + b3 * totalSDES:DPSS_S

   # Index of moderated mediation
   index.mod.med := a*b3

'

fit <- sem(model, data = d2_c)



tbl_orders <- c("cprime", 
                "b1",
                "b2",
                "b3",
                "index.mod.med")


table06 <- fit %>%  
  parameterEstimates(standardized = TRUE) %>%  
  filter(label %in% tbl_orders) %>% 
  select(-c(lhs, op, rhs, std.lv, std.nox, std.all)) %>% 
  arrange(factor(label, levels = tbl_orders)) %>%  
  flextable::flextable() %>%
  flextable::colformat_double(digits = 3) %>% 
  flextable::color(j = ~est, color = "red")




