library(tidyverse)
library(lavaan)


rawdf2 <- haven::read_dta("./data/Dataset2.dta") 

d2 <- rawdf2 %>% 
  select(DPSS_P, TDSSmor, TDSSsex, TDSSpath, totalSDES, DPSS_S, Negative_Body_Image)


model <- '

    totalSDES ~ a * DPSS_P
    Negative_Body_Image ~ b * totalSDES + cprime * DPSS_P

    indirect := a*b
    total := cprime + (a*b)

'

fit <- sem(model, data = d2)



tbl_orders <- c("total", 
                "a", 
                "cprime", 
                "b",
                "indirect")


table05 <- fit %>%  
  parameterEstimates(standardized = TRUE) %>%  
  filter(label %in% tbl_orders) %>% 
  select(-c(lhs, op, rhs, std.lv, std.nox, std.all)) %>% 
  arrange(factor(label, levels = tbl_orders)) %>%  
  flextable::flextable() %>%
  flextable::colformat_double(digits = 3) %>% 
  flextable::color(j = ~est, color = "red")



