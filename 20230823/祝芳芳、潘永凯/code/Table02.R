library(tidyverse)
library(lavaan)


rawdf1 <- haven::read_dta("./data/Dataset1.dta") 


d1 <- rawdf1 %>% 
  select(DPSS_P, totalSDS, DPSS_S, Negative_Body_Image) 



model <- '

    totalSDS ~ a * DPSS_P
    Negative_Body_Image ~ b * totalSDS + cprime * DPSS_P

    indirect := a*b
    total    := cprime + a*b

'

fit <- sem(model, data = d1)




tbl_orders <- c("total", 
                "a", 
                "cprime", 
                "b",
                "indirect")


table02 <- fit %>%  
  parameterEstimates(standardized = TRUE) %>%  
  filter(label %in% tbl_orders) %>% 
  select(-c(lhs, op, rhs, std.lv, std.nox, std.all)) %>% 
  arrange(factor(label, levels = tbl_orders)) %>%  
  flextable::flextable() %>%
  flextable::colformat_double(digits = 3) %>% 
  flextable::color(j = ~est, color = "red")



