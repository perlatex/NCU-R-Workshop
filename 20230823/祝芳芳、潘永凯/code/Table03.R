library(tidyverse)
library(lavaan)


rawdf1 <- haven::read_dta("./data/Dataset1.dta") 


d1 <- rawdf1 %>% 
  select(DPSS_P, totalSDS, DPSS_S, Negative_Body_Image) 



d1_c <- d1 %>% 
  mutate(
    across(c(totalSDS, DPSS_S), ~.x - mean(.x))
  )


model <- '

    totalSDS ~ a*DPSS_P
    Negative_Body_Image ~ b1 * totalSDS + cprime * DPSS_P + b2 * DPSS_S + b3 * totalSDS:DPSS_S

   # Index of moderated mediation
   index.mod.med := a*b3
'

fit <- sem(model, data = d1_c)



tbl_orders <- c("cprime", 
                "b1",
                "b2",
                "b3",
                "index.mod.med")


table03 <- fit %>%  
  parameterEstimates(standardized = TRUE) %>%  
  filter(label %in% tbl_orders) %>% 
  select(-c(lhs, op, rhs, std.lv, std.nox, std.all)) %>% 
  arrange(factor(label, levels = tbl_orders)) %>%  
  flextable::flextable() %>%
  flextable::colformat_double(digits = 3) %>% 
  flextable::color(j = ~est, color = "red")



