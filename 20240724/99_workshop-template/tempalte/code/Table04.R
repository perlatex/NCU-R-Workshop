library(tidyverse)
library(lavaan)
library(semTools)  

d <- haven::read_sav("./data/data.sav")  


model <- '
  # 测量模型
  style =~ style1 + style2 + style3 + style4 + style5
  emp   =~ emp1   + emp2   + emp3   + emp4   + emp5
  iden  =~ iden1  + iden2  + iden3  + iden4  + iden5
  perf  =~ perf1  + perf2  + perf3  + perf4  + perf5

  # 结构模型
  emp  ~ a1 * style
  iden ~ a2 * style + d21 * emp
  perf ~ b1 * emp + b2 * iden + cprime * style
       
'

fit <- sem(model, data = d)

table04 <- fit %>% 
  lavaanExtra::nice_fit(nice_table = TRUE) %>% 
  flextable::fontsize(size = 9, part = "all") %>% 
  flextable::align(align = "center", part = "body") %>% 
  flextable::valign(valign = "center", part = "body") %>% 
  flextable::set_caption("Goodness of fit index of the structural model.")




