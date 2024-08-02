library(tidyverse)
library(lavaan)
library(lavaanExtra)

d <- haven::read_sav("./data/data.sav")  


m_style <- '
  style =~ style1 + style2 + style3 + style4 + style5
'

m_emp <- '
  emp   =~ emp1   + emp2   + emp3   + emp4   + emp5
'


m_iden <- '
  iden  =~ iden1  + iden2  + iden3  + iden4  + iden5
'

m_perf <- '
  perf  =~ perf1  + perf2  + perf3  + perf4  + perf5
'

table01 <- lst(m_style, m_emp, m_iden, m_perf) %>% 
  purrr::map( ~cfa(.x, data = d) ) %>% 
  lavaanExtra::nice_fit(nice_table = TRUE)  %>% 
  flextable::fontsize(size = 9, part = "all") %>% 
  flextable::align(align = "center", part = "body") %>% 
  flextable::valign(valign = "center", part = "body")


