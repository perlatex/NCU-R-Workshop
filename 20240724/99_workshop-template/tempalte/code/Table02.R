library(tidyverse)
library(lavaan)
library(semTools)  

d <- haven::read_sav("./data/data.sav")  


model <- '

  style =~ style1 + style2 + style3 + style4 + style5
  emp   =~ emp1   + emp2   + emp3   + emp4   + emp5
  iden  =~ iden1  + iden2  + iden3  + iden4  + iden5
  perf  =~ perf1  + perf2  + perf3  + perf4  + perf5
       
'

fit_cfa <- cfa(model, data = d)


alpha <- fit_cfa %>% semTools::reliability()   # Cronbach alpha
CR    <- fit_cfa %>% semTools::compRelSEM()    # CR
AVE   <- fit_cfa %>% semTools::AVE()           # AVE 

dd <- tibble(
  name  = names(CR),
  alpha = alpha[1,],
  CR    = CR,
  AVE   = AVE
) 



table02 <- fit_cfa %>% 
  parameterestimates(standardized = TRUE) %>% 
  filter(op == "=~") %>% 
  mutate(pvalue = gtools::stars.pval(pvalue)) %>% 
  select(lhs, rhs, est, se, z, pvalue, std.all) %>% 
  left_join(dd, by = join_by(lhs == name)) %>% 
  flextable::flextable() %>% 
  flextable::merge_v(j = c("lhs", "alpha", "CR", "AVE")) %>% 
  flextable::valign(j = "lhs", valign = "top") %>% 
  flextable::colformat_double(digits = 3) %>% 
  flextable::autofit() %>% 
  flextable::set_caption("Reliability and validity examination.")


