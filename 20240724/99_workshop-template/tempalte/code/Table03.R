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



AVE <- fit_cfa %>% semTools::AVE()  

table03 <- fit_cfa %>% 
  lavaan::lavInspect(what = "cor.lv") %>% 
  corrr::as_cordf(diagonal = sqrt(AVE)) %>% 
  corrr::shave() %>% 
  flextable::flextable() %>% 
  flextable::colformat_double(digits = 3) %>% 
  flextable::width(width = 3, unit = "cm") %>% 
  flextable::set_caption("The discriminate validity test of latent variables.") %>% 
  flextable::add_footer_lines(
    values = "The square root of the AVE of four latent constructs is given in the diagonal, and the correlation coefficient is given on the below diagonal. The bold values represent the square root of AVE."
  )
