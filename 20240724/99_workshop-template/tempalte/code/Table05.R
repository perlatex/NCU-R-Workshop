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
  emp  ~ H2 * style
  iden ~ H5 * style + H6 * emp
  perf ~ H1 * style + H3 * emp + H7 * iden

  # 效应量
  serial_indirect_eff := H2*H6*H7
  indirect_eff        := H2*H6*H7 + H2*H3 + H5*H7
  direct_eff          := H1
  total_eff           := indirect_eff + direct_eff 
  
'



fit_sem <- sem(model, 
               data      = d, 
               se        = "bootstrap",
               bootstrap = 500,
               mimic     = "Mplus"
)



table05 <- fit_sem %>%  
  parameterEstimates(standardized = TRUE) %>% 
  filter(str_detect(label, "^H")) %>% 
  mutate(path = str_c(rhs, " --> ", lhs)) %>% 
  mutate(Remark = if_else(pvalue > 0.05, "Not supported", "Supported")) %>% 
  select(label, path, std.all, se, pvalue, Remark) %>% 
  mutate(" " = gtools::stars.pval(pvalue), .after = pvalue) %>%
  arrange(label) %>% 
  flextable::flextable() %>% 
  flextable::autofit() %>% 
  flextable::colformat_double(digits = 3) %>% 
  flextable::set_caption("The test results of path relationship.") %>% 
  flextable::autofit() %>% 
  flextable::add_footer_lines(
    values = "Standard error. ***p < 0.001, **p < 0.01, *p < 0.05, +p < 0.10"
  )


