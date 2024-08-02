library(tidyverse)
library(lavaan)
library(semPlot)
library(semptools)

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

  # 效应量
  serial_indirect_eff := a1*d21*b2
  indirect_eff        := a1*d21*b2 + a1*b1 + a2*b2
  direct_eff          := cprime
  total_eff           := a1*d21*b2 + a1*b1 + a2*b2 + cprime        
'

fit <- sem(model, data = d, se = "bootstrap", bootstrap = 500)




p <- fit %>% 
  semPaths(
    what           = "path", 
    whatLabels     = "std",       
    style          = "ram",    
    label.cex      = 1.2,        
    edge.label.cex = 0.9, 
    edge.color     = "blue",     
    exoVar         = FALSE,
    residScale     = 8,   
    rotation       = 2,
    nCharEdges     = 0,
    nCharNodes     = 0,          
    residuals      = TRUE,       
    intercepts     = FALSE,       
    groups         = "latents", 
    pastel         = TRUE, 
    borders        = TRUE,
    border.color   = "white",
    DoNotPlot      = TRUE,
    mar            = c(5, 8, 5, 8)
  ) 





factor_layout <- matrix(
  c(NA,    "emp",  NA, "iden", NA,
    "style",  NA,  NA,  NA,  "perf"), 
  byrow = TRUE, 2, 5)



factor_point_to <- matrix(
  c(NA,    "left",  NA, "right", NA,
    "left",  NA,  NA, NA,  "right"), 
  byrow = TRUE, 2, 5)



indicator_push <- c(
  style = 2.8,
  emp   = 2.5,
  iden  = 2.5,
  perf  = 2.8
)


indicator_spread <- c(
  style = 1.2,
  emp   = 1.2,
  iden  = 1.2,
  perf  = 1.2
)




p %>% 
  set_sem_layout(
    factor_layout    = factor_layout,
    factor_point_to  = factor_point_to,
    indicator_push   = indicator_push,
    indicator_spread = indicator_spread
  ) %>% 
  semptools::mark_sig(fit) %>% 
  set_curve( c("perf ~ emp"  = 0) ) %>%
  rotate_resid( c(perf = 180, emp = 0, iden = 0) ) %>% 
  plot()
