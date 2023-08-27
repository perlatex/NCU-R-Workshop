
library(tidyverse)
library(gtsummary)

rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>% 
  drop_na()


library(psych)

d4 <- rawdf %>% 
  select(ex1:ex18) 


fit_efa4 <- d4 %>% 
  fa(nfactors = 4, 
     rotate   = "varimax", 
     fm       = "pa", 
     scores   = TRUE, 
     e.values = TRUE, 
     values   = TRUE)


d4_efa_eigenvalue <- tibble(
  items = c("romance_and_relax",  
           "chance_encounter", 
           "sense_of_loss", 
           "aberration"),
  eigenvalue = fit_efa4$e.values[1:4]
)
  

d4_efa_variance <- fit_efa4$Vaccounted %>%  
  as_tibble() %>%  
  slice(2) %>% 
  set_names(c("romance_and_relax",
              "chance_encounter", 
              "sense_of_loss", 
              "aberration")) %>% 
  pivot_longer(cols = everything(), 
               names_to = "items",
               values_to = "variance explained")



d4_efa_loadings <- fit_efa4$loadings %>%  
  unclass() %>%  
  as.data.frame() %>%  
  rownames_to_column("items") %>%  
  rowwise() %>% 
  mutate(loadings_efa = max(c_across(-items))) %>% 
  ungroup() %>% 
  select(items, loadings_efa)




model_cfa <- '

  romance_and_relax =~ ex1 + ex11 + ex14 + ex2 + ex5 + ex15 + ex4
  chance_encounter  =~ ex10 + ex9 + ex8 + ex6 + ex7
  sense_of_loss     =~ ex17 + ex18 + ex16 
  aberration        =~ ex3 + ex13 + ex12
  
'

fit_cfa4 <- cfa(model = model_cfa, data = d4) 




d4_cfa_loadings <- 
  parameterestimates(fit_cfa4, standardized = TRUE) %>%  
  filter(op == "=~") %>% 
  select(items = rhs, loading_cfa = std.all)





## 合并成大表

CR    <- semTools::compRelSEM(fit_cfa4)
AVE   <- semTools::AVE(fit_cfa4)

d4_cfa_CR_AVE <- tibble(
  items = names(CR),
  CR   = CR,
  AVE  = AVE
 )


pairs <- readxl::read_excel("./data/variables.xlsx") %>% 
  select(colname, long_name) %>% 
  deframe()


target4 <- c(
  "romance_and_relax",
  "ex1",
  "ex14",
  "ex11",
  "ex2",
  "ex15",
  "ex5",
  "ex4",
  
  "chance_encounter", 
  "ex10",
  "ex9",
  "ex8",
  "ex6",
  "ex7",
  
  "sense_of_loss", 
  "ex16",
  "ex17",
  "ex18",
  
  "aberration",
  "ex3",
  "ex13",
  "ex12"
  )


table04 <- d4_cfa_loadings %>%
  left_join(d4_efa_loadings, by = join_by(items)) %>% 
  dplyr::add_row(tibble(items = c("romance_and_relax",
                                  "chance_encounter", 
                                  "sense_of_loss", 
                                  "aberration")), .before = 1) %>%  
  left_join(d4_efa_eigenvalue, by = join_by(items)) %>% 
  left_join(d4_efa_variance, by = join_by(items)) %>% 
  left_join(d4_cfa_CR_AVE, by = join_by(items)) %>% 
  arrange(factor(items, levels = target4)) %>% 
  mutate(items = recode(items, !!!pairs)) %>% 
  relocate(loadings_efa, .before = CR) %>% 
  
  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "l", 
    caption = 'EFA and CFA results of liminal experience.',
    linesep = "",
    digits = 3
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  ) 








