library(tidyverse)
library(gtsummary)

rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>% 
  drop_na()


da <- rawdf %>% 
  select(-gender, -age, -edu, -income, -time)

  
table_appendixr <- da %>% 
  psych::describe() %>% 
  as.data.frame() %>% 
  rownames_to_column("Variables/items") %>% 
  select(`Variables/items`, mean, sd, skew, kurtosis) %>% 
  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "l", 
    caption = 'A Tale of Two Tables.',
    linesep = "",
    digits = 3
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("hold_position"),
    position      = "center"
  ) 

