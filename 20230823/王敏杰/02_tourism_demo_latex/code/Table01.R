library(tidyverse)
library(scales)

rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>% 
  drop_na() %>% 
  filter(gender %in% c(1, 2)) %>% 
  filter(edu %in% c(1, 2, 3)) 
  

d1 <- rawdf %>% 
  select(gender, age, edu, income, time) %>% 
  mutate(gender = factor(gender, 
                         labels = c(
                           "1" = "male",
                           "2" = "female"
                         )),
         age    = factor(age, 
                         labels = c(
                           "1" = "under 20",
                           "2" = "21-30",
                           "3" = "31-40",
                           "4" = "above 41"
                         )),
         edu    = factor(edu,
                         labels = c(
                           "1" = "senior high school and below",
                           "2" = "university",
                           "3" = "master and above"
                         )),
         income = factor(income,
                         labels = c(
                           "1" = "under 1500",
                           "2" = "1500-3000",
                           "3" = "3001-5000",
                           "4" = "5001-7500",
                           "5" = "above 7501"
                         )),
         time   = factor(time,
                         labels = c(
                           "1" = "1",
                           "2" = "2-3",
                           "3" = "4-5",
                           "4" = "> 6"
                         )))





table01 <- d1 %>% 
  pivot_longer(cols = everything()) %>% 
  summarise(Frequency = n(), .by = c(name, value) ) %>% 
  mutate(
    Percent = scales::label_percent(accuracy = 0.01)(Frequency/sum(Frequency)), 
    .by = name
  ) %>% 
  arrange(name, value) %>% 
  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "lllr", 
    caption = 'Sample profile.',
    linesep = "",
    digits = 3
  ) %>%
  
  kableExtra::column_spec(1, bold = T, width = "5em") %>%
  kableExtra::collapse_rows(
    columns = 1, 
    valign  = "top",
    latex_hline = "none"
  ) %>% 
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  ) 
  


