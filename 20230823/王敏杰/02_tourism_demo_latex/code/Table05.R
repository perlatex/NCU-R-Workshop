library(tidyverse)
library(gtsummary)

rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>%
  drop_na()


d5 <- rawdf %>%
  select(ex1:ex18)



library(lavaan)

model_cfa <- "

  romance_and_relax =~ ex1 + ex11 + ex14 + ex2 + ex5 + ex15 + ex4
  chance_encounter  =~ ex10 + ex9 + ex8 + ex6 + ex7
  sense_of_loss     =~ ex17 + ex18 + ex16 
  aberration        =~ ex3 + ex13 + ex12
  
"

fit_cfa5 <- cfa(model = model_cfa, data = d5)



d5 <- rawdf %>%
  transmute(
    romance_and_relax = ex1 + ex11 + ex14 + ex2 + ex5 + ex15 + ex4,
    chance_encounter  = ex10 + ex9 + ex8 + ex6 + ex7,
    sense_of_loss     = ex17 + ex18 + ex16,
    aberration        = ex3 + ex13 + ex12
  )


table05 <- d5 %>%
  summarise(
    across(everything(), list(Mean = mean, SD = sd))
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(Mean|SD)"
  ) %>%
  left_join(
    d5 %>%
      corrr::correlate(
        diagonal = semTools::AVE(fit_cfa5) %>% sqrt()
      ) %>%
      corrr::shave(),
    by = join_by("variable" == "term")
  ) %>%
  set_names(c(" ", "Mean", "SD", "1", "2", "3", "4")) %>%
  knitr::kable(
    "latex",
    booktabs = TRUE,
    align = "l",
    caption = "Discriminant validity test of sub-dimensions of liminal experience (AVE test).",
    linesep = "",
    digits = 3
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  ) %>%
  kableExtra::footnote("The bold diagonal elements are square roots of AVE for each construct. Below diagonal elements are the correlations between constructs.", threeparttable = T)





