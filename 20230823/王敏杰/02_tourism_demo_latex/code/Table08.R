library(tidyverse)
library(gtsummary)

rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>% 
  drop_na()


d7 <- rawdf %>% 
  rowwise() %>% 
  mutate(
    Ambient  = mean(c_across(num_range("Ambient", 1:4))),
    Space    = mean(c_across(num_range("Space",   1:6))),
    Signs    = mean(c_across(num_range("Signs",   1:6)))
  ) %>% 
  ungroup() %>% 
  mutate(
    Experience_seeking       = es1  + es2,
    Thrill_adventure_seeking = tas1 + tas2,
    Disinhibition            = dis1 + dis2,
    romance_and_relax        = ex1  + ex11 + ex14 + ex2 + ex5 + ex15 + ex4,
    chance_encounter         = ex10 + ex9  + ex8  + ex6 + ex7,
    sense_of_loss            = ex17 + ex18 + ex16,
    aberration               = ex3  + ex13 + ex12
  ) 


library(lavaan)

model <- '

  Physical_tourscape  =~ Ambient + Space + Signs

  Social_tourscape    =~ Social1 + Social2 + Social3 + Social4 + Social5

  Socially_symbolic_tourscape  =~ Symbol1 + Symbol2 + Symbol3

  Natural_tourscape  =~ Nature1 + Nature2 + Nature3

  Emotional_arousal =~ em1 + em2 + em3 + em4 + em5

  Sensation_seeking =~ Experience_seeking + Thrill_adventure_seeking + Disinhibition

  Destination_familiarity =~ fm1 + fm2 + fm3 

  Liminal_experience  =~ romance_and_relax + chance_encounter + sense_of_loss + aberration


'

fit_cfa7 <- cfa(model, data = d7)



d8 <- rawdf %>% 
  rowwise() %>% 
  mutate(
    Ambient  = mean(c_across(num_range("Ambient", 1:4))),
    Space    = mean(c_across(num_range("Space", 1:6))),
    Signs    = mean(c_across(num_range("Signs", 1:6)))
  ) %>% 
  ungroup() %>% 
  transmute(
    Physical_tourscape  =  Ambient  + Space + Signs,
    Social_tourscape    = Social1 + Social2 + Social3 + Social4 + Social5,
    Socially_symbolic_tourscape  = Symbol1 + Symbol2 + Symbol3,
    Natural_tourscape  = Nature1 + Nature2 + Nature3,
    Emotional_arousal = em1 + em2 + em3 + em4 + em5,
    Sensation_seeking = es1 + es2 + tas1 + tas2 + dis1 + dis2,
    Destination_familiarity = fm1 + fm2 + fm3,
    Liminal_experience = ex1 + ex11 + ex14 + ex2 + ex5 + ex15 + ex4 +
                         ex10 + ex9 + ex8 + ex6 + ex7 +
                         ex17 + ex18 + ex16 +
                         ex3 + ex13 + ex12
 )




table08 <- d8 %>% 
  summarise(
    across(everything(), list(Mean = mean, SD = sd))
  ) %>% 
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(Mean|SD)"
  ) %>% 

  left_join(
    d8 %>% 
      corrr::correlate(
        diagonal = semTools::AVE(fit_cfa7) %>% sqrt()
      ) %>% 
      corrr::shave(),
    by = join_by("variable" == "term")
  ) %>%  
  set_names(c(" ", "Mean", "SD", "1", "2", "3", "4", "5", "6", "7", "8")) %>% 

  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "l", 
    caption = 'Discriminant validity test of all constructs (AVE test).',
    linesep = "",
    digits = 3
  ) %>%
  
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  ) %>%
  kableExtra::footnote("The bold diagonal elements are square roots of AVE for each construct. Below diagonal elements are the correlations between constructs.", threeparttable = T)







