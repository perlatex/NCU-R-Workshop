library(tidyverse)
library(gtsummary)

rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>% 
  drop_na()



d11 <- rawdf %>% 
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

  # latent variables
  Physical  =~ Ambient + Space + Signs
  Social    =~ Social1 + Social2 + Social3 + Social4 + Social5
  Symbolic  =~ Symbol1 + Symbol2 + Symbol3
  Natural   =~ Nature1 + Nature2 + Nature3

  Emotional_arousal =~ em1 + em2 + em3 + em4 + em5

  Sensation_seeking =~ Experience_seeking + Thrill_adventure_seeking + Disinhibition

  Destination_familiarity =~ fm1 + fm2 + fm3 

  Liminal_experience  =~ romance_and_relax + chance_encounter + sense_of_loss + aberration


  # regressions
  Emotional_arousal ~  H1a * Physical + H2a * Social + H3a *Symbolic +  H4a * Natural
  Liminal_experience  ~ H1b * Physical +  H2b * Social + H3b * Symbolic + H4b * Natural + H5 * Emotional_arousal + H7 *Sensation_seeking + H8 * Destination_familiarity

  # define
  MH6a  := H5 * H1a
  MH6b  := H5 * H2a
  MH6c  := H5 * H3a
  MH6d  := H5 * H4a

'

fit_sem11 <- sem(model, data = d11)




table11 <- fit_sem11 %>% 
  parameterEstimates() %>% 
  filter(op %in% c("~", ":=")) %>% 
  filter(str_detect(label, "^M")) %>% 
  mutate(
    `Independent variable` = str_c(c("Physical", "Social", "Socially symbolic", "Natural"), " tourscape"),
    Mediator = c("Emotional arousal", NA, NA, NA),
    `Dependent variable` = c("Liminal experience", NA, NA, NA)
  ) %>% 
  select(label, `Independent variable`, Mediator, `Dependent variable`,
         #est, z, pvalue, 
         ci.lower, ci.upper) %>% 
  mutate(label = str_remove(label, "M")) %>% 
  knitr::kable(
    "latex", 
    booktabs = TRUE, 
    align = "llllrr", 
    caption = 'Mediating effect (bootstrap =2000).',
    linesep = "",
    digits = 3
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  ) 


