library(tidyverse)
library(gtsummary)

rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>%
  drop_na()


d2 <- rawdf %>%
  rowwise() %>%
  mutate(
    Ambient  = mean(c_across(num_range("Ambient", 1:4))),
    Space    = mean(c_across(num_range("Space", 1:6))),
    Signs    = mean(c_across(num_range("Signs", 1:6)))
  ) %>%
  ungroup() %>%
  select(
    Ambient, Space, Signs,
    starts_with("Social"),
    starts_with("Symbol"),
    starts_with("Nature")
  )


library(lavaan)

model <- "

  Physical_tourscape =~ Ambient + Space + Signs 
  Social_tourscape   =~ Social1 + Social2 + Social3 + Social4 + Social5
  Socially_symbolic_tourscape =~ Symbol1 + Symbol2 + Symbol3
  Natural_tourscape =~ Nature1 + Nature2 + Nature3

"

fit_cfa2 <- cfa(model, data = d2)


d3 <- rawdf %>%
  rowwise() %>%
  mutate(
    Ambient  = mean(c_across(num_range("Ambient", 1:4))),
    Space    = mean(c_across(num_range("Space", 1:6))),
    Signs    = mean(c_across(num_range("Signs", 1:6)))
  ) %>%
  ungroup() %>%
  transmute(
    Physical_tourscape = Ambient + Space + Signs,
    Social_tourscape = Social1 + Social2 + Social3 + Social4 + Social5,
    Socially_symbolic_tourscape = Symbol1 + Symbol2 + Symbol3,
    Natural_tourscape = Nature1 + Nature2 + Nature3
  )


table03 <- d3 %>%
  summarise(
    across(everything(), list(Mean = mean, SD = sd))
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.*)_(Mean|SD)"
  ) %>%
  left_join(
    d3 %>%
      corrr::correlate(
        diagonal = semTools::AVE(fit_cfa2) %>% sqrt()
      ) %>%
      corrr::shave(),
    by = join_by("variable" == "term")
  ) %>%
  set_names(c(" ", "Mean", "SD", "1", "2", "3", "4")) %>%
  
  knitr::kable(
    "latex",
    booktabs = TRUE,
    align = "l",
    caption = "Discriminant validity test of tourscape.",
    linesep = "",
    digits = 3
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  ) %>%
  kableExtra::footnote("The bold diagonal elements are square roots of AVE for each construct. Below diagonal elements are the correlations between constructs.", threeparttable = T)
