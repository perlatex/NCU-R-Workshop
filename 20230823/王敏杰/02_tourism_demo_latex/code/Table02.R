library(tidyverse)
library(gtsummary)
options(knitr.kable.NA = '')


rawdf <- haven::read_sav("./data/Liminal_Experience.sav") %>% 
  drop_na()

pairs <- readxl::read_excel("./data/variables.xlsx") %>% 
  select(colname, long_name) %>% 
  deframe()


d2 <- rawdf %>% 
  rowwise() %>% 
  mutate(
    Ambient  = mean(c_across(num_range("Ambient", 1:4))),
    Space    = mean(c_across(num_range("Space",   1:6))),
    Signs    = mean(c_across(num_range("Signs",   1:6)))
  ) %>% 
  ungroup() %>% 
  select(
    Ambient, Space, Signs,
    starts_with("Social"),
    starts_with("Symbol"),
    starts_with("Nature")
  )


library(lavaan)

model <- '

  Physical_tourscape =~ Ambient + Space + Signs 
  Social_tourscape   =~ Social1 + Social2 + Social3 + Social4 + Social5
  Socially_symbolic_tourscape =~ Symbol1 + Symbol2 + Symbol3
  Natural_tourscape =~ Nature1 + Nature2 + Nature3

'

fit_cfa2 <- cfa(model, data = d2)


CR   <- semTools::compRelSEM(fit_cfa2)
AVE  <- semTools::AVE(fit_cfa2)

d2_CR_AVE <- tibble(
  items = names(CR),
  CR    = CR,
  AVE   = AVE
 )



target <- c(
  "Physical_tourscape", "Ambient", "Space", "Signs" ,
  "Social_tourscape", "Social1", "Social2", "Social3", "Social4", "Social5",
  "Socially_symbolic_tourscape", "Symbol1" , "Symbol2" , "Symbol3",
  "Natural_tourscape", "Nature1" , "Nature2" , "Nature3"
  )


table02 <- parameterestimates(fit_cfa2, standardized = TRUE) %>% 
  filter(op == "=~") %>% 
  select(items = rhs, Loading = std.all) %>% 
  dplyr::add_row(tibble(items = c("Physical_tourscape",
                                  "Social_tourscape", 
                                  "Socially_symbolic_tourscape", 
                                  "Natural_tourscape")), .before = 1) %>%  
  left_join(d2_CR_AVE, by = join_by(items)) %>% 
  arrange(factor(items, levels = target)) %>% 
  mutate(items = recode(items, !!!pairs)) %>% 
  
  knitr::kable(
      "latex", 
      booktabs = TRUE, 
      align = "lllc", 
      caption = 'CFA results of tourscape.',
      linesep = "",
      digits = 3
  ) %>%
  kableExtra::kable_styling(
    latex_options = c("HOLD_position"),
    position      = "center"
  ) 
