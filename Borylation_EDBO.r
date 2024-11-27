library(tidyverse)

df <- read_csv('pred_Borylation_Updated (2).csv')

df$Yield_predicted_mean <- scale(df$Yield_predicted_mean, center = FALSE, scale = max(df$Yield_predicted_mean)/100)

plot <- ggplot(df, aes(x = Yield_predicted_mean, y = Cost_predicted_mean))+
  geom_point(aes(color = Solvent_Fraction, shape = Co_Solvent, size = Concentration), alpha = 0.8)+
  scale_color_viridis_c()+
  facet_grid(rows = vars(Nucleophile_Equiv), cols = vars(Catalyst_Equiv), labeller = label_both)
plot

test <- df %>% filter(Index == 497)
