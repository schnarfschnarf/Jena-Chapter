library(tidyverse)
library(ggdist)
thou = read.csv("fluxes/fluxes1000_dir10.csv") 

thou$Treatment = factor(thou$Treatment, levels = c("Treatment3","Treatment2","Treatment1"))
thou %>%
  ggplot(aes(x = log2(Plant.Richness), y = log10(tot.flux), 
             color = Treatment,
             group = Plot)) +
  stat_pointinterval(position = position_dodge(),
                     point_interval = "mean_hdi",
                     .width = c(0.75, 0.99),) +
  facet_wrap(~Treatment)+  
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 6.321928),
                     labels = c('1', '2', '4', '8', '16', '60')) +
  coord_cartesian(ylim = c(1.7,3)) +
  scale_fill_manual(values = c("#F3BE6140","#AA422E40","#6C6F8040"), 
                    name = "history treatment",
                    labels = c("soil (+), plant (+)", 
                               "soil (+), plant (--)", 
                               "soil (--), plant (--)")) +
  labs(#title = "Energy flux in the soil invertebrate food-web",
    y = "Community level energy flux log10(J/h\u00b7m\u00b2)",
    x = "Plant richness") +
  scale_color_manual(values = c("#F3BE61","#AA422E","#6C6F80"), 
                     name = "history treatment",
                     labels = c("soil (+), plant (+)", 
                                "soil (+), plant (--)", 
                                "soil (--), plant (--)")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2)))
