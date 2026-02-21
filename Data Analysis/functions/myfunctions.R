theme_custom <- function(base_size = 11, base_family = "sans") {
  half_line <- base_size / 2
  theme(
    line = element_line(color = "black", linewidth = .5,
                        linetype = 1, lineend = "butt"),
    rect = element_rect(fill = "white", color = NA),  # no panel box by default
    text = element_text(family = base_family, face = "plain",
                        color = "black", size = base_size,
                        lineheight = .9, hjust = .5, vjust = .5),
    
    # axis lines only on left/bottom
    axis.line = element_line(color = "black", linewidth = .5),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    
    axis.text = element_text(size = base_size * 0.9, color = "gray30"),
    axis.text.x = element_text(margin = margin(t = .8 * half_line / 2), vjust = 1),
    axis.text.y = element_text(margin = margin(r = .8 * half_line / 2), hjust = 1),
    
    axis.ticks = element_line(color = "gray30", linewidth = .5),
    axis.ticks.length = unit(half_line / 2, "pt"),
    
    axis.title.x = element_text(margin = margin(t = half_line),
                                vjust = 1, size = base_size * 1.1),
    axis.title.y = element_text(angle = 90, vjust = 1,
                                margin = margin(r = half_line),
                                size = base_size * 1.1),
    
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(base_size, "pt"),
    
    strip.background = element_rect(fill = "white", color = "gray40"),
    strip.text = element_text(color = "black", size = base_size * 0.9),
    
    plot.background = element_rect(color = NA, fill = "white"),
    plot.title = element_text(size = base_size * 1.4, hjust = 0,  # left-align is common
                              vjust = 1,
                              margin = margin(b = half_line * 1.2)),
    plot.title.position = "plot",
    plot.subtitle = element_text(size = base_size * 1.1,
                                 hjust = 0, vjust = 1,
                                 margin = margin(b = half_line * .9)),
    plot.caption = element_text(size = rel(0.8), hjust = 1, vjust = 1,
                                margin = margin(t = half_line * .9)),
    plot.caption.position = "plot",
    plot.margin = margin(1,1,1,1),
    complete = TRUE
  )
}
