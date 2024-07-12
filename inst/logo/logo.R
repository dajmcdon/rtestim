library(tidyverse)
library(rtestim)
library(hexSticker)

out <- estimate_rt(cancovid$incident_cases, x = cancovid$date)

out <- with(can_tvar, data.frame(
  Rt = c(Rt),
  lam = rep(factor(seq_along(lambda)), each = nrow(Rt)),
  x = rep(x, ncol(Rt))
))


p <- ggplot(out, aes(x, Rt, group = lam)) +
  geom_line(linewidth = .1, color = "#5480be") +
  coord_cartesian(ylim = c(0.5, 2), xlim = c(ymd("2021-08-01", "2022-08-31"))) +
  theme_void() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1, color = "#5480be")

# showtext::font_add_google("Montserrat")
# put it where usethis::use_logo() would
sticker(
  p,
  package = "rtestim", filename = "man/figures/logo.png",
  s_x = .9, s_y = .9, s_width = 2, s_height = 1.8,
  p_size = 20, p_x = 1.25, p_family = "Montserrat", p_color = "white",
  h_fill = "#041E42", h_color = "#5480be", h_size = 1.5,
)
