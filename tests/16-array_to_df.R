# Test the edge correction in K
library(devtools)
load_all(".")
library(spatstat)

x0 <- x <- rmpoispp(rep(100, 5)*2)
r <- seq(.1, .3, l=10)
t0 <- system.time(  a <- pcf_cross_all_box(x, r)  )

d <- array_to_df(a, r)

library(ggplot2)
d |>
  ggplot() +
  geom_line(aes(r, value)) +
  facet_grid(i~j)


d |>
  dplyr::mutate(ij = paste0(pmin(i,j), "-", pmax(i,j)),
                dir = ifelse(i < j, "i < j", "j >= i")) |>
  ggplot() +
  geom_line(aes(r, value, col = dir, linetype = dir)) +
  facet_wrap(~ij)
