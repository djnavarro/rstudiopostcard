
library(Rcpp)
library(tidyverse)
library(ambient)
library(flametree)  
library(ggforce)  
library(paletteer)
library(here)

sourceCpp(here("source", "stepping_stone.cpp")) 
output <- here("images", "rstudio_01.png")


# generate stepping stone background --------------------------------------

cat("generating image...\n")

# parameters
seed_ss <- 339
shades <- 1000
grains_wide <- 525
grains_high <- 740

palette <- paletteer_c(
  palette = "ggthemes::Classic Blue", 
  n = shades
)

# seed for RNG
set.seed(seed_ss)

# create long grid for the raster with appropriate aspect ratio
ar <- grains_high / grains_wide
raster <- long_grid(
  x = seq(0, 1,  length.out = grains_wide), 
  y = seq(0, ar, length.out = grains_high)
)

# initialise raster using worley noise
raster$base <- fracture(
  noise = gen_worley,
  fractal = fbm,
  octaves = 5, 
  frequency = 3,
  value = "distance2",
  seed = seed_ss, 
  x = raster$x,
  y = raster$y
) 

# convert base raster image to integer index
raster$base <- as.integer(ceiling(normalise(raster$base) * shades))

# convert to matrix and run stepping stone automaton
ss <- matrix(raster$base, grains_wide, grains_high)
ss <- t(timestep(ss, 200))

# read colours off the ss matrix 
raster$shade <- palette[ss]





# generate flametree ------------------------------------------------------

seed_ft <- 1
set.seed(seed_ft)

# create the tree
ftree <- flametree_grow(
  time = 12,
  seed = seed_ft,
  angle = c(-2:4) * 10,
  scale = c(.6, .8, .9)
)

# compute aspect ratio of generated tree
ar2 <- with(ftree, (max(coord_y) - min(coord_y))/(max(coord_x) - min(coord_x)))

# scale the tree image to fit the postcard
ftree <- ftree %>%
  mutate(
    coord_x = normalise(coord_x, to = range(raster$x)),
    coord_y = normalise(coord_y, to = range(raster$y))
  ) %>%
  mutate(
    coord_x = coord_x * min(1, ar/ar2), 
    coord_y = coord_y * min(1, ar2/ar)
  ) %>% 
  mutate(
    coord_x = .05 + coord_x * .9,
    coord_y = coord_y * .9
  )


# "leaf" coordinates are at terminal locations (id_step = 2) 
# on the terminal branches (id_leaf == TRUE) in the tree
vleaf <- ftree %>% 
  filter(id_leaf == TRUE, id_step == 2)




# generate raindrops ------------------------------------------------------

seed_rd <- 1
set.seed(seed_rd)

# rain
rain <- expand_grid(
  x = seq(.01, .99, length.out = 50),
  y = seq(.01 * 1.41, .99 * 1.41, length.out = round(50 * 1.41))
) %>% 
  mutate(droplet = runif(n()) < .1) %>%  
  group_by(x) %>%
  mutate(droplet = droplet | lag(droplet) | lag(droplet, 2)) %>%
  ungroup()

rain1 <- rain %>%
  filter(droplet == TRUE)
  
rain2 <- rain %>%
  group_by(x) %>%
  filter(droplet, !lag(droplet) | y == min(y)) %>%
  ungroup()



# render the image --------------------------------------------------------

cat("rendering image...\n")

pic <- ggplot(
  data = raster,
  mapping = aes(x, y, fill = shade)
) + 
  
  # background is the raster object
  geom_raster() + 
  
  # underlay the rain sparkle
  geom_point(
    data = rain2, 
    mapping = aes(x, y),
    size = 1,
    alpha = 1,
    color = "white",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # add the rain drops
  geom_point(
    data = rain1, 
    mapping = aes(x, y),
    size = .5,
    alpha = 1,
    color = "black",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # tree trunk is drawn using geom_bezier from ggforce
  geom_bezier(
    data = ftree, 
    mapping = aes(
      x = coord_x,
      y = coord_y,
      group = id_path,
      size = .3 + seg_wid * 3
    ), 
    lineend = "round", 
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # overlay the leaves
  geom_point(
    data = vleaf,
    mapping = aes(
      x = coord_x,
      y = coord_y
    ),
    size = 1,
    color = "white",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) + 
  
  # plot settings
  scale_fill_identity() + 
  scale_size_identity() + 
  coord_equal() + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() + 
  NULL

# export image
ggsave(
  filename = output,
  plot = pic,
  width = grains_wide / 150,
  height = grains_high / 150,
  dpi = 1200
)





