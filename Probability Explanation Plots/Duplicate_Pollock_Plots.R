#####################
### Load Packages ###
#####################

library(ggplot2)
library(magick)
library(magrittr)
library(tmvtnorm)
library(dplyr)

#####################
### Simulate Data ###
#####################

## Two species independent of each other.

### Species the first. Pr(z>0)=0.76

sp1 <- data.frame(points = rnorm(1000000,
                                 mean = 0.5,
                                 sd = 1))

sp1_PrOcc_marginal <- pnorm(q = 0,
                            mean = 0.5,
                            sd = 1,
                            lower.tail = FALSE)

sp1_density <- data.frame(x = density(sp1$points)$x,
                          y = density(sp1$points)$y)

### Species the second. Pr(z>0)=0.38

sp2 <- data.frame(points = rnorm(1000000,
                                 mean = -1,
                                 sd = 1))

sp2_PrOcc_marginal <- pnorm(q = 0,
                            mean = -1,
                            sd = 1,
                            lower.tail = FALSE)

sp2_density <- data.frame(x = density(sp2$points)$x,
                          y = density(sp2$points)$y)

## Two species with correlations

mean <- c(0.5, -1)

m_spp <- rmvnorm(100000,
                 mean = mean,
                 sigma = matrix(c(1, 0.75, 0.75, 1),
                                2,2))
colnames(m_spp) <- c("x","y")

p_1_1 <- ptmvnorm(mean = mean,
                  sigma = matrix(c(1, 0.75, 0.75, 1),
                                 2,2),
                  lowerx = c(0, 0),
                  upperx = c(Inf, Inf),
                  lower = c(-Inf, -Inf),
                  upper = c(Inf, Inf))

p_1_0 <- ptmvnorm(mean = mean,
                  sigma = matrix(c(1, 0.75, 0.75, 1),
                                 2,2),
                  lowerx = c(0, -Inf),
                  upperx = c(Inf, 0),
                  lower = c(-Inf, -Inf),
                  upper = c(Inf, Inf))

p_0_1 <- ptmvnorm(mean = mean,
                  sigma = matrix(c(1, 0.75, 0.75, 1),
                                 2,2),
                  lowerx = c(-Inf, 0),
                  upperx = c(0, Inf),
                  lower = c(-Inf, -Inf),
                  upper = c(Inf, Inf))

p_0_0 <- ptmvnorm(mean = mean,
                  sigma = matrix(c(1, 0.75, 0.75, 1),
                                 2,2),
                  lowerx = c(-Inf, -Inf),
                  upperx = c(0, 0),
                  lower = c(-Inf, -Inf),
                  upper = c(Inf, Inf))

## Two species with correlations

p_1_1_c <- ptmvnorm(mean = mean,
                    sigma = matrix(c(1, 0.75, 0.75, 1),
                                   2,2),
                    lowerx = c(-Inf, 0),
                    upperx = c(Inf, Inf),
                    lower = c(0, -Inf),
                    upper = c(Inf, Inf))

p_1_0_c <- ptmvnorm(mean = mean,
                    sigma = matrix(c(1, 0.75, 0.75, 1),
                                   2,2),
                    lowerx = c(-Inf, -Inf),
                    upperx = c(Inf, 0),
                    lower = c(0, -Inf),
                    upper = c(Inf, Inf))

p_0_1_c <- ptmvnorm(mean = mean,
                  sigma = matrix(c(1, 0.75, 0.75, 1),
                                 2,2),
                  lowerx = c(-Inf, 0),
                  upperx = c(0, Inf),
                  lower = c(0, -Inf),
                  upper = c(Inf, Inf))

p_0_0_c <- ptmvnorm(mean = mean,
                  sigma = matrix(c(1, 0.75, 0.75, 1),
                                 2,2),
                  lowerx = c(-Inf, -Inf),
                  upperx = c(0, 0),
                  lower = c(0, -Inf),
                  upper = c(Inf, Inf))


#############
### Plots ###
#############

## Load Animal pictures

Hyena <- image_read(paste0("C:/Users/wilko/Pictures/Animal Vector Images/",
                           "faebd7bdd15aa020e0a4d4c3b24ba05d/Hyena.png"))
Hyena <- image_scale(Hyena,"150x150")

Hyena <- image_flop(Hyena)

Hyena.C <- image_read(paste0("C:/Users/wilko/Pictures/Animal Vector Images/",
                             "faebd7bdd15aa020e0a4d4c3b24ba05d/Hyena - Cross.png"))
Hyena.C <- image_scale(Hyena.C,"150x150")

Hyena.C <- image_flop(Hyena.C)

Wildebeest <- image_read(paste0("C:/Users/wilko/Pictures/Animal Vector Images/",
                                 "faebd7bdd15aa020e0a4d4c3b24ba05d/Wildebeest.png"))
Wildebeest <- image_scale(Wildebeest,"150x150")

Wildebeest.C <- image_read(paste0("C:/Users/wilko/Pictures/Animal Vector Images/",
                                  "faebd7bdd15aa020e0a4d4c3b24ba05d/Wildebeest - Cross.png"))
Wildebeest.C <- image_scale(Wildebeest.C,"150x150")

## Single Species

### Set up axis labels

axis_df <- data.frame(x = c(-4, -2, 0.2, 2, 4, -0.3, -0.3, -0.3, -0.3),
                      y = c(-0.015, -0.015, -0.015, -0.015, -0.015, 0.1, 0.2, 0.3, 0.4),
                      z = c(-4, -2, 0, 2, 4, 0.1, 0.2, 0.3, 0.4))

### Wildebeest

p_1 <- 1 - pnorm(0, 0.5, 1)

a <- image_graph(width = 1200, height = 800, antialias = TRUE, res = 300)

ggplot(sp1, aes(points)) +
  geom_ribbon(data = subset(sp1_density, x>=0),
              aes(x = x,
                  ymax = y),
              ymin = 0,
              fill = "grey",
              alpha = 0.8) +
  stat_density(aes(x = points),
               geom = "line",
               size = 0.2) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data = axis_df, aes(x = x, y = y, label = z), size = 3) +
  annotate("text", x = 5, y = -0.02, label = "z[A]", parse = TRUE, size = 5) +
  annotate("text", x = 1, y = 0.1, label = round(p_1, 2), parse = TRUE, size = 4) +
  annotate("text", x = -4, y = 0.4, label = "a)", parse = FALSE, size = 5) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(-5,5)

dev.off()

b <- image_composite(a, Wildebeest, offset = "+1020+590")
c <- image_composite(b, Wildebeest.C, offset = "+60+590")
c

image_write(c, "spp_a.pdf", format = "pdf",density = 300)

### Hyena

p_1 <- 1 - pnorm(0, -1, 1)

a <- image_graph(width = 1200, height = 800, antialias = TRUE, res = 300)

ggplot(sp2, aes(points)) +
  geom_ribbon(data = subset(sp2_density, x>=0),
              aes(x = x,
                  ymax = y),
              ymin = 0,
              fill = "grey",
              alpha = 0.8) +
  stat_density(aes(x = points),
               geom = "line",
               size = 0.2) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(data = axis_df, aes(x = x, y = y, label = z), size = 3) +
  annotate("text", x = 5, y = -0.02, label = "z[B]", parse = TRUE, size = 5) +
  annotate("text", x = 0.5, y = 0.05, label = round(p_1, 2), parse = TRUE, size = 4) +
  annotate("text", x = -4, y = 0.4, label = "b)", parse = FALSE, size = 5) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(-5,5)

dev.off()

b <- image_composite(a, Hyena, offset = "+1010+570")
c <- image_composite(b, Hyena.C, offset = "+60+570")
c

image_write(c, "spp_b.pdf", format = "pdf",density = 300)

## Multi-species

### Set up axis labels

axis_df <- data.frame(x = c(-4, -2, 0, 2, 4, -0.20, -0.20, -0.15, -0.15),
                      y = c(-0.15, -0.15, 0, -0.15, -0.15, -4, -2, 2, 4),
                      z = c(-4, -2, 0, 2, 4, -4, -2, 2, 4))
### Plot

d <- image_graph(width = 1200, height = 1200, antialias = TRUE, res = 300)

ggplot(as.data.frame(m_spp), aes(x = x,
                                 y = y)) +
  stat_ellipse(type = "norm",
               level = 0.90,
               geom = "polygon",
               fill = "#e6e6e6",
               colour = "#262626") +
  stat_ellipse(type = "norm",
               level = 0.7,
               geom = "polygon",
               fill = "#cccccc",
               colour = "#262626") +
  stat_ellipse(type = "norm",
               level = 0.5,
               geom = "polygon",
               fill = "#a6a6a6",
               colour = "#262626") +
  stat_ellipse(type = "norm",
               level = 0.3,
               geom = "polygon",
               fill = "#737373",
               colour = "#262626") +
  stat_ellipse(type = "norm",
               level = 0.1,
               geom = "polygon",
               fill = "#404040",
               colour = "#262626") +
  stat_ellipse(type = "norm",
               level = 0.001,
               geom = "polygon",
               fill = "#1a1a1a",
               colour = "#0d0d0d") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ylim(-4,4) +
  xlim(-4,4) +
  geom_text(data = axis_df, aes(label = z), size = 3) +
  annotate("text", x = 4, y = -0.6, label = "z[A]", parse = TRUE, size = 4) +
  annotate("text", x = -0.5, y = 4, label = "z[B]", parse = TRUE, size = 4) +
  annotate("text", x = 3.8, y = 3.8, label = round(p_1_1, 2), parse = TRUE, size = 4) +
  annotate("text", x = 3.8, y = -3.8, label = round(p_1_0, 2), parse = TRUE, size = 4) +
  annotate("text", x = -3.8, y = 3.8, label = round(p_0_1, 3), parse = TRUE, size = 4) +
  annotate("text", x = -3.8, y = -3.8, label = round(p_0_0, 2), parse = TRUE, size = 4)

dev.off()

e <- image_composite(d, Hyena, offset = "+620+40")
f <- image_composite(e, Hyena.C, offset = "+620+1050")
g <- image_composite(f, Wildebeest, offset = "+1020+480")
h <- image_composite(g, Wildebeest.C, offset = "+40+480")
h

image_write(h, "multi_spp.pdf", format = "pdf", density = 300)

## Multi-species Conditional

### Set up axis labels

axis_df <- data.frame(x = c(-4, -2, 0, 2, 4, -0.20, -0.20, -0.15, -0.15),
                      y = c(-0.15, -0.15, 0, -0.15, -0.15, -4, -2, 2, 4),
                      z = c(-4, -2, 0, 2, 4, -4, -2, 2, 4))

### New function for truncated ellipse
### From here: https://stackoverflow.com/questions/45761468/constrict-ggplot-ellips-to-realistic-possible-values

StatClipEllipse <- ggproto("StatClipEllipse",
                           Stat,
                           required_aes = c("x", "y"),
                           compute_group = function(data,
                                                    scales,
                                                    type = "t",
                                                    level = 0.95,
                                                    segments = 51,
                                                    na.rm = FALSE) {
                             xx <- ggplot2:::calculate_ellipse(data = data,
                                                               vars = c("x", "y"),
                                                               type = type,
                                                               level = level,
                                                               segments = segments)
                             xx %>% mutate(x = pmax(x, 0))
                           }
)

stat_clip_ellipse <- function(mapping = NULL, 
                              data = NULL,
                              geom = "path", 
                              position = "identity",
                              ...,
                              type = "t",
                              level = 0.95,
                              segments = 51,
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatClipEllipse,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      type = type,
      level = level,
      segments = segments,
      na.rm = na.rm,
      ...
    )
  )
}

### Plot

d <- image_graph(width = 1200, height = 1200, antialias = TRUE, res = 300)

ggplot(as.data.frame(m_spp), aes(x = x,
                                 y = y)) +
  stat_clip_ellipse(type = "norm",
               level = 0.90,
               geom = "polygon",
               fill = "#e6e6e6",
               colour = "#262626") +
  stat_clip_ellipse(type = "norm",
               level = 0.7,
               geom = "polygon",
               fill = "#cccccc",
               colour = "#262626") +
  stat_clip_ellipse(type = "norm",
               level = 0.5,
               geom = "polygon",
               fill = "#a6a6a6",
               colour = "#262626") +
  stat_clip_ellipse(type = "norm",
               level = 0.3,
               geom = "polygon",
               fill = "#737373",
               colour = "#262626") +
  stat_clip_ellipse(type = "norm",
               level = 0.1,
               geom = "polygon",
               fill = "#404040",
               colour = "#262626") +
  stat_clip_ellipse(type = "norm",
               level = 0.001,
               geom = "polygon",
               fill = "#1a1a1a",
               colour = "#0d0d0d") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ylim(-4,4) +
  xlim(-4,4) +
  geom_text(data = axis_df, aes(label = z), size = 3) +
  annotate("text", x = 4, y = -0.6, label = "z[A]", parse = TRUE, size = 4) +
  annotate("text", x = -0.5, y = 4, label = "z[B]", parse = TRUE, size = 4) +
  annotate("text", x = 3.8, y = 3.8, label = round(p_1_1_c, 2), parse = TRUE, size = 4) +
  annotate("text", x = 3.8, y = -3.8, label = round(p_1_0_c, 2), parse = TRUE, size = 4) +
  annotate("text", x = -3.8, y = 3.8, label = round(p_0_1_c, 3), parse = TRUE, size = 4) +
  annotate("text", x = -3.8, y = -3.8, label = round(p_0_0_c, 2), parse = TRUE, size = 4)

dev.off()

e <- image_composite(d, Hyena, offset = "+620+40")
f <- image_composite(e, Hyena.C, offset = "+620+1050")
g <- image_composite(f, Wildebeest, offset = "+1020+480")
h <- image_composite(g, Wildebeest.C, offset = "+40+480")
h

image_write(h, "multi_spp_cond.pdf", format = "pdf", density = 300)
