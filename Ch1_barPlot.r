library(tidyverse)
ggplot2::ggplot()

## diamonds df
ggplot2::diamonds  # carat/dbl; cut/ord; color/ord; clarity/ord; depth/dbl; table/dbl; price/int; x/dbl; y/dbl; z/dbl

# Statistical Transformations
ggplot(data = diamonds)

## geom_bar()
ggplot(data = diamonds) + 
    geom_bar(mapping = aes(x = cut))  # calls
## stat_count()
ggplot(data = diamonds) + 
    stat_count(mapping = aes(x = cut))

## override the default stat
demo_data <- tribble(
    ~a,      ~b,
    "bar_1", 20, 
    "bar_2", 30, 
    "bar_3", 40
)
ggplot(data = demo_data) + 
    geom_bar(mapping = aes(x = a, y = b), 
    stat = "identity")

## y = proportion
ggplot(data = diamonds) + 
    geom_bar(mapping = aes(x = cut, y = ..prop.., group = 1))  # vs
ggplot(data = diamonds) + 
    geom_bar(mapping = aes(x = cut))

## stat_summary()
ggplot(data = diamonds) + 
    stat_summary(
        mapping = aes(x = cut, y = depth), 
        fun.ymin = min, 
        fun.ymax = max, 
        fun.y = median
    )

# Position Adjustments
ggplot(data = diamonds) + 
    geom_bar(mapping = aes(x = cut))  # vs
ggplot(data = diamonds) + 
    geom_bar(mapping = aes(x = cut, color = cut))  # colored grid, vs
ggplot(data = diamonds) + 
    geom_bar(mapping = aes(x = cut, fill = cut))  # cut cut, vs
ggplot(data = diamonds) + 
    geom_bar(mapping = aes(x = cut, fill = clarity))  # cut clarity -> stack

ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) + 
    geom_bar(alpha = 1/5, position = "identity")  # vs
ggplot(data = diamonds, mapping = aes(x = cut, color = clarity)) + 
    geom_bar(full = NA, position = "identity")  # vs
ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) + 
    geom_bar(position = "fill")  # vs
ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) + 
    geom_bar(position = "dodge")

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy))  # vs
ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy), 
    position = "jitter")

# Coordinate Systems
## geom_boxplot()
ggplot(data = mpg, mapping = aes(x = class, y = hwy)) + 
    geom_boxplot()  # vs
## coord_flip()
ggplot(data = mpg, mapping = aes(x = class, y = hwy)) + 
    geom_boxplot() + 
    coord_flip()

## geom_polygon()
nz <- map_data("nz")
ggplot(nz, aes(long, lat, group = group)) + 
    geom_polygon(fill = "white", color = "black")  # vs
## coord_quickmap()
ggplot(nz, aes(long, lat, group = group)) + 
    geom_polygon(fill = "white", color = "black") + 
    coord_quickmap()

## coord_polar()
bar <- ggplot(data = diamonds) + 
    geom_bar(mapping = aes(x = cut, fill = cut), 
        show.legend = FALSE, 
        width = 1) + 
    theme(aspect.ratio = 1) + 
    labs(x = NULL, y = NULL)
bar
bar + coord_flip()
bar + coord_polar()

## geom_abline()
ggplot(data = mpg, mapping = aes(x = cty, y = hwy)) + 
    geom_point()  # vs
ggplot(data = mpg, mapping = aes(x = cty, y = hwy)) + 
    geom_point() + 
    geom_abline()  # vs
## geom_fixed()
ggplot(data = mpg, mapping = aes(x = cty, y = hwy)) + 
    geom_point() + 
    geom_fixed()
