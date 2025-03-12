library(tidyverse)
ggplot2::ggplot()

# Aesthetic Mappings
## mpg df
ggplot2::mpg  # manufactuer/chr; model/chr; displ/dbl; year/int; cyl/int; trans/chr; drv/chr; cty/int; hwy/int; fl/chr; class/chr

## empty
ggplot(data = mpg)

## add layer of points
ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy))

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy, color = class))

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy, size = class))  # Warning 4 discrete var

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy, alpha = class))

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy, shape = class))

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy), 
        color = "blue")  # mind the () loci

# Facets
## facet_wrap()
ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy)
        ) + 
    facet_wrap(~ class, nrow = 2)

# facet_grid()
ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy)
        ) + 
    facet_grid(drv ~ cyl)

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy)
        ) + 
    facet_grid(. ~ cyl)

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy)
        ) + 
    facet_grid(drv ~ .)

# Geometric Objects
## geom_smooth()
ggplot(data = mpg) + 
    geom_smooth(mapping = aes(x = displ, y = hwy))  # -> 1 line +- sd

ggplot(data = mpg) + 
    geom_smooth(mapping = aes(x = displ, y = hwy, linetype = drv))  # -> 3 lines

## overlay
ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy)) + 
    geom_smooth(mapping = aes(x = displ, y = hwy))  # or
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) + 
    geom_point() + 
    geom_smooth()

ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) + 
    geom_point(mapping = aes(color = class)) + 
    geom_smooth()

ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) + 
    geom_point() +  # no class
    geom_smooth(
        se = FALSE  # no se
    )
## filter()
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) + 
    geom_point(mapping = aes(color = class)) +  # yes class
    geom_smooth(
        data = filter(mpg, class == "subcompact"), 
        se = FALSE
    )

## geom_smooth(std error, se)
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) + 
    geom_point(mapping = aes(color = class)) + 
    geom_smooth(
        data = filter(mpg, class == "subcompact"), 
        se = TRUE  # yes se
    )

ggplot(data = mpg) + 
    geom_smooth(
        mapping = aes(x = displ, y = hwy, color = drv), 
        show.legend = FALSE
    )

## geom_smooth(show.legend)
ggplot(data = mpg) + 
    geom_smooth(
        mapping = aes(x = displ, y = hwy, color = drv), 
        show.legend = TRUE
    )
