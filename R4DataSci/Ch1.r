library(tidyverse)
ggplot2::ggplot()

# mpg df
ggplot2::mpg  # manufactuer/chr; model/chr; displ/dbl; year/int; cyl/int; trans/chr; drv/chr; cty/int; hwy/int; fl/chr; class/chr

ggplot(data = mpg)

ggplot(data = mpg) + 
    geom_point(mapping = aes(x = displ, y = hwy))

