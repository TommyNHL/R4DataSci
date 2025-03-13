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

