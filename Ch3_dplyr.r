# Data Transformation with dplyr. dplyr overrided stats::filter() and stats::lag()
library(nycflights13)  # data
library(tidyverse)

## df flights
nycflights13::flights  # year/int; month/int; day/int; dep_time/int; sched_dep_time/int; dep_delay/dbl; arr_time/int; sched_arr_time/int
View(flights)  # tibble

## data types
# int, dbl, chr, dttm(a date + a time), lgl(bool), fctr(factors/cat. variables), date

# dplyr basics
## Logical Operations, &, |, !, xor()
## Missing Values
NA > 5
10 == NA
NA + 10
NA / 2
NA == NA
x <- NA
y <- NA
x == y
is.na(x)
df <- tibble(x = c(1, NA, 3))
filter(df, x > 1)
filter(df, is.na(x) | x > 1)
## filter()
filter(flights, month == 1, day == 1)  # >, <, >=, <=, ==, !=
jan1 <- filter(flights, month == 1, day == 1)
jan1
(dec25 <- filter(flights, month == 12, day == 25))
sqrt(2) ^ 2 == 2  # -> FALSE
1/49 * 49 == 1  # -> FALSE
near(sqrt(2) ^ 2, 2)  # -> TRUE
near(1/49 * 49, 1)  # -> TRUE
filter(flights, month == 11 | month == 12)  # or
filter(flights, month %in% c(11, 12))
filter(flights, !(arr_delay > 120 | dep_delay > 120))  # or
filter(flights, arr_delay <= 120, dep_delay <= 120)

