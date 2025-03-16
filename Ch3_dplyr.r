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

## arrange(), sort by
arrange(flights, year, month, day)  # ascending order
arrange(flights, desc(arr_delay))  # descending order
df <- tibble(x = c(5, 2, NA))

arrange(df, x)  # NA is sorted at the end
arrange(df, desc(x))  # NA is sorted at the end

## select()
select(flights, year, month, day)  # select col by name
select(flights, year:day)  # select col between year and day (inclusive)
select(flights, -(year:day))  # select all col except those from year to day
select(flights, starts_with("sched"))
select(flights, ends_with("time"))
select(flights, contains("dep"), contains("arr"))
select(flights, matches("(.)\\1"))  # any variables contain repeated char
df <- tibble(x1 = c(1, 2, 3), x2 = c(2, 3, 4), x3 = c(3, 4, 5), x4 = c(4, 5, 6))
select(df, num_range("x", 1:3))  # matches x1, x2, x3, so x4 is unselected
select(flights, time_hour, air_time, everything())  # move col to the start of the df

## rename(), a variant of select()
rename(flights, tail_num = tailnum)  # new_name = old_name

## mutate(), the new variables are kept at the end of df
flights_sml <- select(flights, 
    year:day, 
    ends_with("delay"), 
    distance, air_time)
View(flights_sml)
mutate(flights_sml, 
    gain = arr_delay - dep_delay, 
    speed = distance / air_time * 60)
mutate(flights_sml, 
    gain = arr_delay - dep_delay, 
    hours = air_time / 60, 
    gain_per_hour = gain / hours)

## transmute(), only keep the new variables
transmute(flights, 
    gain = arr_delay - dep_delay, 
    hours = air_time / 60, 
    gain_per_hour = gain / hours)

