x1 <- list(
  c(0.27, 0.37, 0.57, 0.91, 0.20),
  c(0.90, 0.94, 0.66, 0.63, 0.06), 
  c(0.21, 0.18, 0.69, 0.38, 0.77)
)
x2 <- list(
  c(0.50, 0.72, 0.99, 0.38, 0.78), 
  c(0.93, 0.21, 0.65, 0.13, 0.27), 
  c(0.39, 0.01, 0.38, 0.87, 0.34)
)

threshold <- function(x, cutoff = 0.8) x[x > cutoff]
x1 %>% sapply(threshold) %>% str()
x1 %>% lapply(threshold) %>% str()
x2 %>% sapply(threshold) %>% str()
x2 %>% lapply(threshold) %>% str()

# http://r4ds.had.co.nz/relational-data.html#relational-data
install.packages("nycflights13")
install.packages("lubridate")

library(nycflights13)
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate)

# http://r4ds.had.co.nz/dates-and-times.html
flights %>% select(year, month, day, hour, minute)
flights

datetimes <- flights %>% 
  mutate(departure = make_datetime(year = year, month = month, day = day, 
                                   hour = hour, min = minute))
