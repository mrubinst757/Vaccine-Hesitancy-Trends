library(covidcast)
library(tidyverse)

get_data <- function(signal, start_date, end_date) {
  
  signal0 <- covidcast_signal(
    geo_type = "county",
    data_source = "usa-facts", 
    signal = signal, 
    start_day = start_date,
    end_day = start_date
  )
  
  signal1 <- covidcast_signal(
    geo_type = "county",
    data_source = "usa-facts", 
    signal = signal, 
    start_day = end_date,
    end_day = end_date
  )

  dat0 <- signal0 %>%
    select(geo_value, start = value)
  
  dat1 <- signal1 %>%
    select(geo_value, end = value)
  
  dat0 %>%
    left_join(dat1, by = "geo_value") %>%
    mutate(diff = end - start) %>%
    as_tibble()
}

process_data <- function(signal1, signal2, start_date, end_date) {
  deaths <- get_data(signal1, start_date, end_date)
  cases <- get_data(signal2, start_date, end_date)
  
  deaths <- deaths %>%
    select(geo_value, deaths = diff) %>%
    mutate(deaths_prop_cat = ntile(deaths, 4))
  
  cases <- cases %>%
    select(geo_value, cases = diff) %>%
    mutate(cases_prop_cat = ntile(cases, 4))
  
  full_data <- deaths %>%
    left_join(cases, by = "geo_value")
  
  saveRDS(full_data, paste0("../00_RawData/other/covidcast-usafacts-", end_date, ".rds"))
  full_data
}

signal1 <- "deaths_cumulative_prop"
signal2 <- "confirmed_cumulative_prop"

data.1220 <- process_data(signal1, signal2, "2020-11-30", "2020-12-31")
data.0121 <- process_data(signal1, signal2, "2020-12-31", "2021-01-31")
data.0221 <- process_data(signal1, signal2, "2021-01-31", "2021-02-28")
data.0321 <- process_data(signal1, signal2, "2021-02-28", "2021-03-31")
data.0421 <- process_data(signal1, signal2, "2021-03-31", "2021-04-30")

march <- readRDS("../00_RawData/other/covidcast-usafacts-2021-03-31.rds")
feb <- readRDS("../00_RawData/other/covidcast-usafacts-2021-02-28.rds")
