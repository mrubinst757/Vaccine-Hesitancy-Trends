# program: aggregate-microdata.R
# purpose: take daily microdata aggregates, combine into one file
# author: max rubinstein
# date modified: may 5, 2021

# read libraries ----------------------------------------------------------------------------------------
source("01_RPrograms/trends-paper/01_ProcessData/data-process-funs.R")

# create trend data file --------------------------------------------------------------------------
date_start <- c("2021-01-06", "2021-02-01", "2021-03-01", "2021-04-01", "2021-05-01")
date_end   <- c("2021-02-01", "2021-03-01", "2021-04-01", "2021-05-01", "2021-06-01")
date_args  <- transpose(list(date_start, date_end))

trend_data <- list.files("../00_RawData/microdata/", pattern = "race", full.names = TRUE) %>%
  map(~data.table::fread(.x) %>%
        mutate(hesitant = case_when(
          V1 == 1 | V3 %in% c(1:2) ~ 0,
          V3 %in% c(3:4) ~ 1,
          TRUE ~ NA_real_)) %>%
        filter(wave <= 10) %>%
        set_names(gsub("C17a", "C17", names(.)))) %>%
  map2(date_args, ~filter(.x, as.Date(StartDatetime) >= .y[[1]], as.Date(StartDatetime) < .y[[2]])) %>%
  map(basic_process) %>%
  map(~dplyr::select(.x, weight, age_cat, age_cat2, ethnicity_race, region, 
                     educ, gender, emp, emp_out, rep_lead_ntile,
                     start_date = StartDatetime, hesitant, vaccine_intent_all)) %>%
  invoke(rbind, .) %>%
  replace_na(list(rep_lead_ntile = 199, educ = 199, gender = 199, region = 199)) %>%
  mutate(month = lubridate::month(start_date)) %>%
  fastDummies::dummy_cols(c("rep_lead_ntile", "ethnicity_race", "region", "educ", "vaccine_intent_all"))

saveRDS(trend_data, "../00_RawData/processed-data/jan-may-trendsdata.rds")

# create may analytic file -------------------------------------------------------------------------

# april cases/deaths classification
covid_cases <- readRDS("../00_RawData/other/covidcast-usafacts-2021-04-30.rds") %>%
  rename(fips = geo_value) %>%
  mutate_at("fips", as.numeric) %>%
  fastDummies::dummy_columns(c("cases_prop_cat", "deaths_prop_cat"))

# list of new files given an input data (y-m-d format)
may_data <- readRDS("../00_RawData/microdata/may/may-filtered.rds") %>%
  rename(C17 = C17a) %>%
  filter(wave == 10) %>%
  basic_process() %>%
  final_process(covid_cases)

saveRDS(may_data, "../00_RawData/processed-data/may-data-08-16-2021.rds")
