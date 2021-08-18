library(stringr)
source("01_ProcessData/survey-utils.R")
OUTFILE <- "../00_RawData/microdata/may/may-filtered.rds"

## input_data_monthly <- get_survey_df("monthly-files-1/")
input_data_monthly <- get_survey_df("monthly-files-2/")

cat("Got monthly data; GCing to proceed\n")
gc()

input_data_april <- get_survey_df("../00_RawData/microdata/may")
cat("Got May data; GCing to proceed\n")
gc()

## input_data <- input_data_monthly
input_data <- rbind(input_data_monthly, input_data_may)
input_data <- input_data_may

cat("Concatenated; proceeding with coding\n")
gc()

# To filter on the number of selected symptoms, we must know how many were
# possible to select in each survey wave. Note that this includes "None of the
# above", which Qualtrics makes mutually exclusive of other options, so it
# should not be possible to have all selected.
# Covers waves 1 - 10 (wave 9 is a lie)
max_symptoms <- Vectorize(function(wave) {
  waves <- c(15, 15, 16, 17, 19, 19, 19, 19, NA, 19)

  if (wave >= 1 && wave <= length(waves)) {
    return(waves[wave])
  }

  return(NA_real_)
})

# Similarly for comorbidities. Again, "None of the above" is exclusive of other
# options, so users should only be able to select all but one.
max_comorbidities <- Vectorize(function(wave) {
  waves <- c(9, 9, 9, 11, 11, 11, 11, 12, NA, 12)

  if (wave >= 1 && wave <= length(waves)) {
    return(waves[wave])
  }

  return(NA_real_)
})

## Adding check columns. These columns are all coded so that TRUE means the row
## "passes" and FALSE means it is anomalous or suspicious in some way.

# First, some numeric filters. Check that numeric variables are either blank or
# that reasonable values were added.
input_data <- input_data %>%
  mutate(#chk_temperature = is.na(Q40) | between(Q40, 96, 105),
         chk_hh_n_sick = is.na(A2) | between(A2, 0, 30),
         #chk_hh_size = is.na(A2b) | between(A2b, 0, 30),
         #chk_hh_size_sick = is.na(A2) | is.na(A2b) | A2 <= A2b,
         chk_work_contacts = is.na(C10_1_1) | between(C10_1_1, 0, 100),
         chk_shopping_contacts = is.na(C10_2_1) | between(C10_2_1, 0, 100),
         chk_social_contacts = is.na(C10_3_1) | between(C10_3_1, 0, 100),
         chk_other_contacts = is.na(C10_4_1) | between(C10_4_1, 0, 100))

# How many symptoms and comorbidities did the user select, out of how many were
# listed in this wave?
num_selected_symptoms <- str_count(input_data$B2, ",") + 1
max_selected_symptoms <- max_symptoms(input_data$wave)

num_selected_comorbidities <- str_count(input_data$C1, ",") + 1
max_selected_comorbidities <- max_comorbidities(input_data$wave)

# Filter on number of selected symptoms.
input_data <- input_data %>%
  mutate(chk_all_symptoms = is.na(B2) | num_selected_symptoms <= max_selected_symptoms - 1,
         chk_almost_all_symptoms = is.na(B2) | num_selected_symptoms <= max_selected_symptoms - 2)

# Filter on number of selected comorbidities.
input_data <- input_data %>%
  mutate(chk_all_comorbidities = is.na(C1) | num_selected_comorbidities <= max_selected_comorbidities - 1,
         chk_almost_all_comorbidities = is.na(C1) | num_selected_comorbidities <= max_selected_comorbidities - 2)

## Produce output data. We'll let users choose which subsetting they want to
## use, but provide a "chk_all" column that includes all the filters above.
## Filter on "chk_all" to get the most restricted subset.

input_data <- input_data %>%
  mutate(chk_all = (#chk_temperature & 
                      chk_hh_n_sick & 
                      #chk_hh_size &
                      #chk_hh_size_sick & 
                      chk_work_contacts &
                      chk_shopping_contacts & chk_social_contacts &
                      chk_other_contacts & chk_all_symptoms &
                      chk_almost_all_symptoms & chk_all_comorbidities &
                      chk_almost_all_comorbidities))

write_rds(input_data, OUTFILE)
input_data <- readRDS(OUTFILE)

## Diagnostics
# To examine the percent of responses that fail each check:
input_data %>%
  summarize(across(starts_with("chk_"),
                   function(x) {
                     100 * mean(!x)
                   })
            ) %>%
  as.list
