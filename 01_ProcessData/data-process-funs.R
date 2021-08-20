# program: data-process-funs.R
# purpose: process microdata files to add/combine relevant variables for final analysis
# author: max rubinstein
# date modified: august 16, 2021

library(data.table)
library(tidyverse)
library(gtools)
library(usmap)
library(assertthat)

#########################################################################################
############################## data processing functions ################################
#########################################################################################

#' basic_process: takes raw microdata file and creates combined race/ethnicity variable,
#' primary vaccine intent variable, combined hesitancy reasons variable, combined age variable,
#' adds unemployed category to works outside or at home, combines some college with college and
#' less than HS with HS graduation, very worried with somewhat worried about COVID-19, unsure
#' to no for previous flu vaccine, and recodes ever tested positive for COVID variable. also adds
#' region data, urban classifications, Trump vote-share, and governor data (read from external
#' directories)
#'
#' @param data 
#'
#' @return dataframe with processed columns

basic_process <- function(data) {
  # codebook mapping survey variable names to new variable names
  recode <- read_csv("../00_RawData/codebook/discrete-vars.csv") %>%
    mutate_at("value", as.character) %>%
    filter(create_dummy == 1) %>%
    filter(!grepl("all", variable_name)) 
  
  codebook <- read_csv("../00_RawData/codebook/codebook.csv") %>%
    filter(original_file == 1) 

  # rename variables and recode NA values to 199
  replace_vars <- codebook$variable_name
  names(replace_vars) <- paste0("^", codebook$variable, "$")
  
  # urban/rural classifications
  county_urban <- readxl::read_excel("../00_RawData/other/NCHSURCodes2013.xlsx") %>%
    select(fips = `FIPS code`, urban_fips = `2013 code`) 
  
  # nytimes voting data 
  voting <- readRDS("../00_RawData/other/voting-agg.rds") %>%
    mutate_at("fips", as.numeric) %>%
    select(fips, pct_rep_lead, rep_lead_ntile)
  
  # us region data 
  r1 <- .midwest_region; r2 <- .south_region; r3 <- .pacific; r4 <- .mountain; r5 <- .northeast_region
  
  regions <- map2(list(1, 2, 3, 4, 5), 
                  list(r1, r2, r3, r4, r5),
                  ~tibble(region = .x, state = .y)) %>%
    invoke(rbind, .) %>%
    mutate(state = as.numeric(fips(state))) %>%
    bind_rows(
      tibble(
        state = c(60, 64, 66, 68, 69, 70, 72, 74, 78,
                  81, 84, 86, 67, 89, 71, 76, 95, 79)
      ) %>%
        mutate(region = 6)
    )
  
  # 2021 governor data 
  governors <- readxl::read_excel("../00_RawData/other/governors-2021.xlsx") %>%
    select(state = State, governor = Party) %>%
    mutate(state = as.numeric(fips(state))) %>%
    filter(!is.na(state)) %>%
    mutate(governor = ifelse(governor == "Republican", 1, 0))
  
  res <- data %>%
    mutate(mult_resp = ifelse(stringr::str_length(D7) > 1, 1, 0), # dummy for multiple race categories
           vaccine_intent_all = case_when(
             V1 == 1 ~ 1,
             V3 == 1 ~ 2,
             V3 == 2 ~ 3,
             V3 == 3 ~ 4,
             V3 == 4 ~ 5
           ),
           V3 = case_when(
             V1 == 1 ~ 1, 
             TRUE ~ as.numeric(V3)
           ),
           ethnicity_race = case_when(
             D6 == 1 ~ 1, # race/ethnicity recode; hispanic versus non-hispanic [race]
             (is.na(D6) | D6 == 2) & mult_resp == 0 & D7 == 1 ~ 2,
             (is.na(D6) | D6 == 2) & mult_resp == 0 & D7 == 2 ~ 3,
             (is.na(D6) | D6 == 2) & mult_resp == 0 & D7 == 3 ~ 4,
             (is.na(D6) | D6 == 2) & mult_resp == 0 & D7 == 4 ~ 5,
             (is.na(D6) | D6 == 2) & mult_resp == 0 & D7 == 5 ~ 6,
             (is.na(D6) | D6 == 2) & mult_resp == 0 & D7 == 6 ~ 7,
             (is.na(D6) | D6 == 2) & mult_resp == 1 ~ 8),
           ethnicity_race = ifelse(is.na(ethnicity_race), 7, ethnicity_race),
           age_cat2 = case_when(
             D2 %in% c(1, 2) ~ 1,
             D2 %in% c(3, 4, 5) ~ 2,
             D2 %in% c(6, 7) ~ 3,
             is.na(D2) ~ 199),
           D8  = case_when(D8 == 1 ~ 2, TRUE ~ as.numeric(D8)), # combine less than HS with HS graduation
           D8  = case_when(D8 == 4 ~ 3, TRUE ~ as.numeric(D8)), # combine some college with college
           C9  = case_when(C9 == 2 ~ 1, TRUE ~ as.numeric(C9)), # combine very worried with somewhat worried about COVID
           C17 = ifelse(C17 %in% c(3, 4), 2, C17), # combine unsure and no for flu vaccine 
           B11 = case_when(B11 == 3 ~ 2, # add unsure to not tested positive for COVID
                           B8 == 2 ~ 2, # add never tested for covid to no positive test
                           B10a == 1 ~ 1, # add positive test in last 14 days to positive test
                           TRUE ~ as.numeric(B11)),
           D10 = case_when(D9 == 2 ~ 3, TRUE ~ as.numeric(D10)) # add not working category for people employed outside of home
    ) %>%
    mutate(fips = as.numeric(fips)) %>%
    set_names(stringr::str_replace_all(names(.), replace_vars)) %>%
    mutate(fips.length = stringr::str_length(fips) - 4,
           state = as.numeric(stringr::str_sub(fips, start = 1, end = 1 + fips.length))) %>%
    mutate(gender = ifelse(gender == 5, 199, gender)) %>% # combine prefer not to answer with no response
    left_join(regions, by = "state") %>%
    left_join(governors, by = "state") %>%
    left_join(county_urban, by = "fips") %>%
    left_join(voting %>% select(fips, rep_lead_ntile), by = "fips")
  
  if (min(as.Date(res$StartDatetime)) >= as.Date("2021-02-01")) {
    res <- res %>%
      mutate(reasons = case_when(
        vaccine_intent == 2 ~ V5a,
        vaccine_intent == 3 ~ V5b,
        vaccine_intent == 4 ~ V5c
    ))
  }
  return(res)
}

#' final_process: take output from basic process and add additional variables (
#' including elderly risk, high risk health condition indicators), and recodes
#' missing age to people with race information but no age info. Also adds
#' indicator variables for all of the relevant categorical data. Outputs a new
#' version of the codebook and label cross-walk as a side-effect to include
#' age-race interactions, as well as a text file with the number of respondents
#' with present race information recoded to missing due to missing age information
#'
#' @param basic_data output from basic_process
#' @param covid_cases covid case data from prior month
#'
#' @return expanded dataframe with columns noted above
final_process <- function(basic_data, covid_cases) {
  
  # codebook mapping survey variable names to new variable names
  recode <- read_csv("../00_RawData/codebook/discrete-vars.csv") %>%
    mutate_at("value", as.character) %>%
    filter(create_dummy == 1) %>%
    filter(!grepl("all", variable_name)) 
  
  codebook <- read_csv("../00_RawData/codebook/codebook.csv") %>%
    filter(original_file == 1) 
  
  nas <- rep(list(199), length(codebook$variable_name))
  names(nas) <- codebook$variable_name
  
  reason_values <- read_csv("../00_RawData/codebook/hesitancy-reasons.csv") %>%
    mutate(varname = paste0("reason_", value))

  res <- basic_data %>%
    left_join(covid_cases, by = "fips") %>%
    mutate_at(vars(contains("prop_cat"), contains("num_cat")), funs(as.character(.))) %>%
    mutate(hh_eld = case_when(
      hh_eld >= 1 ~ 1,
      hh_eld == 0 ~ 0,
      hh_eld < 0 ~ NA_real_,
      TRUE ~ NA_real_
    )) %>%
    mutate(elderly_risk = case_when(
      age_cat %in% c(6, 7) | hh_eld == 1 ~ 1, # add elderly risk variables
      age_cat %in% c(1:5) & hh_eld == 0 ~ 0,
      TRUE ~ NA_real_)) %>%
    mutate(elderly_risk_det = case_when(
      age_cat %in% c(6, 7) ~ 1,
      age_cat %in% c(1:5) & hh_eld >= 1 ~ 2,
      age_cat %in% c(1:5) & hh_eld == 0 ~ 3,
      TRUE ~ NA_real_)) %>% 
    replace_na(nas) %>%
    mutate(intent_binary = case_when( # generate outcome variables
        vaccine_intent %in% c(1, 2) ~ 1,
        vaccine_intent %in% c(3, 4) ~ 0,
        is.na(vaccine_intent) ~ NA_real_), 
      hesitant = 1 - intent_binary,  # recode health conditions 
      num_conds = stringr::str_count(health,  "\\,") + 1, # number of health conditions variable
      num_conds = case_when(
        health == 199 ~ -1,
        health == 9 ~ 0,
        TRUE ~ num_conds
      ), 
      hbp = ifelse(grepl("4", health), 1, 0),  # high blood pressure indicator
      health2 = ifelse(num_conds == 2, gsub("4\\,|\\,4", "", health), health), 
      health1 = ifelse(num_conds > 1, 21, health), # create multiple conditions category
      health1 = ifelse(num_conds == 2 & hbp == 1, health2, health1), # recode to other if 2 conds and HBP = 1
      healthc = case_when(
        health == 9 ~ 1,
        num_conds >= 1 ~ 2,
        health == 199 ~ 199
      )) %>% 
    select(-hbp, -num_conds, -health) %>%
    rename(health = health1) %>%
    fastDummies::dummy_cols(c("gender", "age_cat", "age_cat2", "ethnicity_race",
                              "educ", "emp_out", "region", "urban_fips", "governor", 
                              "covid_pos_ever", "health", "healthc", "vaccine_flu", 
                              "elderly_risk", "elderly_risk_det",  "behav_worry", 
                              "behav_contact")) %>%
    generate_reason_indicators(reason_values$value) 
  
  # create age-race interaction variables and add to codebook
  codebook <- read_csv("../00_RawData/codebook/codebook.csv") %>%
    add_codebook_interaction("age_race", "Age by Race/Ethnicity") %>%
    add_codebook_interaction("age_race2", "Age (collapsed) by Race/Ethnicity")
  
  fullvars <- read_csv("../00_RawData/codebook/discrete-vars.csv") %>%
    add_xwalk_interaction("ethnicity_race", "age_cat", "age_race") %>%
    add_xwalk_interaction("ethnicity_race", "age_cat2", "age_race2") 
  
  #- quick check on number of NA-age with present race info
  subset1 <- filter(res, !is.na(hesitant), gender != 4)
  races <- filter(fullvars, variable_name == "ethnicity_race")
  repls <- races$label_full
  names(repls) <- races$value
  
  misstab <- table(subset1$ethnicity_race, subset1$age_cat2, useNA = "always") %>%
    as.data.frame() %>%
    filter(Var2 == 199) %>%
    mutate_at("Var1", ~stringr::str_replace_all(., repls)) %>%
    filter(!is.na(Var1)) %>%
    filter(!grepl("Unknown", Var1)) %>%
    as_tibble() 
  
  nummiss <- misstab %>%
    summarize_if(is.numeric, sum)
  
  # write the number of recoded present race to missing (because of missing age)
  write_csv(misstab, "C:/Users/mdrub/Dropbox/Delphi/01_Analysis/race-missing-age.csv")
  
  # recode missing age to missing race
  res <- res %>%
    mutate(age_race = paste0(age_cat, "_", ethnicity_race), # ethnicity by race
           age_race2 = paste0(age_cat2, "_", ethnicity_race), # ethnicity by race
           age_race = ifelse(grepl("199", age_race), "199", age_race), # recode missing values (if any = missing)
           age_race2 = ifelse(grepl("199", age_race2), "199", age_race2)) %>% # recode missing values (if any = missing)
    fastDummies::dummy_cols(select_columns = c("age_race", "age_race2"))
  
  write_csv(codebook, "../00_RawData/codebook/codebook-interactions.csv")
  write_csv(fullvars, "../00_RawData/codebook/discrete-vars-interactions.csv")
  
  return(res)
}

##############################################################################################
################################### helper functions #########################################
##############################################################################################

#' generate_reason_indicators: create dummy variables for hesitancy reasons
#'
#' @param data output from basic_process
#' @param reason_values list of values and value names of hesitancy reasons
#'
#' @return dataframe with indicators for hesitancy reasons
generate_reason_indicators <- function(simple_data, reason_values) {
  generate_indicator <- function(data, reason_value) {
    varname <- paste0("reason_", reason_value)
    value_string <- paste0("^", reason_value, "$|", "\\,", reason_value, "\\,|\\,", 
                           reason_value, "$|^", reason_value, "\\,")
    
    res <- mutate(data, !!varname := ifelse(grepl(value_string, reasons), 1, 0))
    return(res)
  }
  
  for(reason_value in reason_values) {
    simple_data <- generate_indicator(simple_data, reason_value)
  }
  return(simple_data)
}

#' add_xwalk_interaction: adds a list of interactions to the variable value
#' and description cross-walk file 
#'
#' @param label_xwalk cross-walk for variable label values and full descriptions
#' @param var1 first variable in interaction
#' @param var2 second variable in interaction
#' @param newname name for interacted variables
#' @param start defunct
#' @param missing_label what to relabel missing responses
#'
#' @return cross-walk with rows for interactions added
add_xwalk_interaction <- function(label_xwalk, var1, var2, newname, 
                            start = "", missing_label = "No response (either category)") {
  suffix_gen <- function(label_xwalk, var1, var2) {
    t1 <- list()
    for(var in label_xwalk$value[label_xwalk$variable_name == var1]) {
      t1 <- append(t1, map(label_xwalk$value[label_xwalk$variable_name == var2], 
                           ~paste0(.x, "_", var)))
    }
    t1
  }
  
  num_rep1 <- length(label_xwalk$label_abbrev[label_xwalk$variable_name == var1])
  num_rep2 <- length(label_xwalk$label_abbrev[label_xwalk$variable_name == var2])
  
  refs <- label_xwalk %>%
    filter(variable_name %in% c(var1, var2),
           reference == 1)
  
  ref.val <- paste0(refs$value[refs$variable_name == var2], "_",
                    refs$value[refs$variable_name == var1])  
  
  interaction <- tibble(
    value = unlist(suffix_gen(label_xwalk, var1, var2)),
    label_full1 = rep(label_xwalk$label_full[label_xwalk$variable_name == var1], each = num_rep2),
    label_full2 = rep(label_xwalk$label_full[label_xwalk$variable_name == var2], num_rep1),
    label_full = paste0(label_full1, ": ", label_full2),
    variable_name = newname,
    reference = ifelse(value == ref.val, 1, 0),
    new_varname = paste0(variable_name, "_", value),
    create_dummy = 0,
    label_abbrev = "none",
    orig_var = new_varname
  ) %>%
    select(-label_full1, -label_full2)
  
  interaction <- interaction %>%
    filter(!grepl(paste0("199", start), orig_var))
  
  interaction <- interaction %>%
    bind_rows(
      tibble(
        orig_var = paste0(newname, "_199"),
        variable_name = newname,
        value = "199",
        label_abbrev = "none",
        label_full = missing_label,
        new_varname = orig_var,
        reference = 0,
        create_dummy = 0
      )
    )
  
  label_xwalk <- label_xwalk %>%
    mutate(value = as.character(value)) 
  
  bind_rows(label_xwalk, interaction)
}

#' add_codebook_interaction: adds interaction for two variables to the
#' variable codebook
#'
#' @param codebook codebook containing variable name and coded crosswalk
#' @param newname new variable to add to codebook
#' @param description description of variable to add to codebook
#'
#' @return codebook with additional variable
add_codebook_interaction <- function(codebook, newname, description) {
  newrow <- tibble(
    variable_name = newname,
    variable = newname,
    description = description,
    variable_name_full = description,
    original_file = 0,
    create_dummies = 0,
    order = NA
  )
  bind_rows(codebook, newrow)
}

