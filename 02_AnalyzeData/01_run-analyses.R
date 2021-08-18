source("02_AnalyzeData/trends-paper/analysis-funs.R")

# read data ----------------------------------------------------------------------------------
data <- readRDS("../00_RawData/processed-data/may-data-08-16-2021.rds") %>%
  filter(!is.na(intent_binary), gender != 4) %>%
  fastDummies::dummy_cols(c("rep_lead_ntile", "deaths_prop_cat"))

nonhesitant <- data %>%
  filter(hesitant == 0) %>%
  mutate(probably = ifelse(vaccine_intent == 2, 1, 0))

vars <- read_csv("../00_RawData/codebook/discrete-vars.csv") %>%
  select(groups = variable_name, value, label_full) %>%
  mutate_at("value", as.character)

codebook <- read_csv("../00_RawData/codebook/codebook-interactions.csv") 

fullvars <- read_csv("../00_RawData/codebook/discrete-vars-interactions.csv") 

# specify regression adjustment variables ----------------------------------------------------
all_groups <- rev(c("gender", "age_cat", "ethnicity_race", "educ", "emp_out", 
                "region", "urban_fips", "deaths_prop_cat", "rep_lead_ntile", "governor",
                "vaccine_flu", "covid_pos_ever", "health", "healthc", "elderly_risk", 
                "elderly_risk_det", "behav_worry", "behav_contact", "age_race", "age_race2"))

adjustment_set_1 <- c("gender", "age_race", "educ", "region",
                "emp_out", "behav_contact", "deaths_prop_cat", "healthc",
                "elderly_risk", "behav_worry", "vaccine_flu", "covid_pos_ever",
                "urban_fips", "rep_lead_ntile")

adjustment_set_2 <- c("gender", "age_race", "educ", "region",
                "emp_out", "behav_contact", "deaths_prop_cat", "health",
                "elderly_risk", "behav_worry", "vaccine_flu", "covid_pos_ever",
                "urban_fips", "rep_lead_ntile")

regression_variables <- map(list(adjustment_set_1, adjustment_set_2), generate_regression_variables, fullvars)

###############################################################################################
############################# estimate primary table data #####################################
###############################################################################################

# cross-tabs (intent by group percentages and risk ratios)
main_xtabs <- calc_xtabs(samp, all_groups[1:3], fullvars, "hesitant")
sens_xtabs <- calc_xtabs(nonhesitant, all_groups, fullvars, "probably")

saveRDS(main_xtabs, "02_Output/main-tables.rds")
saveRDS(sens_xtabs, "02_Output/sens-tables.rds")

# adjusted risk-ratios by group
main_models <- estimate_models(data, "hesitant", regression_variables, fullvars)
sens_models <- estimate_models(nonhesitant, "probably", regression_variables, fullvars)

all_model_results_v1 <- list(reg_results0 = main_models$output$reg_results0, 
                             reg_results1 = main_models$output$reg_results1,
                             reg_results0_sens = sens_models$output$reg_results0, 
                             reg_results1_sens = sens_models$output$reg_results1)

all_model_results_v2 <- list(reg_results0 = main_models$output$reg_results0.a, 
                             reg_results1 = main_models$output$reg_results1.a,
                             reg_results0_sens = sens_models$output$reg_results0.a, 
                             reg_results1_sens = sens_models$output$reg_results1.a)

write_csv(main_models$table$reg_results, "02_Output/model-results.csv")
write_csv(main_models$table$reg_results.a, "02_Output/model-results-a.csv")
write_csv(sens_models$table$reg_results, "02_Output/model-results-sens.csv")
write_csv(sens_models$table$reg_results.a, "02_Output/model-results-sens-a.csv")
saveRDS(all_model_results_v1, "02_Output/model-results-v1.rds")
saveRDS(all_model_results_v2, "02_Output/model-results-v2.rds")

###############################################################################################
############################# estimate age-race lincoms #######################################
###############################################################################################

# read data --------------------------------------------------------------------------
reg_results   <- readRDS("02_Output/model-results-v1.rds")
reg_results0 <- reg_results$reg_results0; reg_results1 <- reg_results$reg_results1

replacements <- fullvars$label_full
names(replacements) <- fullvars$orig_var

all_races <- filter(fullvars, orig_var == "ethnicity_race")$label_full
all_ages  <- filter(fullvars, orig_var == "D2")$label_full
all_ages  <- all_ages[-length(all_ages)]

# estimate linear combinations and CIs 
reg0_coefs <- reg_results0$model$coefs; reg0_vcov <- reg_results0$model$vcov
reg1_coefs <- reg_results1$model$coefs; reg1_vcov <- reg_results1$model$vcov

lincoms0 <- estimate_lincoms(reg0_coefs, reg0_vcov)
lincoms1 <- estimate_lincoms(reg1_coefs, reg1_vcov)
lincoms_unadj0 <- unadjusted_poisson_lincoms(data, fullvars, "age_race", "hesitant")

write_csv(lincoms0$lincom_white, "02_Output/rr-by-age-race.csv")
write_csv(lincoms0$lincom_young, "02_Output/arr-by-race.csv")
write_csv(lincoms1$lincom_white, "02_Output/rr-by-age-race-unvaccinated.csv")
write_csv(lincoms1$lincom_young, "02_Output/arr-by-race-unvaccinated.csv")
lincoms_unadj0$lincom_white_unadj %>%
  bind_rows(lincoms_unadj0$lincom_young_unadj) %>%
  write_csv("02_Output/age-race-lincoms-unadjusted.csv")

###############################################################################################
################################ estimate time trends #########################################
###############################################################################################

# preliminaries ---------------------------------------------------------------------

# read data
trend_data <- readRDS("../00_RawData/processed-data/jan-may-trendsdata.rds") %>%
  filter(!is.na(hesitant), gender != 4) %>%
  mutate(all = 1) %>%
  nest(-month) %>%
  mutate(ggmonth = lubridate::ydm(paste0("2021010", month)))

change_data <- trend_data %>%
  filter(month %in% c(1, 5)) %>%
  unnest(cols = c(data))

# change variable names/labels for plots & tables
codebook <- read_csv("../00_RawData/codebook/codebook.csv") %>%
  mutate_at("variable_name_full", ~gsub("Race and ethnicity",
                                        "Race and ethnicity (18 - 34 year olds)", .)) %>%
  mutate_at("variable_name_full", ~gsub("County Trump vote total minus Biden vote total among all votes",
                                        "Trump vote share", .))

fullvars <- read_csv("../00_RawData/codebook/discrete-vars.csv") %>%
  mutate_at("label_full", ~gsub("High school, GED, or less", "High school or less", .)) %>%
  mutate_at("label_full", ~gsub("Some college or two year degree", "Some college", .)) %>%
  mutate_at("label_full", ~gsub("Professional degree", "Professional \n degree", .))

# estimate time trends ---------------------------------------------------------------------
intent_vars <- c(grep("vaccine_intent_all_", names(change_data), value = TRUE)[-6])
educ_vars <- grep("educ_", names(change_data), value = TRUE)
vote_vars <- grep("rep_lead_ntile_", names(change_data), value = TRUE)
race_vars <- grep("ethnicity_race_", names(change_data), value = TRUE)
region_vars <- grep("region", names(change_data), value = TRUE)

all_trend <- trend_tables("all", trend_data, change_data) 
vacc_trend <- map(intent_vars, ~trend_tables("all", trend_data, change_data, outcome = .x)) %>%
  map2(c("Already vaccinated", "Yes, definitely", "Yes, probably",
         "No, probably not", "No, definitely not"), ~mutate(.x, label = .y)) %>%
  map(~set_names(.x, gsub("vaccine_intent_all_[1-5]", "hesitant", names(.)))) %>%
  invoke(rbind, .)
educ_trend <- trend_tables(educ_vars, trend_data, change_data)
vote_trend <- trend_tables(vote_vars, trend_data, change_data)
region_trend <- trend_tables(region_vars, trend_data, change_data)
race_trend <- trend_tables(race_vars, 
                           mutate(trend_data, data = map(data, ~filter(.x, age_cat2 == 1))),
                           filter(change_data, age_cat2 == 1))

aggregate_trends <- list(vacc_trend, educ_trend, vote_trend, race_trend, region_trend, all_trend) %>%
  map2(c("vaccine_intent", "educ", "rep_lead_ntile", "ethnicity_race", "region", "all"),
       ~mutate(.x, variable = .y)) %>%
  invoke(rbind, .)

saveRDS(aggregate_trends, "02_Output/trend-paper-hesitancy-trends.rds")
