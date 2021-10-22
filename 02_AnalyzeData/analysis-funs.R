library(dplyr)
library(readr)
library(sandwich)
library(msm)
library(gt)
library(purrr)
library(R.utils)

#' format_cells: formats numeric values from a table to create cells with
#' an estimate (lower CI, upper CI) with specified number of digits and
#' possibly rescaled
#'
#' @param num point estimate
#' @param lci lower confidence interval
#' @param uci upper confidence interval
#' @param digits how many digits to round to
#' @param mult multiplicative scaling of original values
#'
#' @return formatted numbers 
format_cells <- function(num, lci, uci, digits, mult = 1) {
  paste0(format(round(mult*num, digits), nsmall = digits), " (", 
         format(round(mult*lci, digits), nsmall = digits), ", ", 
         format(round(mult*uci, digits), nsmall = digits), ")")
}

#####################################################################################
########################### primary table functions #################################
#####################################################################################

#' generate_regression_variables: takes a categorical variable and retrieves
#' corresponding indicator variables
#'
#' @param adjustment_set list of categorical variables to be used for regression adjustment
#' @param crosswalk crosswalk from categories to variable labels
generate_regression_variables <- function(adjustment_set, crosswalk) {
  crosswalk %>%
    filter(variable_name %in% adjustment_set) %>%
    filter(reference != 1) %>%
    .$new_varname
}


# intent_by_category: calculates cross-tabs with robust standard errors for a given
#'categorical variable
#'
#' @param data analytic dataset
#' @param crosswalk crosswalk from categories to variable labels 
#' @param group categorical variable
#' @param outcome binary outcome (typically hesitancy)
intent_by_category <- function(data, crosswalk, group, outcome) {
  gvars <- crosswalk %>%
    filter(variable_name == !!group) 
  
  totals <- unlist(map(gvars$new_varname, ~nrow(filter(data, !!sym(.x) == 1))))
  
  formula <- map(gvars$new_varname, ~paste0(outcome, " ~ ", .x,  "-1"))
  
  models <- map(formula, ~lm(as.formula(.x), data = data, weights = data$weight))
  robust <- map(models, ~vcovHC(.x, type = "HC0"))
  mcoefs <- unlist(map(models, ~coef(.x)))
  se_est <- unlist(map(robust, ~try(sqrt(.x[1,1]))))
  
  tibble(
    props = mcoefs,
    label = gvars$new_varname,
    se = se_est,
    total = totals
  )
}

# rr_by_category: calculates risk ratios with robust standard errors for a given
#' categorical variable
#' 
#' @param data analytic dataset
#' @param crosswalk variable label and description crosswalk 
#' @param group categorical variable
#' @param outcome binary outcome (typically hesitancy)
#'
#' @return
#' @export
#'
#' @examples
rr_by_category <- function(data, crosswalk, group, outcome) {
  gvars <- crosswalk %>%
    filter(variable_name == !!group) %>%
    filter(reference == 0)
  
  formula <- paste0(outcome, " ~ ", paste0(gvars$new_varname, collapse = "+"))
  
  rr_model <- glm(as.formula(formula), data = data, weights = data$weight, family = poisson)
  rr_cov <- vcovHC(rr_model, type = "HC0")
  coefs <- coef(rr_model)
  coefs.na <- coefs[!is.na(coefs)]
  my_coefs <- map(sprintf("~exp(x%s)", 2:length(coefs.na)), as.formula)
  my_coefs.na <- my_coefs[!is.na(my_coefs)]
  se_ests <- deltamethod(my_coefs.na, coefs.na, rr_cov)
  if (length(grep(TRUE, is.na(coefs))) > 0) {
    se_ests <- insert(se_ests, grep(TRUE, is.na(coefs)) - 1)
  }
  
  tibble(
    rr = exp(coef(rr_model))[-1],
    label = names(coef(rr_model)[-1]),
    se = se_ests
  )
}

#' calc_xtabs: calculate cross-tabs and risk ratios for a list of categorical variables conducted
#' on primary dataset and filtered to only unvaccinated individuals (not used for analysis)
#'
#' @param data analytic file
#' @param groups list of variables to calculate cross-tabs by
#' @param crosswalk variable label and description crosswalk
#' @param outcome binary outcome (typically hesitancy)
calc_xtabs <- function(data, groups, crosswalk, outcome) {
  rd_table_v1 <- map(groups, ~intent_by_category(data, crosswalk, .x, outcome)) 
  rd_table_v2 <- map(groups, ~intent_by_category(filter(data, vaccine_ever != 1), crosswalk, .x, outcome)) 
  rr_table_v1 <- map(groups, ~rr_by_category(data, crosswalk, .x, outcome))
  rr_table_v2 <- map(groups, ~rr_by_category(filter(data, vaccine_ever != 1), crosswalk, .x, outcome)) 
  list(rd_table_v1, rd_table_v2, rr_table_v1, rr_table_v2)
}

#' adjusted_risk_ratios: calculate adjusted risk-ratios using weighted
#' Poisson regression and robust standard errors 
#'
#' @param data analytic file
#' @param outcome binary outcome (typically hesitancy)
#' @param adjustment_variables 
#' @param crosswalk variable category and label cross-walk
adjusted_risk_ratios <- function(data, outcome, adjustment_variables, crosswalk) {
  my_formula <- paste0(outcome, " ~ ", paste0(adjustment_variables, collapse = "+"))
  rr_model <- glm(as.formula(my_formula), family = poisson, data = data, weights = weight)
  coef.exp <- exp(coef(rr_model))
  rr_cov <- vcovHC(rr_model, type = "HC0")
  coef.exp.na <- coef.exp[!is.na(coef.exp)]
  my_coefs <- map(sprintf("~exp(x%s)", 2:length(coef.exp.na)), as.formula)
  se_ests <- deltamethod(my_coefs, coef(rr_model)[!is.na(coef(rr_model))], rr_cov)
  
  reg_results <- tibble(
    rr.adj = coef.exp[-1][!is.na(coef.exp[-1])],
    se.adj = se_ests,
    l95ci = rr.adj - 1.96*se_ests,
    u95ci = rr.adj + 1.96*se_ests,
    ci.adj = format_cells(rr.adj, l95ci, u95ci, 2, 1),
    new_varname = adjustment_variables[!is.na(coef.exp[-1])]
  ) %>%
    left_join(crosswalk, by = c("new_varname")) 
  
  list(table = reg_results, model = list(coefs = coef(rr_model), vcov = rr_cov))
}

estimate_models <- function(data, outcome, regression_variable_list, crosswalk) {
  # table output as CSV file
  table_bind <- function(reg_results0, reg_results1) {
    reg_results0$table %>%
      mutate(sample = "full") %>%
      bind_rows(
        reg_results1$table %>%
          mutate(sample = "unvaccinated")
      )
  }
  
  reg_results0   <- adjusted_risk_ratios(data, outcome, regression_variable_list[[1]], crosswalk)
  reg_results1   <- adjusted_risk_ratios(subset(data, vaccine_ever != 1), outcome, 
                                       regression_variable_list[[1]], crosswalk)
  reg_results0.a <- adjusted_risk_ratios(data, outcome, regression_variable_list[[2]], crosswalk)
  reg_results1.a <- adjusted_risk_ratios(subset(data, vaccine_ever != 1), outcome, 
                                       regression_variable_list[[2]], crosswalk)
  reg_results   <- table_bind(reg_results0, reg_results1)
  reg_results.a <- table_bind(reg_results0.a, reg_results1.a)
  
  list(tables = list(reg_results = reg_results, reg_results.a = reg_results.a),
       output = list(reg_results0 = reg_results0, reg_results1 = reg_results1, 
                     reg_results0.a = reg_results0.a, reg_results1.a = reg_results1.a))
}


######################################################################################
######################## linear combination functions ################################
######################################################################################

# within year linear combinations

#' lincom_within_year: calculates linear combinations to compare risk ratios
#' across racial groups within an age category
#'
#' @param years list of age groups
#' @param ref_cat reference category for racial group
#' @param coefs model coefficients from poisson regression
#' @param vcov covariance matrix from poisson regression
#'
#' @return dataframe of risk ratios among racial categories within
#' an age group
lincom_within_year <- function(years, ref_cat, coefs, vcov) {
  # calculate within year linear combination
  lincom <- function(year_cat, comp_cat, ref_cat, coefs, vcov) {
    coef_names <- names(coefs)
    
    ref_group <- paste0(ref_cat, ": ", year_cat)
    comp_group <- paste0(comp_cat, ": ", year_cat)
    
    ref_num1 <- grep(ref_group, coef_names, fixed = TRUE)
    comp_num1 <- grep(comp_group, coef_names, fixed = TRUE)
    
    rr_1 <- exp(coefs[comp_num1] - coefs[ref_num1])
    
    ref_num <- paste0("x", grep(ref_group, coef_names, fixed = TRUE))
    comp_num <- paste0("x", grep(comp_group, coef_names, fixed = TRUE))
    
    if (comp_num != "x" & ref_num != "x") {
      form <- paste0("~exp(", comp_num, "-", ref_num, ")")
    }
    if(comp_num == "x") {
      ref_num1 <- grep(ref_group, coef_names, fixed = TRUE)
      rr_1 <- exp(coefs[ref_num1])
      
      form <- paste0("~exp(", ref_num, ")")
    }
    if(ref_num == "x") {
      comp_num1 <- grep(comp_group, coef_names, fixed = TRUE)
      rr_1 <- exp(coefs[comp_num1])
      
      form <- paste0("~exp(", comp_num, ")")
    }
    if(!(ref_num == "x" & comp_num == "x")) {
      se_ests <- msm::deltamethod(as.formula(form), coefs, vcov)
      
      final <- tibble(
        year_cat = year_cat,
        ref_cat = ref_cat,
        comp_cat = comp_cat,
        rr = rr_1,
        se = se_ests,
        lci = rr - 1.96*se,
        uci = rr + 1.96*se
      )
    }
    if((ref_num == "x" & comp_num == "x")) {
      
      final <- tibble(
        year_cat = year_cat,
        ref_cat = ref_cat,
        comp_cat = comp_cat,
        rr = 1,
        se = 0,
        lci = rr - 1.96*se,
        uci = rr + 1.96*se
      )
    }
    final
  }
  
  # iterate across all race categories
  map(all_races, ~lincom(years, .x, ref_cat, coefs, vcov)) %>%
    invoke(rbind, .)
}

#' lincom_within_race: calculates linear combinates to compare risk ratios
#' across age categories within a racial group
#'
#' @param race list of racial groups
#' @param ref_cat reference category for age
#' @param coefs model coefficients from poisson regressoin
#' @param vcov covariance matrix from poisson regression
#'
#' @return dataframe of risk ratios among age categories within a
#' list of racial groups
lincom_within_race <- function(race, ref_cat, coefs, vcov) {
  # calculate within racial group linear combination
  lincom <- function(year_ref, year_comp, race, coefs, vcov) {
    coef_names <- names(coefs)
    ref_group <- paste0(race, ": ", year_ref)
    comp_group <- paste0(race, ": ", year_comp)
    
    ref_num1 <- grep(ref_group, coef_names, fixed = TRUE)
    comp_num1 <- grep(comp_group, coef_names, fixed = TRUE)
    
    rr_1 <- exp(coefs[comp_num1] - coefs[ref_num1])
    
    ref_num <- paste0("x", grep(ref_group, coef_names, fixed = TRUE))
    comp_num <- paste0("x", grep(comp_group, coef_names, fixed = TRUE))
    
    if (comp_num != "x" & ref_num != "x") {
      form <- paste0("~exp(", comp_num, "-", ref_num, ")")
    }
    if(comp_num == "x") {
      ref_num1 <- grep(ref_group, coef_names, fixed = TRUE)
      rr_1 <- exp(coefs[ref_num1])
      
      form <- paste0("~exp(", ref_num, ")")
    }
    if(ref_num == "x") {
      comp_num1 <- grep(comp_group, coef_names, fixed = TRUE)
      rr_1 <- exp(coefs[comp_num1])
      
      form <- paste0("~exp(", comp_num, ")")
    }
    if(!(ref_num == "x" & comp_num == "x")) {
      se_ests <- msm::deltamethod(as.formula(form), coefs, vcov)
      
      final <- tibble(
        race = race,
        ref_cat = year_ref,
        comp_cat = year_comp,
        rr = rr_1,
        se = se_ests,
        lci = rr - 1.96*se,
        uci = rr + 1.96*se
      )
    }
    if((ref_num == "x" & comp_num == "x")) {
      
      final <- tibble(
        race = race,
        ref_cat = year_ref,
        comp_cat = year_comp,
        rr = 1,
        se = 0,
        lci = rr - 1.96*se,
        uci = rr + 1.96*se
      )
    }
    final
  }
  
  # iterate across all age categories
  map(all_ages, ~lincom(ref_cat, .x, race, coefs, vcov)) %>%
    invoke(rbind, .)
}

#' estmate_lincoms: calculate linear combinations of age/race interactions
#' from model output to compare risk ratios within age groups and racial
#' categories
#'
#' @param coefs coefficient vector from poisson regression model (with covariate adjustment)
#' @param vcov covariance matrix from poisson regression model (with covariate adjustment)
#'
#' @return list of linear combinations using white and 65-74 years of
#' age as the reference categories
estimate_lincoms <- function(coefs, vcov) {
  coef_names <- stringr::str_replace_all(names(coefs), replacements)
  coef_names <- coef_names[!is.na(coefs)]
  fcoefs <- coefs[!is.na(coefs)]
  names(fcoefs) <- coef_names
  
  lincom_white <- map(all_ages, ~lincom_within_year(.x, "White", fcoefs, vcov)) %>%
    invoke(rbind, .)
  
  lincom_young <- map(all_races, ~lincom_within_race(.x, "65-74 years", fcoefs, vcov)) %>%
    invoke(rbind, .)
  
  list(lincom_white = lincom_white,
       lincom_young = lincom_young)
}

#' unadjusted_poisson_lincoms: calculates the linear combinations for the age-race
#' interaction terms without a covariate adjustment, using white and 65-74 years old
#' as the reference categories 
#'
#' @param data analytic dataset
#' @param crosswalk variable category descriptions 
#' @param group categorical variable
#' @param outcome binary outcome (typically hesitancy)
#'
#' @return unadjusted linear combinations of the model coefficients for
#' age and race interaction terms
unadjusted_poisson_lincoms <- function(data, fullvars, group, outcome) {
  gvars <- crosswalk %>%
    filter(variable_name == !!group) %>%
    filter(reference == 0)
  
  formula <- paste0(outcome, " ~ ", paste0(gvars$new_varname, collapse = "+"))
  
  model <- glm(as.formula(formula), data = data, weights = data$weight, family = poisson)
  
  rr_cov1 <- vcovHC(model, type = "HC0")
  fcoefs1 <- coef(model)[!is.na(coef(model))]
  
  coef_names <- stringr::str_replace_all(names(coef(model)), replacements)
  coef_names <- coef_names[!is.na(coef(model))]
  names(fcoefs1) <- coef_names
  
  # estimate linear combinations and CIs
  lincom_white_unadj <- map(all_ages, ~iter_lincom(.x, "White", fcoefs1, rr_cov1)) %>%
    invoke(rbind, .) %>%
    mutate(group = "By Race")
  
  lincom_young_unadj <- map(all_races, ~iter_lincom1(.x, "65-74 years", fcoefs1, rr_cov1)) %>%
    invoke(rbind, .) %>%
    mutate(group = "By Age")
  
  list(lincom_white_unadj = lincom_white_unadj,
       lincom_young_unadj = lincom_young_unadj)
}

#########################################################################################
############################# time trend analysis functions #############################
#########################################################################################

#' trend_tables: create time trends for specified list of binary variables, 
#' including both monthly estimates and changes between specified months
#' 
#' @param var_list list of binary variables (typically occupations)
#' @param monthly_data dataframe nested by month
#' @param change_data dataframe containing only two months to compare change 
#' over (typically January and May)
#' @param outcome binary outcome (typically vaccine hesitancy)
#'
#' @return dataframe containing the percent within the outcome for each month and the
#' change over the same time perios
trend_tables <- function(var_list, monthly_data, change_data, outcome = "hesitant") {
  percent_changes <- map(var_list, ~calculate_intent(change_data, .x, outcome, change = TRUE)) %>%
    invoke(rbind, .) 
  
  percent_by_month <- monthly_data %>%
    mutate(hesitant = map(data, ~calculate_intent_loop(.x, var_list, outcome))) %>%
    select(month, hesitant) %>%
    mutate(hesitant = map(hesitant, ~invoke(rbind, .x))) %>%
    unnest(cols = c(hesitant)) %>%
    bind_rows(percent_changes) 
  
  return(percent_by_month)
}

#' calculate_intent: calculates the percentage of a binary outcome within a specified strata,
#' or calculates the change in the outcome between two months (outcome typically vaccine
#' intent)
#'
#' @param data unprocessed microdata
#' @param group binary variable that species a subset of the dataframe
#' @param outcome binary outcome
#' @param change boolean indicating whether to calculate an intercept or a change
#' across two months
#'
#' @return dataframe containing proportion and standard errors associated with
#' the outcome within the specified group. if change is true also returns the 
#' estimated percent change and standard error
calculate_intent <- function(data, group, outcome, change = FALSE) {
  data <- data %>%
    filter(!!sym(group) == 1)
  
  if (change == FALSE) {
    formula <- paste0(outcome, " ~ 1")
  }
  
  if (change == TRUE) {
    formula <- paste0(outcome, " ~ factor(month)")
  }
  
  model <- lm(as.formula(formula), data = data, weights = data$weight)
  robust <- vcovHC(model, type = "HC0")
  
  if (change == FALSE) {
    mcoefs <- coef(model)
    se_est <- sqrt(robust[1,1])
    res <- tibble(
      props = mcoefs,
      se = se_est,
      lci = props - 1.96*se,
      uci = props + 1.96*se
    ) %>%
      set_names(paste0(outcome, "_", names(.))) %>%
      mutate(label = group)
  }
  
  if (change == TRUE) {
    mcoefs <- coef(model)[2]
    se_est <- sqrt(robust[2,2])
    
    calc_ratio <- function(model, vcov) {
      coef1 <- coef(model)[1]; coef2 <- coef(model)[2]
      ratio.est <- coef2/coef1
      ratio.var <- ratio.est^2 * (vcov[1, 1]/(coef1^2) + vcov[2, 2] / (coef2^2) - 2 * vcov[1, 2] / (coef1 * coef2))
      c(ratio.est = ratio.est, se.est = sqrt(ratio.var))
    }
    ratios <- calc_ratio(model, robust)
    res <- tibble(
      props = c(mcoefs, ratios[1]),
      se = c(se_est, ratios[2]),
      lci = props - 1.96*se,
      uci = props + 1.96*se
    ) %>%
      set_names(paste0(outcome, "_", names(.))) %>%
      mutate(label = group, month = c(6, 7))
  }
  return(res)
}

#' calculate_intent_loop: iterates calculate_intent function across a list of variables
#'
#' @param data analytic dataset
#' @param var_list list of binary variables
#' @param outcome binary outcome
#'
#' @return list of dataframes of calc_intent output
calculate_intent_loop <- function(data, var_list, outcome) {
  map(var_list, ~calculate_intent(data, .x, outcome))
}

#######################################################################################
########################### hesitancy reason estimates ################################
#######################################################################################

#' @description estimate percent with a specified hesitancy reason
#'
#' @param data analytic dataset
#' @param outcome hesitancy reason (binary)
#'
#' @return dataframe with percent with a reason and standard error
hesitancy_reason_estimates <- function(data, outcome) {
  my_formula <- paste0(outcome, " ~ 1")
  m1 <- lm(as.formula(my_formula), weight = weight, data = data)
  se <- sqrt(vcovHC(m1, type = "HC0")[1,1])
  tibble(
    estimate = 100*coef(m1),
    se = 100*se
  )
}

#' @description create a table of hesitancy reasons labelled by
#' the subgroup analyzed
#'
#' @param data analytic dataset (typically subsetted to an occupation 
#' category of interest)
#' @param label label for dataset
#'
#' @return formatted dataframe of hesitanc reasons
reasons_table <- function(data, label) {
  num = nrow(data)
  
  map(reason_values$varname, ~hesitancy_reason_estimates(data, .x)) %>%
    invoke(rbind, .) %>%
    mutate(Reason = reason_values$description) %>%  
    mutate(l95ci = estimate - 1.96*se, u95ci = estimate + 1.96*se,
           reason_no = format_cells(estimate, l95ci, u95ci, 1, 1)) %>%
    arrange(-estimate) %>%
    select(Reason, reason_no) %>%
    rename(!!paste0(label, " \n N = ", num) := reason_no)
}

#' sample_characteristics: takes analytic file and outputs summary sample 
#' characteristics (in contrast to sample_characteristics_alll), removing 
#' non-response from percentages
#'
#' @param data analytic file
#'
#' @return dataframe containing relevant sample characteristics
sample_characteristics <- function(data) {
  # age stats
  age_dat <- filter(data, age_cat != 199)
  age <- matrixStats::weightedMedian(age_dat$age_cat, age_dat$weight)
  age_75 <- weighted.mean(age_dat$age_cat == 7, age_dat$weight)
  
  # gender stats
  gender_dat <- filter(data, gender != 199)
  male   <- weighted.mean(gender_dat$gender == 1, gender_dat$weight)
  female <- weighted.mean(gender_dat$gender == 2, gender_dat$weight)
  nonb   <- weighted.mean(gender_dat$gender == 3, gender_dat$weight)
  
  # race stats
  race_dat <- filter(data, ethnicity_race != 7)
  hisp  <- weighted.mean(race_dat$ethnicity_race == 1, race_dat$weight) 
  natam <- weighted.mean(race_dat$ethnicity_race == 2, race_dat$weight) 
  asian <- weighted.mean(race_dat$ethnicity_race == 3, race_dat$weight) 
  black <- weighted.mean(race_dat$ethnicity_race == 4, race_dat$weight) 
  pacis <- weighted.mean(race_dat$ethnicity_race == 5, race_dat$weight) 
  white <- weighted.mean(race_dat$ethnicity_race == 6, race_dat$weight) 
  mult   <- weighted.mean(race_dat$ethnicity_race == 8, race_dat$weight) 
  
  #educ stats
  educ_dat <- filter(data, educ != 199)
  lshs <- weighted.mean(educ_dat$educ == 2, educ_dat$weight)
  college <- weighted.mean(educ_dat$educ %in% c(5:8), educ_dat$weight)
  phd <- weighted.mean(educ_dat$educ == 7, educ_dat$weight)
  
  #emp dat
  emp_dat <- filter(data, emp != 199)
  worked <- weighted.mean(emp_dat$emp == 1, emp_dat$weight)
  
  emp_out_dat <- filter(data, emp_out != 199)
  worked_outside <- weighted.mean(emp_out_dat$emp_out == 1, emp_out_dat$weight)
  
  tibble(
    `Median age` = age,
    `75 or older` = 100*age_75,
    `% male` = 100*male,
    `% female` = 100*female,
    `% non-binary` = 100*nonb,
    `% white` = 100*white,
    `% hispanic` = 100*hisp,
    `% black` = 100*black,
    `% asian` = 100*asian,
    `% native american` = 100*natam,
    `% pacific islander` = 100*pacis,
    `% multiracial` = 100*mult,
    `% less than HS` = 100*lshs,
    `% college or more` = 100*college,
    `% phd` = 100*phd,
    `% worked for pay` = 100*worked,
    `% worked outside` = 100*worked_outside
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    set_names(c("Variable", "Value")) %>%
    mutate_at("Value", ~round(., 2))
}

#' flow_table: takes full data and outputs participant flow table (ie how
#' many respondents met inclusion criteria for our analyses)
#'
#' @param data analytic file
#' @param month_num number specifying targeted month
#'
#' @return dataframe containing relevant sample flow

sample_flow <- function(trend_data, month_data) {

  defyes_may <- month_data %>%
    filter(!is.na(hesitant) & gender != 4, vaccine_intent_all %in% c(1:2)) %>%
    nrow()
  
  na_reasons <- month_data %>%
    filter(vaccine_intent_all %in% c(3:5) & gender != 4) %>%
    filter(is.na(reasons)) %>%
    nrow()
  
  reasons_samp <- month_data %>%
    filter(vaccine_intent_all %in% c(3:5) & gender != 4) %>%
    filter(!is.na(reasons)) %>%
    nrow()
  
  excluded <- trend_data %>%
    nest(-month) %>%
    mutate(original = map(data, nrow)) %>%
    mutate(missing_y = map(data, ~mutate(.x, missing_y = is.na(hesitant)) %>%
                             .$missing_y %>% sum())) %>%
    mutate(selfdesc = map(data, ~mutate(.x, selfdesc = ifelse(gender == 4, 1, 0)) %>%
                            filter(!is.na(hesitant)) %>%
                            .$selfdesc %>% sum()),
           final = map(data, ~nrow(filter(.x, gender != 4, !is.na(hesitant))))) %>%
    select(original, missing_y, selfdesc, final)
  
  table_one <- excluded %>%
    set_names(c("Original N", "Missing Intent", "Self-describe", "Final")) %>%
    as_tibble() %>%
    unnest() %>%
    mutate(Month = c("January", "February", "March", "April", "May")) %>%
    cbind(
      tibble(
        `Vaccinated/DefYes` = c(NA, NA, NA, NA, defyes_may),
        `No reasons provided` = c(NA, NA, NA, NA, na_reasons),
        `Reasons sample` = c(NA, NA, NA, NA, reasons_samp)
      )
    ) 
  
  rows <- colnames(table_one)
  cols <- table_one$Month
  
  table_one %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    set_names(c(" ", cols)) %>%
    mutate(` ` = c("Began survey", "Did not report hesitancy", "Self-describe gender", 
                   "Report sample", "Month", "Vaccinated/Definitely Yes", "Did not provide reasons", 
                   "Reasons sample")) %>%
    filter(` ` != "Month") %>%
    gt()
  
  print(table_one)
}

#' sample_chars_all: calculates percentage of observations within each demographic
#' cell for each variable from January through May
#'
#' @param trend_data survey results nested by month
#' @param missing_values option to remove missing values from cell percentages
#'
#' @return dataframe containing cell percentages by demographic category by month
sample_characteristics_all <- function(trend_data, missing_values) {
  
  sample_characteristics_full <- function(data, variable, group_values, missing_values = TRUE) {
    rm_miss <- function(data) {
      data[data[[variable]] != 199]
    }
    if (missing_values) {
      means <- unlist(map(group_values, ~weighted.mean(data[[variable]] == .x, data$weight)))
    }
    if (!missing_values) {
      data <- data %>%
        mutate(ethnicity_race = ifelse(ethnicity_race == 7, 199, ethnicity_race))
      
      means <- unlist(map(group_values, ~weighted.mean(rm_miss(data)[[variable]] == .x, rm_miss(data)$weight)))
    }
    
    tibble(
      variable = variable,
      value = group_values,
      means = means
    )                
  }
  
  iter_trends <- function(data_list, variable, missing_values) {
    map(data_list, ~sample_characteristics_full(.x, variable, unique(.x[[variable]]),
                                                missing_values)) %>%
      map2(c("January", "February", "March", "April", "May"), ~mutate(.x, month = .y)) %>%
      invoke(rbind, .)
  }
  
  map(c("gender", "age_cat", "ethnicity_race", "educ", "emp_out"), 
      ~iter_trends(trend_data$data, .x, missing_values)) %>%
    invoke(rbind, .) %>%
    spread(month, means) %>%
    mutate_at("variable", ~factor(., levels = c("gender", "age_cat", "ethnicity_race", "educ", "emp_out"))) %>%
    arrange(variable) %>%
    mutate_at("value", as.character) %>%
    left_join(codebook, by = c("variable" = "variable_name")) %>%
    left_join(fullvars, by = c("variable" = "variable_name", "value")) %>%
    select(Variable = description, Description = label_full, January, February, March, April, May) %>%
    mutate_if(is.numeric, ~round(100*., 2)) %>%
    select(-Variable)
}



