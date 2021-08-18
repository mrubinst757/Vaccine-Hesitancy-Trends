source("02_AnalyzeData/trends-paper/analysis-funs.R")
source("02_AnalyzeData/trends-paper/table-funs.R")

# specify formatting options (table ordering) ---------------------------------------------------
codebook <- read_csv("../00_RawData/codebook/codebook-interactions.csv")

fullvars <- read_csv("../00_RawData/codebook/discrete-vars-interactions.csv") %>%
  mutate(label_full = ifelse(label_full == "No response (either category)",
                              "Unknown (other or no response): unknown Age", label_full))

final_table1 <- c("gender", "age_cat", "ethnicity_race", "educ", "emp_out",
                  "region", "urban_fips")

final_table2 <- c("deaths_prop_cat",  "rep_lead_ntile", "governor",
                   "vaccine_flu", "covid_pos_ever", "healthc",
                   "elderly_risk", "elderly_risk_det", "behav_worry", "behav_contact")

final_table3 <- "health"

full_groups <- c(final_table1, final_table2, final_table3, "age_race", "age_race2")

main_order <- codebook %>%
  filter(variable_name %in% full_groups) %>%
  mutate(variable_name = factor(variable_name, levels = full_groups)) %>%
  arrange(variable_name) %>%
  .$description

final_order1 <- main_order[1:length(final_table1)]
final_order2 <- main_order[(length(final_table1) + 1):(length(final_table1) + length(final_table2))]
final_order3 <- main_order[length(main_order) - 2]

out_dir <- "C:/Users/mdrub/Box/Vaccine-Intent-Project/trends-paper/"

#########################################################################################
############################### hesitancy tables ########################################
#########################################################################################

main_xtabs <- readRDS("02_Output/main-tables.rds")
sens_xtabs <- readRDS("02_Output/sens-tables.rds")
main_model_results <- read_csv("02_Output/model-results.csv") 
main_model_results.a <- read_csv("02_Output/model-results-a.csv") 
sens_model_results <- read_csv("02_Output/model-results-sens.csv") 
sens_model_results.a <- read_csv("02_Output/model-results-sens-a.csv") 
full_table_main <- process_table(main_xtabs, codebook, fullvars, "Race and Ethnicity|diagnosed") %>%
  final_process(main_model_results, main_model_results.a, ., codebook)
sens_table_main <- process_table(sens_xtabs, codebook, fullvars, "Race and Ethnicity|diagnosed") %>%
  final_process(sens_model_results, sens_model_results.a, ., codebook)

write_rds(full_table_main, "02_Output/full-results.rds")

health_table <- process_table(main_xtabs, codebook, fullvars, "Race and Ethnicity|diagnosed") %>%
  final_process(main_model_results, main_model_results.a, ., codebook) %>%
  select(-ci.adj.0, -ci.adj.1) %>%
  rename(ci.adj.0 = ci.adj.0.a, ci.adj.1 = ci.adj.1.a) %>%
  filter(Characteristic == "Ever diagnosed with \"high-risk\" condition") %>%
  separate(`Proportion (CI).x`, into = c("rd", "ci"), remove = FALSE, sep = "\\(") %>%
  mutate(rd = as.numeric(rd)) %>%
  arrange(rd)

# complete results
gtsave(complete_gt_table(full_table_main, main_order), 
       paste0(out_dir, "complete-results/risk-table.html"))
gtsave(complete_gt_table(full_table_sens, main_order, outcome = "Probably "), 
       paste0(out_dir, "complete-results/risk-table-sens.html"))

# paper tables
gtsave(final_gt_table(full_table_main, final_order1), paste0(out_dir, "tables/table-1.html"))
gtsave(final_gt_table(full_table_main, final_order2), paste0(out_dir, "tables/table-2.html"))
gtsave(final_gt_table(health_table, final_order3), paste0(out_dir, "tables/supp-table-5.html"))
gtsave(final_gt_table(full_table_sens, final_order1, outcome = "Probably "), paste0(out_dir, "tables/supp-table-6.html"))
gtsave(final_gt_table(full_table_sens, final_order2, outcome = "Probably "), paste0(out_dir, "tables/supp-table-7.html"))

#########################################################################################
############################## age/race interactions ####################################
#########################################################################################

# tables ------------------------------------------------------------------------
out_dir <- "C:/Users/mdrub/Box/Vaccine-Intent-Project/trends-paper/"

# read data 
full_results <- readRDS("02_Output/full-results.rds") %>%
  filter(grepl("by Race", Characteristic), !grepl("collapsed", Characteristic)) %>%
  select(Characteristic, Category, N = N.x, `% Sample` = cat_pct1,
         `Hesitant (95% CI)` = `Proportion (CI).x`,
         RR = ci.x)

rr_by_age_race <- read_csv("02_Output/age-race-lincoms-unadjusted.csv")
arr_by_race <- read_csv("02_Output/arr-by-race.csv")
arr_by_age <- read_csv("02_Output/rr-by-age-race.csv") %>%
  filter(ref_cat == "White")

# create race tables stratified by age group 
table_one <- lincom_table_race(rr_by_age_race, arr_by_age, full_results)
gtsave(table_one, paste0(out_dir, "tables/supp-table-3.html"))

# create age tables stratified by race 
table_two <- lincom_table_age(rr_by_age_race, arr_by_race, full_results)
gtsave(table_two, paste0(out_dir, "tables/supp-table-4.html"))

# create figure 3: hesitancy by age and racial group
figure_three_data <- full_table_main %>% 
  filter(grepl("Race/Ethnicity", Characteristic)) %>%
  select(Characteristic, Category, N.x, Hesitancy = `Proportion (CI).x`) %>%
  separate(Hesitancy, c("Hesitancy", "CI"), sep = "\\(") %>%
  separate(CI, c("l95ci", "u95ci"), ",") %>%
  separate(Category, into = c("Race", "Age"), ":", remove = FALSE) %>%
  mutate_at(vars(contains("ci")), ~as.numeric(gsub("\\)", "", .))) %>%
  mutate(Hesitancy = as.numeric(Hesitancy)) %>%
  filter(!grepl("collapsed", Characteristic))

figure3 <- create_figure_three(plot.dat.full)

ggsave(figure3,
       filename = paste0(out_dir, "figures/figure-3.tiff"), 
       height = 5, width = 7, dpi = 300)

######################################################################################################
################################### monthly hesitancy trends #########################################
######################################################################################################
trend_data <- readRDS("02_Output/trend-paper-hesitancy-trends.rds") %>%
  mutate(ggmonth = lubridate::ydm(paste0("2021010", month)))

codebook <- read_csv("../00_RawData/codebook/codebook.csv") %>%
  mutate_at("variable_name_full", ~gsub("Race and ethnicity",
                                        "Race and ethnicity (18 - 34 year olds)", .)) %>%
  mutate_at("variable_name_full", ~gsub("County Trump vote total minus Biden vote total among all votes",
                                        "Trump vote share", .))

fullvars <- read_csv("../00_RawData/codebook/discrete-vars.csv") %>%
  mutate_at("label_full", ~gsub("High school, GED, or less", "High school or less", .)) %>%
  mutate_at("label_full", ~gsub("Some college or two year degree", "Some college", .)) %>%
  mutate_at("label_full", ~gsub("Professional degree", "Professional \n degree", .))

supp_table1 <- trend_data %>%
  filter(grepl("vaccine_intent|all", variable), month != 7) %>%
  mutate_at("label", ~gsub("all", "Hesitant (probably/definitely not)", .)) %>%
  mutate(props = format_cells(hesitant_props, hesitant_lci, hesitant_uci, 1, 100)) %>%
  select(props, month, label) %>%
  mutate_all(~gsub("\\( ", "\\(", .)) %>%
  mutate_all(~gsub("\\,  ", "\\, ", .)) %>%
  spread(month, props) %>%
  set_names(c("Response", "January", "February", "March", "April", "May",
              "Difference (May - January)")) %>%
  mutate_at("Response", ~factor(., levels = c("Already vaccinated", "Yes, definitely", "Yes, probably",
                                           "No, probably not", "No, definitely not", 
                                           "Hesitant (probably/definitely not)"))) %>%
  arrange(Response) %>%
  gt() 

gtsave(supp_table1, paste0(out_dir, "tables/supp-table-1.html"))

# specify color palette for plots and intent values
option1 <- RColorBrewer::brewer.pal(8, "Paired")[c(6,5,3,4,7)]

# figure 1: hesitancy and vaccination status over time
figure1 <- create_figure_one(trend_data)
ggsave(figure1, height = 5, width = 7, dpi = 300,
       filename = paste0(out_dir, "figures/figure-1.tiff"))

# create table version of figure + the trend differences calculated previously
table_orders <- fullvars %>%
  filter(variable_name %in% c("educ", "ethnicity_race", "region", "rep_lead_ntile")) %>%
  .$label_full

supp_table2 <- trend_data %>%
  filter(!variable %in% c("vaccine_intent", "all"), !grepl("199|NA|^region$|region_6|ethnicity_race_7", label),
         month %in% c(1:6)) %>%
  left_join(fullvars %>% select(label_full, new_varname), by = c("label" = "new_varname")) %>%
  left_join(codebook, by = c("variable" = "variable_name")) %>%
  mutate(props = format_cells(hesitant_props, hesitant_lci, hesitant_uci, 1, 100)) %>%
  select(month, Description = label_full, variable_name_full, props) %>%
  spread(month, props) %>%
  set_names(c("Description", "Variable", "January", "February", "March", 
              "April", "May", "Difference (May - January)")) %>%
  mutate_at("Description", ~factor(., levels = table_orders)) %>%
  arrange(Variable, Description) %>%
  mutate_all(~gsub("\\( ", "\\(", .)) %>%
  mutate_all(~gsub("\\,  ", "\\, ", .)) %>%
  mutate_at("Variable", ~factor(., levels = c("Race and ethnicity (18 - 34 year olds)",
                                "Trump vote share", "US Region", "Education level"))) %>%
  arrange(Variable)

group_rows <- map(unique(supp_table2$Variable), ~grep(.x, supp_table2$Variable))
group_rows[[1]] <- c(1:7)

supp_table2_gt <- supp_table2 %>%
  select(-Variable) %>%
  filter(!is.na(label)) %>%
  gt(rowname_col = "Description") 

for(i in 1:length(group_rows)) {
  label <- paste0(unique(supp_table2$Variable)[i], "; percent hesitant (95% CI)")
  supp_table2_gt <- add_gt_rows(supp_table2_gt, label, group_rows[[i]]) 
}

gtsave(supp_table2_gt, paste0(out_dir, "tables/supp-table-2.html"))

# create figure 2: trends by category
plot_data <- trend_data %>%
  filter(!variable %in% c("vaccine_intent", "all"), !grepl("199|NA|^region$|region_6|ethnicity_race_7", label),
         month %in% c(1:5)) %>%
  mutate_at("month", ~month.name[.]) %>%
  mutate_at("month", ~factor(., levels = c("January", "February", "March", 
                                           "April", "May"))) %>%
  left_join(fullvars %>% select(label_full, new_varname), by = c("label" = "new_varname")) %>%
  left_join(codebook, by = c("variable" = "variable_name"))

plot1 <- create_figure_two(plot_data, "ethnicity_race", "Set1", 60)
plot2 <- create_figure_two(plot_data, "rep_lead_ntile", "Set1", 45)
plot3 <- create_figure_two(plot_data, "region", "Set1", 45)
plot4 <- create_figure_two(plot_data, "educ", "Set1", 45)

saveRDS(list(plot1, plot4, plot3, plot2), "02_Output/trends-paper-trend-plots.RDS")
main_plots <- readRDS("02_Output/trends-paper-trend-plots.RDS")

final_plots <- map2(main_plots, c("A", "B", "C", "D"), format_plot)
tiff(paste0(out_dir, "figures/figure-2.tiff"), 
     res = 260, width = 4.8*480, height = 4.6*480*(2/3))
gridExtra::grid.arrange(final_plots[[1]], final_plots[[2]], 
                        final_plots[[3]], final_plots[[4]], 
                        ncol = 2, nrow = 2,
                        heights = c(3, 3), widths = c(4.5, 4.5))
dev.off() # Close the file

######################################################################################################
################################## hesitancy reasons tables ##########################################
######################################################################################################
may_data <- readRDS("../00_RawData/processed-data/may-data-08-16-2021.rds") %>%
  filter(!is.na(intent_binary), gender != 4, !is.na(reasons)) %>%
  select(vaccine_intent, weight, ethnicity_race, starts_with("reason_"))

reason_values <- read_csv("../00_RawData/codebook/hesitancy-reasons.csv") %>%
  mutate(varname = paste0("reason_", value))

# stratify data by reason
defno <- filter(may_data, vaccine_intent == 4) 
probno <- filter(may_data, vaccine_intent == 3)
probyes <- filter(may_data, vaccine_intent == 2)
pyes1 <- filter(may_data, vaccine_intent == 2)

# combine stratified data for marginal totals
marg <- rbind(defno, probno, probyes)
hesitant <- rbind(defno, probno)

# create reasons table for each vaccine intention category
reasons_by_intent <- map2(list(hesitant, defno, probno, probyes), c("Total among hesitant",
                                                         "Definitely not choose to be vaccinated today",
                                                         "Probably not choose to be vaccinated today",
                                                         "Probably choose to be vaccinated today"),
               reasons_table) %>% 
  reduce(left_join, by = "Reason")

order_var <- names(reasons_by_intent)[2]

table_3 <- reasons_by_intent %>% 
  separate(order_var, c("value", "other"), "\\(", remove = FALSE) %>%
  arrange(-as.numeric(value)) %>%
  mutate_all(~gsub("\\( ", "\\(", .)) %>%
  mutate_all(~gsub("\\,  ", "\\, ", .))

other <- filter(table_3, grepl("^Other$", Reason))

table_3 <- table_3 %>%
  filter(!grepl("^Other$", Reason)) %>%
  bind_rows(other)

table_3 %>%
  select(-value, -other) %>%
  gt(rowname_col = "Reasons") %>%
  gtsave(paste0(out_dir, "tables/table-3.html"))

# create reasons table for each race/ethnicity group
hesitant <- split(hesitant, f = hesitant$ethnicity_race)
reasons_by_race <- map2(hesitant, c("Hispanic", "Native American", "Asian",
                                "Black", "Pacific Islander", "White",
                                "Unknown (other or no response)",
                                "Multiracial"), reasons_table) %>% 
  reduce(left_join, by = "Reason")

supp_table_8 <- reasons_by_race %>%
  select(Reason, contains("White"), contains("Hispanic"),
         contains("Black"), contains("Asian"), contains("Native American"),
         contains("Pacific Islander"), contains("Multiracial"),
         contains("Unknown"))

supp_table_8 %>% 
  mutate_all(~gsub("\\( ", "\\(", .)) %>%
  mutate_all(trimws) %>%
  mutate_at("Reason", ~factor(., levels = table_3$Reason)) %>%
  arrange(Reason) %>%
  gt(rowname_col = "Reason") %>%
  gtsave(paste0(out_dir, "tables/supp-table-8.html"))

########################################################################################
############################## sample characteristics ##################################
########################################################################################
trend_data <- readRDS("../00_RawData/processed-data/jan-may-trendsdata.rds") 

may_data <- readRDS("../00_RawData/processed-data/may-data-08-16-2021.rds")

possible_trolls <- readRDS("../00_RawData/processed-data/may-data-08-16-2021.rds") %>%
  filter(!is.na(intent_binary), gender == 4) 

characteristics_table <- sample_characteristics(filter(may_data, !is.na(hesitant), gender != 4))
troll_characteristics <- sample_characteristics(possible_trolls)
flow_table <- sample_flow(trend_data, may_data)

gtsave(gt(flow_table), paste0(out_dir, "tables/sample-flow.html"))
gtsave(gt(characteristics_table), paste0(out_dir, "tables/sample-characteristics.html"))
gtsave(gt(troll_characteristics), paste0(out_dir, "tables/troll-characteristics.html"))
