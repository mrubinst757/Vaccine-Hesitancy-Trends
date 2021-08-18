library(tidyverse)
library(data.table)
library(sandwich)
library(gt)
library(gridExtra)
extrafont::loadfonts(device = "pdf")
extrafont::loadfonts(device = "win")

######################################################################################
############################# create hesitancy tables ################################
######################################################################################

#' process_rdtab: creates formatted table for risk-differences with respect to a binary
#' outcome and a categorical variable
#' @param table_output output from run-analysis cross-tabs function
#' @param crosswalk cross-walk between category values and labels
#' @param codebook crosswalk between variable names and descriptions
#'
#' @return dataframe with risk differences and estimated CIs
process_rdtab <- function(table_output, crosswalk, codebook) {
  table_output %>%
    invoke(rbind, .) %>%
    mutate_at("se", as.numeric) %>%
    left_join(crosswalk, by = c("label" = "new_varname")) %>%
    left_join(codebook, by = "variable_name") %>%
    select(props, se, N = total, reference,
           Category = label_full, Characteristic = description) %>%
    mutate_at(vars("props", "se"), ~100*.) %>%
    select(Characteristic, Category, props, se, N, reference) %>%
    mutate(lci = props - 1.96*se, uci = props + 1.96*se,
           `Proportion (CI)` = format_cells(props, lci, uci, 1, 1)) %>%
    select(-props, -se, -lci, -uci)
}

#' process_rrtab: creates formatted table for risk-ratios with respect to a binary
#' outcome and a categorical variable
#'
#' @param table_output output from run-analysis cross-tabs function
#' @param crosswalk cross-walk between category values and labels
#' @param codebook crosswalk between variable names and descriptions
#'
#' @return dataframe with risk differences and estimated CIs
process_rrtab <- function(tab_output, fullvars, codebook) {
  tab_output %>%
    invoke(rbind, .) %>%
    mutate_at("se", as.numeric) %>%
    left_join(fullvars, by = c("label" = "new_varname")) %>%
    left_join(codebook, by = "variable_name") %>%
    select(rr, se,
           Category = label_full, Characteristic = description) %>%
    select(Characteristic, Category, rr, se) %>%
    mutate(lci = rr - 1.96*se, uci = rr + 1.96*se,
           ci = format_cells(rr, lci, uci, 2, 1)) %>%
    select(-rr, -se, -lci, -uci)
}

# joins the risk-ratio/risk difference results for full sample, unvaccinated sample
#' process_table: calls process_rdtab and process_rrtab on cross-tab output and joins
#' them into one large table
#'
#' @param crosstabs output from run-analysis program generate_xtabs
#' @param crosswalk cross-walk between category values and labels
#' @param codebook crosswalk between variable names and descriptions
#' @param order_cats specified order of categorical variables
#'
#' @return datatable with risk-ratios and risk-differences for binary
#' outcome across multiple categorical variables
process_table <- function(crosstabs, codebook, crosswalk, order_cats) {
  rdv1 <- process_rdtab(crosstabs[[1]], crosswalk, codebook)
  rdv2 <- process_rdtab(crosstabs[[2]], crosswalk, codebook)
  rrtab1 <- process_rrtab(crosstabs[[3]], crosswalk, codebook) 
  rrtab2 <- process_rrtab(crosstabs[[4]], crosswalk, codebook) 
  
  rdtab <- rdv1 %>%
    left_join(rdv2, by = c("Characteristic", "Category", "reference"))
  
  rrtabf <- left_join(rrtab1, rrtab2, by = c("Characteristic", "Category"))
  
  rdtabf <- rdtab %>%
    left_join(rrtabf) %>%
    replace_na(list(ci.x = "1.0 (NA)", ci.y = "1.0 (NA)")) %>%
    group_by(Characteristic) %>%
    mutate(cat_pct1 = format(round(100*N.x/sum(N.x), 1), nsmall = 1), 
           cat_pct2 = format(round(100*N.y/sum(N.y), 1), nsmall = 1)) %>%
    ungroup() %>%
    mutate(reorder = ifelse(grepl(order_cats, Characteristic), 1, 0)) %>%
    nest(-reorder) %>%
    mutate(data = ifelse(reorder == 1, map(data, ~.x %>% group_by(Characteristic) %>% arrange(-N.x)), data)) %>%
    unnest() %>%
    select(-reorder) %>%
    mutate_if(is.character, ~gsub("\\( ", "\\(", .))
}

#' final_process: add regression adjusted risk-ratios to primary cross-tabs 
#'
#' @param reg_results regression results
#' @param reg_results.a regression results for second covariate set
#' @param main_table output from main_table
#' @param codebook codebook linking variable names to descriptions
#'
#' @return dataframe with risk ratios, differences, and adjusted risk ratios
final_process <- function(reg_results, reg_results.a, main_table, codebook) {
  gen_reg_table <- function(results) {
    results %>%
      select(ci.adj, variable_name, label_full, sample) %>%
      spread(sample, ci.adj) %>%
      rename(ci.adj.0 = full, ci.adj.1 = unvaccinated) %>%
      left_join(codebook %>% select(variable_name, description), by = "variable_name") %>%
      rename(Category = label_full, Characteristic = description) %>%
      select(-variable_name) 
  }
  reg.table <- gen_reg_table(reg_results)
  reg.table.a <- gen_reg_table(reg_results.a) %>%
    rename(ci.adj.0.a = ci.adj.0, ci.adj.1.a = ci.adj.1)
  reg.full <- full_join(reg.table, reg.table.a, by = c("Characteristic", "Category")) 
  
  left_join(main_table, reg.full, by = c("Characteristic", "Category")) %>%
    replace_na(list(ci.adj.0 = "-", ci.adj.1 = "-", ci.adj.0.a = "-", ci.adj.1.a = "-")) %>%
    mutate_if(is.character, ~gsub("\\( ", "\\(", .))
}


#' complete_table: format cross-tabs using gt presenting all available results
#'
#' @param table output from process_table
#' @param row_order order for rows 
#' @param outcome binary outcome
#'
#' @return formatted gt table with all results from process_table
complete_gt_table <- function(table, row_order, outcome = "Hesitant ") {
  char.df <- table %>%
    select(Characteristic, Category) %>%
    mutate(id = 1:nrow(.)) %>%
    nest(-Characteristic)  
  
  char.list <- map(char.df$data, ~.x$id)
  group.list <- char.df$Characteristic
  
  if (outcome == "Hesitant ") {
    rr_label <- "RR (95% CI of Hesitancy)"
  }
  if (outcome != "Hesitant ") {
    rr_label <- "RR (95% CI)"
  }
  gt_tbl <- table %>%
    select(-Characteristic) %>%
    gt(rowname_col = "Category") %>%
    tab_spanner(
      label = "Full Data", 
      columns = vars( "N.x", "cat_pct1", "Proportion (CI).x", "ci.x", "ci.adj.0", "ci.adj.0.a")
    ) %>%
    tab_spanner(
      label = "Unvaccinated", 
      columns = vars("N.y", "cat_pct2", "Proportion (CI).y", "ci.y", "ci.adj.1", "ci.adj.1.a")
    ) %>%
    cols_label(`Proportion (CI).x` = paste0(outcome, "(95% CI)"),
               ci.adj.0 = "Adjusted RR (95% CI)",
               ci.adj.1 = "Adjusted RR (95% CI)",
               ci.adj.0.a = "Adjusted RR (Full health cats)",
               ci.adj.1.a = "Adjusted RR (Full health cats)",
               N.x = "N", cat_pct1 = "% Sample",
               ci.x = rr_label, ci.y = rr_label,
               `Proportion (CI).y` = paste0(outcome, "(95% CI)"),
               N.y = "N", cat_pct2 = "% Sample")
  
  gt_tbl <- gt_tbl %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(rows = reference == 1)
    ) %>%
    cols_hide(columns = "reference") %>%
    tab_style(style = list(cell_fill(color = "gray")), 
              locations = cells_body(columns = vars(N.y, `Proportion (CI).y`, ci.y, cat_pct2, ci.adj.1,
                                                    ci.adj.1.a)))
  for(i in 1:length(char.list)) {
    gt_tbl <- add_gt_rows(gt_tbl, group.list[i], char.list[[i]])
  }
  gt_tbl %>%
    row_group_order(groups = row_order)
}

#' final_table: produce tables for paper (excluding sensitivity results)
#'
#' @param table output from process_table
#' @param row_order order for rows
#' @param outcome binary outcome
#' @param hesitant.only specify whether output is among only the hesitant or everyone
#'
#' @return formatting gt table
final_gt_table <- function(table, row_order, outcome = "Hesitant ", hesitant.only = FALSE) {
  char.df <- table %>%
    select(Characteristic, Category) %>%
    filter(Characteristic %in% row_order) %>%
    mutate(id = 1:nrow(.)) %>%
    nest(-Characteristic) 
  
  char.list <- map(char.df$data, ~.x$id)
  group.list <- char.df$Characteristic
  
  if (outcome == "Hesitant ") {
    rr_label <- "RR (95% CI of Hesitancy)"
  }
  if (outcome != "Hesitant ") {
    rr_label <- "RR (95% CI)"
  }
  
  if (hesitant.only == FALSE) {
    table <- filter(table, Characteristic %in% row_order) %>%
      select(Category, reference, N.x, cat_pct1, `Proportion (CI).x`, ci.x, ci.adj.0) 
  }
  
  if (hesitant.only == TRUE) {
    table <- filter(table, Characteristic %in% row_order) %>%
      select(Category, reference, N.x = N.y, 
             cat_pct1 = cat_pct2, 
             `Proportion (CI).x` = `Proportion (CI).y`, 
             ci.x = ci.y, ci.adj.0 = ci.adj.1) 
  }
  gt_tbl <- table %>%
    gt(rowname_col = "Category") %>%
    cols_label(`Proportion (CI).x` = paste0(outcome, "(95% CI)"),
               ci.adj.0 = "Adjusted RR (95% CI)",
               N.x = "N", cat_pct1 = "% Sample",
               ci.x = rr_label)
  
  gt_tbl <- gt_tbl %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(rows = reference == 1)
    ) 
  if (length(group.list) > 1) {
    for(i in 1:length(char.list)) {
      gt_tbl <- add_gt_rows(gt_tbl, group.list[i], char.list[[i]])
    }
    res <- gt_tbl %>%
      cols_hide(columns = "reference") %>%
      row_group_order(groups = row_order) 
  }
  if (length(group.list == 1)) {
    res <- gt_tbl %>%
      cols_hide(columns = "reference") 
  }
  return(res)
}

#' add_gt_rows: inserts rows within a tab group given row ids
#'
#' @param gttable gt table
#' @param group categorical variable
#' @param row_ids row numbers where the table has the specified category from the 
#' categorical variable
#'
#' @return gt table with tabs for categorical variables
add_gt_rows <- function(gttable, group, row_ids) {
  final <- gttable %>%
    tab_row_group(group = group, rows = row_ids)
  final
} 

###################################################################################
########################## age-race linear combination tables #####################
###################################################################################

#' lincom_table_race: generate linear combinations of risk-ratios within
#' each age group by race
#'
#' @param rr_by_age_race linear combinations by age within race
#' @param arr_by_age adjusted linear combinations
#' @param full_results output from final_process
#'
#' @return table with specified linear combinations
lincom_table_race <- function(rr_by_age_race, arr_by_age, full_results) {
  race_unadj <- rr_by_age_race %>%
    filter(group == "By Race") %>%
    rename(Age = year_cat, Race = comp_cat) %>%
    mutate(`RR (95% CI of Hesitancy)` = paste0(format(round(rr, 2), nsmall = 2), " (", 
                                               format(round(lci, 2), nsmall = 2), ", ", 
                                               format(round(uci, 2), nsmall = 2), ")")) %>%
    select(Age, Race, `RR (95% CI of Hesitancy)`) %>%
    mutate(`RR (95% CI of Hesitancy)` = ifelse(`RR (95% CI of Hesitancy)` == "1.00 (1.00, 1.00)", "1.00 (NA)", 
                                               `RR (95% CI of Hesitancy)`)) %>%
    mutate_at("RR (95% CI of Hesitancy)", ~gsub("\\( ", "\\(", .))
  
  race_adj <- arr_by_age %>%
    rename(Age = year_cat, Race = comp_cat) %>%
    mutate(`Adjusted RR (95% CI)` = paste0(format(round(rr, 2), nsmall = 2), " (", 
                                           format(round(lci, 2), nsmall = 2), ", ", 
                                           format(round(uci, 2), nsmall = 2), ")")) %>%
    mutate(`Adjusted RR (95% CI)` = ifelse(`Adjusted RR (95% CI)` == "1.00 (1.00, 1.00)", "1.00 (NA)", 
                                           `Adjusted RR (95% CI)`)) %>% 
    select(Age, Race, `Adjusted RR (95% CI)`) %>%
    mutate_at("Adjusted RR (95% CI)", ~gsub("\\( ", "\\(", .))
  
  table_one <- full_results %>%
    separate(Category, c("Race", "Age"), remove = FALSE, sep = "\\:") %>%
    mutate_if(is.character, trimws) %>%
    left_join(race_unadj, by = c("Age", "Race")) %>%
    left_join(race_adj, by = c("Age", "Race")) %>%
    arrange(Age) %>%
    select(-RR)
  
  groups <- rev(unique(table_one$Age)[1:7])
  
  rows <- rev(map(0:(length(groups)-1), ~.x*7 + c(1:7)))
  
  table_one_gt <- table_one %>%
    filter(!grepl("Unknown", Category)) %>% 
    select(-Age, -Characteristic, -Category) %>%
    gt(rowname_col = "Race") 
  
  for(i in 1:length(groups)) {
    table_one_gt <- add_gt_rows(table_one_gt, groups[[i]], rows[[i]]) 
  }
  table_one_gt
}

#' lincom_table_age: generate linear combinations of risk-ratios within
#' each race group by age
#'
#' @param rr_by_age_race linear combinations by age within race
#' @param arr_by_race adjusted linear combinations
#' @param full_results output from final_process
#'
#' @return table with specified linear combinations
lincom_table_age <- function(rr_by_age_race, arr_by_race, full_results) {
  age_unadj <- rr_by_age_race %>%
    filter(group == "By Age") %>%
    rename(Age = comp_cat, Race = race) %>%
    mutate(`RR (95% CI of Hesitancy)` = paste0(format(round(rr, 2), nsmall = 2), " (", 
                                               format(round(lci, 2), nsmall = 2), ", ", 
                                               format(round(uci, 2), nsmall = 2), ")")) %>%
    select(Age, Race, `RR (95% CI of Hesitancy)`) %>%
    mutate(`RR (95% CI of Hesitancy)` = ifelse(`RR (95% CI of Hesitancy)` == "1.00 (1.00, 1.00)", "1.00 (NA)", 
                                               `RR (95% CI of Hesitancy)`)) %>%
    mutate_at("RR (95% CI of Hesitancy)", ~gsub("\\( ", "\\(", .))
  
  
  age_adj <- arr_by_race %>%
    rename(Age = comp_cat, Race = race) %>%
    mutate(`Adjusted RR (95% CI)` = paste0(format(round(rr, 2), nsmall = 2), " (", 
                                           format(round(lci, 2), nsmall = 2), ", ", 
                                           format(round(uci, 2), nsmall = 2), ")")) %>%
    mutate(`Adjusted RR (95% CI)` = ifelse(`Adjusted RR (95% CI)` == "1.00 (1.00, 1.00)", "1.00 (NA)", 
                                           `Adjusted RR (95% CI)`)) %>% 
    select(Age, Race, `Adjusted RR (95% CI)`) %>%
    mutate_at("Adjusted RR (95% CI)", ~gsub("\\( ", "\\(", .))
  
  table_two <- full_results %>%
    separate(Category, c("Race", "Age"), remove = FALSE, sep = "\\:") %>%
    mutate_if(is.character, trimws) %>%
    left_join(age_unadj, by = c("Age", "Race")) %>%
    left_join(age_adj, by = c("Age", "Race")) %>%
    select(-RR)
  
  groups <- rev(unique(table_two$Race)[1:7])
  rows <- rev(map(0:(length(groups)-1), ~.x*7 + c(1:7)))
  
  table_two_gt <- table_two %>%
    filter(!grepl("Unknown", Category)) %>% 
    select(-Race, -Characteristic, -Category) %>%
    gt(rowname_col = "Age") 
  
  for(i in 1:length(groups)) {
    table_two_gt <- add_gt_rows(table_two_gt, groups[[i]], rows[[i]]) 
  }
  table_two_gt
}

#####################################################################################
################################ create plots #######################################
#####################################################################################

#' create_figure_one: generate figure depicting time trends in hesitancy from
#' january through may
#'
#' @param trend_data data containing the time trend information
#'
#' @return ggplot output
create_figure_one <- function(trend_data) {
  trend_data %>%
    filter(grepl("vaccine_intent", variable), month %in% c(1:5)) %>%
    mutate_at("label", ~factor(., levels = rev(c("Already vaccinated", "Yes, definitely", "Yes, probably",
                                                 "No, probably not", "No, definitely not")))) %>%
    group_by(month) %>%
    arrange(label) %>%
    mutate_at("month", ~month.name[.]) %>%
    mutate_at("month", ~factor(., levels = c("January", "February", "March", "April", "May"))) %>%
    mutate(lci = cumsum(hesitant_lci), uci = cumsum(hesitant_uci)) %>%
    mutate_at("label", ~factor(., levels = (c("Already vaccinated", "Yes, definitely", "Yes, probably",
                                              "No, probably not", "No, definitely not")))) %>%
    ggplot(aes(x = month, y = hesitant_props, fill = label, 
               ymin = lci, ymax = uci)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_errorbar(stat = "identity", width = 0.2) +
    scale_fill_manual(values = rev(option1)) +
    ylab("Vaccine hesitant (weighted %)") +
    xlab("Month") +
    theme_minimal() +
    labs(fill = "Would you take \n the COVID-19 \n vaccine if offered \n today?") +
    theme(text = element_text(size = 14, family = "Arial", color = "Black"),
          axis.text.x = element_text(color = "Black", size = 12),
          axis.text.y = element_text(color = "Black", size = 12))
}

#' create_figure_two: generate figure depicting time trends in hesitancy from
#' january through may for specified subgroups
#'
#' @param plot.data data containing the time trend information
#' @param category categorical variable being plottes
#' @param palette rcolorbrewer palette
#' @param ymax maximum value on y-axis
#'
#' @return ggplot output
create_figure_two <- function(plot.data, category, palette, ymax) {
  my_levels <- fullvars %>%
    filter(variable_name == category) %>%
    .$label_full
  
  shapes <- 5 + length(my_levels) -1
  
  plot.data %>%
    filter(variable == category) %>%
    mutate_at("label_full", ~factor(., levels = my_levels)) %>%
    ggplot(aes(x = ggmonth, y = 100*hesitant_props, color = label_full, 
               group = label_full, shape = label_full)) +
    geom_point() +
    geom_line() +
    scale_color_brewer(palette = palette) +
    facet_wrap(~variable_name_full) +
    ylab("Vaccine hesitant (%)") +
    xlab("Month") +
    scale_x_date(breaks = "1 month", date_labels = "%B") +
    theme_minimal() +
    scale_shape_manual(values=c(5:shapes)) +
    scale_y_continuous(minor_breaks = seq(0, 65, 5), 
                       breaks = seq(0, 65, 10), 
                       limits = c(0, ymax)) +
    guides(color = guide_legend(nrow = 2)) +
    theme(legend.position = "bottom", legend.title = element_blank())
}

#' format_plot: add a letter in the left hand side of a plot for panels
#'
#' @param plot a ggplotobject (from create_figure_two)
#' @param letter letter to put in plot
#'
#' @return formatted ggplot object
format_plot <- function(plot, letter) {
  library(grid)
  arrangeGrob(plot, top = textGrob(letter, x = unit(0, "npc"), 
                                   y = unit(1, "npc"), just=c("left","top"),
                                   gp=gpar(col="black", fontsize=10, fontfamily="Times Roman")))
}

#' create_figure_three: comparison of hesitancy by age and race
#' @param plot_data data containing the time trend information
#'
#' @return ggplot output
create_figure_three <- function(plot_data) {
  race_factor <- c("White", "Hispanic", "Black", "Asian", 
                   "Native American", "Pacific Islander", "Multi-racial",
                   "Unknown \n (other or no response)", "No response (either category)")
  
  age_factor <- c("18-34 years", "35-64 years", "65 or more years", "No response")
  
  plot_data %>%
    filter(Age != " unknown Age") %>%
    mutate_at("Age", trimws) %>% 
    mutate_at("Race", ~gsub("Unknown \\(other or no response\\)", 
                            "Unknown \n \\(other or no response\\)", .)) %>%
    mutate(Race = factor(Race, levels = race_factor)) %>%
    filter(!Race %in% c("Unknown \n (other or no response)", "Multi-racial")) %>%
    ggplot(aes(y = Age, x = Hesitancy, color = Race)) +
    facet_wrap(~Race, ncol = 3, nrow = 2) +
    geom_point() +
    geom_errorbar(aes(xmin = l95ci, xmax = u95ci)) +
    theme_minimal() +
    xlab("Vaccine hesitant (%) \n (95% CI errorbars)") +
    ylab("Age") +
    labs(title = "Race/ethnicity") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none") +
    scale_color_brewer(palette = "Dark2") +
    theme(text = element_text(size = 14, family = "Arial", color = "Black"),
          axis.text.x = element_text(color = "Black", size = 12),
          axis.text.y = element_text(color = "Black", size = 12),
          plot.title = element_text(size = 14),
          strip.text = element_text(color = "Black", size = 12))
}

