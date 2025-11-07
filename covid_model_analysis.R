
# SET WORKING DIRECTORY TO THE FOLDER CONTAINING MODEL CODE AND DATA
# CREATE A FOLDER IN THE WORKING DIRECTORY NAMED "Results"



library(deSolve)
library(ggplot2)
library(socialmixr)
library(readxl)
library(BCEA)
library(janitor)
library(openxlsx)
library(dplyr)
library(tibble)
library(tidytable)
library(fixr)
library(scales)
library(gridExtra)
library(matlib)
library(forcats)
library(geomtextpath)
library(cowplot)
library(stringr)
library(ggbump)
library(ggtext)
library(grid)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(RCurl)

# Set seed and GGplot theme

set.seed(1)
theme_set(theme_light(base_size = 15))
options(scipen = 999)

source("covid_model_run.R")
source("covid_model_odes.R")
source("calculate_variant_param_change.R")
source("covid_model_economic_parameters.R")

# setting model parameters

n_samples_fit <- 500

fitted_param_sample <- read_excel("Data/Model fit results.xlsx", sheet = "param")
fitted_param_sample <- t(fitted_param_sample)
fitted_param_sample <- cbind(rep(0, 13), fitted_param_sample)

beta_fit <- as.numeric(fitted_param_sample[1, -1])
reinf_fit <- as.numeric(fitted_param_sample[2, -1])
tau_fit <- as.numeric(fitted_param_sample[3, -1])
severity_scaling_0_5_fit <- as.numeric(fitted_param_sample[4, -1])
severity_scaling_6_17_fit <- as.numeric(fitted_param_sample[5, -1])
severity_scaling_18_64_fit <- as.numeric(fitted_param_sample[6, -1])
severity_scaling_65_84_fit <- as.numeric(fitted_param_sample[7, -1])
severity_scaling_85_fit <- as.numeric(fitted_param_sample[8, -1])
hosp_mortality_scaling_0_5_fit <- as.numeric(fitted_param_sample[9, -1])
hosp_mortality_scaling_6_17_fit <- as.numeric(fitted_param_sample[10, -1])
hosp_mortality_scaling_18_64_fit <- as.numeric(fitted_param_sample[11, -1])
hosp_mortality_scaling_65_84_fit <- as.numeric(fitted_param_sample[12, -1])
waned_ve_against_infection_fit <- as.numeric(fitted_param_sample[13, -1])

par_fit <- rbind(beta_fit, reinf_fit, tau_fit, severity_scaling_0_5_fit, severity_scaling_6_17_fit, severity_scaling_18_64_fit,
                 severity_scaling_65_84_fit, severity_scaling_85_fit, hosp_mortality_scaling_0_5_fit, hosp_mortality_scaling_6_17_fit,
                 hosp_mortality_scaling_18_64_fit, hosp_mortality_scaling_65_84_fit, waned_ve_against_infection_fit)
rownames(par_fit) <- c("Beta", "Protection against reinfection", "Relative transmission from asymptomatic infections",
                       "Infection severity age scaling 0-5", "Infection severity age scaling 6-17", "Infection severity age scaling 18-64",
                       "Infection severity age scaling 65-84", "Infection severity age scaling 85+", "Hospital mortality age scaling 0-5",
                       "Hospital mortality age scaling 6-17", "Hospital mortality age scaling 18-64", "Hospital mortality age scaling 65-84",
                       "Waned VE against infection scaling")
par_fit <- par_fit[, 1:n_samples_fit]
colnames(par_fit) <- seq(1:n_samples_fit)

# Setting variant and vaccination scenarios and the time horizon and setting up empty lists for the results

vacc_scenarios <- c("Comparator", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "B1", "B2", "B3",
                    "B4", "B5", "B6", "B7", "B8", "B9", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9")
df_yearly_inf_hosp_deaths_lc <- df_hosp_deaths_weekly_summaries <- vector(mode='list', length = length(vacc_scenarios))
variant_scenario <- "No variants"
n_yearly_variants <- 1
variant_every_n_years <- 1
time_horizon <- 100

# Running the model

for (i in 1:length(vacc_scenarios)) {
  
  print(vacc_scenarios[i])
  print(Sys.time())
  results <- covid_model_run_fixed_params_variant_scenarios(par_fit = par_fit, variant_scenario = variant_scenario, n_yearly_variants = n_yearly_variants,
                                                            vaccination_scenario = vacc_scenarios[i], time_horizon = time_horizon)
  df_yearly_inf_hosp_deaths_lc[[i]] <- results[[1]]
  df_hosp_deaths_weekly_summaries[[i]] <- results[[2]]
  
}

save.image("Results/Model run 1.RData")

# Setting discount rate

set_costs_discount_rate <- 0.035
set_qalys_discount_rate <- 0.035

future_years <- 0:(time_horizon - 1)
costs_discount_rate <- 1 / ((1 + set_costs_discount_rate) ^ future_years)
qalys_discount_rate <- 1 / ((1 + set_qalys_discount_rate) ^ future_years)

# Vaccine uptake rate

vacc_uptake_rate <- c(rep(0, 5), rep(6.4, 7), rep(31.7, 4), rep(31.7, 2), rep(31.7, 2), rep(40.3, 5), rep(40.9, 5), 
                      rep(44.0, 5), rep(49.5, 5), rep(57.0, 5), rep(64.8, 5), rep(64.8, 5), rep(64.8, 5), 
                      rep(64.8, 5), rep(64.8, 5), rep(70.1, 5), rep(75.5, 5), rep(75.7, 21)) / 100

# Extracting life years remaining

life_years_remaining <- read_excel("Data/Life expectancy (England).xlsx", sheet = "Cohort males and females")

life_years_remaining <- cbind(life_years_remaining[, -1],
                              rep(life_years_remaining[, ncol(life_years_remaining)], 51))
colnames(life_years_remaining) <- 2024:2123
rownames(life_years_remaining) <- 0:100
life_years_remaining <- t(life_years_remaining)

# Calculating costs and QALYs

total_vacc_cost <- total_inf_cost <- total_hosp_cost <- total_lc_cost <- total_vacc_cost_no_discounting <- 
  total_inf_cost_no_discounting <- total_hosp_cost_no_discounting <- total_lc_cost_no_discounting <- total_lc_person_years <- 
  total_qalys <- total_life_years_lost <- total_infections <- total_hosp <- total_icu <- total_deaths <-
  total_vacc <- total_symp_inf <- mean_yearly_deaths <- mean_yearly_hosp <- cases_averted <- hosp_averted <- deaths_averted <-
  life_years_saved <- long_covid_py_averted <- matrix(NA, ncol = length(vacc_scenarios), nrow = ncol(par_fit))

mean_cases_age_strat <- mean_qalys_by_age <- mean_vacc_yearly_coverage_by_age <- 
  mean_vacc_age_strat <- matrix(NA, ncol = length(vacc_scenarios), nrow = 101)
mean_annual_cases <- matrix(NA, ncol = length(vacc_scenarios), nrow = 100)

for (j in 1:length(vacc_scenarios)) { 
  
  total_cases_age_strat <- total_qalys_by_age <- vacc_yearly_coverage_by_age <- total_vacc_age_strat <-
    matrix(NA, ncol = 101, nrow = ncol(par_fit))
  annual_cases <-  matrix(NA, ncol = 100, nrow = ncol(par_fit))
  
  for (i in 1:ncol(par_fit)) {
    
    # Set dataframe for vaccination scenario
    
    df_results <- df_yearly_inf_hosp_deaths_lc[[j]]
    
    # Calculating infections, symptomatic infections, hospitalisations, ICU admissions, deaths and vaccinations

    total_infections[i, j] <- sum(df_results$inc_inf_yearly[[i]])
    total_hosp[i, j] <- sum(df_results$inc_hosp_yearly[[i]])
    total_icu[i, j] <- sum(icu_rate * colSums(df_results$inc_hosp_yearly[[i]]))
    total_deaths[i, j] <- sum(df_results$inc_deaths_yearly[[i]])
    total_vacc[i, j] <- sum(df_results$inc_vacc_yearly[[i]])
    total_symp_inf[i, j] <- sum(df_results$inc_symp_inf_yearly[[i]])

    # Calculating life years lost

    total_life_years_lost[i, j] <- sum(df_results$inc_deaths_yearly[[i]] * life_years_remaining)

    # Calculating costs

    total_vacc_cost[i, j] <- sum(costs_discount_rate * ((df_results$inc_vacc_yearly[[i]]) * state_costs$Vaccinated))
    total_inf_cost[i, j] <- sum(costs_discount_rate * ((df_results$inc_symp_inf_yearly[[i]]) * rep(state_costs$`Infectious symptomatic`, each = time_horizon)))
    total_hosp_cost[i, j] <- sum(costs_discount_rate * ((df_results$inc_hosp_yearly[[i]]) * rep(state_costs$Hospitalised, each = time_horizon)))

    # Calculating long Covid costs

    n_lc_recovered <- df_results$n_recovered_past_year_yearly[[i]] * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_vaccinated <- df_results$n_recovered_vaccinated_past_year_yearly[[i]] * red_long_covid * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_waned <- df_results$n_recovered_waned_past_year_yearly[[i]] * wan_long_covid * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_total <- n_lc_recovered + n_lc_recovered_vaccinated + n_lc_recovered_waned
    total_lc_cost[i, j] <- sum(costs_discount_rate * (n_lc_recovered_total * state_costs$`Long Covid`))

    # Calculating LC person-years

    total_lc_person_years[i, j] <- sum(n_lc_recovered_total) / 365

    # Calculating QALYs

    total_qalys[i, j] <- sum(qalys_discount_rate * ((df_results$n_healthy_yearly[[i]] -
                                                       n_lc_recovered_total) * rep(state_utilities$Healthy,
                                                                                   each = time_horizon) +
                                                      df_results$n_symptomatic_yearly[[i]] *
                                                      rep(state_utilities$`Infectious symptomatic`, each = time_horizon) +
                                                      df_results$n_hospitalised_yearly[[i]] *
                                                      rep(state_utilities$Hospitalised, each = time_horizon) +
                                                      n_lc_recovered_total * rep(state_utilities$`Long Covid`,
                                                                                 each = time_horizon)))

    # Calculating mean yearly hospital admissions
    
    mean_yearly_hosp[i, j] <- mean(rowSums(df_results$inc_hosp_yearly[[i]]))
    
    # Calculating mean yearly deaths
    
    mean_yearly_deaths[i, j] <- mean(rowSums(df_results$inc_deaths_yearly[[i]]))
    
    # Calculating mean number of vaccinations a person receives each year, counting only the eligible population (ACC. FOR uptake rate)

    vacc_yearly_coverage_by_age[i, ] <- (100 * colMeans(df_results$inc_vacc_yearly[[i]])) / age_dist

    # Calculating number of infections in each age group, and annual number of cases

    total_cases_age_strat[i, ] <- colSums(df_results$inc_inf_yearly[[i]])
    annual_cases[i, ] <- rowSums(df_results$inc_inf_yearly[[i]])

    # Calculating number of vaccinations in each age group

    total_vacc_age_strat[i, ] <- colSums(df_results$inc_vacc_yearly[[i]])

    # Calculating disease averted summaries

    cases_averted[i, j] <- total_infections[i, 1] - total_infections[i, j]
    hosp_averted[i, j] <- total_hosp[i, 1] - total_hosp[i, j]
    deaths_averted[i, j] <- total_deaths[i, 1] - total_deaths[i, j]
    life_years_saved[i, j] <- total_life_years_lost[i, 1] - total_life_years_lost[i, j]
    long_covid_py_averted[i, j] <- total_lc_person_years[i, 1] - total_lc_person_years[i, j]

    # Calculating costs with no discounting

    total_vacc_cost_no_discounting[i, j] <- sum((df_results$inc_vacc_yearly[[i]]) * state_costs$Vaccinated)
    total_inf_cost_no_discounting[i, j] <- sum((df_results$inc_symp_inf_yearly[[i]]) * rep(state_costs$`Infectious symptomatic`,
                                                                                           each = time_horizon))
    total_hosp_cost_no_discounting[i, j] <- sum((df_results$inc_hosp_yearly[[i]]) * rep(state_costs$Hospitalised, each = time_horizon))
    total_lc_cost_no_discounting[i, j] <- sum(n_lc_recovered_total * state_costs$`Long Covid`)
    
  }
  
  mean_cases_age_strat[, j] <- colMeans(total_cases_age_strat)
  mean_vacc_age_strat[, j] <- colMeans(total_vacc_age_strat)
  mean_annual_cases[, j] <- colMeans(annual_cases)
  mean_vacc_yearly_coverage_by_age[, j] <- colMeans(vacc_yearly_coverage_by_age)
  
}

# Calculating mean yearly hospital admissions and deaths, and 95% credible intervals, and max. weekly hospital admissions post-2030

mean_yearly_hosp_lower <- round(apply(mean_yearly_hosp, 2, function(x) quantile(x, probs = 0.025)))
mean_yearly_hosp_upper <- round(apply(mean_yearly_hosp, 2, function(x) quantile(x, probs = 0.975)))
mean_yearly_deaths_lower <- round(apply(mean_yearly_deaths, 2, function(x) quantile(x, probs = 0.025)))
mean_yearly_deaths_upper <- round(apply(mean_yearly_deaths, 2, function(x) quantile(x, probs = 0.975)))

mean_yearly_hosp <- round(colMeans(mean_yearly_hosp))
mean_yearly_deaths <- round(colMeans(mean_yearly_deaths))
max_weekly_hosp_no_vacc <- df_yearly_inf_hosp_deaths_lc[[1]]$max_weekly_hosp_post_2030

# Calculating total costs and collecting results

total_cost <- total_vacc_cost + total_inf_cost + total_hosp_cost + total_lc_cost

colnames(total_vacc_cost) <- colnames(total_inf_cost) <- colnames(total_hosp_cost) <- colnames(total_lc_cost) <- 
  colnames(total_cost) <- colnames(total_lc_person_years) <- colnames(total_qalys) <- colnames(total_life_years_lost) <- 
  colnames(total_infections) <- colnames(total_hosp) <- colnames(total_deaths) <- colnames(total_symp_inf) <- 
  colnames(total_icu) <- colnames(total_vacc) <- vacc_scenarios 

results <- list("Total vaccination cost" = total_vacc_cost,
                "Total infection cost" = total_inf_cost,
                "Total hospitalisation cost" = total_hosp_cost,
                "Total long Covid cost" = total_lc_cost,
                "Total cost" = total_cost,
                "Total long Covid person-years" = total_lc_person_years,
                "Total QALYs" = total_qalys,
                "Total life years lost" = total_life_years_lost,
                "Total infections" = total_infections,
                "Total hospitalisations" = total_hosp,
                "Total deaths" = total_deaths,
                "Total symptomatic infections" = total_symp_inf,
                "Total ICU admissions" = total_icu,
                "Total vaccinations" = total_vacc)

costs_no_discounting <- list("Total vaccination cost" = total_vacc_cost_no_discounting,
                             "Total infection cost" = total_inf_cost_no_discounting,
                             "Total hospitalisation cost" = total_hosp_cost_no_discounting,
                             "Total long Covid cost" = total_lc_cost_no_discounting)

# Disease averted

disease_averted_summary <- data.frame(scenario = vacc_scenarios,
                                      cases_mean = colMeans(cases_averted),
                                      cases_lower = apply(cases_averted, 2, function(x) quantile(x, probs = 0.025)),
                                      cases_upper = apply(cases_averted, 2, function(x) quantile(x, probs = 0.975)),
                                      cases_min = apply(cases_averted, 2, min),
                                      cases_max = apply(cases_averted, 2, max),
                                      hosp_mean = colMeans(hosp_averted),
                                      hosp_lower = apply(hosp_averted, 2, function(x) quantile(x, probs = 0.025)),
                                      hosp_upper = apply(hosp_averted, 2, function(x) quantile(x, probs = 0.975)),
                                      hosp_min = apply(hosp_averted, 2, min),
                                      hosp_max = apply(hosp_averted, 2, max),
                                      deaths_mean = colMeans(deaths_averted),
                                      deaths_lower = apply(deaths_averted, 2, function(x) quantile(x, probs = 0.025)),
                                      deaths_upper = apply(deaths_averted, 2, function(x) quantile(x, probs = 0.975)),
                                      deaths_min = apply(deaths_averted, 2, min),
                                      deaths_max = apply(deaths_averted, 2, max),
                                      life_years_mean = colMeans(life_years_saved),
                                      life_years_lower = apply(life_years_saved, 2, function(x) quantile(x, probs = 0.025)),
                                      life_years_upper = apply(life_years_saved, 2, function(x) quantile(x, probs = 0.975)),
                                      life_years_min = apply(life_years_saved, 2, min),
                                      life_years_max = apply(life_years_saved, 2, max),
                                      long_covid_py_mean = colMeans(long_covid_py_averted),
                                      long_covid_py_lower = apply(long_covid_py_averted, 2, function(x) quantile(x, probs = 0.025)),
                                      long_covid_py_upper = apply(long_covid_py_averted, 2, function(x) quantile(x, probs = 0.975)),
                                      long_covid_py_min = apply(long_covid_py_averted, 2, min),
                                      long_covid_py_max = apply(long_covid_py_averted, 2, max))
disease_averted_summary <- disease_averted_summary[-1,]

# Check for negative values in results

if(is.null(fixr::check_for_negative_values(c(as.numeric(unlist(results))))) == FALSE) {
  stop("Negative values in results")
}

# Plotting vaccination coverage by age

colnames(mean_vacc_yearly_coverage_by_age) <- vacc_scenarios
rownames(mean_vacc_yearly_coverage_by_age) <- 1:101

mean_vacc_yearly_coverage_by_age <- mean_vacc_yearly_coverage_by_age[, c(2, 3, 4, 7, 10, 11, 12, 13, 16, 19, 20, 21, 22, 25, 28)]
mean_vacc_yearly_coverage_by_age[, c(6:10)] <- mean_vacc_yearly_coverage_by_age[, c(6:10)] / 2
mean_vacc_yearly_coverage_by_age[, 11:15] <- mean_vacc_yearly_coverage_by_age[, 11:15] * 2

mean_vacc_yearly_coverage_by_age_longer <- data.frame(mean_vacc_yearly_coverage_by_age) |>
  pivot_longer(everything(),
               names_to = "Strategy",
               values_to = "Coverage")
mean_vacc_yearly_coverage_by_age_longer$Age <- rep(1:101, nrow(mean_vacc_yearly_coverage_by_age_longer) / 101)
mean_vacc_yearly_coverage_by_age_longer$Uptake <- rep(vacc_uptake_rate * 100, nrow(mean_vacc_yearly_coverage_by_age_longer) / 101)

vacc_coverage_by_age_plot <- ggplot(data = mean_vacc_yearly_coverage_by_age_longer) +
  geom_line(aes(x = Age, y = Coverage, colour = Strategy, group = Strategy), linewidth = 1) +
  scale_colour_manual(values = c(brewer.pal(9, "Blues")[c(4:8)], brewer.pal(9, "Greens")[c(4:8)], 
                                 brewer.pal(9, "Reds")[c(4:8)], "grey", "lightgrey")) +
  scale_y_continuous(limits = c(0, 100)) +
  geom_line(aes(x = Age, y = Uptake, linetype = "Vaccine uptake rate"), colour = "black", linewidth = 1.5) +
  labs(y = "Vaccine coverage (%)", colour = "Vaccination strategy", linetype = "") +
  scale_linetype_manual(values = c("Vaccine uptake rate" = "dotted")) +
  guides(colour = guide_legend(order = 1, ncol = 3),  
         linetype = guide_legend(order = 2, override.aes = list(colour = "black"))) +
  theme(legend.position = c(0.12, 0.75))

print(vacc_coverage_by_age_plot)

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Vaccine coverage plot.jpeg"), 
       plot = vacc_coverage_by_age_plot, width = 10, height = 7)

# Plotting mean yearly deaths for each vacc strategy against 2023

df_mean_yearly_deaths <- data.frame(Deaths = mean_yearly_deaths, Strategy = vacc_scenarios)
mean_yearly_deaths_plot <- ggplot(data = df_mean_yearly_deaths, aes(x = Strategy, y = Deaths)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = deaths_data_total["Total"]), colour = "red", linewidth = 1.5, linetype = "dashed") +
  theme(legend.position = "none", 
        axis.text.x = element_text(hjust = 1, angle =  45, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  labs(x = "Vaccination strategy",
       y = "Mean yearly deaths") +
  annotate(geom = "label", x = 2, y = deaths_data_total["Total"], label = "2023", size = 5)

print(mean_yearly_deaths_plot)

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Plot 1. Mean yearly deaths.jpeg"), 
       plot = mean_yearly_deaths_plot, width = 8, height = 5)

# Carrying out cost-effectiveness evaluation

# Run BCEA package

BCEA <- bcea(eff = results$`Total QALYs`, cost = results$`Total cost`, ref = 1, plot = FALSE, 
             interventions = vacc_scenarios, kmax = 30000)
summary_BCEA <- summary(BCEA, wtp = 20000)
BCEA$ICER

# Cost-effectiveness frontier

# Calculating dominant strategies

cef_data_means <- data.frame(icer = BCEA$ICER, 
                             effectiveness = -colMeans(BCEA$delta_e),
                             cost = -colMeans(BCEA$delta_c),
                             strategy = vacc_scenarios[-1],
                             absolute_domination = rep("", length(vacc_scenarios[-1])),
                             extended_domination = rep("", length(vacc_scenarios[-1])),
                             dominant = rep("", length(vacc_scenarios[-1])))

cef_data_means <- cef_data_means[order(cef_data_means$cost), ]

for (i in 2:length(vacc_scenarios[-1])) {
  if (any(cef_data_means$effectiveness[i] < cef_data_means$effectiveness[1:(i-1)])) {
    cef_data_means$absolute_domination[i] <- "Yes"
  }
}

cef_data_means_undominated <- cef_data_means

repeat {
  cef_data_means_undominated_previous <- cef_data_means_undominated
  cef_data_means_undominated <-  cef_data_means |> filter(absolute_domination != "Yes", extended_domination != "Yes")
  cef_data_means_undominated <- cef_data_means_undominated[order(cef_data_means_undominated$effectiveness), ]
  cef_data_means_undominated$icer <- rep(0, nrow(cef_data_means_undominated))
  for (i in 2:nrow(cef_data_means_undominated)) {
    cef_data_means_undominated$icer[i] <- (cef_data_means_undominated$cost[i] - cef_data_means_undominated$cost[i - 1]) /
      (cef_data_means_undominated$effectiveness[i] - 
         cef_data_means_undominated$effectiveness[i - 1])
  }
  for (i in 1:(nrow(cef_data_means_undominated) - 1)) {
    if (cef_data_means_undominated$icer[i] > cef_data_means_undominated$icer[i + 1]) {
      cef_data_means_undominated$extended_domination[i] <- "Yes"
      cef_data_means[cef_data_means$strategy == cef_data_means_undominated$strategy[i], "extended_domination"] <- "Yes"
    }
  }
  cef_data_means_undominated <- cef_data_means_undominated |>
    filter(extended_domination != "Yes")
  if (identical(cef_data_means_undominated, cef_data_means_undominated_previous)) {
    cef_data_means_undominated <- data.frame(cef_data_means_undominated)
    break
  }
}

for (i in 1:length(vacc_scenarios[-1])) {
  if (cef_data_means$absolute_domination[i] == "" & cef_data_means$extended_domination[i] == "") {
    cef_data_means$dominant[i] <- "Yes"
  } else {
    cef_data_means$dominant[i] <- "No"
  }
}

cef_data_dominant <- cef_data_means[cef_data_means$dominant == "Yes",]

cef_plot <- ggplot(data = cef_data_means, aes(x = effectiveness, y = cost, colour = dominant, label = strategy)) +
  geom_point() +
  geom_line(data = cef_data_dominant, aes(x = effectiveness, y = cost, group = 1)) +
  labs(colour = "Dominant strategy?",
       x = "Incremental QALYs",
       y = "Incremental cost (£)") +
  scale_x_continuous(labels = unit_format(unit = "million", scale = 1e-6)) +
  scale_y_continuous(labels = unit_format(unit = "billion", scale = 1e-9)) +
  geom_label_repel(aes(label = strategy), box.padding   = 0.35, point.padding = 0.5, max.overlaps = 40,
                   segment.color = 'grey50', min.segment.length = 0.1, show.legend = FALSE) +
  theme_light(base_size = 13) +
  theme(legend.position = c(0.2, 0.9))

print(cef_plot)

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Plot 2. Cost-effectiveness frontier.jpeg"),
       plot = cef_plot, width = 8, height = 6)

# Check for negative incremental qalys

if(is.null(fixr::check_for_negative_values(-BCEA$delta_e)) == FALSE) {
  stop("Negative incremental qalys")
}

# Calculate 95% confidence intervals for the ICER

icer_simulations <- BCEA$delta_c / BCEA$delta_e         # delta_c are the incremental costs and delta_e incremental qalys for each simulation
icer_95_ci <- apply(icer_simulations, 2, quantile, probs = c(0.025, 0.975))

# Calculating ICER 95% CI

icer_mean <- c(0, round(as.numeric(BCEA$ICER)))
icer_95_ci_lower <- c(0, round(icer_95_ci[1,]))
icer_95_ci_upper <- c(0, round(icer_95_ci[2,]))

# Calculate the incremental net monetary benefit, with a willingness to pay threshold of £20,000 per QALY
# If the INMB is positive then the intervention is cost-effective

inmb <- (-BCEA$delta_e) * 20000 - (-BCEA$delta_c)
inmb <- cbind(rep(0, ncol(par_fit)), inmb)
inmb_mean <- colMeans(inmb)
inmb_95_ci <- apply(inmb, 2, quantile, probs = c(0.025, 0.975))
inmb_95_ci_lower <- round(inmb_95_ci[1,])
inmb_95_ci_upper <- round(inmb_95_ci[2,])

# Plotting CEP for multiple vaccination strategies

incremental_costs <- -BCEA$delta_c |> 
  pivot_longer(cols = vacc_scenarios[2]:vacc_scenarios[length(vacc_scenarios)],
               names_to = "Scenario",
               values_to = "Costs")
incremental_costs <- data.frame(incremental_costs)

incremental_qalys <- -BCEA$delta_e |> 
  pivot_longer(cols = vacc_scenarios[2]:vacc_scenarios[length(vacc_scenarios)],
               names_to = "Scenario",
               values_to = "QALYs")
incremental_qalys <- data.frame(incremental_qalys)

vacc_frequency <- str_sub(incremental_costs$Scenario, 1, 1)
vacc_groups <- str_sub(incremental_costs$Scenario, 2, 2)

inc_costs_qalys <- cbind(vacc_frequency, vacc_groups, incremental_costs, incremental_qalys[,2])
colnames(inc_costs_qalys)[c(1,2,5)] <- c("Frequency", "Groups", "QALYs")

cost_effectiveness_plane <- ggplot(inc_costs_qalys, aes(x = QALYs, y = Costs)) +
  geom_point(size = 0.8) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  scale_y_continuous(labels = unit_format(unit = "B", scale = 1e-9)) +
  geom_abline(aes(intercept = 20000, slope = 20000), color = "orange", linewidth = 0.8) +
  geom_abline(aes(intercept = 30000, slope = 30000), color = "red", linewidth = 0.8) +
  facet_wrap(~Scenario, ncol = 3, scales = 'free') +
  labs(y = "Incremental costs",
       x = "Incremental QALYs") +
  theme(strip.text.x = element_text(size = 10, face = "bold", colour = "black"),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 13))

label_plot <- ggplot() +
  annotate(geom = "richtext", x = 2, y = 2,
           label = " <span style='color:orange;'>Orange</span>: £20,000 threshold<br><span style='color:red;'>Red</span>: £30,000 threshold",
           size = 5, fill = NA, label.color = NA) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", color = "white"))

cost_effectiveness_plane <- cost_effectiveness_plane + label_plot + plot_layout(nrow = 2, heights = c(16,1))

print(cost_effectiveness_plane)

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Plot 3. Cost-effectiveness plane.jpeg"),
       plot = cost_effectiveness_plane, height = 16, width = 10)

# Plot the cost-effectiveness acceptability curve: probability that vaccination is CE for various WTP thresholds

wtp_thresholds <- seq(0, 40000, by = 100)
prob_cost_effective <- matrix(1, nrow = length(wtp_thresholds), ncol = ncol(icer_simulations))

for (i in 1:length(wtp_thresholds)) {
  prob_cost_effective[i, ] <- colSums(icer_simulations < wtp_thresholds[i]) / nrow(icer_simulations)
}

prob_cost_effective_20k <- round(c(0, prob_cost_effective[201,]) * 100)
prob_cost_effective_30k <- round(c(0, prob_cost_effective[301,]) * 100)

prob_cost_effective <- data.frame(prob_cost_effective)
colnames(prob_cost_effective) <- vacc_scenarios[-1]

prob_cost_effective <- prob_cost_effective |> 
  pivot_longer(cols = vacc_scenarios[2]:vacc_scenarios[length(vacc_scenarios)],
               names_to = "Scenario",
               values_to = "Probability")

prob_cost_effective <- data.frame(prob_cost_effective)

prob_cost_effective <- cbind(rep(wtp_thresholds, ncol(icer_simulations)), prob_cost_effective)
colnames(prob_cost_effective)[1] <- "Threshold"

# Probability that the strategy is cost saving:

prob_cost_saving <- c(0, round(colSums(BCEA$delta_c > 0) / nrow(BCEA$delta_c) * 100))

# Plotting CEAC for multiple vaccination strategies

blue_12 <- colorRampPalette(brewer.pal(9, "Blues"))(12)
green_12 <- colorRampPalette(brewer.pal(9, "Greens"))(12)
red_12 <- colorRampPalette(brewer.pal(9, "Reds"))(12)

cost_effectiveness_acceptability_curve <- ggplot(prob_cost_effective, aes(x = Threshold, y = Probability, colour = Scenario)) + 
  geom_smooth(method = "loess", span = 0.1, se = FALSE) +
  labs(colour = "Vaccination strategy") +
  scale_color_manual(values = c(blue_12[4:12], green_12[4:12], red_12[4:12])) +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 12))

print(cost_effectiveness_acceptability_curve)

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Plot 4. Cost effectiveness acceptability curve.jpeg"), 
       plot = cost_effectiveness_acceptability_curve, width = 8, height = 5)

# Calculating means and confidence intervals for all results

results_mean <- results_95_ci <- list()
results_merged <- rep(NA, length(vacc_scenarios))
for (i in 1:length(results)) {
  results_mean[[i]] <- as.data.frame(apply(results[[i]], 2, mean))
  names(results_mean[[i]]) <- paste0(names(results)[i], " mean")
  results_95_ci[[i]] <- apply(results[[i]], 2, quantile, probs = c(0.025, 0.975))
  results_merged <- rbind(results_merged, t(results_mean[[i]]), results_95_ci[[i]])
}
results_merged <- round(results_merged[-1,])

# Merging results with 0% discounted results so first four items are cost of vacc, inf, hosp and lc without discounting to 
# calculate the annual mean incremental results

results_without_with_discounting <- merge.list(costs_no_discounting, results)

incremental_results_mean <- incremental_results_95_ci <- list()
incremental_results_merged <- rep(NA, length(vacc_scenarios))
for (i in 1:length(results_without_with_discounting)) {
  incremental_results_mean[[i]] <- as.data.frame(apply(((results_without_with_discounting[[i]] - 
                                                           results_without_with_discounting[[i]][,1])), 2, mean))
  names(incremental_results_mean[[i]]) <- paste0(names(results_without_with_discounting)[i], " mean")
  incremental_results_95_ci[[i]] <- apply(((results_without_with_discounting[[i]] - 
                                              results_without_with_discounting[[i]][,1])), 2, quantile, probs = c(0.025, 0.975))
  incremental_results_merged <- rbind(incremental_results_merged, t(incremental_results_mean[[i]]), incremental_results_95_ci[[i]])
}

incremental_results_merged <- round(incremental_results_merged[-1,])
names(results_mean) <- names(results_95_ci) <- names(results)

# Plotting ranked strategies for various epidemiological outcomes

df_disease_averted <- disease_averted_summary[c(-4, -5, -7, -8, -13, -14, -16, -17, -22, -23, -25, -26),]
df_disease_averted$rank.cases <- rank(-df_disease_averted$cases_mean)
df_disease_averted$rank.hosp <- rank(-df_disease_averted$hosp_mean)
df_disease_averted$rank.deaths <- rank(-df_disease_averted$deaths_mean)
df_disease_averted$rank.lys <- rank(-df_disease_averted$life_years_mean)
df_disease_averted$rank.lcpy <- rank(-df_disease_averted$long_covid_py_mean)
df_disease_averted_rank <- df_disease_averted[, c(1, 27:31)]

df_disease_averted_rank <- df_disease_averted_rank |>
  pivot_longer(cols = (rank.cases:rank.lcpy),
               names_to = "Outcome",
               values_to = "Rank")
df_disease_averted_rank$Outcome <- c(rep("Cases", nrow(df_disease_averted)), 
                                     rep("Hospitalisations", nrow(df_disease_averted)), 
                                     rep("Deaths", nrow(df_disease_averted)), 
                                     rep("Life years saved", nrow(df_disease_averted)),
                                     rep("Long Covid person-years", nrow(df_disease_averted)))
df_disease_averted_rank$Outcome_f <- factor(df_disease_averted_rank$Outcome, 
                                            levels = c("Cases", "Hospitalisations", "Deaths", "Life years saved",
                                                       "Long Covid person-years"))

df_disease_averted_rank$Names <- rep(1:nrow(df_disease_averted), 5)
labels <- setNames(c(vacc_scenarios[c(-1, -5, -6, -8, -9, -14, -15, -17, -18, -23, -24, -26, -27)]), 1:nrow(df_disease_averted))

ranked_strategies_plot <- ggplot(df_disease_averted_rank, aes(x = Rank, y = Names, colour = Outcome_f, linetype = Outcome_f)) +
  geom_bump(linewidth = 1.5) +
  labs(x = "Rank (1 = best)",
       y = "Vaccination strategy",
       colour = "Outcome",
       linetype = "Outcome") +
  scale_y_continuous(breaks = seq(1, nrow(df_disease_averted), by = 1), labels = labels, n.breaks = nrow(df_disease_averted)) +
  scale_x_continuous(breaks = seq(1, nrow(df_disease_averted), by = 1), labels = 1:nrow(df_disease_averted)) +
  scale_color_manual(values = c((hue_pal()(4))[1], "green4", (hue_pal()(4))[3:4], "navy")) +
  theme(panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.width = unit(0.7, 'cm'), #change legend key width
        legend.position = c(0.18, 0.85)) +
  geom_text(aes(label = Rank), vjust = -0.8, size = 4.5, show.legend = FALSE)

print(ranked_strategies_plot)

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Plot 5. Ranked strategies by outcome.jpeg"), 
       plot = ranked_strategies_plot, width = 8, height = 6)

# Plotting disease averted for each strategy as a panel

disease_averted_summary_long <- disease_averted_summary |> 
  pivot_longer(everything(),
               names_to = c("Outcome_orig", ".value"),
               names_pattern = "(.+)_(.+)")
disease_averted_summary_long <- disease_averted_summary_long[1:((length(vacc_scenarios) - 1) * 5)]
disease_averted_summary_long$Scenario <- rep(vacc_scenarios[-1], nrow(disease_averted_summary_long) / length(vacc_scenarios[-1]))
disease_averted_summary_long$Outcome <- factor(disease_averted_summary_long$Outcome_orig, 
                                               levels = c("cases", "hosp", "deaths", "life_years", "long_covid_py"))

disease_averted_labels <- c(`cases` = "(a) Cases averted",
                            `hosp` = "(b) Hospitalisations averted", 
                            `deaths` = "(c) Deaths averted", 
                            `life_years` = "(d) Life years saved",
                            `long_covid_py` = "(e) Long Covid person-years averted")

disease_averted_panel_plot <- ggplot(data = disease_averted_summary_long, aes(x = Scenario, fill = Outcome)) + 
  geom_boxplot(aes(ymin = min, lower = lower, middle = mean, upper = upper, ymax = max),
               stat = "identity") +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  facet_wrap(~Outcome, scales = "free_y", labeller = as_labeller(disease_averted_labels), nrow = 3) +
  labs(fill = "Outcome",
       x = "Vaccination strategy",
       y = "Outcome averted") +
  scale_fill_manual(values = c(brewer.pal(5, "Set2")), 
                    labels = c("Cases averted", "Hospitalisations averted", "Deaths averted", "Life years saved",
                               "Long Covid person-years
                            averted")) +
  theme(strip.text.x = element_text(size = 13, colour = "black"),
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13))

print(disease_averted_panel_plot)

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Plot 6. Disease averted panel plot.jpeg"),
       plot = disease_averted_panel_plot, height = 9, width = 9)

# Disease averted per vaccination given

df_disease_averted_per_1000_doses <- data.frame(matrix(NA, nrow = 5, ncol = (length(vacc_scenarios) - 1)))

df_disease_averted_per_1000_doses[1,] <- 
  1000 * disease_averted_summary$cases_mean / results_mean$`Total vaccinations`[-1, ] / 
  median(  1000 * disease_averted_summary$cases_mean / results_mean$`Total vaccinations`[-1, ])
df_disease_averted_per_1000_doses[2,] <- 
  1000 * disease_averted_summary$hosp_mean / results_mean$`Total vaccinations`[-1, ] / 
  median(  1000 * disease_averted_summary$hosp_mean / results_mean$`Total vaccinations`[-1, ])
df_disease_averted_per_1000_doses[3,] <- 
  1000 * disease_averted_summary$deaths_mean / results_mean$`Total vaccinations`[-1, ] / 
  median(  1000 * disease_averted_summary$deaths_mean / results_mean$`Total vaccinations`[-1, ])
df_disease_averted_per_1000_doses[4,] <- 
  1000 * disease_averted_summary$life_years_mean / results_mean$`Total vaccinations`[-1, ] / 
  median(  1000 * disease_averted_summary$life_years_mean / results_mean$`Total vaccinations`[-1, ])
df_disease_averted_per_1000_doses[5,] <- 
  1000 * disease_averted_summary$long_covid_py_mean / results_mean$`Total vaccinations`[-1, ] / 
  median(  1000 * disease_averted_summary$long_covid_py_mean / results_mean$`Total vaccinations`[-1, ])

rownames(df_disease_averted_per_1000_doses) <- c("Cases averted per 1000 doses", "Hospital admissions averted per 1000 doses", 
                                                 "Deaths averted per 1000 doses", "Life years saved per 1000 doses", 
                                                 "Long Covid person-years averted per 1000 doses")
colnames(df_disease_averted_per_1000_doses) <- vacc_scenarios[-1]

df_disease_averted_per_1000_doses_longer <- df_disease_averted_per_1000_doses |> 
  pivot_longer(cols = 1:ncol(df_disease_averted_per_1000_doses),
               names_to = "Scenario",
               values_to = "Disease averted")

df_disease_averted_per_1000_doses_longer$Outcome <- rep(c("Cases", "Hospitalisations", 
                                                          "Deaths", "Life years saved", 
                                                          "Long Covid person-years"), (length(vacc_scenarios) - 1))
df_disease_averted_per_1000_doses_longer$Outcome_f <- factor(df_disease_averted_per_1000_doses_longer$Outcome, 
                                                             levels = c("Cases", "Hospitalisations", "Deaths", "Life years saved",
                                                                        "Long Covid person-years"))

df_disease_averted_per_1000_doses_longer <- df_disease_averted_per_1000_doses_longer[1:45, ]

disease_averted_per_1000_doses_plot <- ggplot(data = df_disease_averted_per_1000_doses_longer, 
                                              aes(x = Scenario, y = `Disease averted`, colour = Outcome_f)) +
  geom_point(aes(group = Outcome_f), shape = 1, stroke = 1.5, size = 2) +
  labs(y = "Outcome averted per 1,000 doses (normalised)",
       x = "Vaccination strategy",
       colour = "Outcome") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = c(0.5, 0.8))

print(disease_averted_per_1000_doses_plot)

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Plot 7. Disease averted per 1000 doses plot.jpeg"),
       plot = disease_averted_per_1000_doses_plot, width = 9, height = 6)

# Creating results table

df_disease_averted_per_1000_doses <- cbind(rep(0, 5), df_disease_averted_per_1000_doses)

colnames(df_disease_averted_per_1000_doses) <- colnames(results_merged) <- colnames(incremental_results_merged) <- 
  vacc_scenarios

full_epi_results <- rbind.data.frame(results_merged[22:42,], incremental_results_merged[22:36,], df_disease_averted_per_1000_doses)

full_ce_results <- rbind.data.frame("ICER" = icer_mean, "ICER lower bound" = icer_95_ci_lower, "ICER upper bound" = icer_95_ci_upper, 
                                    "INMB" = inmb_mean, "INMB lower bound" = inmb_95_ci_lower, "INMB upper bound" = inmb_95_ci_upper,
                                    rbind.data.frame(results_merged[1:21,], incremental_results_merged[1:21,]),
                                    "Probability of cost effectiveness, £20k WTP" = prob_cost_effective_20k, 
                                    "Probability of cost effectiveness, £30k WTP" = prob_cost_effective_30k, 
                                    "Probability cost saving" = prob_cost_saving)

colnames(full_epi_results) <- colnames(full_ce_results) <- c("BASELINE, comparator", vacc_scenarios[-1])
model_output <- list("Epi results" = full_epi_results, "CE results" = full_ce_results, Parameters = par_fit)
write.xlsx((model_output), paste0("Results/", as.character(Sys.Date()), " Results table, baseline.xlsx"), rowNames = TRUE)


############################################################################################

# Calculating disease outcomes by age

# Calculating cases by age

mean_cases_age_strat <- matrix(NA, ncol = length(vacc_scenarios), nrow = 101)

for (j in 1:length(vacc_scenarios)) {
  total_cases_age_strat <- matrix(NA, nrow = ncol(par_fit), ncol = 101)
  for (i in 1:ncol(par_fit)) {
    # Set dataframe for vaccination scenario
    df_results <- df_yearly_inf_hosp_deaths_lc[[j]]
    # Calculating number of infections in each age group
    total_cases_age_strat[i, ] <- colSums(df_results$inc_inf_yearly[[i]])
  }
  mean_cases_age_strat[, j] <- colMeans(total_cases_age_strat)
}

mean_cases_averted_age_strat <- (mean_cases_age_strat[, 1] - mean_cases_age_strat)[, -1]
mean_cases_averted_0_17 <- colSums(mean_cases_age_strat[1:18, 1] - mean_cases_age_strat[1:18, ])[-1]
mean_cases_averted_18_49 <- colSums(mean_cases_age_strat[19:50, 1] - mean_cases_age_strat[19:50, ])[-1]
mean_cases_averted_50_64 <- colSums(mean_cases_age_strat[51:65, 1] - mean_cases_age_strat[51:65, ])[-1]
mean_cases_averted_65_79 <- colSums(mean_cases_age_strat[66:80, 1] - mean_cases_age_strat[66:80, ])[-1]
mean_cases_averted_80_ <- colSums(mean_cases_age_strat[81:101, 1] - mean_cases_age_strat[81:101, ])[-1]

cases_averted_by_age <- t(data.frame(mean_cases_averted_0_17 = mean_cases_averted_0_17,
                                     mean_cases_averted_18_49 = mean_cases_averted_18_49,
                                     mean_cases_averted_50_64 = mean_cases_averted_50_64,
                                     mean_cases_averted_65_79 = mean_cases_averted_65_79,
                                     mean_cases_averted_80_ = mean_cases_averted_80_))
colnames(cases_averted_by_age) <- vacc_scenarios[-1]

cases_averted_by_age_longer <- data.frame(cases_averted_by_age) |>
  pivot_longer(everything(),
               names_to = "Strategy",
               values_to = "Cases averted")

cases_averted_by_age_longer$Age_group <- rep(c("0 to 17", "18 to 49", "50 to 64", "65 to 79", "80+"), length(vacc_scenarios[-1]))
cases_averted_by_age_longer_subset <- cases_averted_by_age_longer[grepl("1|2|3|6|9", cases_averted_by_age_longer$Strategy), ]

cases_averted_by_age_plot <- ggplot(data = cases_averted_by_age_longer_subset, 
                                    aes(x = Strategy, y = `Cases averted`, fill = Age_group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(brewer.pal(5, "Set2")))

# Calculating hosp by age

mean_hosp_age_strat <- matrix(NA, ncol = length(vacc_scenarios), nrow = 101)

for (j in 1:length(vacc_scenarios)) {
  total_hosp_age_strat <- matrix(NA, nrow = ncol(par_fit), ncol = 101)
  for (i in 1:ncol(par_fit)) {
    # Set dataframe for vaccination scenario
    df_results <- df_yearly_inf_hosp_deaths_lc[[j]]
    # Calculating number of hosp in each age group
    total_hosp_age_strat[i, ] <- colSums(df_results$inc_hosp_yearly[[i]])
  }
  mean_hosp_age_strat[, j] <- colMeans(total_hosp_age_strat)
}

mean_hosp_averted_age_strat <- (mean_hosp_age_strat[, 1] - mean_hosp_age_strat)[, -1]
mean_hosp_averted_0_17 <- colSums(mean_hosp_age_strat[1:18, 1] - mean_hosp_age_strat[1:18, ])[-1]
mean_hosp_averted_18_49 <- colSums(mean_hosp_age_strat[19:50, 1] - mean_hosp_age_strat[19:50, ])[-1]
mean_hosp_averted_50_64 <- colSums(mean_hosp_age_strat[51:65, 1] - mean_hosp_age_strat[51:65, ])[-1]
mean_hosp_averted_65_79 <- colSums(mean_hosp_age_strat[66:80, 1] - mean_hosp_age_strat[66:80, ])[-1]
mean_hosp_averted_80_ <- colSums(mean_hosp_age_strat[81:101, 1] - mean_hosp_age_strat[81:101, ])[-1]

hosp_averted_by_age <- t(data.frame(mean_hosp_averted_0_17 = mean_hosp_averted_0_17,
                                    mean_hosp_averted_18_49 = mean_hosp_averted_18_49,
                                    mean_hosp_averted_50_64 = mean_hosp_averted_50_64,
                                    mean_hosp_averted_65_79 = mean_hosp_averted_65_79,
                                    mean_hosp_averted_80_ = mean_hosp_averted_80_))
colnames(hosp_averted_by_age) <- vacc_scenarios[-1]

hosp_averted_by_age_longer <- data.frame(hosp_averted_by_age) |>
  pivot_longer(everything(),
               names_to = "Strategy",
               values_to = "Hospitalisations averted")

hosp_averted_by_age_longer$Age_group <- rep(c("0 to 17", "18 to 49", "50 to 64", "65 to 79", "80+"), length(vacc_scenarios[-1]))
hosp_averted_by_age_longer_subset <- hosp_averted_by_age_longer[grepl("1|2|3|6|9", hosp_averted_by_age_longer$Strategy), ]

hosp_averted_by_age_plot <- ggplot(data = hosp_averted_by_age_longer_subset, 
                                   aes(x = Strategy, y = `Hospitalisations averted`, fill = Age_group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(brewer.pal(5, "Set2")))

# Calculating deaths by age

mean_deaths_age_strat <- matrix(NA, ncol = length(vacc_scenarios), nrow = 101)
mean_deaths_yearly <- mean_deaths_yearly_0_20 <- mean_deaths_yearly_21_40 <- mean_deaths_yearly_41_60 <- mean_deaths_yearly_61_80 <-
  mean_deaths_yearly_81_100 <- matrix(NA, ncol = length(vacc_scenarios), nrow = 100)

for (j in 1:length(vacc_scenarios)) {
  total_deaths_age_strat <- matrix(NA, nrow = ncol(par_fit), ncol = 101)
  total_deaths_yearly <- total_deaths_yearly_0_20 <- total_deaths_yearly_21_40 <- total_deaths_yearly_41_60 <- total_deaths_yearly_61_80 <-
    total_deaths_yearly_81_100 <- matrix(NA, nrow = ncol(par_fit), ncol = 100)
  for (i in 1:ncol(par_fit)) {
    # Set dataframe for vaccination scenario
    df_results <- df_yearly_inf_hosp_deaths_lc[[j]]
    # Calculating number of deaths in each age group
    total_deaths_age_strat[i, ] <- colSums(df_results$inc_deaths_yearly[[i]])
    total_deaths_yearly[i, ] <- rowSums(df_results$inc_deaths_yearly[[i]])
    total_deaths_yearly_0_20[i, ] <- rowSums(df_results$inc_deaths_yearly[[i]][, 1:21])
    total_deaths_yearly_21_40[i, ] <- rowSums(df_results$inc_deaths_yearly[[i]][, 22:41])
    total_deaths_yearly_41_60[i, ] <- rowSums(df_results$inc_deaths_yearly[[i]][, 42:61])
    total_deaths_yearly_61_80[i, ] <- rowSums(df_results$inc_deaths_yearly[[i]][, 62:81])
    total_deaths_yearly_81_100[i, ] <- rowSums(df_results$inc_deaths_yearly[[i]][, 82:101])
    
  }
  mean_deaths_age_strat[, j] <- colMeans(total_deaths_age_strat)
  mean_deaths_yearly[, j] <- colMeans(total_deaths_yearly)
  mean_deaths_yearly_0_20[, j] <- colMeans(total_deaths_yearly_0_20)
  mean_deaths_yearly_21_40[, j] <- colMeans(total_deaths_yearly_21_40)
  mean_deaths_yearly_41_60[, j] <- colMeans(total_deaths_yearly_41_60)
  mean_deaths_yearly_61_80[, j] <- colMeans(total_deaths_yearly_61_80)
  mean_deaths_yearly_81_100[, j] <- colMeans(total_deaths_yearly_81_100)
  
}

mean_deaths_yearly_by_age <- cbind(mean_deaths_yearly_0_20, mean_deaths_yearly_21_40, mean_deaths_yearly_41_60, mean_deaths_yearly_61_80,
                                   mean_deaths_yearly_81_100)

mean_deaths_averted_age_strat <- (mean_deaths_age_strat[, 1] - mean_deaths_age_strat)[, -1]
mean_deaths_averted_0_17 <- colSums(mean_deaths_age_strat[1:18, 1] - mean_deaths_age_strat[1:18, ])[-1]
mean_deaths_averted_18_49 <- colSums(mean_deaths_age_strat[19:50, 1] - mean_deaths_age_strat[19:50, ])[-1]
mean_deaths_averted_50_64 <- colSums(mean_deaths_age_strat[51:65, 1] - mean_deaths_age_strat[51:65, ])[-1]
mean_deaths_averted_65_79 <- colSums(mean_deaths_age_strat[66:80, 1] - mean_deaths_age_strat[66:80, ])[-1]
mean_deaths_averted_80_ <- colSums(mean_deaths_age_strat[81:101, 1] - mean_deaths_age_strat[81:101, ])[-1]

deaths_averted_by_age <- t(data.frame(mean_deaths_averted_0_17 = mean_deaths_averted_0_17,
                                      mean_deaths_averted_18_49 = mean_deaths_averted_18_49,
                                      mean_deaths_averted_50_64 = mean_deaths_averted_50_64,
                                      mean_deaths_averted_65_79 = mean_deaths_averted_65_79,
                                      mean_deaths_averted_80_ = mean_deaths_averted_80_))
colnames(deaths_averted_by_age) <- vacc_scenarios[-1]

deaths_averted_by_age_longer <- data.frame(deaths_averted_by_age) |>
  pivot_longer(everything(),
               names_to = "Strategy",
               values_to = "Deaths averted")

deaths_averted_by_age_longer$Age_group <- rep(c("0 to 17", "18 to 49", "50 to 64", "65 to 79", "80+"), length(vacc_scenarios[-1]))
deaths_averted_by_age_longer_subset <- deaths_averted_by_age_longer[grepl("1|2|3|6|9", deaths_averted_by_age_longer$Strategy), ]

deaths_averted_by_age_plot <- ggplot(data = deaths_averted_by_age_longer_subset, 
                                     aes(x = Strategy, y = `Deaths averted`, fill = Age_group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(brewer.pal(5, "Set2")))

# Calculating long Covid by age

mean_long_covid_age_strat <- matrix(NA, ncol = length(vacc_scenarios), nrow = 101)

for (j in 1:length(vacc_scenarios)) {
  
  total_lc_person_years_age_strat <- matrix(NA, nrow = ncol(par_fit), ncol = 101)
  
  for (i in 1:ncol(par_fit)) {
    
    # Set dataframe for vaccination scenario
    
    df_results <- df_yearly_inf_hosp_deaths_lc[[j]]
    
    # Calculating long Covid person-years by age
    
    n_lc_recovered <- df_results$n_recovered_past_year_yearly[[i]] * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_vaccinated <- df_results$n_recovered_vaccinated_past_year_yearly[[i]] * red_long_covid * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_waned <- df_results$n_recovered_waned_past_year_yearly[[i]] * wan_long_covid * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_total <- n_lc_recovered + n_lc_recovered_vaccinated + n_lc_recovered_waned
    
    total_lc_person_years_age_strat[i, ] <- colSums(n_lc_recovered_total) / 365
    
  }
  
  mean_long_covid_age_strat[, j] <- colMeans(total_lc_person_years_age_strat)
  
}

mean_long_covid_averted_age_strat <- (mean_long_covid_age_strat[, 1] - mean_long_covid_age_strat)[, -1]
mean_long_covid_averted_0_17 <- colSums(mean_long_covid_age_strat[1:18, 1] - mean_long_covid_age_strat[1:18, ])[-1]
mean_long_covid_averted_18_49 <- colSums(mean_long_covid_age_strat[19:50, 1] - mean_long_covid_age_strat[19:50, ])[-1]
mean_long_covid_averted_50_64 <- colSums(mean_long_covid_age_strat[51:65, 1] - mean_long_covid_age_strat[51:65, ])[-1]
mean_long_covid_averted_65_79 <- colSums(mean_long_covid_age_strat[66:80, 1] - mean_long_covid_age_strat[66:80, ])[-1]
mean_long_covid_averted_80_ <- colSums(mean_long_covid_age_strat[81:101, 1] - mean_long_covid_age_strat[81:101, ])[-1]

long_covid_averted_by_age <- t(data.frame(mean_long_covid_averted_0_17 = mean_long_covid_averted_0_17,
                                          mean_long_covid_averted_18_49 = mean_long_covid_averted_18_49,
                                          mean_long_covid_averted_50_64 = mean_long_covid_averted_50_64,
                                          mean_long_covid_averted_65_79 = mean_long_covid_averted_65_79,
                                          mean_long_covid_averted_80_ = mean_long_covid_averted_80_))
colnames(long_covid_averted_by_age) <- vacc_scenarios[-1]

long_covid_averted_by_age_longer <- data.frame(long_covid_averted_by_age) |>
  pivot_longer(everything(),
               names_to = "Strategy",
               values_to = "Long Covid averted")

long_covid_averted_by_age_longer$Age_group <- rep(c("0 to 17", "18 to 49", "50 to 64", "65 to 79", "80+"), length(vacc_scenarios[-1]))
long_covid_averted_by_age_longer_subset <- long_covid_averted_by_age_longer[grepl("1|2|3|6|9", long_covid_averted_by_age_longer$Strategy), ]

long_covid_averted_by_age_plot <- ggplot(data = long_covid_averted_by_age_longer_subset, 
                                         aes(x = Strategy, y = `Long Covid averted`, fill = Age_group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(brewer.pal(5, "Set2")))

# Panel plot

disease_averted_by_age_plot <- cases_averted_by_age_plot + labs(y = "Proportion", title = "(a) Cases averted", fill = "Age group") +
  theme(legend.position = "bottom", legend.justification = "left", title = element_text(size = 13)) +
  long_covid_averted_by_age_plot +  labs(y = "Proportion", title = "(b) Long Covid person-years averted") +
  theme(legend.position = "None", title = element_text(size = 13)) +
  hosp_averted_by_age_plot +  labs(y = "Proportion", title = "(c) Hospitalisations averted") +
  theme(legend.position = "None", title = element_text(size = 13)) +
  deaths_averted_by_age_plot +  labs(y = "Proportion", title = "(d) Deaths averted") +
  theme(legend.position = "None", title = element_text(size = 13)) +
  cases_averted_by_age_plot + labs(y = "Total number", title = "(e) Cases averted", fill = "Age group") +
  theme(legend.position = "bottom", legend.justification = "left", title = element_text(size = 13)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  long_covid_averted_by_age_plot +  labs(y = "Total number", title = "(f) Long Covid person-years averted") +
  theme(legend.position = "None", title = element_text(size = 13)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  hosp_averted_by_age_plot +  labs(y = "Total number", title = "(g) Hospitalisations averted") +
  theme(legend.position = "None", title = element_text(size = 13)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  deaths_averted_by_age_plot +  labs(y = "Total number", title = "(h) Deaths averted") +
  theme(legend.position = "None", title = element_text(size = 13)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) +
  plot_layout(nrow = 4, axes = "collect_y")

disease_averted_by_age_plot

ggsave(path = "Results/", filename = paste0(as.character(Sys.Date()), " Disease averted by age plot.jpeg"), 
       plot = disease_averted_by_age_plot, width = 10, height = 15)


