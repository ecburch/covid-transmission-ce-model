
# BEFORE RUNNING CODE, CREATE A FOLDER IN WORKING DIRECTORY NAMED "RESULTS" AND A SUBFOLDER NAMED "SENSITIVITY ANALYSIS"


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
options(scipen=999)

# LOAD THE OUTPUT FROM RUNNING THE MODEL FOR 100 YEARS WITH 500 PARAMETER SETS. I.E. NEED THE "results",
#  "df_yearly_inf_hosp_deaths_lc", and "df_hosp_deaths_weekly_summaries" dataframes.

load("Model run 1.RData")

###############################################################################

# Testing higher vaccine cost (£31 per vaccine dose + £10 admin fee)

state_costs_vaccinated <- 41

# Combining costs and utilities into lists

state_utilities <- list (state_utilities_healthy, state_utilities_inf_symptomatic, state_utilities_hospitalised, 
                         state_utilities_lc, state_utilities_dead)

state_costs <- list (state_costs_healthy, state_costs_vaccinated, state_costs_inf_symptomatic, state_costs_hospitalised, 
                     state_costs_lc, state_costs_dead)

names(state_utilities) <- c("Healthy", "Infectious symptomatic", "Hospitalised", "Long Covid", "Dead")

names(state_costs) <- c("Healthy", "Vaccinated", "Infectious symptomatic", "Hospitalised", "Long Covid", "Dead")


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

life_years_remaining <- read_excel("Data/Life expectancy (England).xlsx", sheet="Cohort males and females")

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

for (j in 1:length(vacc_scenarios)) { 
  
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
    
    total_lc_person_years[i, j] <- sum(n_lc_recovered_total)/365
    
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
    
    # Calculating costs with no discounting
    
    total_vacc_cost_no_discounting[i, j] <- sum((df_results$inc_vacc_yearly[[i]]) * state_costs$Vaccinated)
    total_inf_cost_no_discounting[i, j] <- sum((df_results$inc_symp_inf_yearly[[i]]) * rep(state_costs$`Infectious symptomatic`,
                                                                                           each = time_horizon))
    total_hosp_cost_no_discounting[i, j] <- sum((df_results$inc_hosp_yearly[[i]]) * rep(state_costs$Hospitalised, each = time_horizon))
    total_lc_cost_no_discounting[i, j] <- sum(n_lc_recovered_total * state_costs$`Long Covid`)
    
  }
  
}

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

# Check for negative values in results

if(is.null(fixr::check_for_negative_values(c(as.numeric(unlist(results))))) == FALSE) {
  stop("Negative values in results")
}

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

cef_data_dominant <- cef_data_means[cef_data_means$dominant == "Yes", ]

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

# ggsave(path = "Results/Sensitivity analysis", 
#        filename = paste0(as.character(Sys.Date()), " Plot 2. Cost-effectiveness frontier (£41 vaccine).jpeg"),
#        plot = cef_plot, width = 8, height = 6)

# Check for negative incremental qalys

if(is.null(fixr::check_for_negative_values(-BCEA$delta_e)) == FALSE) {
  stop("Negative incremental qalys")
}

# Calculate 95% confidence intervals for the ICER

icer_simulations <- BCEA$delta_c / BCEA$delta_e         # delta_c are the incremental costs and delta_e incremental qalys for each simulation
icer_95_ci <- apply(icer_simulations, 2, quantile, probs = c(0.025, 0.975))

# Calculating ICER 95% CI

icer_mean <- c(0, round(as.numeric(BCEA$ICER)))
icer_95_ci_lower <- c(0, round(icer_95_ci[1, ]))
icer_95_ci_upper <- c(0, round(icer_95_ci[2, ]))

# Calculate the incremental net monetary benefit, with a willingness to pay threshold of £20,000 per QALY
# If the INMB is positive then the intervention is cost-effective

inmb <- (-BCEA$delta_e) * 20000 - (-BCEA$delta_c)
inmb <- cbind(rep(0, ncol(par_fit)), inmb)
inmb_mean <- colMeans(inmb)
inmb_95_ci <- apply(inmb, 2, quantile, probs = c(0.025, 0.975))
inmb_95_ci_lower <- round(inmb_95_ci[1, ])
inmb_95_ci_upper <- round(inmb_95_ci[2, ])

# Calculating total ICU admissions

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

# ggsave(path = "Results/Sensitivity analysis", 
#        filename = paste0(as.character(Sys.Date()), " Plot 3. Cost-effectiveness plane (£41 vaccine).jpeg"),
#        plot = cost_effectiveness_plane, height = 12, width = 8)

# Plot the cost-effectiveness acceptability curve: probability that vaccination is CE for various WTP thresholds

wtp_thresholds <- seq(0, 40000, by = 100)
prob_cost_effective <- matrix(1, nrow = length(wtp_thresholds), ncol = ncol(icer_simulations))

for (i in 1:length(wtp_thresholds)) {
  prob_cost_effective[i, ] <- colSums(icer_simulations < wtp_thresholds[i]) / nrow(icer_simulations)
}

prob_cost_effective_20k <- round(c(0, prob_cost_effective[201, ]) * 100)
prob_cost_effective_30k <- round(c(0, prob_cost_effective[301, ]) * 100)

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

# ggsave(path = "Results/Sensitivity analysis", 
#        filename = paste0(as.character(Sys.Date()), " Plot 4. Cost effectiveness acceptability curve (£41 vaccine).jpeg"), 
#        plot = cost_effectiveness_acceptability_curve, width = 8, height = 5)

# Calculating means and confidence intervals for all results

results_mean <- results_95_ci <- list()
results_merged <- rep(NA, length(vacc_scenarios))
for (i in 1:length(results)) {
  results_mean[[i]] <- as.data.frame(apply(results[[i]], 2, mean))
  names(results_mean[[i]]) <- paste0(names(results)[i], " mean")
  results_95_ci[[i]] <- apply(results[[i]], 2, quantile, probs = c(0.025, 0.975))
  results_merged <- rbind(results_merged, t(results_mean[[i]]), results_95_ci[[i]])
}
results_merged <- round(results_merged[-1, ])

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

incremental_results_merged <- round(incremental_results_merged[-1, ])
names(results_mean) <- names(results_95_ci) <- names(results)

# Creating results table

colnames(results_merged) <- colnames(incremental_results_merged) <- vacc_scenarios

full_epi_results <- rbind.data.frame(results_merged[22:42, ], incremental_results_merged[22:36, ])

full_ce_results <- rbind.data.frame("ICER" = icer_mean, "ICER lower bound" = icer_95_ci_lower, "ICER upper bound" = icer_95_ci_upper, 
                                    "INMB" = inmb_mean, "INMB lower bound" = inmb_95_ci_lower, "INMB upper bound" = inmb_95_ci_upper,
                                    rbind.data.frame(results_merged[1:21, ], incremental_results_merged[1:21, ]),
                                    "Probability of cost effectiveness, £20k WTP" = prob_cost_effective_20k, 
                                    "Probability of cost effectiveness, £30k WTP" = prob_cost_effective_30k, 
                                    "Probability cost saving" = prob_cost_saving)

colnames(full_epi_results) <- colnames(full_ce_results) <- c("BASELINE, comparator", vacc_scenarios[-1])
model_output <- list("Epi results" = full_epi_results, "CE results" = full_ce_results, Parameters = par_fit)
write.xlsx((model_output), paste0("Results/Sensitivity analysis/", as.character(Sys.Date()), " Results table, £41 vaccine.xlsx"), 
           rowNames = TRUE)


###############################################################################

# Testing 1.5% discounting

# Clear environment then reload model simulation results

rm(list = ls())

load("Model run 1.RData")

# Setting discount rate

set_costs_discount_rate <- 0.015           # 1.5% discounting
set_qalys_discount_rate <- 0.015           # 1.5% discounting

future_years <- 0:(time_horizon - 1)
costs_discount_rate <- 1 / ((1 + set_costs_discount_rate) ^ future_years)
qalys_discount_rate <- 1 / ((1 + set_qalys_discount_rate) ^ future_years)

# Vaccine uptake rate

vacc_uptake_rate <- c(rep(0, 5), rep(6.4, 7), rep(31.7, 4), rep(31.7, 2), rep(31.7, 2), rep(40.3, 5), rep(40.9, 5), 
                      rep(44.0, 5), rep(49.5, 5), rep(57.0, 5), rep(64.8, 5), rep(64.8, 5), rep(64.8, 5), 
                      rep(64.8, 5), rep(64.8, 5), rep(70.1, 5), rep(75.5, 5), rep(75.7, 21)) / 100

# Extracting life years remaining

life_years_remaining <- read_excel("Data/Life expectancy (England).xlsx", sheet="Cohort males and females")

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

for (j in 1:length(vacc_scenarios)) { 
  
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
    
    total_lc_person_years[i, j] <- sum(n_lc_recovered_total)/365
    
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
    
    # Calculating costs with no discounting
    
    total_vacc_cost_no_discounting[i, j] <- sum((df_results$inc_vacc_yearly[[i]]) * state_costs$Vaccinated)
    total_inf_cost_no_discounting[i, j] <- sum((df_results$inc_symp_inf_yearly[[i]]) * rep(state_costs$`Infectious symptomatic`,
                                                                                           each = time_horizon))
    total_hosp_cost_no_discounting[i, j] <- sum((df_results$inc_hosp_yearly[[i]]) * rep(state_costs$Hospitalised, each = time_horizon))
    total_lc_cost_no_discounting[i, j] <- sum(n_lc_recovered_total * state_costs$`Long Covid`)
    
  }
  
}

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

# Check for negative values in results

if(is.null(fixr::check_for_negative_values(c(as.numeric(unlist(results))))) == FALSE) {
  stop("Negative values in results")
}

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

cef_data_dominant <- cef_data_means[cef_data_means$dominant == "Yes", ]

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 2. Cost-effectiveness frontier (1.5% discounting).jpeg"),
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
icer_95_ci_lower <- c(0, round(icer_95_ci[1, ]))
icer_95_ci_upper <- c(0, round(icer_95_ci[2, ]))

# Calculate the incremental net monetary benefit, with a willingness to pay threshold of £20,000 per QALY
# If the INMB is positive then the intervention is cost-effective

inmb <- (-BCEA$delta_e) * 20000 - (-BCEA$delta_c)
inmb <- cbind(rep(0, ncol(par_fit)), inmb)
inmb_mean <- colMeans(inmb)
inmb_95_ci <- apply(inmb, 2, quantile, probs = c(0.025, 0.975))
inmb_95_ci_lower <- round(inmb_95_ci[1, ])
inmb_95_ci_upper <- round(inmb_95_ci[2, ])

# Calculating total ICU admissions

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 3. Cost-effectiveness plane (1.5% discounting).jpeg"),
       plot = cost_effectiveness_plane, height = 12, width = 8)

# Plot the cost-effectiveness acceptability curve: probability that vaccination is CE for various WTP thresholds

wtp_thresholds <- seq(0, 40000, by = 100)
prob_cost_effective <- matrix(1, nrow = length(wtp_thresholds), ncol = ncol(icer_simulations))

for (i in 1:length(wtp_thresholds)) {
  prob_cost_effective[i, ] <- colSums(icer_simulations < wtp_thresholds[i]) / nrow(icer_simulations)
}

prob_cost_effective_20k <- round(c(0, prob_cost_effective[201, ]) * 100)
prob_cost_effective_30k <- round(c(0, prob_cost_effective[301, ]) * 100)

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 4. Cost effectiveness acceptability curve (1.5% discounting).jpeg"), 
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
results_merged <- round(results_merged[-1, ])

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

incremental_results_merged <- round(incremental_results_merged[-1, ])
names(results_mean) <- names(results_95_ci) <- names(results)

# Creating results table

colnames(results_merged) <- colnames(incremental_results_merged) <- vacc_scenarios

full_epi_results <- rbind.data.frame(results_merged[22:42, ], incremental_results_merged[22:36, ])

full_ce_results <- rbind.data.frame("ICER" = icer_mean, "ICER lower bound" = icer_95_ci_lower, "ICER upper bound" = icer_95_ci_upper, 
                                    "INMB" = inmb_mean, "INMB lower bound" = inmb_95_ci_lower, "INMB upper bound" = inmb_95_ci_upper,
                                    rbind.data.frame(results_merged[1:21, ], incremental_results_merged[1:21, ]),
                                    "Probability of cost effectiveness, £20k WTP" = prob_cost_effective_20k, 
                                    "Probability of cost effectiveness, £30k WTP" = prob_cost_effective_30k, 
                                    "Probability cost saving" = prob_cost_saving)

colnames(full_epi_results) <- colnames(full_ce_results) <- c("BASELINE, comparator", vacc_scenarios[-1])
model_output <- list("Epi results" = full_epi_results, "CE results" = full_ce_results, Parameters = par_fit)
write.xlsx((model_output), paste0("Results/Sensitivity analysis/", as.character(Sys.Date()), " Results table, 1.5% discounting.xlsx"), 
           rowNames = TRUE)




###############################################################################

# Testing 10 year horizon

# Clear environment then reload model simulation results

rm(list = ls())

load("Model run 1.RData")

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

life_years_remaining <- read_excel("Data/Life expectancy (England).xlsx", sheet="Cohort males and females")

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

for (j in 1:length(vacc_scenarios)) { 
  
  for (i in 1:ncol(par_fit)) {
    
    # Set dataframe for vaccination scenario
    
    df_results <- df_yearly_inf_hosp_deaths_lc[[j]]
    
    # Calculating infections, symptomatic infections, hospitalisations, ICU admissions, deaths and vaccinations
    
    total_infections[i, j] <- sum(df_results$inc_inf_yearly[[i]][1:10, ])
    total_hosp[i, j] <- sum(df_results$inc_hosp_yearly[[i]][1:10, ])
    total_icu[i, j] <- sum(icu_rate * colSums(df_results$inc_hosp_yearly[[i]][1:10, ]))
    total_deaths[i, j] <- sum(df_results$inc_deaths_yearly[[i]][1:10, ])
    total_vacc[i, j] <- sum(df_results$inc_vacc_yearly[[i]][1:10, ])
    total_symp_inf[i, j] <- sum(df_results$inc_symp_inf_yearly[[i]][1:10, ])
    
    # Calculating life years lost
    
    total_life_years_lost[i, j] <- sum(df_results$inc_deaths_yearly[[i]][1:10, ] * life_years_remaining[1:10, ])
    
    # Calculating costs
    
    total_vacc_cost[i, j] <- sum(costs_discount_rate[1:10] * ((df_results$inc_vacc_yearly[[i]][1:10, ]) * state_costs$Vaccinated))
    total_inf_cost[i, j] <- sum(costs_discount_rate[1:10] * ((df_results$inc_symp_inf_yearly[[i]][1:10, ]) * 
                                                               rep(state_costs$`Infectious symptomatic`, each = 10)))
    total_hosp_cost[i, j] <- sum(costs_discount_rate[1:10] * ((df_results$inc_hosp_yearly[[i]][1:10, ]) * 
                                                                rep(state_costs$Hospitalised, each = 10)))
    
    # Calculating long Covid costs
    
    n_lc_recovered <- df_results$n_recovered_past_year_yearly[[i]] * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_vaccinated <- df_results$n_recovered_vaccinated_past_year_yearly[[i]] * red_long_covid * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_waned <- df_results$n_recovered_waned_past_year_yearly[[i]] * wan_long_covid * rep(long_covid_risk, each = time_horizon)
    n_lc_recovered_total <- n_lc_recovered + n_lc_recovered_vaccinated + n_lc_recovered_waned
    total_lc_cost[i, j] <- sum(costs_discount_rate[1:10] * (n_lc_recovered_total[1:10, ] * state_costs$`Long Covid`))
    
    # Calculating LC person-years
    
    total_lc_person_years[i, j] <- sum(n_lc_recovered_total[1:10, ])/365
    
    # Calculating QALYs
    
    total_qalys[i, j] <- sum(qalys_discount_rate[1:10] * ((df_results$n_healthy_yearly[[i]][1:10, ] -
                                                             n_lc_recovered_total[1:10, ]) * rep(state_utilities$Healthy,
                                                                                                 each = 10) +
                                                            df_results$n_symptomatic_yearly[[i]][1:10, ] *
                                                            rep(state_utilities$`Infectious symptomatic`, each = 10) +
                                                            df_results$n_hospitalised_yearly[[i]][1:10, ] *
                                                            rep(state_utilities$Hospitalised, each = 10) +
                                                            n_lc_recovered_total[1:10, ] * rep(state_utilities$`Long Covid`,
                                                                                               each = 10)))
    
    # Calculating costs with no discounting
    
    total_vacc_cost_no_discounting[i, j] <- sum((df_results$inc_vacc_yearly[[i]]) * state_costs$Vaccinated)
    total_inf_cost_no_discounting[i, j] <- sum((df_results$inc_symp_inf_yearly[[i]]) * rep(state_costs$`Infectious symptomatic`,
                                                                                           each = time_horizon))
    total_hosp_cost_no_discounting[i, j] <- sum((df_results$inc_hosp_yearly[[i]]) * rep(state_costs$Hospitalised, each = time_horizon))
    total_lc_cost_no_discounting[i, j] <- sum(n_lc_recovered_total * state_costs$`Long Covid`)
    
  }
  
}

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

# Check for negative values in results

if(is.null(fixr::check_for_negative_values(c(as.numeric(unlist(results))))) == FALSE) {
  stop("Negative values in results")
}

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

cef_data_dominant <- cef_data_means[cef_data_means$dominant == "Yes", ]

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 2. Cost-effectiveness frontier (10 year time horizon).jpeg"),
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
icer_95_ci_lower <- c(0, round(icer_95_ci[1, ]))
icer_95_ci_upper <- c(0, round(icer_95_ci[2, ]))

# Calculate the incremental net monetary benefit, with a willingness to pay threshold of £20,000 per QALY
# If the INMB is positive then the intervention is cost-effective

inmb <- (-BCEA$delta_e) * 20000 - (-BCEA$delta_c)
inmb <- cbind(rep(0, ncol(par_fit)), inmb)
inmb_mean <- colMeans(inmb)
inmb_95_ci <- apply(inmb, 2, quantile, probs = c(0.025, 0.975))
inmb_95_ci_lower <- round(inmb_95_ci[1, ])
inmb_95_ci_upper <- round(inmb_95_ci[2, ])

# Calculating total ICU admissions

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 3. Cost-effectiveness plane (10 year time horizon).jpeg"),
       plot = cost_effectiveness_plane, height = 12, width = 8)

# Plot the cost-effectiveness acceptability curve: probability that vaccination is CE for various WTP thresholds

wtp_thresholds <- seq(0, 40000, by = 100)
prob_cost_effective <- matrix(1, nrow = length(wtp_thresholds), ncol = ncol(icer_simulations))

for (i in 1:length(wtp_thresholds)) {
  prob_cost_effective[i, ] <- colSums(icer_simulations < wtp_thresholds[i]) / nrow(icer_simulations)
}

prob_cost_effective_20k <- round(c(0, prob_cost_effective[201, ]) * 100)
prob_cost_effective_30k <- round(c(0, prob_cost_effective[301, ]) * 100)

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 4. Cost effectiveness acceptability curve (10 year time horizon).jpeg"), 
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
results_merged <- round(results_merged[-1, ])

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

incremental_results_merged <- round(incremental_results_merged[-1, ])
names(results_mean) <- names(results_95_ci) <- names(results)

# Creating results table

colnames(results_merged) <- colnames(incremental_results_merged) <- vacc_scenarios

full_epi_results <- rbind.data.frame(results_merged[22:42, ], incremental_results_merged[22:36, ])

full_ce_results <- rbind.data.frame("ICER" = icer_mean, "ICER lower bound" = icer_95_ci_lower, "ICER upper bound" = icer_95_ci_upper, 
                                    "INMB" = inmb_mean, "INMB lower bound" = inmb_95_ci_lower, "INMB upper bound" = inmb_95_ci_upper,
                                    rbind.data.frame(results_merged[1:21, ], incremental_results_merged[1:21, ]),
                                    "Probability of cost effectiveness, £20k WTP" = prob_cost_effective_20k, 
                                    "Probability of cost effectiveness, £30k WTP" = prob_cost_effective_30k, 
                                    "Probability cost saving" = prob_cost_saving)

colnames(full_epi_results) <- colnames(full_ce_results) <- c("BASELINE, comparator", vacc_scenarios[-1])
model_output <- list("Epi results" = full_epi_results, "CE results" = full_ce_results, Parameters = par_fit)
write.xlsx((model_output), paste0("Results/Sensitivity analysis/", as.character(Sys.Date()), " Results table, 10 year time horizon.xlsx"), 
           rowNames = TRUE)





###############################################################################

# Testing no effects from long Covid

# Clear environment then reload model simulation results

rm(list = ls())

load("Model run 1.RData")

# Change long Covid cost and QALY loss

state_utilities_lc <- state_utilities_healthy
state_costs_lc <- 0

# Combining costs and utilities into lists

state_utilities <- list (state_utilities_healthy, state_utilities_inf_symptomatic, state_utilities_hospitalised, 
                         state_utilities_lc, state_utilities_dead)

state_costs <- list (state_costs_healthy, state_costs_vaccinated, state_costs_inf_symptomatic, state_costs_hospitalised, 
                     state_costs_lc, state_costs_dead)

names(state_utilities) <- c("Healthy", "Infectious symptomatic", "Hospitalised", "Long Covid", "Dead")

names(state_costs) <- c("Healthy", "Vaccinated", "Infectious symptomatic", "Hospitalised", "Long Covid", "Dead")

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

life_years_remaining <- read_excel("Data/Life expectancy (England).xlsx", sheet="Cohort males and females")

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

for (j in 1:length(vacc_scenarios)) { 
  
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
    
    total_lc_person_years[i, j] <- sum(n_lc_recovered_total)/365
    
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
    
    # Calculating costs with no discounting
    
    total_vacc_cost_no_discounting[i, j] <- sum((df_results$inc_vacc_yearly[[i]]) * state_costs$Vaccinated)
    total_inf_cost_no_discounting[i, j] <- sum((df_results$inc_symp_inf_yearly[[i]]) * rep(state_costs$`Infectious symptomatic`,
                                                                                           each = time_horizon))
    total_hosp_cost_no_discounting[i, j] <- sum((df_results$inc_hosp_yearly[[i]]) * rep(state_costs$Hospitalised, each = time_horizon))
    total_lc_cost_no_discounting[i, j] <- sum(n_lc_recovered_total * state_costs$`Long Covid`)
    
  }
  
}

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

# Check for negative values in results

if(is.null(fixr::check_for_negative_values(c(as.numeric(unlist(results))))) == FALSE) {
  stop("Negative values in results")
}

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

cef_data_dominant <- cef_data_means[cef_data_means$dominant == "Yes", ]

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 2. Cost-effectiveness frontier (No LC effects).jpeg"),
       plot = cef_plot, width = 8, height = 6)

cef_plot_no_lc <- cef_plot
saveRDS(cef_plot_no_lc, file = paste0("cef_plot_no_lc", ".rds"))

# Check for negative incremental qalys

if(is.null(fixr::check_for_negative_values(-BCEA$delta_e)) == FALSE) {
  stop("Negative incremental qalys")
}

# Calculate 95% confidence intervals for the ICER

icer_simulations <- BCEA$delta_c / BCEA$delta_e         # delta_c are the incremental costs and delta_e incremental qalys for each simulation
icer_95_ci <- apply(icer_simulations, 2, quantile, probs = c(0.025, 0.975))

# Calculating ICER 95% CI

icer_mean <- c(0, round(as.numeric(BCEA$ICER)))
icer_95_ci_lower <- c(0, round(icer_95_ci[1, ]))
icer_95_ci_upper <- c(0, round(icer_95_ci[2, ]))

# Calculate the incremental net monetary benefit, with a willingness to pay threshold of £20,000 per QALY
# If the INMB is positive then the intervention is cost-effective

inmb <- (-BCEA$delta_e) * 20000 - (-BCEA$delta_c)
inmb <- cbind(rep(0, ncol(par_fit)), inmb)
inmb_mean <- colMeans(inmb)
inmb_95_ci <- apply(inmb, 2, quantile, probs = c(0.025, 0.975))
inmb_95_ci_lower <- round(inmb_95_ci[1, ])
inmb_95_ci_upper <- round(inmb_95_ci[2, ])

# Calculating total ICU admissions

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 3. Cost-effectiveness plane (No LC effects).jpeg"),
       plot = cost_effectiveness_plane, height = 12, width = 8)

# Plot the cost-effectiveness acceptability curve: probability that vaccination is CE for various WTP thresholds

wtp_thresholds <- seq(0, 40000, by = 100)
prob_cost_effective <- matrix(1, nrow = length(wtp_thresholds), ncol = ncol(icer_simulations))

for (i in 1:length(wtp_thresholds)) {
  prob_cost_effective[i, ] <- colSums(icer_simulations < wtp_thresholds[i]) / nrow(icer_simulations)
}

prob_cost_effective_20k <- round(c(0, prob_cost_effective[201, ]) * 100)
prob_cost_effective_30k <- round(c(0, prob_cost_effective[301, ]) * 100)

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 4. Cost effectiveness acceptability curve (No LC effects).jpeg"), 
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
results_merged <- round(results_merged[-1, ])

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

incremental_results_merged <- round(incremental_results_merged[-1, ])
names(results_mean) <- names(results_95_ci) <- names(results)

# Creating results table

colnames(results_merged) <- colnames(incremental_results_merged) <- vacc_scenarios

full_epi_results <- rbind.data.frame(results_merged[22:42, ], incremental_results_merged[22:36, ])

full_ce_results <- rbind.data.frame("ICER" = icer_mean, "ICER lower bound" = icer_95_ci_lower, "ICER upper bound" = icer_95_ci_upper, 
                                    "INMB" = inmb_mean, "INMB lower bound" = inmb_95_ci_lower, "INMB upper bound" = inmb_95_ci_upper,
                                    rbind.data.frame(results_merged[1:21, ], incremental_results_merged[1:21, ]),
                                    "Probability of cost effectiveness, £20k WTP" = prob_cost_effective_20k, 
                                    "Probability of cost effectiveness, £30k WTP" = prob_cost_effective_30k, 
                                    "Probability cost saving" = prob_cost_saving)

colnames(full_epi_results) <- colnames(full_ce_results) <- c("BASELINE, comparator", vacc_scenarios[-1])
model_output <- list("Epi results" = full_epi_results, "CE results" = full_ce_results, Parameters = par_fit)
write.xlsx((model_output), paste0("Results/Sensitivity analysis/", as.character(Sys.Date()), " Results table, No LC effects.xlsx"), 
           rowNames = TRUE)




###############################################################################

# Testing greater cost and qaly loss from long Covid

# Clear environment then reload model simulation results

rm(list = ls())

load("Model run 1.RData")

# Load cost-effectiveness frontier plot from previous analysis with no effects from long Covid

cef_plot_no_lc <- readRDS("cef_plot_no_lc.rds")

# Change long Covid cost and QALY loss

lc_utility_loss_upper_bound <- 0.13

# Health utility values for each state

state_utilities_lc <- state_utilities_healthy - lc_utility_loss_upper_bound

# Divide all state utilities by 365 to get yearly health utilities instead of daily

state_utilities_lc <- state_utilities_lc / 365

# Upper bound on cost of long Covid

state_costs_lc <- 2.84

# Combining costs and utilities into lists

state_utilities <- list (state_utilities_healthy, state_utilities_inf_symptomatic, state_utilities_hospitalised, 
                         state_utilities_lc, state_utilities_dead)

state_costs <- list (state_costs_healthy, state_costs_vaccinated, state_costs_inf_symptomatic, state_costs_hospitalised, 
                     state_costs_lc, state_costs_dead)

names(state_utilities) <- c("Healthy", "Infectious symptomatic", "Hospitalised", "Long Covid", "Dead")

names(state_costs) <- c("Healthy", "Vaccinated", "Infectious symptomatic", "Hospitalised", "Long Covid", "Dead")

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

life_years_remaining <- read_excel("Data/Life expectancy (England).xlsx", sheet="Cohort males and females")

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

for (j in 1:length(vacc_scenarios)) { 
  
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
    
    total_lc_person_years[i, j] <- sum(n_lc_recovered_total)/365
    
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
    
    # Calculating costs with no discounting
    
    total_vacc_cost_no_discounting[i, j] <- sum((df_results$inc_vacc_yearly[[i]]) * state_costs$Vaccinated)
    total_inf_cost_no_discounting[i, j] <- sum((df_results$inc_symp_inf_yearly[[i]]) * rep(state_costs$`Infectious symptomatic`,
                                                                                           each = time_horizon))
    total_hosp_cost_no_discounting[i, j] <- sum((df_results$inc_hosp_yearly[[i]]) * rep(state_costs$Hospitalised, each = time_horizon))
    total_lc_cost_no_discounting[i, j] <- sum(n_lc_recovered_total * state_costs$`Long Covid`)
    
  }
  
}

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

# Check for negative values in results

if(is.null(fixr::check_for_negative_values(c(as.numeric(unlist(results))))) == FALSE) {
  stop("Negative values in results")
}

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

cef_data_dominant <- cef_data_means[cef_data_means$dominant == "Yes", ]

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 2. Cost-effectiveness frontier (More severe LC).jpeg"),
       plot = cef_plot, width = 8, height = 6)

cef_plot_severe_lc <- cef_plot

# Check for negative incremental qalys

if(is.null(fixr::check_for_negative_values(-BCEA$delta_e)) == FALSE) {
  stop("Negative incremental qalys")
}

# Calculate 95% confidence intervals for the ICER

icer_simulations <- BCEA$delta_c / BCEA$delta_e         # delta_c are the incremental costs and delta_e incremental qalys for each simulation
icer_95_ci <- apply(icer_simulations, 2, quantile, probs = c(0.025, 0.975))

# Calculating ICER 95% CI

icer_mean <- c(0, round(as.numeric(BCEA$ICER)))
icer_95_ci_lower <- c(0, round(icer_95_ci[1, ]))
icer_95_ci_upper <- c(0, round(icer_95_ci[2, ]))

# Calculate the incremental net monetary benefit, with a willingness to pay threshold of £20,000 per QALY
# If the INMB is positive then the intervention is cost-effective

inmb <- (-BCEA$delta_e) * 20000 - (-BCEA$delta_c)
inmb <- cbind(rep(0, ncol(par_fit)), inmb)
inmb_mean <- colMeans(inmb)
inmb_95_ci <- apply(inmb, 2, quantile, probs = c(0.025, 0.975))
inmb_95_ci_lower <- round(inmb_95_ci[1, ])
inmb_95_ci_upper <- round(inmb_95_ci[2, ])

# Calculating total ICU admissions

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 3. Cost-effectiveness plane (More severe LC).jpeg"),
       plot = cost_effectiveness_plane, height = 12, width = 8)

# Plot the cost-effectiveness acceptability curve: probability that vaccination is CE for various WTP thresholds

wtp_thresholds <- seq(0, 40000, by = 100)
prob_cost_effective <- matrix(1, nrow = length(wtp_thresholds), ncol = ncol(icer_simulations))

for (i in 1:length(wtp_thresholds)) {
  prob_cost_effective[i, ] <- colSums(icer_simulations < wtp_thresholds[i]) / nrow(icer_simulations)
}

prob_cost_effective_20k <- round(c(0, prob_cost_effective[201, ]) * 100)
prob_cost_effective_30k <- round(c(0, prob_cost_effective[301, ]) * 100)

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

ggsave(path = "Results/Sensitivity analysis", 
       filename = paste0(as.character(Sys.Date()), " Plot 4. Cost effectiveness acceptability curve (More severe LC).jpeg"), 
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
results_merged <- round(results_merged[-1, ])

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

incremental_results_merged <- round(incremental_results_merged[-1, ])
names(results_mean) <- names(results_95_ci) <- names(results)

# Creating results table

colnames(results_merged) <- colnames(incremental_results_merged) <- vacc_scenarios

full_epi_results <- rbind.data.frame(results_merged[22:42, ], incremental_results_merged[22:36, ])

full_ce_results <- rbind.data.frame("ICER" = icer_mean, "ICER lower bound" = icer_95_ci_lower, "ICER upper bound" = icer_95_ci_upper, 
                                    "INMB" = inmb_mean, "INMB lower bound" = inmb_95_ci_lower, "INMB upper bound" = inmb_95_ci_upper,
                                    rbind.data.frame(results_merged[1:21, ], incremental_results_merged[1:21, ]),
                                    "Probability of cost effectiveness, £20k WTP" = prob_cost_effective_20k, 
                                    "Probability of cost effectiveness, £30k WTP" = prob_cost_effective_30k, 
                                    "Probability cost saving" = prob_cost_saving)

colnames(full_epi_results) <- colnames(full_ce_results) <- c("BASELINE, comparator", vacc_scenarios[-1])
model_output <- list("Epi results" = full_epi_results, "CE results" = full_ce_results, Parameters = par_fit)
write.xlsx((model_output), paste0("Results/Sensitivity analysis/", as.character(Sys.Date()), " Results table, More severe LC.xlsx"), 
           rowNames = TRUE)


##############################################################################################

lc_sa_cef_plots <- cef_plot_no_lc + 
  theme(legend.position = "none") + 
  labs(title = "(a) No cost or QALY loss from long Covid") +
  
  cef_plot_severe_lc + 
  theme(legend.position = "none", axis.title.y = element_blank()) + labs(title = "(b) Increased cost and QALY loss from long Covid") + 
  plot_layout(nrow = 1, axis_titles = "collect")

print(lc_sa_cef_plots)

ggsave(path = "Results/Sensitivity analysis", filename = paste0(as.character(Sys.Date()),
                                                                " Plot ii. Cost-effectiveness frontier (long covid parameters).jpeg"),
       plot = lc_sa_cef_plots, width = 10, height = 5)
