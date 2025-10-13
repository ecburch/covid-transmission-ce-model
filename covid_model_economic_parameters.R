
# Setting up CE model parameters

# Long Covid parameters

# Age stratified prevalence of LC amongst recoveries from the previous 2 years

long_covid_risk <- c(rep(0.0573, 25), rep(0.1045, 10), rep(0.1924, 10), rep(0.2542, 10), rep(0.2442, 15), rep(0.1430, 10),
                     rep(0.0544, 21))

# Reduction in LC risk due to vaccination

red_long_covid <- 0.59

# Waned reduction in LC risk

wan_long_covid <- 0.59

# ICU rate

icu_rate <- c(rep(0.083, 18), rep(0.128, 32), rep(0.185, 15), rep(0.191, 10), rep(0.098, 10), rep(0.017, 16))

# Health utility loss for each disease state

inf_symp_utility_loss <- c(rep(0.11, 16), rep(0.20, 85))
hosp_utility_loss <- 0.43
icu_utility_loss <- 0.55
lc_utility_loss <- 0.069

# Health utility values for each state

state_utilities_healthy <- c(rep(0.898, 16), rep(0.898, 2), rep(0.893, 2), rep(0.8765, 5), rep(0.882, 5), rep(0.892, 5), 
                             rep(0.8585, 5), rep(0.859, 5), rep(0.814, 5), rep(0.817, 5), rep(0.8, 5), rep(0.7895, 5), 
                             rep(0.786, 5), rep(0.7925, 5), rep(0.759, 5), rep(0.7385, 5), rep(0.6965, 5), rep(0.661, 11))
state_utilities_inf_symptomatic <- state_utilities_healthy - inf_symp_utility_loss
state_utilities_hospitalised <- state_utilities_healthy - (((1 - icu_rate) * hosp_utility_loss) + (icu_rate * icu_utility_loss))
state_utilities_lc <- state_utilities_healthy - lc_utility_loss
state_utilities_dead <- 0

# Divide all state utilities by 365 to get yearly health utilities instead of daily

state_utilities_healthy <- state_utilities_healthy / 365
state_utilities_inf_symptomatic <- state_utilities_inf_symptomatic / 365
state_utilities_hospitalised <- state_utilities_hospitalised / 365
state_utilities_lc <- state_utilities_lc / 365
state_utilities_dead <- state_utilities_dead / 365

# Transition costs for each state

state_costs_healthy <- 0
state_costs_vaccinated <- 25
state_costs_inf_symptomatic <- c(rep(29.63, 5), rep(28.15, 7), rep(37.26, 6), rep(38.29, 27), rep(39.02, 10), rep(39.79, 10), 
                                 rep(41.05, 10), rep(42.83, 10), rep(45.20, 16))
state_costs_hospitalised <- c(rep(4190.37, 5), rep(4594.37, 7), rep(10424.62, 6), rep(11283.58, 27), rep(13367.89, 10), 
                              rep(15631.64, 10), rep(17026.23, 10), rep(13595.86, 10), rep(11419.01, 16))
state_costs_lc <- 0.73
state_costs_dead <- 0

# Combining costs and utilities into lists

state_utilities <- list (state_utilities_healthy, state_utilities_inf_symptomatic, state_utilities_hospitalised, 
                         state_utilities_lc, state_utilities_dead)

state_costs <- list (state_costs_healthy, state_costs_vaccinated, state_costs_inf_symptomatic, state_costs_hospitalised, 
                     state_costs_lc, state_costs_dead)

names(state_utilities) <- c("Healthy", "Infectious symptomatic", "Hospitalised", "Long Covid", "Dead")

names(state_costs) <- c("Healthy", "Vaccinated", "Infectious symptomatic", "Hospitalised", "Long Covid", "Dead")


