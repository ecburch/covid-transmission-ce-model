
# SET WORKING DIRECTORY TO THE FOLDER CONTAINING MODEL CODE AND DATA
# CREATE A FOLDER IN THE WORKING DIRECTORY NAMED "Results" AND A SUBFOLDER WITHIN IT NAMED "Model fitting"



# Set the start run time, load packages and set seed

start_time <- Sys.time()

library(deSolve)
library(ggplot2)
library(socialmixr)
library(readxl)
library(EasyABC)
library(tibble)
library(janitor)
library(tidytable)
library(openxlsx)
library(matlib)
library(fixr)
library(utils)
library(stats)
library(forcats)
library(gridExtra)
library(cowplot)
library(patchwork)
library(GGally)
library(scales)

set.seed(1)
theme_set(theme_light(base_size = 15))
options(scipen = 999)

# Model EasyABC run function with fixed parameters set at the mean value.

covid_model_run_easy_ABC_set_seed <- function(par_fit) {
  
  # Create function "covid_model_odes" that computes the values of the model ODEs.
  
  covid_model_odes <- function(t, states, par){
    with(as.list(c(states, par)), {
      
      # Parameters
      
      ageing <- par$Ageing
      lat_period <- par$`Latent period`
      prop_symp <- par$`Proportion symptomatic`
      symp_period <- par$`Symptomatic period`
      rec_rate <- par$`Recovery rate (asymptomatic)`
      prop_hosp <- par$`Proportion hospitalised`
      hosp_period <- par$`Hospital length of stay`
      prop_hosp_mortalities <- par$`Proportion of hospital mortalities`
      vacc_rate <- par$`Vaccination rate`
      vacc_rate_shift <- par$`Vaccination rate shift`
      n_vacc_waves <- par$`Number of vaccination waves`         # No vaccination for under 12s, two vaccinations a year for 12-100 year olds
      vacc_uptake_rate <- par$`Vaccine uptake rate`
      wane_rate <- par$`Wane rate`                      # Immunity wanes after 270 days (nine months)
      red_I <- par$`Reduction in infectiousness`
      red_S <- par$`Reduction in symptomatic infections`
      red_H <- par$`Reduction in hospitalisations`
      red_HM <- par$`Reduction in hospital mortalities`
      wan_I <- par$`Waned reduction in infectiousness`
      wan_S <- par$`Waned reduction in symptomatic infections`
      wan_H <- par$`Waned reduction in hospitalisations`
      wan_HM <- par$`Waned reduction in hospital mortalities`
      reinf <- par$`Protection against reinfection`
      red_T <- par$`Reduction in transmission`
      wan_T <- par$`Waned reduction in transmission`
      tau <- par$`Reduction in asymptomatic infectiousness`
      beta <- par$`Transmission rate`
      num_peaks <- par$`Number of annual peaks in transmission`
      inf_fluctuation <- par$`Amount of fluctuation in infectiousness`
      inf_shift <- par$`Infection shift`
      n_contacts <- par$`Number of contacts per age group`
      total_pop <- par$`Total population size`
      n_age_groups <- par$`Number of age groups`
      all_cause_death_rate <- par$`All-cause death rate`
      birth_rate <- par$`Birth rate`
      
      # States
      
      S = states[1:n_age_groups]                                        # Susceptible
      E = states[(n_age_groups + 1):(2 * n_age_groups)]                     # Exposed
      I = states[(2 * n_age_groups + 1):(3 * n_age_groups)]                   # Infected, symptomatic
      A = states[(3 * n_age_groups + 1):(4 * n_age_groups)]                   # Infected, asymptomatic
      H = states[(4 * n_age_groups + 1):(5 * n_age_groups)]                   # Hospitalised
      R = states[(5 * n_age_groups + 1):(6 * n_age_groups)]                   # Recovered
      D = states[(6 * n_age_groups + 1):(7 * n_age_groups)]                   # Death
      V = states[(7 * n_age_groups + 1):(8 * n_age_groups)]                   # Vaccinated
      Ev = states[(8 * n_age_groups + 1):(9 * n_age_groups)]                  # Exposed-vaccinated
      Iv = states[(9 * n_age_groups + 1):(10 * n_age_groups)]                 # Infected, symptomatic-vaccinated
      Av = states[(10 * n_age_groups + 1):(11 * n_age_groups)]                # Infected, asymptomatic-vaccinated
      Hv = states[(11 * n_age_groups + 1):(12 * n_age_groups)]                # Hospitalised-vaccinated
      Rv = states[(12 * n_age_groups + 1):(13 * n_age_groups)]                # Recovered-vaccinated
      W = states[(13 * n_age_groups + 1):(14 * n_age_groups)]                 # Waned immunity
      Ew = states[(14 * n_age_groups + 1):(15 * n_age_groups)]                # Exposed-waned
      Iw = states[(15 * n_age_groups + 1):(16 * n_age_groups)]                # Infected, symptomatic-waned
      Aw = states[(16 * n_age_groups + 1):(17 * n_age_groups)]                # Infected, asymptomatic-waned
      Hw = states[(17 * n_age_groups + 1):(18 * n_age_groups)]                # Hospitalised-waned
      Rw = states[(18 * n_age_groups + 1):(19 * n_age_groups)]                # Recovered-waned
      II = states[(19 * n_age_groups + 1):(20 * n_age_groups)]                # Combined, cumulative infected states
      HH = states[(20 * n_age_groups + 1):(21 * n_age_groups)]                # Combined, cumulative hospitalised states
      VV = states[(21 * n_age_groups + 1):(22 * n_age_groups)]                # Combined vaccinated states
      VP = states[(22 * n_age_groups + 1):(23 * n_age_groups)]                # Combined vaccinated (primary programme) states
      VB = states[(23 * n_age_groups + 1):(24 * n_age_groups)]                # Combined vaccinated (booster programme) states
      IIS = states[(24 * n_age_groups + 1):(25 * n_age_groups)]               # Combined, cumulative symptomatic infected states
      RR = states[(25 * n_age_groups + 1):(26 * n_age_groups)]                # Cumulative symptomatic recovered unvaccinated state
      RRV = states[(26 * n_age_groups + 1):(27 * n_age_groups)]               # Cumulative symptomatic recovered vaccinated state
      RRW = states[(27 * n_age_groups + 1):(28 * n_age_groups)]               # Cumulative symptomatic recovered waned state
      
      # Specify total population
      
      N <- S + E + I + A + H + R + V + Ev + Iv + Av + Hv + Rv + W + Ew + Iw + Aw + Hw + Rw
      
      # All age compartments shifted up one for ageing
      
      older_S <- c(birth_rate, S)[1:length(S)]
      older_E <- c(0, E)[1:length(E)]
      older_I <- c(0, I)[1:length(I)]
      older_A<- c(0, A)[1:length(A)]
      older_H <- c(0, H)[1:length(H)]
      older_R <- c(0, R)[1:length(R)]
      older_V <- c(0, V)[1:length(V)]
      older_Ev <- c(0, Ev)[1:length(Ev)]
      older_Iv <- c(0, Iv)[1:length(Iv)]
      older_Av <- c(0, Av)[1:length(Av)]
      older_Hv <- c(0, Hv)[1:length(Hv)]
      older_Rv <- c(0, Rv)[1:length(Rv)]
      older_W <- c(0, W)[1:length(W)]
      older_Ew <- c(0, Ew)[1:length(Ew)]
      older_Iw <- c(0, Iw)[1:length(Iw)]
      older_Aw <- c(0, Aw)[1:length(Aw)]
      older_Hw <- c(0, Hw)[1:length(Hw)]
      older_Rw <- c(0, Rw)[1:length(Rw)]
      
      # Equation for frequency-dependent force of infection
      
      lambda <- (beta[t] + inf_fluctuation * (sin(inf_shift + num_peaks * 2 * pi * t / 365))) * 
        n_contacts * ((sum(I) + red_T[t] * sum(Iv) + wan_T[t] * sum(Iw)) + tau * (sum(A) + red_T[t] * sum(Av)
                                                                                  + wan_T[t] * sum(Aw))) / sum(N)
      
      # Overall time-dependent vaccination rate
      
      overall_vacc_rate <- (vacc_rate + (vacc_rate * sin(vacc_rate_shift + n_vacc_waves * 2 * pi * t / 365))) * vacc_uptake_rate
      
      # ODEs for each model state
      
      dD <- prop_hosp_mortalities[t, ] * (1 / hosp_period) * H + red_HM[t] * prop_hosp_mortalities[t, ] * (1 / hosp_period) * Hv + 
        wan_HM[t] * prop_hosp_mortalities[t, ] * (1 / hosp_period) * Hw
      
      dS <- ageing * older_S - overall_vacc_rate * S - lambda * S - ageing * S - all_cause_death_rate * S
      
      dE <- lambda * S - overall_vacc_rate * E - (1 / lat_period) * E - ageing * E + ageing * older_E - all_cause_death_rate * E
      
      dI <- prop_symp[t, ] * (1 / lat_period) * E - (1 / symp_period) * I  - overall_vacc_rate * I - ageing * I + ageing * older_I - 
        all_cause_death_rate * I
      
      dA <- (1 - prop_symp[t, ]) * (1 / lat_period) * E - rec_rate * A - overall_vacc_rate * A - ageing * A + ageing * older_A - 
        all_cause_death_rate * A
      
      dH <- prop_hosp[t, ] * (1 / symp_period) * I - (1 / hosp_period) * H - overall_vacc_rate * H - ageing * H + ageing * older_H - 
        all_cause_death_rate * H
      
      dR <- rec_rate * A + (1 - prop_hosp_mortalities[t, ]) * (1 / hosp_period) * H + (1 - prop_hosp[t, ]) * (1 / symp_period) * I - 
        wane_rate * R - overall_vacc_rate * R - red_I[t] * lambda * R - ageing * R + ageing * older_R - all_cause_death_rate * R
      
      # Vaccinated compartments
      
      dV <- overall_vacc_rate * S + overall_vacc_rate * W - wane_rate * V - red_I[t] * lambda * V - ageing * V + ageing * older_V - 
        all_cause_death_rate * V - overall_vacc_rate * V + overall_vacc_rate * V
      
      dEv <- red_I[t] * lambda * V + overall_vacc_rate * E + overall_vacc_rate * Ew + red_I[t] * lambda * Rv + red_I[t] * lambda * R + 
        wan_I[t] * lambda * Rw - wane_rate * Ev - (1 / lat_period) * Ev - ageing * Ev + ageing * older_Ev - all_cause_death_rate * Ev -
        overall_vacc_rate * Ev + overall_vacc_rate * Ev
      
      dIv <- red_S[t] * prop_symp[t, ] * (1 / lat_period) * Ev - (1 / symp_period) * Iv + overall_vacc_rate * I + overall_vacc_rate * Iw -
        ageing * Iv + ageing * older_Iv - all_cause_death_rate * Iv - overall_vacc_rate * Iv + overall_vacc_rate * Iv
      
      dAv <- (1 - red_S[t] * prop_symp[t, ]) * (1 / lat_period) * Ev - rec_rate * Av + overall_vacc_rate * A + overall_vacc_rate * Aw -
        ageing * Av + ageing * older_Av - all_cause_death_rate * Av - overall_vacc_rate * Av + overall_vacc_rate * Av
      
      dHv <- red_H[t] * prop_hosp[t, ] * (1 / symp_period) * Iv - (1 / hosp_period) * Hv + overall_vacc_rate * H + overall_vacc_rate * Hw  - 
        ageing * Hv + ageing * older_Hv - all_cause_death_rate * Hv - overall_vacc_rate * Hv + overall_vacc_rate * Hv
      
      dRv <- rec_rate * Av + (1 - red_HM[t] * prop_hosp_mortalities[t, ]) * (1 / hosp_period) * Hv + 
        (1 - (red_H[t] * prop_hosp[t, ])) * (1 / symp_period) * Iv + rec_rate * Aw + 
        (1 - wan_HM[t] * prop_hosp_mortalities[t, ]) * (1 / hosp_period) * Hw + 
        (1 - (wan_H[t] * prop_hosp[t, ])) * (1 / symp_period) * Iw + overall_vacc_rate * R + overall_vacc_rate * Rw -
        wane_rate * Rv - red_I[t] * lambda * Rv - ageing * Rv + ageing * older_Rv - all_cause_death_rate * Rv -
        overall_vacc_rate * Rv + overall_vacc_rate * Rv
      
      # Waned compartments
      
      dW <- wane_rate * V - overall_vacc_rate * W - wan_I[t] * lambda * W - ageing * W + ageing * older_W - all_cause_death_rate * W
      
      dEw <- wan_I[t] * lambda * W + wane_rate * Ev - 
        overall_vacc_rate * Ew - (1 / lat_period) * Ew - ageing * Ew + ageing * older_Ew - all_cause_death_rate * Ew
      
      dIw <- wan_S[t] * prop_symp[t, ] * (1 / lat_period) * Ew - (1 / symp_period) * Iw - overall_vacc_rate * Iw - ageing * Iw + 
        ageing * older_Iw - all_cause_death_rate * Iw
      
      dAw <- (1 - wan_S[t] * prop_symp[t, ]) * (1 / lat_period) * Ew - rec_rate * Aw - overall_vacc_rate * Aw - ageing * Aw + 
        ageing * older_Aw - all_cause_death_rate * Aw
      
      dHw <- wan_H[t] * prop_hosp[t, ] * (1 / symp_period) * Iw - (1 / hosp_period) * Hw - overall_vacc_rate * Hw - ageing * Hw + 
        ageing * older_Hw - all_cause_death_rate * Hw
      
      dRw <- wane_rate * R + wane_rate * Rv - overall_vacc_rate * Rw - wan_I[t] * lambda * Rw - ageing * Rw + ageing * older_Rw - 
        all_cause_death_rate * Rw
      
      # Combined infected, hospitalised, vaccinated, infected symptomatic and recovered states
      
      dII <- (1 / lat_period) * (E + Ev + Ew)
      
      dHH <- prop_hosp[t, ] * (1 / symp_period) * (I + red_H[t] * Iv + wan_H[t] * Iw)
      
      dVV <- overall_vacc_rate * (S + E + I + A + H + R + V + Ev + Iv + Av + Hv + Rv + W + Ew + Iw + Aw + Hw + Rw)
      
      dVP <- overall_vacc_rate * (S + E + I + A + H + R)
      
      dVB <- overall_vacc_rate * (V + Ev + Iv + Av + Hv + Rv + W + Ew + Iw + Aw + Hw + Rw)
      
      dIIS <- prop_symp[t, ] * (1 / lat_period) * (E + red_S[t] * Ev + wan_S[t] * Ew)
      
      dRR <- (1 - prop_hosp_mortalities[t, ]) * (1 / hosp_period) * H + (1 - prop_hosp[t, ]) * (1 / symp_period) * I
      
      dRRV <- (1 - red_HM[t] * prop_hosp_mortalities[t, ]) * (1 / hosp_period) * Hv + 
        (1 - (red_H[t] * prop_hosp[t, ])) * (1 / symp_period) * Iv
      
      dRRW <- (1 - wan_HM[t] * prop_hosp_mortalities[t, ]) * (1 / hosp_period) * Hw + 
        (1 - (wan_H[t] * prop_hosp[t, ])) * (1 / symp_period) * Iw
      
      
      return(list(c(dS, dE, dI, dA, dH, dR, dD, dV, dEv, dIv, dAv, dHv, dRv, dW, dEw, dIw, dAw, dHw, dRw, 
                    dII, dHH, dVV, dVP, dVB, dIIS, dRR, dRRV, dRRW)))	
    }
    )
  }
  
  # Function to calculate the new vaccine effectiveness, taking in the initial vaccine efficacy and the factor change
  
  calculate_new_VE <- function(VE_init, factor_change) {
    new_VE <- c(VE_init * factor_change, VE_init * (factor_change) ^ 2, VE_init * (factor_change) ^ 3, 
                VE_init * (factor_change) ^ 4)
  }
  
  # Function to calculate the new VE reduction parameters (Red_I) from the current reduction parameters and the desired VE
  
  calculate_reduction <- function(current_reduction, desired_VE){
    reduction <- (1 - desired_VE / 100) / current_reduction
  }
  
  # Function to calculate VE over the course of a year with varying number of new variants
  
  calculate_param_change <- function(n_yearly_variants, param_variant_change) {
    if (n_yearly_variants == 1) {
      protection <- c(rep(param_variant_change[1], ceiling((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))), 
                      seq(param_variant_change[1], param_variant_change[2], length.out = 30), 
                      rep(param_variant_change[2], floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
    }  else if (n_yearly_variants == 2) {
      protection <- c(rep(param_variant_change[1], ceiling((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))), 
                      seq(param_variant_change[1], param_variant_change[2], length.out = 30))
      protection <- append(protection, rep(param_variant_change[2],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[2], param_variant_change[3], length.out = 30))
      protection <- append(protection, rep(param_variant_change[3],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
    }  else if (n_yearly_variants == 3) {
      protection <- c(rep(param_variant_change[1], ceiling((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))), 
                      seq(param_variant_change[1], param_variant_change[2], length.out = 30))
      protection <- append(protection, rep(param_variant_change[2],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[2], param_variant_change[3], length.out = 30))
      protection <- append(protection, rep(param_variant_change[3],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[3], param_variant_change[4], length.out = 30))
      protection <- append(protection, rep(param_variant_change[4],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
    }  else if (n_yearly_variants == 4) {
      protection <- c(rep(param_variant_change[1], ceiling((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))), 
                      seq(param_variant_change[1], param_variant_change[2], length.out = 30))
      protection <- append(protection, rep(param_variant_change[2],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[2], param_variant_change[3], length.out = 30))
      protection <- append(protection, rep(param_variant_change[3],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[3], param_variant_change[4], length.out = 30))
      protection <- append(protection, rep(param_variant_change[4],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[4], param_variant_change[5], length.out = 30))
      protection <- append(protection, rep(param_variant_change[5],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
    } 
    protection <- ifelse(protection > 1, 1, protection)
    return (protection)
  }
  
  calculate_variant_parameters <- function(reinf_init, red_I_init, VE_S_init, VE_H_init, VE_HM_init, VE_T_init, 
                                           wan_I_init, wan_VE_S_init, wan_VE_H_init, wan_VE_HM_init, wan_VE_T_init, 
                                           beta_init, protection_factor, beta_factor, n_yearly_variants) {
    
    # Calculating vaccine effectiveness for the param change
    
    VE_reinf <- calculate_new_VE(100 * (1 - reinf_init), protection_factor)
    VE_I <- calculate_new_VE(100 * (1 - red_I_init), protection_factor)
    VE_S <- calculate_new_VE(VE_S_init, protection_factor)
    VE_H <- calculate_new_VE(VE_H_init, protection_factor)
    VE_HM <- calculate_new_VE(VE_HM_init, protection_factor)
    VE_T <- calculate_new_VE(VE_T_init, protection_factor)
    wan_VE_I <- calculate_new_VE(100 * (1 - wan_I_init), protection_factor)
    wan_VE_S <- calculate_new_VE(wan_VE_S_init, protection_factor)
    wan_VE_H <- calculate_new_VE(wan_VE_H_init, protection_factor)
    wan_VE_HM <- calculate_new_VE(wan_VE_HM_init, protection_factor)
    wan_VE_T <- calculate_new_VE(wan_VE_T_init, protection_factor)
    beta_variant_change <- calculate_new_VE(beta_init, beta_factor)
    
    # Calculating reduction parameters from the VE
    
    reinf_variant_change <- c(reinf_init, 1 - (VE_reinf / 100))
    red_I_variant_change <- c(red_I_init, 1 - (VE_I / 100))
    red_S_variant_change <- calculate_reduction(current_reduction = red_I_variant_change, desired_VE = c(VE_S_init, VE_S))
    red_H_variant_change <- calculate_reduction(current_reduction = 1 - (c(VE_S_init, VE_S) / 100), desired_VE = c(VE_H_init, VE_H))
    red_HM_variant_change <- calculate_reduction(current_reduction = 1 - (c(VE_H_init, VE_H) / 100), desired_VE = c(VE_HM_init, VE_HM))
    red_T_variant_change <- calculate_reduction(current_reduction = 1, desired_VE = c(VE_T_init, VE_T))
    wan_I_variant_change <- c(wan_I_init, 1 - (wan_VE_I / 100))
    wan_S_variant_change <- calculate_reduction(current_reduction = wan_I_variant_change, desired_VE = c(wan_VE_S_init, wan_VE_S))
    wan_H_variant_change <- calculate_reduction(current_reduction = 1 - (c(wan_VE_S_init, wan_VE_S) / 100), desired_VE = c(wan_VE_H_init, wan_VE_H))
    wan_HM_variant_change <- calculate_reduction(current_reduction = 1 - (c(wan_VE_H_init, wan_VE_H) / 100), desired_VE = c(wan_VE_HM_init, wan_VE_HM))
    wan_T_variant_change <- calculate_reduction(current_reduction = 1, desired_VE = c(wan_VE_T_init, wan_VE_T))
    beta_variant_change <- c(beta_init, beta_variant_change)
    
    # Calculating variant parameters
    
    reinf_variants <- calculate_param_change(n_yearly_variants, reinf_variant_change)
    red_I_variants <- calculate_param_change(n_yearly_variants, red_I_variant_change)
    red_S_variants <- calculate_param_change(n_yearly_variants, red_S_variant_change)
    red_H_variants <- calculate_param_change(n_yearly_variants, red_H_variant_change)
    red_HM_variants <- calculate_param_change(n_yearly_variants, red_HM_variant_change)
    red_T_variants <- calculate_param_change(n_yearly_variants, red_T_variant_change)
    wan_I_variants <- calculate_param_change(n_yearly_variants, wan_I_variant_change)
    wan_S_variants <- calculate_param_change(n_yearly_variants, wan_S_variant_change)
    wan_H_variants <- calculate_param_change(n_yearly_variants, wan_H_variant_change)
    wan_HM_variants <- calculate_param_change(n_yearly_variants, wan_HM_variant_change)
    wan_T_variants <- calculate_param_change(n_yearly_variants, wan_T_variant_change)
    beta_variants <- calculate_param_change(n_yearly_variants, beta_variant_change)
    
    return (list(reinf_variants, red_I_variants, red_S_variants, red_H_variants, red_HM_variants, red_T_variants, 
                 wan_I_variants, wan_S_variants, wan_H_variants, wan_HM_variants, wan_T_variants, beta_variants))
  }
  
  calculate_severity_parameters_scenario_3 <- function (severity_init, severity_factor_scenario_3) {
    severity_scenario_3 <- c(rep(severity_init, ceiling((365 - 30) / 2)), seq(severity_init, severity_init * severity_factor_scenario_3, 
                                                                              length.out = 30), rep(severity_init * severity_factor_scenario_3,  
                                                                                                    floor((365 - 30) / 2)))
    severity_scenario_3 <- ifelse(severity_scenario_3 > 1, 1, severity_scenario_3)
    return (severity_scenario_3)
  }
  
  # Set seed
  
  set.seed(par_fit[1])
  
  # Set variant scenario
  
  variant_scenario <- "No variants"
  n_yearly_variants <- 1
  vaccination_scenario <- "B9"      # arbitrary choice as the future simlulation time is 0 years
  
  # Parameters
  
  # Sizes of each age group in the population of England, using the output of the population burn-in function for 500 years
  
  n_age_groups <- 101
  age_dist <- c(696785, 696646, 696507, 696367, 696228, 696159, 696089, 696019, 695950, 695880, 695811, 695741, 695671, 
                695602, 695532, 695324, 695115, 694907, 694698, 694490, 694212, 693935, 693657, 693380, 693103, 692826, 
                692549, 692272, 691995, 691718, 691234, 690751, 690268, 689785, 689302, 688614, 687926, 687238, 686552, 685866, 
                684839, 683813, 682789, 681766, 680745, 679183, 677624, 676069, 674518, 672970, 670690, 668417, 666152, 663895, 
                661645, 658354, 655078, 651819, 648576, 645350, 640482, 635651, 630856, 626098, 621376, 613947, 606607, 599355, 
                592189, 585109, 574031, 563162, 552498, 542037, 531774, 515535, 499791, 484528, 469732, 455387, 430301, 406596, 
                384197, 363033, 343034, 310185, 280482, 253623, 229337, 207376, 169134, 137945, 112507, 91760, 74839, 61038, 
                49782, 40602, 33115, 27008, 27008)
  total_pop <- total_original_pop <- sum(age_dist)
  age_dist_proportion <- age_dist / total_pop
  
  # On average there are 630000 deaths during the burn-in period. Add these to the initial population so that these are not lost
  # over the burn-in.
  
  age_dist <- age_dist + 630000 * age_dist_proportion
  age_dist <<- age_dist
  total_pop <- sum(age_dist)
  
  # Time horizon for burn-in period
  
  n_years_0 <- 25                  # Run the model initially for 25 years
  SIMTIME_0 <- n_years_0 * 365
  times_0 <- seq(1, SIMTIME_0, by = 1)
  
  # Time horizon for fitting period (2023)
  
  n_years_fit <- 1
  SIMTIME_fit <- n_years_fit * 365
  times_fit <- seq(1, SIMTIME_fit, by = 1)
  
  # Model run time horizon
  
  n_years_future <- 0
  SIMTIME_future <- n_years_future * 365
  if (SIMTIME_future > 0) {
    times_future <- seq(1, SIMTIME_future, by = 1)
  }
  
  SIMTIME <- SIMTIME_fit + SIMTIME_future
  times <- seq(1, SIMTIME, by = 1)
  
  # Vaccination rates for fitting period
  
  vacc_rate_4 <- rep(0, 65)                                          # Vaccination rate for 0 to 64 year olds (no vaccination)
  vacc_rate_5 <- rep(1 / 365, 10)                                      # Vaccination rate for 65 to 74 year olds
  vacc_rate_6 <- rep(2 / 365, 26)                                      # Vaccination rate for 75 to 100 year olds
  
  # Set vacc rates in burn-in period to same as fitting period
  
  vacc_rate_1 <- vacc_rate_4
  vacc_rate_2 <- vacc_rate_5
  vacc_rate_3 <- vacc_rate_6
  
  # Number of vaccination waves for burn-in and fitting period
  
  n_vacc_waves_fitting <- c(rep(0, 65), rep(1, 10), rep(2, 26))      # Two boosters for 75+, one booster for 65-74, no vaccination for under 64s
  
  n_vacc_waves_burn_in <- n_vacc_waves_fitting
  
  # Setting up the vaccination scenario
  
  vacc_rate_scenarios <- list(Comparator <- rep(0, 101),
                              A1 <- c(rep(0, 5), rep(1 / 365, 96)),
                              A2 <- c(rep(0, 18), rep(1 / 365, 83)),
                              A3 <- c(rep(0, 50), rep(1 / 365, 51)),
                              A4 <- c(rep(0, 55), rep(1 / 365, 46)),
                              A5 <- c(rep(0, 60), rep(1 / 365, 41)),
                              A6 <- c(rep(0, 65), rep(1 / 365, 36)),
                              A7 <- c(rep(0, 70), rep(1 / 365, 31)),
                              A8 <- c(rep(0, 75), rep(1 / 365, 26)),
                              A9 <- c(rep(0, 80), rep(1 / 365, 21)),
                              B1 <- c(rep(0, 5), rep(2 / 365, 96)),
                              B2 <- c(rep(0, 18), rep(2 / 365, 83)),
                              B3 <- c(rep(0, 50), rep(2 / 365, 51)),
                              B4 <- c(rep(0, 55), rep(2 / 365, 46)),
                              B5 <- c(rep(0, 60), rep(2 / 365, 41)),
                              B6 <- c(rep(0, 65), rep(2 / 365, 36)),
                              B7 <- c(rep(0, 70), rep(2 / 365, 31)),
                              B8 <- c(rep(0, 75), rep(2 / 365, 26)),
                              B9 <- c(rep(0, 80), rep(2 / 365, 21)),
                              C1 <- c(rep(0, 5), rep(1 / 730, 96)),
                              C2 <- c(rep(0, 18), rep(1 / 730, 83)),
                              C3 <- c(rep(0, 50), rep(1 / 730, 51)),
                              C4 <- c(rep(0, 55), rep(1 / 730, 46)),
                              C5 <- c(rep(0, 60), rep(1 / 730, 41)),
                              C6 <- c(rep(0, 65), rep(1 / 730, 36)),
                              C7 <- c(rep(0, 70), rep(1 / 730, 31)),
                              C8 <- c(rep(0, 75), rep(1 / 730, 26)),
                              C9 <- c(rep(0, 80), rep(1 / 730, 21)),
                              "fitting_period" <- c(vacc_rate_4, vacc_rate_5, vacc_rate_6))
  
  n_vacc_waves_scenarios <- list(Comparator <- rep(0, 101),
                                 A1 <- c(rep(0, 5), rep(1, 96)),
                                 A2 <- c(rep(0, 18), rep(1, 83)),
                                 A3 <- c(rep(0, 50), rep(1, 51)),
                                 A4 <- c(rep(0, 55), rep(1, 46)),
                                 A5 <- c(rep(0, 60), rep(1, 41)),
                                 A6 <- c(rep(0, 65), rep(1, 36)),
                                 A7 <- c(rep(0, 70), rep(1, 31)),
                                 A8 <- c(rep(0, 75), rep(1, 26)),
                                 A9 <- c(rep(0, 80), rep(1, 21)),
                                 B1 <- c(rep(0, 5), rep(2, 96)),
                                 B2 <- c(rep(0, 18), rep(2, 83)),
                                 B3 <- c(rep(0, 50), rep(2, 51)),
                                 B4 <- c(rep(0, 55), rep(2, 46)),
                                 B5 <- c(rep(0, 60), rep(2, 41)),
                                 B6 <- c(rep(0, 65), rep(2, 36)),
                                 B7 <- c(rep(0, 70), rep(2, 31)),
                                 B8 <- c(rep(0, 75), rep(2, 26)),
                                 B9 <- c(rep(0, 80), rep(2, 21)),
                                 C1 <- c(rep(0, 5), rep(0.5, 96)),
                                 C2 <- c(rep(0, 18), rep(0.5, 83)),
                                 C3 <- c(rep(0, 50), rep(0.5, 51)),
                                 C4 <- c(rep(0, 55), rep(0.5, 46)),
                                 C5 <- c(rep(0, 60), rep(0.5, 41)),
                                 C6 <- c(rep(0, 65), rep(0.5, 36)),
                                 C7 <- c(rep(0, 70), rep(0.5, 31)),
                                 C8 <- c(rep(0, 75), rep(0.5, 26)),
                                 C9 <- c(rep(0, 80), rep(0.5, 21)),
                                 "fitting_period" <- n_vacc_waves_fitting)
  
  vacc_rate_shift_scenarios <- list(Comparator <- 2 * pi / 3,        # Vaccination peaks ~ 10th November
                                    A1 <- 2 * pi / 3,
                                    A2 <- 2 * pi / 3,
                                    A3 <- 2 * pi / 3,
                                    A4 <- 2 * pi / 3,
                                    A5 <- 2 * pi / 3,
                                    A6 <- 2 * pi / 3,
                                    A7 <- 2 * pi / 3,
                                    A8 <- 2 * pi / 3,
                                    A9 <- 2 * pi / 3,
                                    B1 <- pi,                        # Vaccination peaks ~ 23 October and 21st April
                                    B2 <- pi,
                                    B3 <- pi,
                                    B4 <- pi,
                                    B5 <- pi,
                                    B6 <- pi,
                                    B7 <- pi,
                                    B8 <- pi,
                                    B9 <- pi,
                                    C1 <- 2 * pi / 3,
                                    C2 <- 2 * pi / 3,
                                    C3 <- 2 * pi / 3,
                                    C4 <- 2 * pi / 3,
                                    C5 <- 2 * pi / 3,
                                    C6 <- 2 * pi / 3,
                                    C7 <- 2 * pi / 3,
                                    C8 <- 2 * pi / 3,
                                    C9 <- 2 * pi / 3,
                                    "fitting_period" <- pi)
  
  names(vacc_rate_scenarios) <- names(n_vacc_waves_scenarios) <- names(vacc_rate_shift_scenarios) <- 
    c("Comparator", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", 
      "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "fitting_period")
  
  vacc_rate_future <- as.numeric(unlist(vacc_rate_scenarios[vaccination_scenario]))
  n_vacc_waves_future <- as.numeric(unlist(n_vacc_waves_scenarios[vaccination_scenario]))
  vacc_rate_shift_future <- as.numeric(unlist(vacc_rate_shift_scenarios[vaccination_scenario]))
  
  # Contact rates based on earlier bootstrapping
  
  contacts <- c(0.85, 0.96, 1.07, 1.18, 1.29, 1.39, 1.48, 1.55, 1.62, 1.66, 1.7, 1.71, 1.72, 1.72, 1.7, 1.68, 1.64, 1.6, 
                1.56, 1.51, 1.46, 1.42, 1.38, 1.35, 1.32, 1.3, 1.28, 1.27, 1.26, 1.25, 1.25, 1.25, 1.25, 1.26, 1.26, 1.27, 
                1.27, 1.28, 1.28, 1.27, 1.27, 1.27, 1.27, 1.26, 1.25, 1.24, 1.22, 1.2, 1.17, 1.15, 1.12, 1.1, 1.08, 1.06, 
                1.04, 1.02, 1.01, 0.99, 0.97, 0.95, 0.93, 0.91, 0.89, 0.87, 0.84, 0.82, 0.79, 0.76, 0.73, 0.7, 0.68, 0.65, 
                0.63, 0.61, 0.6, 0.58, 0.58, 0.57, 0.57, 0.56, 0.56, 0.56, 0.56, 0.56, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 
                0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57)
  
  # Calculating parameter factor changes according to variant scenario
  
  if (variant_scenario == "Base case" | variant_scenario == "Scenario 2" | variant_scenario == "Scenario 3") {
    protection_factor <- 0.92
    beta_factor <- 1
  } else if (variant_scenario == "Scenario 4") {
    n_yearly_variants <- 3
    protection_factor <- 0.95
    beta_factor <- 1
  } else {
    print("UNRECOGNISED SCENARIO")
    protection_factor <- 1
    beta_factor <- 1
  }
  
  protection_factor_scenario_2 <- 0.92
  beta_factor_scenario_2 <- 1.10
  prop_symp_factor_scenario_3 <- 1
  prop_hosp_factor_scenario_3 <- 1.80
  prop_hosp_mortalities_factor_scenario_3 <- 1.80
  
  VE_S_init <- 65
  VE_H_init <- 90
  VE_HM_init <- 90
  VE_T_init <- 68
  wan_VE_H_init <- 60
  wan_VE_HM_init <- 60
  wan_VE_T_init <- 55
  
  # Infection severity rates
  
  inf_symp_rate_init <- 0.803
  inf_hosp_rate_init <- 0.32
  hosp_mort_rate_init <- 0.499
  inf_symp_age_scaling <- c(rep(0.664, 19), rep(0.846, 41), rep(1, 41))
  age_scaling <- c(rep(0.00513, 6), rep(0.00563, 12), rep(0.074, 47), rep(0.488, 20), rep(1, 16))
  
  inf_symp_rate_base_case <- rep(inf_symp_rate_init, max(SIMTIME_future, SIMTIME_0))
  inf_hosp_rate_base_case <- rep(inf_hosp_rate_init, max(SIMTIME_future, SIMTIME_0))
  hosp_mort_rate_base_case <- rep(hosp_mort_rate_init, max(SIMTIME_future, SIMTIME_0))
  
  # Setting up increased severity in variant scenario 3
  
  if (variant_scenario == "Scenario 3") {
    
    fifth_element <- ((seq(1:floor(n_years_future / 5)) * 5) - 1) * 365
    
    # Setting new severity rates for variant scenario 3
    
    inf_symp_rate_scenario_3 <- calculate_severity_parameters_scenario_3(inf_symp_rate_init, prop_symp_factor_scenario_3)
    inf_hosp_rate_scenario_3 <- calculate_severity_parameters_scenario_3(inf_hosp_rate_init, prop_hosp_factor_scenario_3)
    hosp_mort_rate_scenario_3 <- calculate_severity_parameters_scenario_3(hosp_mort_rate_init, prop_hosp_mortalities_factor_scenario_3)
    
    # Setting severity rates over the time horizon to the base case
    
    inf_symp_rate <- inf_symp_rate_base_case
    inf_hosp_rate <- inf_hosp_rate_base_case
    hosp_mort_rate <- hosp_mort_rate_base_case
    
    # Replace the severity rates for every fifth year with the scenario 3 variant severity rates (only one variant in every fifth year)
    
    for (j in 1:length(fifth_element)) {
      inf_symp_rate[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- inf_symp_rate_scenario_3
      inf_hosp_rate[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- inf_hosp_rate_scenario_3
      hosp_mort_rate[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- hosp_mort_rate_scenario_3
    }
    
  } else {
    
    # If variant scenario is not scenario 3
    
    inf_symp_rate <- inf_symp_rate_base_case
    inf_hosp_rate <- inf_hosp_rate_base_case
    hosp_mort_rate <- hosp_mort_rate_base_case
    
  }
  
  # Set fixed parameters
  
  par <- list(ageing <- n_age_groups / (365 * 101),
              lat_period <- 3.1,
              prop_symp <- outer(inf_symp_rate, inf_symp_age_scaling, "*"),
              symp_period <- 9,
              rec_rate <- 1 / symp_period,
              prop_hosp <- outer(inf_hosp_rate, age_scaling, "*"),
              hosp_period <- c(rep(1.8, 5), rep(2.8, 7), rep(5, 6), rep(6.4, 32), rep(9, 15), rep(10.8, 10), rep(11.2, 10), rep(11.4, 16)),
              prop_hosp_mortalities <- outer(hosp_mort_rate, age_scaling, "*"),
              vacc_rate <- c(vacc_rate_1, vacc_rate_2, vacc_rate_3),
              vacc_rate_shift <- pi,
              n_vacc_waves <- 0,         # No vaccination for under 12s, two vaccinations a year for 12-100 year olds
              vacc_uptake_rate <- c(rep(0, 5), rep(6.4, 7), rep(31.7, 4), rep(31.7, 2), rep(31.7, 2), rep(40.3, 5), rep(40.9, 5), 
                                    rep(44.0, 5), rep(49.5, 5), rep(57.0, 5), rep(64.8, 5), rep(64.8, 5), rep(64.8, 5), 
                                    rep(64.8, 5), rep(64.8, 5), rep(70.1, 5), rep(75.5, 5), rep(75.7, 21)) / 100,
              # Using uptake rates for 0-49 year olds from third dose uptake, 50-64 from autumn 2022 uptake and 65+ from autumn 2023 uptake
              wane_rate <- 2 / 365,                      # Immunity wanes after six months
              red_I <- 0,
              red_S <- 0,
              red_H <- 0,
              red_HM <- 0,
              wan_I <- 0,
              wan_S <- 0,
              wan_H <- 0,
              wan_HM <- 0,
              reinf <- 0,
              red_T <- 0,
              wan_T <- 0,
              tau <- 0,
              beta <- 0,
              num_peaks <- 2,
              inf_shift <- pi / 4,
              inf_fluctuation <- 0.125,
              n_contacts <- contacts,
              total_pop <- total_pop,
              n_age_groups <- n_age_groups,
              all_cause_death_rate <- c(4.2, rep(0.2, 4), rep(0.1, 5), rep(0.1, 5), rep(0.3, 5), rep(0.4, 5), rep(0.4, 5), rep(0.7, 5), rep(1.0, 5),
                                        rep(1.5, 5), rep(2.3, 5), rep(3.4, 5), rep(5.0, 5), rep(7.6, 5), rep(12.1, 5), rep(19.3, 5),
                                        rep(31.5, 5), rep(58.3, 5), rep(105.9, 5), rep(226.1, 10), 0) / 365000,          
              # Using yearly mortality rates per 1000 people from 2023 from the ONS, so dividing by 365*1000 (https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsregisteredinenglandandwalesseriesdrreferencetables)
              birth_rate <- 700000)                               # Number of projected yearly births in England proportional to the UK
  
  names(par) <- c("Ageing", "Latent period", "Proportion symptomatic", "Symptomatic period", "Recovery rate (asymptomatic)", "Proportion hospitalised",
                  "Hospital length of stay", "Proportion of hospital mortalities", "Vaccination rate", "Vaccination rate shift",
                  "Number of vaccination waves", "Vaccine uptake rate", "Wane rate", "Reduction in infectiousness", "Reduction in symptomatic infections",
                  "Reduction in hospitalisations", "Reduction in hospital mortalities", "Waned reduction in infectiousness", 
                  "Waned reduction in symptomatic infections", "Waned reduction in hospitalisations", "Waned reduction in hospital mortalities", 
                  "Protection against reinfection", "Reduction in transmission", "Waned reduction in transmission", 
                  "Reduction in asymptomatic infectiousness", "Transmission rate", "Number of annual peaks in transmission",
                  "Infection shift", "Amount of fluctuation in infectiousness", "Number of contacts per age group", "Total population size",
                  "Number of age groups", "All-cause death rate", "Birth rate")
  
  all_cause_death_rate <<- par$all_cause_death_rate
  
  # Initial state values for the ODE system, followed by the time sequence over which the ODEs are solved
  
  # Set initial number infected to the number of people testing positive in England week ending 28 December 2022, distributed
  # according to contact rates (https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/bulletins/coronaviruscovid19infectionsurveypilot/6january2023)
  
  n_cases <- 2463000
  I_init <- round(n_cases * (contacts/sum(contacts)))
  
  # Set initial number vaccinated to the percentage of the population in each age group that had received a third dose on 31 Dec 2022 - the number infected
  
  V_init <- round(c(rep(0, 5), rep(0.002, 7), rep(0.011, 4), rep(0.072, 2), rep(0.239, 2), rep(0.400, 5), rep(0.406, 5), 
                    rep(0.437, 5), rep(0.492, 5), rep(0.568, 5), rep(0.646, 5), rep(0.737, 5), rep(0.794, 5), rep(0.836, 5), 
                    rep(0.874, 5), rep(0.909, 5), rep(0.936, 5), rep(0.941, 21)) * (age_dist - I_init), 2)
  
  # Set initial number waned to the percentage in each age group that had received a second dose on 31 Dec 2022 and had not received a third - the number infected
  
  W_init <- round(age_dist * c(rep(0, 5), rep(0.056, 7), rep(0.300, 4), rep(0.400, 2), rep(0.333, 2), rep(0.244, 5), rep(0.225, 5), 
                               rep(0.210, 5), rep(0.192, 5), rep(0.168, 5), rep(0.144, 5), rep(0.110, 5), rep(0.085, 5), 
                               rep(0.065, 5), rep(0.045, 5), rep(0.029, 5), rep(0.020, 5), rep(0.018, 21)), 2)
  
  n_recovered <- 19037525                    # Approx cumulative number of first episodes in England as of 31 December 2022
  
  E_init <- rep(0, n_age_groups)
  A_init <- 0 * age_dist
  H_init <- 0 * age_dist
  R_init <- round(pmin(n_recovered * (contacts/sum(contacts)), age_dist - I_init - V_init - W_init), 2)      
  # Number of recovered in each age group equal to minimum of estimated number recovered and the remaining population 
  # in that age group that's not infected, vaccinated, or waned
  D_init <- rep(0, n_age_groups)
  Ev_init <- rep(0, n_age_groups)
  Iv_init <- 0 * age_dist
  Av_init <- 0 * age_dist
  Hv_init <- rep(0, n_age_groups)
  Rv_init <- 0 * age_dist
  Ew_init <- rep(0, n_age_groups)
  Iw_init <- 0 * age_dist
  Aw_init <- 0 * age_dist
  Hw_init <- rep(0, n_age_groups)
  Rw_init <- 0 * age_dist
  
  S_init <- round(age_dist - (E_init + I_init + A_init + H_init + R_init + D_init + V_init + Ev_init + Iv_init + Av_init + 
                                Hv_init + Rv_init + W_init + Ew_init + Iw_init + Aw_init + Hw_init + Rw_init), 2)
  
  II_init <- rep(0, n_age_groups)
  HH_init <- rep(0, n_age_groups)
  VV_init <- rep(0, n_age_groups)
  VP_init <- rep(0, n_age_groups)
  VB_init <- rep(0, n_age_groups)
  IIS_init <- rep(0, n_age_groups)
  RR_init <- rep(0, n_age_groups)
  RRV_init <- rep(0, n_age_groups)
  RRW_init <- rep(0, n_age_groups)
  
  # Model burn-in input
  
  input_init <- c(S_init, E_init, I_init, A_init, H_init, R_init, D_init, V_init, Ev_init, Iv_init, Av_init, Hv_init, 
                  Rv_init, W_init, Ew_init, Iw_init, Aw_init, Hw_init, Rw_init, II_init, HH_init, VV_init, VP_init, 
                  VB_init, IIS_init, RR_init, RRV_init, RRW_init)
  
  # Check for negative values in model initial input and checking population is constant
  
  combined_pop_compartments_init <- S_init + E_init + I_init + A_init + H_init + R_init + D_init + V_init + 
    Ev_init + Iv_init + Av_init + Hv_init + Rv_init + W_init + Ew_init + 
    Iw_init + Aw_init + Hw_init + Rw_init
  
  if(is.null(fixr::check_for_negative_values(input_init)) == FALSE) {
    stop("Negative values in input_init")
  } else if (round(sum(combined_pop_compartments_init)) != round(total_pop)) {
    stop("Total population in input_init not equal to initial population")
  }
  # 
  # # Initialising hospitalisation counting, dates, parameter sets and hospital data set vectors
  # 
  # inc_hosp <- inc_deaths <- inc_hosp_0_5 <- inc_deaths_0_5 <- inc_hosp_6_17 <- inc_deaths_6_17 <- inc_hosp_18_64 <- 
  #   inc_deaths_18_64 <- inc_hosp_65_84 <- inc_deaths_65_84 <- inc_hosp_85_ <- inc_deaths_85_ <- inc_vacc <- matrix(nrow = length(times) - 1, ncol = ncol(par_fit))
  # cum_hosp <- cum_deaths <- matrix(nrow = length(times), ncol = ncol(par_fit))
  # total_inf_future <- total_vacc_future <- total_hosp_future <- total_deaths_future <- total_inf_future_0_5 <- total_inf_future_6_17 <-
  #   total_inf_future_18_64 <- total_inf_future_65_84 <- total_inf_future_85_ <- total_hosp_future_0_5 <- total_hosp_future_6_17 <-
  #   total_hosp_future_18_64 <- total_hosp_future_65_84 <- total_hosp_future_85_<- total_deaths_future_0_5 <- total_deaths_future_6_17 <-
  #   total_deaths_future_18_64 <- total_deaths_future_65_84 <- total_deaths_future_85_ <- total_icu_future_0_5 <-
  #   total_icu_future_6_17 <- total_icu_future_18_64 <- total_icu_future_65_84 <- total_icu_future_85_ <-life_years_lost_future <-
  #   rep(NA, ncol(par_fit))
  # 
  # # Initialising counting for number in each age group in each disease state every year in future horizon for counting costs/QALYs
  # 
  # n_healthy_yearly <- n_symptomatic_yearly <- n_hospitalised_yearly <- n_dead_yearly <-  n_recovered_past_year_yearly <- 
  #   n_recovered_vaccinated_past_year_yearly <- n_recovered_waned_past_year_yearly <- inc_inf_yearly <- inc_hosp_yearly <- 
  #   inc_deaths_yearly <- inc_vacc_yearly <- inc_symp_inf_yearly <- inc_rec_yearly <- vector(mode='list', length = ncol(par_fit))
  # 
  # # Setting dates vector
  # 
  # dates <- seq(as.Date("2022-12-31"), by = "weeks", length = floor(SIMTIME / 7))
  # breaks.vec <- seq(min(dates), max(dates), by = "2 months")
  # 
  # parameter_set <- rep(paste("Parameter set", "1"), floor(SIMTIME / 7))
  # if (ncol(par_fit) > 1) {
  #   for (i in 2:ncol(par_fit)) {
  #     parameter_set <- append(parameter_set, rep(paste("Parameter set", i), floor(SIMTIME / 7)))
  #   }
  # }
  # 
  # # Extracting hospital admissions data from excel file
  # 
  # hosp_data <- read_excel("Data/Weekly hospital admissions 2023 (England).xlsx", sheet = "2023")
  # deaths_data <- read_excel("Data/Weekly deaths 2023 (England).xlsx", sheet = "2023")
  # 
  # # Daily hospital admissions for all age groups
  # 
  # hosp_data_daily <- as.data.frame(hosp_data)
  # hosp_data_daily <- hosp_data_daily[-c(1), -c(8,9)]
  # hosp_data_daily[, 1] <- excel_numeric_to_date(as.numeric(as.character(hosp_data_daily[, 1])), date_system = "modern")
  # colnames(hosp_data_daily) <- c("Date", "Total", "0-5", "6-17", "18-64", "65-84", "85+")
  # 
  # hosp_data_total <<- colSums(hosp_data_daily[, 2:7])
  # 
  # # Daily deaths for all age groups
  # 
  # deaths_data_weekly <- as.data.frame(deaths_data[-53,])
  # colnames(deaths_data_weekly) <- c("Date", "Total", "0-5", "6-17", "18-64", "65-84", "85+")
  # deaths_data_weekly <<- deaths_data_weekly
  # 
  # deaths_data_total <<- colSums(deaths_data_weekly[, 2:7])
  # 
  # # Summing daily admissions to give weekly admissions for all groups
  # 
  # hosp_data_weekly <<- rowsum(hosp_data_daily[1:364,], rep(1:52, each = 7))
  # 
  # # Calculating number of deaths during the burn-in
  # 
  # burn_in_deaths <- rep(0, ncol(par_fit))
  # 
  # # Run model for each parameter set
  # 
  # for (i in 1:ncol(par_fit)) {
  #   
  #   print(paste("Parameter set", i))
  #   
  # Setting parameters
  
  par[["Reduction in asymptomatic infectiousness"]] <- tau <- par_fit[4]
  
  age_scaling <- c(rep(par_fit[5], 6), rep(par_fit[6], 12),
                   rep(par_fit[7], 47), rep(par_fit[8], 20), 
                   rep(par_fit[9], 16))
  
  par[["Proportion symptomatic"]] <- prop_symp <- outer(inf_symp_rate_base_case[1:SIMTIME_0], inf_symp_age_scaling, "*")
  
  par[["Proportion hospitalised"]] <- prop_hosp <- outer(inf_hosp_rate_base_case[1:SIMTIME_0], age_scaling, "*")
  
  par[["Proportion of hospital mortalities"]] <- prop_hosp_mortalities <- 
    outer(hosp_mort_rate_base_case[1:SIMTIME_0], age_scaling * c(rep(par_fit[10], 6), 
                                                                 rep(par_fit[11], 12),
                                                                 rep(par_fit[12], 47), 
                                                                 rep(par_fit[13], 20), 
                                                                 rep(1, 16)), "*")
  
  n_vacc_waves <- par[["Number of vaccination waves"]] <- n_vacc_waves_burn_in         # For burn-in period: no vaccination for under 12s, two vaccinations a year for 12-100 year olds
  
  vacc_rate <- par[["Vaccination rate"]] <- c(vacc_rate_1, vacc_rate_2, vacc_rate_3)       # Set original vaccination rates for burn-in period
  
  vacc_rate_shift <- par[["Vaccination rate shift"]] <- pi                             # Shift for 2 vacc peaks a year so vacc peaks in Apr-Jun and Sep-Dec
  
  # Setting fitted parameters according to the variant scenario
  
  par[["Transmission rate"]] <- beta <- rep(par_fit[2], max(SIMTIME_0, SIMTIME_fit))
  par[["Protection against reinfection"]] <- reinf <- rep(par_fit[3], max(SIMTIME_0, SIMTIME_fit))
  
  # Setting initial values for fitted parameters
  
  reinf_init <- par_fit[3]
  beta_init <- par_fit[2]
  
  # Setting VE against infection to equal natural immunity against infection
  
  red_I_init <- reinf_init
  VE_I_init <- 100 * (1 - red_I_init)
  
  # Setting waned VE against infection to fitted parameter * VE against infection
  
  wan_VE_I_init <- VE_I_init * par_fit[14]
  wan_I_init <- 1 - (wan_VE_I_init / 100)
  red_I <- par[["Reduction in infectiousness"]] <- rep(red_I_init, max(SIMTIME_0, SIMTIME_fit))
  wan_I <- par[["Waned reduction in infectiousness"]] <- rep(wan_I_init, max(SIMTIME_0, SIMTIME_fit))
  
  # Setting waned VE against symptoms to equal 10% or at least the waned VE against infection
  
  wan_VE_S_init <- max(10, wan_VE_I_init)
  
  # Calculating vaccine effectiveness initial values from fitted parameters
  
  red_S_init <- calculate_reduction(current_reduction = red_I_init, desired_VE = VE_S_init)
  red_H_init <- calculate_reduction(current_reduction = 1 - (VE_S_init / 100), desired_VE = VE_H_init)
  red_HM_init <- calculate_reduction(current_reduction = 1 - (VE_H_init / 100), desired_VE = VE_HM_init)
  red_T_init <- calculate_reduction(current_reduction = 1, desired_VE = VE_T_init)
  wan_S_init <- calculate_reduction(current_reduction = wan_I_init, desired_VE = wan_VE_S_init)
  wan_H_init <- calculate_reduction(current_reduction = 1 - (wan_VE_S_init / 100), desired_VE = wan_VE_H_init)
  wan_HM_init <- calculate_reduction(current_reduction = 1 - (wan_VE_H_init / 100), desired_VE = wan_VE_HM_init)
  wan_T_init <- calculate_reduction(current_reduction = 1, desired_VE = wan_VE_T_init)
  
  # Setting vaccine effectiveness to the initial value for the burn-in period
  
  red_S <- par[["Reduction in symptomatic infections"]] <- rep(red_S_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period
  red_H <- par[["Reduction in hospitalisations"]] <- rep(red_H_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period
  red_HM <- par[["Reduction in hospital mortalities"]] <- rep(red_HM_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period
  red_T <- par[["Reduction in transmission"]] <- rep(red_T_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period
  wan_S <- par[["Waned reduction in symptomatic infections"]] <- rep(wan_S_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period/fitting period, whichever is longer
  wan_H <- par[["Waned reduction in hospitalisations"]] <- rep(wan_H_init, max(SIMTIME_0, SIMTIME_fit))           # Set waned vaccine effectiveness to the initial value for the burn-in period
  wan_HM <- par[["Waned reduction in hospital mortalities"]] <- rep(wan_HM_init, max(SIMTIME_0, SIMTIME_fit))           # Set waned vaccine effectiveness to the initial value for the burn-in period
  wan_T <- par[["Waned reduction in transmission"]] <- rep(wan_T_init, max(SIMTIME_0, SIMTIME_fit))           # Set waned vaccine effectiveness to the initial value for the burn-in period
  
  # Data frame containing solutions to the ODEs over the specified time sequence (burn-in period)
  
  out_0 <- as.data.frame(deSolve::ode(y = input_init, times = times_0, func = covid_model_odes, parms = par, method = "euler"))
  
  # Label the output columns for each compartment
  
  S_0 <- out_0[, 2:(n_age_groups + 1)]
  E_0 <- out_0[, (n_age_groups + 2):(2 * n_age_groups + 1)]
  I_0 <- out_0[, (2 * n_age_groups + 2):(3 * n_age_groups + 1)]
  A_0 <- out_0[, (3 * n_age_groups + 2):(4 * n_age_groups + 1)]
  H_0 <- out_0[, (4 * n_age_groups + 2):(5 * n_age_groups + 1)]
  R_0 <- out_0[, (5 * n_age_groups + 2):(6 * n_age_groups + 1)]
  D_0 <- out_0[, (6 * n_age_groups + 2):(7 * n_age_groups + 1)]
  V_0 <- out_0[, (7 * n_age_groups + 2):(8 * n_age_groups + 1)]
  Ev_0 <- out_0[, (8 * n_age_groups + 2):(9 * n_age_groups + 1)]
  Iv_0 <- out_0[, (9 * n_age_groups + 2):(10 * n_age_groups + 1)]
  Av_0 <- out_0[, (10 * n_age_groups + 2):(11 * n_age_groups + 1)]
  Hv_0 <- out_0[, (11 * n_age_groups + 2):(12 * n_age_groups + 1)]
  Rv_0 <- out_0[, (12 * n_age_groups + 2):(13 * n_age_groups + 1)]
  W_0 <- out_0[, (13 * n_age_groups + 2):(14 * n_age_groups + 1)]
  Ew_0 <- out_0[, (14 * n_age_groups + 2):(15 * n_age_groups + 1)]
  Iw_0 <- out_0[, (15 * n_age_groups + 2):(16 * n_age_groups + 1)]
  Aw_0 <- out_0[, (16 * n_age_groups + 2):(17 * n_age_groups + 1)]
  Hw_0 <- out_0[, (17 * n_age_groups + 2):(18 * n_age_groups + 1)]
  Rw_0 <- out_0[, (18 * n_age_groups + 2):(19 * n_age_groups + 1)]
  II_0 <- out_0[, (19 * n_age_groups + 2):(20 * n_age_groups + 1)]
  HH_0 <- out_0[, (20 * n_age_groups + 2):(21 * n_age_groups + 1)]
  VV_0 <- out_0[, (21 * n_age_groups + 2):(22 * n_age_groups + 1)]
  VP_0 <- out_0[, (22 * n_age_groups + 2):(23 * n_age_groups + 1)]
  VB_0 <- out_0[, (23 * n_age_groups + 2):(24 * n_age_groups + 1)]
  IIS_0 <- out_0[, (24 * n_age_groups + 2):(25 * n_age_groups + 1)]
  RR_0 <- out_0[, (25 * n_age_groups + 2):(26 * n_age_groups + 1)]
  RRV_0 <- out_0[, (26 * n_age_groups + 2):(27 * n_age_groups + 1)]
  RRW_0 <- out_0[, (27 * n_age_groups + 2):(28 * n_age_groups + 1)]
  
  # Check for negative values in model burn-in output
  
  if(is.null(fixr::check_for_negative_values(out_0)) == FALSE) {
    stop("Negative values in out_0")
  }
  
  # Run the model again with the previous output as the input
  
  input_fit <- as.numeric(cbind(S_0, E_0, I_0, A_0, H_0, R_0, D_0, V_0, Ev_0, Iv_0, Av_0, Hv_0, Rv_0, W_0, Ew_0, Iw_0,
                                Aw_0, Hw_0, Rw_0, II_0, HH_0, VV_0, VP_0, VB_0, IIS_0, RR_0, RRV_0, RRW_0)[SIMTIME_0, ])
  
  # Reset the Covid deaths compartments in the input to 0
  
  input_fit[607:707] <- rep(0, 101)
  
  # burn_in_deaths[i] <- sum(D_0[nrow(D_0), ])
  
  # Replenish Covid deaths from the burn-in and reintroduce to the susceptible population with the same age distribution as the original
  
  # input_fit <- input_fit + c(sum(D_0[456, ]) * age_dist_proportion, rep(0, (length(input_fit) - 101)))
  
  # Second model run parameters (fitting period)
  
  n_vacc_waves <- par[["Number of vaccination waves"]] <- n_vacc_waves_fitting      # Two boosters for 75+, one booster for 50-74, no vaccination for under 50
  vacc_rate <- par[["Vaccination rate"]] <- c(vacc_rate_4, vacc_rate_5, vacc_rate_6)       # Set new vaccination rates
  
  # Reset severity parameters so matrix has same number of rows as SIMTIME_fit
  
  par[["Proportion symptomatic"]] <- prop_symp <- outer(inf_symp_rate_base_case[1:SIMTIME_fit], inf_symp_age_scaling, "*")
  
  par[["Proportion hospitalised"]] <- prop_hosp <- outer(inf_hosp_rate_base_case[1:SIMTIME_fit], age_scaling, "*")
  
  par[["Proportion of hospital mortalities"]] <- prop_hosp_mortalities <- 
    outer(hosp_mort_rate_base_case[1:SIMTIME_fit], age_scaling * c(rep(par_fit[10], 6), 
                                                                   rep(par_fit[11], 12),
                                                                   rep(par_fit[12], 47), 
                                                                   rep(par_fit[13], 20), 
                                                                   rep(1, 16)), "*")
  
  # Run the model again,  beginning from 31/12/2022, for the year 2023
  
  out_fit <- as.data.frame(deSolve::ode(y = input_fit, times = times_fit, func = covid_model_odes, parms = par, method = "euler"))
  
  # Label the output columns for each compartment
  
  S_fit <- out_fit[, 2:(n_age_groups + 1)]
  E_fit <- out_fit[, (n_age_groups + 2):(2 * n_age_groups + 1)]
  I_fit <- out_fit[, (2 * n_age_groups + 2):(3 * n_age_groups + 1)]
  A_fit <- out_fit[, (3 * n_age_groups + 2):(4 * n_age_groups + 1)]
  H_fit <- out_fit[, (4 * n_age_groups + 2):(5 * n_age_groups + 1)]
  R_fit <- out_fit[, (5 * n_age_groups + 2):(6 * n_age_groups + 1)]
  D_fit <- out_fit[, (6 * n_age_groups + 2):(7 * n_age_groups + 1)]
  V_fit <- out_fit[, (7 * n_age_groups + 2):(8 * n_age_groups + 1)]
  Ev_fit <- out_fit[, (8 * n_age_groups + 2):(9 * n_age_groups + 1)]
  Iv_fit <- out_fit[, (9 * n_age_groups + 2):(10 * n_age_groups + 1)]
  Av_fit <- out_fit[, (10 * n_age_groups + 2):(11 * n_age_groups + 1)]
  Hv_fit <- out_fit[, (11 * n_age_groups + 2):(12 * n_age_groups + 1)]
  Rv_fit <- out_fit[, (12 * n_age_groups + 2):(13 * n_age_groups + 1)]
  W_fit <- out_fit[, (13 * n_age_groups + 2):(14 * n_age_groups + 1)]
  Ew_fit <- out_fit[, (14 * n_age_groups + 2):(15 * n_age_groups + 1)]
  Iw_fit <- out_fit[, (15 * n_age_groups + 2):(16 * n_age_groups + 1)]
  Aw_fit <- out_fit[, (16 * n_age_groups + 2):(17 * n_age_groups + 1)]
  Hw_fit <- out_fit[, (17 * n_age_groups + 2):(18 * n_age_groups + 1)]
  Rw_fit <- out_fit[, (18 * n_age_groups + 2):(19 * n_age_groups + 1)]
  II_fit <- out_fit[, (19 * n_age_groups + 2):(20 * n_age_groups + 1)]
  HH_fit <- out_fit[, (20 * n_age_groups + 2):(21 * n_age_groups + 1)]
  VV_fit <- out_fit[, (21 * n_age_groups + 2):(22 * n_age_groups + 1)]
  VP_fit <- out_fit[, (22 * n_age_groups + 2):(23 * n_age_groups + 1)]
  VB_fit <- out_fit[, (23 * n_age_groups + 2):(24 * n_age_groups + 1)]
  IIS_fit <- out_fit[, (24 * n_age_groups + 2):(25 * n_age_groups + 1)]
  RR_fit <- out_fit[, (25 * n_age_groups + 2):(26 * n_age_groups + 1)]
  RRV_fit <- out_fit[, (26 * n_age_groups + 2):(27 * n_age_groups + 1)]
  RRW_fit <- out_fit[, (27 * n_age_groups + 2):(28 * n_age_groups + 1)]
  
  # Check for negative values in model fitting period output
  
  if(is.null(fixr::check_for_negative_values(out_fit)) == FALSE) {
    stop("Negative values in out_fit")
  }
  
  # Run the model again with the previous output as the input, only if SIMTIME_future > 0
  
  if (SIMTIME_future > 0) {
    
    input <- as.numeric((cbind(S_fit, E_fit, I_fit, A_fit, H_fit, R_fit, D_fit, V_fit, Ev_fit, Iv_fit, Av_fit, Hv_fit, Rv_fit, W_fit, Ew_fit, Iw_fit,
                               Aw_fit, Hw_fit, Rw_fit, II_fit, HH_fit, VV_fit, VP_fit, VB_fit, IIS_fit, RR_fit, RRV_fit, RRW_fit))[SIMTIME_fit, ])
    
    # Calculating parameters for base case/no variant scenario
    
    variant_parameters <- calculate_variant_parameters(reinf_init, red_I_init, VE_S_init, VE_H_init, VE_HM_init, 
                                                       VE_T_init, wan_I_init, wan_VE_S_init, wan_VE_H_init, 
                                                       wan_VE_HM_init, wan_VE_T_init, beta_init, 
                                                       protection_factor, beta_factor, n_yearly_variants)
    
    names(variant_parameters) <- c("reinf_variants", "red_I_variants", "red_S_variants", "red_H_variants", 
                                   "red_HM_variants", "red_T_variants", "wan_I_variants", "wan_S_variants", 
                                   "wan_H_variants", "wan_HM_variants", "wan_T_variants", "beta_variants")
    
    # Setting base case variant parameters with new variants once every "variant_every_n_years"
    
    reinf_variants <- rep(reinf_init, SIMTIME_future)
    red_I_variants <- rep(red_I_init, SIMTIME_future)
    red_S_variants <- rep(red_S_init, SIMTIME_future)
    red_H_variants <- rep(red_H_init, SIMTIME_future)
    red_HM_variants <- rep(red_HM_init, SIMTIME_future)
    red_T_variants <- rep(red_T_init, SIMTIME_future)
    wan_I_variants <- rep(wan_I_init, SIMTIME_future)
    wan_S_variants <- rep(wan_S_init, SIMTIME_future)
    wan_H_variants <- rep(wan_H_init, SIMTIME_future)
    wan_HM_variants <- rep(wan_HM_init, SIMTIME_future)
    wan_T_variants <- rep(wan_T_init, SIMTIME_future)
    beta_variants <- rep(beta_init, SIMTIME_future)
    
    nth_element <- ((seq(1:floor(n_years_future / variant_every_n_years)) * variant_every_n_years) - 1) * 365
    
    for (j in 1:length(nth_element)) {
      reinf_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$reinf_variants
      red_I_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_I_variants
      red_S_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_S_variants
      red_H_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_H_variants
      red_HM_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_HM_variants
      red_T_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_T_variants
      wan_I_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_I_variants
      wan_S_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_S_variants
      wan_H_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_H_variants
      wan_HM_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_HM_variants
      wan_T_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_T_variants
      beta_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$beta_variants
    }
    
    # Setting variant parameters for scenario 2
    
    if (variant_scenario == "Scenario 2") {
      
      fifth_element <- ((seq(1:floor(n_years_future / 5)) * 5) - 1) * 365
      
      # Setting parameters in case of variant scenario 2 every five years
      
      variant_parameters_scenario_2 <- calculate_variant_parameters(reinf_init, red_I_init, VE_S_init, VE_H_init, 
                                                                    VE_HM_init, VE_T_init, wan_I_init, 
                                                                    wan_VE_S_init, wan_VE_H_init, wan_VE_HM_init, 
                                                                    wan_VE_T_init, beta_init, 
                                                                    protection_factor_scenario_2, beta_factor_scenario_2, 
                                                                    n_yearly_variants = 1)
      
      names(variant_parameters_scenario_2) <- c("reinf_variants_scenario_2", "red_I_variants_scenario_2", 
                                                "red_S_variants_scenario_2", "red_H_variants_scenario_2", 
                                                "red_HM_variants_scenario_2", "red_T_variants_scenario_2", 
                                                "wan_I_variants_scenario_2", "wan_S_variants_scenario_2", 
                                                "wan_H_variants_scenario_2", "wan_HM_variants_scenario_2", 
                                                "wan_T_variants_scenario_2", "beta_variants_scenario_2")
      
      # Replace the new variant of every fifth year with the scenario 2 variant (only one variant in every fifth year)
      
      for (j in 1:length(fifth_element)) {
        reinf_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$reinf_variants_scenario_2
        red_I_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_I_variants_scenario_2
        red_S_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_S_variants_scenario_2
        red_H_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_H_variants_scenario_2
        red_HM_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_HM_variants_scenario_2
        red_T_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_T_variants_scenario_2
        wan_I_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_I_variants_scenario_2
        wan_S_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_S_variants_scenario_2
        wan_H_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_H_variants_scenario_2
        wan_HM_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_HM_variants_scenario_2
        wan_T_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_T_variants_scenario_2
        beta_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$beta_variants_scenario_2
      }
    }
    
    # Setting variant parameters
    
    reinf <- par[["Protection against reinfection"]] <- reinf_variants
    red_I <- par[["Reduction in infectiousness"]] <- red_I_variants
    red_S <- par[["Reduction in symptomatic infections"]] <- red_S_variants                  # Set vaccine effectiveness to change according to new variants
    red_H <- par[["Reduction in hospitalisations"]] <- red_H_variants                  # Set vaccine effectiveness to change according to new variants
    red_HM <- par[["Reduction in hospital mortalities"]] <- red_HM_variants                  # Set vaccine effectiveness to change according to new variants
    red_T <- par[["Reduction in transmission"]] <- red_T_variants                  # Set vaccine effectiveness to change according to new variants
    wan_I <- par[["Waned reduction in infectiousness"]] <- wan_I_variants
    wan_S <- par[["Waned reduction in symptomatic infections"]] <- wan_S_variants                  # Set vaccine effectiveness to change according to new variant
    wan_H <- par[["Waned reduction in hospitalisations"]] <- wan_H_variants                  # Set vaccine effectiveness to change according to new variant
    wan_HM <- par[["Waned reduction in hospital mortalities"]] <- wan_HM_variants                  # Set vaccine effectiveness to change according to new variants
    wan_T <- par[["Waned reduction in transmission"]] <- wan_T_variants                  # Set vaccine effectiveness to change according to new variants
    beta <- par[["Transmission rate"]] <- beta_variants
    
    # Setting severity parameters for variant scenario
    
    prop_symp <- outer(inf_symp_rate, inf_symp_age_scaling, "*")
    
    prop_hosp <- outer(inf_hosp_rate, age_scaling, "*")
    
    prop_hosp_mortalities <- outer(hosp_mort_rate, age_scaling * c(rep(par_fit["Hospital mortality age scaling 0-5", i], 6), 
                                                                   rep(par_fit["Hospital mortality age scaling 6-17", i], 12),
                                                                   rep(par_fit["Hospital mortality age scaling 18-64", i], 47), 
                                                                   rep(par_fit["Hospital mortality age scaling 65-84", i], 20), 
                                                                   rep(1, 16)), "*")  
    
    # If any proportions are greater than 1, set to 1
    
    par[["Proportion symptomatic"]] <- prop_symp <- ifelse(prop_symp > 1, 1, prop_symp)
    
    par[["Proportion hospitalised"]] <- prop_hosp <- ifelse(prop_hosp > 1, 1, prop_hosp)
    
    par[["Proportion of hospital mortalities"]] <- prop_hosp_mortalities <-
      ifelse(prop_hosp_mortalities > 1, 1 , prop_hosp_mortalities)
    
    # Set vaccination scenario
    
    n_vacc_waves <- par[["Number of vaccination waves"]] <- n_vacc_waves_future
    
    vacc_rate <- par[["Vaccination rate"]] <- vacc_rate_future
    
    vacc_rate_shift <- par[["Vaccination rate shift"]] <- vacc_rate_shift_future                 # Vacc rate shift depending on vaccination scenario
    
    # Run the model again,  beginning from the year 2024
    
    out_future <- as.data.frame(deSolve::ode(y = input, times = times_future, func = covid_model_odes, parms = par, method = "euler"))
    
    # Merge 2023 output with the final model run output
    
    out <- rbind(out_fit, out_future)
    
    # Label all output columns
    
    S_future <- out_future[, 2:(n_age_groups + 1)]
    E_future <- out_future[, (n_age_groups + 2):(2 * n_age_groups + 1)]
    I_future <- out_future[, (2 * n_age_groups + 2):(3 * n_age_groups + 1)]
    A_future <- out_future[, (3 * n_age_groups + 2):(4 * n_age_groups + 1)]
    H_future <- out_future[, (4 * n_age_groups + 2):(5 * n_age_groups + 1)]
    R_future <- out_future[, (5 * n_age_groups + 2):(6 * n_age_groups + 1)]
    D_future <- out_future[, (6 * n_age_groups + 2):(7 * n_age_groups + 1)]
    V_future <- out_future[, (7 * n_age_groups + 2):(8 * n_age_groups + 1)]
    Ev_future <- out_future[, (8 * n_age_groups + 2):(9 * n_age_groups + 1)]
    Iv_future <- out_future[, (9 * n_age_groups + 2):(10 * n_age_groups + 1)]
    Av_future <- out_future[, (10 * n_age_groups + 2):(11 * n_age_groups + 1)]
    Hv_future <- out_future[, (11 * n_age_groups + 2):(12 * n_age_groups + 1)]
    Rv_future <- out_future[, (12 * n_age_groups + 2):(13 * n_age_groups + 1)]
    W_future <- out_future[, (13 * n_age_groups + 2):(14 * n_age_groups + 1)]
    Ew_future <- out_future[, (14 * n_age_groups + 2):(15 * n_age_groups + 1)]
    Iw_future <- out_future[, (15 * n_age_groups + 2):(16 * n_age_groups + 1)]
    Aw_future <- out_future[, (16 * n_age_groups + 2):(17 * n_age_groups + 1)]
    Hw_future <- out_future[, (17 * n_age_groups + 2):(18 * n_age_groups + 1)]
    Rw_future <- out_future[, (18 * n_age_groups + 2):(19 * n_age_groups + 1)]
    II_future <- out_future[, (19 * n_age_groups + 2):(20 * n_age_groups + 1)]
    HH_future <- out_future[, (20 * n_age_groups + 2):(21 * n_age_groups + 1)]
    VV_future <- out_future[, (21 * n_age_groups + 2):(22 * n_age_groups + 1)]
    VP_future <- out_future[, (22 * n_age_groups + 2):(23 * n_age_groups + 1)]
    VB_future <- out_future[, (23 * n_age_groups + 2):(24 * n_age_groups + 1)]
    IIS_future <- out_future[, (24 * n_age_groups + 2):(25 * n_age_groups + 1)]
    RR_future <- out_future[, (25 * n_age_groups + 2):(26 * n_age_groups + 1)]
    RRV_future <- out_future[, (26 * n_age_groups + 2):(27 * n_age_groups + 1)]
    RRW_future <- out_future[, (27 * n_age_groups + 2):(28 * n_age_groups + 1)]
    
    # Check for negative values in model future output
    
    if(is.null(fixr::check_for_negative_values(out_future)) == FALSE) {
      stop("Negative values in out_future")
    }
    
    # Calculating number of symptomatic recoveries over the past year for each time point in SIMTIME_future from the future and model fit runs
    
    RR_fit_future <- rbind(RR_fit, RR_future)
    RRV_fit_future <- rbind(RRV_fit, RRV_future)
    RRW_fit_future <- rbind(RRW_fit, RRW_future)
    
    inc_rec_unvacc_fit_future <- apply(RR_fit_future, 2, diff)           # Daily number of unvaccinated recoveries over the fitting and future period
    inc_rec_vacc_fit_future <- apply(RRV_fit_future, 2, diff)     # Daily number of vaccinated recoveries over the fitting and future period
    inc_rec_wan_fit_future <- apply(RRW_fit_future, 2, diff)     # Daily number of waned recoveries over the fitting and future period
    
    cumulative_n_recovered_unvacc_past_year <- apply(inc_rec_unvacc_fit_future, 2, cumsum)
    cumulative_n_recovered_unvacc_past_year_extended <- rbind(rep(0, n_age_groups), cumulative_n_recovered_unvacc_past_year)
    n_recovered_unvacc_past_year <- cumulative_n_recovered_unvacc_past_year[365:nrow(inc_rec_unvacc_fit_future), ] - 
      cumulative_n_recovered_unvacc_past_year_extended[1:(nrow(inc_rec_unvacc_fit_future) - 364), ]
    
    cumulative_n_recovered_vacc_past_year <- apply(inc_rec_vacc_fit_future, 2, cumsum)
    cumulative_n_recovered_vacc_past_year_extended <- rbind(rep(0, n_age_groups), cumulative_n_recovered_vacc_past_year)
    n_recovered_vacc_past_year <- cumulative_n_recovered_vacc_past_year[365:nrow(inc_rec_vacc_fit_future), ] - 
      cumulative_n_recovered_vacc_past_year_extended[1:(nrow(inc_rec_vacc_fit_future) - 364), ]
    
    cumulative_n_recovered_wan_past_year <- apply(inc_rec_wan_fit_future, 2, cumsum)
    cumulative_n_recovered_wan_past_year_extended <- rbind(rep(0, n_age_groups), cumulative_n_recovered_wan_past_year)
    n_recovered_wan_past_year <- cumulative_n_recovered_wan_past_year[365:nrow(inc_rec_wan_fit_future), ] - 
      cumulative_n_recovered_wan_past_year_extended[1:(nrow(inc_rec_wan_fit_future) - 364), ]
    
    n_recovered_past_year_yearly[[i]] <- rowsum(n_recovered_unvacc_past_year, rep(1:n_years_future, each = 365))
    n_recovered_vaccinated_past_year_yearly[[i]] <- rowsum(n_recovered_vacc_past_year, rep(1:n_years_future, each = 365))
    n_recovered_waned_past_year_yearly[[i]] <- rowsum(n_recovered_wan_past_year, rep(1:n_years_future, each = 365))
    
    # Calculating number in each compartment in each age group every year over the time horizon for calculating QALYs later
    
    n_healthy_yearly[[i]] <- rowsum((S_future + E_future + A_future + V_future + Ev_future + Av_future + 
                                       W_future + Ew_future + Aw_future + R_future +
                                       Rv_future + Rw_future), rep(1:n_years_future, each = 365))
    n_symptomatic_yearly[[i]] <- rowsum((I_future + Iv_future + Iw_future), rep(1:n_years_future, each = 365))
    n_hospitalised_yearly[[i]] <- rowsum((H_future + Hv_future + Hw_future), rep(1:n_years_future, each = 365))
    n_dead_yearly[[i]] <- rowsum((as.matrix(D_future) - matrix(as.numeric(D_fit[nrow(D_fit), ]), nrow = nrow(D_future), 
                                                               ncol = ncol(D_fit), byrow = TRUE)), 
                                 rep(1:n_years_future, each = 365))
    
    # Calculating incremental infections, hospitalisations, deaths vaccinations and symptomatic infections for the future model run
    
    inc_inf_future <- apply(II_future, 2, diff)
    inc_hosp_future <- apply(HH_future, 2, diff)
    inc_deaths_future <- apply(D_future, 2, diff)
    inc_vacc_future <- apply(VV_future, 2, diff)
    inc_symp_inf_future <- apply(IIS_future, 2, diff)
    inc_rec_future <- apply(RR_future + RRV_future + RRW_future, 2, diff)
    
    # # Calculating incremental infections, hospitalisations, deaths, vaccinations and symptomatic infections each year in each age group for calculating costs later
    # 
    # inc_inf_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_inf_future), rep(1:n_years_future, each = 365))
    # inc_hosp_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_hosp_future), rep(1:n_years_future, each = 365))
    # inc_deaths_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_deaths_future), rep(1:n_years_future, each = 365))
    # inc_vacc_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_vacc_future), rep(1:n_years_future, each = 365))
    # inc_symp_inf_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_symp_inf_future), rep(1:n_years_future, each = 365))
    # inc_rec_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_rec_future), rep(1:n_years_future, each = 365))
    # 
    # # Calculating total cases, hospital admissions, deaths and vaccinations given over the future simulation
    # 
    # total_inf_future[i] <- sum(II_future[SIMTIME_future, ] - II_fit[SIMTIME_fit, ])
    # total_vacc_future[i] <- sum(VV_future[SIMTIME_future, ] - VV_fit[SIMTIME_fit, ])
    # total_hosp_future[i] <- sum(HH_future[SIMTIME_future, ] - HH_fit[SIMTIME_fit, ])
    # total_deaths_future[i] <- sum(D_future[SIMTIME_future, ] - D_fit[SIMTIME_fit, ])
    # 
    # # Calculating total cases in each age group
    # 
    # total_inf_future_0_5[i] <- sum(II_future[SIMTIME_future, 0:6] - II_fit[SIMTIME_fit, 0:6])
    # total_inf_future_6_17[i] <- sum(II_future[SIMTIME_future, 7:18] - II_fit[SIMTIME_fit, 7:18])
    # total_inf_future_18_64[i] <- sum(II_future[SIMTIME_future, 19:65] - II_fit[SIMTIME_fit, 19:65])
    # total_inf_future_65_84[i] <- sum(II_future[SIMTIME_future, 66:85] - II_fit[SIMTIME_fit, 66:85])
    # total_inf_future_85_[i] <- sum(II_future[SIMTIME_future, 86:101] - II_fit[SIMTIME_fit, 86:101])
    # 
    # # Calculating total hospitalisations in each age group
    # 
    # total_hosp_future_0_5[i] <- sum(HH_future[SIMTIME_future, 0:6] - HH_fit[SIMTIME_fit, 0:6])
    # total_hosp_future_6_17[i] <- sum(HH_future[SIMTIME_future, 7:18] - HH_fit[SIMTIME_fit, 7:18])
    # total_hosp_future_18_64[i] <- sum(HH_future[SIMTIME_future, 19:65] - HH_fit[SIMTIME_fit, 19:65])
    # total_hosp_future_65_84[i] <- sum(HH_future[SIMTIME_future, 66:85] - HH_fit[SIMTIME_fit, 66:85])
    # total_hosp_future_85_[i] <- sum(HH_future[SIMTIME_future, 86:101] - HH_fit[SIMTIME_fit, 86:101])
    # 
    # # Calculating total deaths in each age group
    # 
    # total_deaths_future_0_5[i] <- sum(D_future[SIMTIME_future, 0:6] - D_fit[SIMTIME_fit, 0:6])
    # total_deaths_future_6_17[i] <- sum(D_future[SIMTIME_future, 7:18] - D_fit[SIMTIME_fit, 7:18])
    # total_deaths_future_18_64[i] <- sum(D_future[SIMTIME_future, 19:65] - D_fit[SIMTIME_fit, 19:65])
    # total_deaths_future_65_84[i] <- sum(D_future[SIMTIME_future, 66:85] - D_fit[SIMTIME_fit, 66:85])
    # total_deaths_future_85_[i] <- sum(D_future[SIMTIME_future, 86:101] - D_fit[SIMTIME_fit, 86:101])
    # 
    # # Calculating total ICU admissions in each age group
    # 
    # total_icu_future_0_5[i] <- sum(icu_rate[1:6] * (HH_future[SIMTIME_future, 1:6] - HH_fit[SIMTIME_fit, 1:6]))
    # total_icu_future_6_17[i] <- sum(icu_rate[7:18] * (HH_future[SIMTIME_future, 7:18] - HH_fit[SIMTIME_fit, 7:18]))
    # total_icu_future_18_64[i] <- sum(icu_rate[19:65] * (HH_future[SIMTIME_future, 19:65] - HH_fit[SIMTIME_fit, 19:65]))
    # total_icu_future_65_84[i] <- sum(icu_rate[66:85] * (HH_future[SIMTIME_future, 66:85] - HH_fit[SIMTIME_fit, 66:85]))
    # total_icu_future_85_[i] <- sum(icu_rate[86:101] * (HH_future[SIMTIME_future, 86:101] - HH_fit[SIMTIME_fit, 86:101]))
    
  } else {
    
    out <-  out_fit
    
  }
  
  # Label the output column for those in hospital and those that have died, separating into age groups
  
  HH <- out[(20 * n_age_groups + 2):(21 * n_age_groups + 1)]
  HH_0_5 <- out[, (20 * n_age_groups + 2):(20 * n_age_groups + 7)]
  HH_6_17 <- out[, (20 * n_age_groups + 8):(20 * n_age_groups + 19)]
  HH_18_64 <- out[, (20 * n_age_groups + 20):(20 * n_age_groups + 66)]
  HH_65_84 <- out[, (20 * n_age_groups + 67):(20 * n_age_groups + 86)]
  HH_85_ <- out[, (20 * n_age_groups + 87):(21 * n_age_groups + 1)]
  
  DD <- out[(6 * n_age_groups + 2):(7 * n_age_groups + 1)]
  DD_0_5 <- out[, (6 * n_age_groups + 2):(6 * n_age_groups + 7)]
  DD_6_17 <- out[, (6 * n_age_groups + 8):(6 * n_age_groups + 19)]
  DD_18_64 <- out[, (6 * n_age_groups + 20):(6 * n_age_groups + 66)]
  DD_65_84 <- out[, (6 * n_age_groups + 67):(6 * n_age_groups + 86)]
  DD_85_ <- out[, (6 * n_age_groups + 87):(7 * n_age_groups + 1)]
  
  names(HH_0_5) <- names(DD_0_5) <- seq(0, 5)
  names(HH_6_17) <- names(DD_6_17) <- seq(6, 17)
  names(HH_18_64) <- names(DD_18_64) <- seq(18, 64)
  names(HH_65_84) <- names(DD_65_84) <- seq(65, 84)
  names(HH_85_) <- names(DD_85_) <- seq(85, 100)
  
  #   # Calculating cumulative hospital admissions and deaths for all age groups and incident hospital admissions and deaths per age group
  #   
  #   cum_hosp[, i] <- rowSums(HH)
  #   inc_hosp[, i] <- diff(rowSums(HH))
  #   
  #   inc_hosp_0_5[, i] <- diff(rowSums(HH_0_5))
  #   inc_hosp_6_17[, i] <- diff(rowSums(HH_6_17))
  #   inc_hosp_18_64[, i] <- diff(rowSums(HH_18_64))
  #   inc_hosp_65_84[, i] <- diff(rowSums(HH_65_84))
  #   inc_hosp_85_[, i] <- diff(rowSums(HH_85_))
  #   
  #   cum_deaths[, i] <- rowSums(DD)
  #   inc_deaths[, i] <- diff(rowSums(DD))
  #   
  #   inc_deaths_0_5[, i] <- diff(rowSums(DD_0_5))
  #   inc_deaths_6_17[, i] <- diff(rowSums(DD_6_17))
  #   inc_deaths_18_64[, i] <- diff(rowSums(DD_18_64))
  #   inc_deaths_65_84[, i] <- diff(rowSums(DD_65_84))
  #   inc_deaths_85_[, i] <- diff(rowSums(DD_85_))
  #   
  # }
  # 
  # print(mean(burn_in_deaths))
  
  # Calculating total number of hospital admissions per age group in 2023
  
  total_hosp <- rowSums(HH)[SIMTIME] - rowSums(HH)[1]
  total_hosp_0_5 <- rowSums(HH_0_5)[SIMTIME] - rowSums(HH_0_5)[1]
  total_hosp_6_17 <- rowSums(HH_6_17)[SIMTIME] - rowSums(HH_6_17)[1]
  total_hosp_18_64 <- rowSums(HH_18_64)[SIMTIME] - rowSums(HH_18_64)[1]
  total_hosp_65_84 <- rowSums(HH_65_84)[SIMTIME] - rowSums(HH_65_84)[1]
  total_hosp_85_ <- rowSums(HH_85_)[SIMTIME] - rowSums(HH_85_)[1]
  
  # Calculating total deaths per age group in 2023
  
  total_deaths <- rowSums(DD)[SIMTIME] - rowSums(DD)[1]
  total_deaths_0_5 <- rowSums(DD_0_5)[SIMTIME] - rowSums(DD_0_5)[1]
  total_deaths_6_17 <- rowSums(DD_6_17)[SIMTIME] - rowSums(DD_6_17)[1]
  total_deaths_18_64 <- rowSums(DD_18_64)[SIMTIME] - rowSums(DD_18_64)[1]
  total_deaths_65_84 <- rowSums(DD_65_84)[SIMTIME] - rowSums(DD_65_84)[1]
  total_deaths_85_ <- rowSums(DD_85_)[SIMTIME] - rowSums(DD_85_)[1]
  
  total_hosp_deaths_all_ages <- c(total_hosp_0_5, total_hosp_6_17, total_hosp_18_64, total_hosp_65_84, total_hosp_85_, 
                                  total_deaths_0_5, total_deaths_6_17, total_deaths_18_64, total_deaths_65_84, total_deaths_85_)
  
  names(total_hosp_deaths_all_ages) <- c("Hosp 0-5", "Hosp 6-17", "Hosp 18-64", "Hosp 65-84", "Hosp 85+", 
                                         "Deaths 0-5", "Deaths 6-17", "Deaths 18-64", "Deaths 65-84", "Deaths 85+")
  
  return(total_hosp_deaths_all_ages)
  
}

# Model run function using a sample of parameters from the posterior distribution of the model fit

covid_model_run_fixed_params_age_dist <- function(par_fit) {
  
  # Create function "covid_model_odes" that computes the values of the model ODEs.
  
  covid_model_odes <- function(t, states, par){
    with(as.list(c(states, par)), {
      
      # Parameters
      
      ageing <- par$Ageing
      lat_period <- par$`Latent period`
      prop_symp <- par$`Proportion symptomatic`
      symp_period <- par$`Symptomatic period`
      rec_rate <- par$`Recovery rate (asymptomatic)`
      prop_hosp <- par$`Proportion hospitalised`
      hosp_period <- par$`Hospital length of stay`
      prop_hosp_mortalities <- par$`Proportion of hospital mortalities`
      vacc_rate <- par$`Vaccination rate`
      vacc_rate_shift <- par$`Vaccination rate shift`
      n_vacc_waves <- par$`Number of vaccination waves`         # No vaccination for under 12s, two vaccinations a year for 12-100 year olds
      vacc_uptake_rate <- par$`Vaccine uptake rate`
      wane_rate <- par$`Wane rate`                      # Immunity wanes after 270 days (nine months)
      red_I <- par$`Reduction in infectiousness`
      red_S <- par$`Reduction in symptomatic infections`
      red_H <- par$`Reduction in hospitalisations`
      red_HM <- par$`Reduction in hospital mortalities`
      wan_I <- par$`Waned reduction in infectiousness`
      wan_S <- par$`Waned reduction in symptomatic infections`
      wan_H <- par$`Waned reduction in hospitalisations`
      wan_HM <- par$`Waned reduction in hospital mortalities`
      reinf <- par$`Protection against reinfection`
      red_T <- par$`Reduction in transmission`
      wan_T <- par$`Waned reduction in transmission`
      tau <- par$`Reduction in asymptomatic infectiousness`
      beta <- par$`Transmission rate`
      num_peaks <- par$`Number of annual peaks in transmission`
      inf_fluctuation <- par$`Amount of fluctuation in infectiousness`
      inf_shift <- par$`Infection shift`
      n_contacts <- par$`Number of contacts per age group`
      total_pop <- par$`Total population size`
      n_age_groups <- par$`Number of age groups`
      all_cause_death_rate <- par$`All-cause death rate`
      birth_rate <- par$`Birth rate`
      
      # States
      
      S = states[1:n_age_groups]                                        # Susceptible
      E = states[(n_age_groups + 1):(2 * n_age_groups)]                     # Exposed
      I = states[(2 * n_age_groups + 1):(3 * n_age_groups)]                   # Infected, symptomatic
      A = states[(3 * n_age_groups + 1):(4 * n_age_groups)]                   # Infected, asymptomatic
      H = states[(4 * n_age_groups + 1):(5 * n_age_groups)]                   # Hospitalised
      R = states[(5 * n_age_groups + 1):(6 * n_age_groups)]                   # Recovered
      D = states[(6 * n_age_groups + 1):(7 * n_age_groups)]                   # Death
      V = states[(7 * n_age_groups + 1):(8 * n_age_groups)]                   # Vaccinated
      Ev = states[(8 * n_age_groups + 1):(9 * n_age_groups)]                  # Exposed-vaccinated
      Iv = states[(9 * n_age_groups + 1):(10 * n_age_groups)]                 # Infected, symptomatic-vaccinated
      Av = states[(10 * n_age_groups + 1):(11 * n_age_groups)]                # Infected, asymptomatic-vaccinated
      Hv = states[(11 * n_age_groups + 1):(12 * n_age_groups)]                # Hospitalised-vaccinated
      Rv = states[(12 * n_age_groups + 1):(13 * n_age_groups)]                # Recovered-vaccinated
      W = states[(13 * n_age_groups + 1):(14 * n_age_groups)]                 # Waned immunity
      Ew = states[(14 * n_age_groups + 1):(15 * n_age_groups)]                # Exposed-waned
      Iw = states[(15 * n_age_groups + 1):(16 * n_age_groups)]                # Infected, symptomatic-waned
      Aw = states[(16 * n_age_groups + 1):(17 * n_age_groups)]                # Infected, asymptomatic-waned
      Hw = states[(17 * n_age_groups + 1):(18 * n_age_groups)]                # Hospitalised-waned
      Rw = states[(18 * n_age_groups + 1):(19 * n_age_groups)]                # Recovered-waned
      II = states[(19 * n_age_groups + 1):(20 * n_age_groups)]                # Combined, cumulative infected states
      HH = states[(20 * n_age_groups + 1):(21 * n_age_groups)]                # Combined, cumulative hospitalised states
      VV = states[(21 * n_age_groups + 1):(22 * n_age_groups)]                # Combined vaccinated states
      VP = states[(22 * n_age_groups + 1):(23 * n_age_groups)]                # Combined vaccinated (primary programme) states
      VB = states[(23 * n_age_groups + 1):(24 * n_age_groups)]                # Combined vaccinated (booster programme) states
      IIS = states[(24 * n_age_groups + 1):(25 * n_age_groups)]               # Combined, cumulative symptomatic infected states
      RR = states[(25 * n_age_groups + 1):(26 * n_age_groups)]                # Cumulative symptomatic recovered unvaccinated state
      RRV = states[(26 * n_age_groups + 1):(27 * n_age_groups)]               # Cumulative symptomatic recovered vaccinated state
      RRW = states[(27 * n_age_groups + 1):(28 * n_age_groups)]               # Cumulative symptomatic recovered waned state
      
      # Specify total population
      
      N <- S + E + I + A + H + R + V + Ev + Iv + Av + Hv + Rv + W + Ew + Iw + Aw + Hw + Rw
      
      # All age compartments shifted up one for ageing
      
      older_S <- c(birth_rate, S)[1:length(S)]
      older_E <- c(0, E)[1:length(E)]
      older_I <- c(0, I)[1:length(I)]
      older_A<- c(0, A)[1:length(A)]
      older_H <- c(0, H)[1:length(H)]
      older_R <- c(0, R)[1:length(R)]
      older_V <- c(0, V)[1:length(V)]
      older_Ev <- c(0, Ev)[1:length(Ev)]
      older_Iv <- c(0, Iv)[1:length(Iv)]
      older_Av <- c(0, Av)[1:length(Av)]
      older_Hv <- c(0, Hv)[1:length(Hv)]
      older_Rv <- c(0, Rv)[1:length(Rv)]
      older_W <- c(0, W)[1:length(W)]
      older_Ew <- c(0, Ew)[1:length(Ew)]
      older_Iw <- c(0, Iw)[1:length(Iw)]
      older_Aw <- c(0, Aw)[1:length(Aw)]
      older_Hw <- c(0, Hw)[1:length(Hw)]
      older_Rw <- c(0, Rw)[1:length(Rw)]
      
      # Equation for frequency-dependent force of infection
      
      lambda <- (beta[t] + inf_fluctuation * (sin(inf_shift + num_peaks * 2 * pi * t / 365))) * 
        n_contacts * ((sum(I) + red_T[t] * sum(Iv) + wan_T[t] * sum(Iw)) + tau * (sum(A) + red_T[t] * sum(Av)
                                                                                  + wan_T[t] * sum(Aw))) / sum(N)
      
      # Overall time-dependent vaccination rate
      
      overall_vacc_rate <- (vacc_rate + (vacc_rate * sin(vacc_rate_shift + n_vacc_waves * 2 * pi * t / 365))) * vacc_uptake_rate
      
      # ODEs for each model state
      
      dD <- prop_hosp_mortalities[t, ] * (1 / hosp_period) * H + red_HM[t] * prop_hosp_mortalities[t, ] * (1 / hosp_period) * Hv + 
        wan_HM[t] * prop_hosp_mortalities[t, ] * (1 / hosp_period) * Hw
      
      dS <- ageing * older_S - overall_vacc_rate * S - lambda * S - ageing * S - all_cause_death_rate * S
      
      dE <- lambda * S - overall_vacc_rate * E - (1 / lat_period) * E - ageing * E + ageing * older_E - all_cause_death_rate * E
      
      dI <- prop_symp[t, ] * (1 / lat_period) * E - (1 / symp_period) * I  - overall_vacc_rate * I - ageing * I + ageing * older_I - 
        all_cause_death_rate * I
      
      dA <- (1 - prop_symp[t, ]) * (1 / lat_period) * E - rec_rate * A - overall_vacc_rate * A - ageing * A + ageing * older_A - 
        all_cause_death_rate * A
      
      dH <- prop_hosp[t, ] * (1 / symp_period) * I - (1 / hosp_period) * H - overall_vacc_rate * H - ageing * H + ageing * older_H - 
        all_cause_death_rate * H
      
      dR <- rec_rate * A + (1 - prop_hosp_mortalities[t, ]) * (1 / hosp_period) * H + (1 - prop_hosp[t, ]) * (1 / symp_period) * I - 
        wane_rate * R - overall_vacc_rate * R - red_I[t] * lambda * R - ageing * R + ageing * older_R - all_cause_death_rate * R
      
      # Vaccinated compartments
      
      dV <- overall_vacc_rate * S + overall_vacc_rate * W - wane_rate * V - red_I[t] * lambda * V - ageing * V + ageing * older_V - 
        all_cause_death_rate * V - overall_vacc_rate * V + overall_vacc_rate * V
      
      dEv <- red_I[t] * lambda * V + overall_vacc_rate * E + overall_vacc_rate * Ew + red_I[t] * lambda * Rv + red_I[t] * lambda * R + 
        wan_I[t] * lambda * Rw - wane_rate * Ev - (1 / lat_period) * Ev - ageing * Ev + ageing * older_Ev - all_cause_death_rate * Ev -
        overall_vacc_rate * Ev + overall_vacc_rate * Ev
      
      dIv <- red_S[t] * prop_symp[t, ] * (1 / lat_period) * Ev - (1 / symp_period) * Iv + overall_vacc_rate * I + overall_vacc_rate * Iw -
        ageing * Iv + ageing * older_Iv - all_cause_death_rate * Iv - overall_vacc_rate * Iv + overall_vacc_rate * Iv
      
      dAv <- (1 - red_S[t] * prop_symp[t, ]) * (1 / lat_period) * Ev - rec_rate * Av + overall_vacc_rate * A + overall_vacc_rate * Aw -
        ageing * Av + ageing * older_Av - all_cause_death_rate * Av - overall_vacc_rate * Av + overall_vacc_rate * Av
      
      dHv <- red_H[t] * prop_hosp[t, ] * (1 / symp_period) * Iv - (1 / hosp_period) * Hv + overall_vacc_rate * H + overall_vacc_rate * Hw  - 
        ageing * Hv + ageing * older_Hv - all_cause_death_rate * Hv - overall_vacc_rate * Hv + overall_vacc_rate * Hv
      
      dRv <- rec_rate * Av + (1 - red_HM[t] * prop_hosp_mortalities[t, ]) * (1 / hosp_period) * Hv + 
        (1 - (red_H[t] * prop_hosp[t, ])) * (1 / symp_period) * Iv + rec_rate * Aw + 
        (1 - wan_HM[t] * prop_hosp_mortalities[t, ]) * (1 / hosp_period) * Hw + 
        (1 - (wan_H[t] * prop_hosp[t, ])) * (1 / symp_period) * Iw + overall_vacc_rate * R + overall_vacc_rate * Rw -
        wane_rate * Rv - red_I[t] * lambda * Rv - ageing * Rv + ageing * older_Rv - all_cause_death_rate * Rv -
        overall_vacc_rate * Rv + overall_vacc_rate * Rv
      
      # Waned compartments
      
      dW <- wane_rate * V - overall_vacc_rate * W - wan_I[t] * lambda * W - ageing * W + ageing * older_W - all_cause_death_rate * W
      
      dEw <- wan_I[t] * lambda * W + wane_rate * Ev - 
        overall_vacc_rate * Ew - (1 / lat_period) * Ew - ageing * Ew + ageing * older_Ew - all_cause_death_rate * Ew
      
      dIw <- wan_S[t] * prop_symp[t, ] * (1 / lat_period) * Ew - (1 / symp_period) * Iw - overall_vacc_rate * Iw - ageing * Iw + 
        ageing * older_Iw - all_cause_death_rate * Iw
      
      dAw <- (1 - wan_S[t] * prop_symp[t, ]) * (1 / lat_period) * Ew - rec_rate * Aw - overall_vacc_rate * Aw - ageing * Aw + 
        ageing * older_Aw - all_cause_death_rate * Aw
      
      dHw <- wan_H[t] * prop_hosp[t, ] * (1 / symp_period) * Iw - (1 / hosp_period) * Hw - overall_vacc_rate * Hw - ageing * Hw + 
        ageing * older_Hw - all_cause_death_rate * Hw
      
      dRw <- wane_rate * R + wane_rate * Rv - overall_vacc_rate * Rw - wan_I[t] * lambda * Rw - ageing * Rw + ageing * older_Rw - 
        all_cause_death_rate * Rw
      
      # Combined infected, hospitalised, vaccinated, infected symptomatic and recovered states
      
      dII <- (1 / lat_period) * (E + Ev + Ew)
      
      dHH <- prop_hosp[t, ] * (1 / symp_period) * (I + red_H[t] * Iv + wan_H[t] * Iw)
      
      dVV <- overall_vacc_rate * (S + E + I + A + H + R + V + Ev + Iv + Av + Hv + Rv + W + Ew + Iw + Aw + Hw + Rw)
      
      dVP <- overall_vacc_rate * (S + E + I + A + H + R)
      
      dVB <- overall_vacc_rate * (V + Ev + Iv + Av + Hv + Rv + W + Ew + Iw + Aw + Hw + Rw)
      
      dIIS <- prop_symp[t, ] * (1 / lat_period) * (E + red_S[t] * Ev + wan_S[t] * Ew)
      
      dRR <- (1 - prop_hosp_mortalities[t, ]) * (1 / hosp_period) * H + (1 - prop_hosp[t, ]) * (1 / symp_period) * I
      
      dRRV <- (1 - red_HM[t] * prop_hosp_mortalities[t, ]) * (1 / hosp_period) * Hv + 
        (1 - (red_H[t] * prop_hosp[t, ])) * (1 / symp_period) * Iv
      
      dRRW <- (1 - wan_HM[t] * prop_hosp_mortalities[t, ]) * (1 / hosp_period) * Hw + 
        (1 - (wan_H[t] * prop_hosp[t, ])) * (1 / symp_period) * Iw
      
      
      return(list(c(dS, dE, dI, dA, dH, dR, dD, dV, dEv, dIv, dAv, dHv, dRv, dW, dEw, dIw, dAw, dHw, dRw, 
                    dII, dHH, dVV, dVP, dVB, dIIS, dRR, dRRV, dRRW)))	
    }
    )
  }
  
  
  # Function to calculate the new vaccine effectiveness, taking in the initial vaccine efficacy and the factor change
  
  calculate_new_VE <- function(VE_init, factor_change) {
    new_VE <- c(VE_init * factor_change, VE_init * (factor_change) ^ 2, VE_init * (factor_change) ^ 3, 
                VE_init * (factor_change) ^ 4)
  }
  
  # Function to calculate the new VE reduction parameters (Red_I) from the current reduction parameters and the desired VE
  
  calculate_reduction <- function(current_reduction, desired_VE){
    reduction <- (1 - desired_VE / 100) / current_reduction
  }
  
  # Function to calculate VE over the course of a year with varying number of new variants
  
  calculate_param_change <- function(n_yearly_variants, param_variant_change) {
    if (n_yearly_variants == 1) {
      protection <- c(rep(param_variant_change[1], ceiling((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))), 
                      seq(param_variant_change[1], param_variant_change[2], length.out = 30), 
                      rep(param_variant_change[2], floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
    }  else if (n_yearly_variants == 2) {
      protection <- c(rep(param_variant_change[1], ceiling((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))), 
                      seq(param_variant_change[1], param_variant_change[2], length.out = 30))
      protection <- append(protection, rep(param_variant_change[2],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[2], param_variant_change[3], length.out = 30))
      protection <- append(protection, rep(param_variant_change[3],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
    }  else if (n_yearly_variants == 3) {
      protection <- c(rep(param_variant_change[1], ceiling((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))), 
                      seq(param_variant_change[1], param_variant_change[2], length.out = 30))
      protection <- append(protection, rep(param_variant_change[2],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[2], param_variant_change[3], length.out = 30))
      protection <- append(protection, rep(param_variant_change[3],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[3], param_variant_change[4], length.out = 30))
      protection <- append(protection, rep(param_variant_change[4],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
    }  else if (n_yearly_variants == 4) {
      protection <- c(rep(param_variant_change[1], ceiling((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))), 
                      seq(param_variant_change[1], param_variant_change[2], length.out = 30))
      protection <- append(protection, rep(param_variant_change[2],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[2], param_variant_change[3], length.out = 30))
      protection <- append(protection, rep(param_variant_change[3],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[3], param_variant_change[4], length.out = 30))
      protection <- append(protection, rep(param_variant_change[4],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
      protection <- append(protection, seq(param_variant_change[4], param_variant_change[5], length.out = 30))
      protection <- append(protection, rep(param_variant_change[5],  floor((365 - 30 * n_yearly_variants) / (n_yearly_variants + 1))))
    } 
    protection <- ifelse(protection > 1, 1, protection)
    return (protection)
  }
  
  calculate_variant_parameters <- function(reinf_init, red_I_init, VE_S_init, VE_H_init, VE_HM_init, VE_T_init, 
                                           wan_I_init, wan_VE_S_init, wan_VE_H_init, wan_VE_HM_init, wan_VE_T_init, 
                                           beta_init, protection_factor, beta_factor, n_yearly_variants) {
    
    # Calculating vaccine effectiveness for the param change
    
    VE_reinf <- calculate_new_VE(100 * (1 - reinf_init), protection_factor)
    VE_I <- calculate_new_VE(100 * (1 - red_I_init), protection_factor)
    VE_S <- calculate_new_VE(VE_S_init, protection_factor)
    VE_H <- calculate_new_VE(VE_H_init, protection_factor)
    VE_HM <- calculate_new_VE(VE_HM_init, protection_factor)
    VE_T <- calculate_new_VE(VE_T_init, protection_factor)
    wan_VE_I <- calculate_new_VE(100 * (1 - wan_I_init), protection_factor)
    wan_VE_S <- calculate_new_VE(wan_VE_S_init, protection_factor)
    wan_VE_H <- calculate_new_VE(wan_VE_H_init, protection_factor)
    wan_VE_HM <- calculate_new_VE(wan_VE_HM_init, protection_factor)
    wan_VE_T <- calculate_new_VE(wan_VE_T_init, protection_factor)
    beta_variant_change <- calculate_new_VE(beta_init, beta_factor)
    
    # Calculating reduction parameters from the VE
    
    reinf_variant_change <- c(reinf_init, 1 - (VE_reinf / 100))
    red_I_variant_change <- c(red_I_init, 1 - (VE_I / 100))
    red_S_variant_change <- calculate_reduction(current_reduction = red_I_variant_change, desired_VE = c(VE_S_init, VE_S))
    red_H_variant_change <- calculate_reduction(current_reduction = 1 - (c(VE_S_init, VE_S) / 100), desired_VE = c(VE_H_init, VE_H))
    red_HM_variant_change <- calculate_reduction(current_reduction = 1 - (c(VE_H_init, VE_H) / 100), desired_VE = c(VE_HM_init, VE_HM))
    red_T_variant_change <- calculate_reduction(current_reduction = 1, desired_VE = c(VE_T_init, VE_T))
    wan_I_variant_change <- c(wan_I_init, 1 - (wan_VE_I / 100))
    wan_S_variant_change <- calculate_reduction(current_reduction = wan_I_variant_change, desired_VE = c(wan_VE_S_init, wan_VE_S))
    wan_H_variant_change <- calculate_reduction(current_reduction = 1 - (c(wan_VE_S_init, wan_VE_S) / 100), desired_VE = c(wan_VE_H_init, wan_VE_H))
    wan_HM_variant_change <- calculate_reduction(current_reduction = 1 - (c(wan_VE_H_init, wan_VE_H) / 100), desired_VE = c(wan_VE_HM_init, wan_VE_HM))
    wan_T_variant_change <- calculate_reduction(current_reduction = 1, desired_VE = c(wan_VE_T_init, wan_VE_T))
    beta_variant_change <- c(beta_init, beta_variant_change)
    
    # Calculating variant parameters
    
    reinf_variants <- calculate_param_change(n_yearly_variants, reinf_variant_change)
    red_I_variants <- calculate_param_change(n_yearly_variants, red_I_variant_change)
    red_S_variants <- calculate_param_change(n_yearly_variants, red_S_variant_change)
    red_H_variants <- calculate_param_change(n_yearly_variants, red_H_variant_change)
    red_HM_variants <- calculate_param_change(n_yearly_variants, red_HM_variant_change)
    red_T_variants <- calculate_param_change(n_yearly_variants, red_T_variant_change)
    wan_I_variants <- calculate_param_change(n_yearly_variants, wan_I_variant_change)
    wan_S_variants <- calculate_param_change(n_yearly_variants, wan_S_variant_change)
    wan_H_variants <- calculate_param_change(n_yearly_variants, wan_H_variant_change)
    wan_HM_variants <- calculate_param_change(n_yearly_variants, wan_HM_variant_change)
    wan_T_variants <- calculate_param_change(n_yearly_variants, wan_T_variant_change)
    beta_variants <- calculate_param_change(n_yearly_variants, beta_variant_change)
    
    return (list(reinf_variants, red_I_variants, red_S_variants, red_H_variants, red_HM_variants, red_T_variants, 
                 wan_I_variants, wan_S_variants, wan_H_variants, wan_HM_variants, wan_T_variants, beta_variants))
  }
  
  calculate_severity_parameters_scenario_3 <- function (severity_init, severity_factor_scenario_3) {
    severity_scenario_3 <- c(rep(severity_init, ceiling((365 - 30) / 2)), seq(severity_init, severity_init * severity_factor_scenario_3, 
                                                                              length.out = 30), rep(severity_init * severity_factor_scenario_3,  
                                                                                                    floor((365 - 30) / 2)))
    severity_scenario_3 <- ifelse(severity_scenario_3 > 1, 1, severity_scenario_3)
    return (severity_scenario_3)
  }
  
  
  # Set variant scenario
  
  variant_scenario <- "No variants"
  n_yearly_variants <- 1
  vaccination_scenario <- "B9"      # arbitrary choice as the future simlulation time is 0 years
  
  # Parameters
  
  # Sizes of each age group in the population of England, using the output of the population burn-in function for 500 years
  
  n_age_groups <- 101
  age_dist <- c(696785, 696646, 696507, 696367, 696228, 696159, 696089, 696019, 695950, 695880, 695811, 695741, 695671, 
                695602, 695532, 695324, 695115, 694907, 694698, 694490, 694212, 693935, 693657, 693380, 693103, 692826, 
                692549, 692272, 691995, 691718, 691234, 690751, 690268, 689785, 689302, 688614, 687926, 687238, 686552, 685866, 
                684839, 683813, 682789, 681766, 680745, 679183, 677624, 676069, 674518, 672970, 670690, 668417, 666152, 663895, 
                661645, 658354, 655078, 651819, 648576, 645350, 640482, 635651, 630856, 626098, 621376, 613947, 606607, 599355, 
                592189, 585109, 574031, 563162, 552498, 542037, 531774, 515535, 499791, 484528, 469732, 455387, 430301, 406596, 
                384197, 363033, 343034, 310185, 280482, 253623, 229337, 207376, 169134, 137945, 112507, 91760, 74839, 61038, 
                49782, 40602, 33115, 27008, 27008)
  total_pop <- total_original_pop <- sum(age_dist)
  age_dist_proportion <- age_dist / total_pop
  
  # On average there are 630000 deaths during the burn-in period. Add these to the initial population so that these are not lost
  # over the burn-in.
  
  age_dist <- age_dist + 630000 * age_dist_proportion
  age_dist <<- age_dist
  total_pop <- sum(age_dist)
  
  # Time horizon for burn-in period
  
  n_years_0 <- 25                  # Run the model initially for 25 years
  SIMTIME_0 <- n_years_0 * 365
  times_0 <- seq(1, SIMTIME_0, by = 1)
  
  # Time horizon for fitting period (2023)
  
  n_years_fit <- 1
  SIMTIME_fit <- n_years_fit * 365
  times_fit <- seq(1, SIMTIME_fit, by = 1)
  
  # Model run time horizon
  
  n_years_future <- 0
  SIMTIME_future <- n_years_future * 365
  if (SIMTIME_future > 0) {
    times_future <- seq(1, SIMTIME_future, by = 1)
  }
  
  SIMTIME <- SIMTIME_fit + SIMTIME_future
  times <- seq(1, SIMTIME, by = 1)
  
  # Vaccination rates for fitting period
  
  vacc_rate_4 <- rep(0, 65)                                          # Vaccination rate for 0 to 64 year olds (no vaccination)
  vacc_rate_5 <- rep(1 / 365, 10)                                      # Vaccination rate for 65 to 74 year olds
  vacc_rate_6 <- rep(2 / 365, 26)                                      # Vaccination rate for 75 to 100 year olds
  
  # Set vacc rates in burn-in period to same as fitting period
  
  vacc_rate_1 <- vacc_rate_4
  vacc_rate_2 <- vacc_rate_5
  vacc_rate_3 <- vacc_rate_6
  
  # Number of vaccination waves for burn-in and fitting period
  
  n_vacc_waves_fitting <- c(rep(0, 65), rep(1, 10), rep(2, 26))      # Two boosters for 75+, one booster for 65-74, no vaccination for under 64s
  
  n_vacc_waves_burn_in <- n_vacc_waves_fitting
  
  # Setting up the vaccination scenario
  
  vacc_rate_scenarios <- list(Comparator <- rep(0, 101),
                              A1 <- c(rep(0, 5), rep(1 / 365, 96)),
                              A2 <- c(rep(0, 18), rep(1 / 365, 83)),
                              A3 <- c(rep(0, 50), rep(1 / 365, 51)),
                              A4 <- c(rep(0, 55), rep(1 / 365, 46)),
                              A5 <- c(rep(0, 60), rep(1 / 365, 41)),
                              A6 <- c(rep(0, 65), rep(1 / 365, 36)),
                              A7 <- c(rep(0, 70), rep(1 / 365, 31)),
                              A8 <- c(rep(0, 75), rep(1 / 365, 26)),
                              A9 <- c(rep(0, 80), rep(1 / 365, 21)),
                              B1 <- c(rep(0, 5), rep(2 / 365, 96)),
                              B2 <- c(rep(0, 18), rep(2 / 365, 83)),
                              B3 <- c(rep(0, 50), rep(2 / 365, 51)),
                              B4 <- c(rep(0, 55), rep(2 / 365, 46)),
                              B5 <- c(rep(0, 60), rep(2 / 365, 41)),
                              B6 <- c(rep(0, 65), rep(2 / 365, 36)),
                              B7 <- c(rep(0, 70), rep(2 / 365, 31)),
                              B8 <- c(rep(0, 75), rep(2 / 365, 26)),
                              B9 <- c(rep(0, 80), rep(2 / 365, 21)),
                              C1 <- c(rep(0, 5), rep(1 / 730, 96)),
                              C2 <- c(rep(0, 18), rep(1 / 730, 83)),
                              C3 <- c(rep(0, 50), rep(1 / 730, 51)),
                              C4 <- c(rep(0, 55), rep(1 / 730, 46)),
                              C5 <- c(rep(0, 60), rep(1 / 730, 41)),
                              C6 <- c(rep(0, 65), rep(1 / 730, 36)),
                              C7 <- c(rep(0, 70), rep(1 / 730, 31)),
                              C8 <- c(rep(0, 75), rep(1 / 730, 26)),
                              C9 <- c(rep(0, 80), rep(1 / 730, 21)),
                              "fitting_period" <- c(vacc_rate_4, vacc_rate_5, vacc_rate_6))
  
  n_vacc_waves_scenarios <- list(Comparator <- rep(0, 101),
                                 A1 <- c(rep(0, 5), rep(1, 96)),
                                 A2 <- c(rep(0, 18), rep(1, 83)),
                                 A3 <- c(rep(0, 50), rep(1, 51)),
                                 A4 <- c(rep(0, 55), rep(1, 46)),
                                 A5 <- c(rep(0, 60), rep(1, 41)),
                                 A6 <- c(rep(0, 65), rep(1, 36)),
                                 A7 <- c(rep(0, 70), rep(1, 31)),
                                 A8 <- c(rep(0, 75), rep(1, 26)),
                                 A9 <- c(rep(0, 80), rep(1, 21)),
                                 B1 <- c(rep(0, 5), rep(2, 96)),
                                 B2 <- c(rep(0, 18), rep(2, 83)),
                                 B3 <- c(rep(0, 50), rep(2, 51)),
                                 B4 <- c(rep(0, 55), rep(2, 46)),
                                 B5 <- c(rep(0, 60), rep(2, 41)),
                                 B6 <- c(rep(0, 65), rep(2, 36)),
                                 B7 <- c(rep(0, 70), rep(2, 31)),
                                 B8 <- c(rep(0, 75), rep(2, 26)),
                                 B9 <- c(rep(0, 80), rep(2, 21)),
                                 C1 <- c(rep(0, 5), rep(0.5, 96)),
                                 C2 <- c(rep(0, 18), rep(0.5, 83)),
                                 C3 <- c(rep(0, 50), rep(0.5, 51)),
                                 C4 <- c(rep(0, 55), rep(0.5, 46)),
                                 C5 <- c(rep(0, 60), rep(0.5, 41)),
                                 C6 <- c(rep(0, 65), rep(0.5, 36)),
                                 C7 <- c(rep(0, 70), rep(0.5, 31)),
                                 C8 <- c(rep(0, 75), rep(0.5, 26)),
                                 C9 <- c(rep(0, 80), rep(0.5, 21)),
                                 "fitting_period" <- n_vacc_waves_fitting)
  
  vacc_rate_shift_scenarios <- list(Comparator <- 2 * pi / 3,        # Vaccination peaks ~ 10th November
                                    A1 <- 2 * pi / 3,
                                    A2 <- 2 * pi / 3,
                                    A3 <- 2 * pi / 3,
                                    A4 <- 2 * pi / 3,
                                    A5 <- 2 * pi / 3,
                                    A6 <- 2 * pi / 3,
                                    A7 <- 2 * pi / 3,
                                    A8 <- 2 * pi / 3,
                                    A9 <- 2 * pi / 3,
                                    B1 <- pi,                        # Vaccination peaks ~ 23 October and 21st April
                                    B2 <- pi,
                                    B3 <- pi,
                                    B4 <- pi,
                                    B5 <- pi,
                                    B6 <- pi,
                                    B7 <- pi,
                                    B8 <- pi,
                                    B9 <- pi,
                                    C1 <- 2 * pi / 3,
                                    C2 <- 2 * pi / 3,
                                    C3 <- 2 * pi / 3,
                                    C4 <- 2 * pi / 3,
                                    C5 <- 2 * pi / 3,
                                    C6 <- 2 * pi / 3,
                                    C7 <- 2 * pi / 3,
                                    C8 <- 2 * pi / 3,
                                    C9 <- 2 * pi / 3,
                                    "fitting_period" <- pi)
  
  names(vacc_rate_scenarios) <- names(n_vacc_waves_scenarios) <- names(vacc_rate_shift_scenarios) <- 
    c("Comparator", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", 
      "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "fitting_period")
  
  vacc_rate_future <- as.numeric(unlist(vacc_rate_scenarios[vaccination_scenario]))
  n_vacc_waves_future <- as.numeric(unlist(n_vacc_waves_scenarios[vaccination_scenario]))
  vacc_rate_shift_future <- as.numeric(unlist(vacc_rate_shift_scenarios[vaccination_scenario]))
  
  # Contact rates based on earlier bootstrapping
  
  contacts <- c(0.85, 0.96, 1.07, 1.18, 1.29, 1.39, 1.48, 1.55, 1.62, 1.66, 1.7, 1.71, 1.72, 1.72, 1.7, 1.68, 1.64, 1.6, 
                1.56, 1.51, 1.46, 1.42, 1.38, 1.35, 1.32, 1.3, 1.28, 1.27, 1.26, 1.25, 1.25, 1.25, 1.25, 1.26, 1.26, 1.27, 
                1.27, 1.28, 1.28, 1.27, 1.27, 1.27, 1.27, 1.26, 1.25, 1.24, 1.22, 1.2, 1.17, 1.15, 1.12, 1.1, 1.08, 1.06, 
                1.04, 1.02, 1.01, 0.99, 0.97, 0.95, 0.93, 0.91, 0.89, 0.87, 0.84, 0.82, 0.79, 0.76, 0.73, 0.7, 0.68, 0.65, 
                0.63, 0.61, 0.6, 0.58, 0.58, 0.57, 0.57, 0.56, 0.56, 0.56, 0.56, 0.56, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 
                0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57)
  
  # Calculating parameter factor changes according to variant scenario
  
  if (variant_scenario == "Base case" | variant_scenario == "Scenario 2" | variant_scenario == "Scenario 3") {
    protection_factor <- 0.92
    beta_factor <- 1
  } else if (variant_scenario == "Scenario 4") {
    n_yearly_variants <- 3
    protection_factor <- 0.95
    beta_factor <- 1
  } else {
    print("UNRECOGNISED SCENARIO")
    protection_factor <- 1
    beta_factor <- 1
  }
  
  protection_factor_scenario_2 <- 0.92
  beta_factor_scenario_2 <- 1.10
  prop_symp_factor_scenario_3 <- 1
  prop_hosp_factor_scenario_3 <- 1.80
  prop_hosp_mortalities_factor_scenario_3 <- 1.80
  
  VE_S_init <- 65
  VE_H_init <- 90
  VE_HM_init <- 90
  VE_T_init <- 68
  wan_VE_H_init <- 60
  wan_VE_HM_init <- 60
  wan_VE_T_init <- 55
  
  # Infection severity rates
  
  inf_symp_rate_init <- 0.803
  inf_hosp_rate_init <- 0.32
  hosp_mort_rate_init <- 0.499
  inf_symp_age_scaling <- c(rep(0.664, 19), rep(0.846, 41), rep(1, 41))
  age_scaling <- c(rep(0.00513, 6), rep(0.00563, 12), rep(0.074, 47), rep(0.488, 20), rep(1, 16))
  
  inf_symp_rate_base_case <- rep(inf_symp_rate_init, max(SIMTIME_future, SIMTIME_0))
  inf_hosp_rate_base_case <- rep(inf_hosp_rate_init, max(SIMTIME_future, SIMTIME_0))
  hosp_mort_rate_base_case <- rep(hosp_mort_rate_init, max(SIMTIME_future, SIMTIME_0))
  
  # Setting up increased severity in variant scenario 3
  
  if (variant_scenario == "Scenario 3") {
    
    fifth_element <- ((seq(1:floor(n_years_future / 5)) * 5) - 1) * 365
    
    # Setting new severity rates for variant scenario 3
    
    inf_symp_rate_scenario_3 <- calculate_severity_parameters_scenario_3(inf_symp_rate_init, prop_symp_factor_scenario_3)
    inf_hosp_rate_scenario_3 <- calculate_severity_parameters_scenario_3(inf_hosp_rate_init, prop_hosp_factor_scenario_3)
    hosp_mort_rate_scenario_3 <- calculate_severity_parameters_scenario_3(hosp_mort_rate_init, prop_hosp_mortalities_factor_scenario_3)
    
    # Setting severity rates over the time horizon to the base case
    
    inf_symp_rate <- inf_symp_rate_base_case
    inf_hosp_rate <- inf_hosp_rate_base_case
    hosp_mort_rate <- hosp_mort_rate_base_case
    
    # Replace the severity rates for every fifth year with the scenario 3 variant severity rates (only one variant in every fifth year)
    
    for (j in 1:length(fifth_element)) {
      inf_symp_rate[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- inf_symp_rate_scenario_3
      inf_hosp_rate[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- inf_hosp_rate_scenario_3
      hosp_mort_rate[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- hosp_mort_rate_scenario_3
    }
    
  } else {
    
    # If variant scenario is not scenario 3
    
    inf_symp_rate <- inf_symp_rate_base_case
    inf_hosp_rate <- inf_hosp_rate_base_case
    hosp_mort_rate <- hosp_mort_rate_base_case
    
  }
  
  # Set fixed parameters
  
  par <- list(ageing <- n_age_groups / (365 * 101),
              lat_period <- 3.1,
              prop_symp <- outer(inf_symp_rate, inf_symp_age_scaling, "*"),
              symp_period <- 9,
              rec_rate <- 1 / symp_period,
              prop_hosp <- outer(inf_hosp_rate, age_scaling, "*"),
              hosp_period <- c(rep(1.8, 5), rep(2.8, 7), rep(5, 6), rep(6.4, 32), rep(9, 15), rep(10.8, 10), rep(11.2, 10), rep(11.4, 16)),
              prop_hosp_mortalities <- outer(hosp_mort_rate, age_scaling, "*"),
              vacc_rate <- c(vacc_rate_1, vacc_rate_2, vacc_rate_3),
              vacc_rate_shift <- pi,
              n_vacc_waves <- 0,         # No vaccination for under 12s, two vaccinations a year for 12-100 year olds
              vacc_uptake_rate <- c(rep(0, 5), rep(6.4, 7), rep(31.7, 4), rep(31.7, 2), rep(31.7, 2), rep(40.3, 5), rep(40.9, 5), 
                                    rep(44.0, 5), rep(49.5, 5), rep(57.0, 5), rep(64.8, 5), rep(64.8, 5), rep(64.8, 5), 
                                    rep(64.8, 5), rep(64.8, 5), rep(70.1, 5), rep(75.5, 5), rep(75.7, 21)) / 100,
              # Using uptake rates for 0-49 year olds from third dose uptake, 50-64 from autumn 2022 uptake and 65+ from autumn 2023 uptake
              wane_rate <- 2 / 365,                      # Immunity wanes after six months
              red_I <- 0,
              red_S <- 0,
              red_H <- 0,
              red_HM <- 0,
              wan_I <- 0,
              wan_S <- 0,
              wan_H <- 0,
              wan_HM <- 0,
              reinf <- 0,
              red_T <- 0,
              wan_T <- 0,
              tau <- 0,
              beta <- 0,
              num_peaks <- 2,
              inf_shift <- pi / 4,
              inf_fluctuation <- 0.125,
              n_contacts <- contacts,
              total_pop <- total_pop,
              n_age_groups <- n_age_groups,
              all_cause_death_rate <- c(4.2, rep(0.2, 4), rep(0.1, 5), rep(0.1, 5), rep(0.3, 5), rep(0.4, 5), rep(0.4, 5), rep(0.7, 5), rep(1.0, 5),
                                        rep(1.5, 5), rep(2.3, 5), rep(3.4, 5), rep(5.0, 5), rep(7.6, 5), rep(12.1, 5), rep(19.3, 5),
                                        rep(31.5, 5), rep(58.3, 5), rep(105.9, 5), rep(226.1, 10), 0) / 365000,          
              # Using yearly mortality rates per 1000 people from 2023 from the ONS, so dividing by 365*1000 (https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsregisteredinenglandandwalesseriesdrreferencetables)
              birth_rate <- 700000)                               # Number of projected yearly births in England proportional to the UK
  
  names(par) <- c("Ageing", "Latent period", "Proportion symptomatic", "Symptomatic period", "Recovery rate (asymptomatic)", "Proportion hospitalised",
                  "Hospital length of stay", "Proportion of hospital mortalities", "Vaccination rate", "Vaccination rate shift",
                  "Number of vaccination waves", "Vaccine uptake rate", "Wane rate", "Reduction in infectiousness", "Reduction in symptomatic infections",
                  "Reduction in hospitalisations", "Reduction in hospital mortalities", "Waned reduction in infectiousness", 
                  "Waned reduction in symptomatic infections", "Waned reduction in hospitalisations", "Waned reduction in hospital mortalities", 
                  "Protection against reinfection", "Reduction in transmission", "Waned reduction in transmission", 
                  "Reduction in asymptomatic infectiousness", "Transmission rate", "Number of annual peaks in transmission",
                  "Infection shift", "Amount of fluctuation in infectiousness", "Number of contacts per age group", "Total population size",
                  "Number of age groups", "All-cause death rate", "Birth rate")
  
  all_cause_death_rate <<- par$all_cause_death_rate
  
  # Initial state values for the ODE system, followed by the time sequence over which the ODEs are solved
  
  # Set initial number infected to the number of people testing positive in England week ending 28 December 2022, distributed
  # according to contact rates (https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/bulletins/coronaviruscovid19infectionsurveypilot/6january2023)
  
  n_cases <- 2463000
  I_init <- round(n_cases * (contacts/sum(contacts)))
  
  # Set initial number vaccinated to the percentage of the population in each age group that had received a third dose on 31 Dec 2022 - the number infected
  
  V_init <- round(c(rep(0, 5), rep(0.002, 7), rep(0.011, 4), rep(0.072, 2), rep(0.239, 2), rep(0.400, 5), rep(0.406, 5), 
                    rep(0.437, 5), rep(0.492, 5), rep(0.568, 5), rep(0.646, 5), rep(0.737, 5), rep(0.794, 5), rep(0.836, 5), 
                    rep(0.874, 5), rep(0.909, 5), rep(0.936, 5), rep(0.941, 21)) * (age_dist - I_init), 2)
  
  # Set initial number waned to the percentage in each age group that had received a second dose on 31 Dec 2022 and had not received a third - the number infected
  
  W_init <- round(age_dist * c(rep(0, 5), rep(0.056, 7), rep(0.300, 4), rep(0.400, 2), rep(0.333, 2), rep(0.244, 5), rep(0.225, 5), 
                               rep(0.210, 5), rep(0.192, 5), rep(0.168, 5), rep(0.144, 5), rep(0.110, 5), rep(0.085, 5), 
                               rep(0.065, 5), rep(0.045, 5), rep(0.029, 5), rep(0.020, 5), rep(0.018, 21)), 2)
  
  n_recovered <- 19037525                    # Approx cumulative number of first episodes in England as of 31 December 2022
  
  E_init <- rep(0, n_age_groups)
  A_init <- 0 * age_dist
  H_init <- 0 * age_dist
  R_init <- round(pmin(n_recovered * (contacts/sum(contacts)), age_dist - I_init - V_init - W_init), 2)      
  # Number of recovered in each age group equal to minimum of estimated number recovered and the remaining population 
  # in that age group that's not infected, vaccinated, or waned
  D_init <- rep(0, n_age_groups)
  Ev_init <- rep(0, n_age_groups)
  Iv_init <- 0 * age_dist
  Av_init <- 0 * age_dist
  Hv_init <- rep(0, n_age_groups)
  Rv_init <- 0 * age_dist
  Ew_init <- rep(0, n_age_groups)
  Iw_init <- 0 * age_dist
  Aw_init <- 0 * age_dist
  Hw_init <- rep(0, n_age_groups)
  Rw_init <- 0 * age_dist
  
  S_init <- round(age_dist - (E_init + I_init + A_init + H_init + R_init + D_init + V_init + Ev_init + Iv_init + Av_init + 
                                Hv_init + Rv_init + W_init + Ew_init + Iw_init + Aw_init + Hw_init + Rw_init), 2)
  
  II_init <- rep(0, n_age_groups)
  HH_init <- rep(0, n_age_groups)
  VV_init <- rep(0, n_age_groups)
  VP_init <- rep(0, n_age_groups)
  VB_init <- rep(0, n_age_groups)
  IIS_init <- rep(0, n_age_groups)
  RR_init <- rep(0, n_age_groups)
  RRV_init <- rep(0, n_age_groups)
  RRW_init <- rep(0, n_age_groups)
  
  # Model burn-in input
  
  input_init <- c(S_init, E_init, I_init, A_init, H_init, R_init, D_init, V_init, Ev_init, Iv_init, Av_init, Hv_init, 
                  Rv_init, W_init, Ew_init, Iw_init, Aw_init, Hw_init, Rw_init, II_init, HH_init, VV_init, VP_init, 
                  VB_init, IIS_init, RR_init, RRV_init, RRW_init)
  
  # Check for negative values in model initial input and checking population is constant
  
  combined_pop_compartments_init <- S_init + E_init + I_init + A_init + H_init + R_init + D_init + V_init + 
    Ev_init + Iv_init + Av_init + Hv_init + Rv_init + W_init + Ew_init + 
    Iw_init + Aw_init + Hw_init + Rw_init
  
  if(is.null(fixr::check_for_negative_values(input_init)) == FALSE) {
    stop("Negative values in input_init")
  } else if (round(sum(combined_pop_compartments_init)) != round(total_pop)) {
    stop("Total population in input_init not equal to initial population")
  }
  
  # Initialising hospitalisation counting, dates, parameter sets and hospital data set vectors
  
  inc_hosp <- inc_deaths <- inc_hosp_0_5 <- inc_deaths_0_5 <- inc_hosp_6_17 <- inc_deaths_6_17 <- inc_hosp_18_64 <- 
    inc_deaths_18_64 <- inc_hosp_65_84 <- inc_deaths_65_84 <- inc_hosp_85_ <- inc_deaths_85_ <- inc_vacc <- matrix(nrow = length(times) - 1, ncol = ncol(par_fit))
  cum_hosp <- cum_deaths <- matrix(nrow = length(times), ncol = ncol(par_fit))
  total_inf_future <- total_vacc_future <- total_hosp_future <- total_deaths_future <- total_inf_future_0_5 <- total_inf_future_6_17 <-
    total_inf_future_18_64 <- total_inf_future_65_84 <- total_inf_future_85_ <- total_hosp_future_0_5 <- total_hosp_future_6_17 <-
    total_hosp_future_18_64 <- total_hosp_future_65_84 <- total_hosp_future_85_<- total_deaths_future_0_5 <- total_deaths_future_6_17 <-
    total_deaths_future_18_64 <- total_deaths_future_65_84 <- total_deaths_future_85_ <- total_icu_future_0_5 <-
    total_icu_future_6_17 <- total_icu_future_18_64 <- total_icu_future_65_84 <- total_icu_future_85_ <-life_years_lost_future <-
    rep(NA, ncol(par_fit))
  
  # Initialising counting for number in each age group in each disease state every year in future horizon for counting costs/QALYs
  
  n_healthy_yearly <- n_symptomatic_yearly <- n_hospitalised_yearly <- n_dead_yearly <-  n_recovered_past_year_yearly <- 
    n_recovered_vaccinated_past_year_yearly <- n_recovered_waned_past_year_yearly <- inc_inf_yearly <- inc_hosp_yearly <- 
    inc_deaths_yearly <- inc_vacc_yearly <- inc_symp_inf_yearly <- inc_rec_yearly <- vector(mode='list', length = ncol(par_fit))
  
  # Setting dates vector
  
  dates <- seq(as.Date("2022-12-31"), by = "weeks", length = floor(SIMTIME / 7))
  breaks.vec <- seq(min(dates), max(dates), by = "2 months")
  
  parameter_set <- rep(paste("Parameter set", "1"), floor(SIMTIME / 7))
  if (ncol(par_fit) > 1) {
    for (i in 2:ncol(par_fit)) {
      parameter_set <- append(parameter_set, rep(paste("Parameter set", i), floor(SIMTIME / 7)))
    }
  }
  
  # Extracting hospital admissions data from excel file
  
  hosp_data <- read_excel("Data/Weekly hospital admissions 2023 (England).xlsx", sheet = "2023")
  deaths_data <- read_excel("Data/Weekly deaths 2023 (England).xlsx", sheet = "2023")
  
  # Daily hospital admissions for all age groups
  
  hosp_data_daily <- as.data.frame(hosp_data)
  hosp_data_daily <- hosp_data_daily[-c(1), -c(8,9)]
  hosp_data_daily[, 1] <- excel_numeric_to_date(as.numeric(as.character(hosp_data_daily[, 1])), date_system = "modern")
  colnames(hosp_data_daily) <- c("Date", "Total", "0-5", "6-17", "18-64", "65-84", "85+")
  
  hosp_data_total <<- colSums(hosp_data_daily[, 2:7])
  
  # Daily deaths for all age groups
  
  deaths_data_weekly <- as.data.frame(deaths_data[-53,])
  colnames(deaths_data_weekly) <- c("Date", "Total", "0-5", "6-17", "18-64", "65-84", "85+")
  deaths_data_weekly <<- deaths_data_weekly
  
  deaths_data_total <<- colSums(deaths_data_weekly[, 2:7])
  
  # Summing daily admissions to give weekly admissions for all groups
  
  hosp_data_weekly <<- rowsum(hosp_data_daily[1:364,], rep(1:52, each = 7))
  
  # Calculating number of deaths during the burn-in
  
  burn_in_deaths <- rep(0, ncol(par_fit))
  
  # Run model for each parameter set
  
  for (i in 1:ncol(par_fit)) {
    
    print(paste("Parameter set", i))
    
    # Setting parameters
    
    par[["Reduction in asymptomatic infectiousness"]] <- tau <- par_fit["Relative transmission from asymptomatic infections", i]
    
    age_scaling <- c(rep(par_fit["Infection severity age scaling 0-5", i], 6), rep(par_fit["Infection severity age scaling 6-17", i], 12),
                     rep(par_fit["Infection severity age scaling 18-64", i], 47), rep(par_fit["Infection severity age scaling 65-84", i], 20), 
                     rep(par_fit["Infection severity age scaling 85+", i], 16))
    
    par[["Proportion symptomatic"]] <- prop_symp <- outer(inf_symp_rate_base_case[1:SIMTIME_0], inf_symp_age_scaling, "*")
    
    par[["Proportion hospitalised"]] <- prop_hosp <- outer(inf_hosp_rate_base_case[1:SIMTIME_0], age_scaling, "*")
    
    par[["Proportion of hospital mortalities"]] <- prop_hosp_mortalities <- 
      outer(hosp_mort_rate_base_case[1:SIMTIME_0], age_scaling * c(rep(par_fit["Hospital mortality age scaling 0-5", i], 6), 
                                                                   rep(par_fit["Hospital mortality age scaling 6-17", i], 12),
                                                                   rep(par_fit["Hospital mortality age scaling 18-64", i], 47), 
                                                                   rep(par_fit["Hospital mortality age scaling 65-84", i], 20), 
                                                                   rep(1, 16)), "*")
    
    n_vacc_waves <- par[["Number of vaccination waves"]] <- n_vacc_waves_burn_in         # For burn-in period: no vaccination for under 12s, two vaccinations a year for 12-100 year olds
    
    vacc_rate <- par[["Vaccination rate"]] <- c(vacc_rate_1, vacc_rate_2, vacc_rate_3)       # Set original vaccination rates for burn-in period
    
    vacc_rate_shift <- par[["Vaccination rate shift"]] <- pi                             # Shift for 2 vacc peaks a year so vacc peaks in Apr-Jun and Sep-Dec
    
    # Setting fitted parameters according to the variant scenario
    
    par[["Transmission rate"]] <- beta <- rep(par_fit["Beta", i], max(SIMTIME_0, SIMTIME_fit))
    par[["Protection against reinfection"]] <- reinf <- rep(par_fit["Protection against reinfection", i], max(SIMTIME_0, SIMTIME_fit))
    
    # Setting initial values for fitted parameters
    
    reinf_init <- par_fit["Protection against reinfection", i]
    beta_init <- par_fit["Beta", i]
    
    # Setting VE against infection to equal natural immunity against infection
    
    red_I_init <- reinf_init
    VE_I_init <- 100 * (1 - red_I_init)
    
    # Setting waned VE against infection to fitted parameter * VE against infection
    
    wan_VE_I_init <- VE_I_init * par_fit["Waned VE against infection scaling", i]
    wan_I_init <- 1 - (wan_VE_I_init / 100)
    red_I <- par[["Reduction in infectiousness"]] <- rep(red_I_init, max(SIMTIME_0, SIMTIME_fit))
    wan_I <- par[["Waned reduction in infectiousness"]] <- rep(wan_I_init, max(SIMTIME_0, SIMTIME_fit))
    
    # Setting waned VE against symptoms to equal 10% or at least the waned VE against infection
    
    wan_VE_S_init <- max(10, wan_VE_I_init)
    
    # Calculating vaccine effectiveness initial values from fitted parameters
    
    red_S_init <- calculate_reduction(current_reduction = red_I_init, desired_VE = VE_S_init)
    red_H_init <- calculate_reduction(current_reduction = 1 - (VE_S_init / 100), desired_VE = VE_H_init)
    red_HM_init <- calculate_reduction(current_reduction = 1 - (VE_H_init / 100), desired_VE = VE_HM_init)
    red_T_init <- calculate_reduction(current_reduction = 1, desired_VE = VE_T_init)
    wan_S_init <- calculate_reduction(current_reduction = wan_I_init, desired_VE = wan_VE_S_init)
    wan_H_init <- calculate_reduction(current_reduction = 1 - (wan_VE_S_init / 100), desired_VE = wan_VE_H_init)
    wan_HM_init <- calculate_reduction(current_reduction = 1 - (wan_VE_H_init / 100), desired_VE = wan_VE_HM_init)
    wan_T_init <- calculate_reduction(current_reduction = 1, desired_VE = wan_VE_T_init)
    
    # Setting vaccine effectiveness to the initial value for the burn-in period
    
    red_S <- par[["Reduction in symptomatic infections"]] <- rep(red_S_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period
    red_H <- par[["Reduction in hospitalisations"]] <- rep(red_H_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period
    red_HM <- par[["Reduction in hospital mortalities"]] <- rep(red_HM_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period
    red_T <- par[["Reduction in transmission"]] <- rep(red_T_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period
    wan_S <- par[["Waned reduction in symptomatic infections"]] <- rep(wan_S_init, max(SIMTIME_0, SIMTIME_fit))           # Set vaccine effectiveness to the initial value for the burn-in period/fitting period, whichever is longer
    wan_H <- par[["Waned reduction in hospitalisations"]] <- rep(wan_H_init, max(SIMTIME_0, SIMTIME_fit))           # Set waned vaccine effectiveness to the initial value for the burn-in period
    wan_HM <- par[["Waned reduction in hospital mortalities"]] <- rep(wan_HM_init, max(SIMTIME_0, SIMTIME_fit))           # Set waned vaccine effectiveness to the initial value for the burn-in period
    wan_T <- par[["Waned reduction in transmission"]] <- rep(wan_T_init, max(SIMTIME_0, SIMTIME_fit))           # Set waned vaccine effectiveness to the initial value for the burn-in period
    
    # Data frame containing solutions to the ODEs over the specified time sequence (burn-in period)
    
    out_0 <- as.data.frame(deSolve::ode(y = input_init, times = times_0, func = covid_model_odes, parms = par, method = "euler"))
    
    # Label the output columns for each compartment
    
    S_0 <- out_0[, 2:(n_age_groups + 1)]
    E_0 <- out_0[, (n_age_groups + 2):(2 * n_age_groups + 1)]
    I_0 <- out_0[, (2 * n_age_groups + 2):(3 * n_age_groups + 1)]
    A_0 <- out_0[, (3 * n_age_groups + 2):(4 * n_age_groups + 1)]
    H_0 <- out_0[, (4 * n_age_groups + 2):(5 * n_age_groups + 1)]
    R_0 <- out_0[, (5 * n_age_groups + 2):(6 * n_age_groups + 1)]
    D_0 <- out_0[, (6 * n_age_groups + 2):(7 * n_age_groups + 1)]
    V_0 <- out_0[, (7 * n_age_groups + 2):(8 * n_age_groups + 1)]
    Ev_0 <- out_0[, (8 * n_age_groups + 2):(9 * n_age_groups + 1)]
    Iv_0 <- out_0[, (9 * n_age_groups + 2):(10 * n_age_groups + 1)]
    Av_0 <- out_0[, (10 * n_age_groups + 2):(11 * n_age_groups + 1)]
    Hv_0 <- out_0[, (11 * n_age_groups + 2):(12 * n_age_groups + 1)]
    Rv_0 <- out_0[, (12 * n_age_groups + 2):(13 * n_age_groups + 1)]
    W_0 <- out_0[, (13 * n_age_groups + 2):(14 * n_age_groups + 1)]
    Ew_0 <- out_0[, (14 * n_age_groups + 2):(15 * n_age_groups + 1)]
    Iw_0 <- out_0[, (15 * n_age_groups + 2):(16 * n_age_groups + 1)]
    Aw_0 <- out_0[, (16 * n_age_groups + 2):(17 * n_age_groups + 1)]
    Hw_0 <- out_0[, (17 * n_age_groups + 2):(18 * n_age_groups + 1)]
    Rw_0 <- out_0[, (18 * n_age_groups + 2):(19 * n_age_groups + 1)]
    II_0 <- out_0[, (19 * n_age_groups + 2):(20 * n_age_groups + 1)]
    HH_0 <- out_0[, (20 * n_age_groups + 2):(21 * n_age_groups + 1)]
    VV_0 <- out_0[, (21 * n_age_groups + 2):(22 * n_age_groups + 1)]
    VP_0 <- out_0[, (22 * n_age_groups + 2):(23 * n_age_groups + 1)]
    VB_0 <- out_0[, (23 * n_age_groups + 2):(24 * n_age_groups + 1)]
    IIS_0 <- out_0[, (24 * n_age_groups + 2):(25 * n_age_groups + 1)]
    RR_0 <- out_0[, (25 * n_age_groups + 2):(26 * n_age_groups + 1)]
    RRV_0 <- out_0[, (26 * n_age_groups + 2):(27 * n_age_groups + 1)]
    RRW_0 <- out_0[, (27 * n_age_groups + 2):(28 * n_age_groups + 1)]
    
    # Check for negative values in model burn-in output
    
    if(is.null(fixr::check_for_negative_values(out_0)) == FALSE) {
      stop("Negative values in out_0")
    }
    
    # Run the model again with the previous output as the input
    
    input_fit <- as.numeric(cbind(S_0, E_0, I_0, A_0, H_0, R_0, D_0, V_0, Ev_0, Iv_0, Av_0, Hv_0, Rv_0, W_0, Ew_0, Iw_0,
                                  Aw_0, Hw_0, Rw_0, II_0, HH_0, VV_0, VP_0, VB_0, IIS_0, RR_0, RRV_0, RRW_0)[SIMTIME_0, ])
    
    # Reset the Covid deaths compartments in the input to 0
    
    input_fit[607:707] <- rep(0, 101)
    
    burn_in_deaths[i] <- sum(D_0[nrow(D_0), ])
    
    # Second model run parameters (fitting period)
    
    n_vacc_waves <- par[["Number of vaccination waves"]] <- n_vacc_waves_fitting      # Two boosters for 75+, one booster for 50-74, no vaccination for under 50
    vacc_rate <- par[["Vaccination rate"]] <- c(vacc_rate_4, vacc_rate_5, vacc_rate_6)       # Set new vaccination rates
    
    # Reset severity parameters so matrix has same number of rows as SIMTIME_fit
    
    par[["Proportion symptomatic"]] <- prop_symp <- outer(inf_symp_rate_base_case[1:SIMTIME_fit], inf_symp_age_scaling, "*")
    
    par[["Proportion hospitalised"]] <- prop_hosp <- outer(inf_hosp_rate_base_case[1:SIMTIME_fit], age_scaling, "*")
    
    par[["Proportion of hospital mortalities"]] <- prop_hosp_mortalities <- 
      outer(hosp_mort_rate_base_case[1:SIMTIME_fit], age_scaling * c(rep(par_fit["Hospital mortality age scaling 0-5", i], 6), 
                                                                     rep(par_fit["Hospital mortality age scaling 6-17", i], 12),
                                                                     rep(par_fit["Hospital mortality age scaling 18-64", i], 47), 
                                                                     rep(par_fit["Hospital mortality age scaling 65-84", i], 20), 
                                                                     rep(1, 16)), "*")
    
    # Run the model again,  beginning from 31/12/2022, for the year 2023
    
    out_fit <- as.data.frame(deSolve::ode(y = input_fit, times = times_fit, func = covid_model_odes, parms = par, method = "euler"))
    
    # Label the output columns for each compartment
    
    S_fit <- out_fit[, 2:(n_age_groups + 1)]
    E_fit <- out_fit[, (n_age_groups + 2):(2 * n_age_groups + 1)]
    I_fit <- out_fit[, (2 * n_age_groups + 2):(3 * n_age_groups + 1)]
    A_fit <- out_fit[, (3 * n_age_groups + 2):(4 * n_age_groups + 1)]
    H_fit <- out_fit[, (4 * n_age_groups + 2):(5 * n_age_groups + 1)]
    R_fit <- out_fit[, (5 * n_age_groups + 2):(6 * n_age_groups + 1)]
    D_fit <- out_fit[, (6 * n_age_groups + 2):(7 * n_age_groups + 1)]
    V_fit <- out_fit[, (7 * n_age_groups + 2):(8 * n_age_groups + 1)]
    Ev_fit <- out_fit[, (8 * n_age_groups + 2):(9 * n_age_groups + 1)]
    Iv_fit <- out_fit[, (9 * n_age_groups + 2):(10 * n_age_groups + 1)]
    Av_fit <- out_fit[, (10 * n_age_groups + 2):(11 * n_age_groups + 1)]
    Hv_fit <- out_fit[, (11 * n_age_groups + 2):(12 * n_age_groups + 1)]
    Rv_fit <- out_fit[, (12 * n_age_groups + 2):(13 * n_age_groups + 1)]
    W_fit <- out_fit[, (13 * n_age_groups + 2):(14 * n_age_groups + 1)]
    Ew_fit <- out_fit[, (14 * n_age_groups + 2):(15 * n_age_groups + 1)]
    Iw_fit <- out_fit[, (15 * n_age_groups + 2):(16 * n_age_groups + 1)]
    Aw_fit <- out_fit[, (16 * n_age_groups + 2):(17 * n_age_groups + 1)]
    Hw_fit <- out_fit[, (17 * n_age_groups + 2):(18 * n_age_groups + 1)]
    Rw_fit <- out_fit[, (18 * n_age_groups + 2):(19 * n_age_groups + 1)]
    II_fit <- out_fit[, (19 * n_age_groups + 2):(20 * n_age_groups + 1)]
    HH_fit <- out_fit[, (20 * n_age_groups + 2):(21 * n_age_groups + 1)]
    VV_fit <- out_fit[, (21 * n_age_groups + 2):(22 * n_age_groups + 1)]
    VP_fit <- out_fit[, (22 * n_age_groups + 2):(23 * n_age_groups + 1)]
    VB_fit <- out_fit[, (23 * n_age_groups + 2):(24 * n_age_groups + 1)]
    IIS_fit <- out_fit[, (24 * n_age_groups + 2):(25 * n_age_groups + 1)]
    RR_fit <- out_fit[, (25 * n_age_groups + 2):(26 * n_age_groups + 1)]
    RRV_fit <- out_fit[, (26 * n_age_groups + 2):(27 * n_age_groups + 1)]
    RRW_fit <- out_fit[, (27 * n_age_groups + 2):(28 * n_age_groups + 1)]
    
    # Check for negative values in model fitting period output
    
    if(is.null(fixr::check_for_negative_values(out_fit)) == FALSE) {
      stop("Negative values in out_fit")
    }
    
    # Run the model again with the previous output as the input, only if SIMTIME_future > 0
    
    if (SIMTIME_future > 0) {
      
      input <- as.numeric((cbind(S_fit, E_fit, I_fit, A_fit, H_fit, R_fit, D_fit, V_fit, Ev_fit, Iv_fit, Av_fit, Hv_fit, Rv_fit, W_fit, Ew_fit, Iw_fit,
                                 Aw_fit, Hw_fit, Rw_fit, II_fit, HH_fit, VV_fit, VP_fit, VB_fit, IIS_fit, RR_fit, RRV_fit, RRW_fit))[SIMTIME_fit, ])
      
      # Calculating parameters for base case/no variant scenario
      
      variant_parameters <- calculate_variant_parameters(reinf_init, red_I_init, VE_S_init, VE_H_init, VE_HM_init, 
                                                         VE_T_init, wan_I_init, wan_VE_S_init, wan_VE_H_init, 
                                                         wan_VE_HM_init, wan_VE_T_init, beta_init, 
                                                         protection_factor, beta_factor, n_yearly_variants)
      
      names(variant_parameters) <- c("reinf_variants", "red_I_variants", "red_S_variants", "red_H_variants", 
                                     "red_HM_variants", "red_T_variants", "wan_I_variants", "wan_S_variants", 
                                     "wan_H_variants", "wan_HM_variants", "wan_T_variants", "beta_variants")
      
      # Setting base case variant parameters with new variants once every "variant_every_n_years"
      
      reinf_variants <- rep(reinf_init, SIMTIME_future)
      red_I_variants <- rep(red_I_init, SIMTIME_future)
      red_S_variants <- rep(red_S_init, SIMTIME_future)
      red_H_variants <- rep(red_H_init, SIMTIME_future)
      red_HM_variants <- rep(red_HM_init, SIMTIME_future)
      red_T_variants <- rep(red_T_init, SIMTIME_future)
      wan_I_variants <- rep(wan_I_init, SIMTIME_future)
      wan_S_variants <- rep(wan_S_init, SIMTIME_future)
      wan_H_variants <- rep(wan_H_init, SIMTIME_future)
      wan_HM_variants <- rep(wan_HM_init, SIMTIME_future)
      wan_T_variants <- rep(wan_T_init, SIMTIME_future)
      beta_variants <- rep(beta_init, SIMTIME_future)
      
      nth_element <- ((seq(1:floor(n_years_future / variant_every_n_years)) * variant_every_n_years) - 1) * 365
      
      for (j in 1:length(nth_element)) {
        reinf_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$reinf_variants
        red_I_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_I_variants
        red_S_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_S_variants
        red_H_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_H_variants
        red_HM_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_HM_variants
        red_T_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$red_T_variants
        wan_I_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_I_variants
        wan_S_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_S_variants
        wan_H_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_H_variants
        wan_HM_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_HM_variants
        wan_T_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$wan_T_variants
        beta_variants[(nth_element[j] + 1):(nth_element[j] + 365)] <- variant_parameters$beta_variants
      }
      
      # Setting variant parameters for scenario 2
      
      if (variant_scenario == "Scenario 2") {
        
        fifth_element <- ((seq(1:floor(n_years_future / 5)) * 5) - 1) * 365
        
        # Setting parameters in case of variant scenario 2 every five years
        
        variant_parameters_scenario_2 <- calculate_variant_parameters(reinf_init, red_I_init, VE_S_init, VE_H_init, 
                                                                      VE_HM_init, VE_T_init, wan_I_init, 
                                                                      wan_VE_S_init, wan_VE_H_init, wan_VE_HM_init, 
                                                                      wan_VE_T_init, beta_init, 
                                                                      protection_factor_scenario_2, beta_factor_scenario_2, 
                                                                      n_yearly_variants = 1)
        
        names(variant_parameters_scenario_2) <- c("reinf_variants_scenario_2", "red_I_variants_scenario_2", 
                                                  "red_S_variants_scenario_2", "red_H_variants_scenario_2", 
                                                  "red_HM_variants_scenario_2", "red_T_variants_scenario_2", 
                                                  "wan_I_variants_scenario_2", "wan_S_variants_scenario_2", 
                                                  "wan_H_variants_scenario_2", "wan_HM_variants_scenario_2", 
                                                  "wan_T_variants_scenario_2", "beta_variants_scenario_2")
        
        # Replace the new variant of every fifth year with the scenario 2 variant (only one variant in every fifth year)
        
        for (j in 1:length(fifth_element)) {
          reinf_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$reinf_variants_scenario_2
          red_I_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_I_variants_scenario_2
          red_S_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_S_variants_scenario_2
          red_H_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_H_variants_scenario_2
          red_HM_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_HM_variants_scenario_2
          red_T_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$red_T_variants_scenario_2
          wan_I_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_I_variants_scenario_2
          wan_S_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_S_variants_scenario_2
          wan_H_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_H_variants_scenario_2
          wan_HM_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_HM_variants_scenario_2
          wan_T_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$wan_T_variants_scenario_2
          beta_variants[(fifth_element[j] + 1):(fifth_element[j] + 365)] <- variant_parameters_scenario_2$beta_variants_scenario_2
        }
      }
      
      # Setting variant parameters
      
      reinf <- par[["Protection against reinfection"]] <- reinf_variants
      red_I <- par[["Reduction in infectiousness"]] <- red_I_variants
      red_S <- par[["Reduction in symptomatic infections"]] <- red_S_variants                  # Set vaccine effectiveness to change according to new variants
      red_H <- par[["Reduction in hospitalisations"]] <- red_H_variants                  # Set vaccine effectiveness to change according to new variants
      red_HM <- par[["Reduction in hospital mortalities"]] <- red_HM_variants                  # Set vaccine effectiveness to change according to new variants
      red_T <- par[["Reduction in transmission"]] <- red_T_variants                  # Set vaccine effectiveness to change according to new variants
      wan_I <- par[["Waned reduction in infectiousness"]] <- wan_I_variants
      wan_S <- par[["Waned reduction in symptomatic infections"]] <- wan_S_variants                  # Set vaccine effectiveness to change according to new variant
      wan_H <- par[["Waned reduction in hospitalisations"]] <- wan_H_variants                  # Set vaccine effectiveness to change according to new variant
      wan_HM <- par[["Waned reduction in hospital mortalities"]] <- wan_HM_variants                  # Set vaccine effectiveness to change according to new variants
      wan_T <- par[["Waned reduction in transmission"]] <- wan_T_variants                  # Set vaccine effectiveness to change according to new variants
      beta <- par[["Transmission rate"]] <- beta_variants
      
      # Setting severity parameters for variant scenario
      
      prop_symp <- outer(inf_symp_rate, inf_symp_age_scaling, "*")
      
      prop_hosp <- outer(inf_hosp_rate, age_scaling, "*")
      
      prop_hosp_mortalities <- outer(hosp_mort_rate, age_scaling * c(rep(par_fit["Hospital mortality age scaling 0-5", i], 6), 
                                                                     rep(par_fit["Hospital mortality age scaling 6-17", i], 12),
                                                                     rep(par_fit["Hospital mortality age scaling 18-64", i], 47), 
                                                                     rep(par_fit["Hospital mortality age scaling 65-84", i], 20), 
                                                                     rep(1, 16)), "*")  
      
      # If any proportions are greater than 1, set to 1
      
      par[["Proportion symptomatic"]] <- prop_symp <- ifelse(prop_symp > 1, 1, prop_symp)
      
      par[["Proportion hospitalised"]] <- prop_hosp <- ifelse(prop_hosp > 1, 1, prop_hosp)
      
      par[["Proportion of hospital mortalities"]] <- prop_hosp_mortalities <-
        ifelse(prop_hosp_mortalities > 1, 1 , prop_hosp_mortalities)
      
      # Set vaccination scenario
      
      n_vacc_waves <- par[["Number of vaccination waves"]] <- n_vacc_waves_future
      
      vacc_rate <- par[["Vaccination rate"]] <- vacc_rate_future
      
      vacc_rate_shift <- par[["Vaccination rate shift"]] <- vacc_rate_shift_future                 # Vacc rate shift depending on vaccination scenario
      
      # Run the model again,  beginning from the year 2024
      
      out_future <- as.data.frame(deSolve::ode(y = input, times = times_future, func = covid_model_odes, parms = par, method = "euler"))
      
      # Merge 2023 output with the final model run output
      
      out <- rbind(out_fit, out_future)
      
      # Label all output columns
      
      S_future <- out_future[, 2:(n_age_groups + 1)]
      E_future <- out_future[, (n_age_groups + 2):(2 * n_age_groups + 1)]
      I_future <- out_future[, (2 * n_age_groups + 2):(3 * n_age_groups + 1)]
      A_future <- out_future[, (3 * n_age_groups + 2):(4 * n_age_groups + 1)]
      H_future <- out_future[, (4 * n_age_groups + 2):(5 * n_age_groups + 1)]
      R_future <- out_future[, (5 * n_age_groups + 2):(6 * n_age_groups + 1)]
      D_future <- out_future[, (6 * n_age_groups + 2):(7 * n_age_groups + 1)]
      V_future <- out_future[, (7 * n_age_groups + 2):(8 * n_age_groups + 1)]
      Ev_future <- out_future[, (8 * n_age_groups + 2):(9 * n_age_groups + 1)]
      Iv_future <- out_future[, (9 * n_age_groups + 2):(10 * n_age_groups + 1)]
      Av_future <- out_future[, (10 * n_age_groups + 2):(11 * n_age_groups + 1)]
      Hv_future <- out_future[, (11 * n_age_groups + 2):(12 * n_age_groups + 1)]
      Rv_future <- out_future[, (12 * n_age_groups + 2):(13 * n_age_groups + 1)]
      W_future <- out_future[, (13 * n_age_groups + 2):(14 * n_age_groups + 1)]
      Ew_future <- out_future[, (14 * n_age_groups + 2):(15 * n_age_groups + 1)]
      Iw_future <- out_future[, (15 * n_age_groups + 2):(16 * n_age_groups + 1)]
      Aw_future <- out_future[, (16 * n_age_groups + 2):(17 * n_age_groups + 1)]
      Hw_future <- out_future[, (17 * n_age_groups + 2):(18 * n_age_groups + 1)]
      Rw_future <- out_future[, (18 * n_age_groups + 2):(19 * n_age_groups + 1)]
      II_future <- out_future[, (19 * n_age_groups + 2):(20 * n_age_groups + 1)]
      HH_future <- out_future[, (20 * n_age_groups + 2):(21 * n_age_groups + 1)]
      VV_future <- out_future[, (21 * n_age_groups + 2):(22 * n_age_groups + 1)]
      VP_future <- out_future[, (22 * n_age_groups + 2):(23 * n_age_groups + 1)]
      VB_future <- out_future[, (23 * n_age_groups + 2):(24 * n_age_groups + 1)]
      IIS_future <- out_future[, (24 * n_age_groups + 2):(25 * n_age_groups + 1)]
      RR_future <- out_future[, (25 * n_age_groups + 2):(26 * n_age_groups + 1)]
      RRV_future <- out_future[, (26 * n_age_groups + 2):(27 * n_age_groups + 1)]
      RRW_future <- out_future[, (27 * n_age_groups + 2):(28 * n_age_groups + 1)]
      
      # Check for negative values in model future output
      
      if(is.null(fixr::check_for_negative_values(out_future)) == FALSE) {
        stop("Negative values in out_future")
      }
      
      # Calculating number of symptomatic recoveries over the past year for each time point in SIMTIME_future from the future and model fit runs
      
      RR_fit_future <- rbind(RR_fit, RR_future)
      RRV_fit_future <- rbind(RRV_fit, RRV_future)
      RRW_fit_future <- rbind(RRW_fit, RRW_future)
      
      inc_rec_unvacc_fit_future <- apply(RR_fit_future, 2, diff)           # Daily number of unvaccinated recoveries over the fitting and future period
      inc_rec_vacc_fit_future <- apply(RRV_fit_future, 2, diff)     # Daily number of vaccinated recoveries over the fitting and future period
      inc_rec_wan_fit_future <- apply(RRW_fit_future, 2, diff)     # Daily number of waned recoveries over the fitting and future period
      
      cumulative_n_recovered_unvacc_past_year <- apply(inc_rec_unvacc_fit_future, 2, cumsum)
      cumulative_n_recovered_unvacc_past_year_extended <- rbind(rep(0, n_age_groups), cumulative_n_recovered_unvacc_past_year)
      n_recovered_unvacc_past_year <- cumulative_n_recovered_unvacc_past_year[365:nrow(inc_rec_unvacc_fit_future), ] - 
        cumulative_n_recovered_unvacc_past_year_extended[1:(nrow(inc_rec_unvacc_fit_future) - 364), ]
      
      cumulative_n_recovered_vacc_past_year <- apply(inc_rec_vacc_fit_future, 2, cumsum)
      cumulative_n_recovered_vacc_past_year_extended <- rbind(rep(0, n_age_groups), cumulative_n_recovered_vacc_past_year)
      n_recovered_vacc_past_year <- cumulative_n_recovered_vacc_past_year[365:nrow(inc_rec_vacc_fit_future), ] - 
        cumulative_n_recovered_vacc_past_year_extended[1:(nrow(inc_rec_vacc_fit_future) - 364), ]
      
      cumulative_n_recovered_wan_past_year <- apply(inc_rec_wan_fit_future, 2, cumsum)
      cumulative_n_recovered_wan_past_year_extended <- rbind(rep(0, n_age_groups), cumulative_n_recovered_wan_past_year)
      n_recovered_wan_past_year <- cumulative_n_recovered_wan_past_year[365:nrow(inc_rec_wan_fit_future), ] - 
        cumulative_n_recovered_wan_past_year_extended[1:(nrow(inc_rec_wan_fit_future) - 364), ]
      
      n_recovered_past_year_yearly[[i]] <- rowsum(n_recovered_unvacc_past_year, rep(1:n_years_future, each = 365))
      n_recovered_vaccinated_past_year_yearly[[i]] <- rowsum(n_recovered_vacc_past_year, rep(1:n_years_future, each = 365))
      n_recovered_waned_past_year_yearly[[i]] <- rowsum(n_recovered_wan_past_year, rep(1:n_years_future, each = 365))
      
      # Calculating number in each compartment in each age group every year over the time horizon for calculating QALYs later
      
      n_healthy_yearly[[i]] <- rowsum((S_future + E_future + A_future + V_future + Ev_future + Av_future + 
                                         W_future + Ew_future + Aw_future + R_future +
                                         Rv_future + Rw_future), rep(1:n_years_future, each = 365))
      n_symptomatic_yearly[[i]] <- rowsum((I_future + Iv_future + Iw_future), rep(1:n_years_future, each = 365))
      n_hospitalised_yearly[[i]] <- rowsum((H_future + Hv_future + Hw_future), rep(1:n_years_future, each = 365))
      n_dead_yearly[[i]] <- rowsum((as.matrix(D_future) - matrix(as.numeric(D_fit[nrow(D_fit), ]), nrow = nrow(D_future), 
                                                                 ncol = ncol(D_fit), byrow = TRUE)), 
                                   rep(1:n_years_future, each = 365))
      
      # Calculating incremental infections, hospitalisations, deaths vaccinations and symptomatic infections for the future model run
      
      inc_inf_future <- apply(II_future, 2, diff)
      inc_hosp_future <- apply(HH_future, 2, diff)
      inc_deaths_future <- apply(D_future, 2, diff)
      inc_vacc_future <- apply(VV_future, 2, diff)
      inc_symp_inf_future <- apply(IIS_future, 2, diff)
      inc_rec_future <- apply(RR_future + RRV_future + RRW_future, 2, diff)
      
      # Calculating incremental infections, hospitalisations, deaths, vaccinations and symptomatic infections each year in each age group for calculating costs later
      
      inc_inf_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_inf_future), rep(1:n_years_future, each = 365))
      inc_hosp_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_hosp_future), rep(1:n_years_future, each = 365))
      inc_deaths_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_deaths_future), rep(1:n_years_future, each = 365))
      inc_vacc_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_vacc_future), rep(1:n_years_future, each = 365))
      inc_symp_inf_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_symp_inf_future), rep(1:n_years_future, each = 365))
      inc_rec_yearly[[i]] <- rowsum(rbind(rep(0, 101), inc_rec_future), rep(1:n_years_future, each = 365))
      
      # Calculating total cases, hospital admissions, deaths and vaccinations given over the future simulation
      
      total_inf_future[i] <- sum(II_future[SIMTIME_future, ] - II_fit[SIMTIME_fit, ])
      total_vacc_future[i] <- sum(VV_future[SIMTIME_future, ] - VV_fit[SIMTIME_fit, ])
      total_hosp_future[i] <- sum(HH_future[SIMTIME_future, ] - HH_fit[SIMTIME_fit, ])
      total_deaths_future[i] <- sum(D_future[SIMTIME_future, ] - D_fit[SIMTIME_fit, ])
      
      # Calculating total cases in each age group
      
      total_inf_future_0_5[i] <- sum(II_future[SIMTIME_future, 0:6] - II_fit[SIMTIME_fit, 0:6])
      total_inf_future_6_17[i] <- sum(II_future[SIMTIME_future, 7:18] - II_fit[SIMTIME_fit, 7:18])
      total_inf_future_18_64[i] <- sum(II_future[SIMTIME_future, 19:65] - II_fit[SIMTIME_fit, 19:65])
      total_inf_future_65_84[i] <- sum(II_future[SIMTIME_future, 66:85] - II_fit[SIMTIME_fit, 66:85])
      total_inf_future_85_[i] <- sum(II_future[SIMTIME_future, 86:101] - II_fit[SIMTIME_fit, 86:101])
      
      # Calculating total hospitalisations in each age group
      
      total_hosp_future_0_5[i] <- sum(HH_future[SIMTIME_future, 0:6] - HH_fit[SIMTIME_fit, 0:6])
      total_hosp_future_6_17[i] <- sum(HH_future[SIMTIME_future, 7:18] - HH_fit[SIMTIME_fit, 7:18])
      total_hosp_future_18_64[i] <- sum(HH_future[SIMTIME_future, 19:65] - HH_fit[SIMTIME_fit, 19:65])
      total_hosp_future_65_84[i] <- sum(HH_future[SIMTIME_future, 66:85] - HH_fit[SIMTIME_fit, 66:85])
      total_hosp_future_85_[i] <- sum(HH_future[SIMTIME_future, 86:101] - HH_fit[SIMTIME_fit, 86:101])
      
      # Calculating total deaths in each age group
      
      total_deaths_future_0_5[i] <- sum(D_future[SIMTIME_future, 0:6] - D_fit[SIMTIME_fit, 0:6])
      total_deaths_future_6_17[i] <- sum(D_future[SIMTIME_future, 7:18] - D_fit[SIMTIME_fit, 7:18])
      total_deaths_future_18_64[i] <- sum(D_future[SIMTIME_future, 19:65] - D_fit[SIMTIME_fit, 19:65])
      total_deaths_future_65_84[i] <- sum(D_future[SIMTIME_future, 66:85] - D_fit[SIMTIME_fit, 66:85])
      total_deaths_future_85_[i] <- sum(D_future[SIMTIME_future, 86:101] - D_fit[SIMTIME_fit, 86:101])
      
      # Calculating total ICU admissions in each age group
      
      total_icu_future_0_5[i] <- sum(icu_rate[1:6] * (HH_future[SIMTIME_future, 1:6] - HH_fit[SIMTIME_fit, 1:6]))
      total_icu_future_6_17[i] <- sum(icu_rate[7:18] * (HH_future[SIMTIME_future, 7:18] - HH_fit[SIMTIME_fit, 7:18]))
      total_icu_future_18_64[i] <- sum(icu_rate[19:65] * (HH_future[SIMTIME_future, 19:65] - HH_fit[SIMTIME_fit, 19:65]))
      total_icu_future_65_84[i] <- sum(icu_rate[66:85] * (HH_future[SIMTIME_future, 66:85] - HH_fit[SIMTIME_fit, 66:85]))
      total_icu_future_85_[i] <- sum(icu_rate[86:101] * (HH_future[SIMTIME_future, 86:101] - HH_fit[SIMTIME_fit, 86:101]))
      
    } else {
      
      out <-  out_fit
      
    }
    
    # Label the output column for those in hospital and those that have died, separating into age groups
    
    HH <- out[(20 * n_age_groups + 2):(21 * n_age_groups + 1)]
    HH_0_5 <- out[, (20 * n_age_groups + 2):(20 * n_age_groups + 7)]
    HH_6_17 <- out[, (20 * n_age_groups + 8):(20 * n_age_groups + 19)]
    HH_18_64 <- out[, (20 * n_age_groups + 20):(20 * n_age_groups + 66)]
    HH_65_84 <- out[, (20 * n_age_groups + 67):(20 * n_age_groups + 86)]
    HH_85_ <- out[, (20 * n_age_groups + 87):(21 * n_age_groups + 1)]
    
    DD <- out[(6 * n_age_groups + 2):(7 * n_age_groups + 1)]
    DD_0_5 <- out[, (6 * n_age_groups + 2):(6 * n_age_groups + 7)]
    DD_6_17 <- out[, (6 * n_age_groups + 8):(6 * n_age_groups + 19)]
    DD_18_64 <- out[, (6 * n_age_groups + 20):(6 * n_age_groups + 66)]
    DD_65_84 <- out[, (6 * n_age_groups + 67):(6 * n_age_groups + 86)]
    DD_85_ <- out[, (6 * n_age_groups + 87):(7 * n_age_groups + 1)]
    
    names(HH_0_5) <- names(DD_0_5) <- seq(0, 5)
    names(HH_6_17) <- names(DD_6_17) <- seq(6, 17)
    names(HH_18_64) <- names(DD_18_64) <- seq(18, 64)
    names(HH_65_84) <- names(DD_65_84) <- seq(65, 84)
    names(HH_85_) <- names(DD_85_) <- seq(85, 100)
    
    # Calculating cumulative hospital admissions and deaths for all age groups and incident hospital admissions and deaths per age group
    
    cum_hosp[, i] <- rowSums(HH)
    inc_hosp[, i] <- diff(rowSums(HH))
    
    inc_hosp_0_5[, i] <- diff(rowSums(HH_0_5))
    inc_hosp_6_17[, i] <- diff(rowSums(HH_6_17))
    inc_hosp_18_64[, i] <- diff(rowSums(HH_18_64))
    inc_hosp_65_84[, i] <- diff(rowSums(HH_65_84))
    inc_hosp_85_[, i] <- diff(rowSums(HH_85_))
    
    cum_deaths[, i] <- rowSums(DD)
    inc_deaths[, i] <- diff(rowSums(DD))
    
    inc_deaths_0_5[, i] <- diff(rowSums(DD_0_5))
    inc_deaths_6_17[, i] <- diff(rowSums(DD_6_17))
    inc_deaths_18_64[, i] <- diff(rowSums(DD_18_64))
    inc_deaths_65_84[, i] <- diff(rowSums(DD_65_84))
    inc_deaths_85_[, i] <- diff(rowSums(DD_85_))
    
  }
  
  print(mean(burn_in_deaths))
  
  # Calculating total number of hospital admissions and deaths per age group
  
  total_hosp <- colSums(inc_hosp)
  total_hosp_0_5 <- colSums(inc_hosp_0_5)
  total_hosp_6_17 <- colSums(inc_hosp_6_17)
  total_hosp_18_64 <- colSums(inc_hosp_18_64)
  total_hosp_65_84 <- colSums(inc_hosp_65_84)
  total_hosp_85_ <- colSums(inc_hosp_85_)
  
  total_deaths <- colSums(inc_deaths)
  total_deaths_0_5 <- colSums(inc_deaths_0_5)
  total_deaths_6_17 <- colSums(inc_deaths_6_17)
  total_deaths_18_64 <- colSums(inc_deaths_18_64)
  total_deaths_65_84 <- colSums(inc_deaths_65_84)
  total_deaths_85_ <- colSums(inc_deaths_85_)
  
  # creating a data frame containing the age distribution of hospital admissions and another for deaths for the data and model runs
  
  total_hosp_all_ages <- cbind(total_hosp_0_5, total_hosp_6_17, total_hosp_18_64, total_hosp_65_84, total_hosp_85_)
  total_deaths_all_ages <- cbind(total_deaths_0_5, total_deaths_6_17, total_deaths_18_64, total_deaths_65_84, total_deaths_85_)
  hosp_age_distribution <- rbind(hosp_data_total[2:6], total_hosp_all_ages)
  hosp_age_distribution <- data.frame(hosp_age_distribution)
  deaths_age_distribution <- rbind(deaths_data_total[2:6], total_deaths_all_ages)
  deaths_age_distribution <- data.frame(deaths_age_distribution)
  
  row_names_parameter_set <- seq(1, ncol(par_fit), 1)
  row_names_parameter_set_rep <- rep(c("Data", row_names_parameter_set), 5)
  colnames(hosp_age_distribution) <- colnames(deaths_age_distribution)  <-c("AgeGroup1", "AgeGroup2", "AgeGroup3", "AgeGroup4", "AgeGroup5")
  
  hosp_age_distribution_tidy <- hosp_age_distribution |>
    pivot_longer(cols = AgeGroup1:AgeGroup5,
                 names_to = "AgeGroup",
                 values_to = "Value")
  hosp_age_distribution_tidy <- cbind(row_names_parameter_set_rep, hosp_age_distribution_tidy)
  
  deaths_age_distribution_tidy <- deaths_age_distribution |>
    pivot_longer(cols = AgeGroup1:AgeGroup5,
                 names_to = "AgeGroup",
                 values_to = "Value")
  deaths_age_distribution_tidy <- cbind(row_names_parameter_set_rep, deaths_age_distribution_tidy)
  
  colnames(deaths_age_distribution_tidy) <- colnames(hosp_age_distribution_tidy) <- c("Source", "AgeGroup", "Value")
  
  # Summing every seven rows of the incident hospital admissions and deaths to calculate weekly admissions and deaths per age group
  
  if (SIMTIME%%7 == 0) {
    n_weeks <- (SIMTIME - 7) / 7
  } else {
    n_weeks <- floor(SIMTIME / 7)
  }
  
  inc_hosp_weekly <- rowsum(inc_hosp[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_hosp_0_5_weekly <- rowsum(inc_hosp_0_5[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_hosp_6_17_weekly <- rowsum(inc_hosp_6_17[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_hosp_18_64_weekly <- rowsum(inc_hosp_18_64[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_hosp_65_84_weekly <- rowsum(inc_hosp_65_84[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_hosp_85_weekly <- rowsum(inc_hosp_85_[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  
  inc_deaths_weekly <- rowsum(inc_deaths[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_deaths_0_5_weekly <- rowsum(inc_deaths_0_5[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_deaths_6_17_weekly <- rowsum(inc_deaths_6_17[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_deaths_18_64_weekly <- rowsum(inc_deaths_18_64[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_deaths_65_84_weekly <- rowsum(inc_deaths_65_84[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  inc_deaths_85_weekly <- rowsum(inc_deaths_85_[1:(n_weeks * 7), ], rep(1:n_weeks, each = 7))
  
  mean_weekly_hosp <- mean(inc_hosp_weekly)
  mean_weekly_deaths <- mean(inc_deaths_weekly)
  
  # Extending the matrix to have the same dimensions as the hospital admissions and death matrices for the total model run
  
  empty_matrix <- matrix("N/A", nrow = n_weeks - nrow(hosp_data_weekly), ncol = ncol(hosp_data_weekly))
  colnames(empty_matrix) <- c("Date", "Total", "0-5", "6-17", "18-64", "65-84", "85+")
  hosp_data_weekly_extended <- rbind(hosp_data_weekly, empty_matrix)
  deaths_data_weekly_extended <- rbind(deaths_data_weekly, empty_matrix)
  
  suppressWarnings({
    hosp_data_all_ages_weekly <- as.numeric(hosp_data_weekly_extended[, 2])
    hosp_data_0_5_weekly <- as.numeric(hosp_data_weekly_extended[, 3])
    hosp_data_6_17_weekly <- as.numeric(hosp_data_weekly_extended[, 4])
    hosp_data_18_64_weekly <- as.numeric(hosp_data_weekly_extended[, 5])
    hosp_data_65_84_weekly <- as.numeric(hosp_data_weekly_extended[, 6])
    hosp_data_85_weekly <- as.numeric(hosp_data_weekly_extended[, 7])
    deaths_data_all_ages_weekly <- as.numeric(deaths_data_weekly_extended[, 2])
    deaths_data_0_5_weekly <- as.numeric(deaths_data_weekly_extended[, 3])
    deaths_data_6_17_weekly <- as.numeric(deaths_data_weekly_extended[, 4])
    deaths_data_18_64_weekly <- as.numeric(deaths_data_weekly_extended[, 5])
    deaths_data_65_84_weekly <- as.numeric(deaths_data_weekly_extended[, 6])
    deaths_data_85_weekly <- as.numeric(deaths_data_weekly_extended[, 7])
  })
  
  # Appending the hospital admissions and deaths data to the model hospital admissions and deaths to create one long vector for each age group
  
  if (SIMTIME%%7 == 0) {
    parameter_set_hosp <- rep(paste("Parameter set", "1"), floor((SIMTIME) / 7) - 1)
    if (ncol(par_fit) > 1) {
      for (i in 2:ncol(par_fit)) {
        parameter_set_hosp <- append(parameter_set_hosp, rep(paste("Parameter set", i), floor((SIMTIME) / 7) - 1))
      }
    }
    parameter_set_hosp <- append(parameter_set_hosp, rep("England data", n_weeks))
  } else {
    parameter_set_hosp <- append(parameter_set, rep("England data", n_weeks))
  }
  
  hosp_weekly <- append(round(inc_hosp_weekly), hosp_data_all_ages_weekly)
  hosp_0_5_weekly <- append(round(inc_hosp_0_5_weekly), hosp_data_0_5_weekly)
  hosp_6_17_weekly <- append(round(inc_hosp_6_17_weekly), hosp_data_6_17_weekly)
  hosp_18_64_weekly <- append(round(inc_hosp_18_64_weekly), hosp_data_18_64_weekly)
  hosp_65_84_weekly <- append(round(inc_hosp_65_84_weekly), hosp_data_65_84_weekly)
  hosp_85_weekly <- append(round(inc_hosp_85_weekly), hosp_data_85_weekly)
  
  hosp_date <- rep(dates[1:n_weeks], ncol(par_fit) + 1)
  hosp_df_weekly <- data.frame(hosp_date,
                               hosp_weekly <- hosp_weekly,
                               parameter_set_hosp <- parameter_set_hosp)
  
  deaths_weekly <- append(round(inc_deaths_weekly), deaths_data_all_ages_weekly)
  deaths_0_5_weekly <- append(round(inc_deaths_0_5_weekly), deaths_data_0_5_weekly)
  deaths_6_17_weekly <- append(round(inc_deaths_6_17_weekly), deaths_data_6_17_weekly)
  deaths_18_64_weekly <- append(round(inc_deaths_18_64_weekly), deaths_data_18_64_weekly)
  deaths_65_84_weekly <- append(round(inc_deaths_65_84_weekly), deaths_data_65_84_weekly)
  deaths_85_weekly <- append(round(inc_deaths_85_weekly), deaths_data_85_weekly)
  
  deaths_date <- rep(dates[1:n_weeks], ncol(par_fit) + 1)
  deaths_df_weekly <- data.frame(deaths_date,
                                 deaths_weekly <- deaths_weekly,
                                 parameter_set_deaths <- parameter_set_hosp)
  
  # Reordering the dataframes to make sure "England data" is last
  
  hosp_df_weekly$parameter_set_hosp <- fct_relevel(hosp_df_weekly$parameter_set_hosp, "England data", after = Inf)
  deaths_df_weekly$parameter_set_deaths <- fct_relevel(deaths_df_weekly$parameter_set_deaths, "England data", after = Inf)
  
  # Calculating the weekly mean and 95th percentiles of the hospital admissions and deaths for each parameter set
  
  inc_hosp_weekly_means <- rowMeans(inc_hosp_weekly)
  inc_hosp_weekly_95_percentiles <- apply(inc_hosp_weekly, 1, quantile, c(0.025, 0.975))
  
  inc_deaths_weekly_means <- rowMeans(inc_deaths_weekly)
  inc_deaths_weekly_95_percentiles <- apply(inc_deaths_weekly, 1, quantile, c(0.025, 0.975))
  
  # Appending the hospital admissions and deaths data to the model hospital admissions and deaths mean and 95% percentiles 
  # to create one long vector
  
  hosp_weekly_means <- append(round(inc_hosp_weekly_means), hosp_data_all_ages_weekly)
  
  hosp_df_weekly_summaries <- data.frame(dates,
                                         inc_hosp_weekly_means <- inc_hosp_weekly_means,
                                         inc_hosp_weekly_95_percentile_lower <- inc_hosp_weekly_95_percentiles[1, ],
                                         inc_hosp_weekly_95_percentile_upper <- inc_hosp_weekly_95_percentiles[2, ],
                                         hosp_data_all_ages_weekly <- hosp_data_all_ages_weekly)
  
  a3 <- ggplot(hosp_df_weekly_summaries, aes(x = dates , y = inc_hosp_weekly_means)) +
    labs(title = "(a) Hospital admissions",
         x = "Date",
         y = "Weekly hospital admissions") +
    geom_line(data = hosp_df_weekly_summaries, aes(x = dates , y = hosp_data_all_ages_weekly), colour = "red", linewidth = 0.8) +
    geom_ribbon(data = hosp_df_weekly_summaries,
                aes(x = dates, min = inc_hosp_weekly_95_percentile_lower, max = inc_hosp_weekly_95_percentile_upper),
                fill = "grey", alpha = 0.5) +
    geom_line(data = hosp_df_weekly_summaries, aes(x = dates , y = inc_hosp_weekly_means), linewidth = 0.8) +
    geom_hline(yintercept = max(hosp_data_weekly$Total), colour = "red", linetype = "dashed", linewidth = 1) +
    annotate("label", x = as.Date(min(deaths_date) + 20), y = max(hosp_data_weekly$Total), label = "Max", size = 6, colour = "black") +
    geom_hline(yintercept = min(hosp_data_weekly$Total), colour = "red", linetype = "dashed", linewidth = 1) +
    annotate("label", x = as.Date(min(deaths_date) + 20), y = min(hosp_data_weekly$Total), label = "Min", size = 6, colour = "black") +
    scale_x_continuous(breaks = seq(min(deaths_date), max(deaths_date), by = 30.41), 
                       labels = c("Jan", "Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
                                  "Oct", "Nov", "Dec"), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 16))
  
  deaths_weekly_means <- append(round(inc_deaths_weekly_means), deaths_data_all_ages_weekly)
  
  deaths_df_weekly_summaries <- data.frame(dates,
                                           inc_deaths_weekly_means <- inc_deaths_weekly_means,
                                           inc_deaths_weekly_95_percentile_lower <- inc_deaths_weekly_95_percentiles[1, ],
                                           inc_deaths_weekly_95_percentile_upper <- inc_deaths_weekly_95_percentiles[2, ],
                                           deaths_data_all_ages_weekly <- deaths_data_all_ages_weekly)
  
  b3 <- ggplot(deaths_df_weekly_summaries, aes(x = dates , y = inc_deaths_weekly_means)) +
    labs(title = "(b) Deaths",
         x = "Date",
         y = "Weekly deaths") +
    geom_line(data = deaths_df_weekly_summaries, aes(x = dates , y = deaths_data_all_ages_weekly), colour = "red", linewidth = 0.8) +
    geom_ribbon(data = deaths_df_weekly_summaries,
                aes(x = dates, min = inc_deaths_weekly_95_percentile_lower, max = inc_deaths_weekly_95_percentile_upper),
                fill = "grey", alpha = 0.5) +
    geom_line(data = deaths_df_weekly_summaries, aes(x = dates , y = inc_deaths_weekly_means), linewidth = 0.8) +
    geom_hline(yintercept = max(deaths_data_weekly$Total[4:52]), colour = "red", linetype = "dashed", linewidth = 1) +
    annotate("label", x = as.Date(min(deaths_date) + 20), y = max(deaths_data_weekly$Total[4:52]), label = "Max.", size = 6, colour = "black") +
    geom_hline(yintercept = min(deaths_data_weekly$Total), colour = "red", linetype = "dashed", linewidth = 1) +
    annotate("label", x = as.Date(min(deaths_date) + 20), y = min(deaths_data_weekly$Total), label = "Min.", size = 6, colour = "black") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16)) +
    scale_x_continuous(breaks = seq(min(deaths_date), max(deaths_date), by = 30.41),
                       labels = c("Jan", "Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
                                  "Oct", "Nov", "Dec"),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, NA))
  
  c2 <- grid.arrange(a3, b3, nrow = 2)
  
  print(a3)
  
  print(b3)
  
  ggsave(path = "Results/Model fitting/", 
         filename = paste0(as.character(Sys.Date()), "Plot 10X. Weekly hosp and deaths in 2023 (mean and CI).jpeg"), 
         plot = c2, width = 12, height = 10)
  
  # Calculating max. weekly hospital admissions after 2030
  
  if (n_years_future == 100) {
    max_weekly_hosp_post_2030 <- max(rowMeans(inc_hosp_weekly[(313:nrow(inc_hosp_weekly)),]))
    mean_weekly_hosp_post_2030 <- mean(rowMeans(inc_hosp_weekly[(313:nrow(inc_hosp_weekly)),]))
  } else {
    max_weekly_hosp_post_2030 <- NA
    mean_weekly_hosp_post_2030 <- NA
  }
  
  # Creating a dataframe of the number in each disease state each year in each age group to calculate costs/QALYs
  
  df_yearly_inf_hosp_deaths_lc <- list(n_healthy_yearly = n_healthy_yearly,
                                       n_symptomatic_yearly = n_symptomatic_yearly,
                                       n_hospitalised_yearly = n_hospitalised_yearly,
                                       n_dead_yearly = n_dead_yearly,
                                       n_recovered_past_year_yearly = n_recovered_past_year_yearly,
                                       n_recovered_vaccinated_past_year_yearly = n_recovered_vaccinated_past_year_yearly,
                                       n_recovered_waned_past_year_yearly = n_recovered_waned_past_year_yearly,
                                       inc_inf_yearly = inc_inf_yearly,
                                       inc_hosp_yearly = inc_hosp_yearly,
                                       inc_deaths_yearly = inc_deaths_yearly,
                                       inc_vacc_yearly = inc_vacc_yearly,
                                       inc_symp_inf_yearly = inc_symp_inf_yearly,
                                       inc_rec_yearly = inc_rec_yearly,
                                       max_weekly_hosp_post_2030 = max_weekly_hosp_post_2030,
                                       mean_weekly_hosp_post_2030 = mean_weekly_hosp_post_2030)
  
  df_hosp_deaths_weekly_summaries <- data.frame(hosp_df_weekly_summaries = hosp_df_weekly_summaries, deaths_df_weekly_summaries = deaths_df_weekly_summaries)
  
  output <- list(df_yearly_inf_hosp_deaths_lc, df_hosp_deaths_weekly_summaries)
  
  return(output)
  
}


# Set up prior distributions

model_priors_29 <- list(c("unif", 0.25, 0.55), c("unif", 0.3, 0.99999), c("unif", 0.3, 0.99999), c("unif", 0.015, 0.036), 
                        c("unif", 0.005, 0.011), c("unif", 0.03, 0.07), c("unif", 0.27, 0.61), c("unif", 1, 1.70), 
                        c("unif", 0.22, 0.5), c("unif", 2, 7), c("unif", 1, 7), c("unif", 0.5, 3.5), c("unif", 0.001, 0.8))

# Set target summary statistic

data_2023 <- c(2227, 684, 11701, 20488, 10988, 10, 12, 1071, 6845, 7246)
names(data_2023) <- c("Hosp 0-5", "Hosp 6-17", "Hosp 18-64", "Hosp 65-84", "Hosp 85+", 
                      "Deaths 0-5", "Deaths 6-17", "Deaths 18-64", "Deaths 65-84", "Deaths 85+")

# Using ABC sequential

print ("Fitting 10000 parameter sets using model priors_29, with alpha = default and p_acc_min = 0.025.")

abc_seq <- ABC_sequential(method = "Lenormand", model = covid_model_run_easy_ABC_set_seed, prior = model_priors_29,
                          nb_simul = 200, summary_stat_target = data_2023, n_cluster = 24, p_acc_min = 0.025,
                          use_seed = TRUE, seed_count = 1)

save.image("Results/Model fitting/Model fit 1X.RData")

retained_params <- abc_seq$param
sum_stats <- abc_seq$stats

retained_params <- retained_params[,c(1, 2, 13, 3:12)]

param_names <- colnames(retained_params) <- c("Beta", "Reduction in infection rate (boosted immunity)", 
                                              "Residual immunity infection rate scaling parameter",
                                              "Relative infectiousness of an asymptomatic infection", 
                                              "Hospitalisation rate age scaling 0-5", "Hospitalisation rate age scaling 6-17", 
                                              "Hospitalisation rate age scaling 18-64", "Hospitalisation rate age scaling 65-84", 
                                              "Hospitalisation rate age scaling 85+", "Mortality rate age scaling 0-5", 
                                              "Mortality rate age scaling 6-17", "Mortality rate age scaling 18-64", 
                                              "Mortality rate age scaling 65-84")
colnames(sum_stats) <- names(data_2023)

# Plotting posterior distributions

histogram_data <- data.frame(retained_params) |>
  pivot_longer(cols = (1:ncol(retained_params)),
               names_to = "Parameter",
               values_to = "Value")
histogram_data$Number <- rep(1:13, each = nrow(retained_params))
hist_labels <- c(`1` = "(a) Beta", 
                 `2` = "(b) Reduction in infection
rate (boosted immunity)", 
                 `3` = "(c) Residual immunity
infection rate
scaling parameter",
                 `4` = "(d) Relative infectiousness
of an
asymptomatic infection", 
                 `5` = "(e) Hospitalisation rate
age scaling 0-5", 
                 `6` = "(f) Hospitalisation rate
age scaling 6-17", 
                 `7` = "(g) Hospitalisation rate
age scaling 18-64", 
                 `8` = "(h) Hospitalisation rate
age scaling 65-84", 
                 `9` = "(i) Hospitalisation rate
age scaling 85+", 
                 `10` = "(j) Mortality rate
age scaling 0-5", 
                 `11` = "(k) Mortality rate
age scaling 6-17", 
                 `12` = "(l) Mortality rate
age scaling 18-64", 
                 `13` = "(m) Mortality rate
age scaling 65-84")
histograms <- ggplot(histogram_data, aes(Value, after_stat(ndensity))) + 
  theme_light(base_size = 12) +
  geom_histogram(bins = 50, alpha = 0.9, col = "white", fill = "navy") + 
  facet_wrap(~Number, scales = 'free_x', nrow = 5, labeller = as_labeller(hist_labels)) +
  labs(y = "Density",
       x = "Parameter value") +
  theme(strip.text.x = element_text(size = 10, colour = "black"),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 12))

print(histograms)

ggsave(path = "Results/Model fitting/", 
       filename = paste0(as.character(Sys.Date()), "Plot 1X. Histograms.jpeg"), 
       plot = histograms, height = 10)

# Creating a pairs plot

pairs_plot <- ggpairs(retained_params, diag = list(continuous = "blankDiag"), lower = list(continuous = "density"), 
                      upper = "blank", switch = "y", columnLabels = 1:13) +
  labs(caption = "") +
  theme(plot.caption = element_text(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.text = element_text(size = 11))

print(pairs_plot)

ggsave(path = "Results/Model fitting/", 
       filename = paste0(as.character(Sys.Date()), "Plot 2X. Pairs.jpeg"), 
       plot = pairs_plot)

# Plotting goodness of fit

GOF <- apply(sum_stats, 1, function (x) sum(abs(data_2023 - x)))

density_plots <- list()
for (i in 1:13) {
  x_params <- as.numeric(retained_params[,i])
  data <- data.frame(x_params, GOF)
  density_plots[[i]] <-   ggplot(data, aes(x = x_params, y = GOF) ) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    theme(legend.position="none") +
    labs(x = param_names[i])
}
density_plots <- wrap_plots(density_plots, ncol = 5)

ggsave(path = "Results/Model fitting/", 
       filename = paste0(as.character(Sys.Date()), "Plot 3.X Density GOF plots.pdf"), 
       plot = density_plots)

# Summarising results

summary <- data.frame(nrow = 7, ncol = 13)

for (i in 1:13) {
  summary[1, i] <- signif(mean(retained_params[, i]), 2)
  summary[2, i] <- median(retained_params[, i])
  summary[3, i] <- sd(retained_params[, i])
  summary[4, i] <- min(retained_params[, i])
  summary[5, i] <- max(retained_params[, i])
  summary[6, i] <- signif(quantile(retained_params[, i], probs = c(0.025)), 2)
  summary[7, i] <- signif(quantile(retained_params[, i], probs = c(0.975)), 2)
}

summary_tidy <- summary |>
  pivot_longer(cols = 1:13, 
               names_to = "Parameter", 
               values_to = "Value")

abc_seq_results <- append(abc_seq, append(list(GOF), list(summary_tidy)))
names(abc_seq_results) <- c(names(abc_seq), "GOF", "Parameter summaries")
write.xlsx(abc_seq_results, paste0("Results/Model fitting/", as.character(Sys.Date()), "Model fit results 1X.xlsx"))

# Run model using parameters from the fitting

n_samples_fit <- nrow(abc_seq$param)

beta_fit <- abc_seq$param[, 1]
reinf_fit <- abc_seq$param[, 2]
tau_fit <- abc_seq$param[, 3]
severity_scaling_0_5_fit <- abc_seq$param[, 4]
severity_scaling_6_17_fit <- abc_seq$param[, 5]
severity_scaling_18_64_fit <- abc_seq$param[, 6]
severity_scaling_65_84_fit <- abc_seq$param[, 7]
severity_scaling_85_fit <- abc_seq$param[, 8]
hosp_mortality_scaling_0_5_fit <- abc_seq$param[, 9]
hosp_mortality_scaling_6_17_fit <- abc_seq$param[, 10]
hosp_mortality_scaling_18_64_fit <- abc_seq$param[, 11]
hosp_mortality_scaling_65_84_fit <- abc_seq$param[, 12]
waned_ve_against_infection_fit <- abc_seq$param[, 13]

par_fit <- rbind(beta_fit, reinf_fit, tau_fit, severity_scaling_0_5_fit, severity_scaling_6_17_fit, severity_scaling_18_64_fit, 
                 severity_scaling_65_84_fit, severity_scaling_85_fit, hosp_mortality_scaling_0_5_fit, hosp_mortality_scaling_6_17_fit, 
                 hosp_mortality_scaling_18_64_fit, hosp_mortality_scaling_65_84_fit, waned_ve_against_infection_fit)

rownames(par_fit) <- c("Beta", "Protection against reinfection", "Relative transmission from asymptomatic infections", 
                       "Infection severity age scaling 0-5", "Infection severity age scaling 6-17", "Infection severity age scaling 18-64", 
                       "Infection severity age scaling 65-84", "Infection severity age scaling 85+", "Hospital mortality age scaling 0-5", 
                       "Hospital mortality age scaling 6-17", "Hospital mortality age scaling 18-64", "Hospital mortality age scaling 65-84",
                       "Waned VE against infection scaling")
colnames(par_fit) <- seq(1:n_samples_fit)

results <- covid_model_run_fixed_params_age_dist(par_fit = par_fit)

hosp_age_distribution_tidy <- results[[1]]
deaths_age_distribution_tidy <- results[[2]]

data_2023 <- c(2227, 684, 11701, 20488, 10988, 10, 12, 1071, 6845, 7246)
names(data_2023) <- c("Hosp 0-5", "Hosp 6-17", "Hosp 18-64", "Hosp 65-84", "Hosp 85+", 
                      "Deaths 0-5", "Deaths 6-17", "Deaths 18-64", "Deaths 65-84", "Deaths 85+")

# Plotting mean and confidence intervals for summary stats against data

model_fit_sum_stats_means_hosp <- colMeans(sum_stats[, 1:5])
model_fit_sum_stats_max_hosp <- apply(sum_stats[, 1:5], 2, max)
model_fit_sum_stats_min_hosp <- apply(sum_stats[, 1:5], 2, min)

hosp_age_distribution_fit_mean <- data.frame(c(rbind(model_fit_sum_stats_means_hosp, data_2023[1:5])))
hosp_age_distribution_fit_max <- c(rbind(model_fit_sum_stats_max_hosp, data_2023[1:5]))
hosp_age_distribution_fit_min <- c(rbind(model_fit_sum_stats_min_hosp, data_2023[1:5]))

hosp_age_distribution_fit <- (cbind(rep(c("Model", "Data"), 5), rep(c("HospAgeGroup1", "HospAgeGroup2", 
                                                                      "HospAgeGroup3", "HospAgeGroup4", 
                                                                      "HospAgeGroup5"), 
                                                                    each = 2), 
                                    hosp_age_distribution_fit_mean))

model_fit_sum_stats_means_deaths <- colMeans(sum_stats[, 6:10])
model_fit_sum_stats_max_deaths <- apply(sum_stats[, 6:10], 2, max)
model_fit_sum_stats_min_deaths <- apply(sum_stats[, 6:10], 2, min)

deaths_age_distribution_fit_mean <- data.frame(c(rbind(model_fit_sum_stats_means_deaths, data_2023[6:10])))
deaths_age_distribution_fit_max <- c(rbind(model_fit_sum_stats_max_deaths, data_2023[6:10]))
deaths_age_distribution_fit_min <- c(rbind(model_fit_sum_stats_min_deaths, data_2023[6:10]))

deaths_age_distribution_fit <- (cbind(rep(c("Model", "Data"), 5), rep(c("DeathsAgeGroup1", "DeathsAgeGroup2", 
                                                                        "DeathsAgeGroup3", "DeathsAgeGroup4", 
                                                                        "DeathsAgeGroup5"), 
                                                                      each = 2), 
                                      deaths_age_distribution_fit_mean))

colnames(hosp_age_distribution_fit) <- colnames(deaths_age_distribution_fit) <- c("Source", "Statistic", "Value")

hosp_age_distribution_fit <- mutate(hosp_age_distribution_fit, 
                                    Statistic = factor(Statistic, levels = c("HospAgeGroup1", "HospAgeGroup2",
                                                                             "HospAgeGroup3", "HospAgeGroup4",
                                                                             "HospAgeGroup5")))

deaths_age_distribution_fit <- mutate(deaths_age_distribution_fit, 
                                      Statistic = factor(Statistic, levels = c("DeathsAgeGroup1", "DeathsAgeGroup2", 
                                                                               "DeathsAgeGroup3", "DeathsAgeGroup4", 
                                                                               "DeathsAgeGroup5")))

hosp_sum_stats_plot_mean_CI <- ggplot(hosp_age_distribution_fit, aes(x = Statistic, y = Value, group = Source, color = Source)) + 
  geom_point(size = 1) +
  geom_pointrange(aes(ymin = hosp_age_distribution_fit_min, ymax = hosp_age_distribution_fit_max), linewidth = 0.6) +
  scale_x_discrete(labels = c("0-5", "6-17", "18-64", "65-84", "85+")) +
  scale_y_continuous(limits = c(0, NA)) +
  ggtitle("(a) Hospital admissions in 2023") +
  xlab("Age group") + ylab("Total hospital admissions") +
  theme_light(base_size = 12) +
  theme(legend.position = c(0.15, 0.85))

deaths_sum_stats_plot_mean_CI <- ggplot(deaths_age_distribution_fit, aes(x = Statistic, y = Value, group = Source, color = Source)) + 
  geom_point(size = 1) +
  geom_pointrange(aes(ymin = deaths_age_distribution_fit_min, ymax = deaths_age_distribution_fit_max), linewidth = 0.6) +
  scale_x_discrete(labels = c("0-5", "6-17", "18-64", "65-84", "85+")) +
  scale_y_continuous(limits = c(0, NA)) +
  ggtitle("(b) Deaths in 2023") +
  xlab("Age group") + ylab("Total deaths") +
  theme_light(base_size = 12) +
  theme(legend.position = "none")

sum_stats_plot_mean_CI <- hosp_sum_stats_plot_mean_CI + deaths_sum_stats_plot_mean_CI + plot_layout(nrow = 1, axes = "collect")
sum_stats_plot_mean_CI
ggsave(path = "Results/Model fitting/", 
       filename = paste0(as.character(Sys.Date()), "Plot 8X. Sum stats mean and CI.jpeg"), 
       plot = sum_stats_plot_mean_CI, height = 5, width = 8)

# Drawing boxplots

sum_stats_renamed <- sum_stats
colnames(sum_stats_renamed) <- c("Hosp 1", "Hosp 2", "Hosp 3", "Hosp 4", "Hosp 5", "Deaths 1", "Deaths 2", "Deaths 3", "Deaths 4", "Deaths 5")

sum_stats_hosp_longer <- data.frame(sum_stats_renamed[, 1:5]) |> 
  pivot_longer(cols = Hosp.1:Hosp.5,
               names_to = "Statistic",
               values_to = "Value")

hosp_data <- data.frame(colnames(data.frame(sum_stats_renamed[, 1:5])), data_2023[1:5])
colnames(hosp_data) <- c("Age group", "Value")

hosp_boxplot <- ggplot(sum_stats_hosp_longer, aes(x = Statistic, y = Value)) + 
  geom_boxplot(colour = hue_pal()(2)[2]) +
  geom_point(data = hosp_data, (aes(x = `Age group`, y = `Value`)), shape = 4, size = 5, stroke = 1, colour = hue_pal()(2)[1]) +
  scale_x_discrete(labels = c("0-5", "6-17", "18-64", "65-84", "85+")) +
  ggtitle("(a) Hospital admissions in 2023") +
  xlab("Age group") + ylab("Total hospital admissions") +
  theme_light(base_size = 12)

sum_stats_deaths_longer <- data.frame(sum_stats_renamed[, 6:10]) |> 
  pivot_longer(cols = Deaths.1:Deaths.5,
               names_to = "Statistic",
               values_to = "Value")

deaths_data <- data.frame(colnames(data.frame(sum_stats_renamed[, 6:10])), data_2023[6:10])
colnames(deaths_data) <- c("Age group", "Value")

deaths_boxplot <- ggplot(sum_stats_deaths_longer, aes(x = Statistic, y = Value)) + 
  geom_boxplot(colour = hue_pal()(2)[2]) +
  geom_point(data = deaths_data, (aes(x = `Age group`, y = `Value`)), shape = 4, size = 5, stroke = 1, colour = hue_pal()(2)[1]) +
  scale_x_discrete(labels = c("0-5", "6-17", "18-64", "65-84", "85+")) +
  ggtitle("(b) Deaths in 2023") +
  xlab("Age group") + ylab("Total deaths") +
  theme_light(base_size = 12)

boxplot_mean_CI <- hosp_boxplot + deaths_boxplot + plot_layout(nrow = 1, axes = "collect")
boxplot_mean_CI
ggsave(path = "Results/Model fitting/", 
       filename = paste0(as.character(Sys.Date()), "Plot 8X. Sum stats boxplot.jpeg"), 
       plot = boxplot_mean_CI, height = 5, width = 8)

abc_seq_results <- append(abc_seq_results, results)
names(abc_seq_results) <- c(names(abc_seq), "GOF", "Parameter summaries", "Model run hosp age dist",
                            "Model run deaths age dist")
write.xlsx(abc_seq_results,  paste0("Results/Model fitting/", as.character(Sys.Date()), "Model fit results 2X.xlsx"))

save.image("Results/Model fitting/Model fit 2X.RData")

end_time <- Sys.time()
print(end_time-start_time)

