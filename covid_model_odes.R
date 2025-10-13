# Added compartment for RRW calculating cumulative waned recoveries

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
