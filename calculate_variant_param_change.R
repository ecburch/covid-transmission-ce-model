
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
