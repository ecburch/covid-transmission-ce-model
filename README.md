# covid-transmission-ce-model
This is the R code for an age-structured, deterministic, compartmental model of COVID-19 transmission in England, with a 100-year time horizon and with costs and quality-adjusted life-years (QALYs) incorporated.

To run the code, download all scripts and the data folder, and set the containing folder to be the working directory in R. Create a folder within the working directory named "Results", with subfolders named "Model fitting", "Sensitivity analysis", and "Variant scenarios".

There are nine .R scripts for running the Covid model, which are described below:
1. calculate_variant_param_change.R
   
   This code contains five functions which calculate changes in model parameters due to new variants being introduced. These functions are called upon in the "covid_model_fit.R", "covid_model_run.R", "covid_model_analysis.R", "covid_model_analysis_variants_50_year_simulation_comparator_only.R", and "covid_model_analysis_variants_full_simulation.R" scripts.
2. covid_model_economic_parameters.R
   
   This code contains the parameters required for the economic evaluation, including the risk of and vaccine effectiveness against long Covid, the ICU admission rate, and all state utilities and costs. These parameters are required in the "covid_model_analysis.R", "covid_model_analysis_variants_full_simulation.R", and "covid_model_sensitivity_analysis.R" scripts.
3. covid_model_odes.R
   
   This code contains the series of ordinary differential equations which make up the Covid model. The ODE function is called upon in the "covid_model_fit.R", "covid_model_run.R", "covid_model_analysis.R", "covid_model_analysis_variants_50_year_simulation_comparator_only.R", and "covid_model_analysis_variants_full_simulation.R" scripts.
4. covid_model_fit.R
   
   This is the model fitting code which was run using the Advanced Computing Research Centre at the University of Bristol. The fitting function is set up to run in parallel on 24 clusters.
5. covid_model_run.R
    
   This code contains the function to set up and run the Covid model, given a set of input parameters obtained from the model fitting, the variant scenario, number of yearly variants, vaccination scenario, and time horizon. The function outputs the results from running the model, which are then analysed and interpreted in a separate script. This function is called upon in the "covid_model_analysis.R" and "covid_model_analysis_variants_full_simulation.R" scripts.
6. covid_model_analysis.R
    
   This is the code which runs and carries out the main analysis of the Covid model, using the functions and parameters described above. This script outputs the main results tables and figures. In the model analysis, 500 parameter sets obtained from the posterior distribution of the model fit are used, and a time horizon of 100 years.
7. covid_model_analysis_variants_50_year_simulation_comparator_only.R

   This script runs the Covid model for each of the variant scenarios for 50 years in the comparator vaccination strategy (no future vaccination) to plot the projected weekly hospital admissions and deaths in this scenario for each variant, over 50 years. 
8. covid_model_analysis_variants_full_simulation.R

   This script runs the Covid model for each of the variant scenarios for the full 100-year period over all vaccination strategies, and carries out the full analysis for each variant scenario. This script outputs the results tables and figures for the variant scenario analysis.
9. covid_model_sensitivity_analysis.R

    This script loads the output of the model run from the "covid_model_analysis.R" script and carries out the sensitivity analyses on the results. This script outputs the results tables and figures for the sensitivity analyses. The "covid_model_analysis.R" script must be run first to save the "Model run.RData" file.

This repository also contains four .xslx data files:

1. Life expectancy (England).xlsx
   
    This file contains the 2022 ONS principal projection for period and cohort expectation of life in England. [1] The 2022 projections cover the years 1981 to 2072 and are stratified into male and female life expectancy. In the model, life-years remaining are calculated beginning from the projected life expectancies in 2024 and using the 2072 projections for the years 2072 to 2123. Averages of the male and female projection figures are used, because the model is not stratified by sex. These data are required in the "covid_model_analysis.R", "covid_model_analysis_variants_full_simulation.R", and "covid_model_sensitivity_analysis.R" scripts.
2. Model fit results.xlsx
   
    This file contains the results of fitting the model to data to obtain 10,000 parameter sets. These are used for running the model future simulation.
3. Weekly hospital admissions 2023 (England).xlsx
   
    This file contains the age-stratified daily COVID-19 hospital admissions from 1 January to 31 December 2023 in England. [2] These were extracted from the datasets “Monthly COVID Publication 1 October 2022 up to 31st May 2023” and “Monthly COVID Publication 1 June 2023 up to 31 March 2024”, using the “Admissions - Number of patients admitted with COVID-19 (Last 24hrs)”.
4. Weekly deaths 2023 (England).xlsx
   
    This file contains the age-stratified weekly deaths due to COVID-19 from 1 January to 29 December 2023 in England. [3] To obtain the total number of age-stratified deaths from the data, “Weekly provisional figures on death occurrences involving COVID-19 in England and Wales by sex and age group” (Sheet 5) from the 2023 version of the dataset “Deaths registered weekly in England and Wales” were used. The age groups used here do not match the age groups used in the hospital admissions data. Therefore, for consistency, the age groups are formed by merging data from the other groups. To form the 0-5 age group, 20% of deaths in the 5-9 group are incorporated into the 0-4 group. Similarly, the 6-17 group is made up of 80% of the deaths in the 5-9 age group, all of the deaths in the 10-14 group, and 60% of the deaths in the 15-19 group. The remaining 40% of deaths in the 15-19 group is counted in the 18-64 age group. After this, the deaths for each age group are reported explicitly. These data include deaths in both England and Wales, but are also shown stratified by country in Sheet 9. The proportion of total deaths that occurred in England was calculated to be 93.7%, which was then used to scale the total number of deaths in each age group to those that occurred in England only.


References
1. Office for National Statistics. Expectation of life, principal projection, England (2022 edition). 2025. https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/expectationoflifeprincipalprojectionengland (accessed May 30, 2025).
2. NHS England. Statistics » COVID-19 Hospital Activity. Monthly publication of COVID-19 data. 2024. https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/ (accessed July 25, 2024).
3. Office for National Statistics. Deaths registered weekly in England and Wales, provisional. 2023 edition. Sheet 5: Weekly provisional figures on death occurrences involving COVID-19 in England and Wales by sex and age group, registered 2022 and 2023 n.d. https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales (accessed July 25, 2024).
