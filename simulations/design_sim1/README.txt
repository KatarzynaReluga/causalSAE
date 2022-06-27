README file for the results of simulations

Pdfs

Data_analysis._HTE

Model-based and direct HTE plotted together
HTE1 -- when propensity score estimated using
only data from the sample
HTE1_NoDir -- when propensity score estimated using
only data from the sample, no dir estimates
HTE2  -- propensity score estimated the whole population

Boxplots_HTE_density -- clear, from the first sim

CI_compare -- confidence intervals, simulations, prop score 
estimated using sample data
CI_compare2 -- confidence intervals, simulations, prop score 
estimated using population data

Effects, Effects2 --  true effects from the first simulations


RB_rmse... -- relative bias and rmse from the populations with prop score
estimated from the sample

RB_rmse...2 -- relative bias and rmse from the populations with prop score
estimated from the population

csv files
data_filtered_census.csv
data_filtered_survey.csv
Copied files from 
C:\Users\katar\Documents\Kasia\rok_2021_2022\Causality\Comparative_HTE\EU-SILC_Censimento_SAE_Causal_Inference\Censimento\Filtered

sim1 -- simulations with propensity score estimated using only the sample
sim12 -- simulations with propensity score estimated using only the population

RData
Data_analysis.RData

Dir_EBLUP_tau -- EBLUP and population data to compute propensity score model
Dir_EBLUP_tau2 -- EBLUP and only sample data to compute propensity score model

Dir_RF_tau -- random forest and population data to compute propensity score model
Dir_RF_tau2 -- random forest and only sample data to compute propensity score model

EBLUP_tau_adap -- model based EBLUP HTE with adaptive weights
MQ_tau_adap -- model based MQ HTE with adaptive weigths
RF_tau_adap -- model based RF HTE with adaptive weigths
RFt_tau_adap -- model based RFt HTE with adapttive wights

Design_based_sim1.RData -- results without direct effects
Design_based_sim1corr2.RData -- results witht direct effects and propensity score estimates from
the population
Design_based_sim1corr.RData -- results witht direct effects and propensity score estimates from
the sample
