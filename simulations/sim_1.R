###############################################
# Author: Katarzyna Reluga                    #
###############################################
###########################
# Model based simulations #
###########################

# Generate 100 populations without outliers 
gen_pops <- generate_pop(coef_outcome_p_score = get_default_cov("outcome_p_score"),
                         coef_outcome = get_default_cov("outcome"),
                         coef_p_score = get_default_cov("p_score"),
                         coef_A_params = list(mean = 10, var = 1),
                         Ni_size  = 100,
                         m = 50,
                         no_sim = 100,
                         regression_type = "continuous", 
                         additional_params = list(frac_out = 0.00,
                                                  re = list(var_re = 3,
                                                            frac_out_re = 0.0,
                                                            var_re_out = 9, 
                                                            mean_re_out = 20),
                                                  var_re_treat  = 0.25, 
                                                  params_cont = list(var_e = 6,
                                                                     var_e_out = 150,
                                                                     mean_e_out = 20), 
                                                  disturbance = 5), 
                         start_seed = 1)

samples <- generate_sample(gen_pops, get_index = T)

formula_y = y ~ X1 + Xo1 + A + (1 + A||group)
formula_p_score = A ~ X1 + (1|group)

true_tau_adap <- matrix(0, nrow = 100, ncol = 50)
true_tau_reg <- matrix(0, nrow = 100, ncol = 50)


EBLUP_tau_adap <- matrix(0, nrow = 100, ncol = 50)
EBLUP_tau_reg <- matrix(0, nrow = 100, ncol = 50)

MQ_tau_adap <- matrix(0, nrow = 100, ncol = 50)
MQ_tau_reg <- matrix(0, nrow = 100,  ncol = 50)

RF_tau_adap <- matrix(0, nrow = 100, ncol = 50)
RF_tau_reg <- matrix(0, nrow = 100, ncol = 50)

RFt_tau_adap <- matrix(0, nrow = 100, ncol = 50)
RFt_tau_reg <- matrix(0, nrow = 100, ncol = 50)

for (i in 1:100) {
  
  data_sample <- data.frame(samples[[i]]$samp_data)
  index_sample <- samples[[i]]$index_s
  pop_y_na <- data.frame(gen_pops[[i]])[, -ncol(gen_pops[[i]])]
  
  data_out_of_sample <- gen_pops[[i]][-index_sample, ]
  
  true_tau_adap[i, ] <- calculate_tau(list(gen_pops[[i]]), type = "adaptive")[[1]]$tau
  true_tau_reg[i, ] <- calculate_tau(list(gen_pops[[i]]), type = "regular")[[1]]$tau

  # Mixed models
  ht_EBLUP_adap <- hte(formula_y, 
                       formula_p_score, 
                       type_tau = "adaptive",
                       data_sample, 
                       data_out_of_sample,
                       method_y = "EBLUP",
                       method_p_score =  "EBLUP")
  
  EBLUP_tau_adap[i, ] <- ht_EBLUP_adap$tau_hat[[1]]$tau
  
  ht_EBLUP_reg <- hte(formula_y, 
                      formula_p_score, 
                      type_tau = "regular",
                      data_sample, 
                      data_out_of_sample,
                      method_y = "EBLUP",
                      method_p_score =  "EBLUP")
  EBLUP_tau_reg[i, ] <- ht_EBLUP_reg$tau_hat[[1]]$tau
  
  # MQ 
  ht_MQ_adap <- hte(formula_y, 
                     formula_p_score, 
                     type_tau = "adaptive",
                     data_sample, 
                     data_out_of_sample,
                     method_y = "MQ",
                    method_p_score =  "MQ")
  MQ_tau_adap[i, ] <- ht_MQ_adap$tau_hat[[1]]$tau
  
  ht_MQ_reg <- hte(formula_y, 
                   formula_p_score, 
                   type_tau = "regular",
                   data_sample, 
                   data_out_of_sample,
                   method_y = "MQ",
                   formula_p_score =  "MQ")
  MQ_tau_reg[i, ] <- ht_MQ_reg$tau_hat[[1]]$tau
  
  
  # Random forest
  ht_RF_adap <- hte(formula_y, 
                     formula_p_score, 
                     type_tau = "adaptive",
                     data_sample, 
                     data_out_of_sample,
                     method_y = "RF",
                     formula_p_score =  "RF")
  RF_tau_adap[i, ] <- ht_RF_adap$tau_hat[[1]]$tau
  
  ht_RF_reg <- hte(formula_y, 
                   formula_p_score, 
                   type_tau = "regular",
                   data_sample, 
                   data_out_of_sample,
                   method_y = "RF",
                   formula_p_score =  "RF")
  RF_tau_reg[i, ] <- ht_RF_reg$tau_hat[[1]]$tau
  
  # Tuned random forest
  ht_RFt_adap <- hte(formula_y, 
                     formula_p_score, 
                     type_tau = "adaptive",
                     data_sample, 
                     data_out_of_sample,
                     method_y = "RF",
                     formula_p_score =  "RF", 
                     tune_RF = T)
  RFt_tau_adap[i, ] <- ht_RFt_adap$tau_hat[[1]]$tau
  
  ht_RFt_reg <- hte(formula_y, 
                    formula_p_score, 
                    type_tau = "regular",
                    data_sample, 
                    data_out_of_sample,
                    method_y = "RF",
                    formula_p_score =  "RF")
  RFt_tau_reg[i, ] <- ht_RFt_reg$tau_hat[[1]]$tau
  
  
  
}


ER_EBLUP_adap <- EBLUP_tau_adap - true_tau_adap
RB_EBLUP_adap <- colMeans(ER_EBLUP_adap, na.rm = TRUE)/colMeans(true_tau_adap)
RRMSE_EBLUP_adap <- sqrt(colMeans((ER_EBLUP_adap)^2, na.rm = TRUE))/colMeans(true_tau_adap)

ER_MQ_adap <- MQ_tau_adap - true_tau_adap
RB_MQ_adap <- colMeans(ER_MQ_adap, na.rm = TRUE)/colMeans(true_tau_adap)
RRMSE_MQ_adap <- sqrt(colMeans((ER_MQ_adap)^2, na.rm = TRUE))/colMeans(true_tau_adap)

ER_RF_adap <- RF_tau_adap - true_tau_adap
RB_RF_adap <- colMeans(ER_RF_adap, na.rm = TRUE)/colMeans(true_tau_adap)
RRMSE_RF_adap <- sqrt(colMeans((ER_RF_adap)^2, na.rm = TRUE))/colMeans(true_tau_adap)

ER_RFt_adap <- RFt_tau_adap - true_tau_adap
RB_RFt_adap <- colMeans(ER_RFt_adap, na.rm = TRUE)/colMeans(true_tau_adap)
RRMSE_RFt_adap <- sqrt(colMeans((ER_RFt_adap)^2, na.rm = TRUE))/colMeans(true_tau_adap)

data_RB <- data.frame(RB_EB = RB_EBLUP_adap, 
                      RB_MQ = RB_MQ_adap,
                      RB_RF = RB_RF_adap,
                      RB_RFt = RB_RFt_adap)

boxplot(data_RB)

data_RRMSE <- data.frame(RRMSE_EB = RRMSE_EBLUP_adap, 
                         RRMSE_MQ = RRMSE_MQ_adap,
                         RRMSE_RF = RRMSE_RF_adap,
                         RRMSE_RFt = RRMSE_RFt_adap)

boxplot(data_RRMSE)



ER_EBLUP_reg <- EBLUP_tau_reg - true_tau_reg
RB_EBLUP_reg <- colMeans(ER_EBLUP_reg, na.rm = TRUE)/colMeans(true_tau_reg)
RRMSE_EBLUP_reg <- sqrt(colMeans((ER_EBLUP_reg)^2, na.rm = TRUE))/colMeans(true_tau_reg)

ER_MQ_reg <- MQ_tau_reg - true_tau_reg
RB_MQ_reg <- colMeans(ER_MQ_reg, na.rm = TRUE)/colMeans(true_tau_reg)
RRMSE_MQ_reg <- sqrt(colMeans((ER_MQ_reg)^2, na.rm = TRUE))/colMeans(true_tau_reg)

ER_RF_reg <- RF_tau_reg - true_tau_reg
RB_RF_reg <- colMeans(ER_RF_reg, na.rm = TRUE)/colMeans(true_tau_reg)
RRMSE_RF_reg <- sqrt(colMeans((ER_RF_reg)^2, na.rm = TRUE))/colMeans(true_tau_reg)

ER_RFt_reg <- RFt_tau_reg - true_tau_reg
RB_RFt_reg <- colMeans(ER_RFt_reg, na.rm = TRUE)/colMeans(true_tau_reg)
RRMSE_RFt_reg <- sqrt(colMeans((ER_RFt_reg)^2, na.rm = TRUE))/colMeans(true_tau_reg)



save(tau_true, 
     ER_01, ER_02, ER_03, 
     RB_01,RB_02,RB_03,
     RRMSE_01,RRMSE_02,RRMSE_03,
     file = "Result_01_b.RData")

Direct <- ER_01+tau_true
EBLUP<- ER_02+tau_true
MQ <- ER_03+tau_true


for(oo in 1:50){
  plot(density(rnorm(length(MQ[,oo]),10,1)), col="red", main = c("area",oo))
  lines(density(MQ[,oo]),col=4,lty=4) #blue
  lines(density(tau_true[,oo])) #black
  lines(density(Direct[,oo],na.rm = T),col=2,lty=2) #red
  lines(density(EBLUP[,oo]), col=3, lty=3)#green
}


EBLUP <- ER_EBLUP_adap + true_tau_adap
MQ <- ER_MQ_adap + true_tau_adap
RF <- ER_RF_adap + true_tau_adap
RFt <- ER_RFt_adap + true_tau_adap


for(oo in 1:50){
  plot(density(rnorm(length(MQ[,oo]),10,1)), col="red", main = c("area",oo))
  lines(density(true_tau_adap[,oo])) #black
  lines(density(MQ[,oo]),col=3,lty=4) #blue
  lines(density(EBLUP[,oo]), col=4, lty=3)#green
  lines(density(RF[,oo]), col=5, lty=3)#
  lines(density(RFt[,oo]), col=6, lty=3)#
}



# Sample from population
samples <- generate_sample(gen_pops, get_index = T, ni_size = 50)
data_sample <- data.frame(samples[[1]]$samp_data)
index_sample <- samples[[1]]$index_s
gen_popsy <- data.frame(gen_pops)[,-ncol(gen_pops[[1]])]

gen_popsy1 <- data.frame(gen_pops) %>% filter(group ==1)
subpopulation <- gen_popsy1

formula_y = y ~ X1 + Xo1 + A + (1 + A||group)
formula_p_score = A ~ X1 + (1|group)
# "X1"      "Xo1"     "group"    "A"   "tau_repeat" "p_score"    "re_repeat"  "y" 
data_out_of_sample = gen_popsy[-index_sample,]
# method_y = "EBLUP"
# formula_p_score = "EBLUP"
true_tau_adap <- calculate_tau(gen_pops, type = "adaptive")
true_tau_reg <- calculate_tau(gen_pops, type = "regular")

