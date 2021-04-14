# Table 1 (Sample Size Calculation for Randomized Clinical Trials when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

### Individually Randomized Trials (IRTs)

# Weighting for a fully observed auxiliary categorical variable

# Sample Size Calculation
# Y continuous, g identity
# scenario 1
n_standard(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
           type="continuous",link="identity",sigma_y_sq=0.250625,phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
      mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
      type="continuous",link="identity",
      sigma_11_sq=0.245,sigma_21_sq=0.245,sigma_10_sq=0.245,sigma_20_sq=0.245,
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
        mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
        type="continuous",link="identity",
        sigma_11_sq=0.245,sigma_21_sq=0.245,sigma_10_sq=0.245,sigma_20_sq=0.245,
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
         type="continuous",link="identity",
         sigma_y_sq=0.250625,
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 2
n_standard(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
           type="continuous",link="identity",sigma_y_sq=0.250625,phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
      mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
      type="continuous",link="identity",
      sigma_11_sq=0.145,sigma_21_sq=0.145,sigma_10_sq=0.145,sigma_20_sq=0.145,
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
        mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
        type="continuous",link="identity",
        sigma_11_sq=0.145,sigma_21_sq=0.145,sigma_10_sq=0.145,sigma_20_sq=0.145,
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
         type="continuous",link="identity",
         sigma_y_sq=0.250625,
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 3
n_standard(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
           type="continuous",link="identity",sigma_y_sq=0.250625,phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
      mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
      type="continuous",link="identity",
      sigma_11_sq=0.245,sigma_21_sq=0.145,sigma_10_sq=0.245,sigma_20_sq=0.145,
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
        mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
        type="continuous",link="identity",
        sigma_11_sq=0.245,sigma_21_sq=0.145,sigma_10_sq=0.245,sigma_20_sq=0.145,
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
         type="continuous",link="identity",
         sigma_y_sq=0.250625,
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)

# Y binary, g identity
# scenario 1
n_standard(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
           type="binary",link="identity",phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
      mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
      type="binary",link="identity",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
        mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
        type="binary",link="identity",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
         type="binary",link="identity",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 2
n_standard(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
           type="binary",link="identity",phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
      mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
      type="binary",link="identity",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
        mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
        type="binary",link="identity",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
         type="binary",link="identity",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 3
n_standard(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
           type="binary",link="identity",phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
      mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
      type="binary",link="identity",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
        mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
        type="binary",link="identity",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
         type="binary",link="identity",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)

# Y binary, g logit
# scenario 1
n_standard(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
           type="binary",link="logit",phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
      mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
      type="binary",link="logit",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
        mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
        type="binary",link="logit",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
         type="binary",link="logit",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 2
n_standard(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
           type="binary",link="logit",phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
      mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
      type="binary",link="logit",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
        mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
        type="binary",link="logit",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
         type="binary",link="logit",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 3
n_standard(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
           type="binary",link="logit",phi=0.8)
n_IPW_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
      mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
      type="binary",link="logit",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
        mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
        type="binary",link="logit",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,mu_1=0.35,mu_0=0.25,
         type="binary",link="logit",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)

# Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(20120208)

# Y continuous, g identity
# scenario 1
sim_power_cat(rep=10000,n_standard=1318,n_ipw=1324,type="continuous",link="identity",
          mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
          sigma_11_sq=0.245,sigma_12_sq=0.245,sigma_01_sq=0.245,sigma_02_sq=0.245,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)
# scenario 2
sim_power_cat(rep=10000,n_standard=1318,n_ipw=1214,type="continuous",link="identity",
          mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
          sigma_11_sq=0.145,sigma_12_sq=0.145,sigma_01_sq=0.145,sigma_02_sq=0.145,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)
# scenario 3
sim_power_cat(rep=10000,n_standard=1318,n_ipw=1228,type="continuous",link="identity",
          mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
          sigma_11_sq=0.245,sigma_12_sq=0.145,sigma_01_sq=0.245,sigma_02_sq=0.145,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)

# Y binary, g identity
# scenario 1
sim_power_cat(rep=10000,n_standard=1272,n_ipw=1282,type="binary",link="identity",
          mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
          sigma_11_sq=0.245,sigma_12_sq=0.245,sigma_01_sq=0.245,sigma_02_sq=0.245,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)
# scenario 2
sim_power_cat(rep=10000,n_standard=1272,n_ipw=1186,type="binary",link="identity",
          mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
          sigma_11_sq=0.145,sigma_12_sq=0.145,sigma_01_sq=0.145,sigma_02_sq=0.145,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)
# scenario 3
sim_power_cat(rep=10000,n_standard=1092,n_ipw=1094,type="binary",link="identity",
          mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
          sigma_11_sq=0.245,sigma_12_sq=0.145,sigma_01_sq=0.245,sigma_02_sq=0.145,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)

# Y binary, g logit
# scenario 1
sim_power_cat(rep=10000,n_standard=1290,n_ipw=1300,type="binary",link="logit",
          mu_11=0.55,mu_21=0.4,mu_10=0.45,mu_20=0.3,
          sigma_11_sq=0.245,sigma_12_sq=0.245,sigma_01_sq=0.245,sigma_02_sq=0.245,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)
# scenario 2
sim_power_cat(rep=10000,n_standard=1290,n_ipw=1204,type="binary",link="logit",
          mu_11=0.8,mu_21=0.15,mu_10=0.7,mu_20=0.05,
          sigma_11_sq=0.145,sigma_12_sq=0.145,sigma_01_sq=0.145,sigma_02_sq=0.145,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)
# scenario 3
sim_power_cat(rep=10000,n_standard=1112,n_ipw=1114,type="binary",link="logit",
          mu_11=0.55,mu_21=0.15,mu_10=0.45,mu_20=0.05,
          sigma_11_sq=0.245,sigma_12_sq=0.145,sigma_01_sq=0.245,sigma_02_sq=0.145,
          expit_beta_11=0.7,expit_beta_12=0.9,expit_beta_01=0.75,expit_beta_02=0.85)



