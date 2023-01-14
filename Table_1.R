# Table 1 (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
# Explore when get largest differences in sample sizes between formulas
# and make applicable for CRTs scenarios (binary correlation constraints and var >0 in continuous simulations)
# Add scenario from R shiny example
source('useFUNCTIONS.R')

### Individually Randomized Trials (IRTs)

# Weighting for a fully observed auxiliary categorical variable

# Sample Size Calculation
# Y continuous, g identity
# scenario 1
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           type="continuous",link="identity",sigma_y_sq=0.25,phi=0.8)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
           type="continuous",link="identity",
           sigma_11_sq=0.026,sigma_21_sq=0.294,sigma_10_sq=0.026,sigma_20_sq=0.229,
           pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
           expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
            mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
            type="continuous",link="identity",
            sigma_11_sq=0.026,sigma_21_sq=0.294,sigma_10_sq=0.026,sigma_20_sq=0.229,
            pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
             type="continuous",link="identity",
             sigma_y_sq=0.25,
             pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
             expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 2
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           type="continuous",link="identity",sigma_y_sq=0.25,phi=0.8)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
           type="continuous",link="identity",
           sigma_11_sq=0.41955,sigma_21_sq=0.026,sigma_10_sq=0.46795,sigma_20_sq=0.026,
           pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
           expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
            mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
            type="continuous",link="identity",
            sigma_11_sq=0.41955,sigma_21_sq=0.026,sigma_10_sq=0.46795,sigma_20_sq=0.026,
            pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
             type="continuous",link="identity",
             sigma_y_sq=0.25,
             pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
             expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 3
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.23,mu_0=0.13,
           type="continuous",link="identity",sigma_y_sq=0.0991,phi=0.748)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.23,mu_0=0.13,
           mu_11=0.2,mu_21=0.3,mu_10=0.1,mu_20=0.2,
           type="continuous",link="identity",
           sigma_11_sq=0.01,sigma_21_sq=0.3,sigma_10_sq=0.01,sigma_20_sq=0.3,
           pi_11=0.7,pi_21=0.3,pi_10=0.7,pi_20=0.3,
           expit_beta_11=0.64,expit_beta_10=0.64,expit_beta_21=1,expit_beta_20=1)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.23,mu_0=0.13,
           mu_11=0.2,mu_21=0.3,mu_10=0.1,mu_20=0.2,
           type="continuous",link="identity",
           sigma_11_sq=0.01,sigma_21_sq=0.3,sigma_10_sq=0.01,sigma_20_sq=0.3,
           pi_11=0.7,pi_21=0.3,pi_10=0.7,pi_20=0.3,
           expit_beta_11=0.64,expit_beta_10=0.64,expit_beta_21=1,expit_beta_20=1)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.23,mu_0=0.13,
             type="continuous",link="identity",
             sigma_y_sq=0.0991,
             pi_11=0.7,pi_21=0.3,pi_10=0.7,pi_20=0.3,
             expit_beta_11=0.64,expit_beta_10=0.64,expit_beta_21=1,expit_beta_20=1)
# scenario 4
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.23,mu_0=0.13,
           type="continuous",link="identity",sigma_y_sq=0.0991,phi=0.892)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.23,mu_0=0.13,
           mu_11=0.2,mu_21=0.3,mu_10=0.1,mu_20=0.2,
           type="continuous",link="identity",
           sigma_11_sq=0.01,sigma_21_sq=0.3,sigma_10_sq=0.01,sigma_20_sq=0.3,
           pi_11=0.7,pi_21=0.3,pi_10=0.7,pi_20=0.3,
           expit_beta_11=1,expit_beta_10=1,expit_beta_21=0.64,expit_beta_20=0.64)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.23,mu_0=0.13,
            mu_11=0.2,mu_21=0.3,mu_10=0.1,mu_20=0.2,
            type="continuous",link="identity",
            sigma_11_sq=0.01,sigma_21_sq=0.3,sigma_10_sq=0.01,sigma_20_sq=0.3,
            pi_11=0.7,pi_21=0.3,pi_10=0.7,pi_20=0.3,
            expit_beta_11=1,expit_beta_10=1,expit_beta_21=0.64,expit_beta_20=0.64)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.23,mu_0=0.13,
             type="continuous",link="identity",
             sigma_y_sq=0.0991,
             pi_11=0.7,pi_21=0.3,pi_10=0.7,pi_20=0.3,
             expit_beta_11=1,expit_beta_10=1,expit_beta_21=0.64,expit_beta_20=0.64)


# Y binary, g identity
# scenario 1
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           type="binary",link="identity",phi=0.8)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
           type="binary",link="identity",
           pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
           expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
            mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
            type="binary",link="identity",
            pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
             type="binary",link="identity",
             pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
             expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 2
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           type="binary",link="identity",phi=0.8)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
           type="binary",link="identity",
           pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
           expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
            mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
            type="binary",link="identity",
            pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
             type="binary",link="identity",
             pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
             expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)


# Y binary, g logit
# scenario 1
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           type="binary",link="logit",phi=0.8)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
           type="binary",link="logit",
           pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
           expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
            mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
            type="binary",link="logit",
            pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
             type="binary",link="logit",
             pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
             expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
# scenario 2
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           type="binary",link="logit",phi=0.8)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
           type="binary",link="logit",
           pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
           expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
            mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
            type="binary",link="logit",
            pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
             type="binary",link="logit",
             pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
             expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)


# Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(20220628)

# Y continuous, g identity
# scenario 1
sim_power_cat(rep=10000,kappa=0.5,n_standard=1314,n_IPRW=1150,n_known=1266,n_approx=1328,
              type="continuous",link="identity",
              pi_11=0.5,
              mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
              sigma_11_sq=0.026,sigma_21_sq=0.294,sigma_10_sq=0.026,sigma_20_sq=0.229,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85)
# scenario 2
sim_power_cat(rep=10000,kappa=0.5,n_standard=1314,n_IPRW=1412,n_known=1430,n_approx=1328,
              type="continuous",link="identity",
              pi_11=0.5,
              mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
              sigma_11_sq=0.41955,sigma_21_sq=0.026,sigma_10_sq=0.46795,sigma_20_sq=0.026,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85)
# scenario 3
sim_power_cat(rep=10000,kappa=0.5,n_standard=558,n_IPRW=434,n_known=436,n_approx=582,
              type="continuous",link="identity",
              pi_11=0.7,
              mu_11=0.2,mu_21=0.3,mu_10=0.1,mu_20=0.2,
              sigma_11_sq=0.01,sigma_21_sq=0.3,sigma_10_sq=0.01,sigma_20_sq=0.3,
              expit_beta_11=0.64,expit_beta_21=1,expit_beta_10=0.64,expit_beta_20=1)
# scenario 4
sim_power_cat(rep=10000,kappa=0.5,n_standard=468,n_IPRW=630,n_known=634,n_approx=488,
              type="continuous",link="identity",
              pi_11=0.7,
              mu_11=0.2,mu_21=0.3,mu_10=0.1,mu_20=0.2,
              sigma_11_sq=0.01,sigma_21_sq=0.3,sigma_10_sq=0.01,sigma_20_sq=0.3,
              expit_beta_11=1,expit_beta_21=0.64,expit_beta_10=1,expit_beta_20=0.64)


# Y binary, g identity
# scenario 1
sim_power_cat(rep=10000,kappa=0.5,n_standard=1288,n_IPRW=1164,n_known=1280,n_approx=1300,
              type="binary",link="identity",
              pi_11=0.5,
              mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
              sigma_11_sq=0.09,sigma_21_sq=0.21,sigma_10_sq=0.1275,sigma_20_sq=0.1275,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85)
# scenario 2
sim_power_cat(rep=10000,kappa=0.5,n_standard=1012,n_IPRW=1038,n_known=1056,n_approx=1020,
              type="binary",link="identity",
              pi_11=0.5,
              mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
              sigma_11_sq=0.2356,sigma_21_sq=0.0475,sigma_10_sq=0.2331,sigma_20_sq=0.1924,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85)


# Y binary, g logit
# scenario 1
sim_power_cat(rep=10000,kappa=0.5,n_standard=1306,n_IPRW=1180,n_known=1298,n_approx=1318,
              type="binary",link="logit",
              pi_11=0.5,
              mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
              sigma_11_sq=0.09,sigma_21_sq=0.21,sigma_10_sq=0.1275,sigma_20_sq=0.1275,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85)
# scenario 2
sim_power_cat(rep=10000,kappa=0.5,n_standard=1034,n_IPRW=1068,n_known=1088,n_approx=1044,
              type="binary",link="logit",
              pi_11=0.5,
              mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
              sigma_11_sq=0.2356,sigma_21_sq=0.0475,sigma_10_sq=0.2331,sigma_20_sq=0.1924,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85)


