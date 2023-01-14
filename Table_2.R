# Table 2 (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting
# when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

### Cluster Randomized Trials (CRTs)

# Weighting for a fully observed auxiliary categorical variable

# Sample Size Calculation
# Y continuous, g identity
# scenario 1
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           type="continuous",link="identity",sigma_y_sq=0.25,phi=0.8,
           delta=0.05,m=5)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
      mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
      type="continuous",link="identity",
      sigma_11_sq=0.026,sigma_21_sq=0.294,sigma_10_sq=0.026,sigma_20_sq=0.229,
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
      delta=0.05,m=5)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
        mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
        type="continuous",link="identity",
        sigma_11_sq=0.026,sigma_21_sq=0.294,sigma_10_sq=0.026,sigma_20_sq=0.229,
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
        delta=0.05,m=5)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
         type="continuous",link="identity",
         sigma_y_sq=0.25,
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
         delta=0.05,m=5)
# scenario 2
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           type="continuous",link="identity",sigma_y_sq=0.25,phi=0.8,
           delta=0.05,m=5)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
      mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
      type="continuous",link="identity",
      sigma_11_sq=0.41955,sigma_21_sq=0.026,sigma_10_sq=0.46795,sigma_20_sq=0.026,
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
      delta=0.05,m=5)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
        mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
        type="continuous",link="identity",
        sigma_11_sq=0.41955,sigma_21_sq=0.026,sigma_10_sq=0.46795,sigma_20_sq=0.026,
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
        delta=0.05,m=5)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
         type="continuous",link="identity",
         sigma_y_sq=0.25,
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
         delta=0.05,m=5)


# Y binary, g identity
# scenario 1
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           type="binary",link="identity",phi=0.8,
           delta=0.05,m=5)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
      mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
      type="binary",link="identity",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
      delta=0.05,m=5)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
        mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
        type="binary",link="identity",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
        delta=0.05,m=5)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
         type="binary",link="identity",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
         delta=0.05,m=5)
# scenario 2
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           type="binary",link="identity",phi=0.8,
           delta=0.05,m=5)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
      mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
      type="binary",link="identity",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
      delta=0.05,m=5)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
        mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
        type="binary",link="identity",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
        delta=0.05,m=5)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
         type="binary",link="identity",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
         delta=0.05,m=5)


# Y binary, g logit
# scenario 1
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
           type="binary",link="logit",phi=0.8,
           delta=0.05,m=5)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
      mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
      type="binary",link="logit",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
      delta=0.05,m=5)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
        mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
        type="binary",link="logit",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
        delta=0.05,m=5)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.6,mu_0=0.5,
         type="binary",link="logit",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
         delta=0.05,m=5)
# scenario 2
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
           type="binary",link="logit",phi=0.8,
           delta=0.05,m=5)
n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
      mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
      type="binary",link="logit",
      pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
      expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
      delta=0.05,m=5)
n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
        mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
        type="binary",link="logit",
        pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
        delta=0.05,m=5)
n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.785,mu_0=0.685,
         type="binary",link="logit",
         pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85,
         delta=0.05,m=5)


# Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(20210316)

# Y continuous, g identity
# scenario 1
sim_power_cat(rep=10000,kappa=0.5,n_standard=1524,n_IPRW=1360,n_known=1476,n_approx=1538,
              type="continuous",link="identity",
              mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
              sigma_11_sq=0.026,sigma_21_sq=0.294,sigma_10_sq=0.026,sigma_20_sq=0.229,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
              delta=0.05,m=5)
# scenario 2
sim_power_cat(rep=10000,kappa=0.5,n_standard=1524,n_IPRW=1622,n_known=1640,n_approx=1538,
              type="continuous",link="identity",
              mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
              sigma_11_sq=0.41955,sigma_21_sq=0.026,sigma_10_sq=0.46795,sigma_20_sq=0.026,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
              delta=0.05,m=5)


# Y binary, g identity
# scenario 1
sim_power_cat(rep=10000,kappa=0.5,n_standard=1494,n_IPRW=1370,n_known=1486,n_approx=1506,
              type="binary",link="identity",
              mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
              sigma_11_sq=0.09,sigma_21_sq=0.21,sigma_10_sq=0.1275,sigma_20_sq=0.1275,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
              delta=0.05,m=5)
# scenario 2
sim_power_cat(rep=10000,kappa=0.5,n_standard=1172,n_IPRW=1200,n_known=1216,n_approx=1182,
              type="binary",link="identity",
              mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
              sigma_11_sq=0.2356,sigma_21_sq=0.0475,sigma_10_sq=0.2331,sigma_20_sq=0.1924,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
              delta=0.05,m=5)


# Y binary, g logit
# scenario 1
sim_power_cat(rep=10000,kappa=0.5,n_standard=1514,n_IPRW=1388,n_known=1506,n_approx=1528,
              type="binary",link="logit",
              mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,
              sigma_11_sq=0.09,sigma_21_sq=0.21,sigma_10_sq=0.1275,sigma_20_sq=0.1275,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
              delta=0.05,m=5)
# scenario 2
sim_power_cat(rep=10000,kappa=0.5,n_standard=1200,n_IPRW=1232,n_known=1254,n_approx=1210,
              type="binary",link="logit",
              mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,
              sigma_11_sq=0.2356,sigma_21_sq=0.0475,sigma_10_sq=0.2331,sigma_20_sq=0.1924,
              expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
              delta=0.05,m=5)


