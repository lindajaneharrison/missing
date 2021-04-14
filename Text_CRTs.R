# Text (Sample Size Calculation for Randomized Clinical Trials when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

### Cluster Randomized Trials (CRTs)

# Weighting for a fully observed auxiliary continuous variable

# (X,Y) bivariate normal, g identity
n_standard(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
           type="continuous",link="identity",sigma_y_sq=0.245,phi=0.8,
           delta=0.1,m=10)
n_IPW_cont(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
           mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=0.5,
           beta_01=1.45,beta_11=0.47,beta_00=1.4,beta_10=0.21,
           delta=0.1,m=10)
n_known_cont(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
        mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.5,
        beta_01=1.45,beta_11=0.47,beta_00=1.4,beta_10=0.21,
        delta=0.1,m=10)
n_approx_cont(power=0.9,alpha=0.05,mu_1=0.475,mu_0=0.375,
         mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,
         beta_01=1.45,beta_11=0.47,beta_00=1.4,beta_10=0.21,
         delta=0.1,m=10)

# Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(4376758)

sim_power_cont(rep=10000,n_standard=2214,n_ipw=2224,
          mu_1=0.475,mu_0=0.375,mu_x=0,
          sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.5,
          beta_01=1.45,beta_11=0.47,beta_00=1.4,beta_10=0.21,
          delta=0.1,m=10)


