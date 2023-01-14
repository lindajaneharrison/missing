# Text (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

### Cluster Randomized Trials (CRTs)

# Weighting for a fully observed auxiliary continuous variable

# (X,Y) bivariate normal, g identity
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
           type="continuous",link="identity",sigma_y_sq=0.245,phi=0.8,
           delta=0.05,m=5)
n_IPRW_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
           mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
           beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64,
           delta=0.05,m=5)
n_known_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
        mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
        beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64,
        delta=0.05,m=5)
n_approx_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
         mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,
         beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64,
         delta=0.05,m=5)

# Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(4376758)

sim_power_cont(rep=10000,kappa=0.5,n_standard=1494,n_IPRW=1686,n_known=2042,n_approx=1634,
          mu_1=0.475,mu_0=0.375,mu_x=0,
          sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
          beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64,
          delta=0.05,m=5)


