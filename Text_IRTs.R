# Text (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

### Individually Randomized Trials (IRTs)

# Weighting for a fully observed auxiliary continuous variable

# Pick beta_11, beta_01, beta_10 and beta_00, so phi=0.8
library(fastGHQuad)
x_j <- gaussHermiteData(100)$x 
w_j <- gaussHermiteData(100)$w
sigma_x_sq <- 1
mu_x <- 0
beta_01 <- 1.45
beta_11 <- 0.47
sum(plogis(beta_01+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_11)*w_j)/sqrt(pi)
beta_01 <- 1.39
beta_11 <- 0.11
sum(plogis(beta_01+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_11)*w_j)/sqrt(pi)
beta_00 <- 1.4
beta_10 <- 0.21
sum(plogis(beta_00+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_10)*w_j)/sqrt(pi)

# (X,Y) bivariate normal, g identity
n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
           type="continuous",link="identity",sigma_y_sq=0.245,phi=0.8)
f <- function(beta_10) sum(plogis(beta_00+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_10)*w_j)/sqrt(pi)-0.8
for (beta_00 in c(1.39,seq(1.4,3,0.1))){   
  print("scenario")
  print(beta_00)
  beta_10 <- uniroot(f, c(0,5))$root
  print(beta_10)
  print(n_IPRW_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                    mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
                    beta_01=1.4,beta_11=0.21,beta_00=beta_00,beta_10=beta_10))
  print(n_known_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                    mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
                    beta_01=1.4,beta_11=0.21,beta_00=beta_00,beta_10=beta_10))
  print(n_approx_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                    mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,
                    beta_01=1.4,beta_11=0.21,beta_00=beta_00,beta_10=beta_10))
}

# Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(5376758)

sim_power_cont(rep=10000,kappa=0.5,n_standard=1288,n_IPRW=1480,n_known=1836,n_approx=1428,
               mu_1=0.475,mu_0=0.375,mu_x=0,
               sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
               beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64)





