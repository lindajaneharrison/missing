# R functions (Sample Size Calculation for Randomized Clinical Trials  via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)

# Sample Size Calculation Functions

# Standard
n_standard <- function(power,alpha,kappa,mu_1,mu_0,type,link,sigma_y_sq=NULL,phi,
                       delta=NULL,m=NULL)
{
  if (power<=0|power>=1){stop("Power must be strictly between 0 and 1")}
  if (kappa<=0|kappa>=1){stop("Kappa P(Z_i=1) must be strictly between 0 and 1")}
  if (type=="continuous" & link=="logit"){stop("Can only use logit for binary outcomes")}
  if (type=="continuous" && link=="identity") {
    eta_standard <- sigma_y_sq/(phi*kappa*(1-kappa))
    n_standard <- eta_standard*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  }
  if (type=="binary") {
    if (link=="identity") {
      eta_standard <- (mu_1*(1-mu_1)/kappa+mu_0*(1-mu_0)/(1-kappa))/phi
      n_standard <- eta_standard*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
    }
    if (link=="logit") {
      eta_standard <- (1/(kappa*mu_1*(1-mu_1))+1/((1-kappa)*mu_0*(1-mu_0)))/phi
      n_standard <- eta_standard*(qnorm(power)+qnorm(1-alpha/2))^2/(qlogis(mu_1)-qlogis(mu_0))^2
    }
  }
  if (!is.null(delta)) {
    n_standard <- (1+(m-1)*phi*delta)*n_standard  # CRT
  }
  return(n_standard)
}

# Weighting for a fully observed auxiliary categorical variable
n_IPRW_cat <- function(power,alpha,kappa,mu_1,mu_0,
                  mu_11,mu_10,mu_21,mu_20,
                  type,link,
                  sigma_11_sq=NULL,sigma_10_sq=NULL,sigma_21_sq=NULL,sigma_20_sq=NULL,
                  pi_1,pi_2,
                  expit_beta_11,expit_beta_10,expit_beta_21,expit_beta_20,
                  delta=NULL,m=NULL)
{ 
  if (power<=0|power>=1){stop("Power must be strictly between 0 and 1")}
  if (kappa<=0|kappa>=1){stop("Kappa P(Z_i=1) must be strictly between 0 and 1")}
  if (type=="continuous" & link=="logit"){stop("Can only use logit for binary outcomes")}
  if (mu_1!=pi_1*mu_11+pi_2*mu_21){stop("mu_1 should equal pi_1*mu_11+pi_2*mu_21")}
  if (mu_0!=pi_1*mu_10+pi_2*mu_20){stop("mu_0 should equal pi_1*mu_10+pi_2*mu_20")}
  if (type=="continuous" && link=="identity") {
    eta_IPRW <- (pi_1/kappa)*(sigma_11_sq/expit_beta_11+(mu_11-mu_1)^2) +
                (pi_1/(1-kappa))*(sigma_10_sq/expit_beta_10+(mu_10-mu_0)^2) +
                (pi_2/kappa)*(sigma_21_sq/expit_beta_21+(mu_21-mu_1)^2) +
                (pi_2/(1-kappa))*(sigma_20_sq/expit_beta_20+(mu_20-mu_0)^2) 
    if (!is.null(delta)) {
      sigma_y_sq <- (pi_1*(sigma_11_sq+(mu_11-mu_1)^2)+pi_2*(sigma_21_sq+(mu_21-mu_1)^2))*kappa +
                    (pi_1*(sigma_10_sq+(mu_10-mu_0)^2)+pi_2*(sigma_20_sq+(mu_20-mu_0)^2))*(1-kappa)
      eta_IPRW <- eta_IPRW + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
    }
    n_IPRW <- eta_IPRW*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  }
  if (type=="binary") {
    if (link=="identity") {
      eta_IPRW <- (pi_1/kappa)*(mu_11*(1-mu_11)/expit_beta_11+(mu_11-mu_1)^2) +
                  (pi_1/(1-kappa))*(mu_10*(1-mu_10)/expit_beta_10+(mu_10-mu_0)^2) +
                  (pi_2/kappa)*(mu_21*(1-mu_21)/expit_beta_21+(mu_21-mu_1)^2) +
                  (pi_2/(1-kappa))*(mu_20*(1-mu_20)/expit_beta_20+(mu_20-mu_0)^2) 
      if (!is.null(delta)) {
        sigma_y_sq <- mu_1*(1-mu_1)/kappa+mu_0*(1-mu_0)/(1-kappa)
        eta_IPRW <- eta_IPRW + (m-1)*delta*sigma_y_sq
      }
      n_IPRW <- eta_IPRW*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
    }
    if (link=="logit") {
      eta_IPRW <- (pi_1/kappa)*(mu_11*(1-mu_11)/(expit_beta_11*(mu_1*(1-mu_1))^2)+((mu_11-mu_1)/(mu_1*(1-mu_1)))^2) +
                  (pi_1/(1-kappa))*(mu_10*(1-mu_10)/(expit_beta_10*(mu_0*(1-mu_0))^2)+((mu_10-mu_0)/(mu_0*(1-mu_0)))^2) +
                  (pi_2/kappa)*(mu_21*(1-mu_21)/(expit_beta_21*(mu_1*(1-mu_1))^2)+((mu_21-mu_1)/(mu_1*(1-mu_1)))^2) +
                  (pi_2/(1-kappa))*(mu_20*(1-mu_20)/(expit_beta_20*(mu_0*(1-mu_0))^2)+((mu_20-mu_0)/(mu_0*(1-mu_0)))^2)
      if (!is.null(delta)) {
        sigma_y_sq <- 1/(kappa*mu_1*(1-mu_1))+1/((1-kappa)*mu_0*(1-mu_0))
        eta_IPRW <- eta_IPRW + (m-1)*delta*sigma_y_sq
      }
      n_IPRW <- eta_IPRW*(qnorm(power)+qnorm(1-alpha/2))^2/(qlogis(mu_1)-qlogis(mu_0))^2
    }
  }
  return(n_IPRW)
}

n_known_cat <- function(power,alpha,kappa,mu_1,mu_0,
                    mu_11,mu_10,mu_21,mu_20,
                    type,link,
                    sigma_11_sq=NULL,sigma_10_sq=NULL,sigma_21_sq=NULL,sigma_20_sq=NULL,
                    pi_1,pi_2,
                    expit_beta_11,expit_beta_10,expit_beta_21,expit_beta_20,
                    delta=NULL,m=NULL)
{
  if (power<=0|power>=1){stop("Power must be strictly between 0 and 1")}
  if (kappa<=0|kappa>=1){stop("Kappa P(Z_i=1) must be strictly between 0 and 1")}
  if (type=="continuous" & link=="logit"){stop("Can only use logit for binary outcomes")}
  if (mu_1!=pi_1*mu_11+pi_2*mu_21){stop("mu_1 should equal pi_1*mu_11+pi_2*mu_21")}
  if (mu_0!=pi_1*mu_10+pi_2*mu_20){stop("mu_0 should equal pi_1*mu_10+pi_2*mu_20")}
  if (type=="continuous" && link=="identity") {
    eta_known <- (pi_1/kappa)*((sigma_11_sq+(mu_11-mu_1)^2)/expit_beta_11) +
                 (pi_1/(1-kappa))*((sigma_10_sq+(mu_10-mu_0)^2)/expit_beta_10) +
                 (pi_2/kappa)*((sigma_21_sq+(mu_21-mu_1)^2)/expit_beta_21) +
                 (pi_2/(1-kappa))*((sigma_20_sq+(mu_20-mu_0)^2)/expit_beta_20) 
    if (!is.null(delta)) {
      sigma_y_sq <- (pi_1*(sigma_11_sq+(mu_11-mu_1)^2)+pi_2*(sigma_21_sq+(mu_21-mu_1)^2))*kappa +
                    (pi_1*(sigma_10_sq+(mu_10-mu_0)^2)+pi_2*(sigma_20_sq+(mu_20-mu_0)^2))*(1-kappa)
      eta_known <- eta_known + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
    }
    n_known <- eta_known*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  }
  if (type=="binary") {
    if (link=="identity") {
      eta_known <- (pi_1/kappa)*((mu_11*(1-mu_11)+(mu_11-mu_1)^2)/expit_beta_11) +
                   (pi_1/(1-kappa))*((mu_10*(1-mu_10)+(mu_10-mu_0)^2)/expit_beta_10) +
                   (pi_2/kappa)*((mu_21*(1-mu_21)+(mu_21-mu_1)^2)/expit_beta_21) +
                   (pi_2/(1-kappa))*((mu_20*(1-mu_20)+(mu_20-mu_0)^2)/expit_beta_20) 
      if (!is.null(delta)) {
        sigma_y_sq <- mu_1*(1-mu_1)/kappa+mu_0*(1-mu_0)/(1-kappa)
        eta_known <- eta_known + (m-1)*delta*sigma_y_sq
      }
      n_known <- eta_known*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
    }
    if (link=="logit") {
      eta_known <- (pi_1/kappa)*((mu_11*(1-mu_11)+(mu_11-mu_1)^2)/(expit_beta_11*(mu_1*(1-mu_1))^2)) +
                   (pi_1/(1-kappa))*((mu_10*(1-mu_10)+(mu_10-mu_0)^2)/(expit_beta_10*(mu_0*(1-mu_0))^2)) +
                   (pi_2/kappa)*((mu_21*(1-mu_21)+(mu_21-mu_1)^2)/(expit_beta_21*(mu_1*(1-mu_1))^2)) +
                   (pi_2/(1-kappa))*((mu_20*(1-mu_20)+(mu_20-mu_0)^2)/(expit_beta_20*(mu_0*(1-mu_0))^2))
      if (!is.null(delta)) {
        sigma_y_sq <- 1/(kappa*mu_1*(1-mu_1))+1/((1-kappa)*mu_0*(1-mu_0))
        eta_known <- eta_known + (m-1)*delta*sigma_y_sq
      }
      n_known <- eta_known*(qnorm(power)+qnorm(1-alpha/2))^2/(qlogis(mu_1)-qlogis(mu_0))^2
    }
  }
  return(n_known)
}

n_approx_cat <- function(power,alpha,kappa,mu_1,mu_0,
                     type,link,
                     sigma_y_sq=NULL,
                     pi_1,pi_2,
                     expit_beta_11,expit_beta_10,expit_beta_21,expit_beta_20,
                     delta=NULL,m=NULL)
{
  if (power<=0|power>=1){stop("Power must be strictly between 0 and 1")}
  if (kappa<=0|kappa>=1){stop("Kappa P(Z_i=1) must be strictly between 0 and 1")}
  if (type=="continuous" & link=="logit"){stop("Can only use logit for binary outcomes")}
  if (type=="continuous" && link=="identity") {
    eta_approx <- (pi_1/kappa)*(sigma_y_sq/expit_beta_11) +
                  (pi_1/(1-kappa))*(sigma_y_sq/expit_beta_10) +
                  (pi_2/kappa)*(sigma_y_sq/expit_beta_21) +
                  (pi_2/(1-kappa))*(sigma_y_sq/expit_beta_20) 
    if (!is.null(delta)) {
      eta_approx <- eta_approx + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
    }
    n_approx <- eta_approx*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  }
  if (type=="binary") {
    if (link=="identity") {
      eta_approx <- (pi_1/kappa)*(mu_1*(1-mu_1)/expit_beta_11) +
                    (pi_1/(1-kappa))*(mu_0*(1-mu_0)/expit_beta_10) +
                    (pi_2/kappa)*(mu_1*(1-mu_1)/expit_beta_21) +
                    (pi_2/(1-kappa))*(mu_0*(1-mu_0)/expit_beta_20) 
      if (!is.null(delta)) {
        sigma_y_sq <- mu_1*(1-mu_1)/kappa+mu_0*(1-mu_0)/(1-kappa)
        eta_approx <- eta_approx + (m-1)*delta*sigma_y_sq
      }
      n_approx <- eta_approx*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
    }
    if (link=="logit") {
      eta_approx <- (pi_1/kappa)*((1/(mu_1*(1-mu_1)))/expit_beta_11) + 
                    (pi_1/(1-kappa))*((1/(mu_0*(1-mu_0)))/expit_beta_10) +
                    (pi_2/kappa)*((1/(mu_1*(1-mu_1)))/expit_beta_21) + 
                    (pi_2/(1-kappa))*((1/(mu_0*(1-mu_0)))/expit_beta_20)
      if (!is.null(delta)) {
        sigma_y_sq <- 1/(kappa*mu_1*(1-mu_1))+1/((1-kappa)*mu_0*(1-mu_0))
        eta_approx <- eta_approx + (m-1)*delta*sigma_y_sq
      }
      n_approx <- eta_approx*(qnorm(power)+qnorm(1-alpha/2))^2/(qlogis(mu_1)-qlogis(mu_0))^2
    }
  }
  return(n_approx)
}

# Check input across functions
check_it <- function(phi,
                     pi_1,pi_2,
                     expit_beta_11,expit_beta_21,expit_beta_10,expit_beta_20,
                     sigma_y_sq,
                     sigma_11_sq,sigma_21_sq,sigma_10_sq,sigma_20_sq,
                     mu_11,mu_21,mu_10,mu_20,
                     mu_1,mu_0,
                     type)
{
  if (phi!=pi_1*expit_beta_11+pi_2*expit_beta_21){stop("phi should equal pi_1*expit_beta__11+pi_2*expit_beta__21")}
  if (phi!=pi_1*expit_beta_10+pi_2*expit_beta_20){stop("phi should equal pi_1*expit_beta__10+pi_2*expit_beta__20")}
  if (sigma_y_sq!=pi_1*(sigma_10_sq+(mu_10-mu_0)^2)+pi_2*(sigma_20_sq+(mu_20-mu_0)^2) & type=="continuous"){stop("sigma_y_sq should equal pi_1*(sigma_10_sq+(mu_10-mu_0)^2)+pi_2*(sigma_20_sq+(mu_20-mu_0)^2)")}
  if (sigma_y_sq!=pi_1*(sigma_11_sq+(mu_11-mu_1)^2)+pi_2*(sigma_21_sq+(mu_21-mu_1)^2) & type=="continuous"){stop("sigma_y_sq should equal pi_1*(sigma_11_sq+(mu_11-mu_1)^2)+pi_2*(sigma_21_sq+(mu_21-mu_1)^2)")}
}

# Weighting for a fully observed auxiliary continuous variable
n_IPRW_cont <- function(power,alpha,kappa,mu_1,mu_0,
                       mu_x,sigma_x_sq,sigma_y_sq,rho,
                       beta_11,beta_01,beta_10,beta_00,
                       delta=NULL,m=NULL)
{
  if (power<=0|power>=1){stop("Power must be strictly between 0 and 1")}
  if (kappa<=0|kappa>=1){stop("Kappa P(Z_i=1) must be strictly between 0 and 1")}
  require(fastGHQuad)
  x_j <- gaussHermiteData(100)$x 
  w_j <- gaussHermiteData(100)$w
  C_1 <- c(0,sqrt(sigma_x_sq)) -
         sqrt(2/pi)*c(sum(x_j*plogis(beta_01+beta_11*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j),
                      sum(x_j*(sqrt(2*sigma_x_sq)*x_j+mu_x)*plogis(beta_01+beta_11*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j))
  D_1 <- sqrt(1/pi)*matrix(c(sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*w_j),
                  sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                  sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                  sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)^2*w_j)),
                byrow=T,ncol=2,nrow=2)
  C_0 <- c(0,sqrt(sigma_x_sq)) -
         sqrt(2/pi)*c(sum(x_j*plogis(beta_00+beta_10*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j),
                      sum(x_j*(sqrt(2*sigma_x_sq)*x_j+mu_x)*plogis(beta_00+beta_10*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j))
  D_0 <- sqrt(1/pi)*matrix(c(sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*w_j),
                  sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                  sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                  sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)^2*w_j)),
                byrow=T,ncol=2,nrow=2)
  eta_IPRW <- sigma_y_sq*(1/kappa+exp(-mu_x*beta_11+beta_11^2*sigma_x_sq/2-beta_01)*(1+rho^2*sigma_x_sq*beta_11^2)/kappa+
                          1/(1-kappa)+exp(-mu_x*beta_10+beta_10^2*sigma_x_sq/2-beta_00)*(1+rho^2*sigma_x_sq*beta_10^2)/(1-kappa)-
                            rho^2*t(C_1)%*%solve(D_1)%*%C_1/kappa - rho^2*t(C_0)%*%solve(D_0)%*%C_0/(1-kappa)) 
  if (!is.null(delta)) {
    eta_IPRW <- eta_IPRW + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
  }
  n_IPRW <- eta_IPRW*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  return(as.numeric(n_IPRW))
}

n_known_cont <- function(power,alpha,kappa,mu_1,mu_0,
                    mu_x,sigma_x_sq,sigma_y_sq,rho,
                    beta_11,beta_01,beta_10,beta_00,
                    delta=NULL,m=NULL)
{
  if (power<=0|power>=1){stop("Power must be strictly between 0 and 1")}
  if (kappa<=0|kappa>=1){stop("Kappa P(Z_i=1) must be strictly between 0 and 1")}
  eta_known <- sigma_y_sq*(1/kappa+exp(-mu_x*beta_11+beta_11^2*sigma_x_sq/2-beta_01)*(1+rho^2*sigma_x_sq*beta_11^2)/kappa+
                           1/(1-kappa)+exp(-mu_x*beta_10+beta_10^2*sigma_x_sq/2-beta_00)*(1+rho^2*sigma_x_sq*beta_10^2)/(1-kappa))
  if (!is.null(delta)) {
    eta_known <- eta_known + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
  }
  n_known <- eta_known*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  return(n_known)
}

n_approx_cont <- function(power,alpha,kappa,mu_1,mu_0,
                     mu_x,sigma_x_sq,sigma_y_sq,
                     beta_11,beta_01,beta_10,beta_00,
                     delta=NULL,m=NULL)
{
  if (power<=0|power>=1){stop("Power must be strictly between 0 and 1")}
  if (kappa<=0|kappa>=1){stop("Kappa P(Z_i=1) must be strictly between 0 and 1")}
  eta_approx <- sigma_y_sq*(1/kappa+exp(-mu_x*beta_11+beta_11^2*sigma_x_sq/2-beta_01)/kappa+
                            1/(1-kappa)+exp(-mu_x*beta_10+beta_10^2*sigma_x_sq/2-beta_00)/(1-kappa))
  if (!is.null(delta)) {
    eta_approx <- eta_approx + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
  }
  n_approx <- eta_approx*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  return(n_approx)
}



