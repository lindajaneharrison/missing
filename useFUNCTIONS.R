# R functions (Sample Size Calculation for Randomized Clinical Trials  via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)

# Sample Size Calculation Functions

# Standard
n_standard <- function(power,alpha,kappa,mu_1,mu_0,type,link,sigma_y_sq=NULL,phi,
                       delta=NULL,m=NULL)
{
  if (type=="continuous" && link=="identity") {
    eta_standard <- sigma_y_sq/(phi*kappa*(1-kappa))
    n_standard <- eta_standard*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  }
  if (type=="binary") {
    if (!is.null(sigma_y_sq)) {print("sigma_y_sq should not be specified when type is binary")}
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
                  pi_11,pi_21,pi_10,pi_20,
                  expit_beta_11,expit_beta_10,expit_beta_21,expit_beta_20,
                  delta=NULL,m=NULL)
{
  if (type=="continuous" && link=="identity") {
    eta_IPRW <- (pi_11/kappa)*(sigma_11_sq/expit_beta_11+(mu_11-mu_1)^2) +
                (pi_10/(1-kappa))*(sigma_10_sq/expit_beta_10+(mu_10-mu_0)^2) +
                (pi_21/kappa)*(sigma_21_sq/expit_beta_21+(mu_21-mu_1)^2) +
                (pi_20/(1-kappa))*(sigma_20_sq/expit_beta_20+(mu_20-mu_0)^2) 
    if (!is.null(delta)) {
      sigma_y_sq <- (pi_11*(sigma_11_sq+(mu_11-mu_1)^2)+pi_21*(sigma_21_sq+(mu_21-mu_1)^2))*kappa +
                    (pi_10*(sigma_10_sq+(mu_10-mu_0)^2)+pi_20*(sigma_20_sq+(mu_20-mu_0)^2))*(1-kappa)
      eta_IPRW <- eta_IPRW + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
    }
    n_IPRW <- eta_IPRW*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  }
  if (type=="binary") {
    if (!is.null(sigma_11_sq)) {print("sigma_sq variables should not be specified when type is binary")}
    if (link=="identity") {
      eta_IPRW <- (pi_11/kappa)*(mu_11*(1-mu_11)/expit_beta_11+(mu_11-mu_1)^2) +
                  (pi_10/(1-kappa))*(mu_10*(1-mu_10)/expit_beta_10+(mu_10-mu_0)^2) +
                  (pi_21/kappa)*(mu_21*(1-mu_21)/expit_beta_21+(mu_21-mu_1)^2) +
                  (pi_20/(1-kappa))*(mu_20*(1-mu_20)/expit_beta_20+(mu_20-mu_0)^2) 
      if (!is.null(delta)) {
        sigma_y_sq <- mu_1*(1-mu_1)/kappa+mu_0*(1-mu_0)/(1-kappa)
        eta_IPRW <- eta_IPRW + (m-1)*delta*sigma_y_sq
      }
      n_IPRW <- eta_IPRW*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
    }
    if (link=="logit") {
      eta_IPRW <- (pi_11/kappa)*(mu_11*(1-mu_11)/(expit_beta_11*(mu_1*(1-mu_1))^2)+((mu_11-mu_1)/(mu_1*(1-mu_1)))^2) +
                  (pi_10/(1-kappa))*(mu_10*(1-mu_10)/(expit_beta_10*(mu_0*(1-mu_0))^2)+((mu_10-mu_0)/(mu_0*(1-mu_0)))^2) +
                  (pi_21/kappa)*(mu_21*(1-mu_21)/(expit_beta_21*(mu_1*(1-mu_1))^2)+((mu_21-mu_1)/(mu_1*(1-mu_1)))^2) +
                  (pi_20/(1-kappa))*(mu_20*(1-mu_20)/(expit_beta_20*(mu_0*(1-mu_0))^2)+((mu_20-mu_0)/(mu_0*(1-mu_0)))^2)
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
                    pi_11,pi_21,pi_10,pi_20,
                    expit_beta_11,expit_beta_10,expit_beta_21,expit_beta_20,
                    delta=NULL,m=NULL)
{
  if (type=="continuous" && link=="identity") {
    eta_known <- (pi_11/kappa)*((sigma_11_sq+(mu_11-mu_1)^2)/expit_beta_11) +
                 (pi_10/(1-kappa))*((sigma_10_sq+(mu_10-mu_0)^2)/expit_beta_10) +
                 (pi_21/kappa)*((sigma_21_sq+(mu_21-mu_1)^2)/expit_beta_21) +
                 (pi_20/(1-kappa))*((sigma_20_sq+(mu_20-mu_0)^2)/expit_beta_20) 
    if (!is.null(delta)) {
      sigma_y_sq <- (pi_11*(sigma_11_sq+(mu_11-mu_1)^2)+pi_21*(sigma_21_sq+(mu_21-mu_1)^2))*kappa +
                    (pi_10*(sigma_10_sq+(mu_10-mu_0)^2)+pi_20*(sigma_20_sq+(mu_20-mu_0)^2))*(1-kappa)
      eta_known <- eta_known + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
    }
    n_known <- eta_known*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  }
  if (type=="binary") {
    if (!is.null(sigma_11_sq)) {print("sigma_sq variables should not be specified when type is binary")}
    if (link=="identity") {
      eta_known <- (pi_11/kappa)*((mu_11*(1-mu_11)+(mu_11-mu_1)^2)/expit_beta_11) +
                   (pi_10/(1-kappa))*((mu_10*(1-mu_10)+(mu_10-mu_0)^2)/expit_beta_10) +
                   (pi_21/kappa)*((mu_21*(1-mu_21)+(mu_21-mu_1)^2)/expit_beta_21) +
                   (pi_20/(1-kappa))*((mu_20*(1-mu_20)+(mu_20-mu_0)^2)/expit_beta_20) 
      if (!is.null(delta)) {
        sigma_y_sq <- mu_1*(1-mu_1)/kappa+mu_0*(1-mu_0)/(1-kappa)
        eta_known <- eta_known + (m-1)*delta*sigma_y_sq
      }
      n_known <- eta_known*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
    }
    if (link=="logit") {
      eta_known <- (pi_11/kappa)*((mu_11*(1-mu_11)+(mu_11-mu_1)^2)/(expit_beta_11*(mu_1*(1-mu_1))^2)) +
                   (pi_10/(1-kappa))*((mu_10*(1-mu_10)+(mu_10-mu_0)^2)/(expit_beta_10*(mu_0*(1-mu_0))^2)) +
                   (pi_21/kappa)*((mu_21*(1-mu_21)+(mu_21-mu_1)^2)/(expit_beta_21*(mu_1*(1-mu_1))^2)) +
                   (pi_20/(1-kappa))*((mu_20*(1-mu_20)+(mu_20-mu_0)^2)/(expit_beta_20*(mu_0*(1-mu_0))^2))
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
                     pi_11,pi_21,pi_10,pi_20,
                     expit_beta_11,expit_beta_10,expit_beta_21,expit_beta_20,
                     delta=NULL,m=NULL)
{
  if (type=="continuous" && link=="identity") {
    eta_approx <- (pi_11/kappa)*(sigma_y_sq/expit_beta_11) +
                  (pi_10/(1-kappa))*(sigma_y_sq/expit_beta_10) +
                  (pi_21/kappa)*(sigma_y_sq/expit_beta_21) +
                  (pi_20/(1-kappa))*(sigma_y_sq/expit_beta_20) 
    if (!is.null(delta)) {
      eta_approx <- eta_approx + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
    }
    n_approx <- eta_approx*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  }
  if (type=="binary") {
    if (!is.null(sigma_y_sq)) {print("sigma_y_sq should not be specified when type is binary")}
    if (link=="identity") {
      eta_approx <- (pi_11/kappa)*(mu_1*(1-mu_1)/expit_beta_11) +
                    (pi_10/(1-kappa))*(mu_0*(1-mu_0)/expit_beta_10) +
                    (pi_21/kappa)*(mu_1*(1-mu_1)/expit_beta_21) +
                    (pi_20/(1-kappa))*(mu_0*(1-mu_0)/expit_beta_20) 
      if (!is.null(delta)) {
        sigma_y_sq <- mu_1*(1-mu_1)/kappa+mu_0*(1-mu_0)/(1-kappa)
        eta_approx <- eta_approx + (m-1)*delta*sigma_y_sq
      }
      n_approx <- eta_approx*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
    }
    if (link=="logit") {
      eta_approx <- (pi_11/kappa)*((1/(mu_1*(1-mu_1)))/expit_beta_11) + 
                    (pi_10/(1-kappa))*((1/(mu_0*(1-mu_0)))/expit_beta_10) +
                    (pi_21/kappa)*((1/(mu_1*(1-mu_1)))/expit_beta_21) + 
                    (pi_20/(1-kappa))*((1/(mu_0*(1-mu_0)))/expit_beta_20)
      if (!is.null(delta)) {
        sigma_y_sq <- 1/(kappa*mu_1*(1-mu_1))+1/((1-kappa)*mu_0*(1-mu_0))
        eta_approx <- eta_approx + (m-1)*delta*sigma_y_sq
      }
      n_approx <- eta_approx*(qnorm(power)+qnorm(1-alpha/2))^2/(qlogis(mu_1)-qlogis(mu_0))^2
    }
  }
  return(n_approx)
}

# Weighting for a fully observed auxiliary continuous variable
n_IPRW_cont <- function(power,alpha,kappa,mu_1,mu_0,
                       mu_x,sigma_x_sq,sigma_y_sq,rho,
                       beta_11,beta_01,beta_10,beta_00,
                       delta=NULL,m=NULL)
{
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
  eta_approx <- sigma_y_sq*(1/kappa+exp(-mu_x*beta_11+beta_11^2*sigma_x_sq/2-beta_01)/kappa+
                            1/(1-kappa)+exp(-mu_x*beta_10+beta_10^2*sigma_x_sq/2-beta_00)/(1-kappa))
  if (!is.null(delta)) {
    eta_approx <- eta_approx + (m-1)*delta*sigma_y_sq/(kappa*(1-kappa))
  }
  n_approx <- eta_approx*(qnorm(power)+qnorm(1-alpha/2))^2/(mu_1-mu_0)^2
  return(n_approx)
}


# Simulation functions

# Weighting for a fully observed auxiliary categorical variable
sim_power_cat <- function(rep,n_standard,n_IPRW,n_known,n_approx,
                          type,link,
                          pi_11,
                          mu_11,mu_21,mu_10,mu_20,
                          sigma_11_sq,sigma_21_sq,sigma_10_sq,sigma_20_sq,
                          expit_beta_11,expit_beta_21,expit_beta_10,expit_beta_20,
                          delta=NULL,m=NULL,kappa)
{
  para <- foreach(i=1:rep, .combine=rbind) %dorng% {
    # n_standard
    Z <- c(rep(0,ceiling(n_standard*(1-kappa))),rep(1,ceiling(n_standard*kappa)))
    X <- rbinom(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa)),1,1-pi_11)+1
    R <- NULL
    R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
    R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
    R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
    R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
    Y <- NULL
    if (type=="continuous" && is.null(delta)) {
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
    }
    if (type=="binary" && is.null(delta)) {
      Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
      Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
      Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
      Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
    }
    if (type=="continuous" && !is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_standard*(1-kappa))/m)),each=m)[1:ceiling(n_standard*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_standard*kappa)/m)),each=m)[1:ceiling(n_standard*kappa)] + ceiling(ceiling(n_standard*(1-kappa))/m)
      mu_1 <- (mu_11+mu_21)/2
      mu_0 <- (mu_10+mu_20)/2
      sigma_y_sq <- (0.5*(sigma_11_sq+(mu_11-mu_1)^2)+0.5*(sigma_21_sq+(mu_21-mu_1)^2) +
                     0.5*(sigma_10_sq+(mu_10-mu_0)^2)+0.5*(sigma_20_sq+(mu_20-mu_0)^2))/2 
      cluster_effect <- rnorm(ceiling(ceiling(n_standard*kappa)/m)+ceiling(ceiling(n_standard*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
      zeta <- NULL
      zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_standard*(1-kappa))/m)],each=m)[1:ceiling(n_standard*(1-kappa))]
      zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_standard*kappa)/m)+1):(ceiling(ceiling(n_standard*kappa)/m)+ceiling(ceiling(n_standard*(1-kappa))/m))],each=m)[1:ceiling(n_standard*kappa)]
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==1] 
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==2] 
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==1] 
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==2] 
    }
    # method of Qaqish, Biometrika 2003
    if (type=="binary" && !is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_standard*(1-kappa))/m)),each=m)[1:ceiling(n_standard*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_standard*kappa)/m)),each=m)[1:ceiling(n_standard*kappa)] + ceiling(ceiling(n_standard*(1-kappa))/m)
      mu_vec <- NULL
      mu_vec[Z==1 & X==1] <- mu_11
      mu_vec[Z==1 & X==2] <- mu_21
      mu_vec[Z==0 & X==1] <- mu_10
      mu_vec[Z==0 & X==2] <- mu_20 
      for (k in 1:(ceiling(ceiling(n_standard*kappa)/m)+ceiling(ceiling(n_standard*(1-kappa))/m))) {
        mu_k <- mu_vec[K==k]
        Y_k <- NULL
        Y_k[1] <- rbinom(1,1,mu_k[1]) 
        if (sum(K==k)>1) {
          for (j in 2:sum(K==k)) {
            j1 <- j-1
            res <- Y_k[1:j1] - mu_k[1:j1] # residuals
            v_part <- sqrt(mu_k[j]*(1-mu_k[j]))/sqrt(mu_k[1:j1]*(1-mu_k[1:j1])) # variance part (v_ii/v_jj)^(1/2)
            lambda <- mu_k[j] + sum(res*(delta/(1+(j-2)*delta))*v_part) # cond.mean
            lambda <- min(lambda,1)
            lambda <- max(lambda,0)
            Y_k[j] <- rbinom(1,1,lambda)
          }
        }
        Y <- c(Y,Y_k)
      }
    }
    # Equation 3
    mu_hat_1_eqn3 <- sum(R*Z*Y)/sum(R*Z)
    mu_hat_0_eqn3 <- sum(R*(1-Z)*Y)/sum(R*(1-Z))
    if (link=="identity" && is.null(delta)) {
      est_eqn3 <- mu_hat_1_eqn3 - mu_hat_0_eqn3 
      var_eqn3 <- sum(R*Z*(Y-mu_hat_1_eqn3)^2)/(sum(R*Z)^2) + 
                  sum(R*(1-Z)*(Y-mu_hat_0_eqn3)^2)/(sum(R*(1-Z))^2)
    }
    if (link=="logit" && is.null(delta)) {
      est_eqn3 <- qlogis(mu_hat_1_eqn3) - qlogis(mu_hat_0_eqn3)
      var_eqn3 <- sum(R*Z*(Y-mu_hat_1_eqn3)^2)/((sum(R*Z)*mu_hat_1_eqn3*(1-mu_hat_1_eqn3))^2) + 
                  sum(R*(1-Z)*(Y-mu_hat_0_eqn3)^2)/((sum(R*(1-Z))*mu_hat_0_eqn3*(1-mu_hat_0_eqn3))^2)
    }
    if (link=="identity" && !is.null(delta)) {
      est_eqn3 <- mu_hat_1_eqn3 - mu_hat_0_eqn3 
      out1 <- NULL
      out0 <- NULL
      for (k in 1:max(K)) {
        out1[k] <- sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn3))^2
        out0[k] <- sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn3))^2
      }
      var_eqn3 <- sum(out1)/(sum(R*Z)^2) + 
                  sum(out0)/(sum(R*(1-Z))^2)
    }
    if (link=="logit" && !is.null(delta)) {
      est_eqn3 <- qlogis(mu_hat_1_eqn3) - qlogis(mu_hat_0_eqn3)
      out1 <- NULL
      out0 <- NULL
      for (k in 1:max(K)) {
        out1[k] <- sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn3))^2
        out0[k] <- sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn3))^2
      }
      var_eqn3 <- sum(out1)/((sum(R*Z)*mu_hat_1_eqn3*(1-mu_hat_1_eqn3))^2) + 
                  sum(out0)/((sum(R*(1-Z))*mu_hat_0_eqn3*(1-mu_hat_0_eqn3))^2)
    }
    wald_eqn3 <- est_eqn3/sqrt(var_eqn3)
    reject_eqn3 <- wald_eqn3>qnorm(0.975)
    # Equation 4
    X1 <- ifelse(X==1,1,0)
    X2 <- ifelse(X==2,1,0)
    model_1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],family="binomial")
    model_0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],family="binomial")
    e <- NULL
    e[Z==1 & X==1] <- plogis(model_1$coefficients[1])
    e[Z==1 & X==2] <- plogis(model_1$coefficients[2])
    e[Z==0 & X==1] <- plogis(model_0$coefficients[1])
    e[Z==0 & X==2] <- plogis(model_0$coefficients[2])
    mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
    mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
    if (link=="identity" && is.null(delta)) {
      est_standard_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
      A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa)))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa)))
    }
    if (link=="logit" && is.null(delta)) {
      est_standard_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
      A_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa)))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_standard*kappa)+ceiling(n_standard*(1-kappa)))
    }
    if (link=="identity" && !is.null(delta)) {
      est_standard_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
      A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    if (link=="logit" && !is.null(delta)) {
      est_standard_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
      A_C_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    wald_eqn4 <- est_standard_eqn4/sqrt(var_eqn4)
    reject_standard_eqn4 <- wald_eqn4>qnorm(0.975)
    # n_IPRW
    Z <- c(rep(0,ceiling(n_IPRW*(1-kappa))),rep(1,ceiling(n_IPRW*kappa)))
    X <- rbinom(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa),1,1-pi_11)+1
    R <- NULL
    R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
    R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
    R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
    R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
    Y <- NULL
    if (type=="continuous" && is.null(delta)) {
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
    }
    if (type=="binary" && is.null(delta)) {
      Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
      Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
      Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
      Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
    }
    if (type=="continuous" && !is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_IPRW*(1-kappa))/m)),each=m)[1:ceiling(n_IPRW*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_IPRW*kappa)/m)),each=m)[1:ceiling(n_IPRW*kappa)] + ceiling(ceiling(n_IPRW*(1-kappa))/m)
      mu_1 <- (mu_11+mu_21)/2
      mu_0 <- (mu_10+mu_20)/2
      sigma_y_sq <- (0.5*(sigma_11_sq+(mu_11-mu_1)^2)+0.5*(sigma_21_sq+(mu_21-mu_1)^2) +
                       0.5*(sigma_10_sq+(mu_10-mu_0)^2)+0.5*(sigma_20_sq+(mu_20-mu_0)^2))/2 
      cluster_effect <- rnorm(ceiling(ceiling(n_IPRW*kappa)/m)+ceiling(ceiling(n_IPRW*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
      zeta <- NULL
      zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_IPRW*(1-kappa))/m)],each=m)[1:ceiling(n_IPRW*(1-kappa))]
      zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_IPRW*(1-kappa))/m)+1):(ceiling(ceiling(n_IPRW*kappa)/m)+ceiling(ceiling(n_IPRW*(1-kappa))/m))],each=m)[1:ceiling(n_IPRW*kappa)]
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==1] 
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==2] 
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==1] 
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==2] 
    }
    # method of Qaqish, Biometrika 2003
    if (type=="binary" && !is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_IPRW*(1-kappa))/m)),each=m)[1:ceiling(n_IPRW*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_IPRW*kappa)/m)),each=m)[1:ceiling(n_IPRW*kappa)] + ceiling(ceiling(n_IPRW*(1-kappa))/m)
      mu_vec <- NULL
      mu_vec[Z==1 & X==1] <- mu_11
      mu_vec[Z==1 & X==2] <- mu_21
      mu_vec[Z==0 & X==1] <- mu_10
      mu_vec[Z==0 & X==2] <- mu_20 
      for (k in 1:(ceiling(ceiling(n_IPRW*kappa)/m)+ceiling(ceiling(n_IPRW*(1-kappa))/m))) {
        mu_k <- mu_vec[K==k]
        Y_k <- NULL
        Y_k[1] <- rbinom(1,1,mu_k[1]) 
        if (sum(K==k)>1) {
          for (j in 2:sum(K==k)) {
            j1 <- j-1
            res <- Y_k[1:j1] - mu_k[1:j1] # residuals
            v_part <- sqrt(mu_k[j]*(1-mu_k[j]))/sqrt(mu_k[1:j1]*(1-mu_k[1:j1])) # variance part (v_ii/v_jj)^(1/2)
            lambda <- mu_k[j] + sum(res*(delta/(1+(j-2)*delta))*v_part) # cond.mean
            lambda <- min(lambda,1)
            lambda <- max(lambda,0)
            Y_k[j] <- rbinom(1,1,lambda)
          }
        }
        Y <- c(Y,Y_k)
      }
    }
    # Equation 4
    X1 <- ifelse(X==1,1,0)
    X2 <- ifelse(X==2,1,0)
    model_1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],family="binomial")
    model_0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],family="binomial")
    e <- NULL
    e[Z==1 & X==1] <- plogis(model_1$coefficients[1])
    e[Z==1 & X==2] <- plogis(model_1$coefficients[2])
    e[Z==0 & X==1] <- plogis(model_0$coefficients[1])
    e[Z==0 & X==2] <- plogis(model_0$coefficients[2])
    mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
    mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
    if (link=="identity" && is.null(delta)) {
      est_IPRW_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
      A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
    }
    if (link=="logit" && is.null(delta)) {
      est_IPRW_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
      A_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
    }
    if (link=="identity" && !is.null(delta)) {
      est_IPRW_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
      A_C_hat <- matrix(
                 c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
                   0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
                   0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
                   0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
                   0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
                   0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
                 nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    if (link=="logit" && !is.null(delta)) {
      est_IPRW_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
      A_C_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    wald_eqn4 <- est_IPRW_eqn4/sqrt(var_eqn4)
    reject_IPRW_eqn4 <- wald_eqn4>qnorm(0.975)
    # n_known
    Z <- c(rep(0,ceiling(n_known*(1-kappa))),rep(1,ceiling(n_known*kappa)))
    X <- rbinom(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)),1,1-pi_11)+1
    R <- NULL
    R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
    R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
    R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
    R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
    Y <- NULL
    if (type=="continuous" && is.null(delta)) {
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
    }
    if (type=="binary" && is.null(delta)) {
      Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
      Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
      Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
      Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
    }
    if (type=="continuous" && !is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_known*(1-kappa))/m)),each=m)[1:ceiling(n_known*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_known*kappa)/m)),each=m)[1:ceiling(n_known*kappa)] + ceiling(ceiling(n_known*(1-kappa))/m)
      mu_1 <- (mu_11+mu_21)/2
      mu_0 <- (mu_10+mu_20)/2
      sigma_y_sq <- (0.5*(sigma_11_sq+(mu_11-mu_1)^2)+0.5*(sigma_21_sq+(mu_21-mu_1)^2) +
                       0.5*(sigma_10_sq+(mu_10-mu_0)^2)+0.5*(sigma_20_sq+(mu_20-mu_0)^2))/2 
      cluster_effect <- rnorm(ceiling(ceiling(n_known*kappa)/m)+ceiling(ceiling(n_known*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
      zeta <- NULL
      zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_known*(1-kappa))/m)],each=m)[1:ceiling(n_known*(1-kappa))]
      zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_known*(1-kappa))/m)+1):(ceiling(ceiling(n_known*kappa)/m)+ceiling(ceiling(n_known*(1-kappa))/m))],each=m)[1:ceiling(n_known*kappa)]
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==1] 
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==2] 
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==1] 
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==2] 
    }
    # method of Qaqish, Biometrika 2003
    if (type=="binary" && !is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_known*(1-kappa))/m)),each=m)[1:ceiling(n_known*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_known*kappa)/m)),each=m)[1:ceiling(n_known*kappa)] + ceiling(ceiling(n_known*(1-kappa))/m)
      mu_vec <- NULL
      mu_vec[Z==1 & X==1] <- mu_11
      mu_vec[Z==1 & X==2] <- mu_21
      mu_vec[Z==0 & X==1] <- mu_10
      mu_vec[Z==0 & X==2] <- mu_20 
      for (k in 1:(ceiling(ceiling(n_known*kappa)/m)+ceiling(ceiling(n_known*(1-kappa))/m))) {
        mu_k <- mu_vec[K==k]
        Y_k <- NULL
        Y_k[1] <- rbinom(1,1,mu_k[1]) 
        if (sum(K==k)>1) {
          for (j in 2:sum(K==k)) {
            j1 <- j-1
            res <- Y_k[1:j1] - mu_k[1:j1] # residuals
            v_part <- sqrt(mu_k[j]*(1-mu_k[j]))/sqrt(mu_k[1:j1]*(1-mu_k[1:j1])) # variance part (v_ii/v_jj)^(1/2)
            lambda <- mu_k[j] + sum(res*(delta/(1+(j-2)*delta))*v_part) # cond.mean
            lambda <- min(lambda,1)
            lambda <- max(lambda,0)
            Y_k[j] <- rbinom(1,1,lambda)
          }
        }
        Y <- c(Y,Y_k)
      }
    }
    # Equation 4
    X1 <- ifelse(X==1,1,0)
    X2 <- ifelse(X==2,1,0)
    model_1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],family="binomial")
    model_0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],family="binomial")
    e <- NULL
    e[Z==1 & X==1] <- plogis(model_1$coefficients[1])
    e[Z==1 & X==2] <- plogis(model_1$coefficients[2])
    e[Z==0 & X==1] <- plogis(model_0$coefficients[1])
    e[Z==0 & X==2] <- plogis(model_0$coefficients[2])
    mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
    mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
    if (link=="identity" && is.null(delta)) {
      est_known_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
      A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))
    }
    if (link=="logit" && is.null(delta)) {
      est_known_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
      A_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))
    }
    if (link=="identity" && !is.null(delta)) {
      est_known_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
      A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    if (link=="logit" && !is.null(delta)) {
      est_known_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
      A_C_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    wald_eqn4 <- est_known_eqn4/sqrt(var_eqn4)
    reject_known_eqn4 <- wald_eqn4>qnorm(0.975)
    # n_approx
    Z <- c(rep(0,ceiling(n_approx*(1-kappa))),rep(1,ceiling(n_approx*kappa)))
    X <- rbinom(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)),1,1-pi_11)+1
    R <- NULL
    R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
    R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
    R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
    R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
    Y <- NULL
    if (type=="continuous" && is.null(delta)) {
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
    }
    if (type=="binary" && is.null(delta)) {
      Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
      Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
      Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
      Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
    }
    if (type=="continuous" && !is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_approx*(1-kappa))/m)),each=m)[1:ceiling(n_approx*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_approx*kappa)/m)),each=m)[1:ceiling(n_approx*kappa)] + ceiling(ceiling(n_approx*(1-kappa))/m)
      mu_1 <- (mu_11+mu_21)/2
      mu_0 <- (mu_10+mu_20)/2
      sigma_y_sq <- (0.5*(sigma_11_sq+(mu_11-mu_1)^2)+0.5*(sigma_21_sq+(mu_21-mu_1)^2) +
                       0.5*(sigma_10_sq+(mu_10-mu_0)^2)+0.5*(sigma_20_sq+(mu_20-mu_0)^2))/2 
      cluster_effect <- rnorm(ceiling(ceiling(n_approx*kappa)/m)+ceiling(ceiling(n_approx*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
      zeta <- NULL
      zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_approx*(1-kappa))/m)],each=m)[1:ceiling(n_approx*(1-kappa))]
      zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_approx*(1-kappa))/m)+1):(ceiling(ceiling(n_approx*kappa)/m)+ceiling(ceiling(n_approx*(1-kappa))/m))],each=m)[1:ceiling(n_approx*kappa)]
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==1] 
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==2] 
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==1] 
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==2] 
    }
    # method of Qaqish, Biometrika 2003
    if (type=="binary" && !is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_approx*(1-kappa))/m)),each=m)[1:ceiling(n_approx*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_approx*kappa)/m)),each=m)[1:ceiling(n_approx*kappa)] + ceiling(ceiling(n_approx*(1-kappa))/m)
      mu_vec <- NULL
      mu_vec[Z==1 & X==1] <- mu_11
      mu_vec[Z==1 & X==2] <- mu_21
      mu_vec[Z==0 & X==1] <- mu_10
      mu_vec[Z==0 & X==2] <- mu_20 
      for (k in 1:(ceiling(ceiling(n_approx*kappa)/m)+ceiling(ceiling(n_approx*(1-kappa))/m))) {
        mu_k <- mu_vec[K==k]
        Y_k <- NULL
        Y_k[1] <- rbinom(1,1,mu_k[1]) 
        if (sum(K==k)>1) {
          for (j in 2:sum(K==k)) {
            j1 <- j-1
            res <- Y_k[1:j1] - mu_k[1:j1] # residuals
            v_part <- sqrt(mu_k[j]*(1-mu_k[j]))/sqrt(mu_k[1:j1]*(1-mu_k[1:j1])) # variance part (v_ii/v_jj)^(1/2)
            lambda <- mu_k[j] + sum(res*(delta/(1+(j-2)*delta))*v_part) # cond.mean
            lambda <- min(lambda,1)
            lambda <- max(lambda,0)
            Y_k[j] <- rbinom(1,1,lambda)
          }
        }
        Y <- c(Y,Y_k)
      }
    }
    # Equation 4
    X1 <- ifelse(X==1,1,0)
    X2 <- ifelse(X==2,1,0)
    model_1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],family="binomial")
    model_0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],family="binomial")
    e <- NULL
    e[Z==1 & X==1] <- plogis(model_1$coefficients[1])
    e[Z==1 & X==2] <- plogis(model_1$coefficients[2])
    e[Z==0 & X==1] <- plogis(model_0$coefficients[1])
    e[Z==0 & X==2] <- plogis(model_0$coefficients[2])
    mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
    mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
    if (link=="identity" && is.null(delta)) {
      est_approx_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
      A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
    }
    if (link=="logit" && is.null(delta)) {
      est_approx_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
      A_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
    }
    if (link=="identity" && !is.null(delta)) {
      est_approx_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
      A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    if (link=="logit" && !is.null(delta)) {
      est_approx_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
      A_C_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    wald_eqn4 <- est_approx_eqn4/sqrt(var_eqn4)
    reject_approx_eqn4 <- wald_eqn4>qnorm(0.975)
    return(c(est_eqn3,est_standard_eqn4,est_IPRW_eqn4,est_known_eqn4,est_approx_eqn4,
             reject_eqn3,reject_standard_eqn4,reject_IPRW_eqn4,reject_known_eqn4,reject_approx_eqn4))
  }
  results <- colMeans(para)
  names(results) <- c("g(mu_hat_1)-g(mu_hat_0) (eqn 3)",
                      "g(mu_hat_1)-g(mu_hat_0) (standard eqn 4)",
                      "g(mu_hat_1)-g(mu_hat_0) (IPRW eqn 4)",
                      "g(mu_hat_1)-g(mu_hat_0) (known eqn 4)",
                      "g(mu_hat_1)-g(mu_hat_0) (approx eqn 4)",
                      "empirical power (eqn 3)",
                      "empirical power (standard eqn 4)",
                      "empirical power (IPRW eqn 4)",
                      "empirical power (known eqn 4)",
                      "empirical power (approx eqn 4)")
  return(results)
}

# Weighting for a fully observed auxiliary categorical variable - pilot data
sim_power_cat_pilot <- function(rep,type,link,
                          pi_11,
                          mu_11,mu_21,mu_10,mu_20,
                          sigma_11_sq,sigma_21_sq,sigma_10_sq,sigma_20_sq,
                          expit_beta_11,expit_beta_21,expit_beta_10,expit_beta_20,
                          delta=NULL,m=NULL,kappa,n_pilot)
{
  para <- foreach(i=1:rep, .combine=rbind) %dorng% {
    # Simulate pilot data
    Z <- c(rep(0,ceiling(n_pilot*(1-kappa))),rep(1,ceiling(n_pilot*kappa)))
    X <- rbinom(ceiling(n_pilot*kappa)+ceiling(n_pilot*(1-kappa)),1,1-pi_11)+1
    R <- NULL
    R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
    R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
    R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
    R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
    Y <- NULL
    if (type=="continuous") {
      Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
      Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
      Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
      Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
    }
    if (type=="binary") {
      Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
      Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
      Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
      Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
    }
    # Estimate parameters from pilot study
    X <- X-1
    pi_11_out <- sum(X & Z)/sum(Z)
    pi_21_out <- sum((1-X) & Z)/sum(Z)
    pi_10_out <- sum(X & (1-Z))/sum(1-Z)
    pi_20_out <- sum((1-X) & (1-Z))/sum(1-Z)
    mu_11_out <- mean(Y[X & Z & R],na.rm=TRUE)
    mu_21_out <- mean(Y[(1-X) & Z & R],na.rm=TRUE)
    mu_10_out <- mean(Y[X & (1-Z) & R],na.rm=TRUE)
    mu_20_out <- mean(Y[(1-X) & (1-Z) & R],na.rm=TRUE)
    sigma_11_sq_out <- var(Y[X & Z & R],na.rm=TRUE)*(sum(X & Z & R)-1)/sum(X & Z & R)
    sigma_21_sq_out <- var(Y[(1-X) & Z & R],na.rm=TRUE)*(sum((1-X) & Z & R)-1)/sum((1-X) & Z & R)
    sigma_10_sq_out <- var(Y[X & (1-Z) & R],na.rm=TRUE)*(sum(X & (1-Z) & R)-1)/sum(X & (1-Z) & R)
    sigma_20_sq_out <- var(Y[(1-X) & (1-Z) & R],na.rm=TRUE)*(sum((1-X) & (1-Z) & R)-1)/sum((1-X) & (1-Z) & R)
    expit_beta_11_out <- mean(R[X & Z])
    expit_beta_21_out <- mean(R[(1-X) & Z])
    expit_beta_10_out <- mean(R[X & (1-Z)])
    expit_beta_20_out <- mean(R[(1-X) & (1-Z)])
    X1 <- ifelse(X==1,1,0)
    X2 <- ifelse(X==0,1,0)
    miss_model1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],data=as.data.frame(cbind(R,Z,X1,X2)),family=binomial)
    miss_model0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],data=as.data.frame(cbind(R,Z,X1,X2)),family=binomial)
    e <- NULL
    e[X==1 & Z==1] <- plogis(miss_model1$coefficients[1])
    e[X==0 & Z==1] <- plogis(miss_model1$coefficients[2])
    e[Z==0 & Z==0] <- plogis(miss_model0$coefficients[1])
    e[X==0 & Z==0] <- plogis(miss_model0$coefficients[2])
    mu_1_hat <- sum(Z*R*Y/e,na.rm=TRUE)/sum(Z*R/e)
    mu_0_hat <- sum((1-Z)*R*Y/e,na.rm=TRUE)/sum((1-Z)*R/(e))  
    sigma_y_sq_out <- (pi_11_out*(sigma_11_sq_out+(mu_11_out-mu_1_hat)^2)+pi_21_out*(sigma_21_sq_out+(mu_21_out-mu_1_hat)^2))*kappa +
                      (pi_10_out*(sigma_10_sq_out+(mu_10_out-mu_0_hat)^2)+pi_20_out*(sigma_20_sq_out+(mu_20_out-mu_0_hat)^2))*(1-kappa)
    phi_out <- mean(R)
    # Calculate sample size from pilot study estimates
    n_standard_out <- n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1_hat,mu_0=mu_0_hat,
                      type=type,link=link,sigma_y_sq=sigma_y_sq_out,phi=phi_out)
    n_IPRW <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1_hat,mu_0=mu_0_hat,
                      mu_11=mu_11_out,mu_21=mu_21_out,mu_10=mu_10_out,mu_20=mu_20_out,
                      type=type,link=link,
                      sigma_11_sq=sigma_11_sq_out,sigma_21_sq=sigma_21_sq_out,sigma_10_sq=sigma_10_sq_out,sigma_20_sq=sigma_20_sq_out,
                      pi_11=pi_11_out,pi_21=pi_21_out,pi_10=pi_10_out,pi_20=pi_20_out,
                      expit_beta_11=expit_beta_11_out,expit_beta_10=expit_beta_10_out,expit_beta_21=expit_beta_21_out,expit_beta_20=expit_beta_20_out)
    n_known <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1_hat,mu_0=mu_0_hat,
                       mu_11=mu_11_out,mu_21=mu_21_out,mu_10=mu_10_out,mu_20=mu_20_out,
                       type=type,link=link,
                       sigma_11_sq=sigma_11_sq_out,sigma_21_sq=sigma_21_sq_out,sigma_10_sq=sigma_10_sq_out,sigma_20_sq=sigma_20_sq_out,
                       pi_11=pi_11_out,pi_21=pi_21_out,pi_10=pi_10_out,pi_20=pi_20_out,
                       expit_beta_11=expit_beta_11_out,expit_beta_10=expit_beta_10_out,expit_beta_21=expit_beta_21_out,expit_beta_20=expit_beta_20_out)
    n_approx <- n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1_hat,mu_0=mu_0_hat,
                        type=type,link=link,
                        sigma_y_sq=sigma_y_sq_out,
                        pi_11=pi_11_out,pi_21=pi_21_out,pi_10=pi_10_out,pi_20=pi_20_out,
                        expit_beta_11=expit_beta_11_out,expit_beta_10=expit_beta_10_out,expit_beta_21=expit_beta_21_out,expit_beta_20=expit_beta_20_out)
    mu_1 <- pi_11*mu_11 + (1-pi_11)*mu_21
    mu_0 <- pi_11*mu_10 + (1-pi_11)*mu_20
    if (link=="identity") {
      n_standard_out <- (((mu_1_hat-mu_0_hat)^2)/((mu_1-mu_0)^2))*n_standard_out
      n_IPRW <- (((mu_1_hat-mu_0_hat)^2)/((mu_1-mu_0)^2))*n_IPRW
      n_known <- (((mu_1_hat-mu_0_hat)^2)/((mu_1-mu_0)^2))*n_known
      n_approx <- (((mu_1_hat-mu_0_hat)^2)/((mu_1-mu_0)^2))*n_approx
    }
    if (link=="logit") {
      n_standard_out <- (((qlogis(mu_1_hat)-qlogis(mu_0_hat))^2)/((qlogis(mu_1)-qlogis(mu_0))^2))*n_standard_out
      n_IPRW <- (((qlogis(mu_1_hat)-qlogis(mu_0_hat))^2)/((qlogis(mu_1)-qlogis(mu_0))^2))*n_IPRW
      n_known <- (((qlogis(mu_1_hat)-qlogis(mu_0_hat))^2)/((qlogis(mu_1)-qlogis(mu_0))^2))*n_known
      n_approx <- (((qlogis(mu_1_hat)-qlogis(mu_0_hat))^2)/((qlogis(mu_1)-qlogis(mu_0))^2))*n_approx
    }
    # n_standard
    if (!is.na(n_standard_out) & !is.nan(n_standard_out)) {
      Z <- c(rep(0,ceiling(n_standard_out*(1-kappa))),rep(1,ceiling(n_standard_out*kappa)))
      X <- rbinom(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)),1,1-pi_11)+1
      R <- NULL
      R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
      R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
      R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
      R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
      Y <- NULL
      if (type=="continuous" && is.null(delta)) {
        Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
        Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
        Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
        Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
      }
      if (type=="binary" && is.null(delta)) {
        Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
        Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
        Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
        Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
      }
      if (type=="continuous" && !is.null(delta)) {
        K <- NULL
        K[Z==0] <- rep(seq(1,ceiling(ceiling(n_standard_out*(1-kappa))/m)),each=m)[1:ceiling(n_standard_out*(1-kappa))]
        K[Z==1] <- rep(seq(1,ceiling(ceiling(n_standard_out*kappa)/m)),each=m)[1:ceiling(n_standard_out*kappa)] + ceiling(ceiling(n_standard_out*(1-kappa))/m)
        mu_1 <- (mu_11+mu_21)/2
        mu_0 <- (mu_10+mu_20)/2
        sigma_y_sq <- (0.5*(sigma_11_sq+(mu_11-mu_1)^2)+0.5*(sigma_21_sq+(mu_21-mu_1)^2) +
                       0.5*(sigma_10_sq+(mu_10-mu_0)^2)+0.5*(sigma_20_sq+(mu_20-mu_0)^2))/2 
        cluster_effect <- rnorm(ceiling(ceiling(n_standard_out*kappa)/m)+ceiling(ceiling(n_standard_out*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
        zeta <- NULL
        zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_standard_out*(1-kappa))/m)],each=m)[1:ceiling(n_standard_out*(1-kappa))]
        zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_standard_out*kappa)/m)+1):(ceiling(ceiling(n_standard_out*kappa)/m)+ceiling(ceiling(n_standard_out*(1-kappa))/m))],each=m)[1:ceiling(n_standard_out*kappa)]
        Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==1] 
        Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==2] 
        Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==1] 
        Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==2] 
      }
      # method of Qaqish, Biometrika 2003
      if (type=="binary" && !is.null(delta)) {
        K <- NULL
        K[Z==0] <- rep(seq(1,ceiling(ceiling(n_standard_out*(1-kappa))/m)),each=m)[1:ceiling(n_standard_out*(1-kappa))]
        K[Z==1] <- rep(seq(1,ceiling(ceiling(n_standard_out*kappa)/m)),each=m)[1:ceiling(n_standard_out*kappa)] + ceiling(ceiling(n_standard_out*(1-kappa))/m)
        mu_vec <- NULL
        mu_vec[Z==1 & X==1] <- mu_11
        mu_vec[Z==1 & X==2] <- mu_21
        mu_vec[Z==0 & X==1] <- mu_10
        mu_vec[Z==0 & X==2] <- mu_20 
        for (k in 1:(ceiling(ceiling(n_standard_out*kappa)/m)+ceiling(ceiling(n_standard_out*(1-kappa))/m))) {
          mu_k <- mu_vec[K==k]
          Y_k <- NULL
          Y_k[1] <- rbinom(1,1,mu_k[1]) 
          if (sum(K==k)>1) {
            for (j in 2:sum(K==k)) {
              j1 <- j-1
              res <- Y_k[1:j1] - mu_k[1:j1] # residuals
              v_part <- sqrt(mu_k[j]*(1-mu_k[j]))/sqrt(mu_k[1:j1]*(1-mu_k[1:j1])) # variance part (v_ii/v_jj)^(1/2)
              lambda <- mu_k[j] + sum(res*(delta/(1+(j-2)*delta))*v_part) # cond.mean
              lambda <- min(lambda,1)
              lambda <- max(lambda,0)
              Y_k[j] <- rbinom(1,1,lambda)
            }
          }
          Y <- c(Y,Y_k)
        }
      }
      # Equation 3
      mu_hat_1_eqn3 <- sum(R*Z*Y)/sum(R*Z)
      mu_hat_0_eqn3 <- sum(R*(1-Z)*Y)/sum(R*(1-Z))
      if (link=="identity" && is.null(delta)) {
        est_eqn3 <- mu_hat_1_eqn3 - mu_hat_0_eqn3 
        var_eqn3 <- sum(R*Z*(Y-mu_hat_1_eqn3)^2)/(sum(R*Z)^2) + 
          sum(R*(1-Z)*(Y-mu_hat_0_eqn3)^2)/(sum(R*(1-Z))^2)
      }
      if (link=="logit" && is.null(delta)) {
        est_eqn3 <- qlogis(mu_hat_1_eqn3) - qlogis(mu_hat_0_eqn3)
        var_eqn3 <- sum(R*Z*(Y-mu_hat_1_eqn3)^2)/((sum(R*Z)*mu_hat_1_eqn3*(1-mu_hat_1_eqn3))^2) + 
          sum(R*(1-Z)*(Y-mu_hat_0_eqn3)^2)/((sum(R*(1-Z))*mu_hat_0_eqn3*(1-mu_hat_0_eqn3))^2)
      }
      if (link=="identity" && !is.null(delta)) {
        est_eqn3 <- mu_hat_1_eqn3 - mu_hat_0_eqn3 
        out1 <- NULL
        out0 <- NULL
        for (k in 1:max(K)) {
          out1[k] <- sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn3))^2
          out0[k] <- sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn3))^2
        }
        var_eqn3 <- sum(out1)/(sum(R*Z)^2) + 
          sum(out0)/(sum(R*(1-Z))^2)
      }
      if (link=="logit" && !is.null(delta)) {
        est_eqn3 <- qlogis(mu_hat_1_eqn3) - qlogis(mu_hat_0_eqn3)
        out1 <- NULL
        out0 <- NULL
        for (k in 1:max(K)) {
          out1[k] <- sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn3))^2
          out0[k] <- sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn3))^2
        }
        var_eqn3 <- sum(out1)/((sum(R*Z)*mu_hat_1_eqn3*(1-mu_hat_1_eqn3))^2) + 
          sum(out0)/((sum(R*(1-Z))*mu_hat_0_eqn3*(1-mu_hat_0_eqn3))^2)
      }
      wald_eqn3 <- est_eqn3/sqrt(var_eqn3)
      reject_eqn3 <- wald_eqn3>qnorm(0.975)
      # Equation 4
      X1 <- ifelse(X==1,1,0)
      X2 <- ifelse(X==2,1,0)
      model_1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],family="binomial")
      model_0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],family="binomial")
      e <- NULL
      e[Z==1 & X==1] <- plogis(model_1$coefficients[1])
      e[Z==1 & X==2] <- plogis(model_1$coefficients[2])
      e[Z==0 & X==1] <- plogis(model_0$coefficients[1])
      e[Z==0 & X==2] <- plogis(model_0$coefficients[2])
      mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
      mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
      if (link=="identity" && is.null(delta)) {
        est_standard_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
        A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
        B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
        B_hat <- B_comp/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)))
        var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)))
      }
      if (link=="logit" && is.null(delta)) {
        est_standard_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
        A_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
        B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
        B_hat <- B_comp/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)))
        var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)))
      }
      if (link=="identity" && !is.null(delta)) {
        est_standard_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
        A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
        B_C_comp <- matrix(0,ncol=6,nrow=6)
        for (k in 1:max(K)) {
          u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
          outer <- u_k %*% t(u_k)
          B_C_comp <- B_C_comp + outer
        }
        B_C_hat <- B_C_comp/max(K)
        var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
      }
      if (link=="logit" && !is.null(delta)) {
        est_standard_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
        A_C_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
        B_C_comp <- matrix(0,ncol=6,nrow=6)
        for (k in 1:max(K)) {
          u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
          outer <- u_k %*% t(u_k)
          B_C_comp <- B_C_comp + outer
        }
        B_C_hat <- B_C_comp/max(K)
        var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
      }
      wald_eqn4 <- est_standard_eqn4/sqrt(var_eqn4)
      reject_standard_eqn4 <- wald_eqn4>qnorm(0.975)
    }
    if (is.na(n_standard_out)|is.nan(n_standard_out)) {  
      est_eqn3 <- NA
      est_standard_eqn4 <- NA
      reject_eqn3 <- FALSE
      reject_standard_eqn4 <- FALSE
    }
    # n_IPRW
    if (!is.na(n_IPRW) & !is.nan(n_IPRW)) {
      Z <- c(rep(0,ceiling(n_IPRW*(1-kappa))),rep(1,ceiling(n_IPRW*kappa)))
      X <- rbinom(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa),1,1-pi_11)+1
      R <- NULL
      R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
      R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
      R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
      R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
      Y <- NULL
      if (type=="continuous" && is.null(delta)) {
        Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
        Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
        Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
        Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
      }
      if (type=="binary" && is.null(delta)) {
        Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
        Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
        Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
        Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
      }
      if (type=="continuous" && !is.null(delta)) {
        K <- NULL
        K[Z==0] <- rep(seq(1,ceiling(ceiling(n_IPRW*(1-kappa))/m)),each=m)[1:ceiling(n_IPRW*(1-kappa))]
        K[Z==1] <- rep(seq(1,ceiling(ceiling(n_IPRW*kappa)/m)),each=m)[1:ceiling(n_IPRW*kappa)] + ceiling(ceiling(n_IPRW*(1-kappa))/m)
        mu_1 <- (mu_11+mu_21)/2
        mu_0 <- (mu_10+mu_20)/2
        sigma_y_sq <- (0.5*(sigma_11_sq+(mu_11-mu_1)^2)+0.5*(sigma_21_sq+(mu_21-mu_1)^2) +
                       0.5*(sigma_10_sq+(mu_10-mu_0)^2)+0.5*(sigma_20_sq+(mu_20-mu_0)^2))/2 
        cluster_effect <- rnorm(ceiling(ceiling(n_IPRW*kappa)/m)+ceiling(ceiling(n_IPRW*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
        zeta <- NULL
        zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_IPRW*(1-kappa))/m)],each=m)[1:ceiling(n_IPRW*(1-kappa))]
        zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_IPRW*(1-kappa))/m)+1):(ceiling(ceiling(n_IPRW*kappa)/m)+ceiling(ceiling(n_IPRW*(1-kappa))/m))],each=m)[1:ceiling(n_IPRW*kappa)]
        Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==1] 
        Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==2] 
        Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==1] 
        Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==2] 
      }
      # method of Qaqish, Biometrika 2003
      if (type=="binary" && !is.null(delta)) {
        K <- NULL
        K[Z==0] <- rep(seq(1,ceiling(ceiling(n_IPRW*(1-kappa))/m)),each=m)[1:ceiling(n_IPRW*(1-kappa))]
        K[Z==1] <- rep(seq(1,ceiling(ceiling(n_IPRW*kappa)/m)),each=m)[1:ceiling(n_IPRW*kappa)] + ceiling(ceiling(n_IPRW*(1-kappa))/m)
        mu_vec <- NULL
        mu_vec[Z==1 & X==1] <- mu_11
        mu_vec[Z==1 & X==2] <- mu_21
        mu_vec[Z==0 & X==1] <- mu_10
        mu_vec[Z==0 & X==2] <- mu_20 
        for (k in 1:(ceiling(ceiling(n_IPRW*kappa)/m)+ceiling(ceiling(n_IPRW*(1-kappa))/m))) {
          mu_k <- mu_vec[K==k]
          Y_k <- NULL
          Y_k[1] <- rbinom(1,1,mu_k[1]) 
          if (sum(K==k)>1) {
            for (j in 2:sum(K==k)) {
              j1 <- j-1
              res <- Y_k[1:j1] - mu_k[1:j1] # residuals
              v_part <- sqrt(mu_k[j]*(1-mu_k[j]))/sqrt(mu_k[1:j1]*(1-mu_k[1:j1])) # variance part (v_ii/v_jj)^(1/2)
              lambda <- mu_k[j] + sum(res*(delta/(1+(j-2)*delta))*v_part) # cond.mean
              lambda <- min(lambda,1)
              lambda <- max(lambda,0)
              Y_k[j] <- rbinom(1,1,lambda)
            }
          }
          Y <- c(Y,Y_k)
        }
      }
      # Equation 4
      X1 <- ifelse(X==1,1,0)
      X2 <- ifelse(X==2,1,0)
      model_1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],family="binomial")
      model_0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],family="binomial")
      e <- NULL
      e[Z==1 & X==1] <- plogis(model_1$coefficients[1])
      e[Z==1 & X==2] <- plogis(model_1$coefficients[2])
      e[Z==0 & X==1] <- plogis(model_0$coefficients[1])
      e[Z==0 & X==2] <- plogis(model_0$coefficients[2])
      mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
      mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
      if (link=="identity" && is.null(delta)) {
        est_IPRW_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
        A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))),
        nrow=6,ncol=6,byrow=T)
        B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
        B_hat <- B_comp/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
        var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
      }
      if (link=="logit" && is.null(delta)) {
        est_IPRW_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
        A_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa)),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))),
        nrow=6,ncol=6,byrow=T)
        B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
        B_hat <- B_comp/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
        var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
      }
      if (link=="identity" && !is.null(delta)) {
        est_IPRW_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
        A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
        B_C_comp <- matrix(0,ncol=6,nrow=6)
        for (k in 1:max(K)) {
          u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
          outer <- u_k %*% t(u_k)
          B_C_comp <- B_C_comp + outer
        }
        B_C_hat <- B_C_comp/max(K)
        var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
      }
      if (link=="logit" && !is.null(delta)) {
        est_IPRW_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
        A_C_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
        B_C_comp <- matrix(0,ncol=6,nrow=6)
        for (k in 1:max(K)) {
          u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
          outer <- u_k %*% t(u_k)
          B_C_comp <- B_C_comp + outer
        }
        B_C_hat <- B_C_comp/max(K)
        var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
      }
      wald_eqn4 <- est_IPRW_eqn4/sqrt(var_eqn4)
      reject_IPRW_eqn4 <- wald_eqn4>qnorm(0.975)
    }
    if (is.na(n_IPRW)|is.nan(n_IPRW)) {
      est_IPRW_eqn4 <- NA
      reject_IPRW_eqn4 <- FALSE
    }
    # n_known
    if (!is.na(n_known) & !is.nan(n_known)) {
      Z <- c(rep(0,ceiling(n_known*(1-kappa))),rep(1,ceiling(n_known*kappa)))
      X <- rbinom(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)),1,1-pi_11)+1
      R <- NULL
      R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
      R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
      R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
      R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
      Y <- NULL
      if (type=="continuous" && is.null(delta)) {
        Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
        Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
        Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
        Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
      }
      if (type=="binary" && is.null(delta)) {
        Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
        Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
        Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
        Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
      }
      if (type=="continuous" && !is.null(delta)) {
        K <- NULL
        K[Z==0] <- rep(seq(1,ceiling(ceiling(n_known*(1-kappa))/m)),each=m)[1:ceiling(n_known*(1-kappa))]
        K[Z==1] <- rep(seq(1,ceiling(ceiling(n_known*kappa)/m)),each=m)[1:ceiling(n_known*kappa)] + ceiling(ceiling(n_known*(1-kappa))/m)
        mu_1 <- (mu_11+mu_21)/2
        mu_0 <- (mu_10+mu_20)/2
        sigma_y_sq <- (0.5*(sigma_11_sq+(mu_11-mu_1)^2)+0.5*(sigma_21_sq+(mu_21-mu_1)^2) +
                       0.5*(sigma_10_sq+(mu_10-mu_0)^2)+0.5*(sigma_20_sq+(mu_20-mu_0)^2))/2 
        cluster_effect <- rnorm(ceiling(ceiling(n_known*kappa)/m)+ceiling(ceiling(n_known*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
        zeta <- NULL
        zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_known*(1-kappa))/m)],each=m)[1:ceiling(n_known*(1-kappa))]
        zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_known*(1-kappa))/m)+1):(ceiling(ceiling(n_known*kappa)/m)+ceiling(ceiling(n_known*(1-kappa))/m))],each=m)[1:ceiling(n_known*kappa)]
        Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==1] 
        Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==2] 
        Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==1] 
        Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==2] 
      }
      # method of Qaqish, Biometrika 2003
      if (type=="binary" && !is.null(delta)) {
        K <- NULL
        K[Z==0] <- rep(seq(1,ceiling(ceiling(n_known*(1-kappa))/m)),each=m)[1:ceiling(n_known*(1-kappa))]
        K[Z==1] <- rep(seq(1,ceiling(ceiling(n_known*kappa)/m)),each=m)[1:ceiling(n_known*kappa)] + ceiling(ceiling(n_known*(1-kappa))/m)
        mu_vec <- NULL
        mu_vec[Z==1 & X==1] <- mu_11
        mu_vec[Z==1 & X==2] <- mu_21
        mu_vec[Z==0 & X==1] <- mu_10
        mu_vec[Z==0 & X==2] <- mu_20 
        for (k in 1:(ceiling(ceiling(n_known*kappa)/m)+ceiling(ceiling(n_known*(1-kappa))/m))) {
          mu_k <- mu_vec[K==k]
          Y_k <- NULL
          Y_k[1] <- rbinom(1,1,mu_k[1]) 
          if (sum(K==k)>1) {
            for (j in 2:sum(K==k)) {
              j1 <- j-1
              res <- Y_k[1:j1] - mu_k[1:j1] # residuals
              v_part <- sqrt(mu_k[j]*(1-mu_k[j]))/sqrt(mu_k[1:j1]*(1-mu_k[1:j1])) # variance part (v_ii/v_jj)^(1/2)
              lambda <- mu_k[j] + sum(res*(delta/(1+(j-2)*delta))*v_part) # cond.mean
              lambda <- min(lambda,1)
              lambda <- max(lambda,0)
              Y_k[j] <- rbinom(1,1,lambda)
            }
          }
          Y <- c(Y,Y_k)
        }
      }
      # Equation 4
      X1 <- ifelse(X==1,1,0)
      X2 <- ifelse(X==2,1,0)
      model_1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],family="binomial")
      model_0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],family="binomial")
      e <- NULL
      e[Z==1 & X==1] <- plogis(model_1$coefficients[1])
      e[Z==1 & X==2] <- plogis(model_1$coefficients[2])
      e[Z==0 & X==1] <- plogis(model_0$coefficients[1])
      e[Z==0 & X==2] <- plogis(model_0$coefficients[2])
      mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
      mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
      if (link=="identity" && is.null(delta)) {
        est_known_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
        A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
        B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
        B_hat <- B_comp/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))
        var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))
      }
      if (link=="logit" && is.null(delta)) {
        est_known_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
        A_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
        B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
        B_hat <- B_comp/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))
        var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)))
      }
      if (link=="identity" && !is.null(delta)) {
        est_known_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
        A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
        B_C_comp <- matrix(0,ncol=6,nrow=6)
        for (k in 1:max(K)) {
          u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
          outer <- u_k %*% t(u_k)
          B_C_comp <- B_C_comp + outer
        }
        B_C_hat <- B_C_comp/max(K)
        var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
      }
      if (link=="logit" && !is.null(delta)) {
        est_known_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
        A_C_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
        B_C_comp <- matrix(0,ncol=6,nrow=6)
        for (k in 1:max(K)) {
          u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
          outer <- u_k %*% t(u_k)
          B_C_comp <- B_C_comp + outer
        }
        B_C_hat <- B_C_comp/max(K)
        var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
      }
      wald_eqn4 <- est_known_eqn4/sqrt(var_eqn4)
      reject_known_eqn4 <- wald_eqn4>qnorm(0.975)
    }
    if (is.na(n_known)| is.nan(n_known)) {
      est_known_eqn4 <- NA
      reject_known_eqn4 <- FALSE
    }
    # n_approx
    if (!is.na(n_approx) & !is.nan(n_approx)) {
      Z <- c(rep(0,ceiling(n_approx*(1-kappa))),rep(1,ceiling(n_approx*kappa)))
      X <- rbinom(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)),1,1-pi_11)+1
      R <- NULL
      R[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,expit_beta_11)
      R[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,expit_beta_21)
      R[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,expit_beta_10)
      R[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,expit_beta_20)
      Y <- NULL
      if (type=="continuous" && is.null(delta)) {
        Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq))
        Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq))
        Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq))
        Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq))
      }
      if (type=="binary" && is.null(delta)) {
        Y[Z==1 & X==1] <- rbinom(sum(Z==1 & X==1),1,mu_11)
        Y[Z==1 & X==2] <- rbinom(sum(Z==1 & X==2),1,mu_21)
        Y[Z==0 & X==1] <- rbinom(sum(Z==0 & X==1),1,mu_10)
        Y[Z==0 & X==2] <- rbinom(sum(Z==0 & X==2),1,mu_20)
      }
      if (type=="continuous" && !is.null(delta)) {
        K <- NULL
        K[Z==0] <- rep(seq(1,ceiling(ceiling(n_approx*(1-kappa))/m)),each=m)[1:ceiling(n_approx*(1-kappa))]
        K[Z==1] <- rep(seq(1,ceiling(ceiling(n_approx*kappa)/m)),each=m)[1:ceiling(n_approx*kappa)] + ceiling(ceiling(n_approx*(1-kappa))/m)
        mu_1 <- (mu_11+mu_21)/2
        mu_0 <- (mu_10+mu_20)/2
        sigma_y_sq <- (0.5*(sigma_11_sq+(mu_11-mu_1)^2)+0.5*(sigma_21_sq+(mu_21-mu_1)^2) +
                       0.5*(sigma_10_sq+(mu_10-mu_0)^2)+0.5*(sigma_20_sq+(mu_20-mu_0)^2))/2 
        cluster_effect <- rnorm(ceiling(ceiling(n_approx*kappa)/m)+ceiling(ceiling(n_approx*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
        zeta <- NULL
        zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_approx*(1-kappa))/m)],each=m)[1:ceiling(n_approx*(1-kappa))]
        zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_approx*(1-kappa))/m)+1):(ceiling(ceiling(n_approx*kappa)/m)+ceiling(ceiling(n_approx*(1-kappa))/m))],each=m)[1:ceiling(n_approx*kappa)]
        Y[Z==1 & X==1] <- rnorm(sum(Z==1 & X==1),mu_11,sqrt(sigma_11_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==1] 
        Y[Z==1 & X==2] <- rnorm(sum(Z==1 & X==2),mu_21,sqrt(sigma_21_sq-delta*sigma_y_sq)) + zeta[Z==1 & X==2] 
        Y[Z==0 & X==1] <- rnorm(sum(Z==0 & X==1),mu_10,sqrt(sigma_10_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==1] 
        Y[Z==0 & X==2] <- rnorm(sum(Z==0 & X==2),mu_20,sqrt(sigma_20_sq-delta*sigma_y_sq)) + zeta[Z==0 & X==2] 
      }
      # method of Qaqish, Biometrika 2003
      if (type=="binary" && !is.null(delta)) {
        K <- NULL
        K[Z==0] <- rep(seq(1,ceiling(ceiling(n_approx*(1-kappa))/m)),each=m)[1:ceiling(n_approx*(1-kappa))]
        K[Z==1] <- rep(seq(1,ceiling(ceiling(n_approx*kappa)/m)),each=m)[1:ceiling(n_approx*kappa)] + ceiling(ceiling(n_approx*(1-kappa))/m)
        mu_vec <- NULL
        mu_vec[Z==1 & X==1] <- mu_11
        mu_vec[Z==1 & X==2] <- mu_21
        mu_vec[Z==0 & X==1] <- mu_10
        mu_vec[Z==0 & X==2] <- mu_20 
        for (k in 1:(ceiling(ceiling(n_approx*kappa)/m)+ceiling(ceiling(n_approx*(1-kappa))/m))) {
          mu_k <- mu_vec[K==k]
          Y_k <- NULL
          Y_k[1] <- rbinom(1,1,mu_k[1]) 
          if (sum(K==k)>1) {
            for (j in 2:sum(K==k)) {
              j1 <- j-1
              res <- Y_k[1:j1] - mu_k[1:j1] # residuals
              v_part <- sqrt(mu_k[j]*(1-mu_k[j]))/sqrt(mu_k[1:j1]*(1-mu_k[1:j1])) # variance part (v_ii/v_jj)^(1/2)
              lambda <- mu_k[j] + sum(res*(delta/(1+(j-2)*delta))*v_part) # cond.mean
              lambda <- min(lambda,1)
              lambda <- max(lambda,0)
              Y_k[j] <- rbinom(1,1,lambda)
            }
          }
          Y <- c(Y,Y_k)
        }
      }
      # Equation 4
      X1 <- ifelse(X==1,1,0)
      X2 <- ifelse(X==2,1,0)
      model_1 <- glm(R[Z==1]~0+X1[Z==1]+X2[Z==1],family="binomial")
      model_0 <- glm(R[Z==0]~0+X1[Z==0]+X2[Z==0],family="binomial")
      e <- NULL
      e[Z==1 & X==1] <- plogis(model_1$coefficients[1])
      e[Z==1 & X==2] <- plogis(model_1$coefficients[2])
      e[Z==0 & X==1] <- plogis(model_0$coefficients[1])
      e[Z==0 & X==2] <- plogis(model_0$coefficients[2])
      mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
      mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
      if (link=="identity" && is.null(delta)) {
        est_approx_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
        A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
        B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
        B_hat <- B_comp/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
        var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
      }
      if (link=="logit" && is.null(delta)) {
        est_approx_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
        A_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),
          0,0,sum(Z*e*(1-e)*X1)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa))),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))),
        nrow=6,ncol=6,byrow=T)
        B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X1/e),0,sum(Z*(R-e)^2*X1),0,0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X2/e),0,0,sum(Z*(R-e)^2*X2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X1/e),0,0,sum((1-Z)*(R-e)^2*X1),0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X2/e),0,0,0,sum((1-Z)*(R-e)^2*X2)),
                       nrow=6,ncol=6,byrow=T)
        B_hat <- B_comp/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
        var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
      }
      if (link=="identity" && !is.null(delta)) {
        est_approx_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4
        A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
        B_C_comp <- matrix(0,ncol=6,nrow=6)
        for (k in 1:max(K)) {
          u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
          outer <- u_k %*% t(u_k)
          B_C_comp <- B_C_comp + outer
        }
        B_C_hat <- B_C_comp/max(K)
        var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
      }
      if (link=="logit" && !is.null(delta)) {
        est_approx_eqn4 <- qlogis(mu_hat_1_eqn4) - qlogis(mu_hat_0_eqn4)
        A_C_hat <- matrix(
        c(sum(R*Z*mu_hat_1_eqn4*(1-mu_hat_1_eqn4)/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X1/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X2/e)/max(K),0,0,
          0,sum(R*(1-Z)*mu_hat_0_eqn4*(1-mu_hat_0_eqn4)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X1/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X2/e)/max(K),
          0,0,sum(Z*e*(1-e)*X1)/max(K),0,0,0,
          0,0,0,sum(Z*e*(1-e)*X2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e)*X1)/max(K),0,
          0,0,0,0,0,sum((1-Z)*e*(1-e)*X2)/max(K)),
        nrow=6,ncol=6,byrow=T)
        B_C_comp <- matrix(0,ncol=6,nrow=6)
        for (k in 1:max(K)) {
          u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*X1[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X2[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X1[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X2[K==k]*(R[K==k]-e[K==k])))
          outer <- u_k %*% t(u_k)
          B_C_comp <- B_C_comp + outer
        }
        B_C_hat <- B_C_comp/max(K)
        var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
        var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
      }
      wald_eqn4 <- est_approx_eqn4/sqrt(var_eqn4)
      reject_approx_eqn4 <- wald_eqn4>qnorm(0.975)
    }
    if (is.na(n_approx)|is.nan(n_approx)) {
      est_approx_eqn4 <- NA
      reject_approx_eqn4 <- FALSE
    }
    return(c(n_standard_out,n_IPRW,n_known,n_approx,
             est_eqn3,est_standard_eqn4,est_IPRW_eqn4,est_known_eqn4,est_approx_eqn4,
             reject_eqn3,reject_standard_eqn4,reject_IPRW_eqn4,reject_known_eqn4,reject_approx_eqn4))
  }
  results <- colMeans(para,na.rm=TRUE)
  names(results) <- c("n_standard","n_IPRW","n_known","n_approx",
                      "g(mu_hat_1)-g(mu_hat_0) (eqn 3)",
                      "g(mu_hat_1)-g(mu_hat_0) (standard eqn 4)",
                      "g(mu_hat_1)-g(mu_hat_0) (IPRW eqn 4)",
                      "g(mu_hat_1)-g(mu_hat_0) (known eqn 4)",
                      "g(mu_hat_1)-g(mu_hat_0) (approx eqn 4)",
                      "empirical power (eqn 3)",
                      "empirical power (standard eqn 4)",
                      "empirical power (IPRW eqn 4)",
                      "empirical power (known eqn 4)",
                      "empirical power (approx eqn 4)")
  return(rbind(para,results))
}

# Weighting for a fully observed auxiliary continuous variable
sim_power_cont <- function(rep,n_standard,n_IPRW,n_known,n_approx,
                           mu_1,mu_0,mu_x,
                           sigma_x_sq,sigma_y_sq,rho,
                           beta_11,beta_01,beta_10,beta_00,
                           delta=NULL,m=NULL,kappa)
{
  para <- foreach(i=1:rep, .combine=rbind) %dorng% {
    # n_standard
    Z <- c(rep(0,ceiling(n_standard*(1-kappa))),rep(1,ceiling(n_standard*kappa)))
    X <- rnorm(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa),0,1)
    R <- NULL
    R[Z==1] <- rbinom(sum(Z==1),1,plogis(beta_01+beta_11*X[Z==1]))
    R[Z==0] <- rbinom(sum(Z==0),1,plogis(beta_00+beta_10*X[Z==0]))
    Y <- NULL
    if (is.null(delta)) {
      Y[Z==1] <- rnorm(sum(Z==1),mu_1+rho*sqrt(sigma_y_sq)*X[Z==1],sqrt(sigma_y_sq*(1-rho^2)))
      Y[Z==0] <- rnorm(sum(Z==0),mu_0+rho*sqrt(sigma_y_sq)*X[Z==0],sqrt(sigma_y_sq*(1-rho^2)))
    }
    if (!is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_standard*(1-kappa))/m)),each=m)[1:ceiling(n_standard*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_standard*kappa)/m)),each=m)[1:ceiling(n_standard*kappa)] + ceiling(ceiling(n_standard*(1-kappa))/m)
      cluster_effect <- rnorm(ceiling(ceiling(n_standard*(1-kappa))/m)+ceiling(ceiling(n_standard*kappa)/m),0,sqrt(delta*sigma_y_sq))
      zeta <- NULL
      zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_standard*(1-kappa))/m)],each=m)[1:ceiling(n_standard*(1-kappa))]
      zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_standard*kappa)/m)+1):(ceiling(ceiling(n_standard*kappa)/m)+ceiling(ceiling(n_standard*(1-kappa)/m)))],each=m)[1:ceiling(n_standard*kappa)]
      Y[Z==1] <- rnorm(sum(Z==1),mu_1+rho*sqrt(sigma_y_sq)*X[Z==1],sqrt(sigma_y_sq*(1-delta-rho^2))) + zeta[Z==1] 
      Y[Z==0] <- rnorm(sum(Z==0),mu_0+rho*sqrt(sigma_y_sq)*X[Z==0],sqrt(sigma_y_sq*(1-delta-rho^2))) + zeta[Z==0] 
    }
    # Equation 3
    mu_hat_1_eqn3 <- sum(R*Z*Y)/sum(R*Z)
    mu_hat_0_eqn3 <- sum(R*(1-Z)*Y)/sum(R*(1-Z))
    est_eqn3 <- mu_hat_1_eqn3 - mu_hat_0_eqn3 
    if (is.null(delta)) {
      var_eqn3 <- sum(R*Z*(Y-mu_hat_1_eqn3)^2)/(sum(R*Z)^2) + 
                  sum(R*(1-Z)*(Y-mu_hat_0_eqn3)^2)/(sum(R*(1-Z))^2)
    }
    if (!is.null(delta)) {
      out1 <- NULL
      out0 <- NULL
      for (k in 1:max(K)) {
        out1[k] <- sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn3))^2
        out0[k] <- sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn3))^2
      }
      var_eqn3 <- sum(out1)/(sum(R*Z)^2) + 
                  sum(out0)/(sum(R*(1-Z))^2)
    }
    wald_eqn3 <- est_eqn3/sqrt(var_eqn3)
    reject_eqn3 <- wald_eqn3>qnorm(0.975)
    # Equation 4
    model_1 <- glm(R[Z==1]~X[Z==1],family="binomial")
    model_0 <- glm(R[Z==0]~X[Z==0],family="binomial")
    e <- NULL
    e[Z==1] <- plogis(model_1$coefficients[1]+model_1$coefficients[2]*X[Z==1])
    e[Z==0] <- plogis(model_0$coefficients[1]+model_0$coefficients[2]*X[Z==0])
    mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
    mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
    est_standard_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4 
    if (is.null(delta)) {
      A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)/e)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X/e)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)/e)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X/e)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),
          0,0,sum(Z*e*(1-e))/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),sum(Z*e*(1-e)*X)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),0,0,
          0,0,sum(Z*e*(1-e)*X)/(ceiling(n_standard*(1-kappa)+ceiling(n_standard*kappa))),sum(Z*e*(1-e)*X^2)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e))/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),sum((1-Z)*e*(1-e)*X)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),
          0,0,0,0,sum((1-Z)*e*(1-e)*X)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa)),sum((1-Z)*e*(1-e)*X^2)/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)/e),0,sum(Z*(R-e)^2),sum(Z*(R-e)^2*X),0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X/e),0,sum(Z*(R-e)^2*X),sum(Z*(R-e)^2*X^2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)/e),0,0,sum((1-Z)*(R-e)^2),sum((1-Z)*(R-e)^2*X),
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X/e),0,0,sum((1-Z)*(R-e)^2*X),sum((1-Z)*(R-e)^2*X^2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_standard*(1-kappa))+ceiling(n_standard*kappa))
    }
    if (!is.null(delta)) {
      A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X/e)/max(K),
          0,0,sum(Z*e*(1-e))/max(K),sum(Z*e*(1-e)*X)/max(K),0,0,
          0,0,sum(Z*e*(1-e)*X)/max(K),sum(Z*e*(1-e)*X^2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e))/max(K),sum((1-Z)*e*(1-e)*X)/max(K),
          0,0,0,0,sum((1-Z)*e*(1-e)*X)/max(K),sum((1-Z)*e*(1-e)*X^2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    wald_eqn4 <- est_standard_eqn4/sqrt(var_eqn4)
    reject_standard_eqn4 <- wald_eqn4>qnorm(0.975)
    # n_IPRW
    Z <- c(rep(0,ceiling(n_IPRW*(1-kappa))),rep(1,ceiling(n_IPRW*kappa)))
    X <- rnorm(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa),0,1)
    R <- NULL
    R[Z==1] <- rbinom(sum(Z==1),1,plogis(beta_01+beta_11*X[Z==1]))
    R[Z==0] <- rbinom(sum(Z==0),1,plogis(beta_00+beta_10*X[Z==0]))
    Y <- NULL
    if (is.null(delta)) {
      Y[Z==1] <- rnorm(sum(Z==1),mu_1+rho*sqrt(sigma_y_sq)*X[Z==1],sqrt(sigma_y_sq*(1-rho^2)))
      Y[Z==0] <- rnorm(sum(Z==0),mu_0+rho*sqrt(sigma_y_sq)*X[Z==0],sqrt(sigma_y_sq*(1-rho^2)))
    }
    if (!is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_IPRW*(1-kappa))/m)),each=m)[1:ceiling(n_IPRW*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_IPRW*kappa)/m)),each=m)[1:ceiling(n_IPRW*kappa)] + ceiling(ceiling(n_IPRW*(1-kappa))/m)
      cluster_effect <- rnorm(ceiling(ceiling(n_IPRW*kappa)/m)+ceiling(ceiling(n_IPRW*(1-kappa)/m)),0,sqrt(delta*sigma_y_sq))
      zeta <- NULL
      zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_IPRW*(1-kappa))/m)],each=m)[1:ceiling(n_IPRW*(1-kappa))]
      zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_IPRW*(1-kappa))/m)+1):(ceiling(ceiling(n_IPRW*kappa)/m)+ceiling(ceiling(n_IPRW*(1-kappa))/m))],each=m)[1:ceiling(n_IPRW*kappa)]
      Y[Z==1] <- rnorm(sum(Z==1),mu_1+rho*sqrt(sigma_y_sq)*X[Z==1],sqrt(sigma_y_sq*(1-delta-rho^2))) + zeta[Z==1] 
      Y[Z==0] <- rnorm(sum(Z==0),mu_0+rho*sqrt(sigma_y_sq)*X[Z==0],sqrt(sigma_y_sq*(1-delta-rho^2))) + zeta[Z==0] 
    }
    # Equation 4
    model_1 <- glm(R[Z==1]~X[Z==1],family="binomial")
    model_0 <- glm(R[Z==0]~X[Z==0],family="binomial")
    e <- NULL
    e[Z==1] <- plogis(model_1$coefficients[1]+model_1$coefficients[2]*X[Z==1])
    e[Z==0] <- plogis(model_0$coefficients[1]+model_0$coefficients[2]*X[Z==0])
    mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
    mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
    est_IPRW_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4 
    if (is.null(delta)) {
      A_comp <- matrix(c(sum(R*Z/e),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X/e),0,0,
                         0,sum(R*(1-Z)/e),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X/e),
                         0,0,sum(Z*e*(1-e)),sum(Z*e*(1-e)*X),0,0,
                         0,0,sum(Z*e*(1-e)*X),sum(Z*e*(1-e)*X^2),0,0,
                         0,0,0,0,sum((1-Z)*e*(1-e)),sum((1-Z)*e*(1-e)*X),
                         0,0,0,0,sum((1-Z)*e*(1-e)*X),sum((1-Z)*e*(1-e)*X^2)),
                      nrow=6,ncol=6,byrow=T)
      A_hat <- A_comp/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)/e),0,sum(Z*(R-e)^2),sum(Z*(R-e)^2*X),0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X/e),0,sum(Z*(R-e)^2*X),sum(Z*(R-e)^2*X^2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)/e),0,0,sum((1-Z)*(R-e)^2),sum((1-Z)*(R-e)^2*X),
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X/e),0,0,sum((1-Z)*(R-e)^2*X),sum((1-Z)*(R-e)^2*X^2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_IPRW*(1-kappa))+ceiling(n_IPRW*kappa))
    }
    if (!is.null(delta)) {
      A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X/e)/max(K),
          0,0,sum(Z*e*(1-e))/max(K),sum(Z*e*(1-e)*X)/max(K),0,0,
          0,0,sum(Z*e*(1-e)*X)/max(K),sum(Z*e*(1-e)*X^2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e))/max(K),sum((1-Z)*e*(1-e)*X)/max(K),
          0,0,0,0,sum((1-Z)*e*(1-e)*X)/max(K),sum((1-Z)*e*(1-e)*X^2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    wald_eqn4 <- est_IPRW_eqn4/sqrt(var_eqn4)
    reject_IPRW_eqn4 <- wald_eqn4>qnorm(0.975)
    # n_known
    Z <- c(rep(0,ceiling(n_known*(1-kappa))),rep(1,ceiling(n_known*kappa)))
    X <- rnorm(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa),0,1)
    R <- NULL
    R[Z==1] <- rbinom(sum(Z==1),1,plogis(beta_01+beta_11*X[Z==1]))
    R[Z==0] <- rbinom(sum(Z==0),1,plogis(beta_00+beta_10*X[Z==0]))
    Y <- NULL
    if (is.null(delta)) {
      Y[Z==1] <- rnorm(sum(Z==1),mu_1+rho*sqrt(sigma_y_sq)*X[Z==1],sqrt(sigma_y_sq*(1-rho^2)))
      Y[Z==0] <- rnorm(sum(Z==0),mu_0+rho*sqrt(sigma_y_sq)*X[Z==0],sqrt(sigma_y_sq*(1-rho^2)))
    }
    if (!is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_known*(1-kappa))/m)),each=m)[1:ceiling(n_known*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_known*kappa)/m)),each=m)[1:ceiling(n_known*kappa)] + ceiling(ceiling(n_known*(1-kappa))/m)
      cluster_effect <- rnorm(ceiling(ceiling(n_known*kappa)/m)*2,0,sqrt(delta*sigma_y_sq))
      zeta <- NULL
      zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_known*(1-kappa))/m)],each=m)[1:ceiling(n_known*(1-kappa))]
      zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_known*(1-kappa))/m)+1):(ceiling(ceiling(n_known*(1-kappa))/m)+ceiling(ceiling(n_known*kappa)/m))],each=m)[1:ceiling(n_known*kappa)]
      Y[Z==1] <- rnorm(sum(Z==1),mu_1+rho*sqrt(sigma_y_sq)*X[Z==1],sqrt(sigma_y_sq*(1-delta-rho^2))) + zeta[Z==1] 
      Y[Z==0] <- rnorm(sum(Z==0),mu_0+rho*sqrt(sigma_y_sq)*X[Z==0],sqrt(sigma_y_sq*(1-delta-rho^2))) + zeta[Z==0] 
    }
    # Equation 4
    model_1 <- glm(R[Z==1]~X[Z==1],family="binomial")
    model_0 <- glm(R[Z==0]~X[Z==0],family="binomial")
    e <- NULL
    e[Z==1] <- plogis(model_1$coefficients[1]+model_1$coefficients[2]*X[Z==1])
    e[Z==0] <- plogis(model_0$coefficients[1]+model_0$coefficients[2]*X[Z==0])
    mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
    mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
    est_known_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4 
    if (is.null(delta)) {
      A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)/e)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X/e)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)/e)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X/e)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),
          0,0,sum(Z*e*(1-e))/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),sum(Z*e*(1-e)*X)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),0,0,
          0,0,sum(Z*e*(1-e)*X)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),sum(Z*e*(1-e)*X^2)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e))/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),sum((1-Z)*e*(1-e)*X)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),
          0,0,0,0,sum((1-Z)*e*(1-e)*X)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa)),sum((1-Z)*e*(1-e)*X^2)/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)/e),0,sum(Z*(R-e)^2),sum(Z*(R-e)^2*X),0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X/e),0,sum(Z*(R-e)^2*X),sum(Z*(R-e)^2*X^2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)/e),0,0,sum((1-Z)*(R-e)^2),sum((1-Z)*(R-e)^2*X),
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X/e),0,0,sum((1-Z)*(R-e)^2*X),sum((1-Z)*(R-e)^2*X^2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_known*(1-kappa))+ceiling(n_known*kappa))
    }
    if (!is.null(delta)) {
      A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X/e)/max(K),
          0,0,sum(Z*e*(1-e))/max(K),sum(Z*e*(1-e)*X)/max(K),0,0,
          0,0,sum(Z*e*(1-e)*X)/max(K),sum(Z*e*(1-e)*X^2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e))/max(K),sum((1-Z)*e*(1-e)*X)/max(K),
          0,0,0,0,sum((1-Z)*e*(1-e)*X)/max(K),sum((1-Z)*e*(1-e)*X^2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    wald_eqn4 <- est_known_eqn4/sqrt(var_eqn4)
    reject_known_eqn4 <- wald_eqn4>qnorm(0.975)
    # n_approx
    Z <- c(rep(0,ceiling(n_approx*(1-kappa))),rep(1,ceiling(n_approx*kappa)))
    X <- rnorm(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa),0,1)
    R <- NULL
    R[Z==1] <- rbinom(sum(Z==1),1,plogis(beta_01+beta_11*X[Z==1]))
    R[Z==0] <- rbinom(sum(Z==0),1,plogis(beta_00+beta_10*X[Z==0]))
    Y <- NULL
    if (is.null(delta)) {
      Y[Z==1] <- rnorm(sum(Z==1),mu_1+rho*sqrt(sigma_y_sq)*X[Z==1],sqrt(sigma_y_sq*(1-rho^2)))
      Y[Z==0] <- rnorm(sum(Z==0),mu_0+rho*sqrt(sigma_y_sq)*X[Z==0],sqrt(sigma_y_sq*(1-rho^2)))
    }
    if (!is.null(delta)) {
      K <- NULL
      K[Z==0] <- rep(seq(1,ceiling(ceiling(n_approx*(1-kappa))/m)),each=m)[1:ceiling(n_approx*(1-kappa))]
      K[Z==1] <- rep(seq(1,ceiling(ceiling(n_approx*kappa)/m)),each=m)[1:ceiling(n_approx*kappa)] + ceiling(ceiling(n_approx*(1-kappa))/m)
      cluster_effect <- rnorm(ceiling(ceiling(n_approx*kappa)/m)+ceiling(ceiling(n_approx*(1-kappa))/m),0,sqrt(delta*sigma_y_sq))
      zeta <- NULL
      zeta[Z==0] <- rep(cluster_effect[1:ceiling(ceiling(n_approx*(1-kappa))/m)],each=m)[1:ceiling(n_approx*(1-kappa))]
      zeta[Z==1] <- rep(cluster_effect[(ceiling(ceiling(n_approx*(1-kappa))/m)+1):(ceiling(ceiling(n_approx*kappa)/m)+ceiling(ceiling(n_approx*(1-kappa))/m))],each=m)[1:ceiling(n_approx*kappa)]
      Y[Z==1] <- rnorm(sum(Z==1),mu_1+rho*sqrt(sigma_y_sq)*X[Z==1],sqrt(sigma_y_sq*(1-delta-rho^2))) + zeta[Z==1] 
      Y[Z==0] <- rnorm(sum(Z==0),mu_0+rho*sqrt(sigma_y_sq)*X[Z==0],sqrt(sigma_y_sq*(1-delta-rho^2))) + zeta[Z==0] 
    }
    # Equation 4
    model_1 <- glm(R[Z==1]~X[Z==1],family="binomial")
    model_0 <- glm(R[Z==0]~X[Z==0],family="binomial")
    e <- NULL
    e[Z==1] <- plogis(model_1$coefficients[1]+model_1$coefficients[2]*X[Z==1])
    e[Z==0] <- plogis(model_0$coefficients[1]+model_0$coefficients[2]*X[Z==0])
    mu_hat_1_eqn4 <- sum(R*Z*Y/e)/sum(R*Z/e)
    mu_hat_0_eqn4 <- sum(R*(1-Z)*Y/e)/sum(R*(1-Z)/e)
    est_approx_eqn4 <- mu_hat_1_eqn4 - mu_hat_0_eqn4 
    if (is.null(delta)) {
      A_hat <- matrix(
        c(sum(R*Z/e)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)/e)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X/e)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),0,0,
          0,sum(R*(1-Z)/e)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)/e)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X/e)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),
          0,0,sum(Z*e*(1-e))/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),sum(Z*e*(1-e)*X)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),0,0,
          0,0,sum(Z*e*(1-e)*X)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),sum(Z*e*(1-e)*X^2)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e))/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),sum((1-Z)*e*(1-e)*X)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),
          0,0,0,0,sum((1-Z)*e*(1-e)*X)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa)),sum((1-Z)*e*(1-e)*X^2)/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa))),
        nrow=6,ncol=6,byrow=T)
      B_comp <- matrix(c(sum(R*Z*(Y-mu_hat_1_eqn4)^2/(e^2)),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)/e),sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X/e),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)^2/(e^2)),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)/e),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X/e),
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)/e),0,sum(Z*(R-e)^2),sum(Z*(R-e)^2*X),0,0,
                         sum(R*Z*(Y-mu_hat_1_eqn4)*(R-e)*X/e),0,sum(Z*(R-e)^2*X),sum(Z*(R-e)^2*X^2),0,0,
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)/e),0,0,sum((1-Z)*(R-e)^2),sum((1-Z)*(R-e)^2*X),
                         0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(R-e)*X/e),0,0,sum((1-Z)*(R-e)^2*X),sum((1-Z)*(R-e)^2*X^2)),
                       nrow=6,ncol=6,byrow=T)
      B_hat <- B_comp/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa))
      var_comp <- solve(A_hat) %*% B_hat %*% t(solve(A_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/(ceiling(n_approx*(1-kappa))+ceiling(n_approx*kappa))
    }
    if (!is.null(delta)) {
      A_C_hat <- matrix(
        c(sum(R*Z/e)/max(K),0,sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)/e)/max(K),sum(R*Z*(Y-mu_hat_1_eqn4)*(1-e)*X/e)/max(K),0,0,
          0,sum(R*(1-Z)/e)/max(K),0,0,sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)/e)/max(K),sum(R*(1-Z)*(Y-mu_hat_0_eqn4)*(1-e)*X/e)/max(K),
          0,0,sum(Z*e*(1-e))/max(K),sum(Z*e*(1-e)*X)/max(K),0,0,
          0,0,sum(Z*e*(1-e)*X)/max(K),sum(Z*e*(1-e)*X^2)/max(K),0,0,
          0,0,0,0,sum((1-Z)*e*(1-e))/max(K),sum((1-Z)*e*(1-e)*X)/max(K),
          0,0,0,0,sum((1-Z)*e*(1-e)*X)/max(K),sum((1-Z)*e*(1-e)*X^2)/max(K)),
        nrow=6,ncol=6,byrow=T)
      B_C_comp <- matrix(0,ncol=6,nrow=6)
      for (k in 1:max(K)) {
        u_k <- c(sum(R[K==k]*Z[K==k]*(Y[K==k]-mu_hat_1_eqn4)/e[K==k]),
                 sum(R[K==k]*(1-Z[K==k])*(Y[K==k]-mu_hat_0_eqn4)/e[K==k]),
                 sum(Z[K==k]*(R[K==k]-e[K==k])),
                 sum(Z[K==k]*X[K==k]*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*(R[K==k]-e[K==k])),
                 sum((1-Z[K==k])*X[K==k]*(R[K==k]-e[K==k])))
        outer <- u_k %*% t(u_k)
        B_C_comp <- B_C_comp + outer
      }
      B_C_hat <- B_C_comp/max(K)
      var_comp <- solve(A_C_hat) %*% B_C_hat %*% t(solve(A_C_hat))
      var_eqn4 <- (var_comp[1,1] + var_comp[2,2])/max(K)
    }
    wald_eqn4 <- est_approx_eqn4/sqrt(var_eqn4)
    reject_approx_eqn4 <- wald_eqn4>qnorm(0.975)
    return(c(est_eqn3,est_standard_eqn4,est_IPRW_eqn4,est_known_eqn4,est_approx_eqn4,
             reject_eqn3,reject_standard_eqn4,reject_IPRW_eqn4,reject_known_eqn4,reject_approx_eqn4))
  }
  results <- colMeans(para)
  names(results) <- c("mu_hat_1-mu_hat_0 (eqn 3)","mu_hat_1-mu_hat_0 (standard eqn 4)","mu_hat_1-mu_hat_0 (IPRW eqn 4)",
                      "mu_hat_1-mu_hat_0 (kwown eqn 4)","mu_hat_1-mu_hat_0 (approx eqn 4)",
                      "empirical power (eqn 3)","empirical power (standard eqn 4)","empirical power (IPRW eqn 4)",
                      "empirical power (known eqn 4)","empirical power (approx eqn 4)")
  return(results)
}

