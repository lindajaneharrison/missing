# Individually Randomized Trial (IRT) Case Study of the Relative Efficiency of the IPRW versus Standard Approach
# (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)

# R libraries needed
library(haven)
library(sjlabelled)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(fastGHQuad)

rx <- read_sas("NTISA5273LANHAV2016sent/upack/here/rx.sas7bdat",
               "NTISA5273LANHAV2016sent/upack/here/formats.sas7bcat")
cd4_analysis <- read_sas("NTISA5273LANHAV2016sent/upack/here/cd4_analysis.sas7bdat",
                         "NTISA5273LANHAV2016sent/upack/here/formats.sas7bcat")
lipid_analysis <- read_sas("NTISA5273LANHAV2016sent/upack/here/lipid_analysis.sas7bdat",
                           "NTISA5273LANHAV2016sent/upack/here/formats.sas7bcat")
baseline <- read_sas("NTISA5273LANHAV2016sent/upack/here/baseline.sas7bdat",
                     "NTISA5273LANHAV2016sent/upack/here/formats.sas7bcat")

# Randomized treatment
rx_R <- rx[,c(1,3)]
save(rx_R,file="rx_R.Rda")

# CD4 at 96 weeks
cd4_R <- cd4_analysis[cd4_analysis$calwk==96 & cd4_analysis$cd4_sched==1,c(5,19)]
cd4_R <- copy_labels(cd4_R,cd4_analysis)
save(cd4_R,file="cd4_R.Rda")

# Triglyceride at 96 weeks
lipid_R <- lipid_analysis[lipid_analysis$calwk==96 & !is.na(lipid_analysis$calwk) & 
                            !is.na(lipid_analysis$trigl) & lipid_analysis$lipid_sched==1,c(2,20)]
lipid_R <- copy_labels(lipid_R,lipid_analysis)
save(lipid_R,file="lipid_R.Rda")

# Baseline CD4
baseline_cd4_R <- baseline[,c(1,49)]
baseline_cd4_R <- copy_labels(baseline_cd4_R,baseline)
save(baseline_cd4_R,file="baseline_cd4_R.Rda")

# Baseline triglyceride
baseline_lipid_R <- lipid_analysis[,c(12,20)]
baseline_lipid_R <- baseline_lipid_R %>% group_by(ntisid) %>% filter(row_number()==1)
baseline_lipid_R <- copy_labels(baseline_lipid_R,lipid_analysis)
save(baseline_lipid_R,file="baseline_lipid_R.Rda")

# Participants with complete data for baseline CD4 and triglycerides
complete <- merge(baseline_cd4_R,rx_R, by.x="ntisid", by.y="ntisid", all.x=TRUE)
complete <- merge(complete,baseline_lipid_R,by.x="ntisid", by.y="ntisid", all=TRUE)
complete <- merge(complete,cd4_R,by.x="ntisid", by.y="ntisid", all=TRUE)
complete <- merge(complete,lipid_R,by.x="ntisid", by.y="ntisid", all=TRUE)
save(complete,file="complete.Rda")

# Half participants in each arm
kappa <- 0.5

# log(var)  (continuous outcome, weight binary variable)
# CD4 
complete$var <- log(complete$cd4)
complete$var_bl <- log(complete$bs_cd4)
cut_bl <- 100
cut <- 400

sigma_11_sq <- var(complete$var[complete$trt=="A" & complete$var_bl<log(cut_bl)],na.rm=T)
sigma_21_sq <- var(complete$var[complete$trt=="A" & complete$var_bl>=log(cut_bl)],na.rm=T)
sigma_10_sq <- var(complete$var[complete$trt=="B" & complete$var_bl<log(cut_bl)],na.rm=T)
sigma_20_sq <- var(complete$var[complete$trt=="B" & complete$var_bl>=log(cut_bl)],na.rm=T)

mu_11 <- mean(complete$var[complete$trt=="A" & complete$var_bl<log(cut_bl)],na.rm=T)
mu_21 <- mean(complete$var[complete$trt=="A" & complete$var_bl>=log(cut_bl)],na.rm=T)
mu_10 <- mean(complete$var[complete$trt=="B" & complete$var_bl<log(cut_bl)],na.rm=T)
mu_20 <- mean(complete$var[complete$trt=="B" & complete$var_bl>=log(cut_bl)],na.rm=T)  

mu_1 <- mean(complete$var[complete$trt=="A"],na.rm=T)
mu_0 <- mean(complete$var[complete$trt=="B"],na.rm=T)

pi_1 <- sum(complete$var_bl<log(cut_bl),na.rm=T)/sum(!is.na(complete$var_bl),na.rm=T)
pi_2 <- sum(complete$var_bl>=log(cut_bl),na.rm=T)/sum(!is.na(complete$var_bl),na.rm=T)


n_IPRW <- NULL
n_standard <- NULL
expit_beta_21 <- NULL
expit_beta_20 <- NULL
i <- 0
for (expit_beta_11 in seq(0.85,0.95,0.01)) {
  for (expit_beta_10 in seq(0.7,0.9,0.01)) {
    i <- i+1
    
    expit_beta_21[i] <- (0.9-pi_1*expit_beta_11)/pi_2
    expit_beta_20[i] <- (0.7-pi_1*expit_beta_10)/pi_2

    eta_IPRW <- (pi_1/kappa)*(sigma_11_sq/expit_beta_11 + (mu_11-mu_1)^2) +
                (pi_2/kappa)*(sigma_21_sq/expit_beta_21[i] + (mu_21-mu_1)^2) +
                (pi_1/(1-kappa))*(sigma_10_sq/expit_beta_10 + (mu_10-mu_0)^2) +
                (pi_2/(1-kappa))*(sigma_20_sq/expit_beta_20[i] + (mu_20-mu_0)^2) 

    n_IPRW[i] <- eta_IPRW*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 

    phi <- (pi_1*expit_beta_11 + 
            pi_2*expit_beta_21[i] + 
            pi_1*expit_beta_10 + 
            pi_2*expit_beta_20[i])/2 

    eta_standard <- (pi_1/kappa)*(sigma_11_sq + (mu_11-mu_1)^2)/phi +
                    (pi_2/kappa)*(sigma_21_sq + (mu_21-mu_1)^2)/phi +
                    (pi_1/(1-kappa))*(sigma_10_sq + (mu_10-mu_0)^2)/phi +
                    (pi_2/(1-kappa))*(sigma_20_sq + (mu_20-mu_0)^2)/phi 

    n_standard[i] <- eta_standard*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
  }
}

max(n_IPRW/n_standard)
min(n_IPRW/n_standard)
max(expit_beta_20)
min(expit_beta_20)
max(expit_beta_21)
min(expit_beta_21)

expit_beta_11 <- (0.9-pi_2*expit_beta_21)/pi_1
expit_beta_10 <- (0.7-pi_2*expit_beta_20)/pi_1

RE <- n_IPRW/n_standard

df <- as.data.frame(cbind(expit_beta_11,
                          expit_beta_21,
                          expit_beta_10,
                          expit_beta_20,
                          RE))
colnames(df) <- c("x","x1","y","y1","z")

a <- ggplot(df, aes(x, y, z=z)) + stat_contour(aes(colour = factor(stat(level))),linetype=2) 
data<- ggplot_build(a)$data[[1]] 
indices <- setdiff(1:nrow(data), which(duplicated(data$level))) # distinct levels
a1 <- a + 
  geom_text(aes(label=seq(0.985,1.02,by=0.005), z=NULL), size=2, data=data[indices,],nudge_x=-0.0001) +
  xlab(TeX("expit$(\\beta_{11})$")) +
  scale_x_continuous(breaks=c(0.85,0.875,0.9,0.925,0.95),
    sec.axis = sec_axis(~ (0.9-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{21})$"),breaks=c(0.85,0.9,0.95))) +
  ylab(TeX("expit$(\\beta_{10})$")) + 
  scale_y_continuous(breaks=c(0.7,0.75,0.8,0.85,0.9),
                     sec.axis = sec_axis(~ (0.7-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{20})$"))) +
  theme_classic() + theme(legend.position = "none", plot.title=element_text(size=11)) +
  ggtitle(TeX("C. $Y_i$ continuous, $g$ identity, $X_i$ categorical"))


# log(var)>log(cut)  (binary outcome, identity link, weight binary variable)
mu_11 <- mean(complete$var[complete$trt=="A" & complete$var_bl<log(cut_bl)]>log(cut),na.rm=T)
mu_21 <- mean(complete$var[complete$trt=="A" & complete$var_bl>=log(cut_bl)]>log(cut),na.rm=T)
mu_10 <- mean(complete$var[complete$trt=="B" & complete$var_bl<log(cut_bl)]>log(cut),na.rm=T)
mu_20 <- mean(complete$var[complete$trt=="B" & complete$var_bl>=log(cut_bl)]>log(cut),na.rm=T) 

sigma_11_sq <- mu_11*(1-mu_11)
sigma_21_sq <- mu_21*(1-mu_21)
sigma_10_sq <- mu_10*(1-mu_10)
sigma_20_sq <- mu_20*(1-mu_20)

mu_1 <- mean(complete$var[complete$trt=="A"]>log(cut),na.rm=T)
mu_0 <- mean(complete$var[complete$trt=="B"]>log(cut),na.rm=T)

n_IPRW <- NULL
n_standard <- NULL
i <- 0
for (expit_beta_11 in seq(0.85,0.95,0.01)) {
  for (expit_beta_10 in seq(0.7,0.9,0.01)) {
    i <- i+1
    
    expit_beta_21[i] <- (0.9-pi_1*expit_beta_11)/pi_2
    expit_beta_20[i] <- (0.7-pi_1*expit_beta_10)/pi_2
    
    eta_IPRW <- (pi_1/kappa)*(sigma_11_sq/expit_beta_11 + (mu_11-mu_1)^2) +
                (pi_2/kappa)*(sigma_21_sq/expit_beta_21[i] + (mu_21-mu_1)^2) +
                (pi_1/(1-kappa))*(sigma_10_sq/expit_beta_10 + (mu_10-mu_0)^2) +
                (pi_2/(1-kappa))*(sigma_20_sq/expit_beta_20[i] + (mu_20-mu_0)^2) 
    
    n_IPRW[i] <- eta_IPRW*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
    
    phi <- (pi_1*expit_beta_11 + 
            pi_2*expit_beta_21[i] + 
            pi_1*expit_beta_10 + 
            pi_2*expit_beta_20[i])/2 
    
    eta_standard <- (pi_1/kappa)*(sigma_11_sq + (mu_11-mu_1)^2)/phi +
                    (pi_2/kappa)*(sigma_21_sq + (mu_21-mu_1)^2)/phi +
                    (pi_1/(1-kappa))*(sigma_10_sq + (mu_10-mu_0)^2)/phi +
                    (pi_2/(1-kappa))*(sigma_20_sq + (mu_20-mu_0)^2)/phi 
    
    n_standard[i] <- eta_standard*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
  }
}

max(n_IPRW/n_standard)
min(n_IPRW/n_standard)
max(expit_beta_20)
min(expit_beta_20)
max(expit_beta_21)
min(expit_beta_21)

expit_beta_11 <- (0.9-pi_2*expit_beta_21)/pi_1
expit_beta_10 <- (0.7-pi_2*expit_beta_20)/pi_1

RE <- n_IPRW/n_standard

df <- as.data.frame(cbind(expit_beta_11,
                          expit_beta_21,
                          expit_beta_10,
                          expit_beta_20,
                          RE))
colnames(df) <- c("x","x1","y","y1","z")

b <- ggplot(df, aes(x, y, z=z)) + stat_contour(aes(colour = factor(stat(level))),linetype=2) 
data<- ggplot_build(b)$data[[1]] 
indices <- setdiff(1:nrow(data), which(duplicated(data$level))) # distinct levels
b1 <- b + 
  geom_text(aes(label=seq(0.976,0.996,by=0.002), z=NULL), size=2, data=data[indices,],nudge_x=0.0001) +
  xlab(TeX("expit$(\\beta_{11})$")) +
  scale_x_continuous(breaks=c(0.85,0.875,0.9,0.925,0.95),
                     sec.axis = sec_axis(~ (0.9-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{21})$"),breaks=c(0.85,0.9,0.95))) +
  ylab(TeX("expit$(\\beta_{10})$")) + 
  scale_y_continuous(breaks=c(0.7,0.75,0.8,0.85,0.9),
                     sec.axis = sec_axis(~ (0.7-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{20})$"))) +
  theme_classic() + theme(legend.position = "none", plot.title=element_text(size=11)) +
  ggtitle(TeX("A. $Y_i$ binary, $g$ identity, $X_i$ categorical"))


# log(var)>log(cut)  (binary outcome, logit link, weight binary variable)
n_IPRW <- NULL
n_standard <- NULL
i <- 0
for (expit_beta_11 in seq(0.85,0.95,0.01)) {
  for (expit_beta_10 in seq(0.7,0.9,0.01)) {
    i <- i+1
    
    expit_beta_21[i] <- (0.9-pi_1*expit_beta_11)/pi_2
    expit_beta_20[i] <- (0.7-pi_1*expit_beta_10)/pi_2
    
    eta_IPRW <- (pi_1/(kappa*(mu_11*(1-mu_11))^2))*(sigma_11_sq/expit_beta_11 + (mu_11-mu_1)^2) +
                (pi_2/(kappa*(mu_21*(1-mu_21))^2))*(sigma_21_sq/expit_beta_21[i] + (mu_21-mu_1)^2) +
                (pi_1/((1-kappa)*(mu_10*(1-mu_10))^2))*(sigma_10_sq/expit_beta_10 + (mu_10-mu_0)^2) +
                (pi_2/((1-kappa)*(mu_20*(1-mu_20))^2))*(sigma_20_sq/expit_beta_20[i] + (mu_20-mu_0)^2) 
    
    n_IPRW[i] <- eta_IPRW*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
    
    phi <- (pi_1*expit_beta_11 + 
            pi_2*expit_beta_21[i] + 
            pi_1*expit_beta_10 + 
            pi_2*expit_beta_20[i])/2 
    
    eta_standard <- (pi_1/(kappa*(mu_11*(1-mu_11))^2))*(sigma_11_sq + (mu_11-mu_1)^2)/phi +
                    (pi_2/(kappa*(mu_21*(1-mu_21))^2))*(sigma_21_sq + (mu_21-mu_1)^2)/phi +
                    (pi_1/((1-kappa)*(mu_10*(1-mu_10))^2))*(sigma_10_sq + (mu_10-mu_0)^2)/phi +
                    (pi_2/((1-kappa)*(mu_20*(1-mu_20))^2))*(sigma_20_sq + (mu_20-mu_0)^2)/phi 
    
    n_standard[i] <- eta_standard*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
  }
}

max(n_IPRW/n_standard)
min(n_IPRW/n_standard)
max(expit_beta_20)
min(expit_beta_20)
max(expit_beta_21)
min(expit_beta_21)

expit_beta_11 <- (0.9-pi_2*expit_beta_21)/pi_1
expit_beta_10 <- (0.7-pi_2*expit_beta_20)/pi_1

RE <- n_IPRW/n_standard

df <- as.data.frame(cbind(expit_beta_11,
                          expit_beta_21,
                          expit_beta_10,
                          expit_beta_20,
                          RE))
colnames(df) <- c("x","x1","y","y1","z")

c <- ggplot(df, aes(x, y, z=z)) + stat_contour(aes(colour = factor(stat(level))),linetype=2) 
data<- ggplot_build(c)$data[[1]] 
indices <- setdiff(1:nrow(data), which(duplicated(data$level))) # distinct levels
c1 <- c + 
  geom_text(aes(label=seq(0.985,1.01,by=0.005), z=NULL), size=2, data=data[indices,],nudge_x=-0.0001) +
  xlab(TeX("expit$(\\beta_{11})$")) +
  scale_x_continuous(breaks=c(0.85,0.875,0.9,0.925,0.95),
                     sec.axis = sec_axis(~ (0.9-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{21})$"),breaks=c(0.85,0.9,0.95))) +
  ylab(TeX("expit$(\\beta_{10})$")) + 
  scale_y_continuous(breaks=c(0.7,0.75,0.8,0.85,0.9),
                     sec.axis = sec_axis(~ (0.7-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{20})$"))) +
  theme_classic() + theme(legend.position = "none", plot.title=element_text(size=11)) +
  ggtitle(TeX("B. $Y_i$ binary, $g$ logit, $X_i$ categorical"))


# outcome cont., weight cont.
sigma_x_sq <- var(complete$var_bl,na.rm=T)
mu_x <- mean(complete$var_bl,na.rm=T)
sigma_y_sq <- var(complete$var,na.rm=T)
rho <- 0.5
to_opt_trt  <- function(beta_11,x,sigma_x_sq,mu_x) {
  x_j <- gaussHermiteData(100)$x 
  w_j <- gaussHermiteData(100)$w
  return(
    abs(sum(plogis(x+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_11)*w_j)/sqrt(pi)
        -0.9))
}
to_opt_control  <- function(beta_10,x,sigma_x_sq,mu_x) {
  x_j <- gaussHermiteData(100)$x 
  w_j <- gaussHermiteData(100)$w
  return(
    abs(sum(plogis(x+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_10)*w_j)/sqrt(pi)
        -0.7))
}
eta_IPRW <- NULL
beta_00_out <- NULL
beta_10_out <- NULL
beta_01_out <- NULL
beta_11_out <- NULL
x_j <- gaussHermiteData(100)$x 
w_j <- gaussHermiteData(100)$w
i <- 0
for (beta_10 in seq(-1.4,1.4,0.1)) {
  for (beta_11 in seq(-1.4,1.4,0.1)) {
    i <- i+1
    beta_01 <- optimize(interval=c(-10,10),f=to_opt_trt,
                    beta_11=beta_11,sigma_x_sq=sigma_x_sq,mu_x=mu_x)$minimum
    beta_00 <- optimize(interval=c(-10,10),f=to_opt_control,
                    beta_10=beta_10,sigma_x_sq=sigma_x_sq,mu_x=mu_x)$minimum
    C_1 <- c(0,sqrt(sigma_x_sq)/2) -
            (1/(sqrt(2*pi)))*c(sum(x_j*plogis(beta_01+beta_11*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j),
                     sum(x_j*(sqrt(2*sigma_x_sq)*x_j+mu_x)*plogis(beta_01+beta_11*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j))
    D_1 <- matrix(c(sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*w_j),
                sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)^2*w_j)),
              byrow=T,ncol=2,nrow=2)
    C_0 <- c(0,sqrt(sigma_x_sq)/2) -
            (1/(sqrt(2*pi)))*c(sum(x_j*plogis(beta_00+beta_10*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j),
                     sum(x_j*(sqrt(2*sigma_x_sq)*x_j+mu_x)*plogis(beta_00+beta_10*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j))
    D_0 <- matrix(c(sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*w_j),
                sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)^2*w_j)),
              byrow=T,ncol=2,nrow=2)
    eta_IPRW[i] <- sigma_y_sq*(1/kappa+exp(-mu_x*beta_11^2+beta_11^2*sigma_x_sq/2-beta_01)*(1+rho^2*sigma_x_sq*beta_11^2)/kappa+
                               1/(1-kappa)+exp(-mu_x*beta_10^2+beta_10^2*sigma_x_sq/2-beta_00)*(1+rho^2*sigma_x_sq*beta_10^2)/(1-kappa)-
                               rho^2*t(C_1)%*%solve(D_1)%*%C_1/kappa - rho^2*t(C_0)%*%solve(D_0)%*%C_0/(1-kappa)) 
    beta_00_out[i] <- beta_00
    beta_10_out[i] <- beta_10
    beta_01_out[i] <- beta_01
    beta_11_out[i] <- beta_11
  }
}
eta_standard <- sigma_y_sq/(0.8*kappa*(1-kappa))
RE <- eta_IPRW/eta_standard
RE

df <- as.data.frame(cbind(beta_11_out,
                          beta_01_out,
                          beta_10_out,
                          beta_00_out,
                          RE))
colnames(df) <- c("x","x1","y","y1","z")

d <- ggplot(df, aes(x, y, z=z)) + stat_contour(aes(colour = factor(stat(level))),linetype=2) 
data <- ggplot_build(d)$data[[1]] 
indices <- setdiff(1:nrow(data), which(duplicated(data$level))) # distinct levels
d1 <- d + 
  geom_text(aes(label=seq(0.8,1.6,by=0.1), z=NULL), size=2, data=data[indices,]) +
  xlab(TeX("$\\beta_{11}$")) +
  scale_x_continuous(sec.axis = sec_axis(~ exp(.),name=TeX("exp$(\\beta_{11})$"),
                                         breaks=c(0.25,0.5,1,2,4))) +
  ylab(TeX("$\\beta_{10}$")) + 
  scale_y_continuous(sec.axis = sec_axis(~ exp(.),name=TeX("exp$(\\beta_{10})$"),
                                         breaks=c(0.25,0.5,1,2,4))) +
  theme_classic() + theme(legend.position = "none", plot.title=element_text(size=11)) +
  ggtitle(TeX("D. $Y_i$ continuous, $g$ identity, $X_i$ continuous"))

out <- ggarrange(b1,c1,a1,d1,
          nrow=2,ncol=2)
annotate_figure(out,top = text_grob(TeX("CD4, Relative Efficiency: $\\tau_{IPRW}/\\tau_{standard}$"), face = "bold", size = 12))
ggsave("IRT_cd4.pdf", width = 7, height = 7)


# Triglyceride
complete$var <- log(complete$trigl)
complete$var_bl <- log(complete$trigl_bl)
cut_bl <- 150
cut <- 200

sigma_11_sq <- var(complete$var[complete$trt=="A" & complete$var_bl<log(cut_bl)],na.rm=T)
sigma_21_sq <- var(complete$var[complete$trt=="A" & complete$var_bl>=log(cut_bl)],na.rm=T)
sigma_10_sq <- var(complete$var[complete$trt=="B" & complete$var_bl<log(cut_bl)],na.rm=T)
sigma_20_sq <- var(complete$var[complete$trt=="B" & complete$var_bl>=log(cut_bl)],na.rm=T)

mu_11 <- mean(complete$var[complete$trt=="A" & complete$var_bl<log(cut_bl)],na.rm=T)
mu_21 <- mean(complete$var[complete$trt=="A" & complete$var_bl>=log(cut_bl)],na.rm=T)
mu_10 <- mean(complete$var[complete$trt=="B" & complete$var_bl<log(cut_bl)],na.rm=T)
mu_20 <- mean(complete$var[complete$trt=="B" & complete$var_bl>=log(cut_bl)],na.rm=T)  

mu_1 <- mean(complete$var[complete$trt=="A"],na.rm=T)
mu_0 <- mean(complete$var[complete$trt=="B"],na.rm=T)

pi_1 <- sum(complete$var_bl<log(cut_bl),na.rm=T)/sum(!is.na(complete$var_bl),na.rm=T)
pi_2 <- sum(complete$var_bl>=log(cut_bl),na.rm=T)/sum(!is.na(complete$var_bl),na.rm=T)


n_IPRW <- NULL
n_standard <- NULL
expit_beta_21 <- NULL
expit_beta_20 <- NULL
i <- 0
for (expit_beta_11 in seq(0.85,0.95,0.01)) {
  for (expit_beta_10 in seq(0.7,0.9,0.01)) {
    i <- i+1
    
    expit_beta_21[i] <- (0.9-pi_1*expit_beta_11)/pi_2
    expit_beta_20[i] <- (0.7-pi_1*expit_beta_10)/pi_2
    
    eta_IPRW <- (pi_1/kappa)*(sigma_11_sq/expit_beta_11 + (mu_11-mu_1)^2) +
                (pi_2/kappa)*(sigma_21_sq/expit_beta_21[i] + (mu_21-mu_1)^2) +
                (pi_1/(1-kappa))*(sigma_10_sq/expit_beta_10 + (mu_10-mu_0)^2) +
                (pi_2/(1-kappa))*(sigma_20_sq/expit_beta_20[i] + (mu_20-mu_0)^2) 
    
    n_IPRW[i] <- eta_IPRW*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
    
    phi <- (pi_1*expit_beta_11 + 
            pi_2*expit_beta_21[i] + 
            pi_1*expit_beta_10 + 
            pi_2*expit_beta_20[i])/2 
    
    eta_standard <- (pi_1/kappa)*(sigma_11_sq + (mu_11-mu_1)^2)/phi +
                    (pi_2/kappa)*(sigma_21_sq + (mu_21-mu_1)^2)/phi +
                    (pi_1/(1-kappa))*(sigma_10_sq + (mu_10-mu_0)^2)/phi +
                    (pi_2/(1-kappa))*(sigma_20_sq + (mu_20-mu_0)^2)/phi 
    
    n_standard[i] <- eta_standard*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
  }
}

max(n_IPRW/n_standard)
min(n_IPRW/n_standard)
max(expit_beta_20)
min(expit_beta_20)
max(expit_beta_21)
min(expit_beta_21)

expit_beta_11 <- (0.9-pi_2*expit_beta_21)/pi_1
expit_beta_10 <- (0.7-pi_2*expit_beta_20)/pi_1

RE <- n_IPRW/n_standard

df <- as.data.frame(cbind(expit_beta_11,
                          expit_beta_21,
                          expit_beta_10,
                          expit_beta_20,
                          RE))
colnames(df) <- c("x","x1","y","y1","z")

a <- ggplot(df, aes(x, y, z=z)) + stat_contour(aes(colour = factor(stat(level))),linetype=2) 
data<- ggplot_build(a)$data[[1]] 
indices <- setdiff(1:nrow(data), which(duplicated(data$level))) # distinct levels
a1 <- a + 
  geom_text(aes(label=seq(0.96,1.08,by=0.02), z=NULL), size=2, data=data[indices,],nudge_x=-0.0001) +
  xlab(TeX("expit$(\\beta_{11})$")) +
  scale_x_continuous(breaks=c(0.85,0.875,0.9,0.925,0.95),
                     sec.axis = sec_axis(~ (0.9-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{21})$"))) +
  ylab(TeX("expit$(\\beta_{10})$")) + 
  scale_y_continuous(breaks=c(0.7,0.75,0.8,0.85,0.9),
                     sec.axis = sec_axis(~ (0.7-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{20})$"))) +
  theme_classic() + theme(legend.position = "none", plot.title=element_text(size=11)) +
  ggtitle(TeX("C. $Y_i$ continuous, $g$ identity, $X_i$ categorical"))


# log(var)>log(cut)  (binary outcome, identity link, weight binary variable)
mu_11 <- mean(complete$var[complete$trt=="A" & complete$var_bl<log(cut_bl)]>log(cut),na.rm=T)
mu_21 <- mean(complete$var[complete$trt=="A" & complete$var_bl>=log(cut_bl)]>log(cut),na.rm=T)
mu_10 <- mean(complete$var[complete$trt=="B" & complete$var_bl<log(cut_bl)]>log(cut),na.rm=T)
mu_20 <- mean(complete$var[complete$trt=="B" & complete$var_bl>=log(cut_bl)]>log(cut),na.rm=T) 

sigma_11_sq <- mu_11*(1-mu_11)
sigma_21_sq <- mu_21*(1-mu_21)
sigma_10_sq <- mu_10*(1-mu_10)
sigma_20_sq <- mu_20*(1-mu_20)

mu_1 <- mean(complete$var[complete$trt=="A"]>log(cut),na.rm=T)
mu_0 <- mean(complete$var[complete$trt=="B"]>log(cut),na.rm=T)

n_IPRW <- NULL
n_standard <- NULL
i <- 0
for (expit_beta_11 in seq(0.85,0.95,0.01)) {
  for (expit_beta_10 in seq(0.7,0.9,0.01)) {
    i <- i+1
    
    expit_beta_21[i] <- (0.9-pi_1*expit_beta_11)/pi_2
    expit_beta_20[i] <- (0.7-pi_1*expit_beta_10)/pi_2
    
    eta_IPRW <- (pi_1/kappa)*(sigma_11_sq/expit_beta_11 + (mu_11-mu_1)^2) +
                (pi_2/kappa)*(sigma_21_sq/expit_beta_21[i] + (mu_21-mu_1)^2) +
                (pi_1/(1-kappa))*(sigma_10_sq/expit_beta_10 + (mu_10-mu_0)^2) +
                (pi_2/(1-kappa))*(sigma_20_sq/expit_beta_20[i] + (mu_20-mu_0)^2) 
    
    n_IPRW[i] <- eta_IPRW*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
    
    phi <- (pi_1*expit_beta_11 + 
            pi_2*expit_beta_21[i] + 
            pi_1*expit_beta_10 + 
            pi_2*expit_beta_20[i])/2 
    
    eta_standard <- (pi_1/kappa)*(sigma_11_sq + (mu_11-mu_1)^2)/phi +
                    (pi_2/kappa)*(sigma_21_sq + (mu_21-mu_1)^2)/phi +
                    (pi_1/(1-kappa))*(sigma_10_sq + (mu_10-mu_0)^2)/phi +
                    (pi_2/(1-kappa))*(sigma_20_sq + (mu_20-mu_0)^2)/phi 
    
    n_standard[i] <- eta_standard*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
  }
}

max(n_IPRW/n_standard)
min(n_IPRW/n_standard)
max(expit_beta_20)
min(expit_beta_20)
max(expit_beta_21)
min(expit_beta_21)

expit_beta_11 <- (0.9-pi_2*expit_beta_21)/pi_1
expit_beta_10 <- (0.7-pi_2*expit_beta_20)/pi_1

RE <- n_IPRW/n_standard

df <- as.data.frame(cbind(expit_beta_11,
                          expit_beta_21,
                          expit_beta_10,
                          expit_beta_20,
                          RE))
colnames(df) <- c("x","x1","y","y1","z")

b <- ggplot(df, aes(x, y, z=z)) + stat_contour(aes(colour = factor(stat(level))),linetype=2) 
data<- ggplot_build(b)$data[[1]] 
indices <- setdiff(1:nrow(data), which(duplicated(data$level))) # distinct levels
b1 <- b + 
  geom_text(aes(label=seq(0.98,1.14,by=0.02), z=NULL), size=2, data=data[indices,],nudge_x=0.0001) +
  xlab(TeX("expit$(\\beta_{11})$")) +
  scale_x_continuous(breaks=c(0.85,0.875,0.9,0.925,0.95),
                     sec.axis = sec_axis(~ (0.9-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{21})$"))) +
  ylab(TeX("expit$(\\beta_{10})$")) + 
  scale_y_continuous(breaks=c(0.7,0.75,0.8,0.85,0.9),
                     sec.axis = sec_axis(~ (0.7-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{20})$"))) +
  theme_classic() + theme(legend.position = "none", plot.title=element_text(size=11)) +
  ggtitle(TeX("C. $Y_i$ binary, $g$ identity, $X_i$ categorical"))

# log(var)>log(cut)  (binary outcome, logit link, weight binary variable)
n_IPRW <- NULL
n_standard <- NULL
i <- 0
for (expit_beta_11 in seq(0.85,0.95,0.01)) {
  for (expit_beta_10 in seq(0.7,0.9,0.01)) {
    i <- i+1
    
    expit_beta_21[i] <- (0.9-pi_1*expit_beta_11)/pi_2
    expit_beta_20[i] <- (0.7-pi_1*expit_beta_10)/pi_2
    
    eta_IPRW <- (pi_1/(kappa*(mu_11*(1-mu_11))^2))*(sigma_11_sq/expit_beta_11 + (mu_11-mu_1)^2) +
                (pi_2/(kappa*(mu_21*(1-mu_21))^2))*(sigma_21_sq/expit_beta_21[i] + (mu_21-mu_1)^2) +
                (pi_1/((1-kappa)*(mu_10*(1-mu_10))^2))*(sigma_10_sq/expit_beta_10 + (mu_10-mu_0)^2) +
                (pi_2/((1-kappa)*(mu_20*(1-mu_20))^2))*(sigma_20_sq/expit_beta_20[i] + (mu_20-mu_0)^2) 
    
    n_IPRW[i] <- eta_IPRW*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
    
    phi <- (pi_1*expit_beta_11 + 
            pi_2*expit_beta_21[i] + 
            pi_1*expit_beta_10 + 
            pi_2*expit_beta_20[i])/2 
    
    eta_standard <- (pi_1/(kappa*(mu_11*(1-mu_11))^2))*(sigma_11_sq + (mu_11-mu_1)^2)/phi +
                    (pi_2/(kappa*(mu_21*(1-mu_21))^2))*(sigma_21_sq + (mu_21-mu_1)^2)/phi +
                    (pi_1/((1-kappa)*(mu_10*(1-mu_10))^2))*(sigma_10_sq + (mu_10-mu_0)^2)/phi +
                    (pi_2/((1-kappa)*(mu_20*(1-mu_20))^2))*(sigma_20_sq + (mu_20-mu_0)^2)/phi 
    
    n_standard[i] <- eta_standard*(qnorm(0.9)+qnorm(0.975))^2/(10^2) 
  }
}

max(n_IPRW/n_standard)
min(n_IPRW/n_standard)
max(expit_beta_20)
min(expit_beta_20)
max(expit_beta_21)
min(expit_beta_21)

expit_beta_11 <- (0.9-pi_2*expit_beta_21)/pi_1
expit_beta_10 <- (0.7-pi_2*expit_beta_20)/pi_1

RE <- n_IPRW/n_standard

df <- as.data.frame(cbind(expit_beta_11,
                          expit_beta_21,
                          expit_beta_10,
                          expit_beta_20,
                          RE))
colnames(df) <- c("x","x1","y","y1","z")

c <- ggplot(df, aes(x, y, z=z)) + stat_contour(aes(colour = factor(stat(level))),linetype=2) 
data<- ggplot_build(c)$data[[1]] 
indices <- setdiff(1:nrow(data), which(duplicated(data$level))) # distinct levels
c1 <- c + 
  geom_text(aes(label=seq(0.975,1.03,by=0.005), z=NULL), size=2, data=data[indices,],nudge_x=-0.0001) +
  xlab(TeX("expit$(\\beta_{11})$")) +
  scale_x_continuous(breaks=c(0.85,0.875,0.9,0.925,0.95),
                     sec.axis = sec_axis(~ (0.9-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{21})$"))) +
  ylab(TeX("expit$(\\beta_{10})$")) + 
  scale_y_continuous(breaks=c(0.7,0.75,0.8,0.85,0.9),
                     sec.axis = sec_axis(~ (0.7-pi_2*.)/pi_1,name=TeX("expit$(\\beta_{20})$"))) +
  theme_classic() + theme(legend.position = "none", plot.title=element_text(size=11)) +
  ggtitle(TeX("B. $Y_i$ binary, $g$ logit, $X_i$ categorical"))


# outcome cont., weight cont.
sigma_x_sq <- var(complete$var_bl,na.rm=T)
mu_x <- mean(complete$var_bl,na.rm=T)
sigma_y_sq <- var(complete$var,na.rm=T)
rho <- 0.5
to_opt_trt  <- function(beta_11,x,sigma_x_sq,mu_x) {
  x_j <- gaussHermiteData(100)$x 
  w_j <- gaussHermiteData(100)$w
  return(
    abs(sum(plogis(x+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_11)*w_j)/sqrt(pi)
        -0.9))
}
to_opt_control  <- function(beta_10,x,sigma_x_sq,mu_x) {
  x_j <- gaussHermiteData(100)$x 
  w_j <- gaussHermiteData(100)$w
  return(
    abs(sum(plogis(x+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_10)*w_j)/sqrt(pi)
        -0.7))
}
eta_IPRW <- NULL
beta_00_out <- NULL
beta_10_out <- NULL
beta_01_out <- NULL
beta_11_out <- NULL
x_j <- gaussHermiteData(100)$x 
w_j <- gaussHermiteData(100)$w
i <- 0
for (beta_10 in seq(-1.4,1.4,0.1)) {
  for (beta_11 in seq(-1.4,1.4,0.1)) {
    i <- i+1
    beta_01 <- optimize(interval=c(-10,10),f=to_opt_trt,
                        beta_11=beta_11,sigma_x_sq=sigma_x_sq,mu_x=mu_x)$minimum
    beta_00 <- optimize(interval=c(-10,10),f=to_opt_control,
                        beta_10=beta_10,sigma_x_sq=sigma_x_sq,mu_x=mu_x)$minimum
    C_1 <- c(0,sqrt(sigma_x_sq)/2) -
      (1/(sqrt(2*pi)))*c(sum(x_j*plogis(beta_01+beta_11*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j),
                         sum(x_j*(sqrt(2*sigma_x_sq)*x_j+mu_x)*plogis(beta_01+beta_11*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j))
    D_1 <- matrix(c(sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*w_j),
                    sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                    sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                    sum(plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_01+beta_11*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)^2*w_j)),
                  byrow=T,ncol=2,nrow=2)
    C_0 <- c(0,sqrt(sigma_x_sq)/2) -
      (1/(sqrt(2*pi)))*c(sum(x_j*plogis(beta_00+beta_10*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j),
                         sum(x_j*(sqrt(2*sigma_x_sq)*x_j+mu_x)*plogis(beta_00+beta_10*(x_j*sqrt(2*sigma_x_sq)+mu_x))*w_j))
    D_0 <- matrix(c(sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*w_j),
                    sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                    sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)*w_j),
                    sum(plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x))*(1-plogis(beta_00+beta_10*(sqrt(2*sigma_x_sq)*x_j+mu_x)))*(sqrt(2*sigma_x_sq)*x_j+mu_x)^2*w_j)),
                  byrow=T,ncol=2,nrow=2)
    eta_IPRW[i] <- sigma_y_sq*(1/kappa+exp(-mu_x*beta_11^2+beta_11^2*sigma_x_sq/2-beta_01)*(1+rho^2*sigma_x_sq*beta_11^2)/kappa+
                               1/(1-kappa)+exp(-mu_x*beta_10^2+beta_10^2*sigma_x_sq/2-beta_00)*(1+rho^2*sigma_x_sq*beta_10^2)/(1-kappa)-
                               rho^2*t(C_1)%*%solve(D_1)%*%C_1/kappa - rho^2*t(C_0)%*%solve(D_0)%*%C_0/(1-kappa)) 
    beta_00_out[i] <- beta_00
    beta_10_out[i] <- beta_10
    beta_01_out[i] <- beta_01
    beta_11_out[i] <- beta_11
  }
}
eta_standard <- sigma_y_sq/(0.8*kappa*(1-kappa))
RE <- eta_IPRW/eta_standard
RE

df <- as.data.frame(cbind(beta_11_out,
                          beta_01_out,
                          beta_10_out,
                          beta_00_out,
                          RE))
colnames(df) <- c("x","x1","y","y1","z")

d <- ggplot(df, aes(x, y, z=z)) + stat_contour(aes(colour = factor(stat(level))),linetype=2) 
data<- ggplot_build(d)$data[[1]] 
indices <- setdiff(1:nrow(data), which(duplicated(data$level))) # distinct levels
d1 <- d + 
  geom_text(aes(label=seq(0.8,1.5,by=0.1), z=NULL), size=2, data=data[indices,]) +
  xlab(TeX("$\\beta_{11}$")) +
  scale_x_continuous(sec.axis = sec_axis(~ exp(.),name=TeX("exp$(\\beta_{11})$"),
                                         breaks=c(0.25,0.5,1,2,4))) +
  ylab(TeX("$\\beta_{10}$")) + 
  scale_y_continuous(sec.axis = sec_axis(~ exp(.),name=TeX("exp$(\\beta_{10})$"),
                                         breaks=c(0.25,0.5,1,2,4))) +
  theme_classic() + theme(legend.position = "none", plot.title=element_text(size=11)) +
  ggtitle(TeX("D. $Y_i$ continuous, $g$ identity, $X_i$ continuous"))

out <- ggarrange(b1,c1,a1,d1,
                 nrow=2,ncol=2)
annotate_figure(out,top = text_grob(TeX("Triglyceride, Relative Efficiency: $\\tau_{IPRW}/\\tau_{standard}$"), face = "bold", size = 12))
ggsave("IRT_triglyceride.pdf", width = 7, height = 7)



