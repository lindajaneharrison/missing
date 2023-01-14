
# R libraries you may need
library(haven)
library(tidyverse)
library(latex2exp)
library(ggpubr)

# Read in pilot CRT data
zest <- read_stata("ZEST.dta")

# Limit to standard and delivery arm
zest_2arm <- zest[zest$study_arm %in% c("Delivery","Standard"),]

# Get Y, Z, R and X
zest_2arm$Y <- zest_2arm$hivtest_1mo_1M
zest_2arm$Z <- zest_2arm$arm
zest_2arm$R <- 1
zest_2arm$R[is.na(zest_2arm$Y)] <- 0
zest_2arm$X <- zest_2arm$B5B01_1
zest_2arm$X1 <- ifelse(zest_2arm$B5B01_1==0,1,0)
zest_2arm$X2 <- ifelse(zest_2arm$B5B01_1==1,1,0)

# Make more missing data if healthcare at a community clinic
set.seed(20220830)
zest_2arm$rand_number <- runif(dim(zest_2arm)[1])
zest_2arm$R[zest_2arm$rand_number<0.36 & zest_2arm$Z==1 & zest_2arm$X2==1] <- 0
zest_2arm$R[zest_2arm$rand_number<0.33 & zest_2arm$Z==0 & zest_2arm$X2==1] <- 0
zest_2arm$Y[zest_2arm$R==0] <- NA

# Healthcare at a community clinic (missingness model)
miss_model1 <- glm(R[arm==1]~0+X1[arm==1]+X2[arm==1],data=zest_2arm,family=binomial)
miss_model0 <- glm(R[arm==0]~0+X1[arm==0]+X2[arm==0],data=zest_2arm,family=binomial)
zest_2arm$e <- plogis(miss_model1$coefficients[1])
zest_2arm$e[zest_2arm$B5B01_1==1 & zest_2arm$Z==1] <- plogis(miss_model1$coefficients[2])
zest_2arm$e[zest_2arm$Z==0] <- plogis(miss_model0$coefficients[1])
zest_2arm$e[zest_2arm$B5B01_1==1 & zest_2arm$Z==0] <- plogis(miss_model0$coefficients[2])
zest_2arm$w <- 1/zest_2arm$e
zest_2arm$CLUSTER <- zest_2arm$peer

# Keep relevant data
keep <- as.data.frame(cbind(zest_2arm$Y,zest_2arm$Z, zest_2arm$R, 
                            zest_2arm$X,zest_2arm$CLUSTER,zest_2arm$e,zest_2arm$w,
                            zest_2arm$rand_number))
colnames(keep) <- c("Y","Z","R","X","CLUSTER","e","w","rand")

# New ID for clusters (peer educators)
first <- zest_2arm %>% group_by(peer) %>% filter(row_number()==1)
first$new_id <- seq(1,106,1)
first <- as.data.frame(cbind(first$peer,first$new_id))
colnames(first) <- c("peer","new_id")
keep <- merge(keep,first, by.x="CLUSTER", by.y="peer")

# Estimate mean and treatment difference (working independence)
mu_1_hat <- sum(keep$Z*keep$R*keep$Y/keep$e,na.rm=TRUE)/sum(keep$Z*keep$R/keep$e)
mu_0_hat <- sum((1-keep$Z)*keep$R*keep$Y/keep$e,na.rm=TRUE)/sum((1-keep$Z)*keep$R/(keep$e))  
lOR <- qlogis(mu_1_hat)-qlogis(mu_0_hat)

# Estimate weighted correlation parameter (GEE 1.5 framework)
keep$res <- (keep$Y - mu_1_hat)*keep$w
keep$res[keep$Z==0] <- (keep$Y[keep$Z==0] - mu_0_hat)*keep$w[keep$Z==0]
keep$var <- mu_1_hat*(1-mu_1_hat)
keep$var[keep$Z==0] <- mu_0_hat*(1-mu_0_hat)
keep$factor <- keep$res/sqrt(keep$var)
clust_sum <- rep(NA,max(keep$new_id))
clust_size <- rep(NA,max(keep$new_id))
for (i in 1:max(keep$new_id)){
  clust <- keep[keep$new_id==i,]
  clust <- clust[!is.na(clust$factor),]
  if (dim(clust)[1]>1)
    {clust_sum[i] <- sum(diag(clust$factor) %*% matrix(1,nrow=dim(clust)[1],ncol=dim(clust)[1]) %*% diag(clust$factor) -
                      diag(clust$factor) %*% diag(clust$factor))
     clust_size[i] <- sum(diag(clust$w) %*% matrix(1,nrow=dim(clust)[1],ncol=dim(clust)[1]) %*% diag(clust$w) -
                       diag(clust$w) %*% diag(clust$w))}
  if (dim(clust)[1]==1)
    {clust_sum[i] <- sum(clust$factor %*% matrix(1,nrow=dim(clust)[1],ncol=dim(clust)[1]) %*% clust$factor -
                      clust$factor %*% clust$factor)
     clust_size[i] <- sum(clust$w %*% matrix(1,nrow=dim(clust)[1],ncol=dim(clust)[1]) %*% clust$w -
                       clust$w %*% clust$w)}
}
delta <- sum(clust_sum)/sum(clust_size)

# Estimate all other parameters for tau_IPRW
pi_1 <- sum(keep$X)/length(keep$Z)
pi_2 <- sum((1-keep$X))/length(keep$Z)
mu_11 <- mean(keep$Y[keep$X & keep$Z],na.rm=TRUE)
mu_21 <- mean(keep$Y[(1-keep$X) & keep$Z],na.rm=TRUE)
mu_10 <- mean(keep$Y[keep$X & (1-keep$Z)],na.rm=TRUE)
mu_20 <- mean(keep$Y[(1-keep$X) & (1-keep$Z)],na.rm=TRUE)
expit_beta_11 <- mean(keep$R[keep$X & keep$Z])
expit_beta_21 <- mean(keep$R[(1-keep$X) & keep$Z])
expit_beta_10 <- mean(keep$R[keep$X & (1-keep$Z)])
expit_beta_20 <- mean(keep$R[(1-keep$X) & (1-keep$Z)])
c(pi_1,pi_2,expit_beta_11,expit_beta_21,expit_beta_10,expit_beta_20,
  mu_11,mu_21,mu_10,mu_20,mu_1_hat,mu_0_hat,delta,lOR)

# Proportion observed in pilot trial
kappa <- 1/2
pilot_prop <- kappa*(expit_beta_11*pi_1 + expit_beta_21*pi_2)+(1-kappa)*(expit_beta_10*pi_1 + expit_beta_20*pi_2)
pilot_prop

# Perform the sample size calculation for the comfirmatory cluster randomized trial
# IPRW approach
# with rounding to match R shiny app
round(c(pi_1,pi_2,expit_beta_11,expit_beta_21,expit_beta_10,expit_beta_20,
        mu_11,mu_21,mu_10,mu_20,mu_1_hat,mu_0_hat,delta,lOR),2)
tau_IPRW <- (0.67/(kappa*(0.9532*(1-0.9532))^2))*(0.94*(1-0.94)/0.61+(0.94-0.9532)^2) +
  (0.33/(kappa*(0.9532*(1-0.9532))^2))*(0.98*(1-0.98)/0.96+(0.98-0.9532)^2) +
  (0.67/((1-kappa)*(0.8797*(1-0.8797))^2))*(0.85*(1-0.85)/0.57+(0.85-0.8797)^2) +
  (0.33/((1-kappa)*(0.8797*(1-0.8797))^2))*(0.94*(1-0.94)/0.97+(0.94-0.8797)^2) +
  (6-1)*0.36*((1/(kappa*0.9532*(1-0.9532)))+(1/((1-kappa)*0.8797*(1-0.8797))))
tau_IPRW
sample_size <- tau_IPRW*(qnorm(0.9)+qnorm(0.975))^2/((qlogis(0.9532)-qlogis(0.8797))^2)
# Number of individuals
ceiling(sample_size)
# Number of clusters
2*ceiling((sample_size/2)/6)

# Standard approach
tau <- 1/(kappa*0.9532*(1-0.9532)) + 1/((1-kappa)*0.8797*(1-0.8797))
tau_standard <- ((1/0.71375) + (6-1)*0.36)*tau
tau_standard
sample_size_standard <- tau_standard*(qnorm(0.9)+qnorm(0.975))^2/((qlogis(0.9532)-qlogis(0.8797))^2)
# Number of individuals
ceiling(sample_size_standard)
# Number of clusters
2*ceiling((sample_size_standard/2)/6)







