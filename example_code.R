
library(haven)
library(tidyverse)

zest <- read_stata("ZEST.dta")

# Limit to standard and delivery arm
zest_2arm <- zest[zest$study_arm %in% c("Delivery","Standard"),]

# Healthcare at a community clinic (missingness model)
zest_2arm$Y <- zest_2arm$hivtest_1mo_1M
zest_2arm$Z <- zest_2arm$arm
zest_2arm$R <- 1
zest_2arm$R[is.na(zest_2arm$Y)] <- 0
zest_2arm$X <- zest_2arm$B5B01_1
zest_2arm$X1 <- ifelse(zest_2arm$B5B01_1==0,1,0)
zest_2arm$X2 <- ifelse(zest_2arm$B5B01_1==1,1,0)
miss_model1 <- glm(R[arm==1]~0+X1[arm==1]+X2[arm==1],data=zest_2arm,family=binomial)
miss_model0 <- glm(R[arm==0]~0+X1[arm==0]+X2[arm==0],data=zest_2arm,family=binomial)
zest_2arm$e <- plogis(miss_model1$coefficients[1])
zest_2arm$e[zest_2arm$B5B01_1==1 & zest_2arm$Z==1] <- plogis(miss_model1$coefficients[2])
zest_2arm$e[zest_2arm$Z==0] <- plogis(miss_model0$coefficients[1])
zest_2arm$e[zest_2arm$B5B01_1==1 & zest_2arm$Z==0] <- plogis(miss_model0$coefficients[2])
zest_2arm$w <- 1/zest_2arm$e
zest_2arm$CLUSTER <- zest_2arm$peer
keep <- as.data.frame(cbind(zest_2arm$Y,zest_2arm$Z, zest_2arm$R, 
                            zest_2arm$X,zest_2arm$CLUSTER,zest_2arm$e,zest_2arm$w))
colnames(keep) <- c("Y","Z","R","X","CLUSTER","e","w")

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
  clust_sum[i] <- sum(diag(clust$factor) %*% matrix(1,nrow=dim(clust)[1],ncol=dim(clust)[1]) %*% diag(clust$factor) -
                      diag(clust$factor) %*% diag(clust$factor))
  clust_size[i] <- sum(diag(clust$w) %*% matrix(1,nrow=dim(clust)[1],ncol=dim(clust)[1]) %*% diag(clust$w) -
                       diag(clust$w) %*% diag(clust$w))
}
delta <- sum(clust_sum)/sum(clust_size)

# Estimate all other parameters for eta_IPRW
pi_11 <- sum(keep$X & keep$Z)/sum(keep$Z)
pi_21 <- sum((1-keep$X) & keep$Z)/sum(keep$Z)
pi_10 <- sum(keep$X & (1-keep$Z))/sum(1-keep$Z)
pi_20 <- sum((1-keep$X) & (1-keep$Z))/sum(1-keep$Z)
mu_11 <- mean(keep$Y[keep$X & keep$Z],na.rm=TRUE)
mu_21 <- mean(keep$Y[(1-keep$X) & keep$Z],na.rm=TRUE)
mu_10 <- mean(keep$Y[keep$X & (1-keep$Z)],na.rm=TRUE)
mu_20 <- mean(keep$Y[(1-keep$X) & (1-keep$Z)],na.rm=TRUE)
expit_beta_11 <- mean(keep$R[keep$X & keep$Z])
expit_beta_21 <- mean(keep$R[(1-keep$X) & keep$Z])
expit_beta_10 <- mean(keep$R[keep$X & (1-keep$Z)])
expit_beta_20 <- mean(keep$R[(1-keep$X) & (1-keep$Z)])
c(pi_11,pi_21,pi_10,pi_20,expit_beta_11,expit_beta_21,expit_beta_10,expit_beta_20,
  mu_11,mu_21,mu_10,mu_20,mu_1_hat,mu_0_hat,delta,lOR)

# Perform the sample size calculation for the comfirmatory cluster randomized trial
kappa <- 1/2
eta_IPRW <- (pi_11/(kappa*(mu_1_hat*(1-mu_1_hat))^2))*(mu_11*(1-mu_11)/expit_beta_11+(mu_11-mu_1_hat)^2) +
            (pi_21/(kappa*(mu_1_hat*(1-mu_1_hat))^2))*(mu_21*(1-mu_21)/expit_beta_21+(mu_21-mu_1_hat)^2) +
            (pi_10/((1-kappa)*(mu_0_hat*(1-mu_0_hat))^2))*(mu_10*(1-mu_10)/expit_beta_10+(mu_10-mu_0_hat)^2) +
            (pi_20/((1-kappa)*(mu_0_hat*(1-mu_0_hat))^2))*(mu_20*(1-mu_20)/expit_beta_20+(mu_20-mu_0_hat)^2) +
            (6-1)*delta*((1/(kappa*mu_1_hat*(1-mu_1_hat)))+(1/((1-kappa)*mu_0_hat*(1-mu_0_hat))))
eta_IPRW
sample_size <- eta_IPRW*(qnorm(0.9)+qnorm(0.975))^2/((qlogis(mu_1_hat)-qlogis(mu_0_hat))^2)
# Number of individuals
sample_size
# Number of clusters
sample_size/6












