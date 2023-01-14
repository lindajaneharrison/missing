# Figure 4 (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
# Pilot simulation
# change size of pilot study
source('useFUNCTIONS.R')

### Individually Randomized Trials (IRTs)

# Weighting for a fully observed auxiliary categorical variable

# Sample Size Calculation and Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(20220318)

# Other R packages
library(ggplot2)
library(ggbreak)
library(latex2exp)
library(aplot)
library(cowplot)
library(patchwork)


# Y continuous, g identity
# scenario 1
combine <- NULL
rep <- 10000
for (n_pilot in c(50,100,250,500)){
  out_cont_1 <- sim_power_cat_pilot(rep=rep,kappa=0.5,
                                    type="continuous",link="identity",
                                    mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,pi_11=0.5,
                                    sigma_11_sq=0.026,sigma_21_sq=0.294,sigma_10_sq=0.026,sigma_20_sq=0.229,
                                    expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
                                    n_pilot=n_pilot)
  results <- 100*tail(out_cont_1,n=1)
  results[1:9] <- tail(out_cont_1,n=1)[1:9]
  results[15] <- n_pilot
  results[16] <- median(out_cont_1[1:rep,1],na.rm=TRUE)
  results[17] <- median(out_cont_1[1:rep,2],na.rm=TRUE)
  results[18] <- median(out_cont_1[1:rep,3],na.rm=TRUE)
  results[19] <- median(out_cont_1[1:rep,4],na.rm=TRUE)
  combine <- rbind(combine,results)
}
colnames(combine) <- c(colnames(out_cont_1),"n_pilot","n_standard","n_IPRW","n_known","n_approx")
combine_cont_id_1 <- data.frame(combine)
# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC79A7")
nice.shapes <- c(16,17,15,3,4)
fig_cont_id_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..eqn.3.,color="iv. standard (complete case)",shape="iv. standard (complete case)")) + 
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(0,10,20)) +
  scale_y_break(c(103, 1000),scales=2.5,ticklabels=c(1000,2000,3000,4000,5000,6000)) + scale_y_break(c(65,25),ticklabels=c(70,80,90,100)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power      n",color="") + ggtitle("Scenario 1") + theme(plot.title = element_text(size=10))
fig_cont_id_1 
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)
fig_cont_id_n_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  scale_y_continuous(breaks=c(1000,1250,1500,1750),lim=c(900,1750)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y=TeX("$\\bar{n}$"),color="") + ggtitle(TeX("$Y_i$ continuous, $g$ identity")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=1750, label = "Scenario 1", size=3)
fig_cont_id_n_1 
fig_cont_id_power_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_cont_id_1,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(75,80,85,90,95,100), lim=c(74,105)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power",color="") + ggtitle(TeX("$Y_i$ continuous, $g$ identity")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=105, label = "Scenario 1", size=3)
fig_cont_id_power_1

# scenario 2
combine <- NULL
for (n_pilot in c(50,100,250,500)){
  out_cont_2 <- sim_power_cat_pilot(rep=rep,kappa=0.5,
                                    type="continuous",link="identity",
                                    mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,pi_11=0.5,
                                    sigma_11_sq=0.41955,sigma_21_sq=0.026,sigma_10_sq=0.46795,sigma_20_sq=0.026,
                                    expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
                                    n_pilot=n_pilot)
  results <- 100*tail(out_cont_2,n=1)
  results[1:9] <- tail(out_cont_2,n=1)[1:9]
  results[15] <- n_pilot
  results[16] <- median(out_cont_2[1:rep,1],na.rm=TRUE)
  results[17] <- median(out_cont_2[1:rep,2],na.rm=TRUE)
  results[18] <- median(out_cont_2[1:rep,3],na.rm=TRUE)
  results[19] <- median(out_cont_2[1:rep,4],na.rm=TRUE)
  combine <- rbind(combine,results)
}
colnames(combine) <- c(colnames(out_cont_2),"n_pilot","n_standard","n_IPRW","n_known","n_approx")
combine_cont_id_2 <- data.frame(combine)
# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC79A7")
nice.shapes <- c(16,17,15,3,4)
fig_cont_id_2 <- ggplot() + 
  theme_classic() +
  geom_hline(aes(yintercept=8500),color="white") +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..eqn.3.,color="iv. standard (complete case)",shape="iv. standard (complete case)")) + 
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  geom_hline(aes(yintercept=36),color="white") +
  scale_y_continuous(breaks=c(50,70,90)) +
  scale_y_break(c(99, 900),scales=1.5,ticklabels=c(2000,4000,6000,8000)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power      n",color="") + ggtitle("Scenario 2") + theme(plot.title = element_text(size=10))
fig_cont_id_2
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)
fig_cont_id_n_2 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  scale_y_continuous(breaks=c(1000,1250,1500,1750),lim=c(900,1750)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y=TeX("$\\bar{n}$"),color="") + ggtitle(TeX("$Y_i$ continuous, $g$ identity")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=1750, label = "Scenario 2", size=3)
fig_cont_id_n_2
fig_cont_id_power_2 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_cont_id_2,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(75,80,85,90,95,100), lim=c(74,105)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power",color="") + ggtitle(TeX("$Y_i$ continuous, $g$ identity")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=105, label = "Scenario 2", size=3)
fig_cont_id_power_2


# Y binary, g identity
# scenario 1
combine <- NULL
for (n_pilot in c(50,100,250,500)){
  out_bin_1 <- sim_power_cat_pilot(rep=rep,kappa=0.5,
                                    type="binary",link="identity",
                                    mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,pi_11=0.5,
                                    sigma_11_sq=0.09,sigma_21_sq=0.21,sigma_10_sq=0.1275,sigma_20_sq=0.1275,
                                    expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
                                    n_pilot=n_pilot)
  results <- 100*tail(out_bin_1,n=1)
  results[1:9] <- tail(out_bin_1,n=1)[1:9]
  results[15] <- n_pilot
  results[16] <- median(out_bin_1[1:rep,1],na.rm=TRUE)
  results[17] <- median(out_bin_1[1:rep,2],na.rm=TRUE)
  results[18] <- median(out_bin_1[1:rep,3],na.rm=TRUE)
  results[19] <- median(out_bin_1[1:rep,4],na.rm=TRUE)
  combine <- rbind(combine,results)
}
colnames(combine) <- c(colnames(out_bin_1),"n_pilot","n_standard","n_IPRW","n_known","n_approx")
combine_bin_id_1 <- data.frame(combine)
# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC79A7")
nice.shapes <- c(16,17,15,3,4)
fig_bin_id_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..eqn.3.,color="iv. standard (complete case)",shape="iv. standard (complete case)")) + 
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(0,10,20)) +
  scale_y_break(c(103, 1000),scales=2.5,ticklabels=c(1000,2000,3000,4000,5000,6000)) + scale_y_break(c(25,65),ticklabels=c(70,80,90,100)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power      n",color="")  + ggtitle("Scenario 3") + theme(plot.title = element_text(size=10))
fig_bin_id_1
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)
fig_bin_id_n_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  scale_y_continuous(breaks=c(1000,1250,1500,1750),lim=c(900,1750)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y=TeX("$\\bar{n}$"),color="") + ggtitle(TeX("$Y_i$ binary, $g$ identity")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=1750, label = "Scenario 1", size=3)
fig_bin_id_n_1 
fig_bin_id_power_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_id_1,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(75,80,85,90,95,100), lim=c(74,105)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power",color="") + ggtitle(TeX("$Y_i$ binary, $g$ identity")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=105, label = "Scenario 1", size=3)
fig_bin_id_power_1

# scenario 2
combine <- NULL
for (n_pilot in c(50,100,250,500)){
  out_bin_2 <- sim_power_cat_pilot(rep=rep,kappa=0.5,
                                   type="binary",link="identity",
                                   mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,pi_11=0.5,
                                   sigma_11_sq=0.2356,sigma_21_sq=0.0475,sigma_10_sq=0.2331,sigma_20_sq=0.1924,
                                   expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
                                   n_pilot=n_pilot)
  results <- 100*tail(out_bin_2,n=1)
  results[1:9] <- tail(out_bin_2,n=1)[1:9]
  results[15] <- n_pilot
  results[16] <- median(out_bin_2[1:rep,1],na.rm=TRUE)
  results[17] <- median(out_bin_2[1:rep,2],na.rm=TRUE)
  results[18] <- median(out_bin_2[1:rep,3],na.rm=TRUE)
  results[19] <- median(out_bin_2[1:rep,4],na.rm=TRUE)
  combine <- rbind(combine,results)
}
colnames(combine) <- c(colnames(out_bin_2),"n_pilot","n_standard","n_IPRW","n_known","n_approx")
combine_bin_id_2 <- data.frame(combine)
# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC79A7")
nice.shapes <- c(16,17,15,3,4)
fig_bin_id_2 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..eqn.3.,color="iv. standard (complete case)",shape="iv. standard (complete case)")) + 
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  geom_hline(aes(yintercept=36),color="white") +
  geom_hline(aes(yintercept=8500),color="white") +
  scale_y_continuous(breaks=c(50,70,90)) +
  scale_y_break(c(100, 900),scales=1.5,ticklabels=c(2000,4000,6000,8000)) + 
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_pilot$"),y="power      n",color="") + ggtitle("Scenario 4") + theme(plot.title = element_text(size=10))
fig_bin_id_2
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)
fig_bin_id_n_2 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  scale_y_continuous(breaks=c(1000,1250,1500,1750),lim=c(900,1750)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y=TeX("$\\bar{n}$"),color="") + ggtitle(TeX("$Y_i$ binary, $g$ identity")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=1750, label = "Scenario 2", size=3)
fig_bin_id_n_2 
fig_bin_id_power_2 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_id_2,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(75,80,85,90,95,100), lim=c(74,105)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power",color="") + ggtitle(TeX("$Y_i$ binary, $g$ identity")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=105, label = "Scenario 2", size=3)
fig_bin_id_power_2


# Y binary, g logit
# scenario 1
combine <- NULL
for (n_pilot in c(50,100,250,500)){
  out_bin_1 <- sim_power_cat_pilot(rep=rep,kappa=0.5,
                                   type="binary",link="logit",
                                   mu_11=0.9,mu_21=0.3,mu_10=0.15,mu_20=0.85,pi_11=0.5,
                                   sigma_11_sq=0.09,sigma_21_sq=0.21,sigma_10_sq=0.1275,sigma_20_sq=0.1275,
                                   expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
                                   n_pilot=n_pilot)
  results <- 100*tail(out_bin_1,n=1)
  results[1:9] <- tail(out_bin_1,n=1)[1:9]
  results[15] <- n_pilot
  results[16] <- median(out_bin_1[1:rep,1],na.rm=TRUE)
  results[17] <- median(out_bin_1[1:rep,2],na.rm=TRUE)
  results[18] <- median(out_bin_1[1:rep,3],na.rm=TRUE)
  results[19] <- median(out_bin_1[1:rep,4],na.rm=TRUE)
  combine <- rbind(combine,results)
}
colnames(combine) <- c(colnames(out_bin_1),"n_pilot","n_standard","n_IPRW","n_known","n_approx")
combine_bin_logit_1 <- data.frame(combine)
# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC79A7")
nice.shapes <- c(16,17,15,3,4)
fig_bin_logit_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..eqn.3.,color="iv. standard (complete case)",shape="iv. standard (complete case)")) + 
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(0,10,20)) +
  scale_y_break(c(103, 1000),scales=2.5,ticklabels=c(1000,2000,3000,4000,5000,6000)) + scale_y_break(c(25,65),ticklabels=c(70,80,90,100)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_pilot$"),y="power      n",color="") + ggtitle("Scenario 5") + theme(plot.title = element_text(size=10))
fig_bin_logit_1
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)
fig_bin_logit_n_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  scale_y_continuous(breaks=c(1000,1250,1500,1750),lim=c(900,1750)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y=TeX("$\\bar{n}$"),color="") + ggtitle(TeX("$Y_i$ binary, $g$ logit")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=1750, label = "Scenario 1", size=3)
fig_bin_logit_n_1 
fig_bin_logit_power_1 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_logit_1,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(75,80,85,90,95,100), lim=c(74,105)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power",color="") + ggtitle(TeX("$Y_i$ binary, $g$ logit")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=105, label = "Scenario 1", size=3)
fig_bin_logit_power_1

# scenario 2
combine <- NULL
for (n_pilot in c(50,100,250,500)){
  out_bin_2 <- sim_power_cat_pilot(rep=rep,kappa=0.5,
                                   type="binary",link="logit",
                                   mu_11=0.62,mu_21=0.95,mu_10=0.63,mu_20=0.74,pi_11=0.5,
                                   sigma_11_sq=0.2356,sigma_21_sq=0.0475,sigma_10_sq=0.2331,sigma_20_sq=0.1924,
                                   expit_beta_11=0.7,expit_beta_21=0.9,expit_beta_10=0.75,expit_beta_20=0.85,
                                   n_pilot=n_pilot)
  results <- 100*tail(out_bin_2,n=1)
  results[1:9] <- tail(out_bin_2,n=1)[1:9]
  results[15] <- n_pilot
  results[16] <- median(out_bin_2[1:rep,1],na.rm=TRUE)
  results[17] <- median(out_bin_2[1:rep,2],na.rm=TRUE)
  results[18] <- median(out_bin_2[1:rep,3],na.rm=TRUE)
  results[19] <- median(out_bin_2[1:rep,4],na.rm=TRUE)
  combine <- rbind(combine,results)
}
colnames(combine) <- c(colnames(out_bin_2),"n_pilot","n_standard","n_IPRW","n_known","n_approx")
combine_bin_logit_2 <- data.frame(combine)
# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC79A7")
nice.shapes <- c(16,17,15,3,4)
fig_bin_logit_2 <- ggplot() + 
  theme_classic() +
  geom_hline(aes(yintercept=8500),color="white") +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..eqn.3.,color="iv. standard (complete case)",shape="iv. standard (complete case)")) + 
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  geom_hline(aes(yintercept=36),color="white") +
  scale_y_continuous(breaks=c(50,70,90)) +
  scale_y_break(c(100, 900),scales=1.5,ticklabels=c(2000,4000,6000,8000)) + 
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_pilot$"),y="power      n",color="") + ggtitle("Scenario 6") + theme(plot.title = element_text(size=10))
fig_bin_logit_2
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)
fig_bin_logit_n_2 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  scale_y_continuous(breaks=c(1000,1250,1500,1750),lim=c(900,1750)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y=TeX("$\\bar{n}$"),color="") + ggtitle(TeX("$Y_i$ binary, $g$ logit")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=1750, label = "Scenario 2", size=3)
fig_bin_logit_n_2 
fig_bin_logit_power_2 <- ggplot() + 
  theme_classic() +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine_bin_logit_2,aes(x=n_pilot,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(75,80,85,90,95,100), lim=c(74,105)) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\n_{pilot}$"),y="power",color="") + ggtitle(TeX("$Y_i$ binary, $g$ logit")) + 
  theme(plot.title = element_text(size=10)) +
  annotate("text", x = 400, y=105, label = "Scenario 2", size=3)
fig_bin_logit_power_2

legend <- cowplot::get_legend(fig_cont_id_1 + theme(legend.position = "bottom"))
out_plot <- aplot::plot_list(fig_cont_id_1 + theme(legend.position = "none"),
                 fig_cont_id_2 + theme(legend.position = "none"),
                 fig_bin_id_1 + theme(legend.position = "none"), 
                 fig_bin_id_2 + theme(legend.position = "none"), 
                 fig_bin_logit_1 + theme(legend.position = "none"),
                 fig_bin_logit_2 + theme(legend.position = "none"),
                 ncol = 2, nrow = 3)
aplot::plot_list(out_plot,legend,ncol=1,heights=c(3,0.1))
ggsave("Fig_IRT_cat_pilot_size.pdf")

legend <- cowplot::get_legend(fig_cont_id_n_1 + theme(legend.position = "bottom"))
out_plot <- aplot::plot_list(fig_cont_id_n_1 + theme(legend.position = "none"),
                             fig_cont_id_n_2 + theme(legend.position = "none"),
                             fig_bin_id_n_1 + theme(legend.position = "none"), 
                             fig_bin_id_n_2 + theme(legend.position = "none"), 
                             fig_bin_logit_n_1 + theme(legend.position = "none"),
                             fig_bin_logit_n_2 + theme(legend.position = "none"),
                             ncol = 2, nrow = 3, labels="(a)")
aplot::plot_list(out_plot,legend,ncol=1,heights=c(3,0.1))
ggsave("Fig_IRT_cat_n_pilot_size.pdf",width=5.75,height=8)

legend <- cowplot::get_legend(fig_cont_id_power_1 + theme(legend.position = "bottom"))
out_plot <- aplot::plot_list(fig_cont_id_power_1 + theme(legend.position = "none"),
                             fig_cont_id_power_2 + theme(legend.position = "none"),
                             fig_bin_id_power_1 + theme(legend.position = "none"), 
                             fig_bin_id_power_2 + theme(legend.position = "none"), 
                             fig_bin_logit_power_1 + theme(legend.position = "none"),
                             fig_bin_logit_power_2 + theme(legend.position = "none"),
                             ncol = 2, nrow = 3, labels="(b)")
aplot::plot_list(out_plot,legend,ncol=1,heights=c(3,0.1))
ggsave("Fig_IRT_cat_power_pilot_size.pdf",width=5.75,height=8)



