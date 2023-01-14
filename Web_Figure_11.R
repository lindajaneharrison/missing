# Web Figure 11 (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

### Cluster Randomized Trials (CRTs)

# Weighting for a fully observed auxiliary continuous variable

# Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(4376758)

# Other R packages
library(ggplot2)
library(ggbreak)
library(latex2exp)
library(aplot)
library(cowplot)
library(patchwork)
library(fastGHQuad)

# (X,Y) bivariate normal, g identity
combine <- NULL
rep <- 10000
kappa <- 0.5
for (delta in c(0.03,0.05,0.07)){
  print(delta)
  n_standard_out <- n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                               type="continuous",link="identity",sigma_y_sq=0.245,phi=0.8,
                               delta=delta,m=5)
  n_IPRW <- n_IPRW_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                        mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
                        beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64,
                        delta=delta,m=5)
  n_known <- n_known_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                          mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
                          beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64,
                          delta=delta,m=5)
  n_approx <- n_approx_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                            mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,
                            beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64,
                            delta=delta,m=5)
  sim_out <- sim_power_cont(rep=rep,kappa=0.5,n_standard=n_standard_out,n_IPRW=n_IPRW,n_known=n_known,n_approx=n_approx,
                            mu_1=0.475,mu_0=0.375,mu_x=0,
                            sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
                            beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64,
                            delta=delta,m=5)
  n_out <- c(delta,
             ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)),
             ceiling(n_IPRW*kappa)+ceiling(n_IPRW*(1-kappa)),
             ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)),
             ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
  names(n_out) <- c("delta","n_standard","n_IPRW","n_known","n_approx")
  results <- c(n_out,sim_out[1:5],sim_out[6:10]*100)
  combine <- rbind(combine,results)
}
combine <- data.frame(combine)

# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)

fig_n <- ggplot() + 
  theme_classic() +
  geom_point(data=combine,aes(x=delta,y=n_standard,color="iv. C-standard",shape="iv. C-standard")) + 
  geom_point(data=combine,aes(x=delta,y=n_known,color="ii. C-known",shape="ii. C-known")) + 
  geom_point(data=combine,aes(x=delta,y=n_approx,color="iii. C-approx",shape="iii. C-approx")) +
  geom_point(data=combine,aes(x=delta,y=n_IPRW,color="i. C-IPRW",shape="i. C-IPRW")) + 
  scale_y_continuous(breaks=c(1000,1500,2000,2500),lim=c(1000,2600)) +
  scale_x_continuous(breaks=c(0.03,0.05,0.07)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\delta$"),y="n",color="") + theme(plot.title = element_text(size=10))
fig_n 
fig_power <- ggplot() + 
  theme_classic() +
  geom_point(data=combine,aes(x=delta,y=empirical.power..standard.eqn.4.,color="iv. C-standard",shape="iv. C-standard")) + 
  geom_point(data=combine,aes(x=delta,y=empirical.power..known.eqn.4.,color="ii. C-known",shape="ii. C-known")) +
  geom_point(data=combine,aes(x=delta,y=empirical.power..approx.eqn.4.,color="iii. C-approx",shape="iii. C-approx")) +
  geom_point(data=combine,aes(x=delta,y=empirical.power..IPRW.eqn.4.,color="i. C-IPRW",shape="i. C-IPRW")) +
  scale_y_continuous(breaks=c(75,80,85,90,95,100), lim=c(74,100)) +
  scale_x_continuous(breaks=c(0.03,0.05,0.07)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$\\delta$"),y="power",color="") + theme(plot.title = element_text(size=10))
fig_power


legend <- cowplot::get_legend(fig_n + theme(legend.position = "bottom"))
out_plot <- aplot::plot_list(fig_n + theme(legend.position = "none"),
                             ncol = 1, nrow = 1, labels="(a)")
aplot::plot_list(out_plot,legend,ncol=1,heights=c(3,0.1))
ggsave("Fig_CRT_cont_n_delta.pdf",scale=0.5)

legend <- cowplot::get_legend(fig_power + theme(legend.position = "bottom"))
out_plot <- aplot::plot_list(fig_power + theme(legend.position = "none"),
                             ncol = 1, nrow = 1, labels="(b)")
aplot::plot_list(out_plot,legend,ncol=1,heights=c(3,0.1))
ggsave("Fig_CRT_cont_power_delta.pdf",scale=0.5)


