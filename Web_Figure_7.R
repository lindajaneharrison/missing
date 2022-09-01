# Web Figure 7 (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

### Individually Randomized Trials (IRTs)

# Weighting for a fully observed auxiliary continuous variable

# Simulations
# set-up parallel and reproducible
library(doParallel)
library(doRNG)
numCores <- detectCores()
registerDoParallel(numCores)
set.seed(5376758)

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
x_j <- gaussHermiteData(100)$x 
w_j <- gaussHermiteData(100)$w
sigma_x_sq <- 1
mu_x <- 0
kappa <- 0.5
for (miss_int in c(0.75,0.775,0.8)){
  print(miss_int)
  miss_con <- 1.6-miss_int
  f1 <- function(beta_11) sum(plogis(1.4+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_11)*w_j)/sqrt(pi)-miss_int
  f0 <- function(beta_10) sum(plogis(2+(sqrt(2)*sqrt(sigma_x_sq)*x_j+mu_x)*beta_10)*w_j)/sqrt(pi)-miss_con
  beta_11 <- uniroot(f1, c(0,100))$root
  beta_10 <- uniroot(f0, c(0,100))$root
  print(beta_11)
  print(beta_10)
  n_standard_out <- n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
           type="continuous",link="identity",sigma_y_sq=0.245,phi=0.8)
  n_IPRW <- n_IPRW_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                    mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
                    beta_01=1.4,beta_11=beta_11,beta_00=2,beta_10=beta_10)
  n_known <- n_known_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                    mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
                    beta_01=1.4,beta_11=beta_11,beta_00=2,beta_10=beta_10)
  n_approx <- n_approx_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.475,mu_0=0.375,
                    mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,
                    beta_01=1.4,beta_11=beta_11,beta_00=2,beta_10=beta_10)
  sim_out <- sim_power_cont(rep=rep,kappa=0.5,n_standard=n_standard_out,n_IPRW=n_IPRW,n_known=n_known,n_approx=n_approx,
               mu_1=0.475,mu_0=0.375,mu_x=0,
               sigma_x_sq=1,sigma_y_sq=0.245,rho=-0.75,
               beta_01=1.4,beta_11=beta_11,beta_00=2,beta_10=beta_10)
  n_out <- c(miss_int,
           ceiling(n_standard_out*kappa)+ceiling(n_standard_out*(1-kappa)),
           ceiling(n_IPRW*kappa)+ceiling(n_IPRW*(1-kappa)),
           ceiling(n_known*kappa)+ceiling(n_known*(1-kappa)),
           ceiling(n_approx*kappa)+ceiling(n_approx*(1-kappa)))
  names(n_out) <- c("miss_int","n_standard","n_IPRW","n_known","n_approx")
  results <- c(n_out,sim_out[1:5],sim_out[6:10]*100)
  combine <- rbind(combine,results)
}
combine <- data.frame(combine)

# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)

fig_n <- ggplot() + 
  theme_classic() +
  geom_point(data=combine,aes(x=miss_int,y=n_standard,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine,aes(x=miss_int,y=n_known,color="ii. known",shape="ii. known")) + 
  geom_point(data=combine,aes(x=miss_int,y=n_approx,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine,aes(x=miss_int,y=n_IPRW,color="i. IPRW",shape="i. IPRW")) + 
  scale_y_continuous(breaks=c(1000,1500,2000),lim=c(1000,2000)) +
  scale_x_continuous(breaks=c(0.75,0.775,0.8)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$P(R_i|Z_i=1)$"),y="n",color="") + theme(plot.title = element_text(size=10))
fig_n 
fig_power <- ggplot() + 
  theme_classic() +
  geom_point(data=combine,aes(x=miss_int,y=empirical.power..standard.eqn.4.,color="iv. standard",shape="iv. standard")) + 
  geom_point(data=combine,aes(x=miss_int,y=empirical.power..known.eqn.4.,color="ii. known",shape="ii. known")) +
  geom_point(data=combine,aes(x=miss_int,y=empirical.power..approx.eqn.4.,color="iii. approx",shape="iii. approx")) +
  geom_point(data=combine,aes(x=miss_int,y=empirical.power..IPRW.eqn.4.,color="i. IPRW",shape="i. IPRW")) +
  scale_y_continuous(breaks=c(75,80,85,90,95,100), lim=c(74,100)) +
  scale_x_continuous(breaks=c(0.75,0.775,0.8)) +
  scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
  guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes)), shape=FALSE) +
  labs(x = TeX("$P(R_i|Z_i=1)$"),y="power",color="") + theme(plot.title = element_text(size=10))
fig_power


legend <- cowplot::get_legend(fig_n + theme(legend.position = "bottom"))
out_plot <- aplot::plot_list(fig_n + theme(legend.position = "none"),
                             ncol = 1, nrow = 1, labels="(a)")
aplot::plot_list(out_plot,legend,ncol=1,heights=c(3,0.1))
ggsave("Fig_IRT_cont_n_miss_int.pdf",scale=0.5)

legend <- cowplot::get_legend(fig_power + theme(legend.position = "bottom"))
out_plot <- aplot::plot_list(fig_power + theme(legend.position = "none"),
                             ncol = 1, nrow = 1, labels="(b)")
aplot::plot_list(out_plot,legend,ncol=1,heights=c(3,0.1))
ggsave("Fig_IRT_cont_power_miss_int.pdf",scale=0.5)




