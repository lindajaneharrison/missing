# Figure 2 (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

# Weighting for a fully observed auxiliary continuous variable

library(ggplot2)
library(latex2exp)
library(ggpubr)

cbPalette_web <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7", "#F0E442", "#0072B2", "#D55E00")
nice_lines_web <- c(1,2,3,4)

cbPalette <- c("#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#0072B2", "#D55E00")
nice_lines <- c(1,2,4)

n_standard_out <- n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                             type="continuous",link="identity",sigma_y_sq=0.245,phi=0.8)
n_approx_out <- n_approx_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                         mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,
                         beta_01=1.4,beta_11=0.21,beta_00=1.45,beta_10=0.47)

n_IPRW_out <- rep(NA,length(seq(-0.99,0.99,0.01)))
n_known_out <- rep(NA,length(seq(-0.99,0.99,0.01)))
i <- 0
for (rho in seq(-0.99,0.99,0.01)) {
  i <- i+1
  n_IPRW_out[i] <- n_IPRW_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                            mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=rho,
                            beta_01=1.4,beta_11=0.21,beta_00=1.45,beta_10=0.47)
  n_known_out[i] <- n_known_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                            mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=rho,
                            beta_01=1.4,beta_11=0.21,beta_00=1.45,beta_10=0.47)
}

rho <- seq(-0.99,0.99,0.01)
df <- as.data.frame(cbind(rho,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

a_web <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=rho,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=rho,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_x_continuous(breaks=c(-0.9,-0.5,0,0.5,0.9)) + 
  scale_y_continuous(limits=c(1000,2100),breaks=c(1000,1200,1400,1600,1800,2000)) + 
  labs(x = TeX("$\\rho$"),y="n", color="",
       caption=TeX("$\\beta_{01}=1.4$, $\\beta_{11}=0.21$, $\\beta_{00}=1.45$, $\\beta_{10}=0.47$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette_web) + 
  scale_linetype_manual(values=nice_lines_web) + 
  guides(color = guide_legend(override.aes = list(linetype = nice_lines_web)), linetype=FALSE)

a <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=rho,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=rho,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iii. standard",linetype="iii. standard")) +
  scale_x_continuous(breaks=c(-0.9,-0.5,0,0.5,0.9)) + 
  scale_y_continuous(limits=c(1000,2100),breaks=c(1000,1200,1400,1600,1800,2000)) + 
  labs(x = TeX("$\\rho$"),y="n", color="",
       caption=TeX("$\\beta_{01}=1.4$, $\\beta_{11}=0.21$, $\\beta_{00}=1.45$, $\\beta_{10}=0.47$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)


n_approx_out <- n_approx_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                         mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,
                         beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64)

n_IPRW_out <- rep(NA,length(seq(-0.99,0.99,0.01)))
n_known_out <- rep(NA,length(seq(-0.99,0.99,0.01)))
i <- 0
for (rho in seq(-0.99,0.99,0.01)) {
  i <- i+1
  n_IPRW_out[i] <- n_IPRW_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                        mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=rho,
                        beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64)
  n_known_out[i] <- n_known_cont(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                            mu_x=0,sigma_x_sq=1,sigma_y_sq=0.245,rho=rho,
                            beta_01=1.4,beta_11=0.21,beta_00=2,beta_10=1.64)
}

rho <- seq(-0.99,0.99,0.01)
df <- as.data.frame(cbind(rho,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

b_web <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=rho,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=rho,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_x_continuous(breaks=c(-0.9,-0.5,0,0.5,0.9)) + 
  scale_y_continuous(limits=c(1000,2100),breaks=c(1000,1200,1400,1600,1800,2000)) + 
  labs(x = TeX("$\\rho$"),y="n", color="",
       caption=TeX("$\\beta_{01}=1.4$, $\\beta_{11}=0.21$, $\\beta_{00}=2$, $\\beta_{10}=1.64$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette_web) + 
  scale_linetype_manual(values=nice_lines_web) + 
  guides(color = guide_legend(override.aes = list(linetype = nice_lines_web)), linetype=FALSE)

b <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=rho,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=rho,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iii. standard",linetype="iii. standard")) +
  scale_x_continuous(breaks=c(-0.9,-0.5,0,0.5,0.9)) + 
  scale_y_continuous(limits=c(1000,2100),breaks=c(1000, 1200,1400,1600,1800,2000)) + 
  labs(x = TeX("$\\rho$"),y="n", color="",
       caption=TeX("$\\beta_{01}=1.4$, $\\beta_{11}=0.21$, $\\beta_{00}=2$, $\\beta_{10}=1.64$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)


ggarrange(b,a,common.legend = TRUE, legend = "right")
ggsave("Figure_2.pdf", width = 7, height = 3.5)

ggarrange(b_web,a_web,common.legend = TRUE, legend = "right")
ggsave("Figure_Web2.pdf", width = 7, height = 3.5)


