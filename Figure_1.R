# Figure 1 (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
source('useFUNCTIONS.R')

# Weighting for a fully observed auxiliary categorical variable

# Y binary, g identity
library(ggplot2)
library(latex2exp)
library(ggpubr)

cbPalette_web <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7", "#F0E442", "#0072B2", "#D55E00")
nice_lines_web <- c(1,2,3,4)

cbPalette <- c("#E69F00", "#56B4E9","#CC79A7", "#F0E442", "#0072B2", "#D55E00")
nice_lines <- c(1,2,4)

n_standard_out <- n_standard(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                             type="binary",link="identity",phi=0.8)
n_approx_out <- n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                         type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)

n_IPRW_out <- rep(NA,length(seq(0.1,0.75,0.001)))
n_known_out <- rep(NA,length(seq(0.1,0.75,0.001)))
i <- 0
for (mu_11 in seq(0.1,0.75,0.001)) {
  i <- i+1
  mu_10 <- 0.35
  mu_1 <- 0.375
  mu_0 <- 0.275
  mu_21 <- 2*mu_1-mu_11
  mu_20 <- 2*mu_0-mu_10
  n_IPRW_out[i] <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                        mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                        type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
  n_known_out[i] <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                            mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                            type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
}

mu_11 <- seq(0.1,0.75,0.001)
df <- as.data.frame(cbind(mu_11,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

a_web <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=mu_11,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=mu_11,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_x_continuous(limits=c(0.1,0.75),breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_y_continuous(limits=c(1070,1205),breaks=c(1075,1100,1125,1150,1175,1200)) + 
  labs(x = TeX("$\\mu_{11}$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.35$, expit$(\\beta_{11})=0.7$, expit$(\\beta_{21})=0.9$, expit$(\\beta_{10})=0.75$, expit$(\\beta_{20})=0.85$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette_web) + 
  scale_linetype_manual(values=nice_lines_web) + 
  guides(color = guide_legend(override.aes = list(linetype = nice_lines_web)), linetype=FALSE)

a <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=mu_11,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=mu_11,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iii. standard",linetype="iii. standard")) +
  scale_x_continuous(limits=c(0.1,0.75),breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_y_continuous(limits=c(1070,1205),breaks=c(1075,1100,1125,1150,1175,1200)) + 
  labs(x = TeX("$\\mu_{11}$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.35$, expit$(\\beta_{11})=0.7$, expit$(\\beta_{21})=0.9$, expit$(\\beta_{10})=0.75$, expit$(\\beta_{20})=0.85$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)

n_IPRW_out <- rep(NA,length(seq(0.1,0.75,0.001)))
n_known_out <- rep(NA,length(seq(0.1,0.75,0.001)))
i <- 0
for (mu_11 in seq(0.1,0.75,0.001)) {
  i <- i+1
  mu_10 <- 0.5
  mu_1 <- 0.375
  mu_0 <- 0.275
  mu_21 <- 2*mu_1-mu_11
  mu_20 <- 2*mu_0-mu_10
  n_IPRW_out[i] <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                        mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                        type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
  n_known_out[i] <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                            mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                            type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
}

mu_11 <- seq(0.1,0.75,0.001)
df <- as.data.frame(cbind(mu_11,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

b_web <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=mu_11,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=mu_11,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_x_continuous(limits=c(0.1,0.75),breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_y_continuous(limits=c(1070,1205),breaks=c(1075,1100,1125,1150,1175,1200)) + 
  labs(x = TeX("$\\mu_{11}$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.5$, expit$(\\beta_{11})=0.7$, expit$(\\beta_{21})=0.9$, expit$(\\beta_{10})=0.75$, expit$(\\beta_{20})=0.85$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette_web) + 
  scale_linetype_manual(values=nice_lines_web) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines_web)), linetype=FALSE)

b <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=mu_11,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=mu_11,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iii. standard",linetype="iii. standard")) +
  scale_x_continuous(limits=c(0.1,0.75),breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_y_continuous(limits=c(1070,1205),breaks=c(1075,1100,1125,1150,1175,1200)) + 
  labs(x = TeX("$\\mu_{11}$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.5$, expit$(\\beta_{11})=0.7$, expit$(\\beta_{21})=0.9$, expit$(\\beta_{10})=0.75$, expit$(\\beta_{20})=0.85$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)


n_approx_out <- n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                         type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                         expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)

n_IPRW_out <- rep(NA,length(seq(0.1,0.75,0.001)))
n_known_out <- rep(NA,length(seq(0.1,0.75,0.001)))
i <- 0
for (mu_11 in seq(0.1,0.75,0.001)) {
  i <- i+1
  mu_10 <- 0.35
  mu_1 <- 0.375
  mu_0 <- 0.275
  mu_21 <- 2*mu_1-mu_11
  mu_20 <- 2*mu_0-mu_10
  n_IPRW_out[i] <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                        mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                        type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                        expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)
  n_known_out[i] <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                            mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                            type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                            expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)
}

mu_11 <- seq(0.1,0.75,0.001)
df <- as.data.frame(cbind(mu_11,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

c_web <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=mu_11,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=mu_11,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_x_continuous(limits=c(0.1,0.75),breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_y_continuous(limits=c(995,1180),breaks=c(1000,1025,1050,1075,1100,1125,1150,1175)) + 
  labs(x = TeX("$\\mu_{11}$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.35$, expit$(\\beta_{11})=0.9$, expit$(\\beta_{21})=0.7$, expit$(\\beta_{10})=0.85$, expit$(\\beta_{20})=0.75$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette_web) + 
  scale_linetype_manual(values=nice_lines_web) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines_web)), linetype=FALSE)

c <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=mu_11,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=mu_11,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iii. standard",linetype="iii. standard")) +
  scale_x_continuous(limits=c(0.1,0.75),breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_y_continuous(limits=c(995,1180),breaks=c(1000,1025,1050,1075,1100,1125,1150,1175)) + 
  labs(x = TeX("$\\mu_{11}$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.35$, expit$(\\beta_{11})=0.9$, expit$(\\beta_{21})=0.7$, expit$(\\beta_{10})=0.85$, expit$(\\beta_{20})=0.75$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)


n_IPRW_out <- rep(NA,length(seq(0.1,0.75,0.001)))
n_known_out <- rep(NA,length(seq(0.1,0.75,0.001)))
i <- 0
for (mu_11 in seq(0.1,0.75,0.001)) {
  i <- i+1
  mu_10 <- 0.5
  mu_1 <- 0.375
  mu_0 <- 0.275
  mu_21 <- 2*mu_1-mu_11
  mu_20 <- 2*mu_0-mu_10
  n_IPRW_out[i] <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                        mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                        type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                        expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)
  n_known_out[i] <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                            mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                            type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                            expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)
}

mu_11 <- seq(0.1,0.75,0.001)
df <- as.data.frame(cbind(mu_11,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

d_web <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=mu_11,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=mu_11,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_x_continuous(limits=c(0.1,0.75),breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_y_continuous(limits=c(995,1180),breaks=c(1000,1025,1050,1075,1100,1125,1150,1175)) + 
  labs(x = TeX("$\\mu_{11}$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.5$, expit$(\\beta_{11})=0.9$, expit$(\\beta_{21})=0.7$, expit$(\\beta_{10})=0.85$, expit$(\\beta_{20})=0.75$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette_web) + 
  scale_linetype_manual(values=nice_lines_web) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines_web)), linetype=FALSE)

d <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=mu_11,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=mu_11,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iii. standard",linetype="iii. standard")) +
  scale_x_continuous(limits=c(0.1,0.75),breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_y_continuous(limits=c(995,1180),breaks=c(1000,1025,1050,1075,1100,1125,1150,1175)) + 
  labs(x = TeX("$\\mu_{11}$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.5$, expit$(\\beta_{11})=0.9$, expit$(\\beta_{21})=0.7$, expit$(\\beta_{10})=0.85$, expit$(\\beta_{20})=0.75$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)

ggarrange(a_web,b_web,c_web,d_web,common.legend = TRUE, legend = "right")
ggsave("Figure_Web1.pdf", width = 7, height = 5)

ggarrange(a,b,c,d,common.legend = TRUE, legend = "right")
ggsave("Figure_1.pdf", width = 7, height = 5)



########
n_approx_out <- n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                         type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                         expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)

n_IPRW_out <- rep(NA,length(seq(0.1,0.75,0.001)))
n_known_out <- rep(NA,length(seq(0.1,0.75,0.001)))
i <- 0
for (mu_11 in seq(0.1,0.75,0.001)) {
  i <- i+1
  mu_10 <- 0.35
  mu_1 <- 0.375
  mu_0 <- 0.275
  mu_21 <- 2*mu_1-mu_11
  mu_20 <- 2*mu_0-mu_10
  n_IPRW_out[i] <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                        mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                        type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
  n_known_out[i] <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                            mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                            type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
}

mu_11 <- seq(0.1,0.75,0.001)
rho <- (mu_1-mu_11)/sqrt(mu_1*(1-mu_1))
df <- as.data.frame(cbind(mu_11,rho,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

a <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=rho,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=rho,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_y_continuous(limits=c(1070,1205),breaks=c(1075,1100,1125,1150,1175,1200)) + 
  labs(x = TeX("$\\rho$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.35$, expit$(\\beta_{11})=0.7$, expit$(\\beta_{21})=0.9$, expit$(\\beta_{10})=0.75$, expit$(\\beta_{20})=0.85$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)

n_IPRW_out <- rep(NA,length(seq(0.1,0.75,0.001)))
n_known_out <- rep(NA,length(seq(0.1,0.75,0.001)))
i <- 0
for (mu_11 in seq(0.1,0.75,0.001)) {
  i <- i+1
  mu_10 <- 0.5
  mu_1 <- 0.375
  mu_0 <- 0.275
  mu_21 <- 2*mu_1-mu_11
  mu_20 <- 2*mu_0-mu_10
  n_IPRW_out[i] <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                        mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                        type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                        expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
  n_known_out[i] <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                            mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                            type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                            expit_beta_11=0.7,expit_beta_10=0.75,expit_beta_21=0.9,expit_beta_20=0.85)
}

mu_11 <- seq(0.1,0.75,0.001)
rho <- (mu_1-mu_11)/sqrt(mu_1*(1-mu_1))
df <- as.data.frame(cbind(mu_11,rho,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

b <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=rho,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=rho,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_y_continuous(limits=c(1070,1205),breaks=c(1075,1100,1125,1150,1175,1200)) + 
  labs(x = TeX("$\\rho$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.5$, expit$(\\beta_{11})=0.7$, expit$(\\beta_{21})=0.9$, expit$(\\beta_{10})=0.75$, expit$(\\beta_{20})=0.85$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)


n_approx_out <- n_approx_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=0.375,mu_0=0.275,
                         type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                         expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)

n_IPRW_out <- rep(NA,length(seq(0.1,0.75,0.001)))
n_known_out <- rep(NA,length(seq(0.1,0.75,0.001)))
i <- 0
for (mu_11 in seq(0.1,0.75,0.001)) {
  i <- i+1
  mu_10 <- 0.35
  mu_1 <- 0.375
  mu_0 <- 0.275
  mu_21 <- 2*mu_1-mu_11
  mu_20 <- 2*mu_0-mu_10
  n_IPRW_out[i] <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                        mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                        type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                        expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)
  n_known_out[i] <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                            mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                            type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                            expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)
}

mu_11 <- seq(0.1,0.75,0.001)
rho <- (mu_1-mu_11)/sqrt(mu_1*(1-mu_1))
df <- as.data.frame(cbind(mu_11,rho,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

c <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=rho,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=rho,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_y_continuous(limits=c(995,1180),breaks=c(1000,1025,1050,1075,1100,1125,1150,1175)) + 
  labs(x = TeX("$\\rho$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.35$, expit$(\\beta_{11})=0.9$, expit$(\\beta_{21})=0.7$, expit$(\\beta_{10})=0.85$, expit$(\\beta_{20})=0.75$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)

n_IPRW_out <- rep(NA,length(seq(0.1,0.75,0.001)))
n_known_out <- rep(NA,length(seq(0.1,0.75,0.001)))
i <- 0
for (mu_11 in seq(0.1,0.75,0.001)) {
  i <- i+1
  mu_10 <- 0.5
  mu_1 <- 0.375
  mu_0 <- 0.275
  mu_21 <- 2*mu_1-mu_11
  mu_20 <- 2*mu_0-mu_10
  n_IPRW_out[i] <- n_IPRW_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                        mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                        type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                        expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)
  n_known_out[i] <- n_known_cat(power=0.9,alpha=0.05,kappa=0.5,mu_1=mu_1,mu_0=mu_0,
                            mu_11=mu_11,mu_21=mu_21,mu_10=mu_10,mu_20=mu_20,
                            type="binary",link="identity",pi_11=0.5,pi_21=0.5,pi_10=0.5,pi_20=0.5,
                            expit_beta_11=0.9,expit_beta_10=0.85,expit_beta_21=0.7,expit_beta_20=0.75)
}

mu_11 <- seq(0.1,0.75,0.001)
rho <- (mu_1-mu_11)/sqrt(mu_1*(1-mu_1))
df <- as.data.frame(cbind(mu_11,rho,n_standard_out,n_IPRW_out,n_known_out,n_approx_out))

d <- ggplot() + 
  theme_classic() +
  geom_line(data=df,aes(x=rho,y=n_IPRW_out,color="i. IPRW",linetype="i. IPRW")) +
  geom_line(data=df,aes(x=rho,y=n_known_out,color="ii. known",linetype="ii. known")) +
  geom_hline(data=df,aes(yintercept=n_approx_out,color="iii. approx",linetype="iii. approx")) +
  geom_hline(data=df,aes(yintercept=n_standard_out,color="iv. standard",linetype="iv. standard")) +
  scale_y_continuous(limits=c(995,1180),breaks=c(1000,1025,1050,1075,1100,1125,1150,1175)) + 
  labs(x = TeX("$\\rho$"),y="n", color="",
       caption=TeX("$\\mu_{10}=0.5$, expit$(\\beta_{11})=0.9$, expit$(\\beta_{21})=0.7$, expit$(\\beta_{10})=0.85$, expit$(\\beta_{20})=0.75$")) +
  theme(axis.title.y=element_text(angle = 0), plot.caption=element_text(size=6)) +
  scale_colour_manual(values=cbPalette) + 
  scale_linetype_manual(values=nice_lines) + 
  theme(axis.title.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(linetype = nice_lines)), linetype=FALSE)

ggarrange(a,b,c,d,common.legend = TRUE, legend = "right")
ggsave("Figure_1_rho.pdf", width = 7, height = 5)

