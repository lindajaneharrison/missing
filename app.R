# Shiny App (Sample Size Calculation for Randomized Clinical Trials via Inverse Probability of Response Weighting 
# when Outcome Data are Missing at Random)
source('shinyFUNCTIONS.R')

# R libraries needed
library(shiny)
library(shinyjs)
library(ggplot2)
library(latex2exp)
library(fastGHQuad)

# Gauss-Hermite 100 quadrature points
x_j <- gaussHermiteData(100)$x 
w_j <- gaussHermiteData(100)$w

# Define UI app ----
ui <- 
navbarPage("Sample size calculation for randomized clinical trials via inverse probability of response weighting when outcome data are missing at random",
tabPanel("weighting by a categorical variable",
fluidPage(
  
  shinyjs::useShinyjs(),
  
  tags$style("[type='number'] {font-size:11px;height:20px;}"),
  # Sidebar panel for inputs ----
  fluidRow(
    
    # Input: Selector for variables ----
    column(4,
    sliderInput("power", 
                 "Power [\\(1-\\beta\\)]",
                 min=0.5,max=1,step=0.01,value=0.9,width='100%',ticks=FALSE),
    sliderInput("alpha", 
                 "Type I error [\\(\\alpha\\)]",
                 min=0.01,max=0.1,step=0.01,value=0.05,width='100%',ticks=FALSE),
    radioButtons("type", 
                 "Outcome type",
                 choices = list("continuous" = "continuous", "binary" = "binary"),selected = "continuous",
                 width='100%'),
    radioButtons("link", 
                 "Link function",
                 choices = list("identity" = "identity", "logit" = "logit"),selected = "identity",
                 width='100%'),
    numericInput("mu_11", 
                 "Mean outcome, group 1, intervention [\\(\\mu_{11}\\)]",
                 value=0.9,width='100%'),
    numericInput("mu_21", 
                 "Mean outcome, group 2, intervention [\\(\\mu_{21}\\)]",
                 value=0.3,width='100%'),
    numericInput("mu_10", 
                 "Mean outcome, group 1, control [\\(\\mu_{10}\\)]",
                 value=0.15,width='100%'),
    numericInput("mu_20", 
                 "Mean outcome, group 2, control [\\(\\mu_{20}\\)]",
                 value=0.85,width='100%'),
    numericInput("sigma_11_sq", 
                 "Variance outcome, group 1, intervention [\\(\\sigma_{11}^2\\)]",
                 min=0,value=0.026,width='100%'),
    numericInput("sigma_21_sq", 
                 "Variance outcome, group 2, intervention [\\(\\sigma_{21}^2\\)]",
                 min=0,value=0.294,width='100%'),
    numericInput("sigma_10_sq", 
                 "Variance outcome, group 1, control [\\(\\sigma_{10}^2\\)]",
                 min=0,value=0.026,width='100%'),
    numericInput("sigma_20_sq", 
                 "Variance outcome, group 2, control [\\(\\sigma_{20}^2\\)]",
                 min=0,value=0.229,width='100%')),
    column(4,
    sliderInput("kappa", 
                       "Proportion randomized to intervention [\\(\\kappa\\)]",
                       min=0,max=1,value=0.5,width='100%',ticks=FALSE),
    sliderInput("pi_1", 
                 "Proportion in group 1 [\\(\\pi_1\\)] (weighting by a binary categorical variable)",
                 min=0,max=1,step=0.01,value=0.5,width='100%',ticks=FALSE),
    sliderInput("expit_beta_11", 
                 "Proportion with observed outcome, group 1, intervention [expit(\\(\\beta_{11}\\))]",
                 min=0,max=1,step=0.01,value=0.7,width='100%',ticks=FALSE),
    sliderInput("expit_beta_21", 
                 "Proportion with observed outcome, group 2, intervention [expit(\\(\\beta_{21}\\))]",
                 min=0,max=1,step=0.01,value=0.9,width='100%',ticks=FALSE),
    sliderInput("expit_beta_10", 
                 "Proportion with observed outcome, group 1, control group [expit(\\(\\beta_{10}\\))]",
                 min=0,max=1,step=0.01,value=0.75,width='100%',ticks=FALSE),
    sliderInput("expit_beta_20", 
                 "Proportion with observed outcome, group 2, control group [expit(\\(\\beta_{20}\\))]",
                 min=0,max=1,step=0.01,value=0.85,width='100%',ticks=FALSE),
    withMathJax(uiOutput("standard_inputs",width='100%')),
    withMathJax(uiOutput("standard_inputs_extra",width='100%')),
    radioButtons("CRT", 
                 "Individually or Cluster Randomized Trial",
                 choices = list("individually" = "IRT", "cluster" = "CRT"),selected = "IRT",
                 width='100%'),
    sliderInput("m", 
                "cluster size [m]",
                min=0,max=200,step=1,value=5,width='100%',ticks=FALSE),
    sliderInput("delta", 
                "intercluster correlation [\\(\\delta\\)]",
                min=0,max=1,step=0.01,value=0.05,width='100%',ticks=FALSE)),
  
  # Main panel for displaying outputs ----
  column(4,
    plotOutput("plot",width='100%'),
    withMathJax(uiOutput("n_IPRW_text",width='100%')),
    verbatimTextOutput("n_IPRW"),
    withMathJax(uiOutput("n_known_text",width='100%')),
    verbatimTextOutput("n_known"),
    withMathJax(uiOutput("n_approx_text",width='100%')),
    verbatimTextOutput("n_approx"),
    withMathJax(uiOutput("n_standard_text",width='100%')),
    verbatimTextOutput("n_standard"),
    withMathJax(uiOutput("ratio_text",width='100%')),
    verbatimTextOutput("ratio")
  )
  
  )
)
),
tabPanel("weighting by a continuous variable",
         fluidPage(
           tags$style("[type='number'] {font-size:11px;height:20px;}"),
           # Sidebar panel for inputs ----
           fluidRow(
             
             # Input: Selector for variables ----
             column(4,
                    sliderInput("power_cont", 
                                "Power [\\(1-\\beta\\)]",
                                min=0.5,max=1,step=0.01,value=0.9,width='100%',ticks=FALSE),
                    sliderInput("alpha_cont", 
                                "Type I error [\\(\\alpha\\)]",
                                min=0.01,max=0.1,step=0.01,value=0.05,width='100%',ticks=FALSE),
                    numericInput("mu_1_cont", 
                                 "Mean outcome, intervention [\\(\\mu_1\\)]",
                                 value=0.475,width='100%'),
                    numericInput("mu_0_cont", 
                                 "Mean outcome, control [\\(\\mu_0\\)]",
                                 value=0.375,width='100%'),
                    numericInput("sigma_y_sq_cont", 
                                 "Variance outcome [\\(\\sigma_y^2\\)]",
                                 min=0,value=0.245,width='100%'),
                    numericInput("mu_x", 
                                 "Mean covariate [\\(\\mu_x\\)]",
                                 value=0,width='100%'),
                    numericInput("sigma_x_sq", 
                                 "Variance covariate [\\(\\sigma_x^2\\)]",
                                 min=0,value=1,width='100%'),
                    sliderInput("rho", 
                                "Correlation covariate and outcome [\\(\\rho\\)]",
                                min=-1,max=1,step=0.05,value=-0.75,width='100%',ticks=FALSE),
                    "Assuming the outcome and covariate are normally distributed"),
             column(4,
                    sliderInput("kappa_cont", 
                                "Proportion randomized to intervention [\\(\\kappa\\)]",
                                min=0,max=1,value=0.5,width='100%',ticks=FALSE),
                    numericInput("beta_01", 
                                 "Log odds of observation when x=0, intervention [\\(\\beta_{01}\\)]",
                                 value=1.4,width='100%'),
                    numericInput("beta_11", 
                                 "Log odds ratio of observation for 1 unit higher x, intervention [\\(\\beta_{11}\\)]",
                                 value=0.21,width='100%'),
                    numericInput("beta_00", 
                                 "Log odds of observation when x=0, control [\\(\\beta_{00}\\)]",
                                 value=2,width='100%'),
                    numericInput("beta_10", 
                                 "Log odds ratio of observation for 1 unit higher x, control [\\(\\beta_{10}\\)]",
                                 value=1.64,width='100%'),
                    withMathJax(uiOutput("phi_cont",width='100%')),
                    withMathJax(uiOutput("beta_0_cont",width='100%')),
                    withMathJax(uiOutput("beta_1_cont",width='100%')),
                    radioButtons("CRT_cont", 
                                 "Individually or Cluster Randomized Trial",
                                 choices = list("individually" = "IRT", "cluster" = "CRT"),selected = "IRT",
                                 width='100%'),
                    sliderInput("m_cont", 
                                "cluster size [m]",
                                min=0,max=200,step=1,value=5,width='100%',ticks=FALSE),
                    sliderInput("delta_cont", 
                                "intercluster correlation [\\(\\delta\\)]",
                                min=0,max=1,step=0.001,value=0.05,width='100%',ticks=FALSE)),
             
             # Main panel for displaying outputs ----
             column(4,
                    plotOutput("plot_cont",width='100%'),
                    withMathJax(uiOutput("n_IPRW_cont_text",width='100%')),
                    verbatimTextOutput("n_IPRW_cont"),
                    withMathJax(uiOutput("n_known_cont_text",width='100%')),
                    verbatimTextOutput("n_known_cont"),
                    withMathJax(uiOutput("n_approx_cont_text",width='100%')),
                    verbatimTextOutput("n_approx_cont"),
                    withMathJax(uiOutput("n_standard_cont_text",width='100%')),
                    verbatimTextOutput("n_standard_cont"),
                    withMathJax(uiOutput("ratio_cont_text",width='100%')),
                    verbatimTextOutput("ratio_cont")
             )
             
           )
         )
)
)

# Colors and shapes for figures
cbPalette <- c("#E69F00", "#56B4E9", "#009E73","#CC79A7")
nice.shapes <- c(16,17,15,3)

# Define server logic ----
server <- function(input, output) {
  
  shinyjs::onclick("advanced2",
                   shinyjs::toggle(id = "advanced2", anim = TRUE))
  
  observeEvent(input$type, {
    if(input$type == "binary"){
      shinyjs::disable("sigma_11_sq")
      shinyjs::disable("sigma_21_sq")
      shinyjs::disable("sigma_10_sq")
      shinyjs::disable("sigma_20_sq")
    } else {
      shinyjs::enable("sigma_11_sq")
      shinyjs::enable("sigma_21_sq")
      shinyjs::enable("sigma_10_sq")
      shinyjs::enable("sigma_20_sq")
    }
  })
  
  observeEvent(input$CRT, {
    if(input$CRT == "IRT"){
      shinyjs::disable("m")
      shinyjs::disable("delta")
    } else {
      shinyjs::enable("m")
      shinyjs::enable("delta")
    }
  })
  
  observeEvent(input$CRT_cont, {
    if(input$CRT_cont == "IRT"){
      shinyjs::disable("m_cont")
      shinyjs::disable("delta_cont")
    } else {
      shinyjs::enable("m_cont")
      shinyjs::enable("delta_cont")
    }
  })
  
  pi_2 <- reactive(1-input$pi_1)
  mu_1 <- reactive(input$pi_1*input$mu_11+pi_2()*input$mu_21)
  mu_0 <- reactive(input$pi_1*input$mu_10+pi_2()*input$mu_20)
  sigma_y1_sq <- reactive(input$pi_1*(input$sigma_11_sq+(input$mu_11-mu_1())^2)+pi_2()*(input$sigma_21_sq+(input$mu_21-mu_1())^2))
  sigma_y0_sq <- reactive(input$pi_1*(input$sigma_10_sq+(input$mu_10-mu_0())^2)+pi_2()*(input$sigma_20_sq+(input$mu_20-mu_0())^2))
  sigma_y_sq <- reactive(if      (input$type=="binary")     {NULL}
                         else if (input$type=="continuous") {(1-input$kappa)*sigma_y1_sq()+input$kappa*sigma_y0_sq()})
  phi <- reactive(input$kappa*(input$pi_1*input$expit_beta_11+pi_2()*input$expit_beta_21)+
                  (1-input$kappa)*(input$pi_1*input$expit_beta_10+pi_2()*input$expit_beta_20))
  
  output$standard_inputs <- renderPrint({ 
    withMathJax(sprintf("$$\\mu_1= %0.2f, \\mu_0= %0.2f, \\phi= %0.2f$$",mu_1(),mu_0(),phi()))
  })
  output$standard_inputs_extra <- renderPrint({ 
    withMathJax(sprintf("$$\\sigma_y^2= %0.2f$$",sigma_y_sq()))
  })
  
  output$n_standard_text <- renderPrint({ if (input$CRT=="CRT") {
                                                  withMathJax("$$n_{standard} [K_{standard}]$$")
                                                }
                                           else { 
                                                  withMathJax("$$n_{standard}$$")
                                                }
                                            })
  n_standard_out <- reactive (if (input$CRT=="CRT") {ceiling(n_standard(power=input$power,alpha=input$alpha,kappa=input$kappa,
                                                         mu_1=mu_1(),mu_0=mu_0(),
                                                         type=input$type,link=input$link,
                                                         sigma_y_sq=sigma_y_sq(),
                                                         phi=phi(),
                                                         m=input$m,delta=input$delta))}
                              else                  {ceiling(n_standard(power=input$power,alpha=input$alpha,kappa=input$kappa,
                                                         mu_1=mu_1(),mu_0=mu_0(),
                                                         type=input$type,link=input$link,
                                                         sigma_y_sq=sigma_y_sq(),
                                                         phi=phi()))})
  K_standard_out <- reactive (if (input$CRT=="CRT")
                             {paste("[", ceiling((n_standard_out()*input$kappa)/input$m)+ceiling((n_standard_out()*(1-input$kappa))/input$m),"]")}
                             )
  output$n_standard <- renderText({ 
     paste(n_standard_out(),K_standard_out())
  })
  
  output$n_IPRW_text <- renderPrint({ if (input$CRT=="CRT") {
    withMathJax("$$n_{IPRW} [K_{IPRW}]$$")
  }
    else { 
    withMathJax("$$n_{IPRW}$$")
    }
  })
  n_IPRW <- reactive(if (input$CRT=="CRT") {ceiling(n_IPRW_cat(power=input$power,alpha=input$alpha,kappa=input$kappa,mu_1=mu_1(),mu_0=mu_0(),
                                            mu_11=input$mu_11,mu_21=input$mu_21,mu_10=input$mu_10,mu_20=input$mu_20,
                                            type=input$type,link=input$link,
                                            sigma_11_sq=input$sigma_11_sq,sigma_21_sq=input$sigma_21_sq,sigma_10_sq=input$sigma_10_sq,sigma_20_sq=input$sigma_20_sq,
                                            pi_1=input$pi_1,pi_2=pi_2(),
                                            expit_beta_11=input$expit_beta_11,expit_beta_10=input$expit_beta_10,expit_beta_21=input$expit_beta_21,expit_beta_20=input$expit_beta_20,
                                            m=input$m,delta=input$delta))}
                    else {ceiling(n_IPRW_cat(power=input$power,alpha=input$alpha,kappa=input$kappa,mu_1=mu_1(),mu_0=mu_0(),
                          mu_11=input$mu_11,mu_21=input$mu_21,mu_10=input$mu_10,mu_20=input$mu_20,
                          type=input$type,link=input$link,
                          sigma_11_sq=input$sigma_11_sq,sigma_21_sq=input$sigma_21_sq,sigma_10_sq=input$sigma_10_sq,sigma_20_sq=input$sigma_20_sq,
                          pi_1=input$pi_1,pi_2=pi_2(),
                          expit_beta_11=input$expit_beta_11,expit_beta_10=input$expit_beta_10,expit_beta_21=input$expit_beta_21,expit_beta_20=input$expit_beta_20))})
  K_IPRW <- reactive (if (input$CRT=="CRT")
           {paste("[", ceiling((n_IPRW()*input$kappa)/input$m)+ceiling((n_IPRW()*(1-input$kappa))/input$m),"]")}
           )
  output$n_IPRW <- renderText({ 
    paste(n_IPRW(),K_IPRW())
  })
  
  output$n_known_text <- renderPrint({ if (input$CRT=="CRT") {
    withMathJax("$$n_{known} [K_{known}]$$")
  }
    else { 
      withMathJax("$$n_{known}$$")
    }
  })
  n_known <- reactive(if (input$CRT=="CRT") {ceiling(n_known_cat(power=input$power,alpha=input$alpha,kappa=input$kappa,mu_1=mu_1(),mu_0=mu_0(),
                                             mu_11=input$mu_11,mu_21=input$mu_21,mu_10=input$mu_10,mu_20=input$mu_20,
                                             type=input$type,link=input$link,
                                             sigma_11_sq=input$sigma_11_sq,sigma_21_sq=input$sigma_21_sq,sigma_10_sq=input$sigma_10_sq,sigma_20_sq=input$sigma_20_sq,
                                             pi_1=input$pi_1,pi_2=pi_2(),
                                             expit_beta_11=input$expit_beta_11,expit_beta_10=input$expit_beta_10,expit_beta_21=input$expit_beta_21,expit_beta_20=input$expit_beta_20,
                                             m=input$m,delta=input$delta))}
                      else {ceiling(n_known_cat(power=input$power,alpha=input$alpha,kappa=input$kappa,mu_1=mu_1(),mu_0=mu_0(),
                            mu_11=input$mu_11,mu_21=input$mu_21,mu_10=input$mu_10,mu_20=input$mu_20,
                            type=input$type,link=input$link,
                            sigma_11_sq=input$sigma_11_sq,sigma_21_sq=input$sigma_21_sq,sigma_10_sq=input$sigma_10_sq,sigma_20_sq=input$sigma_20_sq,
                            pi_1=input$pi_1,pi_2=pi_2(),
                            expit_beta_11=input$expit_beta_11,expit_beta_10=input$expit_beta_10,expit_beta_21=input$expit_beta_21,expit_beta_20=input$expit_beta_20))})
  K_known <- reactive (if (input$CRT=="CRT")
  {paste("[", ceiling((n_known()*input$kappa)/input$m)+ceiling((n_known()*(1-input$kappa))/input$m),"]")}
  )
  output$n_known <- renderText({ 
    paste(n_known(),K_known())
  })
  
  output$n_approx_text <- renderPrint({ if (input$CRT=="CRT") {
    withMathJax("$$n_{approx} [K_{approx}]$$")
  }
    else { 
      withMathJax("$$n_{approx}$$")
    }
  })
  n_approx <-reactive(if (input$CRT=="CRT") {ceiling(n_approx_cat(power=input$power,alpha=input$alpha,kappa=input$kappa,mu_1=mu_1(),mu_0=mu_0(),
                                             type=input$type,link=input$link,
                                             sigma_y_sq=sigma_y_sq(),
                                             pi_1=input$pi_1,pi_2=pi_2(),
                                             expit_beta_11=input$expit_beta_11,expit_beta_10=input$expit_beta_10,expit_beta_21=input$expit_beta_21,expit_beta_20=input$expit_beta_20,
                                             m=input$m,delta=input$delta))}
                      else {ceiling(n_approx_cat(power=input$power,alpha=input$alpha,kappa=input$kappa,mu_1=mu_1(),mu_0=mu_0(),
                                    type=input$type,link=input$link,
                                    sigma_y_sq=sigma_y_sq(),
                                    pi_1=input$pi_1,pi_2=pi_2(),
                                    expit_beta_11=input$expit_beta_11,expit_beta_10=input$expit_beta_10,expit_beta_21=input$expit_beta_21,expit_beta_20=input$expit_beta_20))})
  K_approx <- reactive (if (input$CRT=="CRT")
  {paste("[", ceiling((n_approx()*input$kappa)/input$m)+ceiling((n_approx()*(1-input$kappa))/input$m),"]")}
  )
  output$n_approx <- renderText({ 
    paste(n_approx(),K_approx())
  })

  output$ratio_text <- renderPrint({ 
    withMathJax("$$n_{IPRW}/n_{standard}$$")
  })
  output$ratio <- renderText({round(n_IPRW()/n_standard_out(),3)})
  
  plot <- reactive(ggplot() + 
                   theme_classic() +
                   geom_point(aes(x=input$power,y=n_standard_out(),color="iv. standard",shape="iv. standard",size=1)) + 
                   geom_point(aes(x=input$power,y=n_IPRW(),color="i. IPRW",shape="i. IPRW",size=1)) + 
                   geom_point(aes(x=input$power,y=n_known(),color="ii. known",shape="ii. known",size=1)) +
                   geom_point(aes(x=input$power,y=n_approx(),color="iii. approx",shape="iii. approx",size=1)) + 
                   scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
                   guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes, size=3)), shape=FALSE, size=FALSE) +
                   theme(axis.title.y=element_text(angle = 0),text=element_text(size=14),legend.text = element_text(size = 7.5),legend.position="bottom") +
                   scale_x_continuous(breaks=c(input$power)) +
                   labs(x = TeX("Power ($1-\\beta$)"),y="n",color="")) 
  output$plot <- renderPlot({plot()})
  
  phi_cont <- reactive(input$kappa_cont*sum(plogis(input$beta_01+(sqrt(2)*sqrt(input$sigma_x_sq)*x_j+input$mu_x)*input$beta_11)*w_j)/sqrt(pi)
                      +(1-input$kappa_cont)*sum(plogis(input$beta_00+(sqrt(2)*sqrt(input$sigma_x_sq)*x_j+input$mu_x)*input$beta_10)*w_j)/sqrt(pi))
  
  exp_beta_01 <- reactive(exp(input$beta_01))
  exp_beta_11 <- reactive(exp(input$beta_11))
  exp_beta_00 <- reactive(exp(input$beta_00))
  exp_beta_10 <- reactive(exp(input$beta_10))
  
  output$phi_cont <- renderPrint({ 
    withMathJax(sprintf("$$\\phi= %0.2f$$",
                        phi_cont()))
  })
  output$beta_1_cont <- renderPrint({ 
    withMathJax(sprintf("$$exp(\\beta_{01})= %0.2f, exp(\\beta_{11})= %0.2f$$",
                        exp_beta_01(),exp_beta_11()))
  })
  output$beta_0_cont <- renderPrint({ 
    withMathJax(sprintf("$$exp(\\beta_{00})= %0.2f, exp(\\beta_{10})= %0.2f$$",
                        exp_beta_00(),exp_beta_10()))
  })
  
  output$n_standard_cont_text <- renderPrint({ if (input$CRT_cont=="CRT") {
    withMathJax("$$n_{standard} [K_{standard}]$$")
  }
    else { 
      withMathJax("$$n_{standard}$$")
    }
  })
  n_standard_cont_out <- reactive(if (input$CRT_cont=="CRT") {ceiling(n_standard(power=input$power_cont,alpha=input$alpha_cont,kappa=input$kappa_cont,mu_1=input$mu_1_cont,mu_0=input$mu_0_cont,
                                        type="continuous",link="identity",
                                        sigma_y_sq=input$sigma_y_sq_cont,phi=phi_cont(),
                                        m=input$m_cont,delta=input$delta_cont))}
                                  else {ceiling(n_standard(power=input$power_cont,alpha=input$alpha_cont,kappa=input$kappa_cont,mu_1=input$mu_1_cont,mu_0=input$mu_0_cont,
                                         type="continuous",link="identity",
                                         sigma_y_sq=input$sigma_y_sq_cont,phi=phi_cont()))})
  K_standard_cont_out <- reactive (if (input$CRT_cont=="CRT")
  {paste("[", ceiling((n_standard_cont_out()*input$kappa_cont)/input$m_cont)+ceiling((n_standard_cont_out()*(1-input$kappa_cont))/input$m_cont),"]")}
  )
  output$n_standard_cont <- renderText({ 
    paste(n_standard_cont_out(),K_standard_cont_out())
  })
  
  output$n_IPRW_cont_text <- renderPrint({ if (input$CRT_cont=="CRT") {
    withMathJax("$$n_{IPRW} [K_{IPRW}]$$")
  }
    else { 
      withMathJax("$$n_{IPRW}$$")
    }
  })
  n_IPRW_cont_out <- reactive(if (input$CRT_cont=="CRT") {ceiling(n_IPRW_cont(power=input$power_cont,alpha=input$alpha_cont,kappa=input$kappa_cont,mu_1=input$mu_1_cont,mu_0=input$mu_0_cont,
                                                          mu_x=input$mu_x,sigma_x_sq=input$sigma_x_sq,sigma_y_sq=input$sigma_y_sq_cont,rho=input$rho,
                                                          beta_01=input$beta_01,beta_11=input$beta_11,beta_00=input$beta_00,beta_10=input$beta_10,
                                                          m=input$m_cont,delta=input$delta_cont))}
                              else {ceiling(n_IPRW_cont(power=input$power_cont,alpha=input$alpha_cont,kappa=input$kappa_cont,mu_1=input$mu_1_cont,mu_0=input$mu_0_cont,
                                                        mu_x=input$mu_x,sigma_x_sq=input$sigma_x_sq,sigma_y_sq=input$sigma_y_sq_cont,rho=input$rho,
                                                        beta_01=input$beta_01,beta_11=input$beta_11,beta_00=input$beta_00,beta_10=input$beta_10))})
  K_IPRW_cont_out <- reactive (if (input$CRT_cont=="CRT")
  {paste("[", ceiling((n_IPRW_cont_out()*input$kappa_cont)/input$m_cont)+ceiling((n_IPRW_cont_out()*(1-input$kappa_cont))/input$m_cont),"]")}
  )
  output$n_IPRW_cont <- renderText({ 
    paste(n_IPRW_cont_out(),K_IPRW_cont_out())
  })
  
  output$n_known_cont_text <- renderPrint({ if (input$CRT_cont=="CRT") {
    withMathJax("$$n_{known} [K_{known}]$$")
  }
    else { 
      withMathJax("$$n_{known}$$")
    }
  })
  n_known_cont_out <- reactive(if (input$CRT_cont=="CRT") {ceiling(n_known_cont(power=input$power_cont,alpha=input$alpha_cont,kappa=input$kappa_cont,mu_1=input$mu_1_cont,mu_0=input$mu_0_cont,
                                                           mu_x=input$mu_x,sigma_x_sq=input$sigma_x_sq,sigma_y_sq=input$sigma_y_sq_cont,rho=input$rho,
                                                           beta_01=input$beta_01,beta_11=input$beta_11,beta_00=input$beta_00,beta_10=input$beta_10,
                                                           m=input$m_cont,delta=input$delta_cont))}
                              else {ceiling(n_known_cont(power=input$power_cont,alpha=input$alpha_cont,kappa=input$kappa_cont,mu_1=input$mu_1_cont,mu_0=input$mu_0_cont,
                                    mu_x=input$mu_x,sigma_x_sq=input$sigma_x_sq,sigma_y_sq=input$sigma_y_sq_cont,rho=input$rho,
                                    beta_01=input$beta_01,beta_11=input$beta_11,beta_00=input$beta_00,beta_10=input$beta_10))})
  K_known_cont_out <- reactive (if (input$CRT_cont=="CRT")
  {paste("[", ceiling((n_known_cont_out()*input$kappa_cont)/input$m_cont)+ceiling((n_known_cont_out()*(1-input$kappa_cont))/input$m_cont),"]")}
  )
  output$n_known_cont <- renderText({ 
    paste(n_known_cont_out(),K_known_cont_out())
  })
  
  output$n_approx_cont_text <- renderPrint({ if (input$CRT_cont=="CRT") {
    withMathJax("$$n_{approx} [K_{approx}]$$")
  }
    else { 
      withMathJax("$$n_{approx}$$")
    }
  })
  n_approx_cont_out <- reactive(if (input$CRT_cont=="CRT") {ceiling(n_approx_cont(power=input$power_cont,alpha=input$alpha_cont,kappa=input$kappa_cont,mu_1=input$mu_1_cont,mu_0=input$mu_0_cont,
                                                            mu_x=input$mu_x,sigma_x_sq=input$sigma_x_sq,sigma_y_sq=input$sigma_y_sq_cont,
                                                            beta_01=input$beta_01,beta_11=input$beta_11,beta_00=input$beta_00,beta_10=input$beta_10,
                                                            m=input$m_cont,delta=input$delta_cont))}
                              else {ceiling(n_approx_cont(power=input$power_cont,alpha=input$alpha_cont,kappa=input$kappa_cont,mu_1=input$mu_1_cont,mu_0=input$mu_0_cont,
                                                          mu_x=input$mu_x,sigma_x_sq=input$sigma_x_sq,sigma_y_sq=input$sigma_y_sq_cont,
                                                          beta_01=input$beta_01,beta_11=input$beta_11,beta_00=input$beta_00,beta_10=input$beta_10))})
  K_approx_cont_out <- reactive (if (input$CRT_cont=="CRT")
  {paste("[", ceiling((n_approx_cont_out()*input$kappa_cont)/input$m_cont)+ceiling((n_approx_cont_out()*(1-input$kappa_cont))/input$m_cont),"]")}
  )
  output$n_approx_cont <- renderText({ 
    paste(n_approx_cont_out(),K_approx_cont_out())
  })
  
  output$ratio_cont_text <- renderPrint({ 
    withMathJax("$$n_{IPRW}/n_{standard}$$")
  })
  output$ratio_cont <- renderText({round(n_IPRW_cont_out()/n_standard_cont_out(),3)})
  
  plot_cont <- reactive(ggplot() + 
                     theme_classic() +
                     geom_point(aes(x=input$power_cont,y=n_standard_cont_out(),color="iv. standard",shape="iv. standard",size=1)) + 
                     geom_point(aes(x=input$power_cont,y=n_IPRW_cont_out(),color="i. IPRW",shape="i. IPRW",size=1)) + 
                     geom_point(aes(x=input$power_cont,y=n_known_cont_out(),color="ii. known",shape="ii. known",size=1)) + 
                     geom_point(aes(x=input$power_cont,y=n_approx_cont_out(),color="iii. approx",shape="iii. approx",size=1)) + 
                     scale_colour_manual(values=cbPalette) + scale_shape_manual(values=nice.shapes) +
                     guides(color=guide_legend(override.aes = list(color=cbPalette, shape=nice.shapes, size=3)), shape=FALSE, size=FALSE) +   
                     theme(axis.title.y=element_text(angle = 0),text=element_text(size=14),legend.text = element_text(size = 7.5),legend.position="bottom") +
                     scale_x_continuous(breaks=c(input$power_cont)) +
                     labs(x = TeX("Power ($1-\\beta$)"),y="n",color=""))
  output$plot_cont <- renderPlot({plot_cont()})
  
}

shinyApp(ui, server)

#runApp("~/Documents/year4/Rui/manuscript_3/sim_revision/paper_revision/programs/shinyapp")


