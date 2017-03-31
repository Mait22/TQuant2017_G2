library(shiny)
library(shinythemes)


shinyUI(fluidPage(theme = shinytheme("flatly"),
                  titlePanel("Data application"),
                  sidebarLayout(
                    sidebarPanel("Parameters:",width = 3,
                                 p(""),
                                 p(""),
                                 sliderInput("a_pram", "Value of parameter a:",
                                             min = 1, max = 5, value = 2.5,step = 0.1
                                 ),
                                 
                                 sliderInput("b_pram", "Value of parameter b:",
                                             min = 1, max = 5, value = 2.5,step = 0.1
                                 ),
                                 
                                 sliderInput("c_pram", "Value of parameter c:",
                                             min = 1, max = 5, value = 2.5,step = 0.1
                                 ),
                                 
                                 sliderInput("noise_pram", "Noise proportion:",
                                             min = 0.05, max = 1, value = 0.1,step = 0.05
                                 ),
                                 
                                 # sliderInput("n_pram", "N of datapoints to generate:",
                                 #             min = 5, max = 2000, value = 500,step = 50
                                 # ),
                                 
                                 sliderInput("n_sims", "N of simulation to run:",
                                             min = 5, max = 20, value = 10,step = 1
                                 ),
                                
                                 p(""),
                                 actionButton("go","Compute the paper results"),
                                 p("")
                    ),
                    
                    

                    
                    
                    mainPanel(
                      # tags$style(type="text/css",
                      #            ".shiny-output-error { visibility: hidden; }",
                      #            ".shiny-output-error:before { visibility: hidden; }"),
                      
                      
                      tabsetPanel(
                        
                        tabPanel("Description of the application",
                                 includeHTML("description.html")
                        ),
                        
                        
                      
                        
                        tabPanel("Forgetting curves and fitted models",
                                 plotOutput("plot.m1",width = 1000, height = 600),
                                 p(""),
                                 plotOutput("plot.m2",width = 1000, height = 600),
                                 p(""),
                                 plotOutput("plot.m3",width = 1000, height = 600)
                        ),
                        
                        tabPanel("Explore the fit indexies from simulations (raw data)",
                                 uiOutput("aic.sim.title"),
                                 uiOutput("aic.sim"),
                                 uiOutput("aic.sim.exp"),
                                 p(""),
                                 uiOutput("bic.sim.title"),
                                 uiOutput("bic.sim"),
                                 uiOutput("bic.sim.exp"),
                                 p(""),
                                 uiOutput("rmse.sim.title"),
                                 uiOutput("rmse.sim"),
                                 uiOutput("rmse.sim.exp"),
                                 p(""),
                                 uiOutput("pvaf.sim.title"),
                                 uiOutput("pvaf.sim"),
                                 uiOutput("pvaf.sim.exp")
                                 
                                 # DT::dataTableOutput('bic'),
                                 # p(""),
                                 # DT::dataTableOutput('rmse'),
                                 # p(""),
                                 # DT::dataTableOutput('pvaf')
                        ),
                        
                        tabPanel("Generalizability - assessed by cross-validation procedure",
                                 
                                 sliderInput("n",
                                             "Number of samples:",
                                             min = 50,
                                             max = 800,
                                             value = 150),
                                 p(""),
                                 sliderInput("bins",
                                             "Complexity:",
                                             min = 1,
                                             max = 26,
                                             value = 1),
                                 p(""),
                                 actionButton("upN", "Run cross-validation"),
                                 plotOutput("distPlot"),
                                 p(""),
                                 uiOutput("comp.text"),
                                 plotOutput("fixPlot")
                        )
                        
                        
                        
                        
                        
                        
                        
                      )))))


