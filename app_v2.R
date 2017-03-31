#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Cairo)
options(shiny.usecairo=T)
source("app_functions.R")
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Am I number or am I color?"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("n",
                     "Number of samples:",
                     min = 50,
                     max = 800,
                     value = 150),
         actionButton("upN", "Change sample size"),
         sliderInput("bins",
                     "Complexity:",
                     min = 1,
                     max = 26,
                     value = 1),
         h3("da model"),
        verbatimTextOutput("modx")
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot"),
         #verbatimTextOutput("outout"),
         plotOutput("fixPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  myData = dataGen(n=150)
  observeEvent(input$upN, {
    myData = dataGen(n=trunc(input$n))
    # temp <- values$df_data
    # temp = c(temp, shiPro(3,input$name,input$descr))
    #ds$val[[nrow(ds$val)+1]] <- shiPro(3,input$name,input$descr)
    # ds$val = rbind(ds$val,data.frame("id"=nrow(ds$val)+1, "name"=input$name, "descr"=input$descr))
    # write.table(ds$val, "14.txt",row.names = F)
  })
  #nSlider = reactive({dataGen(n=trunc(input$n))})

  
  
  parList = c("x1", "x2", "x3", "I(x1^2)","I(x1^3)",
              "x4","x1*x2","I(x2^2)","I(x2^3)",
              "x1*x3","I(x3^2)","I(x3^3)",
              "x2*x3")
  parList = c(parList,paste0("I(x1^",4:16,")")
              )
  newmula = "y~1"
  xModels = c(newmula)
  for (i in 1:26){
    newmula = paste(newmula, parList[i],sep=" + ")
    xModels = c(xModels,newmula)
  }
  
  # #poly only
  # parList = paste0("I(x1^",2:16,")")
  # newmula = "y~1+x1"
  # xModels = c(newmula)
  # for (i in 1:15){
  #   newmula = paste(newmula, parList[i],sep=" + ")
  #   xModels = c(xModels,newmula)
  # }
  
  
  
  output$modx <- renderPrint({xModels[input$bins]})
  output$distPlot <- renderPlot({
    input$upN
    myn <- isolate(input$n)
      # # generate bins based on input$bins from ui.R
      # x    <- faithful[, 2] 
      # bins <- seq(min(x), max(x), length.out = input$bins + 1)
      # 
      # # draw the histogram with the specified number of bins
      # hist(x, breaks = bins, col = 'darkgray', border = 'white')
     jaVal(as.formula(xModels[input$bins]), myData,folds = 5,repetitions = 10,plot=T)
   })
  ########################
  # All together
  ########################
  # vRMSEtrain = c()
  # vRMSEtest = c()
  # vR2train = c()
  # vR2test = c()
  # for(i in 1:length(xModels)){
  #   temp = jaVal(as.formula(xModels[i]), dataGen(n=trunc(150)),folds = 5,repetitions = 10,plot=F)
  #   vRMSEtrain = c(vRMSEtrain, temp$trainRmse)
  #   vRMSEtest = c(vRMSEtest, temp$testRmse)
  #   vR2train = c(vR2train, temp$trainR2)
  #   vR2test = c(vR2test, temp$testR2)
  # }
  # myYlim = c(
  #   min(min(vRMSEtrain,na.rm=T),min(vRMSEtest,na.rm=T)),
  #   max(max(vRMSEtrain,na.rm=T),max(vRMSEtest,na.rm=T))
  # )
  # myYlim = c(0,1)
  # print(myYlim)
  # output$fixPlot = renderPlot({
  #   plot(0, type="n",xlim=c(0,length(xModels)), ylim=myYlim)
  #   lines(vR2train,col="red")
  #   lines(vR2test,col="blue")
  # })
  
  output$fixPlot = renderPlot({
   # multiGen(input$n,150,xModels)
    input$upN
    myn <- isolate(input$n)
    multiGen(100,myn,xModels[1:18])
  })
  
  output$outout <- renderPrint({
    #print(vRMSEtest)
    #jaVal(as.formula(xModels[input$bins]), dataGen(n=trunc(input$n)),folds = 5,repetitions = 10,plot=F)$trainR
    
  })
}




# Run the application 
shinyApp(ui = ui, server = server)