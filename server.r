source("functions.R")

library(Cairo)
library(DT)
library(shiny)
options(shiny.usecairo=TRUE)







#Shiny server function
shinyServer(function(input,output,session){
  
  
  main.results <- eventReactive(input$go,{
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait, generating the data and computing the results", value = 0)
    
    results <- paper_data_gen(a = input$a_pram,
                              b = input$b_pram,
                              c = input$c_pram,
                              n.total = 5,
                              noise.prop = input$noise_pram)
    
    return(results)
  })
  
  
  graphs <- eventReactive(input$go,{
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait, generating the graphs", value = 0)
    
    results <- paper_data_fitter(main.results())
    
    return(results)
  })
  
  
  output$plot.m1 <- renderPlot({graphs()@data.1.comb.plot})
  output$plot.m2 <- renderPlot({graphs()@data.2.comb.plot})  
  output$plot.m3 <- renderPlot({graphs()@data.3.comb.plot})  
  
  
  
  sims <- eventReactive(input$go,{
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait, running simulations", value = 0)
    
    results <- simulator_aicbic(n.sims = input$n_sims,
                                n.total = 5,
                                noise.prop = input$noise_pram,
                                a = input$a_pram,
                                b = input$b_pram,
                                c = input$c_pram
                                )
    return(results)
  })
  
  
  aic.sim <- eventReactive(input$go,{
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait, running simulations and generating neat output tables!", value = 0)
    
    results <- data_table_output(sims()@raw.results.aic, colors.set = c("white","#4CA24C"))
    return(results)
  })
  
  output$aic.sim <- renderUI(aic.sim())
  
  
  
  bic.sim <- eventReactive(input$go,{
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait, running simulations and generating neat output tables!", value = 0)
    
    results <- data_table_output(sims()@raw.results.bic, colors.set = c("white","#4CA24C"))
    return(results)
  })
  
  output$bic.sim <- renderUI(bic.sim())
  
  
  rmse.sim <- eventReactive(input$go,{
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait, running simulations and generating neat output tables!", value = 0)
    
    results <- data_table_output(sims()@raw.results.rmse, colors.set = c("white","#FFA500"))
    return(results)
  })
  
  output$rmse.sim <- renderUI(rmse.sim())
  
  
  pvaf.sim <- eventReactive(input$go,{
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait, running simulations and generating neat output tables!", value = 0)
    
    results <- data_table_output(sims()@raw.results.pvaf, colors.set = c("#FFA500","white"))
    return(results)
  })
  
  output$pvaf.sim <- renderUI(pvaf.sim())
  
  # output$aic = DT::renderDataTable(sims()@raw.results.aic,
  #                                  server = TRUE, 
  #                                  options = list(dom = 'C<"clear">lfrtip'))
  # 
  # output$bic = DT::renderDataTable(sims()@raw.results.bic,
  #                                  server = TRUE, 
  #                                  options = list(dom = 'C<"clear">lfrtip'))
  # 
  # output$rmse = DT::renderDataTable(sims()@raw.results.rmse,
  #                                  server = TRUE, 
  #                                  options = list(dom = 'C<"clear">lfrtip'))
  # 
  # output$pvaf = DT::renderDataTable(sims()@raw.results.pvaf,
  #                                  server = TRUE, 
  #                                  options = list(dom = 'C<"clear">lfrtip'))
  
  
  
  
  
  aic.sim.exp <- eventReactive(input$go,{
    return(paste('<p align="justify"><font size="4"> Note:  <b> Model fit assessed by AIC. FULLER green colour of the cell indicates better fit.</b> <br></font></p>',sep = ""))
    
  })
  
  bic.sim.exp <- eventReactive(input$go,{
    return(paste('<p align="justify"><font size="4"> Note:  <b> Model fit assessed by BIC. FULLER green colour of the cell indicates better fit.</b> <br></font></p>',sep = ""))
    
  })
  
  rmse.sim.exp <- eventReactive(input$go,{
    return(paste('<p align="justify"><font size="4"> Note:  <b> Model fit assessed by RMSE. FULLER yellow colour of the cell indicates better fit.</b><br> </font></p>',sep = ""))
    
  })
  
  pvaf.sim.exp <- eventReactive(input$go,{
    return(paste('<p align="justify"><font size="4"> Note:  <b> Model fit assessed by PVAF. FULLER yellow colour of the cell indicates better fit.</b> <br></font></p>',sep = ""))
    
  })
  
  
  output$aic.sim.exp <- renderUI({
      HTML(aic.sim.exp())
  })
  
  output$bic.sim.exp <- renderUI({
    HTML(bic.sim.exp())
  })
  
  output$rmse.sim.exp <- renderUI({
    HTML(rmse.sim.exp())
  })
  
  output$pvaf.sim.exp <- renderUI({
    HTML(pvaf.sim.exp())
  })
  
  
  aic.sim.title <- eventReactive(input$go,{
    return(paste('<p align="justify"><font size="11"><br><b>Model fit assessed by AIC </b></font></p>',sep = ""))
    
  })
  
  bic.sim.title <- eventReactive(input$go,{
    return(paste('<p align="justify"><font size="11"><br><b>Model fit assessed by BIC</b></font></p>',sep = ""))
    
  })
  
  rmse.sim.title <- eventReactive(input$go,{
    return(paste('<p align="justify"><font size="11"><br><b>Model fit assessed by RMSE</b></font></p>',sep = ""))
    
  })
  
  pvaf.sim.title <- eventReactive(input$go,{
    return(paste('<p align="justify"><font size="11"><br><b>Model fit assessed by PVAF</b></font></p>',sep = ""))
    
  })
  
  
  output$aic.sim.title <- renderUI({
    HTML(aic.sim.title())
  })
  
  output$bic.sim.title <- renderUI({
    HTML(bic.sim.title())
  })
  
  output$rmse.sim.title <- renderUI({
    HTML(rmse.sim.title())
  })
  
  output$pvaf.sim.title <- renderUI({
    HTML(pvaf.sim.title())
  })
  
  
  
  
  
  
  
  
  
  
  
  #Johann's and Javier's app
  
  myData <- eventReactive(input$upN,{
    return(dataGen(n=trunc(input$n)))
  })
  
  
  distPlot  <- eventReactive(input$upN,{
    
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
    
  return(jaVal(as.formula(xModels[input$bins]), myData(),folds = 5,repetitions = 10,plot=T))  
  })
  
  output$distPlot <-  renderPlot({distPlot()})
  
  
  
  
  fixPlot   <- eventReactive(input$upN,{
    if(input$bins <4){
      return(NULL)
    }

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

    return(multiGen(100,input$n,xModels[1:input$bins]))
  })

  output$fixPlot <-  renderPlot({fixPlot()})
  
  
  
  
  comp.text <- eventReactive(input$upN,{
    if(input$bins < 4){
      return(paste('<p align="justify"><font size="6"><br><b>You need to have at least 4 levels of complexity to see complexity and generalizability trade-off plot</b></font></p>',sep = ""))
    }
    if(input$bins >= 4){
      return(NULL)
    }
  })
  
  
  output$comp.text <- renderUI({
    HTML(comp.text())
  })
  
  
  
  
  
})
  
  













