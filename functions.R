library(ggplot2)
library(ggthemes)
library(ggiraph)
library(ReporteRs)

noise_helper <- function(noise.prop, df, col){
  if(noise.prop == 0){
    return(df)
  }
  if(noise.prop > 0){
    for(r in c(1:(dim(df)[1]))){
      val <- min(c((df[r,col]*(1-noise.prop) + runif(1)*noise.prop),1))
      df[r,col] <- val
    }
  }
  return(df)
  
}


make_prec2 <- function(df, columns) {
  w_data <- df
  n_rows <- dim(w_data)[1]
  for (r in 1:n_rows){
    for (c in columns){
      w_data[r,c] <- round(w_data[r,c]*100,digits = 0)
    }
  }
  
  for (r in 1:n_rows){
    for (c in columns){
      if(w_data[r,c] != 0){
        w_data[r,c] <- as.character(w_data[r,c])
      }
      
      if(w_data[r,c] == 0){
        w_data[r,c] <- ""
      }
    }
  }
  
  for (r in 1:n_rows){
    for (c in columns){
      if(w_data[r,c] != ""){
        w_data[r,c] <- paste((w_data[r,c]),"%",sep = "")
      }
    }
  }
  
  return(w_data)
}



#DATA GENERATION OBJECT
setClass("generated.data",
         slots = c(a = "numeric",
                   b = "numeric",
                   c = "numeric",
                   n.total = "numeric",
                   
                   noise.prop = "numeric",
                   
                   #Data generated w/o noise
                   data.m1.prob.r = "data.frame",
                   data.m2.prob.r = "data.frame",
                   data.m3.prob.r = "data.frame",
                   
                   data.comb.prob.r = "data.frame",
                   
                   data.m1.freq.r = "numeric",
                   data.m2.freq.r= "numeric",
                   data.m3.freq.r = "numeric",
                   
                   pram.estimates.m1.r = "numeric",
                   pram.estimates.m2.r = "numeric",
                   pram.estimates.m3.r = "numeric",
                   
                   
                   #Data generated w/ noise
                   data.m1.prob.n = "data.frame",
                   data.m2.prob.n = "data.frame",
                   data.m3.prob.n = "data.frame",
                   
                   data.comb.prob.n = "data.frame",
                   
                   data.m1.freq.n = "numeric",
                   data.m2.freq.n = "numeric",
                   data.m3.freq.n = "numeric",
                   
                   
                   pram.estimates.m1.n = "numeric",
                   pram.estimates.m2.n = "numeric",
                   pram.estimates.m3.n = "numeric"
                   
                   
                   )
         )





paper_data_gen <- function(a,
                            b,
                            c,
                            n.total,
                            noise.prop
                            ){
  
  
  return_object <- new("generated.data")
  
  return_object@n.total <- n.total
  return_object@a <- a
  return_object@b <- b
  return_object@c <- c
  
  #Defining the data generating functions
  m.1 <- function(param, t){
    (1 + t)^-param[1]
  }
  
  m.2 <- function(param, t){
    (param[2] + t)^-param[1]
  }
  
  m.3 <- function(param, t){
    (1 + param[2]*t)^-param[1]
  }
  
  
  #Timepoint and runs per sample
  n <- 50
  t <- seq(0.1, 8.1, by=2)
  
  
  
  
  #Generating frequencies
  yi.1.r <- sapply(t, 
                 function(x){
                   rbinom(1, n, m.1(a, x))})
  
  yi.2.r <- sapply(t,
                 function(x){
                   rbinom(1, n, m.2(c(a, b), x))
                 })
  
  yi.3.r <- sapply(t,
                 function(x){
                   rbinom(1, n, m.3(c(a, b), x))
                 })
  
  return_object@data.m1.freq.r <- yi.1.r
  return_object@data.m2.freq.r <- yi.2.r
  return_object@data.m3.freq.r <- yi.3.r
  
  

  
  
  
  ### Likelihoods
  loglike.1 <- function(param, data, t){
    prob <- m.1(param, t)
    prob <- c(prob, 1-prob)
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  
  loglike.2 <- function(param, data, t){
    prob  <- sapply(t, 
                    function(x){
                      m.2(param, x)})
    prob <- c(prob, 1-prob) 
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  
  loglike.3 <- function(param, data, t){
    prob <- m.3(param, t)
    prob <- c(prob, 1-prob) 
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  


  

  
  ## maximising likelihood
  max.1.r <-  nlm(loglike.1, data=c(yi.1.r, n-yi.1.r), t=t, p=c(0.8))
  
  max.2.r <- nlm(loglike.2, data=c(yi.2.r, n-yi.2.r), t=t, p=c(.8, 0.99))
  
  max.3.r <- nlm(loglike.3, data=c(yi.3.r, n-yi.3.r), t=t, p=c(0.8, 0.5))
  
  return_object@pram.estimates.m1.r <- max.1.r$estimate
  return_object@pram.estimates.m2.r <- max.2.r$estimate
  return_object@pram.estimates.m3.r <- max.3.r$estimate
  
  
  #Data from model 1
  data1.r <- data.frame(ID=seq_len(n.total),
                      time=rep(t, each=n.total/length(t)))
  data1.r$prob.m1 <- sapply(data1.r$time,
                          function(x){
                            m.1(max.1.r$estimate, x)})
  
  return_object@data.m1.prob.r <- data1.r
  data1.n <- noise_helper(df = data1.r,noise.prop = noise.prop, col = "prob.m1")
  return_object@data.m1.prob.n <- data1.n


  # #Data from model 2
  data2.r <- data.frame(ID=seq_len(n.total),
                     time=rep(t, each=n.total/length(t)))
  data2.r$prob.m2 <- sapply(data2.r$time,
                         function(x){
                           m.2(max.2.r$estimate, x)})

  return_object@data.m2.prob.r <- data2.r
  data2.n <- noise_helper(df = data2.r,noise.prop = noise.prop, col = "prob.m2")
  return_object@data.m2.prob.n <- data2.n


  #Data from model 3
  data3.r <- data.frame(ID=seq_len(n.total),
                     time=rep(t, each=n.total/length(t)))
  data3.r$prob.m3 <- sapply(data3.r$time,
                         function(x){
                           m.3(max.3.r$estimate, x)})
  
  return_object@data.m3.prob.r <- data3.r
  data3.n <- noise_helper(df = data3.r,noise.prop = noise.prop, col = "prob.m3")
  return_object@data.m3.prob.n <- data3.n


  #Data combined
  data.c.r <- data1.r
  data.c.r <- cbind(data.c.r,data2.r[,3],data3.r[,3])
  return_object@data.comb.prob.r <- data.c.r
  
  data.c.n <- data1.n
  data.c.n <- cbind(data.c.n,data2.n[,3],data3.n[,3])
  return_object@data.comb.prob.n <- data.c.n
  
  
  
  
  ## Estimating parameters from data
  #Generating frequencies
  yi.1.n <- sapply(return_object@data.m1.prob.n$prob.m1, function(x){rbinom(1, n, x)
                    
                  })
  
  yi.2.n <- sapply(return_object@data.m2.prob.n$prob.m2, function(x){rbinom(1, n, x)
    
                  })
  
  yi.3.n <- sapply(return_object@data.m3.prob.n$prob.m3, function(x){rbinom(1, n, x)
    
                  })
  
  return_object@data.m1.freq.n <- yi.1.n
  return_object@data.m2.freq.n <- yi.2.n
  return_object@data.m3.freq.n <- yi.3.n
  
  
  ## maximising likelihood
  max.1.n <- nlm(loglike.1, data=c(yi.1.n, n-yi.1.n), t=data1.n$time, p=c(0.5))
  
  max.2.n <- nlm(loglike.2, data=c(yi.2.n, n-yi.2.n), t=data2.n$time, p=c(0.9, 0.5))
  
  max.3.n <- nlm(loglike.3, data=c(yi.3.n, n-yi.3.n), t=data3.n$time, p=c(0.5, 0.5))
  
  
  return_object@pram.estimates.m1.n <- max.1.n$estimate
  return_object@pram.estimates.m2.n <- max.2.n$estimate
  return_object@pram.estimates.m3.n <- max.3.n$estimate
  
  
  

  
  return(return_object)
  
}


#test <- paper_data_gen(1,2,3,100,noise.prop = 0.1)





curve_generator <- function(param, model = 1){
  
  t = seq(from = 0.1,to = 10,by = 0.1)
                  
  return_results <- as.data.frame(matrix(rep(NA,times = 200),nrow = 100,ncol = 2))
  names(return_results) <- c("time","prob")
  
  m.1 <- function(param, t){
    return((1 + t)^-param[1])
  }
  
  m.2 <- function(param, t){
    return((param[2] + t)^-param[1])
  }
  
  m.3 <- function(param, t){
    return((1 + param[2]*t)^-param[1])
  }
  
  if(model == 1){
    for(i in c(1:100)){
      return_results[i,"time"] <- t[i]
      return_results[i,"prob"] <- m.1(param,t[i])
    }
  }
  
  if(model == 2){
    for(i in c(1:100)){
      return_results[i,"time"] <- t[i]
      return_results[i,"prob"] <- m.2(param,t[i])
    }
  }
  
  if(model == 3){
    for(i in c(1:100)){
      return_results[i,"time"] <- t[i]
      return_results[i,"prob"] <- m.3(param,t[i])
    }
  }

  
  return(return_results)
}
  


setOldClass("gg")
  
#DATA FITTING
setClass("fitted.prams",
         slots = c(data.1.m1.pram.n = "numeric",
                   data.1.m2.pram.n = "numeric",
                   data.1.m3.pram.n = "numeric",
                   
                   data.2.m1.pram.n = "numeric",
                   data.2.m2.pram.n = "numeric",
                   data.2.m3.pram.n = "numeric",
                   
                   data.3.m1.pram.n = "numeric",
                   data.3.m2.pram.n = "numeric",
                   data.3.m3.pram.n = "numeric",
                   
                   
                   data.1.m1.curve = "data.frame",
                   data.1.m2.curve = "data.frame",
                   data.1.m3.curve = "data.frame",
                   
                   data.2.m1.curve = "data.frame",
                   data.2.m2.curve = "data.frame",
                   data.2.m3.curve = "data.frame",
                   
                   data.3.m1.curve = "data.frame",
                   data.3.m2.curve = "data.frame",
                   data.3.m3.curve = "data.frame",
                   
                   data.1.comb.curve = "data.frame",
                   data.2.comb.curve = "data.frame",
                   data.3.comb.curve = "data.frame",
                   
                   data.1.comb.plot = "gg",
                   data.2.comb.plot = "gg",
                   data.3.comb.plot = "gg",
                   
                   data.1.m1.fit.n = "list",
                   data.1.m2.fit.n = "list",
                   data.1.m3.fir.n = "list",
                   
                   data.2.m1.fit.n = "list",
                   data.2.m2.fit.n = "list",
                   data.2.m3.fit.n = "list",
                   
                   data.3.m1.fit.n = "list",
                   data.3.m2.fit.n = "list",
                   data.3.m3.fit.n = "list",
                   
                   data.1.freq.comb = "data.frame",
                   data.1.freq.plot = "gg"
                   )
)


paper_data_fitter <- function(data.object){
  
  
  #Return object
  return_object <- new("fitted.prams")
  
  #Timepoint and runs per sample
  n <- 50
  t <- seq(0.1, 8.1, by=2)
  
  ### Likelihoods
  loglike.1 <- function(param, data, t){
    prob <- m.1(param, t)
    prob <- c(prob, 1-prob)
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  
  loglike.2 <- function(param, data, t){
    prob  <- sapply(t, 
                    function(x){
                      m.2(param, x)})
    prob <- c(prob, 1-prob) 
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  
  loglike.3 <- function(param, data, t){
    prob <- m.3(param, t)
    prob <- c(prob, 1-prob) 
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  
  #Defining the data generating functions
  m.1 <- function(param, t){
    (1 + t)^-param[1]
  }
  
  m.2 <- function(param, t){
    (param[2] + t)^-param[1]
  }
  
  m.3 <- function(param, t){
    (1 + param[2]*t)^-param[1]
  }
  
  
  #Loading data
  dat.1 <- data.object@data.m1.freq.n 
  dat.2 <- data.object@data.m2.freq.n
  dat.3 <- data.object@data.m3.freq.n
  
  #Estimating pramd
  return_object@data.1.m1.pram.n <- nlm(loglike.1, data=c(dat.1, n-dat.1), t=rep(t, each=length(dat.1)/length(t)), p=c(0.8))$estimate
  return_object@data.1.m2.pram.n <- nlm(loglike.2, data=c(dat.1, n-dat.1), t=rep(t, each=length(dat.1)/length(t)), p=c(0.8, 0.99))$estimate
  return_object@data.1.m3.pram.n <- nlm(loglike.3, data=c(dat.1, n-dat.1), t=rep(t, each=length(dat.1)/length(t)), p=c(0.8, 0.99))$estimate
  
  return_object@data.2.m1.pram.n <- nlm(loglike.1, data=c(dat.2, n-dat.2), t=rep(t, each=length(dat.2)/length(t)), p=c(0.8))$estimate
  return_object@data.2.m2.pram.n <- nlm(loglike.2, data=c(dat.2, n-dat.2), t=rep(t, each=length(dat.2)/length(t)), p=c(0.8, 0.99))$estimate
  return_object@data.2.m3.pram.n <- nlm(loglike.3, data=c(dat.2, n-dat.2), t=rep(t, each=length(dat.2)/length(t)), p=c(0.8, 0.99))$estimate
  
  return_object@data.3.m1.pram.n <- nlm(loglike.1, data=c(dat.3, n-dat.3), t=rep(t, each=length(dat.3)/length(t)), p=c(0.8))$estimate
  return_object@data.3.m2.pram.n <- nlm(loglike.2, data=c(dat.3, n-dat.3), t=rep(t, each=length(dat.3)/length(t)), p=c(0.8, 0.99))$estimate
  return_object@data.3.m3.pram.n <- nlm(loglike.3, data=c(dat.3, n-dat.3), t=rep(t, each=length(dat.3)/length(t)), p=c(0.8, 0.99))$estimate
  
  
  
  #Generating curves
  return_object@data.1.m1.curve = "data.frame" <- curve_generator(param = return_object@data.1.m1.pram.n, model = 1)
  return_object@data.1.m1.curve[,"model"] <- rep("Model 1", times = dim(return_object@data.1.m1.curve)[1])
  return_object@data.1.m2.curve = "data.frame" <- curve_generator(param = return_object@data.1.m2.pram.n, model = 2)
  return_object@data.1.m2.curve[,"model"] <- rep("Model 2", times = dim(return_object@data.1.m2.curve)[1])
  return_object@data.1.m3.curve = "data.frame" <- curve_generator(param = return_object@data.1.m3.pram.n, model = 3)
  return_object@data.1.m3.curve[,"model"] <- rep("Model 3", times = dim(return_object@data.1.m3.curve)[1])
  
  return_object@data.2.m1.curve = "data.frame" <- curve_generator(param = return_object@data.2.m1.pram.n, model = 1)
  return_object@data.2.m1.curve[,"model"] <- rep("Model 1", times = dim(return_object@data.2.m1.curve)[1])
  return_object@data.2.m2.curve = "data.frame" <- curve_generator(param = return_object@data.2.m2.pram.n, model = 2)
  return_object@data.2.m2.curve[,"model"] <- rep("Model 2", times = dim(return_object@data.2.m2.curve)[1])
  return_object@data.2.m3.curve = "data.frame" <- curve_generator(param = return_object@data.2.m3.pram.n, model = 3)
  return_object@data.2.m3.curve[,"model"] <- rep("Model 3", times = dim(return_object@data.2.m3.curve)[1])
  
  return_object@data.3.m1.curve = "data.frame" <- curve_generator(param = return_object@data.3.m1.pram.n, model = 1)
  return_object@data.3.m1.curve[,"model"] <- rep("Model 1", times = dim(return_object@data.3.m1.curve)[1])
  return_object@data.3.m2.curve = "data.frame" <- curve_generator(param = return_object@data.3.m2.pram.n, model = 2)
  return_object@data.3.m2.curve[,"model"] <- rep("Model 2", times = dim(return_object@data.3.m2.curve)[1])
  return_object@data.3.m3.curve = "data.frame" <- curve_generator(param = return_object@data.3.m3.pram.n, model = 3)
  return_object@data.3.m3.curve[,"model"] <- rep("Model 3", times = dim(return_object@data.3.m3.curve)[1])
  
  return_object@data.1.comb.curve  <- rbind(return_object@data.1.m1.curve,return_object@data.1.m2.curve,return_object@data.1.m3.curve)
  return_object@data.2.comb.curve  <- rbind(return_object@data.2.m1.curve,return_object@data.2.m2.curve,return_object@data.2.m3.curve)
  return_object@data.3.comb.curve  <- rbind(return_object@data.3.m1.curve,return_object@data.3.m2.curve,return_object@data.3.m3.curve)
  return_object@data.1.comb.curve[,"ID"] <- c(1:dim( return_object@data.1.comb.curve)[1])
  return_object@data.2.comb.curve[,"ID"] <- c(1:dim( return_object@data.1.comb.curve)[1])
  return_object@data.3.comb.curve[,"ID"] <- c(1:dim( return_object@data.1.comb.curve)[1])
  
  
  
  #Generating frequencies
  
  freq.template <- as.data.frame(matrix(rep(NA, times = 5*3*3),nrow = 15,ncol = 3))
  names(freq.template) <- c("time","successes","model")
  freq.template[,"time"] <- rep(t, times = 3)
  freq.template[,"model"] <- rep(c("Model 1","Model 1","Model 1","Model 1","Model 1",
                                   "Model 2","Model 2","Model 2","Model 2","Model 2",
                                   "Model 3","Model 3","Model 3","Model 3","Model 3"), times = 1)
  
  freq.template[,"successes"] <- c(dat.1[seq(1,length(dat.1), by=(data.object@n.total/5))]/50,
                                   dat.2[seq(1,length(dat.2), by=(data.object@n.total/5))]/50,
                                   dat.3[seq(1,length(dat.3), by=(data.object@n.total/5))]/50)
  
  return_object@data.1.freq.comb <- freq.template
  
  
  
  return_object@data.1.comb.curve[,"successes"] <- rep(NA, times = dim(return_object@data.1.comb.curve)[1])
  temp <- freq.template[freq.template[,"model"] == "Model 1",]
  for(i in c(1:dim(temp)[1])){
    return_object@data.1.comb.curve[(return_object@data.1.comb.curve[,"time"] == temp[i,"time"] & 
                                    return_object@data.1.comb.curve[,"model"] == temp[i,"model"]),"successes"] <- temp[i,"successes"]
  }
  
  
  return_object@data.2.comb.curve[,"successes"] <- rep(NA, times = dim(return_object@data.2.comb.curve)[1])
  temp <- freq.template[freq.template[,"model"] == "Model 2",]
  for(i in c(1:dim(temp)[1])){
    return_object@data.2.comb.curve[(return_object@data.2.comb.curve[,"time"] == temp[i,"time"] & 
                                       return_object@data.2.comb.curve[,"model"] == temp[i,"model"]),"successes"] <- temp[i,"successes"]
  }
  
  
  return_object@data.3.comb.curve[,"successes"] <- rep(NA, times = dim(return_object@data.3.comb.curve)[1])
  temp <- freq.template[freq.template[,"model"] == "Model 3",]
  for(i in c(1:dim(temp)[1])){
    return_object@data.3.comb.curve[(return_object@data.3.comb.curve[,"time"] == temp[i,"time"] & 
                                       return_object@data.3.comb.curve[,"model"] == temp[i,"model"]),"successes"] <- temp[i,"successes"]
  }
  
  
  
  
  return_object@data.1.freq.plot <- ggplot(return_object@data.1.freq.comb, aes(x=time, y=successes, colour=model, group=model))+
    geom_point(size=5, aes(shape = model), position = position_jitter(w = 0.1, h = 0)) + 
    theme(legend.position="bottom") + 
    scale_y_continuous(labels = scales::percent,limits = c(0,1),breaks = seq(0,1,0.2))+
    
    theme(axis.title.y = element_text(colour="black", size=4),
          axis.text.y  = element_text(vjust=0, size=11))+
    theme(axis.title.y = element_blank())+
    
    ggtitle("Realized successes by three models") +
    theme(plot.title = element_text(lineheight=.8, face="bold"))+
    theme(plot.title = element_text(hjust = 1, vjust=2.5))+
    
    theme_fivethirtyeight(base_size= 18)
  
  
  
  
  return_object@data.1.comb.plot <-   ggplot(return_object@data.1.comb.curve, aes(x=time, y=prob, colour=model, group=model))+
                                              geom_line(size=1.5) + 
                                              theme(legend.position="bottom") + 
                                              scale_y_continuous(labels = scales::percent,limits = c(0,1),breaks = seq(0,1,0.2))+
                                              geom_point(aes(y= successes, x = time),size=5, shape = 17) + 
                                              
                                              theme(axis.title.y = element_text(colour="black", size=4),
                                                    axis.text.y  = element_text(vjust=0, size=11))+
                                              theme(axis.title.y = element_blank())+
                                              
                                              ggtitle("Probabilities of success based on data generated by model 1") +
                                              theme(plot.title = element_text(lineheight=.8, face="bold"))+
                                              theme(plot.title = element_text(hjust = 1, vjust=2.5))+
                                              
                                              theme_fivethirtyeight(base_size= 18) 
  
  return_object@data.2.comb.plot <-   ggplot(return_object@data.2.comb.curve, aes(x=time, y=prob, colour=model, group=model))+
                                              geom_line(size=1.5) + 
                                              theme(legend.position="bottom") + 
                                              scale_y_continuous(labels = scales::percent,limits = c(0,1),breaks = seq(0,1,0.2))+
                                              geom_point(aes(y= successes, x = time),size=5, shape = 17) + 
                                              
                                              theme(axis.title.y = element_text(colour="black", size=4),
                                                    axis.text.y  = element_text(vjust=0, size=11))+
                                              theme(axis.title.y = element_blank())+
                                              
                                              ggtitle("Probabilities of success based on data generated by model 2") +
                                              theme(plot.title = element_text(lineheight=.8, face="bold"))+
                                              theme(plot.title = element_text(hjust = 1, vjust=2.5))+
                                              
                                              theme_fivethirtyeight(base_size= 18) 
  
  return_object@data.3.comb.plot <-   ggplot(return_object@data.3.comb.curve, aes(x=time, y=prob, colour=model, group=model))+
                                              geom_line(size=1.5) + 
                                              theme(legend.position="bottom") + 
                                              scale_y_continuous(labels = scales::percent,limits = c(0,1),breaks = seq(0,1,0.2))+
                                              geom_point(aes(y= successes, x = time),size=5, shape = 17) + 
                                              
                                              theme(axis.title.y = element_text(colour="black", size=4),
                                                    axis.text.y  = element_text(vjust=0, size=11))+
                                              theme(axis.title.y = element_blank())+
                                              
                                              ggtitle("Probabilities of success based on data generated by model 3") +
                                              theme(plot.title = element_text(lineheight=.8, face="bold"))+
                                              theme(plot.title = element_text(hjust = 1, vjust=2.5))+
                                              
                                              theme_fivethirtyeight(base_size= 18)


  return(return_object)

}


#test22 <- paper_data_fitter(test)










#DATA FITTING
setClass("simulator.aicbic",
         slots = c(n.runs = "numeric",
                   raw.results.aic = "data.frame",
                   raw.results.bic = "data.frame",
                   raw.results.rmse = "data.frame",
                   raw.results.pvaf = "data.frame",
                   agregated.results.aic = "data.frame",
                   agregated.results.bic = "data.frame",
                   agregated.results.rmse = "data.frame",
                   agregated.results.pvaf = "data.frame"
                   
         )
)



simulator_aicbic <- function(n.sims, 
                             n.total,
                             noise.prop,
                             a,
                             b,
                             c
                             ){
  
  return_object <- new("simulator.aicbic")
  
  #Timepoint and runs per sample
  n <- 50
  t <- seq(0.1, 8.1, by=2)
  
  
  #Defining the data generating functions
  m.1 <- function(param, t){
    (1 + t)^-param[1]
  }
  
  m.2 <- function(param, t){
    (param[2] + t)^-param[1]
  }
  
  m.3 <- function(param, t){
    (1 + param[2]*t)^-param[1]
  }
  
  
  #Generating frequencies
  yi.1.r <- sapply(t, 
                   function(x){
                     rbinom(1, n, m.1(a, x))})
  
  yi.2.r <- sapply(t,
                   function(x){
                     rbinom(1, n, m.2(c(a, b), x))
                   })
  
  yi.3.r <- sapply(t,
                   function(x){
                     rbinom(1, n, m.3(c(a, b), x))
                   })
  
  
  
  ### Likelihoods
  loglike.1 <- function(param, data, t){
    prob <- m.1(param, t)
    prob <- c(prob, 1-prob)
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  
  loglike.2 <- function(param, data, t){
    prob  <- sapply(t, 
                    function(x){
                      m.2(param, x)})
    prob <- c(prob, 1-prob) 
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  
  loglike.3 <- function(param, data, t){
    prob <- m.3(param, t)
    prob <- c(prob, 1-prob) 
    -sum(sapply(data, 
                function(x)
                {lchoose(n, x)}) + data*log(prob))
  }
  
  
  
  
  ## maximising likelihood
  max.1.r <-  nlm(loglike.1, data=c(yi.1.r, n-yi.1.r), t=t, p=c(0.8))
  
  max.2.r <- nlm(loglike.2, data=c(yi.2.r, n-yi.2.r), t=t, p=c(.8, 0.99))
  
  max.3.r <- nlm(loglike.3, data=c(yi.3.r, n-yi.3.r), t=t, p=c(0.8, 0.5))
  
  
  data1.r <- data.frame(ID=seq_len(n.total),
                        time=rep(t, each=n.total/length(t)))
  data1.r$prob.m1 <- sapply(data1.r$time,
                            function(x){
                              m.1(max.1.r$estimate, x)})
  
  data2.r <- data.frame(ID=seq_len(n.total),
                        time=rep(t, each=n.total/length(t)))
  data2.r$prob.m2 <- sapply(data2.r$time,
                            function(x){
                              m.2(max.2.r$estimate, x)})
  
  data3.r <- data.frame(ID=seq_len(n.total),
                        time=rep(t, each=n.total/length(t)))
  data3.r$prob.m3 <- sapply(data3.r$time,
                            function(x){
                              m.3(max.3.r$estimate, x)})
  
  
  
  bic.new <- function(loglike, nparam, sample){
    -2*loglike + nparam*log(sample)
  }
  
  aic.new <- function(loglike, nparam){
    -2*loglike + 2*nparam
  }
  
  rmse <- function(obs, pred){
    sqrt(sum((obs - pred)^2) / length(obs))
    ## or:
    # sqrt(mean((obs - pred)^2))
  }
  
  pvaf <- function(obs, pred){
    1 - (sum((obs - pred)^2) / sum((obs - mean(obs))^2))
  }
  
  
  raw.results.aic <- as.data.frame(matrix(rep(NA, times = 10*n.sims),nrow = n.sims,ncol = 10))
  names(raw.results.aic) <- c("Simulation run", 
                          "Data from M1 fitted by M1 (AIC)",
                          "Data from M1 fitted by M2 (AIC)",
                          "Data from M1 fitted by M3 (AIC)",
                          
                          "Data from M2 fitted by M1 (AIC)",
                          "Data from M2 fitted by M2 (AIC)",
                          "Data from M2 fitted by M3 (AIC)",
                          
                          "Data from M3 fitted by M1 (AIC)",
                          "Data from M3 fitted by M2 (AIC)",
                          "Data from M3 fitted by M3 (AIC)"
                          )
  
  raw.results.bic <- as.data.frame(matrix(rep(NA, times = 10*n.sims),nrow = n.sims,ncol = 10))
  names(raw.results.bic) <- c("Simulation run", 
                              "Data from M1 fitted by M1 (BIC)",
                              "Data from M1 fitted by M2 (BIC)",
                              "Data from M1 fitted by M3 (BIC)",
                              
                              "Data from M2 fitted by M1 (BIC)",
                              "Data from M2 fitted by M2 (BIC)",
                              "Data from M2 fitted by M3 (BIC)",
                              
                              "Data from M3 fitted by M1 (BIC)",
                              "Data from M3 fitted by M2 (BIC)",
                              "Data from M3 fitted by M3 (BIC)"
  )
  
  raw.results.rmse <- as.data.frame(matrix(rep(NA, times = 10*n.sims),nrow = n.sims,ncol = 10))
  names(raw.results.rmse) <- c("Simulation run", 
                              "Data from M1 fitted by M1 (RMSE)",
                              "Data from M1 fitted by M2 (RMSE)",
                              "Data from M1 fitted by M3 (RMSE)",
                              
                              "Data from M2 fitted by M1 (RMSE)",
                              "Data from M2 fitted by M2 (RMSE)",
                              "Data from M2 fitted by M3 (RMSE)",
                              
                              "Data from M3 fitted by M1 (RMSE)",
                              "Data from M3 fitted by M2 (RMSE)",
                              "Data from M3 fitted by M3 (RMSE)"
  )
  
  raw.results.pvaf <- as.data.frame(matrix(rep(NA, times = 10*n.sims),nrow = n.sims,ncol = 10))
  names(raw.results.pvaf) <- c("Simulation run", 
                               "Data from M1 fitted by M1 (PVAF)",
                               "Data from M1 fitted by M2 (PVAF)",
                               "Data from M1 fitted by M3 (PVAF)",
                               
                               "Data from M2 fitted by M1 (PVAF)",
                               "Data from M2 fitted by M2 (PVAF)",
                               "Data from M2 fitted by M3 (PVAF)",
                               
                               "Data from M3 fitted by M1 (PVAF)",
                               "Data from M3 fitted by M2 (PVAF)",
                               "Data from M3 fitted by M3 (PVAF)"
  )
    
  ## SIMULATION GOES HERE
  for(i in c(1:n.sims)){
    
    # #Data from model 1
    data1.n <- noise_helper(df = data1.r,noise.prop = noise.prop, col = "prob.m1")
    # #Data from model 2
    data2.n <- noise_helper(df = data2.r,noise.prop = noise.prop, col = "prob.m2")
    #Data from model 3
    data3.n <- noise_helper(df = data3.r,noise.prop = noise.prop, col = "prob.m3")

  
    
    ## Estimating parameters from data
    #Generating frequencies
    yi.1.n <- sapply(data1.n$prob.m1, function(x){rbinom(1, n, x)
      
    })
    max.1.m1.n <- nlm(loglike.1, data=c(yi.1.n, n-yi.1.n), t=data1.n$time, p=c(0.8))
    max.1.m2.n <- nlm(loglike.2, data=c(yi.1.n, n-yi.1.n), t=data1.n$time, p=c(.8, 0.99))
    max.1.m3.n <- nlm(loglike.3, data=c(yi.1.n, n-yi.1.n), t=data1.n$time, p=c(0.8, 0.5))
    
    #PRED 1-1
    max.1.m1.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data1.n$time)))
    max.1.m1.n.pred$prob.m1  <- sapply(max.1.m1.n.pred$time,
                              function(x){
                                m.1(max.1.m1.n$estimate, x)})
    max.1.m1.n.pred <- sapply(max.1.m1.n.pred$prob.m1, function(x){rbinom(1, n, x)})

    #PRED 1-2
    max.1.m2.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data2.n$time)))
    max.1.m2.n.pred$prob.m2  <- sapply(max.1.m2.n.pred$time,
                                       function(x){
                                         m.2(max.1.m2.n$estimate, x)})
    max.1.m2.n.pred <- sapply(max.1.m2.n.pred$prob.m2, function(x){rbinom(1, n, x)})

    #PRED 1-3
    max.1.m3.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data3.n$time)))
    max.1.m3.n.pred$prob.m3  <- sapply(max.1.m3.n.pred$time,
                                       function(x){
                                         m.3(max.1.m3.n$estimate, x)})
    max.1.m3.n.pred <- sapply(max.1.m3.n.pred$prob.m3, function(x){rbinom(1, n, x)})

    
    raw.results.aic[i,"Simulation run"] <- i  
    raw.results.aic[i,"Data from M1 fitted by M1 (AIC)"] <-  aic.new(loglike = loglike.1(max.1.m1.n$estimate,
                                                                                     data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                 nparam = length(max.1.m1.n$estimate))
    
    raw.results.aic[i,"Data from M1 fitted by M2 (AIC)"] <-  aic.new(loglike = loglike.2(max.1.m2.n$estimate,
                                                                                     data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                 nparam = length(max.1.m2.n$estimate))
    
    raw.results.aic[i,"Data from M1 fitted by M3 (AIC)"] <-  aic.new(loglike = loglike.3(max.1.m3.n$estimate,
                                                                                     data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                 nparam = length(max.1.m3.n$estimate))
    
    
    raw.results.bic[i,"Simulation run"] <- i  
    raw.results.bic[i,"Data from M1 fitted by M1 (BIC)"] <-  bic.new(loglike = loglike.1(max.1.m1.n$estimate,
                                                                                         data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                     nparam = length(max.1.m1.n$estimate),sample = length(yi.1.n))
    
    raw.results.bic[i,"Data from M1 fitted by M2 (BIC)"] <-  bic.new(loglike = loglike.2(max.1.m2.n$estimate,
                                                                                         data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                     nparam = length(max.1.m2.n$estimate),sample = length(yi.1.n))
    
    raw.results.bic[i,"Data from M1 fitted by M3 (BIC)"] <-  bic.new(loglike = loglike.3(max.1.m3.n$estimate,
                                                                                         data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                     nparam = length(max.1.m3.n$estimate),sample = length(yi.1.n))
    
    
    raw.results.rmse[i,"Simulation run"] <- i
    raw.results.rmse[i,"Data from M1 fitted by M1 (RMSE)"] <-  rmse(obs = yi.1.n,max.1.m1.n.pred)

    raw.results.rmse[i,"Data from M1 fitted by M2 (RMSE)"] <-  rmse(obs = yi.1.n,max.1.m2.n.pred)

    raw.results.rmse[i,"Data from M1 fitted by M3 (RMSE)"] <-  rmse(obs = yi.1.n,max.1.m3.n.pred)
    
    
    raw.results.pvaf[i,"Simulation run"] <- i
    raw.results.pvaf[i,"Data from M1 fitted by M1 (PVAF)"] <-  pvaf(obs = yi.1.n,max.1.m1.n.pred)
    
    raw.results.pvaf[i,"Data from M1 fitted by M2 (PVAF)"] <-  pvaf(obs = yi.1.n,max.1.m2.n.pred)
    
    raw.results.pvaf[i,"Data from M1 fitted by M3 (PVAF)"] <-  pvaf(obs = yi.1.n,max.1.m3.n.pred)
      
      
      
      
    yi.2.n <- sapply(data2.n$prob.m2, function(x){rbinom(1, n, x)})
    
    max.2.m1.n <-  nlm(loglike.1, data=c(yi.2.n, n-yi.2.n), t=data2.n$time, p=c(0.8))
    max.2.m2.n <- nlm(loglike.2, data=c(yi.2.n, n-yi.2.n), t=data2.n$time, p=c(.8, 0.99))
    max.2.m3.n <- nlm(loglike.3, data=c(yi.2.n, n-yi.2.n), t=data2.n$time, p=c(0.8, 0.5))
    
    #PRED 2-1
    max.2.m1.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data1.n$time)))
    max.2.m1.n.pred$prob.m1  <- sapply(max.2.m1.n.pred$time,
                                       function(x){
                                         m.1(max.2.m1.n$estimate, x)})
    max.2.m1.n.pred <- sapply(max.2.m1.n.pred$prob.m1, function(x){rbinom(1, n, x)})
    
    # #PRED 2-2
    max.2.m2.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data2.n$time)))
    max.2.m2.n.pred$prob.m2  <- sapply(max.2.m2.n.pred$time,
                                       function(x){
                                         m.2(max.2.m2.n$estimate, x)})
    max.2.m2.n.pred <- sapply(max.2.m2.n.pred$prob.m2, function(x){rbinom(1, n, x)})
    
    #PRED 2-3
    max.2.m3.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data3.n$time)))
    max.2.m3.n.pred$prob.m3  <- sapply(max.2.m3.n.pred$time,
                                       function(x){
                                         m.3(max.2.m3.n$estimate, x)})
    max.2.m3.n.pred <- sapply(max.2.m3.n.pred$prob.m3, function(x){rbinom(1, n, x)})
    
    
    raw.results.aic[i,"Data from M2 fitted by M1 (AIC)"] <-  aic.new(loglike = loglike.1(max.2.m1.n$estimate,
                                                                                     data=c(yi.2.n, n-yi.2.n), t = data2.n$time),
                                                                 nparam = length(max.2.m1.n$estimate))
    
    raw.results.aic[i,"Data from M2 fitted by M2 (AIC)"] <-  aic.new(loglike = loglike.2(max.2.m2.n$estimate,
                                                                                     data=c(yi.2.n, n-yi.2.n), t = data2.n$time),
                                                                 nparam = length(max.2.m2.n$estimate))
    
    raw.results.aic[i,"Data from M2 fitted by M3 (AIC)"] <-  aic.new(loglike = loglike.3(max.2.m3.n$estimate,
                                                                                     data=c(yi.2.n, n-yi.2.n), t = data2.n$time),
                                                                 nparam = length(max.2.m3.n$estimate))
    
    
    raw.results.bic[i,"Data from M2 fitted by M1 (BIC)"] <-  bic.new(loglike = loglike.1(max.2.m1.n$estimate,
                                                                                         data=c(yi.2.n, n-yi.2.n), t = data2.n$time),
                                                                     nparam = length(max.2.m1.n$estimate),sample = length(yi.2.n))
    
    raw.results.bic[i,"Data from M2 fitted by M2 (BIC)"] <-  bic.new(loglike = loglike.2(max.2.m2.n$estimate,
                                                                                         data=c(yi.2.n, n-yi.2.n), t = data2.n$time),
                                                                     nparam = length(max.2.m2.n$estimate),sample = length(yi.2.n))
    
    raw.results.bic[i,"Data from M2 fitted by M3 (BIC)"] <-  bic.new(loglike = loglike.3(max.2.m3.n$estimate,
                                                                                         data=c(yi.2.n, n-yi.2.n), t = data2.n$time),
                                                                     nparam = length(max.2.m3.n$estimate),sample = length(yi.2.n))
    
    
    raw.results.rmse[i,"Data from M2 fitted by M1 (RMSE)"] <-  rmse(obs = yi.2.n,max.2.m1.n.pred)
    
    raw.results.rmse[i,"Data from M2 fitted by M2 (RMSE)"] <-  rmse(obs = yi.2.n,max.2.m2.n.pred)
    
    raw.results.rmse[i,"Data from M2 fitted by M3 (RMSE)"] <-  rmse(obs = yi.2.n,max.2.m3.n.pred)
    
    
    raw.results.pvaf[i,"Data from M2 fitted by M1 (PVAF)"] <-  pvaf(obs = yi.2.n,max.2.m1.n.pred)
    
    raw.results.pvaf[i,"Data from M2 fitted by M2 (PVAF)"] <-  pvaf(obs = yi.2.n,max.2.m2.n.pred)
    
    raw.results.pvaf[i,"Data from M2 fitted by M3 (PVAF)"] <-  pvaf(obs = yi.2.n,max.2.m3.n.pred)
    
    yi.3.n <- sapply(data3.n$prob.m3, function(x){rbinom(1, n, x)    
    })
    
    max.3.m1.n <-  nlm(loglike.1, data=c(yi.3.n, n-yi.3.n), t=data3.n$time, p=c(0.8))
    max.3.m2.n <- nlm(loglike.2, data=c(yi.3.n, n-yi.3.n), t=data3.n$time, p=c(.8, 0.99))
    max.3.m3.n <- nlm(loglike.3, data=c(yi.3.n, n-yi.3.n), t=data3.n$time, p=c(0.8, 0.5))
    
    #PRED 3-1
    max.3.m1.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data1.n$time)))
    max.3.m1.n.pred$prob.m1  <- sapply(max.3.m1.n.pred$time,
                                       function(x){
                                         m.1(max.3.m1.n$estimate, x)})
    max.3.m1.n.pred <- sapply(max.3.m1.n.pred$prob.m1, function(x){rbinom(1, n, x)})
    
    #PRED 3-2
    max.3.m2.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data2.n$time)))
    max.3.m2.n.pred$prob.m2  <- sapply(max.3.m2.n.pred$time,
                                       function(x){
                                         m.2(max.3.m2.n$estimate, x)})
    max.3.m2.n.pred <- sapply(max.3.m2.n.pred$prob.m2, function(x){rbinom(1, n, x)})
    
    #PRED 3-3
    max.3.m3.n.pred <- data.frame(ID=seq_len(n.total),
                                  time=rep(t, each=n.total/length(data3.n$time)))
    max.3.m3.n.pred$prob.m3  <- sapply(max.3.m3.n.pred$time,
                                       function(x){
                                         m.3(max.3.m3.n$estimate, x)})
    max.3.m3.n.pred <- sapply(max.3.m3.n.pred$prob.m3, function(x){rbinom(1, n, x)})

    
    raw.results.aic[i,"Data from M3 fitted by M1 (AIC)"] <-  aic.new(loglike = loglike.1(max.3.m1.n$estimate,
                                                                                     data=c(yi.3.n, n-yi.3.n), t = data3.n$time),
                                                                 nparam = length(max.3.m1.n$estimate))
    
    raw.results.aic[i,"Data from M3 fitted by M2 (AIC)"] <-  aic.new(loglike = loglike.2(max.3.m2.n$estimate,
                                                                                     data=c(yi.3.n, n-yi.3.n), t = data3.n$time),
                                                                 nparam = length(max.3.m2.n$estimate))
    
    raw.results.aic[i,"Data from M3 fitted by M3 (AIC)"] <-  aic.new(loglike = loglike.3(max.3.m3.n$estimate,
                                                                                     data=c(yi.3.n, n-yi.3.n), t = data3.n$time),
                                                                 nparam = length(max.3.m3.n$estimate))
    
    
    raw.results.bic[i,"Data from M3 fitted by M1 (BIC)"] <-  bic.new(loglike = loglike.1(max.3.m1.n$estimate,
                                                                                         data=c(yi.3.n, n-yi.3.n), t = data3.n$time),
                                                                     nparam = length(max.3.m1.n$estimate),sample = length(yi.3.n))
    
    raw.results.bic[i,"Data from M3 fitted by M2 (BIC)"] <-  bic.new(loglike = loglike.2(max.3.m2.n$estimate,
                                                                                         data=c(yi.3.n, n-yi.3.n), t = data3.n$time),
                                                                     nparam = length(max.3.m2.n$estimate),sample = length(yi.3.n))
    
    raw.results.bic[i,"Data from M3 fitted by M3 (BIC)"] <-  bic.new(loglike = loglike.3(max.3.m3.n$estimate,
                                                                                         data=c(yi.3.n, n-yi.3.n), t = data3.n$time),
                                                                     nparam = length(max.3.m3.n$estimate),sample = length(yi.3.n))
    
    
    raw.results.rmse[i,"Data from M3 fitted by M1 (RMSE)"] <-  rmse(obs = yi.3.n,max.3.m1.n.pred)
    
    raw.results.rmse[i,"Data from M3 fitted by M2 (RMSE)"] <-  rmse(obs = yi.3.n,max.3.m2.n.pred)
    
    raw.results.rmse[i,"Data from M3 fitted by M3 (RMSE)"] <-  rmse(obs = yi.3.n,max.3.m3.n.pred)
    
    
    raw.results.pvaf[i,"Data from M3 fitted by M1 (PVAF)"] <-  pvaf(obs = yi.3.n,max.3.m1.n.pred)
    
    raw.results.pvaf[i,"Data from M3 fitted by M2 (PVAF)"] <-  pvaf(obs = yi.3.n,max.3.m2.n.pred)
    
    raw.results.pvaf[i,"Data from M3 fitted by M3 (PVAF)"] <-  pvaf(obs = yi.3.n,max.3.m3.n.pred)
    
  }

  return_object@raw.results.aic <- raw.results.aic
  return_object@raw.results.bic <- raw.results.bic
  return_object@raw.results.rmse <- raw.results.rmse
  return_object@raw.results.pvaf <- raw.results.pvaf
  
  
    
  return(return_object)
  

}





#test33 <- simulator_aicbic(n.sims = 10,n.total = 100,noise.prop = 0.1,a = 1,b = 2,c = 3)



data_table_output <- function(df,
                              col.width.1 = 2.5,
                              col.width.o = 3.5,
                              colors.set = c("white","#FFBA07")){
  
  
  rows <- dim(df)[1]
  
  
  word.table <- FlexTable(df)
  word.table <- setFlexTableWidths(word.table, widths = c(col.width.1,rep(col.width.o,9)))
  
  word.table[, 1:10, to = 'header'] = textProperties(color = "black", 
                                                     font.size = 14, 
                                                     font.weight = "bold", 
                                                     font.family = "Calibri" )
  word.table[, 1:10, to = 'header'] = cellProperties(vertical.align = "bottom" )
  word.table[, 1:10, to = 'header',] = parProperties( text.align = 'center' )
  
  word.table[, 2:10, to = 'header'] = cellProperties(vertical.align = "bottom", text.direction = "btlr" )
  
  word.table[, 1:01, to = 'header',] = parProperties( text.align = 'center' )
  
  word.table[, 1:10, to = 'header', side = "left"] <- borderProperties(style = 'none' )
  word.table[, 1:10, to = 'header', side = "right"] <- borderProperties(style = 'none' )
  word.table[, 1:10, to = 'header', side = 'bottom'] <- borderProperties( width=1, style = 'solid' )
  word.table[, 1,to = "header",  side = "left"] <- borderProperties(style = 'none',width=1)
  word.table[, 10,to = "header",  side = "left"] <- borderProperties(style = 'none',width=1)
  
  
  
  word.table[,1:10, side = "left"] <- borderProperties(style = 'none' )
  word.table[,1:10, side = "right"] <- borderProperties(style = 'none' )
  word.table[,1:10, side = 'bottom'] <- borderProperties(style = "solid", width=1)
  
  word.table[, 2:10,] = parProperties( text.align = 'center' )
  word.table[, 1,] = textProperties(color = "black", 
                                    font.size = 14, 
                                    font.family = "Calibri" )
  word.table[, 1:10,] = textProperties(color = "black", 
                                       font.size = 14, 
                                       font.family = "Calibri" )
  
  
  
  
  
  
  
  
  
  
  uni.q <- unlist(df[,2:10])
  uni.q <- uni.q[is.na(uni.q) == FALSE]
  uni.q <- unique(uni.q)
  uni.q <- uni.q[order(uni.q,na.last = TRUE,decreasing = TRUE)]
  
  colfunc <- colorRampPalette(colors.set)
  pallette <- colfunc(length(uni.q))
  pallette <- data.frame(Value=uni.q,Colour=pallette)
  pallette[,2] <- as.character(pallette[,2])
  
  
  for(r in c(1:rows)){
    for(c in c(2:10)){
      
      colour <- pallette[pallette[,"Value"] == df[r,c],"Colour"]
      word.table = setFlexTableBackgroundColors(word.table, i = r, j = c,colors = colour)
    }
  }
  
  
  x <- as.html(word.table)
  return(div(HTML(x), class = "shiny-html-output"))
  
  
}







  


################################################################################
################################################################################
## 1.- Muestra de Training y Muestra de Testing
################################################################################
################################################################################

rmse.j <- function(error){
  sqrt(mean(error^2))
}

multiGen = function(nIteration=100, nSample = 150,xModels){
  ml.out <- list()
  pred.lm.testing <- list()
  trainlist <- list()
  testinglist <- list()
  
  res = lapply(seq_len(nIteration),function(i){
    xx.train <- dataGen(nSample,plot=F,generator = 4)
    xx.testing <- dataGen(nSample,plot=F,generator = 4)
    n.all = nrow ( xx.train )
    trainingv <- c()
    testingv <- c()
    dftv <- c()
    valList = lapply(xModels,function(mod){
      model.out = lm(mod, data=xx.train )
      c(rmse.j(xx.train$y-fitted(model.out)),
        rmse.j(xx.testing$y-predict( model.out, newdata=xx.testing) )
      )
    })
    list(
      "train" = sapply(valList,function(x)x[1]),
      "test" = sapply(valList,function(x)x[2])
    )
    
    
  })
  rmseTr <- do.call(rbind,lapply(res,function(x)x$train))
  rmseTs <- do.call(rbind,lapply(res,function(x)x$test))
  dftv = lapply(xModels,function(mod){
    xx.train <- dataGen(nSample,plot=F,generator = 4)
    summary(lm(mod, data=xx.train ))$df[2]
  })
  
  ##PLOT
  
  plot(dftv,rmseTr[1,],type = "n", xlim = rev(range(dftv)), ylim=c(min(rmseTr),max(rmseTr)*1.2),
       xlab = "Model Complexity (DF)",ylab="Prediction Error",main ="Generalization and model complexity" )
  apply(rmseTr,1, function(x){
    lines(dftv,x,col=rgb(0.3,0.3,1,0.1),lwd=2)
  })
  
  apply(rmseTs,1, function(x){
    lines(dftv,x,col=rgb(1,0.3,0.3,0.1),lwd=2)
  })
  
  averagermseTr <- apply(rmseTr,2,mean)
  averagermseTs <- apply(rmseTs,2,mean)
  
  lines(smooth.spline(unlist(dftv),averagermseTr,spar = 0.25),col = "darkblue", lwd = 3)
  lines(smooth.spline(unlist(dftv),averagermseTs,spar = 0.25),col = "darkred", lwd = 3)
  
  # str(summary(model.out))
  
}
# multiGen(100,100,xModels[1:26])
##########################################################################################
########### J A V A L  C R O S S   V A L I D A T I O N   A L G O R I T H M ###############
##########################################################################################

jaVal = function(formula, data, folds=5, repetitions = 10, fixSamp = T, plot=F){
  
  n.all = nrow(data)
  xModel = as.formula(formula)
  n.cv = folds
  n.times = repetitions
  if(fixSamp){
    train.ind = sample( c(rep(0,n.all*0.3), rep(1,n.all*0.7)) )
  } else{
    train.ind = sample( c(0,1), size=n.all , replace=TRUE, prob=c(0.3, 0.7) )
  }
  xx.train = data[ train.ind == 1 , ]
  xx.test  = data[ train.ind == 0 , ] 
  
  ## Matrices for storing the goodness of fit indices
  w.r2 = matrix ( NA, n.times, n.cv )
  train.r2 = matrix ( NA, n.times, n.cv )
  w.rmse   = matrix ( NA, n.times, n.cv )
  train.rmse = matrix ( NA, n.times, n.cv )
  train.bic = matrix ( NA, n.times, n.cv )
  train.aic = matrix ( NA, n.times, n.cv )
  test.predict = list()
  k=0
  #do the magic
  for ( times in 1:n.times )
  { 
    ## Creates the vector "groups" which contains the group (fold) of each observation 
    
    groups = sample ( rep ( 1:n.cv, length=n.all ) )
    
    ## Cross validation
    for ( ind.cv in 1:n.cv )
    {
      k = k+1
      ## Creating the datasets for training and testing in the fold
      cv.train = data [ groups != ind.cv , ]
      cv.test  = data [ groups == ind.cv , ] 
      
      ## Fitting the linear model in the training subsample
      lm.out = lm(xModel, data=cv.train )
      train.r2 [times,ind.cv] <- summary(lm.out)$r.squared
      train.rmse [times,ind.cv] <- round((summary(lm.out)$sigma/nrow(cv.train))^1/2,3)
      train.aic [times, ind.cv] <- AIC(lm.out)
      train.bic [times, ind.cv] <- BIC(lm.out)
      
      ## Predictions in the testing sample
      pred.lm.cv.testing = predict ( lm.out , newdata=cv.test)
      
      #toDo for bands of confidence plot
      #plot(y~x1,data=data)
      test.predict[[k]] = pred.lm.cv.testing
      #points(cv.test$x1,pred.lm.cv.testing,col="blue")
      #lines(data$x1,all.predict[[k]],col="red")
      
      w.r2 [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$R2
      w.rmse [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$RMSE
      
    }
  }
  
  train.bic <- cbind(round(mean(train.bic),2),sd(as.vector(train.bic)))
  train.aic <- cbind(round(mean(train.aic),2),sd(as.vector(train.aic))) 
  
  training.r2 <-  cbind(round(mean(train.r2),2),sd(as.vector(train.r2)))
  testing.r2 <- cbind(round(mean(w.r2),2),sd(as.vector(w.r2)))
  
  training.rmse <- cbind(round(mean(train.rmse),2),sd(as.vector(train.rmse)))
  testing.rmse <- cbind(round(mean(w.rmse),2),sd(as.vector(w.rmse)))
  
  names(train.bic) = c("BIC Mean","BIC SD")
  names(train.aic) = c("AIC Mean","AIC SD")
  names(training.r2) = names(testing.r2) = c("r2_mean","r2_sd")
  names(training.rmse) = names(testing.rmse) = c("rmse_mean","rmse_sd")
  
  if(plot) {
    # rTrain = cv.train$y - fitted(lm(xModel, data=cv.train))
    # rTest = cv.test$y - pred.lm.cv.testing
    # plot(0,type="n",xlim=c(min(cv.test$y,cv.train$y,na.rm=T),max(cv.test$y,cv.train$y,na.rm=T)),
    #      ylim=c(min(rTrain,rTest),max(rTrain,rTest)))
    # points(cv.train$y,rTrain, col="orange", lwd=2)
    # points(cv.test$y, rTest, col="green",lwd=2)
    
    # layout(matrix(c(1,2,3,4,5,5),nrow=3,byrow=T))
    par(mfrow=c(1,2)) 
    plot(cv.train$y, main = paste("Training | R2:",round(mean(train.r2),2) ),xaxt="n",yaxt="n",xlab="Index",ylab="y",las=1,cex.axis=1.5)
    lines(fitted(lm(xModel, data=cv.train)),col="red")
    plot(cv.test$y, main = paste("Test | R2:",round(mean(testing.r2),2) ),xaxt="n",yaxt="n",xlab="Index",ylab="y",las=1,cex.axis=1.5)
    lines(pred.lm.cv.testing,col="red")
    # plot(y~x1,data = cv.train, main = paste("Training | R2:",round(mean(train.r2),2) ))
    # lines(cv.train$x1,fitted(lm(xModel, data=cv.train)),col="red")
    # plot(y~x1,data = cv.test, main = paste("Test | R2:",round(mean(testing.r2),2) ))
    # lines(cv.test$x1,pred.lm.cv.testing,col="red")
    
    # plot(y~x1,data)
    # bands = do.call(rbind,all.predict)
    # lines(data$x1,apply(bands,2,mean),col="blue")
    # lines(data$x1,apply(bands,2,function(x){
    #   mean(x) + sd(x)
    # }),col="blue")
    # lines(data$x1,apply(bands,2,function(x){
    #   mean(x) - sd(x)
    #   }),col="blue")
    
    #lines(predict(),col="red")
    # DiagPlot(cv.test$y, pred.lm.cv.testing)
    # DiagPlot(cv.train$y,  predict ( lm.out , newdata=cv.train))
    
  }
  
  
  output <- list(train = rbind(training.r2,training.rmse,train.aic,train.bic), test=rbind(testing.r2,testing.rmse))
  rownames(output$train) = c("Adjusted R squared","RMSE","AIC","BIC")
  rownames(output$test) = c("Adjusted R squared","RMSE")
  colnames(output$train) = colnames(output$test) = c("Mean","SD")
  # 
  # return(output)
  list("trainRmse"=round(mean(train.rmse),2),
       "testRmse" =round(mean(w.rmse),2),
       "trainR2"=round(mean(train.r2),2),
       "testR2" =round(mean(w.r2),2)
  )
}

#debug
# jaVal(as.formula(xModels[4]), dataGen(n=100,4),folds = 5,repetitions = 10)
gen4 = function(x1,x2,x3,x4){
  n = length(x1)
  return(   0
            +3*x1
            + 3*x1^2 
            + 1*x2
            + 1*x1*x2 
            + 1*x1*x3 
            + 1*x2*x3 
            + 0.008*x4 
            + 0.01*rnorm(length(x1),0,2*n/10)*sin(x1)
            + 0.01*sin(x2)*runif(n,-sqrt(n),sqrt(n))
            *sample(c(0,1),n,prob = c(0.04,0.96),replace = T)
            
            +rnorm(n = n, mean = 0, sd = 5000)
  )
}


# gen5 = function(x1,x2,x3,x4){
#   n = length(x1)
#   return(   0
#             +3*x1
#             + 3*x1^2
#             + 1*x2
#             + 1*x1*x2 
#             + 1*x1*x3 
#             + 1*x2*x3 
#             + 0.008*x4 
#             + 0.01*rnorm(length(x1),0,2*n/10)*sin(x1)
#             #+ 0.0001*sin(x2)*runif(n,-n,n)
#             #*sample(c(0,1),n,prob = c(0.04,0.96),replace = T)
#             #*n^1.5
#   )
# }

EvalRegr = function ( obs, pred )
{ 
  
  if ( length(obs) != length(pred) ) 
  { stop ( "Vectors must have the same length!" ) }   
  
  ## MSE and RMSE
  mse = mean ( ( obs - pred ) ** 2 )
  rmsex = sqrt ( mse)
  
  ## Standard Error of the MSE
  
  se.mse = sd ( ( obs - pred ) ** 2 ) / sqrt(length(obs))  
  
  ## MAE  
  mae = mean ( abs( obs - pred ) )
  
  ## R2
  R2 = cor ( obs, pred ) ** 2
  
  ## Outout
  list( MSE = mse, RMSE = rmsex, MAE = mae, R2 = R2, SE.MSE = se.mse )
}


dataGen = function(n, generator=4, plot=F) {
  if(generator!=4)stop("This version only supports gen4")
  x1 = sort(rnorm(n,15,30))
  x2 = rnorm(n,20,8)
  x3 = rnorm(n,4,70)
  x4 = rep(runif(4,-10,40),ceiling(n/4))[1:n]
  y = gen4(x1,x2,x3,x4)
  ds = data.frame("y"=y,"x1"=x1,"x2"=x2,"x3"=x3,"x4"=x4)
  if(plot)plot(ds$y,main=paste0("Blackbox | n=",n))
  return(ds)
}
