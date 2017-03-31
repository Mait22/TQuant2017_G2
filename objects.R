library(ggplot2)
library(ggthemes)
library(ggiraph)

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


test <- paper_data_gen(1,2,3,100,noise.prop = 0.1)





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


test22 <- paper_data_fitter(test)










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
  
  raw.results.pVaf <- as.data.frame(matrix(rep(NA, times = 10*n.sims),nrow = n.sims,ncol = 10))
  names(raw.results.pVaf) <- c("Simulation run", 
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
    max.1.m1.n <-  nlm(loglike.1, data=c(yi.1.n, n-yi.1.n), t=data1.n$time, p=c(0.8))
    max.1.m2.n <- nlm(loglike.2, data=c(yi.1.n, n-yi.1.n), t=data1.n$time, p=c(.8, 0.99))
    max.1.m3.n <- nlm(loglike.3, data=c(yi.1.n, n-yi.1.n), t=data1.n$time, p=c(0.8, 0.5))
    
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
    raw.results.rmse[i,"Data from M1 fitted by M1 (BIC)"] <-  rmse(loglike = loglike.1(max.1.m1.n$estimate,
                                                                                         data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                     nparam = length(max.1.m1.n$estimate),sample = length(yi.1.n))
    
    raw.results.rmse[i,"Data from M1 fitted by M2 (BIC)"] <-  rmse(loglike = loglike.2(max.1.m2.n$estimate,
                                                                                         data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                     nparam = length(max.1.m2.n$estimate),sample = length(yi.1.n))
    
    raw.results.rmse[i,"Data from M1 fitted by M3 (BIC)"] <-  rmse(loglike = loglike.3(max.1.m3.n$estimate,
                                                                                         data=c(yi.1.n, n-yi.1.n), t = data1.n$time),
                                                                     nparam = length(max.1.m3.n$estimate),sample = length(yi.1.n))
    
      
      
      
      
    yi.2.n <- sapply(data2.n$prob.m2, function(x){rbinom(1, n, x)})
    
    max.2.m1.n <-  nlm(loglike.1, data=c(yi.2.n, n-yi.2.n), t=data2.n$time, p=c(0.8))
    max.2.m2.n <- nlm(loglike.2, data=c(yi.2.n, n-yi.2.n), t=data2.n$time, p=c(.8, 0.99))
    max.2.m3.n <- nlm(loglike.3, data=c(yi.2.n, n-yi.2.n), t=data2.n$time, p=c(0.8, 0.5))
    
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
    
    
    yi.3.n <- sapply(data3.n$prob.m3, function(x){rbinom(1, n, x)    
    })
    
    max.3.m1.n <-  nlm(loglike.1, data=c(yi.3.n, n-yi.3.n), t=data3.n$time, p=c(0.8))
    max.3.m2.n <- nlm(loglike.2, data=c(yi.3.n, n-yi.3.n), t=data3.n$time, p=c(.8, 0.99))
    max.3.m3.n <- nlm(loglike.3, data=c(yi.3.n, n-yi.3.n), t=data3.n$time, p=c(0.8, 0.5))
    
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
    
  }

  return_object@raw.results.aic <- raw.results.aic
  return_object@raw.results.bic <- raw.results.bic
  
  return(return_object)
  

}





simulator_aicbic(n.sims = 10,n.total = 100,noise.prop = 0.1,a = 1,b = 2,c = 3)





















bic.new <- function(loglike, nparam, sample){
  -2*loglike + nparam*log(sample)
}

rmse <- function(obs, pred){
  sqrt(sum((obs - pred)^2) / length(obs))
  ## or:
  # sqrt(mean((obs - pred)^2))
}

pvaf <- function(obs, pred){
  1 - (sum((obs - pred)^2) / sum((obs - mean(obs))^2))
}

## Example on how to use it,
## dat.1 ist teh dataframe in which frequencies and probabilities as well as t are stored
#liste <- list(list(max.1.n, loglike.1),
#              list(max.2.n, loglike.2),
#              list(max.3.n, loglike.3))
#
#aic.1 <- data.frame(Model=c(1,2,3),
#                    AIC=unlist(
#                      lapply(liste,
#                             function(x){
#                               b.aic(
#                                 x[[2]](param=x[[1]]$estimate,
#                                           data=c(dat.1$frq_obs.m1,
#                                           dat.1$frq_nobs.m1),
#                                           t=dat.1$time),
#                                 length(x[[1]]$estimate))})),
#                    BIC=unlist(
#                      lapply(liste,
#                             function(x){
#                               b.bic(
#                                 x[[2]](param=x[[1]]$estimate,
#                                        data=c(dat.1$frq_obs.m1,
#                                               dat.1$frq_nobs.m1),
#                                        t=dat.1$time),
#                                 length(x[[1]]$estimate),
#                                 length(dat.1$frq_obs.m1))
#                             })
#                    ))




