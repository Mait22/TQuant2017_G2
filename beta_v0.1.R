#setwd("C:\\Users\\Johann\\OneDrive - Università degli Studi di Padova\\_Statistica_e_metodi\\_Corsi e dispense\\TquanT_2017\\app")
EvalRegr = function ( obs, pred )
{ 
  
  if ( length(obs) != length(pred) ) 
  { stop ( "Vectors must have the same length!" ) }   
  
  ## MSE and RMSE
  mse = mean ( ( obs - pred ) ** 2 )
  rmse = sqrt ( mse)
  
  ## Standard Error of the MSE
  
  se.mse = sd ( ( obs - pred ) ** 2 ) / sqrt(length(obs))  
  
  ## MAE  
  mae = mean ( abs( obs - pred ) )
  
  ## R2
  R2 = cor ( obs, pred ) ** 2
  
  ## Outout
  list( MSE = mse, RMSE = rmse, MAE = mae, R2 = R2, SE.MSE = se.mse )
}
gen = function(x1,x2){
  return( -500 + 8.5*x1 +300*x1^2 + 0.2*x1*x2 + 29600*rnorm(length(x1),0,90)*sin(x1) + 3000*sin(x2)*runif(length(x1),-15570,15570)*sample(c(0,1),length(x1),prob = c(0.04,0.96),replace = T) )
}
gen2 = function(x1,x2){
  return( -500 + 8.5*x1 +300*x1^2 + 0.2*x1*x2 + 29600*rnorm(length(x1),0,90)*sin(x1) + 3000*sin(x2)*runif(length(x1),-15570,15570)*sample(c(0,1),length(x1),prob = c(0.04,0.96),replace = T) )
}
gen3 = function(x1,x2){
  n = length(x1)
  return(   0 +
            3*x1 +
            7.2*x1^2 +
            0.008*x1*x2*n +
            0.01*rnorm(length(x1),0,2*n/10)*sin(x1)*n^2 +
            0.01*sin(x2)*runif(n,-n,n)
              #*sample(c(0,1),n,prob = c(0.04,0.96),replace = T)
              *n^1.5
            )
}
gen4 = function(x1,x2,x3,x4){
  n = length(x1)
  return(   0 +
              3*x1 +
              720*x1^2 +
              1*x2+
              1*x1*x2 +
              # 1*x1*x3 +
              # 1*x2*x3 +
              # 0.008*x4 +
              0.01*rnorm(length(x1),0,2*n/10)*sin(x1)*n^2 +
              0.01*sin(x2)*runif(n,-n,n)
            #*sample(c(0,1),n,prob = c(0.04,0.96),replace = T)
            *n^1.5
  )
}

dataGen = function(n, generator = 1,plot=F) {
  x1 = -(n/4):(n*3/4-1)
  x2 = rnorm(n,5,30)
  if(generator == 1){
    y = gen(x1,x2)
  } else if(generator ==2){
    y = gen2(x1,x2)
  }else if(generator ==3){
    y = gen3(x1,x2)
  }else if(generator ==4){
    x1 = rnorm(n,2,30)
    x2 = rnorm(n,20,8)
    x3 = rnorm(n,4,70)
    x4 = rep(runif(4,0,40),trunc(n/4))
    y = gen4(x1,x2,x3,x4)
  }
  ds = data.frame("y"=y,"x1"=x1,"x2"=x2)
  if(plot)plot(ds$y,main=paste0("Blackbox | n=",n))
  return(ds)
}

#check the generated data
par(mfrow=c(1,2))
dataGen(1000,plot=T,generator = 3)
dataGen(1000,plot=T,generator = 3)
par(mfrow=c(1,1))



set.seed(111) ## Se fija una semilla, para reproducir los datos

jaVal = function(formula, data, folds=5, repetitions = 10, fixSamp = T){
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
  ## Matriz para almacenar los errores y los AUCs
  w.r2 = matrix ( NA, n.times, n.cv )
  train.r2 = matrix ( NA, n.times, n.cv )
  train.rmse = matrix ( NA, n.times, n.cv )
  w.rmse   = matrix ( NA, n.times, n.cv )
  train.rmse = matrix ( NA, n.times, n.cv )
  test.predict = list()
  k=0
  #do the magic
  for ( times in 1:n.times )
  { 
    ## Crea el vector "groups" que contiene el grupo (fold) de cada observaci?n  
    groups = sample ( rep ( 1:n.cv, length=n.all ) )
    
    ## Validaci?n Cruzada
    for ( ind.cv in 1:n.cv )
    {
      k = k+1
      ## Creamos los datasets de training y testing en este fold
      cv.train = data [ groups != ind.cv , ]
      cv.test  = data [ groups == ind.cv , ] 
      
      ## Regresi?n Log?stica ajustada con la muestra de Training
      lm.out = lm(xModel, data=cv.train )
      train.r2 [times,ind.cv] <- summary(lm.out)$r.squared
      train.rmse [times,ind.cv] <- round((summary(lm.out)$sigma/nrow(cv.train))^1/2,3)
      
      ## Predicciones en la Muestra de Testing
      pred.lm.cv.testing = predict ( lm.out , newdata=cv.test)
      #toDo for bands of confidence plot
      #test.predict[[k]] = pred.lm.cv.testing 
      
      EvalRegr ( cv.test$y , pred.lm.cv.testing )
      
      w.r2 [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$R2
      w.rmse [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$RMSE

    }
  }
  training.r2 <-  cbind(round(mean(train.r2),2),sd(as.vector(train.r2)))
  testing.r2 <- cbind(round(mean(w.r2),2),sd(as.vector(w.r2)))
  names(training.r2) = names(testing.r2) = c("r2_mean","r2_sd")
  

  par(mfrow=c(1,2))
  plot(cv.train$y, main = "Training")
  lines(fitted(lm(xModel, data=cv.train)),col="red")
  plot(cv.test$y, main = "Test")
  lines(pred.lm.cv.testing,col="red")
  
  
  #lines(predict(),col="red")
  
  
  training.rmse <- cbind(round(mean(train.rmse),2),sd(as.vector(train.rmse)))
  testing.rmse <- cbind(round(mean(w.rmse),2),sd(as.vector(w.rmse)))
  # DiagPlot(cv.test$y, pred.lm.cv.testing)
  # DiagPlot(cv.train$y,  predict ( lm.out , newdata=cv.train))
  return(list(train = training.r2, test=testing.r2))
}


jaVal("y~x1+I(x1^2)+I(x1^3)+I(x1^4)+I(x1^5)+I(x1^6)+I(x1^7)+I(x1^8) + x1*x2 + I(x2^2) + I(x2^3)",
      dataGen(n=100,1),folds = 5,repetitions = 10)

jaVal("y~x1+I(x1^2)+I(x1^3)+I(x1^4)+I(x1^5)+I(x1^6)+I(x1^7)+I(x1^8) ",
      dataGen(n=1000,3),folds = 5,repetitions = 10)
