
################################################################################
################################################################################
## 1.- Muestra de Training y Muestra de Testing
################################################################################
################################################################################

rmse <- function(error){
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
    c(rmse(xx.train$y-fitted(model.out)),
         rmse(xx.testing$y-predict( model.out, newdata=xx.testing) )
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


