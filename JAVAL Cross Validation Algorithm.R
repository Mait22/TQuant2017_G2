
library(caret)

##########################################################################################
########### J A V A L  C R O S S   V A L I D A T I O N   A L G O R I T H M ###############
##########################################################################################

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
      #test.predict[[k]] = pred.lm.cv.testing 
      
      w.r2 [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$R2
      w.rmse [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$RMSE
      
    }
  }
  
  training.r2 <-  cbind(round(mean(train.r2),2),sd(as.vector(train.r2)))
  testing.r2 <- cbind(round(mean(w.r2),2),sd(as.vector(w.r2)))
  
  training.rmse <- cbind(round(mean(train.rmse),2),sd(as.vector(train.rmse)))
  testing.rmse <- cbind(round(mean(w.rmse),2),sd(as.vector(w.rmse)))
  
  names(training.r2) = names(testing.r2) = c("r2_mean","r2_sd")
  names(training.rmse) = names(testing.rmse) = c("rmse_mean","rmse_sd")
  
  par(mfrow=c(1,2))
  plot(cv.train$y, main = "Training")
  lines(fitted(lm(xModel, data=cv.train)),col="red")
  plot(cv.test$y, main = "Test")
  lines(pred.lm.cv.testing,col="red")
  
  
  #lines(predict(),col="red")
  # DiagPlot(cv.test$y, pred.lm.cv.testing)
  # DiagPlot(cv.train$y,  predict ( lm.out , newdata=cv.train))
  
  output <- list(train = rbind(training.r2,training.rmse), test=rbind(testing.r2,testing.rmse))
  rownames(output$train) = rownames(output$test) = c("Adjusted R squared","RMSE")
  colnames(output$train) = colnames(output$test) = c("Mean","SD")
  
  return(output)
}

##################################
############ EXAMPLE #############
##################################

jaVal(formula = "y~1 + x1 + x2 + x3 + I(x1^2) + I(x1^3) + x4 + x1*x2 + I(x2^2) + I(x2^3) + x1*x3 + I(x3^2) + I(x3^3) + x2*x3 + I(x1^4) + I(x1^5) + I(x1^6) + I(x1^7) + I(x1^8) + I(x1^9) + I(x1^10) + I(x1^11) + I(x1^12) + I(x1^13) + I(x1^14) + I(x1^15) + I(x1^16)",
      dataGen(n=100,4),folds = 5,repetitions = 10, fixSamp = T)
