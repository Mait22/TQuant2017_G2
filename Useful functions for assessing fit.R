library(pROC)

################################################################################
## ROC Analysis
################################################################################
     
ROCAnalysis = function ( bin.var, cont.var, last = TRUE, 
                         ci = FALSE, plot = TRUE, ... )
{
  ## Checking dichotomous dependent variables
  if ( length( unique (bin.var ) ) != 2 ) 
  { print ( bin.var ) ; stop ( "Variable is not binomial!" ) }   
  
  ## ROC Analysis
  
  if ( last == TRUE ) { roc.out <- roc( bin.var , cont.var, direction="<" ) }
  else                { roc.out <- roc( bin.var , cont.var, direction=">" ) }

  ## 95% IC AUC 
  w.auc  = pROC::auc (roc.out) 
  if ( ci == TRUE ) { w.ci   = ci.auc (roc.out)  } 

  ## Plot ROC (opcional)
  if ( plot == TRUE )
  {
    dev.new()
    plot( roc.out , legacy.axes = TRUE, col="red" , ... )
    text (0.2, 0.2, paste("AUC = ", round(w.auc ,3 )))
  }

  ## Output
  if ( ci == TRUE ) { list ( AUC = w.auc[1] , CI.AUC = c( w.ci[1] , w.ci[3] ) ) }
  else              { list ( AUC = w.auc[1] ) }
}


################################################################################
## Error Rate
################################################################################

ErrRate = function ( obs, pred )
{ 
  
  if ( length(obs) != length(pred) ) 
  { stop ( "Vectors must have the same length!" ) }   
  
  ## Table 
  tt = table( obs, pred )
   
  ## Proportion of correct and incorrect classifications
  
  accuracy = sum ( as.character(obs) == as.character(pred) ) 
 
  accuracy = accuracy / sum(tt)
  
  err.rate = 1 - accuracy
  
  list( accuracy = accuracy, err.rate = err.rate )
}


################################################################################
## Evaluations methods for regression models (continous DV)
################################################################################

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

################################################################################
## Diagnistic plots for regression models (continous DV)
################################################################################

DiagPlot = function ( obs, pred, ... )
{ 
  
  if ( length(obs) != length(pred) ) 
  { stop ( "Vectors must have the same length!" ) }   
  
  ## Residuals
  res = obs - pred

  ## Graph of the observed versus predicted values
  
  w.min = min ( pred, obs )
  w.max = max ( pred, obs )  
  dev.new()
  plot ( pred, obs, 
         xlim = c( w.min, w.max ), ylim = c( w.min, w.max ),  
         xlab="Predicted", ylab="Observed" , pch=16, ... )
  abline ( b=1, a=0 )          

  ## Graph of the residuals versus predicted values
  
  dev.new()
  plot ( pred, res, 
         xlab="Predicted", ylab="Residuals" , pch=16, ... )
  abline ( h=0 )   
  
}


################################################################################
## 2D graph of the region of decision of a classifier
################################################################################

DecisionPlot = function ( w.points, w.grid, w.pred, ... )
{
   dev.new()
   plot   ( w.grid [w.pred == 0, 1 ] , w.grid [w.pred == 0, 2 ] , pch=16, cex=0.3,
               xlab="X1", ylab="X2", col= "cadetblue1", xlim=c(-3, 4.5 ), ylim=c (-2.5, 3 ), ... )
   points ( w.grid [w.pred == 1, 1 ] , w.grid [w.pred == 1, 2 ] , pch=16, cex=0.3,
               col= "tan1" )

   points ( w.points [ w.points$y == 0 , 1] , w.points [ w.points$y == 0 , 2] , col= "blue", pch=16, cex=0.8 )
   points ( w.points [ w.points$y == 1 , 1] , w.points [ w.points$y == 1 , 2] , col= "red2", pch=16, cex=0.8 )
}


################################################################################
## Variable selection filter based on p-values with binary DV
################################################################################
     
Pval.class = function ( x.var , y.var )
{ 
  
  if ( length(unique(y.var)) != 2 ) stop ( "DV must be binary" )
  
  
  p.val = rep ( NA, ncol(x.var) )
  
  for ( i in 1:ncol(x.var) )
  { 

    ## Continous predictors
    if ( length(unique( x.var[ ,i] )) > 10  & is.factor( x.var[ ,i] ) == FALSE )
    { p.val [i] = t.test ( x.var[ ,i] ~ y.var )$p.value }
    else
    { p.val [i] = chisq.test ( y.var , x.var[ ,i] )$p.value }    
  
  }
  
  ## Output
  list ( p.val = p.val )
}
