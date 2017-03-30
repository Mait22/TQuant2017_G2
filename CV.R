

library(pROC)
library(caret)

xx <- read.csv("badData.csv", sep=" ")
xx <- xx[sample(1:500, size = 150, replace = F),]
dim(xx)

y = xx$y

n.all = nrow ( xx )  ## Número de observaciones


################################################################################
################################################################################
## 1.- Muestra de Training y Muestra de Testing
################################################################################
################################################################################

set.seed(111) ## Se fija una semilla, para reproducir los datos

## train.ind = 1 para la Muestra de Training (70%) 
## train.ind = 0 para la Muestra de Testing  (30%)
train.ind = sample ( c(0,1), size=n.all , replace=TRUE, prob=c(0.3, 0.7) )
table(train.ind)
prop.table ( table(train.ind) )

## Definición de la Muestra de Training y Muestra de Testing

xx.train = xx [ train.ind == 1 , ]
xx.test  = xx [ train.ind == 0 , ] 
dim(xx.train)
dim(xx.test)



################################################################################
################################################################################
## 10 times 10-fold CV
################################################################################
################################################################################

## 10 times 10-fold CV
n.cv = 5
n.times = 10

## Matriz para almacenar los errores y los AUCs

w.r2 = matrix ( NA, n.times, n.cv )
train.r2 = matrix ( NA, n.times, n.cv )
train.rmse = matrix ( NA, n.times, n.cv )
w.rmse   = matrix ( NA, n.times, n.cv )
train.rmse = matrix ( NA, n.times, n.cv )



for ( times in 1:n.times )
{ 
  
  ## Crea el vector "groups" que contiene el grupo (fold) de cada observación  
  groups = sample ( rep ( 1:n.cv, length=n.all ) )
  
  ## Validación Cruzada
  for ( ind.cv in 1:n.cv )
  {
    ## Creamos los datasets de training y testing en este fold
    cv.train = xx [ groups != ind.cv , ]
    cv.test  = xx [ groups == ind.cv , ] 
    
    ## Regresión Logística ajustada con la muestra de Training
    lm.out = lm(y~x1+I(x1^2)+I(x1^3)+I(x1^4)+I(x1^5)+I(x1^6)+I(x1^7)+I(x1^8) + x1*x2 + I(x2^2) + I(x2^3), data=cv.train )
    train.r2 [times,ind.cv] <- summary(lm.out)$r.squared
    train.rmse [times,ind.cv] <- round((summary(lm.out)$sigma/150)^1/2,3)
    
    ## Predicciones en la Muestra de Testing
    pred.lm.cv.testing = predict ( lm.out , newdata=cv.test)
    EvalRegr ( cv.test$y , pred.lm.cv.testing )
    
    w.r2 [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$R2
    w.rmse [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$RMSE
    w.rmse.SE [times,ind.cv] = EvalRegr ( cv.test$y , pred.lm.cv.testing )$SE.MSE
    
  }
}

training.r2 <-  cbind(round(mean(train.r2),2),sd(as.vector(train.r2)))
testing.r2 <- cbind(round(mean(w.r2),2),sd(as.vector(w.r2)))

training.rmse <- cbind(round(mean(train.rmse),2),sd(as.vector(train.rmse)))
testing.rmse <- cbind(round(mean(w.rmse),2),sd(as.vector(w.rmse)))

DiagPlot(cv.test$y, pred.lm.cv.testing)
DiagPlot(cv.train$y,  predict ( lm.out , newdata=cv.train))


################################################################################
################################################################################
## C A R E T  -  C R O S S   V A L I D A T I O N
################################################################################
################################################################################


cv.ctrl <- trainControl(method = 'repeatedCV', number = 5, repeats = 10,
             summaryFunction = "defaultSummary" )

lm.fit <- train( xx.train[-1], xx.train$y,
                 method = 'lm',
                 #trControl = cv.ctrl,
                 metric = 'Rsquared',
                 preProcess = "center", 
                 prob.model = TRUE)


lm.fit$bestTune
lm.fit$results
lm.fit$finalModel
