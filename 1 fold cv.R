
parList = c("x1", "x2", "x3", "I(x1^2)","I(x1^3)",
            "x4","x1*x2","I(x2^2)","I(x2^3)",
            "x1*x3","I(x3^2)","I(x3^3)",
            "x2*x3")
parList = c(parList,paste0("I(x1^",4:16,")")
)
newmula = "y~1"
xModels = c(newmula)
# m0 = lm (y~1,ds)


#i is the slider value for complexity
for (i in 1:26){
  newmula = paste(newmula, parList[i],sep=" + ")
  xModels = c(xModels,newmula)
}

################################################################################
################################################################################
## 1.- Muestra de Training y Muestra de Testing
################################################################################
################################################################################

rmse <- function(error){
  sqrt(mean(error^2))
}




ml.out <- list()
pred.lm.testing <- list()

i = 1
m = 2

trainlist <- list()
testinglist <- list()

for (i in 1:100){

xx.train <- dataGen(150,plot=F,generator = 4)
xx.testing <- dataGen(150,plot=F,generator = 4)
n.all = nrow ( xx.train )

#set.seed(111) 

# ## train.ind = 1 para la Muestra de Training (70%) 
# ## train.ind = 0 para la Muestra de Testing  (30%)
# train.ind = sample ( c(0,1), size=n.all , replace=TRUE, prob=c(0.3, 0.7) )
# table(train.ind)
# prop.table ( table(train.ind) )
# 
# ## Definición de la Muestra de Training y Muestra de Testing
# xx.train = xx [ train.ind == 1 , ]
# xx.test  = xx [ train.ind == 0 , ] 
# dim(xx.train)
# dim(xx.test)


####################################
## Regression

## Regresión Logística ajustada con la muestra de Training

trainingv <- c()
testingv <- c()
dftv <- c()

for (m in 1:length(xModels)){
  
  model.out = lm(xModels[m], data=xx.train )
  tempList = predict ( model.out, newdata=xx.testing)
  dftv = c(dftv,summary(model.out)$df[2])
  
  rmse.training <- rmse(xx.train$y-fitted(model.out))
  rmse.testing <-  rmse(xx.testing$y-tempList )
  
  
  trainingv[m] <- rmse.training
  testingv[m] <- rmse.testing
  
}


trainlist[[i]] <- trainingv
testinglist[[i]] <- testingv

}

rmseTr <- do.call(rbind,trainlist)
rmseTs <- do.call(rbind,testinglist)

plot(dftv,rmseTr[1,],type = "n", xlim = rev(range(dftv)))
apply(rmseTr,1, function(x){
  lines(dftv,x,col=rgb(1,0,0,0.1),lwd=2)
})

apply(rmseTs,1, function(x){
  lines(dftv,x,col=rgb(0,0,1,0.1),lwd=2)
})

averagermseTr <- apply(rmseTr,2,mean)
averagermseTs <- apply(rmseTs,2,mean)

lines(dftv,averagermseTr,col = "red", lwd = 3)
lines(dftv,averagermseTs,col = "blue", lwd = 3)

str(summary(model.out))

