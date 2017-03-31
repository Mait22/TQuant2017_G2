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
            + 0.0001*rnorm(length(x1),0,2*n/10)*sin(x1)*n^2
            + 0.0001*sin(x2)*runif(n,-n,n)
            *sample(c(0,1),n,prob = c(0.04,0.96),replace = T)
            *n^1.5
            +rnorm(n = n, mean = 0, sd = 5000)
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
    x1 = sort(rnorm(n,15,30))
    x2 = rnorm(n,20,8)
    x3 = rnorm(n,4,70)
    x4 = rep(runif(4,-10,40),ceiling(n/4))[1:n]
    y = gen4(x1,x2,x3,x4)
  }
  ds = data.frame("y"=y,"x1"=x1,"x2"=x2, "x3" = x3, "x4" = x4)
  if(plot)plot(ds$y,main=paste0("Blackbox | n=",n))
  return(ds)
}

#check the generated data
par(mfrow=c(1,1))
plot(dataGen(100,plot=T,generator = 4)$y)

