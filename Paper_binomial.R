
# Models
m.b <- function(param, t){
  (param[2] + param[3]*t)^-param[1]
}

m.a <- m.1 <- function(param, t){
  (1 + t)^-param[1]
}

m.2 <- function(param, t){
  (param[2] + t)^-param[1]
}

m.3 <- function(param, t){
  (1 + param[2]*t)^-param[1]
}

# SSE sum(observe - prediction)^2

# SST sum(observe - expected mean)^2

## data generation
t <- seq(0.1, 8.1, by=2)
n <- 50
a <- .4
b <- 5
c <- .2

dat.n <- 1000
dat.n2 <- 500
error <- .9 # making it not homogenous? Drawing it from a Distribution, CHi^2?`

## Frequencies
yi.1 <- sapply(t, 
               function(x){
                 rbinom(1, n, m.1(a, x))})

yi.2 <- sapply(t,
               function(x){
                 rbinom(1, n, m.2(c(a, b), x))
               })
  
yi.3 <- sapply(t,
               function(x){
                 rbinom(1, n, m.3(c(a, b), x))
               })

yi.b <- sapply(t,
               function(x){
                 rbinom(1, n, m.b(c(a, b, c), x))
               })

### Likelihoods
loglike.1 <- function(param, data, t){
  prob <- prob <- sapply(t, 
                         function(x){
                           m.1(param, x)})
  prob <- c(prob, 1-prob)
  -sum(log(sapply(data, 
      function(x)
              {choose(50, x)})) + data*log(prob))
  }

loglike.2 <- function(param, data, t){
  prob  <- sapply(t, 
                  function(x){
                    m.2(param, x)})
  prob <- c(prob, 1-prob) 
  #c((param[2] + t)^-param[1], 1-(param[2] + t)^-param[1])
  -sum(log(sapply(data, 
              function(x)
              {choose(50, x)})) + data*log(prob))
  }

loglike.3 <- function(param, data, t){
  prob  <- sapply(t, 
                  function(x){
                    m.3(param, x)})
  prob <- c(prob, 1-prob) 
  #c((1 + param[2]*t)^-param[1], 1- (1 + param[2]*t)^-param[1])
  -sum(log(sapply(data, 
              function(x)
              {choose(50, x)})) + data*log(prob))
  }

loglike.b <- function(param, data, t){
  prob  <- sapply(t, 
                  function(x){
                    m.b(param, x)})
  prob <- c(prob, 1-prob) 
  # c((param[2] + param[3]*t)^-param[1], 1 - (param[2] + param[3]*t)^-param[1])
  -sum(log(sapply(data, 
              function(x)
              {choose(50, x)})) +  data*log(prob))
}

## Getting the binomial Coef
# sapply(x, function(x){choose(50, x)}))

## maximising likelihood
max.1 <- nlm(loglike.1, data=c(yi.1, n-yi.1), t=t, p=c(0.9))

max.2 <- nlm(loglike.2, data=c(yi.2, n-yi.2), t=t, p=c(0.5, 0.5))

max.3 <- nlm(loglike.3, data=c(yi.3, n-yi.3), t=t, p=c(0.5, 0.5))

max.b <- nlm(loglike.b, data=c(yi.b, n-yi.b), t=t, p=c(0.5, 0.5, 0.5))

### generating Data

dat.m1 <- data.frame(ID=seq_len(dat.n),
                     time=rep(t, each=dat.n/length(t)),
                     prob=sapply(dat.m1$time, 
                                 function(x){
                                   m.1(max.1$estimate, x)}))
dat.m1$m_distr <- n*dat.m1$prob

dat.m1$frq_obs <- sapply(dat.m1$prob, 
                         function(x){
                           rbinom(1, n, x * error)
                           })
dat.m1$frq_nobs <- n - dat.m1$frq_obs
dat.m1$frq_prd <- sapply(dat.m1$prob, 
                        function(x){
                          rbinom(1, n, x)
                        })

## RMSE
# from Tutorial on maximum likelihood estimation
# In Jae Myung
dat.m1$SE <- (dat.m1$frq_obs - dat.m1$frq_prd)^2

SSE.m1.1 <- sum(dat.m1$SE)

RMSE.m1.1 <- sqrt(SSE.m1.1/length(dat.m1$frq_obs))

# #SST.m1.1 <- sum((dat.m1$frq_obs - dat.m1$m_distr)^2) 
# # ?, seems strange.
# #
# # What about
 SST.m1.1 <- sum((dat.m1$frq_obs - mean(dat.m1$frq_obs))^2)
# # this averages over the observeations, like a 0 model.
# # But this could also be achieved by averaging the probability and then computing
# # n*p?
# # yeah close
# mean(dat.m1$prob)*50
# # 30.61481
# mean(dat.m1$frq_obs)
# # 27.164
# # Still, for a deviation of .8, through the sampling error, this seems too large a number

#PVAF:
PVAF.m1.1 <- (1 - (SSE.m1.1/SST.m1.1))

## AIC
b.aic <- function(loglike, nparam){
  -2*loglike + 2*nparam
}
b.bic <- function(loglike, nparam, sample){
  -2*loglike + nparam*log(sample)
}

m1.1.aic <- b.aic(loglike.1(param=max.1$estimate, 
                            data=c(dat.m1$frq_obs,
                                   dat.m1$frq_nobs),
                            t=rep(dat.m1$time)),
                  length(max.1$estimate))

m1.1.bic <- b.bic(loglike.1(param=max.1$estimate, 
                            data=c(dat.m1$frq_obs,
                                   dat.m1$frq_nobs),
                            t=rep(dat.m1$time)),
                  length(max.1$estimate),
                  length(dat.m1$frq_obs))
  
  
