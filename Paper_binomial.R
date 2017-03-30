###############
###
### Model 1,3 seem to work, 2 and b doesnt
###
###############
# Models
m.b <- function(param, t){
  (param[2] + param[3]*t)^(-param[1])
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
yt <- rbinom(5, 50, .8)
yt <- yt[order(yt, decreasing=TRUE)]
n <- 50
a <- 0.4
b <- 2
c <- .4

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
  prob <- m.1(param, t)
  prob <- c(prob, 1-prob)
  -sum(sapply(data, 
              function(x)
              {lchoose(n, x)}) + data*log(prob))
}

# loglike.1 <- function(param, data, t){
#   prob <- sapply(t, 
#                  function(x){
#                    m.1(param, x)})
#   prob <- c(prob, 1-prob)
#   -sum(data*log(prob))
# }

loglike.2 <- function(param, data, t){
  prob  <- m.2(param, t)
  prob <- c(prob, 1-prob) 
  #c((param[2] + t)^-param[1], 1-(param[2] + t)^-param[1])
  -sum(sapply(data, 
                  function(x)
                  {lchoose(n, x)}) + data*log(prob))
}

loglike.3 <- function(param, data, t){
  prob <- m.3(param, t)
  prob <- c(prob, 1-prob) 
  #c((1 + param[2]*t)^-param[1], 1- (1 + param[2]*t)^-param[1])
  -sum(sapply(data, 
                  function(x)
                  {lchoose(n, x)}) + data*log(prob))
}

loglike.b <- function(param, data, t){
  prob  <- m.b(param, t)
  prob <- c(prob, 1-prob) 
  # c((param[2] + param[3]*t)^-param[1], 1 - (param[2] + param[3]*t)^-param[1])
  -sum(sapply(data, 
                  function(x)
                  {lchoose(n, x)}) +  data*log(prob))
}

## Getting the binomial Coef
# sapply(x, function(x){choose(50, x)}))

## maximising likelihood
max.1 <- nlm(loglike.1, data=c(yi.1, n-yi.1), t=t, p=c(0.8))

# optim(loglike.1, data=c(yi.1, n-yi.1), t=t, 
#       p=c(0.8),lower=.1,upper=5)
# 
# m.1(5,t)
# 
# aa=seq(.1,2,.01)
# faa <- sapply(aa, function(x){
#       loglike.1(x, data=c(yi.1, n-yi.1), t=t)})
# plot(aa,faa,type="l")

max.2 <- nlm(loglike.2, data=c(yi.2, n-yi.2), t=t, p=c(.8, 0.99))

max.3 <- nlm(loglike.3, data=c(yi.3, n-yi.3), t=t, p=c(0.5, 0.6))

max.b <- nlm(loglike.b, data=c(yi.b, n-yi.b), t=t, p=c(0.8, 0.99, 0.2))

### generating Data

dat.1 <- data.frame(ID=seq_len(dat.n),
                    time=rep(t, each=dat.n/length(t)))
dat.1$prob.m1 <- sapply(dat.1$time, 
                        function(x){
                          m.1(max.1$estimate, x)})

dat.1$frq_obs.m1 <- sapply(dat.1$prob.m1, 
                           function(x){
                             rbinom(1, n, x * error)
                           })
dat.1$frq_nobs.m1 <- n - dat.1$frq_obs.m1
dat.1$m_distr <- n*dat.1$prob.m1
max.2.n <- nlm(loglike.2, 
               data=c(dat.1$frq_obs.m1, dat.1$frq_nobs.m1),
               t=dat.1$time, p=c(.8, 0.99))

dat.1$prob.m2 <- sapply(dat.1$time,
                        function(x){
                          m.2(max.2.n$estimate, x)})

max.3.n <- nlm(loglike.3, 
               data=c(dat.1$frq_obs.m1, dat.1$frq_nobs.m1),
               t=dat.1$time, p=c(.8, 0.99))
dat.1$prob.m3 <- sapply(dat.1$time,
                        function(x){
                          m.3(max.3.n$estimate, x)})
max.b.n <- nlm(loglike.b, 
               data=c(dat.1$frq_obs.m1, dat.1$frq_nobs.m1),
               t=dat.1$time, p=c(0.8, 0.99, 0.2))
dat.1$prob.mb <- sapply(dat.1$time,
                        function(x){
                          m.b(max.b.n$estimate, x)})

dat.1$frq_prd.m1 <- sapply(dat.1$prob.m1, 
                           function(x){
                             rbinom(1, n, x)
                           })
dat.1$frq_prd.m2 <- sapply(dat.1$prob.m2, 
                           function(x){
                             rbinom(1, n, x)
                           })
dat.1$frq_prd.m3 <- sapply(dat.1$prob.m3, 
                           function(x){
                             rbinom(1, n, x)
                           })
dat.1$frq_prd.mb <- sapply(dat.1$prob.mb, 
                           function(x){
                             rbinom(1, n, x)
                           })

###########################################
## Continue with caution
### There seem to be some kind of error in the computation of the SSE and 
### espacially on SST (old)
###########################################

## RMSE
# from Tutorial on maximum likelihood estimation
# In Jae Myung
dat.1$SE.m1 <- (dat.1$frq_obs.m1 - dat.1$frq_prd.m1)^2
dat.1$SE.m2 <- (dat.1$frq_obs.m1 - dat.1$frq_prd.m2)^2
dat.1$SE.m3 <- (dat.1$frq_obs.m1 - dat.1$frq_prd.m3)^2
dat.1$SE.mb <- (dat.1$frq_obs.m1 - dat.1$frq_prd.mb)^2

RMSE.1 <- data.frame(Model=c(1,2,3, "b"),
                   SSE=c( sum(dat.1$SE.m1),
                          sum(dat.1$SE.m2),
                          sum(dat.1$SE.m3),
                          sum(dat.1$SE.mb)),
                   RMSE=c( sqrt(SSE.m1.1/length(dat.1$frq_obs.m1)),
                           sqrt(SSE.m2.1/length(dat.1$frq_obs.m1)),
                           sqrt(SSE.m3.1/length(dat.1$frq_obs.m1)),
                           sqrt(SSE.mb.1/length(dat.1$frq_obs.m1))),
                   SST=rep(sum((dat.1$frq_obs.m1 - 
                                  mean(dat.1$frq_obs.m1))^2), 4))


# #SST.m1.1 <- sum((dat.1$frq_obs - dat.1$m_distr)^2) 
# # ?, seems strange.
# #
# # What about (yep, this is the way to go, Myung,)
#SST.m1.1 <- sum((dat.1$frq_obs.m1 - mean(dat.1$frq_obs.m1))^2)
#SST.m2.1 <- sum((dat.1$frq_obs.m1 - mean(dat.1$frq_obs.m1))^2)
#SST.m3.1 <- sum((dat.1$frq_obs.m1 - mean(dat.1$frq_obs.m1))^2)
#SST.mb.1 <- sum((dat.1$frq_obs.m1 - mean(dat.1$frq_obs.m1))^2)

#PVAF:
PVAF.1 <- data.frame(Model=c(1,2,3,"b"),
                     PVAF=c((1 - (SSE.m1.1/SST.m1.1)),
                            (1 - (SSE.m2.1/SST.m2.1)),
                            (1 - (SSE.m3.1/SST.m3.1)),
                            (1 - (SSE.mb.1/SST.mb.1))))


## AIC, BIC
b.aic <- function(loglike, nparam){
  -2*loglike + 2*nparam
}
b.bic <- function(loglike, nparam, sample){
  -2*loglike + nparam*log(sample)
}


aic.1 <- data.frame(Model=c(1,2,3),
                    AIC=c(b.aic(loglike.1(param=max.1$estimate, 
                                          data=c(dat.1$frq_obs.m1,
                                                 dat.1$frq_nobs.m1),
                                          t=dat.1$time),
                                length(max.1$estimate)),
                          b.aic(loglike.2(param=max.2$estimate, 
                                          data=c(dat.1$frq_obs.m1,
                                                 dat.1$frq_nobs.m1),
                                          t=dat.1$time),
                                length(max.2$estimate)),
                          b.aic(loglike.3(param=max.3$estimate,
                                          data=c(dat.1$frq_obs.m1,
                                                 dat.1$frq_nobs.m1),
                                          t=dat.1$time),
                                length(max.3$estimate))),
                    BIC=c(b.bic(loglike.1(param=max.1$estimate, 
                                          data=c(dat.1$frq_obs.m1,
                                                 dat.1$frq_nobs.m1),
                                          t=dat.1$time),
                                length(max.1$estimate),
                                length(dat.1$frq_obs.m1)),
                          b.bic(loglike.2(param=max.2$estimate, 
                                          data=c(dat.1$frq_obs.m1,
                                                 dat.1$frq_nobs.m1),
                                          t=dat.1$time),
                                length(max.2$estimate),
                                length(dat.1$frq_obs.m1)),
                          b.bic(loglike.3(param=max.3$estimate,
                                          data=c(dat.1$frq_obs.m1,
                                                 dat.1$frq_nobs.m1),
                                          t=dat.1$time),
                                length(max.3$estimate),
                                length(dat.1$frq_nobs.m1))))
