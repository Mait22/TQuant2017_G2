aic.new <- function(loglike, nparam){
  -2*loglike + 2*nparam
}
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
