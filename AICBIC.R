aic.new <- function(loglike, nparam){
  -2*loglike + 2*nparam
}
bic.new <- function(loglike, nparam, sample){
  -2*loglike + nparam*log(sample)
}

rmse <- function(obs, pred){
  sqrt(sum((obs - pred)^2) / length(obs))
}

pvaf <- function(obs, pred){
  1 - (sum((obs - pred)^2) / sum((obs - mean(obs))^2))
}
