setOldClass("lm")

setClass("polynomial_fitter",
         poly.model = "lm",
         coefs = "numeric"
         )



polynomial_fitter <- function(order, data, response.var, pred.var){
  
  model <- lm(as.formula(paste(response.var,"~","poly(",pred.var,",",order,")",sep = "")),data = data)
  
  return(model)
  
}




#TEST  
test <- polynomial_fitter(order = 3,data = Chocolate, response.var = "Chocolate", pred.var = "Anxiety")

