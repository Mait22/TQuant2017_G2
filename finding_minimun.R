##### Fidning minimum and couting upwards

dat <- expand.grid(data=c(1,2,3),
                   model=c(1,2,3),
                   sample=1:10)
# Order by data and sample
dat <- dat[order(dat$data, dat$sample),]
dat$y <- rbinom(90, 2, .5)

mat.min <- matrix(numeric(9), ncol=3, nrow=3)

row.names(mat.min) <- unique(dat$model)
colnames(mat.min) <- unique(dat$data)

get.min <- function(data, matr){
  for(i in unique(dat$data)){
    for(j in unique(dat$sample)){
      slct <- dat[dat$data==i & dat$sample==j, c("model", "y")]
      matr[, i] <- matr[, i] + slct$y %in%  min(slct$y)
    }
  }
  
  matr.marg <- addmargins(matr)
  
  matr.min <- round(matr/matr.marg[nrow(matr.marg),
                                          ncol(matr.marg)-1], 2)
  return(matr.min)
}

