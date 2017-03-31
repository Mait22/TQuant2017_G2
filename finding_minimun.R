##### Fidning minimum and couting upwards
### Use this reshape for the wide format
dat1.l <- reshape(dat.wide, direction="long",
                  varying = list(2:dim(dat1)[2]))
names(dat1.l)[names(dat1.l) %in% "run"] <- "sample"
names(dat1.l)[names(dat1.l) %in% "bla1.1"] <- "y"
names(dat1.l)[names(dat1.l) %in% "time"] <- "fittedmodel"
dat1.l$model <- dat1.l$fittedmodel %% 3
dat1.l$model[dat1.l$model==0] <- 3
dat1.l$data <- 0
dat1.l$data[dat1.l$fittedmodel %in% c(1,2,3)] <- 1
dat1.l$data[dat1.l$fittedmodel %in% c(4,5,6)] <- 2
dat1.l$data[dat1.l$fittedmodel %in% c(7,8,9)] <- 3

dat1.l <- dat1.l[order(dat1.l$data, dat1.l$sample),]

### Long data
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

