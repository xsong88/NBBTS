######################################################################################################################
# Load libraries and set up cluster for parallel processing
######################################################################################################################
# set seed for simulations 
seed <- 1
file.dir <- NULL

library(foreach)
library(doParallel)
library(parallel)
library(splines)
library(boot)

no_cores <- 5

# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#make sure packages are on each core
clusterEvalQ(cl, library(boot))
clusterEvalQ(cl, library(splines))

######################################################################################################################
# RUN NAIVE AND NONPARAMETRIC SIMULATIONS
######################################################################################################################

# functions needed for naive and nonparametric simulations

# SIMEX dataset creation function
SIMEX.data <- function(data,simex.var.pos,omega,B1){
  counter <- 0
  lambda1 = c(0,seq(1/8,2,length.out=16)) 
  
  W.lambda1.B1 <- matrix(ncol=3,nrow=nrow(data)*length(lambda1)*B1)
  for (i in 1:length(lambda1)) {
    for (j in 1:B1) {
      U <- rnorm(n = nrow(data),mean=0, sd=sqrt(omega))                    
      epsilon <- matrix(U, ncol = 1, nrow = nrow(data))
      marker <- matrix(data[,simex.var.pos],ncol=1)
      new.var = marker + (sqrt(lambda1[i])*epsilon)
      
      # save generated values 
      n <- length(new.var)
      W.lambda1.B1[((counter*n)+1):((counter*n)+n),1] <- new.var
      W.lambda1.B1[((counter*n)+1):((counter*n)+n),2] <- rep(lambda1[i],times=n)
      W.lambda1.B1[((counter*n)+1):((counter*n)+n),3] <- rep(j,times=length(new.var))
      counter <- counter + 1
    }
  } 
  colnames(W.lambda1.B1) <- c("gen.w","lambda2","B2")
  return(W.lambda1.B1)
}

# SIMEX function for coefficients
my.simex1 <- function(data,SIMEXdat, B1,int, l, r){ # estimates, 
  lambda1 = c(0,seq(1/8,2,length.out=16)) 
  n=nrow(data)
  n.coef = 2*(length(int)+3)
  simex.est1 <- NULL 
  s1 <- foreach (i=1:length(lambda1)) %dopar% {
    my.vals <- matrix(NA,B1,n.coef)
    for(j in 1:B1){
      Y <- data$Y
      trt <- data$trt
      new.var1 <- SIMEXdat[(i-1)*n*B1+(j-1)*n+(1:n),1]
      zz <- bs(new.var1,knots=int, Boundary.knots=c(l,r), degree=2)
      new.vals <- glm(Y~zz + (trt*zz), family="binomial")
      v <- as.numeric(summary(new.vals)$coefficients[,1])
      my.vals[j,] <- v
      
    }
    rbind(simex.est1,colMeans(my.vals,na.rm=TRUE))
  }
  simex.est1 <- matrix(unlist(s1), ncol=length(s1[[1]]), byrow=T)
  #####
  
  ext2 <- lm(simex.est1 ~ lambda1 + I(lambda1^2))
  my.est2 <- predict(ext2, newdata = data.frame(lambda1 = -1))
  
  return(my.est2)
}

# SIMEX function for theta
my.simex2 <- function(naive, data, SIMEXdat, B2, int, l, r){ # estimates, 
  
  lambda = c(0,seq(1/8,2,length.out=16)) 
  k=length(lambda)
  n=nrow(data)
  simex.est <- NULL
  simex.est.quad <- NULL
  w.vars = SIMEXdat
  d1 = cbind(data,NA)
  s2 <- foreach(i=1:length(lambda)) %dopar% {
    my.thetas <- NULL
    # my.thetas.quad <- NULL
    my.thetas.quad <- rep(NA,B2)
    for (j in 1:B2) {
      # new.dat <-  w.vars[w.vars[,2]==lambda[i] & w.vars[,3]==j,]
      new.dat <-  w.vars[(i-1)*n*B2+(j-1)*n+(1:n),1]
      d1[,ncol(d1)] <- new.dat
      new.val.quad <- theta.fcn(est.coeffs=naive, data=d1,marker.pos=ncol(d1), int=int, l=l, r=r)
      # my.thetas.quad <- rbind(my.thetas.quad, new.val.quad)
      my.thetas.quad[j] <- new.val.quad
    }
    bb <- rbind(simex.est.quad,mean(my.thetas.quad,na.rm=TRUE))
    return(bb)
  }
  simex.est.all <- matrix(unlist(s2), ncol=length(s2[[1]]), byrow=T)
  simex.est.quad <- matrix(simex.est.all[,1],ncol=1)
  
  
  ext.quad <- lm(simex.est.quad ~ lambda + I(lambda^2))
  my.est.quad <- predict(ext.quad, newdata = data.frame(lambda = -1))
  
  my.est.2 <- my.est.quad
  return(my.est.2)
}

# function to get empirical theta estimates
theta.fcn <- function(est.coeffs,data,marker.pos, int, l, r){
  marker <- data[,marker.pos]
  colnames(data) <- c('Y','trt','marker')
  data=data.frame(data)
  data.0 <- data
  data.0$trt <- 0
  data.1 <- data
  data.1$trt= 1
  b0 <- bs(marker,knots=int, Boundary.knots=c(l,r),degree=2)
  n.coef = length(est.coeffs)/2
  
  eta <- est.coeffs[n.coef+1:n.coef]
  delta.hat <-  -(eta[1]+b0%*%eta[2:n.coef])
  
  marker.neg<-ifelse(delta.hat<0, 1,0)
  p.neg <- mean(delta.hat<0)
  data=cbind(data,marker.neg)
  
  #empirical estimate
  B.neg.emp <- ifelse(length(data$Y[data$trt==1 & data$marker.neg==1])>0 & length(data$Y[data$trt==0 & data$marker.neg==1])>0,
                      mean(data$Y[data$trt==1 & data$marker.neg==1]) - mean(data$Y[data$trt==0 & data$marker.neg==1]),0)
  theta.empirical <- B.neg.emp*p.neg
  
  return(theta.empirical)
}


# simulation to bootstrap
boot.func <- function(data, ind, data.reprod,b, kn, lb,rb, kn.reprod,sdat,sigma.mom,is.original=FALSE){
  library(splines)
  library(doParallel)
  library(foreach)
  
  n <- length(ind)
  n.rep <- nrow(data.reprod)
  if(sum(ind==1:n)==n)
  {
    index <- 1:n.rep
  }else
  {
    index <- sample(1:n.rep,n.rep,replace=TRUE)
  }
  dat <- data[ind,]

  ########################################################################################################
  ## G FUNCTION AND ME ESTIMATE
  ########################################################################################################
  
  options(warn=1)
  
  f1 <- glm(data.reprod[index,1] ~bs(data.reprod[index,2], degree=2, knots=kn.reprod, Boundary.knots=c(min(data.reprod[,2]), max(data.reprod[,2]))),family="gaussian")
  vals <- as.numeric(f1$fitted.values)

  # ESTIMATE MEASUREMENT ERROR
  x <- data.reprod[index,1]-vals
  sigma.mom <- ((sum(x^2))/length(x)) - (mean(x))^2
  
  ########################################################################################################
  ## SIMEX DATASET
  ########################################################################################################
  
  
  #generate SIMEX dataset to use in SIMEX's and to select boundary knots
  if(!is.original)
  {
    sdat <- SIMEX.data(data=dat,simex.var.pos=3,omega=sigma.mom,B1=b)
  }
  
  
  ########################################################################################################
  ## NAIVE MODEL
  ########################################################################################################
  
  options(warn=1)
  hh <- bs(dat$obs_marker, knots=kn, Boundary.knots=c(min(data$obs_marker),max(data$obs_marker)), degree=2)
  Nmod0 <- glm(dat$Y~hh + (dat$trt*hh), family="binomial") 
  Np0 <- as.numeric(summary(Nmod0)$coefficients[,1])
  
  k.naive <- rep(NA,100)
  k.naive[1:length(Np0)] <- Np0
  
  ########################################################################################################
  ## CORRECTED
  ########################################################################################################
  
  estimated.vals <- my.simex1(data=dat, SIMEXdat=sdat, B1=b, int=kn, l=lb, r=rb)
  
  k.corr.quad <- rep(NA,100)
  k.corr.quad[1:length(estimated.vals)] <- estimated.vals
  
  ########################################################################################################
  ## THETAS
  ########################################################################################################
  # get naive  theta estimates
  naive.theta <- theta.fcn(est.coeffs=Np0,data=dat,marker.pos=3, int=kn,l=min(data[,3]), r=max(data[,3]))
  
  # get SIMEX corrected theta estimates 
  corrected.theta.quad.only <- my.simex2(naive=estimated.vals, data=dat, SIMEXdat=sdat,B2=b, int=kn, l=lb,r=rb)
  
  
  ########################################################################################################
  ## OUTPUT
  ########################################################################################################
  
  out1 <- c(naive.theta, corrected.theta.quad.only, lb, rb, length(kn), k.naive, k.corr.quad)
  
  return(out1)
}



# set up parameters, simulation seed, and ensure needed functions are on each core

b=200 # number of SIMEX datasets
rf=100 # number of bootstrap iterations # rf=100 # number of bootstrap iterations
set.seed(seed, kind = "L'Ecuyer-CMRG")
clusterExport(cl, c("SIMEX.data","my.simex1","my.simex2","theta.fcn","boot.func"))


# read in reproducibility data as dat.reprod
#ILLUMINA PLATFORM IS UNOBSERVED X
#AFFY IS observed Z 
dat.reprod <- read.csv(paste0(file.dir,'Reproducibility Data.csv'))
dat.reprod <- data.frame(dat.reprod)

# get measurement error estimate between the two platforms
options(warn=2)
model.AIC <- NULL
for (a in 2:10){
  m <- try(glm(Affy~bs(Illumina,df=a,degree=2),family="gaussian", data=dat.reprod),silent=T)
  if("try-error" %in% class(m)) m <- 100000000000000000000000
  else m <- m$aic
  
  model.AIC <- c(model.AIC, m)
}
mat <- cbind(rep(2:10),model.AIC)
mat <- mat[order(mat[,2],decreasing=F),]
n.rep <- nrow(dat.reprod)
df.g <- mat[1,1]
g.knots <- quantile(dat.reprod[,2], probs=seq(0,1,length.out=df.g), na.rm=T)
if (df.g<=2){knots.g=NULL
}else {knots.g=as.numeric(g.knots[2:(length(g.knots)-1)])}
# knots.g are the internal knots for the g function...


options(warn=1)
f1 <- glm(dat.reprod[,1]~bs(dat.reprod[,2],df=mat[1,1], knots=knots.g,degree=2),family="gaussian") ##7  
vals <- as.numeric(f1$fitted.values)
dd <- cbind(dat.reprod[,2], vals)
dd <- dd[order(dd[,1]),]
# linear model to compare
mod3 <- lm(dat.reprod[,1]~dat.reprod[,2])
vals3 <- as.numeric(mod3$fitted.values)
dd3 <- cbind(dat.reprod[,2], vals3)
dd3 <- dd3[order(dd3[,1]),]


# estimate of the measurement error
extended.bspline <- vals-(dat.reprod[,1])
sigma.mom1 <- ((sum(extended.bspline^2))/length(extended.bspline)) - (mean(extended.bspline))^2
cat("=== Estimate sigma2 =======\n")
cat("sigma2 = ", sigma.mom1, "\n\n")


cat("=== Estimate g function========\n")
cat("knots:", c(min(dat.reprod[,2]),knots.g,max(dat.reprod[,2])),"\n")
cat("gamma:", f1$coefficients,"\n\n")

# read in clinical data file
dat.clinic <- read.csv(paste0(file.dir,'Clinical Data.csv'))

#generate SIMEX dataset to use in SIMEX's and to select boundary knots

sdat <- SIMEX.data(data=dat.clinic,simex.var.pos=3,omega=sigma.mom1,B1=b)

#select boundary knots for the corrected models
lb <- min(sdat[,1])
rb <- max(sdat[,1])
n <- nrow(dat.clinic)

# determine number of knots for all models based on the naive models

options(warn=2)
naive1.AIC <- NULL
# naive2.AIC <- NULL
for (j in 2:10){
  gg <- bs(dat.clinic$obs_marker,df=j, degree=2, Boundary.knots = c(min(dat.clinic$obs_marker),max(dat.clinic$obs_marker)))
  m1 <- try(glm(dat.clinic$Y~gg + (dat.clinic$trt*gg),family="binomial")$aic,silent=T)
  if("try-error" %in% class(m1)) m1 <- 1000000000000000000000000
  naive1.AIC <- c(naive1.AIC, m1)
}
mat.naive <- cbind(rep(2:10),naive1.AIC) # , naive2.AIC
mat.naive <- mat.naive[order(mat.naive[,2],decreasing=F),]

df.naive <- mat.naive[1,1]
naive.knots <- quantile(dat.clinic$obs_marker, probs=seq(0,1,length.out=df.naive), na.rm=T)
if (df.naive<=2){knots.naive=NULL
}else {knots.naive=as.numeric(naive.knots[2:(length(naive.knots)-1)])}
# knots.naive are the internal knots for the naive method and the corrected method.


# Get Estimate
est.original <- boot.func(data=dat.clinic, ind=1:n, data.reprod=dat.reprod, b=b,kn=knots.naive, lb=lb, rb=rb, kn.reprod=knots.g, sdat=sdat,sigma.mom = sigma.mom1,is.original=TRUE)
boundary.points <- c(range(dat.clinic[,3]),lb,rb)
method <- c("naive","SIMEX")
cat("=== Estimate h* function========\n")
index <- 5
for(r in 1:2)
{
  beta.est <- (est.original[index+(1:100)])[!(is.na(est.original[index+(1:100)]))]
  n.beta <- length(beta.est)/2
  beta.est[n.beta+(1:n.beta)] <- beta.est[1:n.beta]+beta.est[n.beta+(1:n.beta)]
  n.beta <- length(beta.est)/2
  beta.est[n.beta+(1:n.beta)] <- beta.est[1:n.beta]+beta.est[n.beta+(1:n.beta)]
  cat("==",method[r],":\n")
  cat("knots:", c(boundary.points[(r-1)*2+1],knots.naive,boundary.points[(r-1)*2+2]),"\n")
  cat("beta:", beta.est,"\n\n")
  index <- index+100
}
cat("\n")

# run bootstrap
boot.results <- suppressWarnings(boot(data=dat.clinic, statistic=boot.func,data.reprod=newdat,R=rf, b=b, kn=knots.naive, lb=lb, rb=rb, kn.reprod=knots.g,sdat=sdat,sigma.mom = sigma.mom1))
stopCluster(cl)
results100011 <- rbind(est.original, boot.results$t)

# output compilation for naive and nonparameteric scenarios
rr <- results100011[, c(1:11, 106:111)]
mad.list <- apply(rr[-1,],2,mad) # median absolute deviation 

# getting estimates based on full dataset
top.rows2 <- rr[1,]

#bootstrap MAD coverage
boot.l=top.rows2-(qnorm(0.975)*mad.list)
boot.u=top.rows2+(qnorm(0.975)*mad.list)
s2.boot1=cbind(boot.l,boot.u)

result <- cbind(top.rows2[1:2],mad.list[1:2],s2.boot1[1:2,])
colnames(result) <- c("Theta.est","se","CI.lower","CI.upper")
rownames(result) <- c("naive","SIMEX")
cat("=== Estimate Theta ========\n")
print(result)
