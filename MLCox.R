

# LOAD PACKAGES

detach("package:mgcv",unload=T)
install.packages("SuperLearner")
install.packages("glmnet")
install.packages("ranger")
install.packages("grf")
library(survival)
library(glmnet)
library(grf)
library(mgcv)
library(ranger)
library(MASS)
library(riskRegression)
library(pec)
library(SuperLearner)
library(gam)

g=function(x){
  res=log(-log(x))
  return(res)
} # Cox-link
g.prime=function(x){1/(x*log(x))}
logit=function(x){log(x/(1-x))} 
expit=function(x){exp(x)/(exp(x)+1)}

# ESTIMATOR USING SURVIVAL RANDOM FORESTS, WITH DICHOTOMOUS EXPOSURE

estimator<-function(dataml,dataest){
  
  sxtime<-c(0,sort(dataest$xtime))
  x = sort(sxtime) 

  fit_rf = survival_forest(X=dataml[,1:(ncol(dataml)-2)],Y=dataml$xtime,D=dataml$status,failure.times = sxtime,num.trees = 2000)
  S = predict(fit_rf,newdata = dataest[,1:(ncol(dataml)-2)])$predictions
  dataml.0 = data.frame(a=0,l=dataest[,2:(ncol(dataml)-2)])
  S0 = predict(fit_rf,newdata = dataml.0)$predictions
  dataml.1 = data.frame(a=1,l=dataest[,2:(ncol(dataml)-2)])
  S1 = predict(fit_rf,newdata = dataml.1)$predictions

  fit_rf = survival_forest(X=dataml[,1:(ncol(dataml)-2)],Y=dataml$xtime,D=1-dataml$status,failure.times = sxtime,num.trees = 2000)
  C = predict(fit_rf,newdata = dataest[,1:(ncol(dataml)-2)])$predictions
  C[,ncol(C)] = C[,ncol(C)-1]
  n=nrow(dataest)

  fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=dataml$a,family=binomial(),SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
  p0 = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)
  q0 = g(S1)*p0 + g(S0)*(1-p0)

# Components

  k = length(x)-1
  delta.tau=c(0,x[-1]-x[1:k])
  delta.tau[(1:(k+1))[apply(S,2,max)>exp(-1/(0.1*n))]]=0
  delta.tau[(1:(k+1))[apply(S*C,2,min)<(1/(0.1*n))]]=0

  dN = matrix(0,ncol=k+1,nrow=n)
  dNc = matrix(0,ncol=ncol(dN),nrow=nrow(dN))
  R = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  Rc = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  s = sort.list(dataest$xtime)
  for (i in 1:n){
    dN[s[i],i+1] = dataest$status[s[i]]
    dNc[s[i],i+1] = 1-dataest$status[s[i]]
    if ((i+2) <= ncol(R)){
      R[s[i],(i+2):ncol(R)] = 1-dataest$status[s[i]]
      Rc[s[i],(i+2):ncol(R)] = dataest$status[s[i]]
    }
  }

# Time-fixed estimand

# Plug-in

  h = (dataest$a-p0)*(g(S)-q0)
  h[is.na(h)] = 0
  h[is.infinite(h)] = 0
  hd = (dataest$a-p0)^2
  gamma = mean(hd*sum(delta.tau))
  theta = mean(apply(t(t(h)*delta.tau),1,sum))/gamma
  plug = theta
  IF1 = apply(t(t(h)*delta.tau),1,sum)/gamma
  IF2 = hd*sum(delta.tau)/gamma

# Augmentation S

  hs = -(dataest$a-p0)*g.prime(S)*S
  hs[is.na(hs)] = 0
  hs[is.infinite(hs)] = 0
  hs[,(1:(k+1))[apply(S,2,max)==1]]=0

  dM=R
  dM.int=dM
  dLam=-log(S[,1])
  tmp=1/(S[,1]*C[,1])
  dM[,1]=tmp*(dN[,1]-R[,1]*Rc[,1]*dLam)
  dM.int[,1]=dM[,1]
  for (j in 2:ncol(dM)){
    dLam=-log(S[,j])+log(S[,j-1])
    tmp=1/(S[,j]*C[,j])
    dM[,j]=tmp*(dN[,j]-R[,j]*Rc[,j]*dLam)
    dM.int[,j]=dM.int[,j-1]+dM[,j]
  }
  hs = hs*dM.int
  IFS = apply(t(t(hs)*delta.tau),1,sum)/gamma

  list(IF1=IF1,IF2=IF2,IFS=IFS)
}

# ESTIMATOR USING SURVIVAL SUPERLEARNER, WITH DICHOTOMOUS EXPOSURE

install.packages("devtools")
library(devtools)
devtools::install_github("tedwestling/survSuperLearner",force=TRUE)
library(survSuperLearner)

estimator<-function(dataml,dataest){
  
  sxtime<-c(0,sort(dataest$xtime))
  x = sort(sxtime) 
  
  event.SL.library <- cens.SL.library <- c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc")
  
  fit_rf <- survSuperLearner(time = dataml$xtime, event = dataml$status, X = dataml[,1:(ncol(dataml)-2)], newX = dataest[,1:(ncol(dataml)-2)], new.times = sxtime, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = TRUE)
  S = fit_rf$event.SL.predict
  dataml.0 = data.frame(a=0,l=dataest[,2:(ncol(dataml)-2)])
  colnames(dataml.0) = colnames(dataml[,1:(ncol(dataml)-2)])
  fit_rf0 = predict.survSuperLearner(fit_rf,newdata=dataml.0,new.times = sxtime,onlySL=TRUE)
  S0 = fit_rf0$event.SL.predict
  dataml.1 = data.frame(a=1,l=dataest[,2:(ncol(dataml)-2)])
  colnames(dataml.1) = colnames(dataml[,1:(ncol(dataml)-2)])
  fit_rf1 = predict.survSuperLearner(fit_rf,newdata=dataml.1,new.times = sxtime,onlySL=TRUE)
  S1 = fit_rf1$event.SL.predict
  
  C = fit_rf$cens.SL.predict
  C[,ncol(C)] = C[,ncol(C)-1]
  n=nrow(dataest)
  
  fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=dataml$a,family=binomial(),SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
  p0 = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)
  q0 = g(S1)*p0 + g(S0)*(1-p0)
  
  # Components
  
  k = length(x)-1
  delta.tau=c(0,x[-1]-x[1:k])
  delta.tau[(1:(k+1))[apply(S,2,max)>exp(-1/(0.1*n))]]=0
  delta.tau[(1:(k+1))[apply(S*C,2,min)<(1/(0.1*n))]]=0
  
  dN = matrix(0,ncol=k+1,nrow=n)
  dNc = matrix(0,ncol=ncol(dN),nrow=nrow(dN))
  R = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  Rc = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  s = sort.list(dataest$xtime)
  for (i in 1:n){
    dN[s[i],i+1] = dataest$status[s[i]]
    dNc[s[i],i+1] = 1-dataest$status[s[i]]
    if ((i+2) <= ncol(R)){
      R[s[i],(i+2):ncol(R)] = 1-dataest$status[s[i]]
      Rc[s[i],(i+2):ncol(R)] = dataest$status[s[i]]
    }
  }
  
  # Plug-in
  
  h = (dataest$a-p0)*(g(S)-q0)
  h[is.na(h)] = 0
  h[is.infinite(h)] = 0
  hd = (dataest$a-p0)^2
  gamma = mean(hd*sum(delta.tau))
  theta = mean(apply(t(t(h)*delta.tau),1,sum))/gamma
  plug = theta
  IF1 = apply(t(t(h)*delta.tau),1,sum)/gamma
  IF2 = hd*sum(delta.tau)/gamma
  
  # Augmentation S
  
  hs = -(dataest$a-p0)*g.prime(S)*S
  hs[is.na(hs)] = 0
  hs[is.infinite(hs)] = 0
  hs[,(1:(k+1))[apply(S,2,max)==1]]=0
  
  dM=R
  dM.int=dM
  dLam=-log(S[,1])
  tmp=1/(S[,1]*C[,1])
  dM[,1]=tmp*(dN[,1]-R[,1]*Rc[,1]*dLam)
  dM.int[,1]=dM[,1]
  for (j in 2:ncol(dM)){
    dLam=-log(S[,j])+log(S[,j-1])
    tmp=1/(S[,j]*C[,j])
    dM[,j]=tmp*(dN[,j]-R[,j]*Rc[,j]*dLam)
    dM.int[,j]=dM.int[,j-1]+dM[,j]
  }
  hs = hs*dM.int
  IFS = apply(t(t(hs)*delta.tau),1,sum)/gamma
  
  list(IF1=IF1,IF2=IF2,IFS=IFS)
}

### ESTIMATOR USING VARIABLE SELECTION, WITH CONTINUOUS EXPOSURE

estimator.cont<-function(dataml,dataest){
  
  sxtime<-c(0,sort(dataest$xtime))
  x = sort(sxtime) 
  
  cv.lasso = cv.glmnet(as.matrix(dataml[,1:(ncol(dataml)-2)]), Surv(dataml$xtime,dataml$status), alpha = 1, family = "cox", nfolds=20)
  model_cox = glmnet(as.matrix(dataml[,1:(ncol(dataml)-2)]), Surv(dataml$xtime,dataml$status), alpha = 1, family = "cox", lambda = cv.lasso$lambda.min)
  selection1 <- which(c(1,model_cox$beta[2:(ncol(dataml)-2)] != 0) != 0) ## always includes the treatment
  
  lhs = "Surv(dataml$xtime, dataml$status)"
  if (selection1[1] == 1){
    rhs = paste(c("a",paste(rep("l.l",length=length(selection1)-1),selection1[-1]-1,sep="")),collapse=" + ")
  }
  if (selection1[1] != 1){
    rhs = paste(paste(rep("l.l",length=length(selection1)),selection1,sep=""),collapse=" + ")
  }
  form = as.formula(paste(lhs, "~", rhs))
  fit_rf = coxph(form, x=T, data = as.data.frame(dataml[,1:(ncol(dataml)-2)]))
  S = predictSurvProb(fit_rf,newdata = data.frame(dataest[,1:(ncol(dataml)-2)]),times = sxtime)
  model_coxf = fit_rf
  
  ### Estimate censoring chances 
  
  cv.lasso = cv.glmnet(as.matrix(dataml[,1:(ncol(dataml)-2)]), Surv(dataml$xtime,1-dataml$status), alpha = 1, family = "cox", nfolds=20)
  model_cox = glmnet(as.matrix(dataml[,1:(ncol(dataml)-2)]), Surv(dataml$xtime,1-dataml$status), alpha = 1, family = "cox", lambda = cv.lasso$lambda.min)
  selection2 <- which(c(1,model_cox$beta[2:(ncol(dataml)-2)] != 0) != 0) ## always includes the treatment
  lhs = "Surv(dataml$xtime, dataml$status)"
  if (selection2[1] == 1){
    rhs = paste(c("a",paste(rep("l.l",length=length(selection2)-1),selection2[-1]-1,sep="")),collapse=" + ")
  }
  if (selection2[1] != 1){
    rhs = paste(paste(rep("l.l",length=length(selection2)),selection2,sep=""),collapse=" + ")
  }
  form = as.formula(paste(lhs, "~", rhs))
  fit_rf = coxph(form, x=T, data = as.data.frame(dataml[,1:(ncol(dataml)-2)]))
  C = predictSurvProb(fit_rf,newdata = data.frame(dataest[,1:(ncol(dataml)-2)]),times = sxtime)
  
  n=nrow(dataest)
  
  fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=dataml$a,SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
  p0 = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)
  
  # Components
  
  k = length(x)-1
  delta.tau=c(0,x[-1]-x[1:k])
  delta.tau[(1:(k+1))[apply(S,2,max)>exp(-1/(0.1*n))]]=0
  delta.tau[(1:(k+1))[apply(S*C,2,min)<(1/(0.1*n))]]=0
  #delta.tau[1:(0.05*k+1)]=0
  
  dN = matrix(0,ncol=k+1,nrow=n)
  dNc = matrix(0,ncol=ncol(dN),nrow=nrow(dN))
  R = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  Rc = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  s = sort.list(dataest$xtime)
  for (i in 1:n){
    dN[s[i],i+1] = dataest$status[s[i]]
    dNc[s[i],i+1] = 1-dataest$status[s[i]]
    if ((i+2) <= ncol(R)){
      R[s[i],(i+2):ncol(R)] = 1-dataest$status[s[i]]
      Rc[s[i],(i+2):ncol(R)] = dataest$status[s[i]]
    }
  }
  
  # Plug-in
  
  hd = (dataest$a-p0)^2
  gamma = mean(hd*sum(delta.tau))
  theta = model_coxf$coef[1]
  IF = (dataest$a-p0)^2*sum(delta.tau)/gamma
  
  # Augmentation S
  
  hs = -(dataest$a-p0)*g.prime(S)*S
  hs[is.na(hs)] = 0
  hs[is.infinite(hs)] = 0
  hs[,(1:(k+1))[apply(S,2,max)==1]]=0
  #hs[,1:(0.05*k+1)]=0
  
  dM=R
  dM.int=dM
  dLam=-log(S[,1])
  tmp=1/(S[,1]*C[,1])
  dM[,1]=tmp*(dN[,1]-R[,1]*Rc[,1]*dLam)
  dM.int[,1]=dM[,1]
  for (j in 2:ncol(dM)){#j=ncol(dM)
    dLam=-log(S[,j])+log(S[,j-1])
    tmp=1/(S[,j]*C[,j])
    dM[,j]=tmp*(dN[,j]-R[,j]*Rc[,j]*dLam)
    dM.int[,j]=dM.int[,j-1]+dM[,j]
  }
  hs = hs*dM.int
  IFS = apply(t(t(hs)*delta.tau),1,sum)/gamma
  
  list(IF=IF,IFS=IFS,theta=theta)
}

### SIMULATION STUDY USING MACHINE LEARNING, DICHOTOMOUS EXPOSURE

ndat = 500
nsim = 1000
Sigma = toeplitz((10:1)/10)

res = matrix(0,nsim,16)

for (index in 1:nsim){
### Data-generating mechanism

  set.seed(index)

  n<-1000
  ndat <- n
  l = mvrnorm(n,rep(0,10),Sigma)
  u1 = sqrt(abs(l[,1]*l[,2]))-sqrt(abs(l[,10]))+cos(l[,5])-cos(l[,5])*cos(l[,6])
  u2 = sqrt(abs(l[,1]*l[,10]))-sqrt(abs(l[,9]))+cos(l[,5])-cos(l[,7])*cos(l[,6])
  a = rbinom(n,1,expit(-2*u1))
  weibA  <- 1.5                  # Weibull shape parameter
  weibB  <- 100                  # Weibull scale parameter
  stime <- rweibull(n,shape = weibA, scale = (weibB*exp(-a-u1/2))^(1/weibA))
  weibB  <- 125                  # Weibull scale parameter
  ctime <- rweibull(n,shape = weibA, scale = (weibB*exp(-a-u2/2))^(1/weibA))
  #ctime <- rweibull(n,shape = weibA, scale = (0.5*weibB*exp(-a-u2/2))^(1/weibA))
  #ctime <- rweibull(n,shape = weibA, scale = 1.5*(weibB*exp(-a-u2/2))^(1/weibA))
  #ctime <- rep(13,n)
  ctime <- ifelse(ctime>12,12,ctime)
  status<-ifelse(ctime<stime,0,1)
  xtime<-pmin(stime,ctime)

### Standard Cox analysis

  mod = summary(coxph(Surv(xtime,status)~a+l))$coef
  res[index,1:2] = mod[1,c(1,3)]
  mod = summary(coxph(Surv(xtime,status)~a+u1))$coef
  res[index,3:4] = mod[1,c(1,3)]
  mod = summary(coxph(Surv(xtime,status)~a+ns(l[,1],df=4)+ns(l[,2],df=4)+ns(l[,3],df=4)+ns(l[,4],df=4)+ns(l[,5],df=4)+ns(l[,6],df=4)+ns(l[,7],df=4)+ns(l[,8],df=4)+ns(l[,9],df=4)+ns(l[,10],df=4)))$coef
  res[index,5:6] = mod[1,c(1,3)]
  mod = summary(coxph(Surv(xtime,status)~a+ns(l[,1],df=3)+ns(l[,2],df=3)+ns(l[,3],df=3)+ns(l[,4],df=3)+ns(l[,5],df=3)+ns(l[,6],df=3)+ns(l[,7],df=3)+ns(l[,8],df=3)+ns(l[,9],df=3)+ns(l[,10],df=3)))$coef
  res[index,7:8] = mod[1,c(1,3)]

# Nuisance parameters

  dataml <- data.frame(a=a,l=l,xtime=xtime,status=status)  
  dataest <- data.frame(a=a,l=l,xtime=xtime,status=status)

  obj = estimator(dataml,dataest)

# Plug
  IF1 = obj$IF1
  IF2 = obj$IF2
  theta = mean(IF1)/mean(IF2)
  res[index,9] = theta
  IF = IF1 - theta*IF2
  res[index,10] = sd(IF)/sqrt(length(IF))

# Estt
  IF3 = IF1+obj$IFS
  theta = mean(IF3)/mean(IF2)
  res[index,11] = theta
  IF = IF3 - theta*IF2
  res[index,12] = sd(IF)/sqrt(length(IF))

  obj1 = estimator(dataml[-(1:(ndat/5)),],dataest[1:(ndat/5),])
  obj2 = estimator(dataml[-((ndat/5+1):(2*ndat/5)),],dataest[(ndat/5+1):(2*ndat/5),])
  obj3 = estimator(dataml[-((2*ndat/5+1):(3*ndat/5)),],dataest[(2*ndat/5+1):(3*ndat/5),])
  obj4 = estimator(dataml[-((3*ndat/5+1):(4*ndat/5)),],dataest[(3*ndat/5+1):(4*ndat/5),])
  obj5 = estimator(dataml[-((4*ndat/5+1):ndat),],dataest[(4*ndat/5+1):ndat,])

# Plug
  IF1 = c(obj1$IF1,obj2$IF1,obj3$IF1,obj4$IF1,obj5$IF1)
  IF2 = c(obj1$IF2,obj2$IF2,obj3$IF2,obj4$IF2,obj5$IF2)
  theta = mean(IF1)/mean(IF2)
  res[index,13] = theta
  IF = IF1 - theta*IF2
  res[index,14] = sd(IF)/sqrt(length(IF))

# Estt
  IF3 = IF1+c(obj1$IFS,obj2$IFS,obj3$IFS,obj4$IFS,obj5$IFS)
  theta = mean(IF3)/mean(IF2)
  res[index,15] = theta
  IF = IF3 - theta*IF2
  res[index,16] = sd(IF)/sqrt(length(IF))

# boxplot(res[1:index,1+2*(0:7)])
# abline(h=1)
# cat(index,res[index,],"\n")
  cat(index,"\n")
} 

r<-1:8
for (i in 1:8){
  a<-2*i-1
  r[i]<-mean((res[,a]-qnorm(.975)*res[,a+1]<1)&(res[,a]+qnorm(.975)*res[,a+1]>1))
}

cbind(apply(res-1,2,mean)[2*(1:8)-1],
      sqrt(n)*apply(res-1,2,mean)[2*(1:8)-1],
      apply(res,2,sd)[2*(1:8)-1],
      apply(res,2,mean)[2*(1:8)],r)

### SIMULATIONS VARIABLE SELECTION, CONTINUOUS EXPOSURE

Sigma = toeplitz((10:1)/10)
n<-200
nsim = 1000
res = matrix(0,nsim,8)

for (index in 1:1000){
  
  set.seed(index)
  
  #  n<-1000
  n<-500
  
  ndat <- n
  l = mvrnorm(n,rep(0,10),Sigma)
  u1 = sqrt(abs(l[,1]*l[,2]))-sqrt(abs(l[,10]))+cos(l[,5])-cos(l[,5])*cos(l[,6])
  u2 = sqrt(abs(l[,1]*l[,10]))-sqrt(abs(l[,9]))+cos(l[,5])-cos(l[,7])*cos(l[,6])
  a = rbinom(n,1,expit(-2*u1))
  weibA  <- 1.5                  # Weibull shape parameter
  weibB  <- 100                  # Weibull scale parameter
  stime <- rweibull(n,shape = weibA, scale = (weibB*exp(-a-u1/2))^(1/weibA))
  weibB  <- 125                  # Weibull scale parameter
  ctime <- rweibull(n,shape = weibA, scale = (weibB*exp(-a-u2/2))^(1/weibA))
  #ctime <- rweibull(n,shape = weibA, scale = (0.5*weibB*exp(-a-u2/2))^(1/weibA))
  #ctime <- rweibull(n,shape = weibA, scale = 1.5*(weibB*exp(-a-u2/2))^(1/weibA))
  ctime <- ifelse(ctime>12,12,ctime)
  status<-ifelse(ctime<stime,0,1)
  xtime<-pmin(stime,ctime)
  
  ### Standard Cox analysis
  
  mod = summary(coxph(Surv(xtime,status)~a+l))$coef
  res[index,1:2] = mod[1,c(1,3)]
  mod = summary(coxph(Surv(xtime,status)~a+u1))$coef
  res[index,3:4] = mod[1,c(1,3)]
  mod = summary(coxph(Surv(xtime,status)~a+ns(l[,1],df=4)+ns(l[,2],df=4)+ns(l[,3],df=4)+ns(l[,4],df=4)+ns(l[,5],df=4)+ns(l[,6],df=4)+ns(l[,7],df=4)+ns(l[,8],df=4)+ns(l[,9],df=4)+ns(l[,10],df=4)))$coef
  res[index,5:6] = mod[1,c(1,3)]
  mod = summary(coxph(Surv(xtime,status)~a+ns(l[,1],df=3)+ns(l[,2],df=3)+ns(l[,3],df=3)+ns(l[,4],df=3)+ns(l[,5],df=3)+ns(l[,6],df=3)+ns(l[,7],df=3)+ns(l[,8],df=3)+ns(l[,9],df=3)+ns(l[,10],df=3)))$coef
  res[index,7:8] = mod[1,c(1,3)]
  
  # Nuisance parameters
  
  dataml <- data.frame(a=a,l=l,xtime=xtime,status=status)  
  dataest <- data.frame(a=a,l=l,xtime=xtime,status=status)
  
  obj = estimator(dataml,dataest)
  
  # Plug
  IF1 = obj$IF1
  IF2 = obj$IF2
  theta = mean(IF1)/mean(IF2)
  res[index,9] = theta
  IF = IF1 - theta*IF2
  res[index,10] = sd(IF)/sqrt(length(IF))
  
  # Estt
  IF3 = IF1+obj$IFS
  theta = mean(IF3)/mean(IF2)
  res[index,11] = theta
  IF = IF3 - theta*IF2
  res[index,12] = sd(IF)/sqrt(length(IF))
  
#  boxplot(res[1:index,1+2*(0:5)])
#  abline(h=1)
  cat(index,"\n")
}

r<-1:6
for (i in 1:6){
  a<-2*i-1
  r[i]<-mean((res[,a]-qnorm(.975)*res[,a+1]<1)&(res[,a]+qnorm(.975)*res[,a+1]>1),na.rm=TRUE)
}

cbind(apply(res-0.5,2,mean,na.rm=T)[2*(1:4)-1],
      sqrt(n)*apply(res-0.5,2,mean,na.rm=T)[2*(1:4)-1],
      apply(res,2,sd,na.rm=T)[2*(1:4)-1],
      apply(res,2,mean,na.rm=T)[2*(1:4)],r)

### DATA ANALYSIS

library(survival)
data(mgus2,package="survival")
attach(mgus2)
detach(mgus)
mgus3=data.frame(xtime=futime,status=death,a=mspike,creat=creat,hgb=hgb,age=age,sex=sex,dxyr=dxyr)
mgus3=mgus3[complete.cases(mgus3),]
detach(mgus2)
attach(mgus3)
l=cbind(creat,hgb,age,sex,dxyr)

# Descriptives

par(mfrow=c(1,2))
plot(survfit(Surv(xtime,status)~1),ylab="Proportion alive",xlab="Months")
plot(survfit(Surv(xtime,1-status)~1),ylab="Proportion uncensored",xlab="Months")

survfit(Surv(xtime,1-status)~1)$surv

boxplot(a)

par(mfrow=c(1,2))
plot(survfit(Surv(mgus2$futime,mgus2$death)~1))
plot(survfit(Surv(mgus2$futime,1-mgus2$death)~1))

survfit(Surv(xtime,1-status)~1)$surv

# Variable selection

summary(coxph(Surv(xtime,status) ~ a))$coef
summary(coxph(Surv(xtime,status) ~ cbind(a,l)))$coef
cv.lasso = cv.glmnet(cbind(a,l), Surv(xtime,status), alpha = 1, family = "cox", nfolds=20)
model_cox = glmnet(cbind(a,l), Surv(xtime,status), alpha = 1, family = "cox", lambda = cv.lasso$lambda.1se)
selection1 <- which(c(1,model_cox$beta[2:6] != 0) != 0) ## always includes the treatment
mod = summary(coxph(Surv(xtime,status) ~ cbind(a,l)[,selection1]))$coef
mod[1,c(1,3)]  #Coef and its SE for Lasso

# Splines

mod = coxph(Surv(xtime,status) ~ a+ns(age,df=4)+ns(hgb,df=4)+ns(creat,df=4)+ns(dxyr,df=4))

summary(coxph(Surv(xtime,1-status) ~ cbind(a,l)))$coef
summary(lm(a~l))$coef

# Assumption-lean analysis

n=length(a)
dataml <- data.frame(a=a,l=l,xtime=xtime,status=status)  
dataest <- data.frame(a=a,l=l,xtime=xtime,status=status)

sxtime<-c(0,sort(dataest$xtime))
x = unique(sort(sxtime)) 

fit_rf = survival_forest(X=dataml[,1:(ncol(dataml)-2)],Y=dataml$xtime,D=dataml$status,failure.times = x,num.trees = 2000)
S = predict(fit_rf,newdata = dataest[,1:(ncol(dataml)-2)])$predictions

fit_rf = survival_forest(X=dataml[,1:(ncol(dataml)-2)],Y=dataml$xtime,D=1-dataml$status,failure.times = x,num.trees = 2000)
C = predict(fit_rf,newdata = dataest[,1:(ncol(dataml)-2)])$predictions

fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=dataml$a,SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
p0 = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)

q = S
for (i in 3:ncol(S)){
  fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=as.numeric(S[,i]),SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
  q[,i] = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)
  cat(i,"\t")
}

# Components

k = length(x)-1
delta.tau=c(0,x[-1]-x[1:k])
delta.tau[(1:(k+1))[apply(S,2,max)>exp(-1/(0.1*n))]]=0
delta.tau[(1:(k+1))[apply(S,2,max)>0.99]]=0
delta.tau[(1:(k+1))[apply(S*C,2,min)<(1/(0.1*n))]]=0

dN = matrix(0,ncol=k+1,nrow=n)
dNc = matrix(0,ncol=ncol(dN),nrow=nrow(dN))
R = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
Rc = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
s = sort.list(dataest$xtime)
dataest$xtime[s[i]]
for (i in 1:n){
  dN[s[i],which(dataest$xtime[s[i]]==x)] = dataest$status[s[i]]
  dNc[s[i],which(dataest$xtime[s[i]]==x)] = 1-dataest$status[s[i]]
  if ((i+2) <= ncol(R)){
    R[s[i],(i+2):ncol(R)] = 1-dataest$status[s[i]]
    Rc[s[i],(i+2):ncol(R)] = dataest$status[s[i]]
  }
}

h = (dataest$a-p0)*(g(S)-q)
h[is.na(h)] = 0
h[is.infinite(h)] = 0
hd = (dataest$a-p0)^2
gamma = mean(hd*sum(delta.tau))
theta = mean(apply(t(t(h)*delta.tau),1,sum))/gamma
plug = theta
theta
IF1 = apply(t(t(h)*delta.tau),1,sum)/gamma
IF2 = hd*sum(delta.tau)/gamma
IF = IF1 - theta*IF2
sd(IF)/sqrt(length(IF))

# Augmentation S

hs = -(dataest$a-p0)*g.prime(S)*S
hs[is.na(hs)] = 0
hs[is.infinite(hs)] = 0
hs[,(1:(k+1))[apply(S,2,max)==1]]=0

dM=R
dM.int=dM
dLam=-log(S[,1])
tmp=1/(S[,1]*C[,1])
dM[,1]=tmp*(dN[,1]-R[,1]*Rc[,1]*dLam)
dM.int[,1]=dM[,1]
for (j in 2:ncol(dM)){
  dLam=-log(S[,j])+log(S[,j-1])
  tmp=1/(S[,j]*C[,j])
  dM[,j]=tmp*(dN[,j]-R[,j]*Rc[,j]*dLam)
  dM.int[,j]=dM.int[,j-1]+dM[,j]
}
hs = hs*dM.int
IFS = apply(t(t(hs)*delta.tau),1,sum,na.rm=T)/gamma

IF3 = IF1+IFS
theta = mean(IF3)/mean(IF2)
theta
IF = IF3 - theta*IF2
sd(IF)/sqrt(length(IF))
boxplot(IFS)

dataml <- data.frame(a=a,l=l,xtime=xtime,status=status)  
dataest <- data.frame(a=a,l=l,xtime=xtime,status=status)
obj<-estimator(dataml,dataest)
theta = mean(obj$IF1)/mean(obj$IF2)
IF = obj$IF1 - theta*obj$IF2
sd(IF)/sqrt(length(IF))
IF3 = obj$IF1+obj$IFS
theta = mean(IF3)/mean(obj$IF2)
IF = IF3 - theta*obj$IF2
sd(IF)/sqrt(length(IF))

# Cross-fitting
# Code for `estimator' below
`
obj1 = estimator(dataml[-(1:268),],dataest[1:268,])
obj2 = estimator(dataml[-(269:535),],dataest[269:535,])
obj3 = estimator(dataml[-(536:803),],dataest[536:803,])
obj4 = estimator(dataml[-(804:1070),],dataest[804:1070,])
obj5 = estimator(dataml[-(1071:1338),],dataest[1071:1338,])

estimator<-function(dataml,dataest){
  
  n = nrow(dataest)
  sxtime<-c(0,sort(dataest$xtime))
  x = unique(sort(sxtime)) 
  
  fit_rf = survival_forest(X=dataml[,1:(ncol(dataml)-2)],Y=dataml$xtime,D=dataml$status,failure.times = x,num.trees = 2000)
  Sh = predict(fit_rf,newdata = dataml[,1:(ncol(dataml)-2)])$predictions
  S = predict(fit_rf,newdata = dataest[,1:(ncol(dataml)-2)])$predictions
  
  fit_rf = survival_forest(X=dataml[,1:(ncol(dataml)-2)],Y=dataml$xtime,D=1-dataml$status,failure.times = x,num.trees = 2000)
  C = predict(fit_rf,newdata = dataest[,1:(ncol(dataml)-2)])$predictions
  
  fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=dataml$a,SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
  p0 = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)
  
  q = S
  t0 = (1:ncol(S))[apply(S,2,min)<1][1]
  for (i in t0:ncol(S)){
    fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=as.numeric(Sh[,i]),SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
    q[,i] = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)
    cat(i,"\t")
  }
  
  # Components
  
  k = length(x)-1
  delta.tau=c(0,x[-1]-x[1:k])
  delta.tau[(1:(k+1))[apply(S,2,max)>exp(-1/(0.1*n))]]=0
  delta.tau[(1:(k+1))[apply(S,2,max)>0.99]]=0
  delta.tau[(1:(k+1))[apply(S*C,2,min)<(1/(0.1*n))]]=0
  
  dN = matrix(0,ncol=k+1,nrow=n)
  dNc = matrix(0,ncol=ncol(dN),nrow=nrow(dN))
  R = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  Rc = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  s = sort.list(dataest$xtime)
  dataest$xtime[s[i]]
  for (i in 1:n){
    dN[s[i],which(dataest$xtime[s[i]]==x)] = dataest$status[s[i]]
    dNc[s[i],which(dataest$xtime[s[i]]==x)] = 1-dataest$status[s[i]]
    if ((i+2) <= ncol(R)){
      R[s[i],(i+2):ncol(R)] = 1-dataest$status[s[i]]
      Rc[s[i],(i+2):ncol(R)] = dataest$status[s[i]]
    }
  }
  
  h = (dataest$a-p0)*(g(S)-q)
  h[is.na(h)] = 0
  h[is.infinite(h)] = 0
  hd = (dataest$a-p0)^2
  gamma = mean(hd*sum(delta.tau))
  theta = mean(apply(t(t(h)*delta.tau),1,sum))/gamma
  IF1 = apply(t(t(h)*delta.tau),1,sum)/gamma
  IF2 = hd*sum(delta.tau)/gamma

  # Augmentation S
  
  hs = -(dataest$a-p0)*g.prime(S)*S
  hs[is.na(hs)] = 0
  hs[is.infinite(hs)] = 0
  hs[,(1:(k+1))[apply(S,2,max)==1]]=0
  
  dM=R
  dM.int=dM
  dLam=-log(S[,1])
  tmp=1/(S[,1]*C[,1])
  dM[,1]=tmp*(dN[,1]-R[,1]*Rc[,1]*dLam)
  dM.int[,1]=dM[,1]
  for (j in 2:ncol(dM)){
    dLam=-log(S[,j])+log(S[,j-1])
    tmp=1/(S[,j]*C[,j])
    dM[,j]=tmp*(dN[,j]-R[,j]*Rc[,j]*dLam)
    dM.int[,j]=dM.int[,j-1]+dM[,j]
  }
  hs = hs*dM.int
  IFS = apply(t(t(hs)*delta.tau),1,sum,na.rm=T)/gamma
  
  list(IF1=IF1,IF2=IF2,IFS=IFS)
}

### With survSuperLearner

estimator<-function(dataml,dataest){
  
  sxtime<-c(0,sort(dataest$xtime))
  x = sort(sxtime) 
  
  event.SL.library <- cens.SL.library <- c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc")
  
  fit_rf <- survSuperLearner(time = dataml$xtime, event = dataml$status, X = dataml[,1:(ncol(dataml)-2)], newX = dataest[,1:(ncol(dataml)-2)], new.times = sxtime, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = TRUE)
  Sh = fit_rf$event.SL.predict
  fit_rf0 = predict.survSuperLearner(fit_rf,newdata=dataest,new.times = sxtime,onlySL=TRUE)
  S = fit_rf0$event.SL.predict
  C = fit_rf0$cens.SL.predict
  C[,ncol(C)] = C[,ncol(C)-1]
  n=nrow(dataest)
  
  fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=dataml$a,SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
  p0 = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)
  
  q = S
  t0 = (1:ncol(S))[apply(S,2,min)<1][1]
  for (i in t0:ncol(S)){
    fit_rf = SuperLearner(X=dataml[,2:(ncol(dataml)-2)],Y=as.numeric(Sh[,i]),SL.library = c("SL.glm","SL.glmnet","SL.ranger","SL.gam"))
    q[,i] = as.numeric(predict(fit_rf,newdata = dataest[,2:(ncol(dataml)-2)])$pred)
    cat(i,"\t")
  }
  
  # Components
  
  k = length(x)-1
  delta.tau=c(0,x[-1]-x[1:k])
  delta.tau[(1:(k+1))[apply(S,2,max)>exp(-1/(0.1*n))]]=0
  delta.tau[(1:(k+1))[apply(S,2,max)>0.99]]=0
  delta.tau[(1:(k+1))[apply(S*C,2,min)<(1/(0.1*n))]]=0
  
  dN = matrix(0,ncol=k+1,nrow=n)
  dNc = matrix(0,ncol=ncol(dN),nrow=nrow(dN))
  R = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  Rc = matrix(1,ncol=ncol(dN),nrow=nrow(dN))
  s = sort.list(dataest$xtime)
  dataest$xtime[s[i]]
  for (i in 1:n){
    dN[s[i],which(dataest$xtime[s[i]]==x)] = dataest$status[s[i]]
    dNc[s[i],which(dataest$xtime[s[i]]==x)] = 1-dataest$status[s[i]]
    if ((i+2) <= ncol(R)){
      R[s[i],(i+2):ncol(R)] = 1-dataest$status[s[i]]
      Rc[s[i],(i+2):ncol(R)] = dataest$status[s[i]]
    }
  }
  
  h = (dataest$a-p0)*(g(S)-q)
  h[is.na(h)] = 0
  h[is.infinite(h)] = 0
  hd = (dataest$a-p0)^2
  gamma = mean(hd*sum(delta.tau))
  theta = mean(apply(t(t(h)*delta.tau),1,sum))/gamma
  IF1 = apply(t(t(h)*delta.tau),1,sum)/gamma
  IF2 = hd*sum(delta.tau)/gamma
  
  # Augmentation S
  
  hs = -(dataest$a-p0)*g.prime(S)*S
  hs[is.na(hs)] = 0
  hs[is.infinite(hs)] = 0
  hs[,(1:(k+1))[apply(S,2,max)==1]]=0
  
  dM=R
  dM.int=dM
  dLam=-log(S[,1])
  tmp=1/(S[,1]*C[,1])
  dM[,1]=tmp*(dN[,1]-R[,1]*Rc[,1]*dLam)
  dM.int[,1]=dM[,1]
  for (j in 2:ncol(dM)){
    dLam=-log(S[,j])+log(S[,j-1])
    tmp=1/(S[,j]*C[,j])
    dM[,j]=tmp*(dN[,j]-R[,j]*Rc[,j]*dLam)
    dM.int[,j]=dM.int[,j-1]+dM[,j]
  }
  hs = hs*dM.int
  IFS = apply(t(t(hs)*delta.tau),1,sum,na.rm=T)/gamma
  
  list(IF1=IF1,IF2=IF2,IFS=IFS)
}
