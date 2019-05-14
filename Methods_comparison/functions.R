# title: Functions
# author: Yiwen (Eva) Wang
# date: May 13th, 2019

## ---------------
## LOGISTIC LASSO
## ---------------

# With linear constraint: $\sum \beta_{j} = 0$ (for $j>1$)
# y: dependent variable, binary, vector of length n
# X: matrix of k covarites (positive values, taxa abundances in counts, proportions, intensities, ...), matrix of dimension n by k

coda_logistic_lasso <- function(y,X,lambda, maxiter=400, maxiter2=50, r=10, 
                                tol=1.e-4, tol2=1.e-6){
  
  library(MASS)   # for the computation of the generalized inverse ginv()
  
  if (!is.numeric(y)){
    y<-as.numeric(y)-1
  }
  
  # Transform initial data 
  ztransformation(X)
  
  #initial values for beta
  #   p<<-ncol(X)   # p: number of covariates after filtering
  #   n<<-nrow(X);  # both defined on ztransformation function
  # beta_ini<-rep(1,p+1)/(p+1)  # uniform values 
  beta_ini<-c(log(mean(y)/(1-mean(y))),rep(1/p,p))   
  # b0 related to mean(y) and uniform values for the other components
  
  #null deviance
  nulldev<-glm(y~1, family=binomial())[[10]]
  
  
  #initialization parameters
  k<-1
  epsilon<-0.1
  beta_ant<-beta_ini
  beta<-beta_ant
  y_ant<-beta_ant
  y_k<-y_ant
  d_k<-y_ant
  dev_explained_ant<-0
  
  
  # Optimization with constraint
  start_iter = Sys.time(); 
  while((abs(epsilon)>tol)&(k<=maxiter)){
    #print(append("loop",k))
    k0<-0
    t_k<-10
    condition<-0.1
    while ((condition>0)&(k0<=maxiter2)){
      k0<-k0+1
      #print(append("k0",k0))
      d_k<-y_ant-t_k*grad_g(y_ant, z, y, n)
      
      # Soft thresholding:
      zproxi<-c(d_k[1],soft_thres(d_k[-1],lambda*t_k))
      
      # Projection:
      zproj<-projection(zproxi,p_c)
      
      # line search condition
      Gt<-(y_ant-zproj)/t_k
      condition<-g(y_ant-t_k*Gt,z, y, n)-g(y_ant,z, y, n)+t_k*t(grad_g(y_ant,z, y, n))%*%Gt-t_k/2*normSqr(Gt)
      #print(append("condition",condition))
      #print(append("t_k",t_k))
      t_k<-t_k/2
    }
    beta<-zproj
    
    y_k<-beta+(k-1)/(k+r-1)*(beta-beta_ant)
    
    
    dev_explained<-1-(nrow(X)*2*g(beta,z, y, n)/nulldev)
    epsilon<-abs(dev_explained-dev_explained_ant)  
    
    
    y_ant<-y_k
    beta_ant<-beta
    dev_explained_ant<-dev_explained
    k<-k+1
    #print(append("epsilon",epsilon))
    
  }
  end_iter = Sys.time(); 
  sprintf("iter time = %f",end_iter-start_iter)
  
  #Projection of the optimal beta to fulfil the constraint sum(beta[j])=0, for j>1 
  indx<-which(abs(zproxi)>tol2)
  if (abs(zproxi[1])>0) indx<-indx[-1]
  c0<-rep(1,(length(indx)))
  c0<-c0/sqrt(normSqr(c0))
  p_c0<-c0%*%t(c0)
  #beta1<-as.numeric(projection(beta[indx],p_c0))
  beta1<-projection(beta[indx],p_c0)
  #beta1
  beta_res<-c(beta[1],rep(0,p))
  beta_res[indx]<-beta1
  
  # print("beta coefficients:")
  # beta_res
  
  #print("taxa with non-zero coeff:")
  # selec<-colnames(z)[abs(beta_res)>0]
  
  #print("beta non-zero coefficients:")
  # beta_res[abs(beta_res)>0]
  
  #print("proportion of explained deviance")
  # dev_explained
  #print("proportion of explained deviance beta_res")
  
  dev_explained_beta_res<-1-(nrow(X)*2*g(beta_res,z, y, n)/nulldev)
  # dev_explained_beta_res
  
  
  results<-list(
    "number of iterations"=k,
    "number of selected taxa"=
      sum(abs(beta_res)>0)-1, 
    "indices of taxa with non-zero coeff"=
      which(abs(beta_res)>0)-1,
    "taxa with non-zero coeff"=
      colnames(z)[abs(beta_res)>0], 
    "beta non-zero coefficients"=
      beta_res[abs(beta_res)>0],
    "proportion of explained deviance"=
      dev_explained_beta_res, 
    "betas"=
      beta_res)
  
  return(results)
  
} # END function coda_logistic_lasso



## ---------------------
## LOGISTIC ELASTIC NET
## ---------------------

##-----------------------------------------------------------------

coda_logistic_elasticNet<-function(y,X,lambda, alpha=0.5, maxiter=1000, maxiter2=50, r=10, 
                                   tol=1.e-4, tol2=1.e-6){
  
  library(MASS)   # for the computation of the generalized inverse ginv()
  
  # Transform initial data 
  ztransformation(X)
  
  #initial values for beta
  #beta_ini<-rep(1,p+1)/(p+1)  # uniform values 
  beta_ini<-c(log(mean(y)/(1-mean(y))),rep(1/p,p))   
  # b0 related to mean(y) and uniform values for the other components
  
  #null deviance
  nulldev<-glm(y~1, family=binomial())[[10]]
  
  
  
  
  
  #initialization parameters
  
  k<-1
  epsilon<-0.1
  beta_ant<-beta_ini
  beta<-beta_ant
  y_ant<-beta_ant
  y_k<-y_ant
  d_k<-y_ant
  dev_explained_ant<-0
  # elastic parameters
  lambda=alpha*lambda;
  gamma=2*(1-alpha)/alpha;
  
  
  # Optimization with constraint
  start_iter = Sys.time()
  while((abs(epsilon)>tol)&(k<=maxiter)){
    #print(append("loop",k))
    k0<-0
    t_k<-10
    condition<-0.1
    while ((condition>0)&(k0<=maxiter2)){
      k0<-k0+1
      #print(append("k0",k0))
      d_k<-y_ant-t_k*grad_g(y_ant, z, y, n)
      
      # Soft thresholding:
      # zproxi<-c(d_k[1],soft_thres(d_k[-1],lambda*t_k))
      # zproxi = zproxi/(1+gamma*lambda);
      zproxi<-c(d_k[1],soft_thres(d_k[-1],lambda*t_k)/(1+gamma*lambda*t_k))
      
      
      # Projection:
      zproj<-projection(zproxi,p_c)
      
      # line search condition
      Gt<-(y_ant-zproj)/t_k
      condition<-g(y_ant-t_k*Gt,z, y, n)-g(y_ant,z, y, n)+t_k*t(grad_g(y_ant,z, y, n))%*%Gt-t_k/2*normSqr(Gt)
      
      #print(append("condition",condition))
      #print(append("t_k",t_k))
      t_k<-t_k/2
    }
    beta<-zproj
    
    y_k<-beta+(k-1)/(k+r-1)*(beta-beta_ant)
    
    
    dev_explained<-1-(nrow(X)*2*g(beta,z, y, n)/nulldev)
    epsilon<-abs(dev_explained-dev_explained_ant)  
    
    
    y_ant<-y_k
    beta_ant<-beta
    dev_explained_ant<-dev_explained
    k<-k+1
    #print(append("epsilon",epsilon))
    
  }
  end_iter = Sys.time()
  sprintf("iter time = %f",end_iter - start_iter)
  
  
  #Projection of the optimal beta to fulfil the constraint sum(beta[j])=0, for j>1 
  indx<-which(abs(zproxi)>tol2)
  if (abs(zproxi[1])>0) indx<-indx[-1]
  c0<-rep(1,(length(indx)))
  c0<-c0/sqrt(normSqr(c0))
  p_c0<-c0%*%t(c0)
  #beta1<-as.numeric(projection(beta[indx],p_c0))
  beta1<-projection(beta[indx],p_c0)
  #beta1
  beta_res<-c(beta[1],rep(0,p))
  beta_res[indx]<-beta1
  
  #print("beta coefficients:")
  # beta_res
  
  #print("taxa with non-zero coeff:")
  # selec<-colnames(z)[abs(beta_res)>0]
  #return(selec)
  
  #print("beta non-zero coefficients:")
  # beta_res[abs(beta_res)>0]
  
  #print("proportion of explained deviance")
  # dev_explained
  #print("proportion of explained deviance beta_res")
  
  dev_explained_beta_res<-1-(nrow(X)*2*g(beta_res,z, y, n)/nulldev)
  # dev_explained_beta_res
  
  
  results<-list("number of iterations"=k,
                "number of selected taxa"=
                  sum(abs(beta_res)>0)-1, 
                "indices of taxa with non-zero coeff"=
                  which(abs(beta_res)>0)-1,
                "name of taxa with non-zero coeff"=
                  colnames(z)[abs(beta_res)>0], 
                "beta non-zero coefficients"=
                  beta_res[abs(beta_res)>0],
                "proportion of explained deviance"=
                  dev_explained_beta_res, 
                "betas"=
                  beta_res)
  
  return(results)
  
} # END function coda_logistic_ElasticNet

## -----------
## rangLambda
## -----------

# It provides a rang of lambda values corresponding to a given number of variables to be selected (numVar).      
# The default initial lambda is lambdaIni=1.

rangLambda <- function(y,X,numVar, lambdaIni =1){
  lambdaB = lambdaIni;
  lambdaA = 0;
  #lambda=0.5*(lambdaB+lambdaA);
  lambda=lambdaIni;
  results <- coda_logistic_lasso(y,X,lambda, maxiter = 100);
  numVarAct = results[[2]];
  if (numVarAct > numVar){
    lambda=lambda+0.5;
    #lambda=lambda+1;
    lambdaB =lambda;
  }else{
    lambdaB =lambda; 
  }
  nIter=1;
  presentLambda = NULL;
  presentnumVar = NULL;
  presentLambda[nIter] = lambda;
  presentnumVar[nIter] = numVarAct;
  numvarA=ncol(X);
  numvarB=numVarAct;
  print(c(lambdaA, lambdaB, numvarA, numvarB))
  diffNvar = abs(numvarB-numvarA);
  while ((diffNvar>1) & (abs(numVarAct-numVar)>0) & (nIter < 6)){
    nIter=nIter+1;
    lambda=0.5*(lambdaB+lambdaA);
    results <- coda_logistic_lasso(y,X,lambda, maxiter = 100);
    numVarAct = results[[2]];
    if (numVarAct < numVar) {
      lambdaB = lambda;
      numvarB = numVarAct;
    }else{
      lambdaA = lambda;
      numvarA = numVarAct;
    }
    presentLambda[nIter] = lambda;
    presentnumVar[nIter] = numVarAct;
    diffNvar = abs(numvarB-numvarA);
    print(c(lambdaA, lambdaB, numvarA, numvarB))
  }
  indx=which(presentnumVar >= 0);
  results = list( 'rang lambdas'=c(lambdaA,lambdaB), 'num selected variables'=c(numvarA,numvarB))
  return(results);
}


## ------------------
## Soft thresholding
## ------------------

# http://www.simonlucey.com/soft-thresholding/

soft_thres<-function(b, lambda){
  x<-rep(0,length(b))
  # Set the threshold
  th = lambda/2; 
  
  #First find elements that are larger than the threshold
  k <- which(b > th)
  x[k] <- b[k] - th 
  
  # Next find elements that are less than abs
  k <-which(abs(b) <= th)
  x[k] <- 0 
  
  # Finally find elements that are less than -th
  k <-which(b < -th)
  x[k] = b[k] + th
  
  return(x)
}

## ----------------
## Other functions
## ----------------

##-----------------------------------------------------------------

norm1<-function(x){
  return(sum(abs(x)))
}

normSqr<-function(x){
  return(sum(x^2))
}


ztransformation<-function(X){
  p<<-ncol(X)   # p: number of covariates after filtering
  n<<-nrow(X);
  # log transformation Z=log(X)
  z<-log(X)
  
  # z=matrix of covariates: add a first column of 1's for beta0
  z<-cbind(rep(1,nrow(z)),z)
  z<-as.matrix(z)
  
  # c=linear constraint sum(betas)=0 (except beta0)
  c<-c(0,rep(1,ncol(X)))
  c<-c/sqrt(sum(c^2))
  c<-as.matrix(c)
  p_c<<-c%*%t(c)
  
  # z transformation (centering) for improvement of optimization
  # this transformation does not affect the estimation since the linear predictor is the same
  z<<-(z-(z%*%p_c))
  colnames(z)<<-c("beta0",colnames(X))
  
  return(z)
  
}

A<-function(x){
  y<-log(1+exp(x))
  return(y)
}


mu_beta<-function(x, Z){
  zetabybeta<-Z%*%x
  res<-rep(1,length(zetabybeta))
  indzb<-which(zetabybeta<=100)
  res[indzb]<-exp(zetabybeta[indzb])/(1+exp(zetabybeta[indzb]))
  return(res)
}

g<-function(x,Z,Y,n){
  res<- (t(Y))%*%Z%*%x
  res<-res-sum(log(1+exp(Z%*%x)))
  res<- as.vector((-res)/n)
  return(res)
}

h<-function(x, lambda){
  lambda*norm1(x)
}

grad_g<-function(x, Z, Y, n){
  res<- (t(Y-mu_beta(x,Z)))%*%Z
  res<- as.vector((-res)/n)
  return(res)
}


# Projection

projection<-function(x, M){
  if (ncol(M)*nrow(M)>0){
    res<-x-ginv(M)%*%M%*%x
  } else {res<-rep(0,length(x))}
  return(res)
}


F<-function(x,t_k,d_k){
  lambda*t_k*norm1(x[-1])+normSqr(d_k-x)/2
}

F2<-function(x,Z,Y,n,lambda){
  g(x,Z,Y,n)+lambda*norm1(x[-1])
}

trapezInteg <- function(x,y) {
  # Compute AUC using trapezoid numerical method 
  n = length(x);
  sumArea = 0;
  for (i in 1:(n-1)){
    h=x[i+1]-x[i];
    if (abs(h) > 1.e-7){
      sumArea = sumArea + 0.5*(y[i+1]+y[i])*h;
    }else{
      sumArea = sumArea;
    }
  }
  return(sumArea)
}


